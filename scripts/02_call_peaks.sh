#!/usr/bin/env bash
# Step 02: call peaks only for validated samples that did not provide peak files.
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
# shellcheck source=../config/tools.sh
source config/tools.sh
JOBS=${1:-results/01_inventory/peak_jobs.tsv}
[[ -s "$JOBS" ]] || { printf 'Missing or empty peak jobs file: %s\n' "$JOBS" >&2; exit 1; }
mkdir -p results/02_peaks/bedgraphs
declare -A IGG_BEDGRAPH_CACHE=()
# Reuse an IgG bedGraph when multiple target libraries reference the same control BAM.

require_tool() {
  command -v "$1" >/dev/null 2>&1 || { printf 'Required tool not found on PATH: %s\n' "$1" >&2; exit 1; }
}

# Split literal tabs without IFS whitespace-collapsing, preserving empty fields.
split_tsv() {
  local rest=$1 field
  TSV_FIELDS=()
  while [[ "$rest" == *$'\t'* ]]; do
    field=${rest%%$'\t'*}
    TSV_FIELDS+=("$field")
    rest=${rest#*$'\t'}
  done
  TSV_FIELDS+=("$rest")
}

# For paired-end data, derive same-chromosome fragments shorter than 1 kb before coverage.
make_bedgraph() (
  local bam=$1 out=$2 tmpdir namesorted fragments genome
  [[ -s "$bam" ]] || { printf 'Missing or empty BAM: %s\n' "$bam" >&2; return 1; }
  tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/h2az-bedgraph.XXXXXX")
  trap 'rm -rf "$tmpdir"' EXIT
  namesorted="$tmpdir/namesorted.bam"
  fragments="$tmpdir/fragments.bed"
  genome="$tmpdir/genome.txt"
  samtools sort -n -@ "$THREADS" -o "$namesorted" "$bam"
  samtools view -H "$bam" | awk -F '\t' '
    $1 == "@SQ" {
      name = seqlen = ""
      for (i = 2; i <= NF; i++) {
        if ($i ~ /^SN:/) name = substr($i, 4)
        if ($i ~ /^LN:/) seqlen = substr($i, 4)
      }
      if (name != "" && seqlen ~ /^[0-9]+$/) print name "\t" seqlen
    }' > "$genome"
  [[ -s "$genome" ]] || { printf 'No @SQ sequence dictionary in BAM: %s\n' "$bam" >&2; return 1; }
  bedtools bamtobed -bedpe -i "$namesorted" \
    | awk 'BEGIN{OFS="\t"} $1 == $4 && $2 < $6 && ($6 - $2) < 1000 {print $1,$2,$6}' \
    | sort -k1,1 -k2,2n -k3,3n > "$fragments"
  bedtools genomecov -bg -i "$fragments" -g "$genome" > "$out"
  [[ -s "$out" ]] || { printf 'Generated empty bedGraph from BAM: %s\n' "$bam" >&2; return 1; }
)

# Header is checked explicitly before processing data records.
IFS= read -r header < "$JOBS" || { printf 'Empty peak jobs file: %s\n' "$JOBS" >&2; exit 1; }
[[ "$header" == $'sample_id\tlayout\tbam\tigg_bam' ]] || { printf 'Unexpected peak_jobs.tsv header: %s\n' "$header" >&2; exit 1; }
line_no=1
while IFS= read -r line <&3 || [[ -n "$line" ]]; do
  ((line_no += 1))
  [[ -n "$line" ]] || continue
  split_tsv "$line"
  [[ ${#TSV_FIELDS[@]} -eq 4 ]] || { printf 'Expected 4 tab-separated fields at %s:%d\n' "$JOBS" "$line_no" >&2; exit 1; }
  sid=${TSV_FIELDS[0]}; layout=${TSV_FIELDS[1]}; bam=${TSV_FIELDS[2]}; igg=${TSV_FIELDS[3]}
  [[ "$sid" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] || { printf 'Invalid sample_id at %s:%d: %s\n' "$JOBS" "$line_no" "$sid" >&2; exit 1; }
  [[ "$layout" == PE || "$layout" == SE ]] || { printf 'Invalid layout for %s: %s\n' "$sid" "$layout" >&2; exit 1; }
  [[ -s "$bam" ]] || { printf 'Missing or empty BAM for %s: %s\n' "$sid" "$bam" >&2; exit 1; }
  # PE libraries use SEACR on fragment coverage; SE libraries use MACS3 on alignments.
  if [[ "$layout" == PE ]]; then
    require_tool samtools
    require_tool bedtools
    require_tool Rscript
    [[ -f "$SEACR" ]] || { printf 'SEACR script not found: %s\n' "$SEACR" >&2; exit 1; }
    SEACR=$(realpath -- "$SEACR")
    target_bg="results/02_peaks/bedgraphs/${sid}.bedgraph"
    make_bedgraph "$bam" "$target_bg"
    prefix="results/02_peaks/$sid"
    peak="${prefix}.stringent.bed"
    rm -f "$peak"
    if [[ -n "$igg" ]]; then
      [[ -s "$igg" ]] || { printf 'Missing or empty IgG BAM for %s: %s\n' "$sid" "$igg" >&2; exit 1; }
      if [[ -z "${IGG_BEDGRAPH_CACHE[$igg]+set}" ]]; then
        igg_key=$(printf '%s' "$igg" | sha256sum | cut -d' ' -f1)
        igg_bg="results/02_peaks/bedgraphs/igg_${igg_key}.bedgraph"
        make_bedgraph "$igg" "$igg_bg"
        IGG_BEDGRAPH_CACHE[$igg]=$igg_bg
      fi
      igg_bg=${IGG_BEDGRAPH_CACHE[$igg]}
      (cd results/02_peaks && bash "$SEACR" "bedgraphs/${sid}.bedgraph" "${igg_bg#results/02_peaks/}" norm stringent "$sid")
    else
      # Without an IgG control, SEACR uses its top 1% signal-block threshold.
      (cd results/02_peaks && bash "$SEACR" "bedgraphs/${sid}.bedgraph" 0.01 non stringent "$sid")
    fi
    
  else
    require_tool "$MACS3"
    mkdir -p results/02_peaks
    peak="results/02_peaks/${sid}_peaks.narrowPeak"
    rm -f "$peak"
    macs_args=(-t "$bam" -f BAM -g mm -q 0.01 -n "$sid" --outdir results/02_peaks)
    if [[ -n "$igg" ]]; then
      [[ -s "$igg" ]] || { printf 'Missing or empty IgG BAM for %s: %s\n' "$sid" "$igg" >&2; exit 1; }
      macs_args+=(-c "$igg")
    else
      # Without an IgG control, MACS3 disables local lambda estimation.
      macs_args+=(--nolambda)
    fi
    "$MACS3" callpeak "${macs_args[@]}"
    
  fi
  [[ -s "$peak" ]] || { printf 'Peak caller produced missing or empty peak file for %s: %s\n' "$sid" "$peak" >&2; exit 1; }
  printf '%s\n' "Peak calling complete: $sid ($layout) -> $peak"
done 3< <(tail -n +2 "$JOBS")
