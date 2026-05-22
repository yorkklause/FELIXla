#!/usr/bin/env bash
# designed by Kai, implemented by codex
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  shapeit_chunks_to_tractor_args.sh --chunks FILE [--chunks FILE ...] [options]
  shapeit_chunks_to_tractor_args.sh --chunks-dir DIR [options]

Convert SHAPEIT5/GLIMPSE-style chunk files to flare_subset_to_tractor_hybrid
argument rows:

  source_phase_vcf  source_flare_vcf  n_ancestries  mac_threshold  out_prefix  region

Options:
  --chunks FILE              Chunk file to read. May be repeated.
  --chunks-dir DIR           Directory with chunks_chr1.txt ... chunks_chr22.txt.
  --phase-template STR       Source phased VCF template. Default: phase/{chrom}.phased.vcf.gz
  --flare-template STR       Source FLARE VCF template. Default: flare/{chrom}.flare.vcf.gz
  --out-prefix-template STR  Output prefix template. Default: hybrid/{chrom}.shapeit4cM.chunk{chunk0}
  --n-ancestries INT         Required converter n_ancestries argument.
  --mac-threshold INT        Required converter mac_threshold argument.
  --region-column core|buffered
                             Use SHAPEIT chunk column 4 (core/no-buffer) or column 3 (buffered).
                             Default: core.
  --chrom-style keep|chr|nochr
                             Keep, add, or remove chr prefix in emitted chrom/region strings.
                             Default: keep.
  --chunk-width INT          Width for {chunk0} and {global_chunk0}. Default: 4.
  --header                   Emit a header row.
  -h, --help                 Show this help.

Template placeholders:
  {chrom} {chunk} {chunk0} {global_chunk} {global_chunk0}
  {shapeit_chunk} {region} {core_region} {buffered_region}

Notes:
  SHAPEIT5 chunk files use column 3 for the buffered region and column 4 for
  the no-buffer/core region. For tractor_hybrid conversion, keep the default
  core region so adjacent chunks do not duplicate variants.
EOF
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

phase_template='phase/{chrom}.phased.vcf.gz'
flare_template='flare/{chrom}.flare.vcf.gz'
out_prefix_template='hybrid/{chrom}.shapeit4cM.chunk{chunk0}'
region_column='core'
chrom_style='keep'
chunk_width=4
header=0
n_ancestries=''
mac_threshold=''
chunks_dir=''
chunks=()

while (($#)); do
  case "$1" in
    --chunks)
      (($# >= 2)) || die '--chunks needs a file'
      chunks+=("$2")
      shift 2
      ;;
    --chunks-dir)
      (($# >= 2)) || die '--chunks-dir needs a directory'
      chunks_dir="$2"
      shift 2
      ;;
    --phase-template)
      (($# >= 2)) || die '--phase-template needs a value'
      phase_template="$2"
      shift 2
      ;;
    --flare-template)
      (($# >= 2)) || die '--flare-template needs a value'
      flare_template="$2"
      shift 2
      ;;
    --out-prefix-template)
      (($# >= 2)) || die '--out-prefix-template needs a value'
      out_prefix_template="$2"
      shift 2
      ;;
    --n-ancestries)
      (($# >= 2)) || die '--n-ancestries needs a value'
      n_ancestries="$2"
      shift 2
      ;;
    --mac-threshold)
      (($# >= 2)) || die '--mac-threshold needs a value'
      mac_threshold="$2"
      shift 2
      ;;
    --region-column)
      (($# >= 2)) || die '--region-column needs core or buffered'
      region_column="$2"
      shift 2
      ;;
    --chrom-style)
      (($# >= 2)) || die '--chrom-style needs keep, chr, or nochr'
      chrom_style="$2"
      shift 2
      ;;
    --chunk-width)
      (($# >= 2)) || die '--chunk-width needs a value'
      chunk_width="$2"
      shift 2
      ;;
    --header)
      header=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

[[ "$region_column" == 'core' || "$region_column" == 'buffered' ]] || die '--region-column must be core or buffered'
[[ "$chrom_style" == 'keep' || "$chrom_style" == 'chr' || "$chrom_style" == 'nochr' ]] || die '--chrom-style must be keep, chr, or nochr'
[[ "$chunk_width" =~ ^[0-9]+$ ]] || die '--chunk-width must be a non-negative integer'
[[ -n "$n_ancestries" ]] || die '--n-ancestries is required'
[[ -n "$mac_threshold" ]] || die '--mac-threshold is required'

if [[ -n "$chunks_dir" ]]; then
  for chr in $(seq 1 22); do
    chunk_file="${chunks_dir%/}/chunks_chr${chr}.txt"
    [[ -f "$chunk_file" ]] || die "missing chunk file: $chunk_file"
    chunks+=("$chunk_file")
  done
fi

((${#chunks[@]} > 0)) || die 'provide --chunks or --chunks-dir'

style_chrom() {
  local chrom="$1"
  local base="${chrom#chr}"
  case "$chrom_style" in
    chr) printf 'chr%s' "$base" ;;
    nochr) printf '%s' "$base" ;;
    *) printf '%s' "$chrom" ;;
  esac
}

style_region() {
  local region="$1"
  [[ "$region" == *:* ]] || die "region lacks chromosome prefix: $region"
  local chrom="${region%%:*}"
  local coords="${region#*:}"
  printf '%s:%s' "$(style_chrom "$chrom")" "$coords"
}

render_template() {
  local value="$1"
  local chrom="$2"
  local chunk="$3"
  local chunk0="$4"
  local global_chunk="$5"
  local global_chunk0="$6"
  local shapeit_chunk="$7"
  local region="$8"
  local core_region="$9"
  local buffered_region="${10}"

  value="${value//\{chrom\}/$chrom}"
  value="${value//\{chunk\}/$chunk}"
  value="${value//\{chunk0\}/$chunk0}"
  value="${value//\{global_chunk\}/$global_chunk}"
  value="${value//\{global_chunk0\}/$global_chunk0}"
  value="${value//\{shapeit_chunk\}/$shapeit_chunk}"
  value="${value//\{region\}/$region}"
  value="${value//\{core_region\}/$core_region}"
  value="${value//\{buffered_region\}/$buffered_region}"
  printf '%s' "$value"
}

if ((header)); then
  printf 'source_phase_vcf\tsource_flare_vcf\tn_ancestries\tmac_threshold\tout_prefix\tregion\n'
fi

global_chunk=0
for chunk_file in "${chunks[@]}"; do
  [[ -f "$chunk_file" ]] || die "chunk file not found: $chunk_file"
  per_chrom_chunk=0
  line_no=0
  while read -r shapeit_chunk chrom buffered_region core_region _; do
    line_no=$((line_no + 1))
    [[ -n "${shapeit_chunk:-}" ]] || continue
    [[ "${shapeit_chunk:0:1}" == '#' ]] && continue
    [[ -n "${chrom:-}" && -n "${buffered_region:-}" && -n "${core_region:-}" ]] || die "$chunk_file:$line_no needs at least four columns"

    per_chrom_chunk=$((per_chrom_chunk + 1))
    global_chunk=$((global_chunk + 1))
    chunk0="$(printf "%0${chunk_width}d" "$per_chrom_chunk")"
    global_chunk0="$(printf "%0${chunk_width}d" "$global_chunk")"
    out_chrom="$(style_chrom "$chrom")"
    out_buffered_region="$(style_region "$buffered_region")"
    out_core_region="$(style_region "$core_region")"
    if [[ "$region_column" == 'buffered' ]]; then
      region="$out_buffered_region"
    else
      region="$out_core_region"
    fi

    source_phase="$(
      render_template "$phase_template" "$out_chrom" "$per_chrom_chunk" "$chunk0" "$global_chunk" "$global_chunk0" "$shapeit_chunk" "$region" "$out_core_region" "$out_buffered_region"
    )"
    source_flare="$(
      render_template "$flare_template" "$out_chrom" "$per_chrom_chunk" "$chunk0" "$global_chunk" "$global_chunk0" "$shapeit_chunk" "$region" "$out_core_region" "$out_buffered_region"
    )"
    out_prefix="$(
      render_template "$out_prefix_template" "$out_chrom" "$per_chrom_chunk" "$chunk0" "$global_chunk" "$global_chunk0" "$shapeit_chunk" "$region" "$out_core_region" "$out_buffered_region"
    )"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$source_phase" \
      "$source_flare" \
      "$n_ancestries" \
      "$mac_threshold" \
      "$out_prefix" \
      "$region"
  done < "$chunk_file"
done
