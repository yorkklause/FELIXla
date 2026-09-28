#!/usr/bin/env bash
# designed by Kai, implemented by codex
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN_DIR="${BIN_DIR:-$ROOT_DIR/bin}"
OUT_DIR="$(mktemp -d "${TMPDIR:-/tmp}/tractor-hybrid-tiny.XXXXXX")"
trap 'rm -rf "$OUT_DIR"' EXIT

find_htslib_tool() {
  local tool="$1"
  local prefix=""
  if command -v "$tool" >/dev/null 2>&1; then
    command -v "$tool"
    return 0
  fi
  prefix="$(pkg-config --variable=prefix htslib 2>/dev/null || true)"
  if [[ -n "$prefix" && -x "$prefix/bin/$tool" ]]; then
    printf '%s\n' "$prefix/bin/$tool"
    return 0
  fi
  return 1
}

python3 - "$BIN_DIR" <<'PY'
import pathlib
import sys

bin_dir = pathlib.Path(sys.argv[1])
unexpected = sorted(
    path.name for path in bin_dir.iterdir()
    if path.name not in {"felixla", "vcf_tbi_chunks"}
)
assert not unexpected, unexpected
assert (bin_dir / "felixla").is_file(), bin_dir
assert (bin_dir / "vcf_tbi_chunks").is_file(), bin_dir
PY

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  2 \
  auto \
  "$OUT_DIR/tiny" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny" \
  "$OUT_DIR/tiny.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/vcf_tbi_chunks" \
  --phase-vcf "$OUT_DIR/tiny.roundtrip.vcf.gz" \
  --chunk-bp 20000 \
  --out "$OUT_DIR/tiny.chunks.tsv" \
  --command-template 'felixla --phase-vcf {phase_vcf_q} --region {region_q} --out chunks/{chrom}.chunk{chrom_chunk0}' \
  --commands-out "$OUT_DIR/tiny.commands.txt" \
  >/dev/null

"$BIN_DIR/vcf_tbi_chunks" \
  --phase-vcf "$OUT_DIR/tiny.roundtrip.vcf.gz" \
  --tbi "$OUT_DIR/tiny.roundtrip.vcf.gz.tbi" \
  --chrom chr2 \
  --chunk-mb 1 \
  --out "$OUT_DIR/tiny.chr2.mb.tsv" \
  >/dev/null

"$BIN_DIR/vcf_tbi_chunks" \
  --phase-vcf "$OUT_DIR/absent.vcf.gz" \
  --tbi "$OUT_DIR/tiny.roundtrip.vcf.gz.tbi" \
  --chunk-mb 1 \
  --out "$OUT_DIR/tiny.detached-index.tsv" \
  >/dev/null

python3 - \
  "$OUT_DIR/tiny.chunks.tsv" \
  "$OUT_DIR/tiny.commands.txt" \
  "$OUT_DIR/tiny.chr2.mb.tsv" \
  "$OUT_DIR/tiny.detached-index.tsv" \
  "$OUT_DIR/tiny.roundtrip.vcf.gz" <<'PY'
import pathlib
import sys

manifest_path = pathlib.Path(sys.argv[1])
commands_path = pathlib.Path(sys.argv[2])
chr2_path = pathlib.Path(sys.argv[3])
detached_path = pathlib.Path(sys.argv[4])
phase_vcf = sys.argv[5]

header = [
    "global_chunk",
    "chrom_chunk",
    "chrom",
    "start",
    "end",
    "region",
    "contig_first_pos",
    "contig_last_pos",
    "contig_records",
]
rows = [line.split("\t") for line in manifest_path.read_text().splitlines()]
assert rows[0] == header, rows[0]
assert rows[1:] == [
    ["1", "1", "chr1", "1", "20001", "chr1:1-20001", "1", "20000", "6"],
    ["2", "1", "chr2", "1", "20001", "chr2:1-20001", "1", "20000", "3"],
], rows

quoted_vcf = "'" + phase_vcf + "'"
assert commands_path.read_text().splitlines() == [
    f"felixla --phase-vcf {quoted_vcf} --region 'chr1:1-20001' --out chunks/chr1.chunk0001",
    f"felixla --phase-vcf {quoted_vcf} --region 'chr2:1-20001' --out chunks/chr2.chunk0001",
]

chr2_rows = [line.split("\t") for line in chr2_path.read_text().splitlines()]
assert chr2_rows[0] == header, chr2_rows[0]
assert chr2_rows[1:] == [
    ["1", "1", "chr2", "1", "1000001", "chr2:1-1000001", "1", "1000000", "3"],
], chr2_rows

detached_rows = [line.split("\t") for line in detached_path.read_text().splitlines()]
assert detached_rows == [
    header,
    ["1", "1", "chr1", "1", "1000001", "chr1:1-1000001", "1", "1000000", "6"],
    ["2", "1", "chr2", "1", "1000001", "chr2:1-1000001", "1", "1000000", "3"],
], detached_rows
PY

python3 - "$ROOT_DIR/tests" "$OUT_DIR/tiny" <<'QUERYCHECK'
import pathlib
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

# Dosages read out of the packed files rather than asked of FELIXla, so a
# bug in the reader cannot cancel a matching bug in the writer.
dosages = felixla_check.read_dosages(pathlib.Path(sys.argv[2]))
assert dosages[("chr1", 160, "A", "T")] == [(2, [1, 1]), (0, [0, 0])], dosages
assert dosages[("chr1", 100, "A", "G")] == [(0, [0, 0]), (1, [0, 1])], dosages
QUERYCHECK

"$BIN_DIR/felixla" --version >"$OUT_DIR/felixla.version.txt"
grep -q "FELIXla CLI v0" "$OUT_DIR/felixla.version.txt"

"$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --export-felixla \
  --out \
  "$OUT_DIR/tiny.felixla_cli" >/dev/null

"$BIN_DIR/felixla" \
  --felixla \
  "$OUT_DIR/tiny.felixla_cli" \
  --export-vcf \
  --out \
  "$OUT_DIR/tiny.felixla_cli.roundtrip.vcf.gz" >/dev/null

python3 - "$ROOT_DIR/tests" "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/tiny.felixla_cli.roundtrip.vcf.gz" <<'ROUNDTRIP'
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

diffs = felixla_check.compare_vcfs(sys.argv[2], sys.argv[3], split_lhs=True)
assert not diffs, "round trip differs:\n  " + "\n  ".join(diffs)
ROUNDTRIP

cat >"$OUT_DIR/keep.samples" <<'EOF'
s2
EOF

cat >"$OUT_DIR/extract.sites.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ignored_sample
chr1	100	extract_id_is_ignored_1	A	C,G	.	PASS	.	GT	0|1
chr1	160	extract_id_is_ignored_2	A	T	.	PASS	.	GT	0|0
EOF

python3 - "$OUT_DIR/extract.sites.vcf" "$OUT_DIR/extract.sites.vcf.gz" <<'PY'
import gzip
import pathlib
import sys

gzip.open(sys.argv[2], "wb").write(pathlib.Path(sys.argv[1]).read_bytes())
PY

"$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --keep "$OUT_DIR/keep.samples" \
  --extract "$OUT_DIR/extract.sites.vcf.gz" \
  --export-felixla \
  --out \
  "$OUT_DIR/tiny.keep_extract" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.keep_extract" \
  "$OUT_DIR/tiny.keep_extract.roundtrip.vcf.gz" >/dev/null

python3 - "$OUT_DIR/tiny.keep_extract" <<'PY'
import gzip
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
meta = pathlib.Path(str(prefix) + ".meta").read_text()
log = pathlib.Path(str(prefix) + ".log").read_text()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
assert "n_samples\t1" in meta, meta
# The filters applied are run facts, so they are recorded in the log.
assert "\nkeep " in log, log
assert "\nextract " in log, log
assert log.rstrip().splitlines()[-1].startswith("completed "), log
assert samples == ["s2"], samples

with gzip.open(str(prefix) + ".roundtrip.vcf.gz", "rt") as fh:
    lines = [line.rstrip() for line in fh]
header = [line for line in lines if line.startswith("#CHROM")][0].split("\t")
assert header[9:] == ["s2"], header
rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
assert [(r[0], r[1], r[2], r[3], r[4], r[9]) for r in rows] == [
    ("chr1", "100", "multi_A_C", "A", "C", "0|0"),
    ("chr1", "100", "multi_A_G", "A", "G", "1|0"),
    ("chr1", "160", "common_A_T", "A", "T", "0|0"),
], rows
PY

BGZIP_BIN="$(find_htslib_tool bgzip || true)"
TABIX_BIN="$(find_htslib_tool tabix || true)"
if [[ -n "$BGZIP_BIN" && -n "$TABIX_BIN" ]]; then
  "$BGZIP_BIN" -c "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
    >"$OUT_DIR/tiny.region_seek.genotypes.vcf.gz"
  "$BGZIP_BIN" -c "$ROOT_DIR/testdata/tiny.flare.region_seek.vcf" \
    >"$OUT_DIR/tiny.region_seek.flare.vcf.gz"
  "$TABIX_BIN" -f -p vcf "$OUT_DIR/tiny.region_seek.genotypes.vcf.gz"
  "$TABIX_BIN" -f -p vcf "$OUT_DIR/tiny.region_seek.flare.vcf.gz"

  "$BIN_DIR/felixla" \
    --phase-vcf "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
    --flare-vcf "$ROOT_DIR/testdata/tiny.flare.region_seek.vcf" \
      --region chr1:150-161 \
    --extract "$OUT_DIR/extract.sites.vcf" \
    --export-felixla \
    --out "$OUT_DIR/tiny.region_seek.sequential" \
    >/dev/null 2>"$OUT_DIR/tiny.region_seek.sequential.err"

  "$BIN_DIR/felixla" \
    --phase-vcf "$OUT_DIR/tiny.region_seek.genotypes.vcf.gz" \
    --flare-vcf "$OUT_DIR/tiny.region_seek.flare.vcf.gz" \
      --region chr1:150-161 \
    --extract "$OUT_DIR/extract.sites.vcf" \
    --export-felixla \
    --out "$OUT_DIR/tiny.region_seek.indexed" \
    >/dev/null 2>"$OUT_DIR/tiny.region_seek.indexed.err"

  grep -q 'Using indexed FLARE reader from region start: chr1:150-2147483647.' \
    "$OUT_DIR/tiny.region_seek.indexed.err"
  for suffix in \
    .common.geno.bin .common.variant.mks .common.variant.idx \
    .rare.carrier.bin .rare.variant.mks .rare.variant.idx \
    .ancblock.bin .ancblock.mks .ancblock.idx .samples; do
    cmp "$OUT_DIR/tiny.region_seek.sequential${suffix}" \
      "$OUT_DIR/tiny.region_seek.indexed${suffix}"
  done

  cat >"$OUT_DIR/region_tail.extract.pvar" <<'EOF'
#CHROM	POS	ID	REF	ALT
chr1	250	ignored	A	T
EOF

  "$BIN_DIR/felixla" \
    --phase-vcf "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
    --flare-vcf "$ROOT_DIR/testdata/tiny.flare.region_seek.vcf" \
      --region chr1:250-251 \
    --extract "$OUT_DIR/region_tail.extract.pvar" \
    --export-felixla \
    --out "$OUT_DIR/tiny.region_tail.sequential" \
    >/dev/null 2>"$OUT_DIR/tiny.region_tail.sequential.err"

  "$BIN_DIR/felixla" \
    --phase-vcf "$OUT_DIR/tiny.region_seek.genotypes.vcf.gz" \
    --flare-vcf "$OUT_DIR/tiny.region_seek.flare.vcf.gz" \
      --region chr1:250-251 \
    --extract "$OUT_DIR/region_tail.extract.pvar" \
    --export-felixla \
    --out "$OUT_DIR/tiny.region_tail.indexed" \
    >/dev/null 2>"$OUT_DIR/tiny.region_tail.indexed.err"

  grep -q 'Using preceding indexed FLARE window: chr1:186-249.' \
    "$OUT_DIR/tiny.region_tail.indexed.err"
  for suffix in \
    .common.geno.bin .common.variant.mks .common.variant.idx \
    .rare.carrier.bin .rare.variant.mks .rare.variant.idx \
    .ancblock.bin .ancblock.mks .ancblock.idx .samples; do
    cmp "$OUT_DIR/tiny.region_tail.sequential${suffix}" \
      "$OUT_DIR/tiny.region_tail.indexed${suffix}"
  done
else
  echo "WARNING: bgzip/tabix unavailable; skipping indexed FLARE region-start test." >&2
fi

cat >"$OUT_DIR/extract.bad_ref.vcf" <<'EOF'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	multi	T	G	.	PASS	.
EOF

if "$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --extract "$OUT_DIR/extract.bad_ref.vcf" \
  --export-felixla \
  --out \
  "$OUT_DIR/tiny.bad_extract_ref" >/dev/null 2>"$OUT_DIR/tiny.bad_extract_ref.err"; then
  echo "expected bad REF in --extract fixture to fail" >&2
  exit 1
fi

grep -q "REF mismatch" "$OUT_DIR/tiny.bad_extract_ref.err"

"$BIN_DIR/felixla" \
  --felixla \
  "$OUT_DIR/tiny.felixla_cli" \
  --region chr1:100-161 \
  --export-felixla \
  --out \
  "$OUT_DIR/tiny.felixla_cli.chr1_100_160" >/dev/null


python3 - "$ROOT_DIR/tests" "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/tiny.roundtrip.vcf.gz" <<'ROUNDTRIP'
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

diffs = felixla_check.compare_vcfs(sys.argv[2], sys.argv[3], split_lhs=True)
assert not diffs, "round trip differs:\n  " + "\n  ".join(diffs)
ROUNDTRIP

"$BIN_DIR/felixla" \
  extract \
  "$OUT_DIR/tiny" \
  chr1:100-161 \
  "$OUT_DIR/tiny.chr1_100_160" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.chr1_100_160" \
  "$OUT_DIR/tiny.chr1_100_160.roundtrip.vcf.gz" >/dev/null

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  2 \
  auto \
  "$OUT_DIR/tiny.direct_region" \
  chr1:100-161 >/dev/null

grep -q $'selected_region\tchr1:100-161' "$OUT_DIR/tiny.direct_region.meta"

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.direct_region" \
  "$OUT_DIR/tiny.direct_region.roundtrip.vcf.gz" >/dev/null

python3 - "$ROOT_DIR/tests" \
  "$OUT_DIR/tiny.chr1_100_160.roundtrip.vcf.gz" \
  "$OUT_DIR/tiny.direct_region.roundtrip.vcf.gz" <<'REGIONCHECK'
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

# Extracting a region from a prefix must match packing that region directly.
diffs = felixla_check.compare_vcfs(sys.argv[2], sys.argv[3])
assert not diffs, "region extract differs:\n  " + "\n  ".join(diffs)
REGIONCHECK

grep -v '^##contig=' "$ROOT_DIR/testdata/tiny.genotypes.vcf" >"$OUT_DIR/tiny.genotypes.no_contig.vcf"
grep -v '^##contig=' "$ROOT_DIR/testdata/tiny.flare.vcf" >"$OUT_DIR/tiny.flare.no_contig.vcf"

"$BIN_DIR/felixla" \
  from-flare \
  "$OUT_DIR/tiny.genotypes.no_contig.vcf" \
  "$OUT_DIR/tiny.flare.no_contig.vcf" \
  2 \
  auto \
  "$OUT_DIR/tiny.no_contig" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.no_contig" \
  "$OUT_DIR/tiny.no_contig.roundtrip.vcf.gz" >/dev/null

python3 - "$ROOT_DIR/tests" \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/tiny.no_contig.roundtrip.vcf.gz" <<'ROUNDTRIP_NO_CONTIG'
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

diffs = felixla_check.compare_vcfs(sys.argv[2], sys.argv[3], split_lhs=True)
assert not diffs, "no_contig round trip differs:\n  " + "\n  ".join(diffs)
ROUNDTRIP_NO_CONTIG

"$BIN_DIR/felixla" \
  from-flare \
  "$OUT_DIR/tiny.genotypes.no_contig.vcf" \
  "$OUT_DIR/tiny.flare.no_contig.vcf" \
  2 \
  auto \
  "$OUT_DIR/tiny.no_contig_region" \
  chr1:100-161 >/dev/null

grep -q $'selected_region\tchr1:100-161' "$OUT_DIR/tiny.no_contig_region.meta"

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/tiny.no_contig_region" \
  "$OUT_DIR/tiny.no_contig_region.roundtrip.vcf.gz" >/dev/null

python3 - "$ROOT_DIR/tests" \
  "$OUT_DIR/tiny.direct_region.roundtrip.vcf.gz" \
  "$OUT_DIR/tiny.no_contig_region.roundtrip.vcf.gz" <<'ROUNDTRIP_NO_CONTIG_REGION'
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

diffs = felixla_check.compare_vcfs(sys.argv[2], sys.argv[3], split_lhs=False)
assert not diffs, "no_contig_region round trip differs:\n  " + "\n  ".join(diffs)
ROUNDTRIP_NO_CONTIG_REGION

"$BIN_DIR/felixla" \
  --phase-vcf \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --rfmix-msp \
  "$ROOT_DIR/testdata/tiny.rfmix.msp.tsv" \
  --export-felixla \
  --out \
  "$OUT_DIR/rfmix" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/rfmix" \
  "$OUT_DIR/rfmix.roundtrip.vcf.gz" >/dev/null

python3 - "$ROOT_DIR/tests" \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$OUT_DIR/rfmix.roundtrip.vcf.gz" <<'ROUNDTRIP_RFMIX'
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

diffs = felixla_check.compare_vcfs(sys.argv[2], sys.argv[3], split_lhs=True)
assert not diffs, "rfmix round trip differs:\n  " + "\n  ".join(diffs)
ROUNDTRIP_RFMIX

"$BIN_DIR/felixla" \
  --tractor-dosage-vcf \
  "$ROOT_DIR/testdata/tiny.tractor_dosage.vcf" \
  --export-felixla \
  --out \
  "$OUT_DIR/dosage" >/dev/null

"$BIN_DIR/felixla" \
  to-vcf \
  "$OUT_DIR/dosage" \
  "$OUT_DIR/dosage.roundtrip.vcf.gz" >/dev/null

python3 - "$OUT_DIR/tiny" <<'PY'
import gzip
import pathlib
import struct
import sys

prefix = pathlib.Path(sys.argv[1])

def read_string(data, off):
    (n,) = struct.unpack_from("<I", data, off)
    off += 4
    value = data[off:off + n].decode()
    off += n
    return value, off

def read_common_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRCMMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        mks_offset = off
        common_index, global_index = struct.unpack_from("<QI", data, off)
        off += 12
        chrom, off = read_string(data, off)
        (pos,) = struct.unpack_from("<q", data, off)
        off += 8
        vid, off = read_string(data, off)
        ref, off = read_string(data, off)
        alt, off = read_string(data, off)
        alt_index, block_id = struct.unpack_from("<II", data, off)
        off += 8
        geno_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        mac, = struct.unpack_from("<I", data, off)
        off += 4
        records.append({
            "mks_offset": mks_offset,
            "common_index": common_index,
            "global_index": global_index,
            "chr": chrom,
            "pos": pos,
            "id": vid,
            "ref": ref,
            "alt": alt,
            "alt_index": alt_index,
            "block_id": block_id,
            "geno_offset": geno_offset,
            "mac": mac,
        })
    return records

def read_rare_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRRAMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        mks_offset = off
        rare_index, global_index = struct.unpack_from("<QI", data, off)
        off += 12
        chrom, off = read_string(data, off)
        (pos,) = struct.unpack_from("<q", data, off)
        off += 8
        vid, off = read_string(data, off)
        ref, off = read_string(data, off)
        alt, off = read_string(data, off)
        (alt_index,) = struct.unpack_from("<I", data, off)
        off += 4
        carrier_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        n_carriers, mac = struct.unpack_from("<II", data, off)
        off += 8
        records.append({
            "mks_offset": mks_offset,
            "rare_index": rare_index,
            "global_index": global_index,
            "chr": chrom,
            "pos": pos,
            "id": vid,
            "ref": ref,
            "alt": alt,
            "alt_index": alt_index,
            "carrier_offset": carrier_offset,
            "n_carriers": n_carriers,
            "mac": mac,
        })
    return records

def read_anc_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRANMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        mks_offset = off
        (block_id,) = struct.unpack_from("<I", data, off)
        off += 4
        chrom, off = read_string(data, off)
        start, end = struct.unpack_from("<qq", data, off)
        off += 16
        anc_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        records.append({
            "mks_offset": mks_offset,
            "block_id": block_id,
            "chr": chrom,
            "start": start,
            "end": end,
            "anc_offset": anc_offset,
        })
    return records

def read_common_offset_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRCMIDX2", data[:8]
    return [
        {
            "common_index": struct.unpack_from("<Q", data, off)[0],
            "global_index": struct.unpack_from("<I", data, off + 8)[0],
            "mks_offset": struct.unpack_from("<Q", data, off + 12)[0],
            "geno_offset": struct.unpack_from("<Q", data, off + 20)[0],
        }
        for off in range(8, len(data), 28)
    ]

def read_rare_offset_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRRAIDX2", data[:8]
    return [
        {
            "rare_index": struct.unpack_from("<Q", data, off)[0],
            "global_index": struct.unpack_from("<I", data, off + 8)[0],
            "mks_offset": struct.unpack_from("<Q", data, off + 12)[0],
            "carrier_offset": struct.unpack_from("<Q", data, off + 20)[0],
            "n_carriers": struct.unpack_from("<I", data, off + 28)[0],
        }
        for off in range(8, len(data), 32)
    ]

def read_anc_offset_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRANIDX2", data[:8]
    return [
        {
            "block_id": struct.unpack_from("<I", data, off)[0],
            "mks_offset": struct.unpack_from("<Q", data, off + 4)[0],
            "anc_offset": struct.unpack_from("<Q", data, off + 12)[0],
        }
        for off in range(8, len(data), 20)
    ]

common_idx = read_common_idx(str(prefix) + ".common.variant.mks")
rare_idx = read_rare_idx(str(prefix) + ".rare.variant.mks")
anc_idx = read_anc_idx(str(prefix) + ".ancblock.mks")
common_offsets = read_common_offset_idx(str(prefix) + ".common.variant.idx")
rare_offsets = read_rare_offset_idx(str(prefix) + ".rare.variant.idx")
anc_offsets = read_anc_offset_idx(str(prefix) + ".ancblock.idx")
meta = pathlib.Path(str(prefix) + ".meta").read_text().strip().splitlines()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
roundtrip_path = pathlib.Path(str(prefix) + ".roundtrip.vcf.gz")
assert roundtrip_path.exists(), roundtrip_path
assert pathlib.Path(str(roundtrip_path) + ".tbi").exists(), str(roundtrip_path) + ".tbi"
with gzip.open(roundtrip_path, "rt") as fh:
    roundtrip_vcf = fh.read().strip().splitlines()

assert len(common_idx) == 1, common_idx
assert len(rare_idx) == 8, rare_idx
assert len(anc_idx) == 2, anc_idx
assert common_offsets == [
    {
        "common_index": r["common_index"],
        "global_index": r["global_index"],
        "mks_offset": r["mks_offset"],
        "geno_offset": r["geno_offset"],
    }
    for r in common_idx
], common_offsets
assert rare_offsets == [
    {
        "rare_index": r["rare_index"],
        "global_index": r["global_index"],
        "mks_offset": r["mks_offset"],
        "carrier_offset": r["carrier_offset"],
        "n_carriers": r["n_carriers"],
    }
    for r in rare_idx
], rare_offsets
assert anc_offsets == [
    {
        "block_id": r["block_id"],
        "mks_offset": r["mks_offset"],
        "anc_offset": r["anc_offset"],
    }
    for r in anc_idx
], anc_offsets
assert "n_samples\t2" in meta, meta
assert "n_words\t1" in meta, meta
assert samples == ["s1", "s2"], samples

assert {k: common_idx[0][k] for k in ("global_index", "chr", "pos", "id", "ref", "alt", "alt_index", "block_id")} == {
    "global_index": 4,
    "chr": "chr1",
    "pos": 160,
    "id": "common_A_T",
    "ref": "A",
    "alt": "T",
    "alt_index": 1,
    "block_id": 0,
}, common_idx[0]

rare_globals = [record["global_index"] for record in rare_idx]
assert rare_globals == [0, 1, 2, 3, 5, 6, 7, 8], rare_globals

assert {k: anc_idx[0][k] for k in ("block_id", "chr", "start", "end")} == {
    "block_id": 0, "chr": "chr1", "start": 1, "end": 250
}, anc_idx[0]
assert {k: anc_idx[1][k] for k in ("block_id", "chr", "start", "end")} == {
    "block_id": 1, "chr": "chr2", "start": 1, "end": 150
}, anc_idx[1]

common_bin = pathlib.Path(str(prefix) + ".common.geno.bin").read_bytes()
assert struct.unpack("<Q", common_bin) == (0b0011,), common_bin

rare_bin = pathlib.Path(str(prefix) + ".rare.carrier.bin").read_bytes()
records = [
    struct.unpack("<II", rare_bin[i:i + 8])
    for i in range(0, len(rare_bin), 8)
]

expected = [
    (0, (1 << 27) | 1),
    (1, (1 << 27) | 1),
    (2, (1 << 27) | 2),
    (3, 0),
    (5, 3),
    (6, (1 << 27) | 0),
    (7, 2),
    (8, (1 << 27) | 1),
]
assert records == expected, records

vcf_records = [
    line.split("\t")
    for line in roundtrip_vcf
    if line and not line.startswith("#")
]
assert len(vcf_records) == 9, vcf_records
assert [(r[0], r[1], r[2], r[3], r[4], r[9], r[10]) for r in vcf_records] == [
    ("chr1", "50", "pre_A_T", "A", "T", "0|1", "0|0"),
    ("chr1", "100", "multi_A_C", "A", "C", "0|1", "0|0"),
    ("chr1", "100", "multi_A_G", "A", "G", "0|0", "1|0"),
    ("chr1", "150", "rare_A_T", "A", "T", "1|0", "0|0"),
    ("chr1", "160", "common_A_T", "A", "T", "1|1", "0|0"),
    ("chr1", "250", "tail1_A_T", "A", "T", "0|0", "0|1"),
    ("chr2", "50", "should_skip_A_T", "A", "T", "1|0", "0|0"),
    ("chr2", "100", "chr2var_A_T", "A", "T", "0|0", "1|0"),
    ("chr2", "150", "tail2_A_T", "A", "T", "0|1", "0|0"),
], vcf_records
PY

python3 - "$OUT_DIR/tiny.chr1_100_160" <<'PY'
import gzip
import pathlib
import struct
import sys

prefix = pathlib.Path(sys.argv[1])
meta = pathlib.Path(str(prefix) + ".meta").read_text()
assert "extracted_region\tchr1:100-161" in meta, meta
assert "source_hybrid_prefix" in meta, meta

with gzip.open(str(prefix) + ".roundtrip.vcf.gz", "rt") as fh:
    rows = [line.rstrip().split("\t") for line in fh if not line.startswith("#")]
assert [(r[0], r[1], r[2], r[3], r[4], r[9], r[10]) for r in rows] == [
    ("chr1", "100", "multi_A_C", "A", "C", "0|1", "0|0"),
    ("chr1", "100", "multi_A_G", "A", "G", "0|0", "1|0"),
    ("chr1", "150", "rare_A_T", "A", "T", "1|0", "0|0"),
    ("chr1", "160", "common_A_T", "A", "T", "1|1", "0|0"),
], rows

rare_bin = pathlib.Path(str(prefix) + ".rare.carrier.bin").read_bytes()
records = [
    struct.unpack("<II", rare_bin[i:i + 8])
    for i in range(0, len(rare_bin), 8)
]
assert [r[0] for r in records] == [0, 1, 2], records
PY

python3 - "$OUT_DIR/dosage" <<'PY'
import gzip
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
log = pathlib.Path(str(prefix) + ".log").read_text()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
assert "\ntractor-dosage-vcf " in log, log
assert "\ndosage-source-note " in log, log
assert log.rstrip().splitlines()[-1].startswith("completed "), log
assert samples == ["s1", "s2"], samples
with gzip.open(str(prefix) + ".roundtrip.vcf.gz", "rt") as fh:
    rows = [line.rstrip().split("\t") for line in fh if not line.startswith("#")]
assert [(r[0], r[1], r[2], r[3], r[4], r[9], r[10]) for r in rows] == [
    ("chr1", "100", "v1", "A", "G", "1|0", "0|1"),
    ("chr1", "160", "v2", "C", "T", "1|1", "1|0"),
    ("chr1", "220", "v3", "G", "A", "1|0", "1|0"),
], rows
PY

python3 - "$OUT_DIR/rfmix" <<'PY'
import pathlib
import sys

prefix = pathlib.Path(sys.argv[1])
meta = pathlib.Path(str(prefix) + ".meta").read_text()
log = pathlib.Path(str(prefix) + ".log").read_text()
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
assert "\nrfmix-msp " in log, log
assert "\nrfmix-msp-interval-note " in log, log
assert log.rstrip().splitlines()[-1].startswith("completed "), log
# The subpopulation order stays in .meta: it is what an ancestry index means.
assert "rfmix_subpopulation_order_codes" in meta, meta
assert samples == ["s1", "s2"], samples
PY

if "$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.missing_gt.vcf" \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  2 \
  auto \
  "$OUT_DIR/missing_gt" >/dev/null 2>"$OUT_DIR/missing_gt.err"; then
  echo "expected missing GT fixture to fail" >&2
  exit 1
fi

grep -q "missing genotype" "$OUT_DIR/missing_gt.err"

if "$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.missing_flare.vcf" \
  2 \
  auto \
  "$OUT_DIR/missing_flare" >/dev/null 2>"$OUT_DIR/missing_flare.err"; then
  echo "expected missing FLARE fixture to fail" >&2
  exit 1
fi

grep -q "missing FORMAT/AN1" "$OUT_DIR/missing_flare.err"

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.flare_swapped_samples.vcf" \
  2 \
  auto \
  "$OUT_DIR/swapped_samples" >/dev/null 2>"$OUT_DIR/swapped_samples.err"

grep -q "Using genotype/FLARE sample intersection" "$OUT_DIR/swapped_samples.err"

python3 - "$ROOT_DIR/tests" "$OUT_DIR/swapped_samples" <<'SWAPPED'
import pathlib
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

prefix = pathlib.Path(sys.argv[2])
samples = pathlib.Path(str(prefix) + ".samples").read_text().strip().splitlines()
assert samples == ["s1", "s2"], samples

# Samples given in a different order by the FLARE file must still be packed
# against the right haplotypes.
dosages = felixla_check.read_dosages(prefix)
assert dosages[("chr1", 100, "A", "C")] == [(1, [0, 1]), (0, [0, 0])], dosages
SWAPPED

"$BIN_DIR/felixla" \
  from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.duplicate_flare.vcf" \
  2 \
  auto \
  "$OUT_DIR/duplicate_flare" >/dev/null

python3 - "$OUT_DIR/duplicate_flare" <<'PY'
import pathlib
import struct
import sys

prefix = pathlib.Path(sys.argv[1])
rare_bin = pathlib.Path(str(prefix) + ".rare.carrier.bin").read_bytes()

def read_string(data, off):
    (n,) = struct.unpack_from("<I", data, off)
    off += 4
    value = data[off:off + n].decode()
    off += n
    return value, off

def read_rare_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRRAMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        rare_index, global_index = struct.unpack_from("<QI", data, off)
        off += 12
        chrom, off = read_string(data, off)
        (pos,) = struct.unpack_from("<q", data, off)
        off += 8
        vid, off = read_string(data, off)
        ref, off = read_string(data, off)
        alt, off = read_string(data, off)
        (alt_index,) = struct.unpack_from("<I", data, off)
        off += 4
        carrier_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        n_carriers, mac = struct.unpack_from("<II", data, off)
        off += 8
        records.append({"global_index": global_index, "n_carriers": n_carriers, "mac": mac})
    return records

def read_anc_idx(path):
    data = pathlib.Path(path).read_bytes()
    assert data[:8] == b"TRANMKS1", data[:8]
    off = 8
    records = []
    while off < len(data):
        (block_id,) = struct.unpack_from("<I", data, off)
        off += 4
        chrom, off = read_string(data, off)
        start, end = struct.unpack_from("<qq", data, off)
        off += 16
        anc_offset, = struct.unpack_from("<Q", data, off)
        off += 8
        records.append({"block_id": block_id, "chr": chrom, "start": start, "end": end})
    return records

anc_idx = read_anc_idx(str(prefix) + ".ancblock.mks")
rare_idx = read_rare_idx(str(prefix) + ".rare.variant.mks")

assert len(anc_idx) == 1, anc_idx
assert {k: anc_idx[0][k] for k in ("block_id", "chr", "start", "end")} == {
    "block_id": 0, "chr": "chr1", "start": 1, "end": 250
}, anc_idx[0]
assert [record["global_index"] for record in rare_idx] == [0, 1, 2, 3, 5], rare_idx

records = [
    struct.unpack("<II", rare_bin[i:i + 8])
    for i in range(0, len(rare_bin), 8)
]
assert records == [
    (0, 1),
    (1, 1),
    (2, 2),
    (3, (1 << 27) | 0),
    (5, (1 << 27) | 3),
], records
PY


# The three ancestry views of the prefix, against values worked out by hand
# from testdata/tiny.flare.vcf: s1 is AFR|EUR then EUR|EUR, s2 is EUR|AFR then
# AFR|AFR, over blocks of 250 and 150 bases.
"$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" --export-lai --out "$OUT_DIR/tiny"
"$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" --export-global-admixture --out "$OUT_DIR/tiny"
"$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" --export-local-admixture --out "$OUT_DIR/tiny"

python3 - "$ROOT_DIR/tests" "$OUT_DIR/tiny" <<'ANCESTRY_EXPORTS'
import pathlib
import sys

sys.path.insert(0, sys.argv[1])
import felixla_check

prefix = pathlib.Path(sys.argv[2])

names, columns, rows = felixla_check.read_lai(str(prefix) + ".lai.gz")
assert names == [("ANC1", "AFR"), ("ANC2", "EUR")], names
assert columns == ["s1_1", "s1_2", "s2_1", "s2_2"], columns
# Codes are one-based, leaving zero for a haplotype no caller labelled.
assert rows == [
    ("chr1", 1, 250, [1, 2, 2, 1]),
    ("chr2", 1, 150, [2, 2, 1, 1]),
], rows

header, table = felixla_check.read_tsv(str(prefix) + ".global.admixture.tsv")
assert header == ["#ID", "AFR", "EUR"], header
# s1 is AFR over 250 of its 800 haplotype bases, EUR over the other 550.
assert table == [
    ["s1", "0.312500", "0.687500"],
    ["s2", "0.687500", "0.312500"],
], table

header, table = felixla_check.read_tsv(str(prefix) + ".local.admixture.tsv")
assert header == ["#CHR", "START", "END", "AFR", "EUR"], header
assert table == [
    ["chr1", "1", "250", "0.500000", "0.500000"],
    ["chr2", "1", "150", "0.500000", "0.500000"],
], table

# Every --out spelling of the same file lands in the same place, so a run is
# never left looking for its output under a name nothing wrote.
ANCESTRY_EXPORTS

for spelling in "$OUT_DIR/alias.lai" "$OUT_DIR/alias.gz" "$OUT_DIR/alias"; do
  "$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" --export-lai --out "$spelling" 2>/dev/null
  test -f "$OUT_DIR/alias.lai.gz"
  rm -f "$OUT_DIR/alias.lai.gz"
done


# --threads reaches only the two paths that can use it. Anywhere else it used
# to be accepted and ignored, which reads as a tuning knob that does not work.
"$BIN_DIR/felixla" --phase-vcf "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --export-felixla --out "$OUT_DIR/threaded" --threads 2 >/dev/null
"$BIN_DIR/felixla" --felixla "$OUT_DIR/threaded" --export-lai \
  --out "$OUT_DIR/threaded" --threads 2 >/dev/null 2>&1

printf '%s\n' "$OUT_DIR/threaded" >"$OUT_DIR/threaded.list"
for inert in \
  "--felixla $OUT_DIR/threaded --export-vcf --out $OUT_DIR/nothreads" \
  "--felixla $OUT_DIR/threaded --export-global-admixture --out $OUT_DIR/nothreads" \
  "--felixla $OUT_DIR/threaded --region chr1:1-1000 --export-felixla --out $OUT_DIR/nothreads" \
  "--merge-list $OUT_DIR/threaded.list --export-felixla --out $OUT_DIR/nothreads" \
  "--phase-vcf $ROOT_DIR/testdata/tiny.genotypes.vcf --rfmix-msp $ROOT_DIR/testdata/tiny.rfmix.msp.tsv --export-felixla --out $OUT_DIR/nothreads" \
  "--tractor-dosage-vcf $ROOT_DIR/testdata/tiny.tractor_dosage.vcf --export-felixla --out $OUT_DIR/nothreads"
do
  # shellcheck disable=SC2086
  if "$BIN_DIR/felixla" $inert --threads 2 >/dev/null 2>"$OUT_DIR/threads.err"; then
    echo "--threads was accepted where it does nothing: $inert" >&2
    exit 1
  fi
  grep -q -- "--threads applies to" "$OUT_DIR/threads.err"
done

# The compatibility subcommands are no longer advertised in --help, but the
# command lines people already have must keep working.
"$BIN_DIR/felixla" --help >"$OUT_DIR/help.txt"
if grep -qi "compatibility" "$OUT_DIR/help.txt"; then
  echo "--help still advertises the compatibility subcommands" >&2
  exit 1
fi
"$BIN_DIR/felixla" from-flare \
  "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  "$ROOT_DIR/testdata/tiny.flare.vcf" \
  2 auto "$OUT_DIR/compat" >/dev/null
"$BIN_DIR/felixla" to-vcf "$OUT_DIR/compat" "$OUT_DIR/compat.vcf.gz" >/dev/null
"$BIN_DIR/felixla" extract "$OUT_DIR/compat" chr1:1-1000 "$OUT_DIR/compat_region" >/dev/null
printf '%s\n' "$OUT_DIR/compat" >"$OUT_DIR/compat.list"
"$BIN_DIR/felixla" concat "$OUT_DIR/compat.list" "$OUT_DIR/compat_merged" >/dev/null


# Selecting alleles out of a packed prefix must land exactly where selecting
# them while packing lands, since the two are different code paths for the
# same filter.
cat >"$OUT_DIR/keep.pvar" <<'SITES'
#CHROM	POS	ID	REF	ALT
chr1	100	.	A	G
chr1	160	.	A	T
chr2	100	.	A	T
SITES

"$BIN_DIR/felixla" --phase-vcf "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
  --flare-vcf "$ROOT_DIR/testdata/tiny.flare.vcf" \
  --extract "$OUT_DIR/keep.pvar" \
  --export-felixla --out "$OUT_DIR/packed_extract" >/dev/null
"$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" \
  --extract "$OUT_DIR/keep.pvar" \
  --export-felixla --out "$OUT_DIR/prefix_extract" >/dev/null

for component in common.geno.bin common.variant.mks common.variant.idx \
                 rare.carrier.bin rare.variant.mks rare.variant.idx \
                 ancblock.bin ancblock.mks ancblock.idx samples; do
  if ! cmp -s "$OUT_DIR/packed_extract.$component" "$OUT_DIR/prefix_extract.$component"; then
    echo "--extract on a prefix differs from --extract while packing: $component" >&2
    exit 1
  fi
done

# --exclude composes the same way. Every block here still covers a retained
# variant, so all of them survive; the case where one does not is below.
"$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" \
  --exclude "$OUT_DIR/keep.pvar" \
  --export-felixla --out "$OUT_DIR/prefix_exclude" >/dev/null
cmp -s "$OUT_DIR/tiny.ancblock.bin" "$OUT_DIR/prefix_exclude.ancblock.bin"
cmp -s "$OUT_DIR/tiny.ancblock.mks" "$OUT_DIR/prefix_exclude.ancblock.mks"

"$BIN_DIR/felixla" --felixla "$OUT_DIR/prefix_exclude" --export-vcf \
  --out "$OUT_DIR/prefix_exclude" >/dev/null
python3 - "$OUT_DIR/prefix_exclude.vcf.gz" <<'EXCLUDED'
import gzip
import sys

kept = []
with gzip.open(sys.argv[1], "rt") as handle:
    for line in handle:
        if not line.startswith("#"):
            fields = line.split("\t")
            kept.append((fields[0], int(fields[1]), fields[3], fields[4]))

# The three excluded alleles are gone; the other six survive, including the
# second ALT of the multiallelic record at chr1:100.
assert ("chr1", 100, "A", "G") not in kept, kept
assert ("chr1", 100, "A", "C") in kept, kept
assert ("chr1", 160, "A", "T") not in kept, kept
assert ("chr2", 100, "A", "T") not in kept, kept
assert len(kept) == 6, kept
EXCLUDED

# A list naming a REF the prefix does not have is a mistake when selecting and
# inert when dropping, exactly as it is while packing.
cat >"$OUT_DIR/badref.pvar" <<'BADREF'
#CHROM	POS	ID	REF	ALT
chr1	160	.	G	T
BADREF
if "$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" --extract "$OUT_DIR/badref.pvar" \
     --export-felixla --out "$OUT_DIR/badref" >/dev/null 2>"$OUT_DIR/badref.err"; then
  echo "--extract accepted a REF that disagrees with the prefix" >&2
  exit 1
fi
grep -q "REF mismatch" "$OUT_DIR/badref.err"
"$BIN_DIR/felixla" --felixla "$OUT_DIR/tiny" --exclude "$OUT_DIR/badref.pvar" \
  --export-felixla --out "$OUT_DIR/badref_drop" >/dev/null
cmp -s "$OUT_DIR/tiny.common.variant.mks" "$OUT_DIR/badref_drop.common.variant.mks"


# The two conventions the README pins down: a region closed at both ends, and
# a variant placed in it by POS alone. The deletion at 195 spans 195-205, so a
# rule that consulted REF length would put it in more than one region.
cat >"$OUT_DIR/boundary.vcf" <<'BOUNDARY_GENO'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	rs1	A	T	.	PASS	.	GT	0|1	1|0
chr1	195	rs2	ACCCCCCCCCC	A	.	PASS	.	GT	0|1	1|0
chr1	200	rs3	G	C	.	PASS	.	GT	1|0	0|1
chr1	205	rs4	T	G	.	PASS	.	GT	0|1	1|0
BOUNDARY_GENO

cat >"$OUT_DIR/boundary.flare.vcf" <<'BOUNDARY_FLARE'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First">
##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second">
##ANCESTRY=<ID=0,Name=AFR>
##ANCESTRY=<ID=1,Name=EUR>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr1	195	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr1	200	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr1	205	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
BOUNDARY_FLARE

"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/boundary.vcf" \
  --flare-vcf "$OUT_DIR/boundary.flare.vcf" \
  --export-felixla --out "$OUT_DIR/boundary" >/dev/null

positions_in_region() {
  "$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --region "$1" \
    --export-felixla --out "$OUT_DIR/region_probe" >/dev/null 2>&1 || { echo ""; return; }
  "$BIN_DIR/felixla" --felixla "$OUT_DIR/region_probe" --export-vcf \
    --out "$OUT_DIR/region_probe" >/dev/null 2>&1
  python3 - "$OUT_DIR/region_probe.vcf.gz" <<'REGION_PROBE'
import gzip
import sys

with gzip.open(sys.argv[1], "rt") as handle:
    print(",".join(line.split("\t")[1] for line in handle if not line.startswith("#")))
REGION_PROBE
}

# Half-open [START, END): START is in, END is not. The variant at 200 is
# outside chr1:100-200 and inside chr1:100-201, which is the whole point.
test "$(positions_in_region chr1:100-201)" = "100,195,200"
test "$(positions_in_region chr1:100-200)" = "100,195"
test "$(positions_in_region chr1:101-200)" = "195"
test "$(positions_in_region chr1:195-196)" = "195"

# Consecutive regions tile: END of one is START of the next, and the variant
# at 195 is claimed by exactly one of them.
test "$(positions_in_region chr1:1-195)" = "100"
test "$(positions_in_region chr1:195-300)" = "195,200,205"

# The deletion spans 195-205 but sits at 195, so a region starting after 195
# does not claim it however far it reaches.
test "$(positions_in_region chr1:196-300)" = "200,205"

# END equal to START covers nothing, which is a mistake rather than an empty
# run that looks like it worked.
if "$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --region chr1:195-195 \
     --export-felixla --out "$OUT_DIR/empty_region" >/dev/null 2>"$OUT_DIR/empty_region.err"; then
  echo "an empty half-open region was accepted" >&2
  exit 1
fi
grep -q "covers nothing" "$OUT_DIR/empty_region.err"

# --extract identifies a variant by CHROM/POS/REF/ALT. The ID column takes no
# part, an allele-swapped entry is refused rather than silently missing, and a
# contig spelled differently from the input matches nothing.
cat >"$OUT_DIR/wrong_id.pvar" <<'WRONG_ID'
#CHROM	POS	ID	REF	ALT
chr1	100	rs_not_this_variant	A	T
WRONG_ID
"$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/wrong_id.pvar" \
  --export-felixla --out "$OUT_DIR/by_coords" >/dev/null
"$BIN_DIR/felixla" --felixla "$OUT_DIR/by_coords" --export-vcf \
  --out "$OUT_DIR/by_coords" >/dev/null 2>&1
python3 - "$OUT_DIR/by_coords.vcf.gz" <<'BY_COORDS'
import gzip
import sys

with gzip.open(sys.argv[1], "rt") as handle:
    rows = [line.split("\t") for line in handle if not line.startswith("#")]
assert len(rows) == 1, rows
assert rows[0][1] == "100", rows
# The ID written out is the source ID with the allele appended, so the rows a
# multiallelic record splits into stay distinguishable. It is a label, not the
# identity that was matched on -- the list said rs_not_this_variant.
assert rows[0][2] == "rs1_A_T", rows
BY_COORDS

cat >"$OUT_DIR/swapped.pvar" <<'SWAPPED_ALLELES'
#CHROM	POS	ID	REF	ALT
chr1	100	.	T	A
SWAPPED_ALLELES
if "$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/swapped.pvar" \
     --export-felixla --out "$OUT_DIR/swapped" >/dev/null 2>"$OUT_DIR/swapped.err"; then
  echo "--extract accepted REF and ALT the wrong way round" >&2
  exit 1
fi
grep -q "REF mismatch" "$OUT_DIR/swapped.err"

# A contig may be spelled with or without the chr prefix, and a list that
# mixes the two at one coordinate still contributes both alleles rather than
# losing whichever sorted second.
cat >"$OUT_DIR/nochr.pvar" <<'NO_CHR_PREFIX'
#CHROM	POS	ID	REF	ALT
1	100	.	A	T
NO_CHR_PREFIX
"$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/nochr.pvar" \
  --export-felixla --out "$OUT_DIR/nochr" >/dev/null
for component in common.geno.bin common.variant.mks rare.carrier.bin rare.variant.mks; do
  if ! cmp -s "$OUT_DIR/by_coords.$component" "$OUT_DIR/nochr.$component"; then
    echo "a list spelled without the chr prefix selected something else: $component" >&2
    exit 1
  fi
done

cat >"$OUT_DIR/mixed_chr.pvar" <<'MIXED_CHR'
#CHROM	POS	ID	REF	ALT
chr1	100	.	A	T
1	200	.	G	C
MIXED_CHR
"$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/mixed_chr.pvar" \
  --export-felixla --out "$OUT_DIR/mixed_chr" >/dev/null
"$BIN_DIR/felixla" --felixla "$OUT_DIR/mixed_chr" --export-vcf \
  --out "$OUT_DIR/mixed_chr" >/dev/null 2>&1
python3 - "$OUT_DIR/mixed_chr.vcf.gz" <<'MIXED_CHR_CHECK'
import gzip
import sys

with gzip.open(sys.argv[1], "rt") as handle:
    rows = [line.split("\t") for line in handle if not line.startswith("#")]
assert [(row[0], int(row[1])) for row in rows] == [("chr1", 100), ("chr1", 200)], rows
MIXED_CHR_CHECK

# Packing with the same list must land in the same place, whichever spelling
# it used, or the two routes disagree about what the list meant.
"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/boundary.vcf" \
  --flare-vcf "$OUT_DIR/boundary.flare.vcf" \
  --extract "$OUT_DIR/mixed_chr.pvar" \
  --export-felixla --out "$OUT_DIR/mixed_chr_packed" >/dev/null
for component in common.geno.bin common.variant.mks common.variant.idx \
                 rare.carrier.bin rare.variant.mks rare.variant.idx \
                 ancblock.bin ancblock.mks ancblock.idx samples; do
  if ! cmp -s "$OUT_DIR/mixed_chr_packed.$component" "$OUT_DIR/mixed_chr.$component"; then
    echo "mixed chr spellings differ between packing and prefix extract: $component" >&2
    exit 1
  fi
done

# The boundary fixture has a single block spanning every variant, so keeping
# one variant keeps that block whole: a block that still covers a retained
# variant is not clipped to it.
cat >"$OUT_DIR/first_only.pvar" <<'FIRST_ONLY'
#CHROM	POS	ID	REF	ALT
chr1	100	.	A	T
FIRST_ONLY
"$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/first_only.pvar" \
  --export-felixla --out "$OUT_DIR/first_only" >/dev/null
cmp -s "$OUT_DIR/boundary.ancblock.bin" "$OUT_DIR/first_only.ancblock.bin"
cmp -s "$OUT_DIR/boundary.ancblock.mks" "$OUT_DIR/first_only.ancblock.mks"
cmp -s "$OUT_DIR/boundary.ancblock.idx" "$OUT_DIR/first_only.ancblock.idx"


# A one-column line is a gnomAD-style variant ID. Either separator, either
# contig spelling, and the result has to be the bytes the equivalent PVAR
# produces -- a second way to name a variant is only worth having if it names
# the same one.
cat >"$OUT_DIR/ids_colon.txt" <<'IDS_COLON'
chr1:100:A:T
chr1:200:G:C
IDS_COLON
cat >"$OUT_DIR/ids_dash.txt" <<'IDS_DASH'
1-100-A-T
1-200-G-C
IDS_DASH
cat >"$OUT_DIR/ids_columns.pvar" <<'IDS_COLUMNS'
#CHROM	POS	ID	REF	ALT
chr1	100	.	A	T
chr1	200	.	G	C
IDS_COLUMNS

for spelling in ids_colon.txt ids_dash.txt ids_columns.pvar; do
  "$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/$spelling" \
    --export-felixla --out "$OUT_DIR/ids_${spelling%%.*}" >/dev/null
done
for variant in ids_colon ids_dash; do
  for component in common.geno.bin common.variant.mks common.variant.idx \
                   rare.carrier.bin rare.variant.mks rare.variant.idx \
                   ancblock.bin ancblock.mks ancblock.idx samples; do
    if ! cmp -s "$OUT_DIR/ids_ids_columns.$component" "$OUT_DIR/ids_$variant.$component"; then
      echo "variant IDs in $variant selected something other than the PVAR: $component" >&2
      exit 1
    fi
  done
done

# The same list has to mean the same thing while packing.
"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/boundary.vcf" \
  --flare-vcf "$OUT_DIR/boundary.flare.vcf" \
  --extract "$OUT_DIR/ids_dash.txt" \
  --export-felixla --out "$OUT_DIR/ids_packed" >/dev/null
for component in common.geno.bin common.variant.mks common.variant.idx \
                 rare.carrier.bin rare.variant.mks rare.variant.idx \
                 ancblock.bin ancblock.mks ancblock.idx samples; do
  if ! cmp -s "$OUT_DIR/ids_packed.$component" "$OUT_DIR/ids_ids_dash.$component"; then
    echo "a variant ID list differs between packing and prefix extract: $component" >&2
    exit 1
  fi
done

# An indel keeps its full alleles in an ID, and the separators inside a contig
# name are not field separators: the split counts three from the right.
printf 'chr1:195:ACCCCCCCCCC:A\n' >"$OUT_DIR/ids_indel.txt"
"$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/ids_indel.txt" \
  --export-felixla --out "$OUT_DIR/ids_indel" >/dev/null
"$BIN_DIR/felixla" --felixla "$OUT_DIR/ids_indel" --export-vcf \
  --out "$OUT_DIR/ids_indel" >/dev/null 2>&1
python3 - "$OUT_DIR/ids_indel.vcf.gz" <<'IDS_INDEL_CHECK'
import gzip
import sys

with gzip.open(sys.argv[1], "rt") as handle:
    rows = [line.split("\t") for line in handle if not line.startswith("#")]
assert [(row[1], row[3], row[4]) for row in rows] == [("195", "ACCCCCCCCCC", "A")], rows
IDS_INDEL_CHECK

# Anything that is not four fields is refused rather than guessed at, which is
# what keeps a file of bare rs IDs from looking like it worked.
for malformed in rs12345 chr1:100 chr1:100:A chr1:100::T; do
  printf '%s\n' "$malformed" >"$OUT_DIR/ids_bad.txt"
  if "$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/ids_bad.txt" \
       --export-felixla --out "$OUT_DIR/ids_bad" >/dev/null 2>&1; then
    echo "a malformed variant ID was accepted: $malformed" >&2
    exit 1
  fi
done


# The separators inside a contig name are not field separators. FELIXla packs
# autosomes only, so such a name never reaches the packed data, but it can
# still appear in a site list -- and the right-to-left split is what keeps the
# list readable instead of fatal. Splitting from the left would read
# HLA-A-100-A-T as position "A-T" and die on the whole file, taking the chr1
# entry with it.
cat >"$OUT_DIR/ids_contig_sep.txt" <<'IDS_CONTIG_SEP'
chr1-100-A-T
HLA-A-100-A-T
IDS_CONTIG_SEP

"$BIN_DIR/felixla" --felixla "$OUT_DIR/boundary" --extract "$OUT_DIR/ids_contig_sep.txt" \
  --export-felixla --out "$OUT_DIR/ids_contig_sep" >/dev/null
"$BIN_DIR/felixla" --felixla "$OUT_DIR/ids_contig_sep" --export-vcf \
  --out "$OUT_DIR/ids_contig_sep" >/dev/null 2>&1
python3 - "$OUT_DIR/ids_contig_sep.vcf.gz" <<'IDS_CONTIG_SEP_CHECK'
import gzip
import sys

with gzip.open(sys.argv[1], "rt") as handle:
    rows = [line.split("\t") for line in handle if not line.startswith("#")]
# The chr1 entry is selected and the HLA-A one names a contig the prefix does
# not have, so it is simply inert.
assert [(row[0], row[1]) for row in rows] == [("chr1", 100)][:1] or \
    [(row[0], int(row[1])) for row in rows] == [("chr1", 100)], rows
IDS_CONTIG_SEP_CHECK

# Packing with the same list must read it the same way.
"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/boundary.vcf" \
  --flare-vcf "$OUT_DIR/boundary.flare.vcf" \
  --extract "$OUT_DIR/ids_contig_sep.txt" \
  --export-felixla --out "$OUT_DIR/ids_contig_sep_packed" >/dev/null
for component in common.geno.bin common.variant.mks common.variant.idx \
                 rare.carrier.bin rare.variant.mks rare.variant.idx \
                 ancblock.bin ancblock.mks ancblock.idx samples; do
  if ! cmp -s "$OUT_DIR/ids_contig_sep_packed.$component" "$OUT_DIR/ids_contig_sep.$component"; then
    echo "a contig name containing a separator differs between the two routes: $component" >&2
    exit 1
  fi
done

# Only the autosomes are supported, and a contig outside 1-22 is refused by
# name rather than packed into a prefix whose ancestry blocks mean nothing.
for refused in chrX chrM chr23 chrUn_GL000220v1; do
  sed -e "s/^chr1\t/$refused\t/" -e "s/ID=chr1/ID=$refused/" \
    "$OUT_DIR/boundary.vcf" >"$OUT_DIR/refused.vcf"
  sed -e "s/^chr1\t/$refused\t/" -e "s/ID=chr1/ID=$refused/" \
    "$OUT_DIR/boundary.flare.vcf" >"$OUT_DIR/refused.flare.vcf"
  if "$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/refused.vcf" \
       --flare-vcf "$OUT_DIR/refused.flare.vcf" \
       --export-felixla --out "$OUT_DIR/refused" \
       >/dev/null 2>"$OUT_DIR/refused.err"; then
    echo "contig $refused was packed, but only autosomes are supported" >&2
    exit 1
  fi
  grep -q "is not supported" "$OUT_DIR/refused.err"
done

# chr22 is an autosome and packs.
sed -e "s/^chr1\t/chr22\t/" -e "s/ID=chr1/ID=chr22/" "$OUT_DIR/boundary.vcf" >"$OUT_DIR/chr22.vcf"
sed -e "s/^chr1\t/chr22\t/" -e "s/ID=chr1/ID=chr22/" "$OUT_DIR/boundary.flare.vcf" >"$OUT_DIR/chr22.flare.vcf"
"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/chr22.vcf" --flare-vcf "$OUT_DIR/chr22.flare.vcf" \
  --export-felixla --out "$OUT_DIR/chr22" >/dev/null

# --n-ancestries is gone; the count comes from the input header.
if "$BIN_DIR/felixla" --phase-vcf "$ROOT_DIR/testdata/tiny.genotypes.vcf" \
     --flare-vcf "$ROOT_DIR/testdata/tiny.flare.vcf" --n-ancestries 2 \
     --export-felixla --out "$OUT_DIR/nanc" >/dev/null 2>"$OUT_DIR/nanc.err"; then
  echo "--n-ancestries was accepted after removal" >&2
  exit 1
fi
grep -q "has been removed" "$OUT_DIR/nanc.err"


# A block that no retained variant falls in is dropped, and the blocks left
# behind may have gaps between them. The fixture gives every variant its own
# block so the two ends can be kept and the middle dropped.
cat >"$OUT_DIR/blocks.vcf" <<'BLOCKS_GENO'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	T	.	PASS	.	GT	0|1	1|0
chr1	200	.	A	T	.	PASS	.	GT	0|1	1|0
chr1	300	.	A	T	.	PASS	.	GT	0|1	1|0
chr1	400	.	A	T	.	PASS	.	GT	0|1	1|0
chr1	500	.	A	T	.	PASS	.	GT	0|1	1|0
BLOCKS_GENO

cat >"$OUT_DIR/blocks.flare.vcf" <<'BLOCKS_FLARE'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First">
##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second">
##ANCESTRY=<ID=0,Name=AFR>
##ANCESTRY=<ID=1,Name=EUR>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr1	200	.	A	C	.	PASS	.	AN1:AN2	1:0	0:1
chr1	300	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr1	400	.	A	C	.	PASS	.	AN1:AN2	1:0	0:1
chr1	500	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
BLOCKS_FLARE

"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/blocks.vcf" \
  --flare-vcf "$OUT_DIR/blocks.flare.vcf" \
  --export-felixla --out "$OUT_DIR/blocks" >/dev/null

cat >"$OUT_DIR/blocks_ends.pvar" <<'BLOCKS_ENDS'
#CHROM	POS	ID	REF	ALT
chr1	100	.	A	T
chr1	500	.	A	T
BLOCKS_ENDS

lai_rows() {
  "$BIN_DIR/felixla" --felixla "$1" --export-lai --out "$1" >/dev/null 2>&1
  python3 - "$1.lai.gz" <<'LAI_ROWS'
import gzip
import sys

with gzip.open(sys.argv[1], "rt") as handle:
    print(" ".join(
        "%s-%s" % (line.split("\t")[1], line.split("\t")[2])
        for line in handle if not line.startswith("#")
    ))
LAI_ROWS
}

test "$(lai_rows "$OUT_DIR/blocks")" = "1-100 101-200 201-300 301-400 401-500"

"$BIN_DIR/felixla" --felixla "$OUT_DIR/blocks" --extract "$OUT_DIR/blocks_ends.pvar" \
  --export-felixla --out "$OUT_DIR/blocks_ends" >/dev/null
# The three middle blocks held only variants the list did not name, so they are
# gone, and what is left is not contiguous.
test "$(lai_rows "$OUT_DIR/blocks_ends")" = "1-100 401-500"
grep -q $'ancestry_blocks\t2' "$OUT_DIR/blocks_ends.meta"

# Packing with the same list has to reach the same place, blocks included.
"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/blocks.vcf" \
  --flare-vcf "$OUT_DIR/blocks.flare.vcf" \
  --extract "$OUT_DIR/blocks_ends.pvar" \
  --export-felixla --out "$OUT_DIR/blocks_ends_packed" >/dev/null
for component in common.geno.bin common.variant.mks common.variant.idx \
                 rare.carrier.bin rare.variant.mks rare.variant.idx \
                 ancblock.bin ancblock.mks ancblock.idx samples; do
  if ! cmp -s "$OUT_DIR/blocks_ends_packed.$component" "$OUT_DIR/blocks_ends.$component"; then
    echo "dropping empty blocks differs between packing and prefix extract: $component" >&2
    exit 1
  fi
done

# A prefix whose blocks have gaps is still a valid prefix: concat revalidates
# all eleven files, and block ids must have stayed consecutive through the
# renumbering.
printf '%s\n' "$OUT_DIR/blocks_ends" >"$OUT_DIR/blocks_ends.list"
"$BIN_DIR/felixla" --merge-list "$OUT_DIR/blocks_ends.list" \
  --export-felixla --out "$OUT_DIR/blocks_ends_merged" >/dev/null
test "$(lai_rows "$OUT_DIR/blocks_ends_merged")" = "1-100 401-500"

# --exclude reaches the same rule from the other side.
"$BIN_DIR/felixla" --felixla "$OUT_DIR/blocks" --exclude "$OUT_DIR/blocks_ends.pvar" \
  --export-felixla --out "$OUT_DIR/blocks_middle" >/dev/null
test "$(lai_rows "$OUT_DIR/blocks_middle")" = "101-200 201-300 301-400"

# Per-sample proportions are summed over whatever the blocks cover, so a
# contig with a hole in it is worth saying out loud. Where a contig starts and
# stops is not: a region extract legitimately holds part of one.
warned() {
  "$BIN_DIR/felixla" --felixla "$1" --export-global-admixture \
    --out "$1" 2>"$1.warn" >/dev/null
  grep -oE "[0-9]+ hole\(s\) totalling [0-9]+ bp, the first at [^.]*" "$1.warn" \
    || echo "no warning"
}

# A hole in the middle, from blocks that held no retained variant.
test "$(warned "$OUT_DIR/blocks_ends")" \
  = "1 hole(s) totalling 300 bp, the first at chr1:101-400"

# The rows themselves are still internally consistent: each sums to one over
# the span that is covered. That is what makes them easy to misread.
python3 - "$OUT_DIR/blocks_ends.global.admixture.tsv" <<'HOLEY_ROWS'
import sys

with open(sys.argv[1]) as handle:
    rows = [line.rstrip("\n").split("\t") for line in handle][1:]
for row in rows:
    total = sum(float(value) for value in row[1:])
    assert abs(total - 1.0) < 1e-5, row
HOLEY_ROWS

# A contig that simply starts late is still one piece, and says nothing. Both
# ways of getting there: blocks dropped off the front, and a --region that
# asked for the rest.
cat >"$OUT_DIR/blocks_last.pvar" <<'BLOCKS_LAST'
#CHROM	POS	ID	REF	ALT
chr1	500	.	A	T
BLOCKS_LAST
"$BIN_DIR/felixla" --felixla "$OUT_DIR/blocks" --extract "$OUT_DIR/blocks_last.pvar" \
  --export-felixla --out "$OUT_DIR/blocks_last" >/dev/null
test "$(lai_rows "$OUT_DIR/blocks_last")" = "401-500"
test "$(warned "$OUT_DIR/blocks_last")" = "no warning"

"$BIN_DIR/felixla" --felixla "$OUT_DIR/blocks" --region chr1:200-451 \
  --export-felixla --out "$OUT_DIR/blocks_region" >/dev/null
test "$(warned "$OUT_DIR/blocks_region")" = "no warning"

# Blocks running unbroken say nothing, or the warning would be noise that gets
# ignored on the run that matters.
test "$(warned "$OUT_DIR/blocks")" = "no warning"

# A contig boundary is not a break: the next contig is expected to start at
# one, not to continue from where the previous one ended.
cat >"$OUT_DIR/twochr.vcf" <<'TWOCHR_GENO'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##contig=<ID=chr2,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	T	.	PASS	.	GT	0|1	1|0
chr1	300	.	A	T	.	PASS	.	GT	0|1	1|0
chr2	1000	.	A	T	.	PASS	.	GT	0|1	1|0
chr2	1200	.	A	T	.	PASS	.	GT	0|1	1|0
TWOCHR_GENO

cat >"$OUT_DIR/twochr.flare.vcf" <<'TWOCHR_FLARE'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##contig=<ID=chr2,length=100000>
##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First">
##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second">
##ANCESTRY=<ID=0,Name=AFR>
##ANCESTRY=<ID=1,Name=EUR>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr1	300	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr2	1000	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
chr2	1200	.	A	C	.	PASS	.	AN1:AN2	0:1	1:0
TWOCHR_FLARE

"$BIN_DIR/felixla" --phase-vcf "$OUT_DIR/twochr.vcf" \
  --flare-vcf "$OUT_DIR/twochr.flare.vcf" \
  --export-felixla --out "$OUT_DIR/twochr" >/dev/null
test "$(lai_rows "$OUT_DIR/twochr")" = "1-300 1-1200"
test "$(warned "$OUT_DIR/twochr")" = "no warning"

# The same boundary, with the second contig starting well past where the first
# one stopped. Packing always starts a contig at one, so this takes a merge of
# a whole chr1 with a region-extracted chr2 -- which is what chunked conversion
# produces. A check that forgot to compare contigs would report a 699 bp hole
# at chr2:301-999 here.
"$BIN_DIR/felixla" --felixla "$OUT_DIR/twochr" --region chr1:1-400 \
  --export-felixla --out "$OUT_DIR/part_chr1" >/dev/null
"$BIN_DIR/felixla" --felixla "$OUT_DIR/twochr" --region chr2:1000-1300 \
  --export-felixla --out "$OUT_DIR/part_chr2" >/dev/null
printf '%s\n%s\n' "$OUT_DIR/part_chr1" "$OUT_DIR/part_chr2" >"$OUT_DIR/parts.list"
"$BIN_DIR/felixla" --merge-list "$OUT_DIR/parts.list" \
  --export-felixla --out "$OUT_DIR/parts" >/dev/null
test "$(lai_rows "$OUT_DIR/parts")" = "1-300 1000-1200"
test "$(warned "$OUT_DIR/parts")" = "no warning"

echo "tiny smoke test passed"
