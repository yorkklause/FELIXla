# Ancestry-aware Packed Backend Prototype

> **Acknowledgement:** designed by Kai, implemented by codex.

This repo starts from the ancestry-aware packed backend design for SAIGE-TRACTOR.
The first implemented piece is a standalone converter:

```bash
make
bin/flare_subset_to_tractor_hybrid \
  genotype.phased.vcf.gz \
  flare.anc.vcf.gz \
  3 \
  512 \
  out/chr1
```

For All of Us Workbench or other Linux x86_64 environments where you do not
want to compile, use the checked-in static binaries:

```bash
export PATH="$PWD/prebuilt/linux-x86_64-static:$PATH"
estimate_mac_threshold testdata/tiny.genotypes.vcf
estimate_mac_threshold testdata/tiny.genotypes.vcf --sample-every 100
```

These prebuilt tools do not need conda or htslib at runtime. They are intended
for local VCF/BCF paths; remote URL/S3/GCS support is disabled in the bundled
htslib to keep the binaries dependency-free.

A tools-only Docker image with the same standalone binaries is published as:

```bash
docker pull kyuan1024/tractor-hybrid-tools:latest
docker run --rm -v "$PWD:/data" -w /data kyuan1024/tractor-hybrid-tools:latest \
  rfmix_msp_to_tractor_hybrid genotype.phased.vcf.gz rfmix.msp.tsv.gz 5 512 out/chr22
```

The forward converters and MAC threshold estimator print progress to stderr.
When the input has a usable `.csi` or `.tbi` index with record statistics, the
tool reports `scanned records / total records` and a percentage; otherwise it
keeps reporting the number of scanned records. For the phased genotype plus
FLARE converter, the percentage is based on genotype VCF records, while the
converted count is the number of emitted split-biallelic hybrid variants.

Chunking is driven by region strings passed as the optional final argument to
`flare_subset_to_tractor_hybrid`. Use a prepared interval file, such as a
SHAPEIT5 4 cM chunk file, and pass each desired `chr:start-end` interval to the
converter. Region coordinates follow VCF/tabix style: they are 1-based and both
ends are inclusive, `[start, end]`. A variant belongs to a chunk when
`start <= POS <= end`, so adjacent chunks should start at the previous chunk's
`end + 1`:

```text
chr1:1-50000000
chr1:50000001-100000000
```

Do not define adjacent chunks as `chr1:1-50000000` and
`chr1:50000000-100000000`, because a variant with `POS=50000000` would be
included in both. This chunk convention is separate from FLARE LAI intervals,
which remain `(previous_lai_pos, current_lai_pos]`.

For an interval file with one tabix-style region in the first column:

```bash
i=0
while read -r region _; do
  case "$region" in ""|\#*) continue ;; esac
  chrom="${region%%:*}"
  i=$((i + 1))
  chunk="$(printf '%04d' "$i")"
  flare_subset_to_tractor_hybrid \
    "phase/${chrom}.phased.vcf.gz" \
    "flare/${chrom}.flare.vcf.gz" \
    5 \
    512 \
    "hybrid/${chrom}.chunk${chunk}" \
    "$region"
done < shapeit5.4cm.regions.txt
```

If your interval file has separate `chrom start end` columns, convert the
desired columns to region strings first:

```bash
awk 'BEGIN{OFS=""} !/^#/ && NF >= 3 {print $1,":",$2,"-",$3}' intervals.tsv > regions.txt
```

You can also call the converter directly for one region:

```bash
flare_subset_to_tractor_hybrid \
  phase/chr22.phased.vcf.gz \
  flare/chr22.flare.vcf.gz \
  5 \
  512 \
  hybrid/chr22.chunk0001 \
  chr22:1-50000000
```

The packed files can be converted back to a split-biallelic VCF:

```bash
bin/tractor_hybrid_to_vcf out/chr1 out/chr1.roundtrip.vcf.gz
```

The reverse converter writes BGZF-compressed VCF and creates a tabix index at
`out/chr1.roundtrip.vcf.gz.tbi`.

Extract a 1-based inclusive region from an existing packed hybrid prefix into a
new packed hybrid prefix:

```bash
bin/tractor_hybrid_extract_region \
  out/chr22 \
  chr22:16000000-17000000 \
  out/chr22.16_17mb
```

The extractor rewrites marker indexes and payload offsets, renumbers selected
split variants from zero, and keeps only ancestry blocks overlapping the
requested region. It also accepts split arguments:
`bin/tractor_hybrid_extract_region in_prefix chr22 16000000 17000000 out_prefix`.

Convert a phased genotype VCF/BCF plus an RFMix `.msp.tsv` local ancestry file:

```bash
bin/rfmix_msp_to_tractor_hybrid \
  genotype.phased.vcf.gz \
  rfmix.msp.tsv.gz \
  5 \
  512 \
  out/chr22
```

The MSP converter expects haplotype columns in VCF sample order, named like
`sample.0` and `sample.1`, and validates the optional
`#Subpopulation order/codes:` line against `n_ancestries`. It accepts exact
chromosome names or simple `chr`/non-`chr` equivalents, for example MSP `22`
with VCF `chr22`. Adjacent MSP rows often share an endpoint; the converter
treats the first interval as including `spos`, then assigns later shared
`epos`/`spos` boundaries to the previous interval. Progress reports include
scanned VCF records, MSP rows consumed, converted split variants, and
rare/common counts.

Convert a SAIGE-TRACTOR-style VCF/BCF that already has hardcall
`DS1..DSk` and `ANC1..ANCk` FORMAT fields:

```bash
bin/tractor_dosage_vcf_to_hybrid \
  tractor.step2.input.vcf.gz \
  3 \
  512 \
  out/chr1
```

This path reconstructs canonical haplotypes from the per-sample ancestry counts
and ancestry-specific ALT dosages. It is equivalent for SAIGE-TRACTOR step2
fields (`ANC#`, `DS#`, `DSALL`) but it is not a recovery of the original
phased haplotype order. The converter requires hardcall integer values: each
sample's `ANC#` values must sum to `2`, each `DS#` value must be an integer in
`0..ANC#`, and fractional imputed dosages are rejected.

Estimate a MAC threshold from a phased genotype VCF/BCF:

```bash
bin/estimate_mac_threshold genotype.phased.vcf.gz
```

The estimator's progress line also reports split ALT variants seen so far, the
maximum observed MAC, elapsed time, ETA when the input record count is known,
and scan rate. `--sample-every N` keeps the full scan/progress denominator but
uses only every Nth VCF record for the MAC distribution; the default is `1`,
which uses all records. `--max-records` caps the denominator when an indexed
input has more records than the scan limit.

The estimator scans split ALT MAC values and compares sparse carrier payload
against dense bitset payload. By default it optimizes storage bytes only. You can
add a query-work term when you want to penalize sparse carrier iteration or dense
word scanning:

```bash
bin/estimate_mac_threshold genotype.phased.vcf.gz \
  --query-weight 1 \
  --dense-word-cost 1 \
  --sparse-carrier-cost 16
```

Compare two phased genotype VCF/BCF files:

```bash
bin/compare_vcfs expected.vcf.gz observed.vcf.gz --split-multiallelic
```

`--split-multiallelic` compares both inputs after logically splitting each ALT
allele, which is useful for checking a packed roundtrip against an original
multi-allelic VCF.

SAIGE/SAIGE-TRACTOR step2 adapter:

- `saige_step2_adapter/TractorHybrid.hpp`
- `saige_step2_adapter/TractorHybrid.cpp`
- `saige_step2_adapter/Main.cpp.patch`
- `saige_step2_adapter/Geno.R.patch`
- `saige_step2_adapter/SAIGE_SPATest_Tractor.R.patch`
- `saige_step2_adapter/SAIGE_SPATest_Marker.R.patch`
- `docker/saigetractor-hybrid/Dockerfile`
- `docker/saigetractor-hybrid/patch_saigetractor_hybrid.py`

Build the local SAIGE-TRACTOR image with tractor_hybrid step2 support:

```bash
docker build --platform linux/amd64 \
  -t saigetractor:1.4.9-tractor-hybrid \
  -f docker/saigetractor-hybrid/Dockerfile .
```

The Docker image also includes the standalone hybrid utilities in `/scripts/bin`
and their source code in `/scripts/src`:

```bash
/scripts/bin/flare_subset_to_tractor_hybrid
/scripts/bin/tractor_hybrid_to_vcf
/scripts/bin/tractor_hybrid_extract_region
/scripts/bin/tractor_dosage_vcf_to_hybrid
/scripts/bin/rfmix_msp_to_tractor_hybrid
/scripts/bin/estimate_mac_threshold
/scripts/bin/compare_vcfs
```

Frozen Docker Hub image for the current stable `.2` line:

```text
kyuan1024/saigetractor:1.4.9-tractor-hybrid.2
sha256:9f872fbd2df8e2ee57a427fbbb28b529fb4165ef823e4894750b17f60910fb8e
```

Run step2 with the packed prefix:

```bash
docker run --rm --platform linux/amd64 \
  -v /path/to/data:/data \
  saigetractor:1.4.9-tractor-hybrid \
  step2_SPAtests.R \
    --tractorHybridPrefix=/data/chr1 \
    --chrom=chr1 \
    --is_admixed=TRUE \
    --number_of_ancestry=3 \
    --markers_per_chunk=1000 \
    --GMMATmodelFile=/data/null.rda \
    --varianceRatioFile=/data/varianceRatio.txt \
    --SAIGEOutputFile=/data/chr1.tractor.out
```

Use `--markers_per_chunk=1000` to process about 1000 variants per step2 chunk.
SAIGE-TRACTOR 1.4.9 enforces `1000` as the minimum for single-variant tests; if
you omit this option the default is typically larger, such as `10000`.

The adapter is designed to be copied into the SAIGE/SAIGE-TRACTOR `src/`
directory. It has been checked against `wzhou88/saigetractor:1.4.9`
(`sha256:c599ccff1f3c46322809f0239b56fd4fc5caa44671cd833c65b7f372c0481799`).
It adds a `tractor_hybrid` genotype backend for step2 only. The hook is narrow:
`Unified_getOneMarker()` gets one new branch for ordinary dosage, and
`Unified_getOneMarker_Admixed()` gets one new branch for SAIGE-TRACTOR fields.
The score-test and TRACTOR internal calculation code can stay unchanged.

`TractorHybridClass::getOneMarker()` returns the ordinary ALT dosage vector
`0/1/2` in SAIGE sample-in-model order. Missing rate is always `0` because the
forward converter rejects missing genotypes and missing LAI calls.

For SAIGE-TRACTOR code paths that need ancestry-specific ALT dosage,
`TractorHybridClass::getOneMarkerAncestry()` also fills `GByAncestry`, an
`n_samples_in_model x n_ancestries` matrix. Rare variants use the ancestry bits
stored directly in each carrier record; common variants combine the dense ALT
haplotype bitset with the cached ancestry block.

For SAIGE-TRACTOR 1.4.9 admixed step2, the main entry point is
`TractorHybridClass::getOneMarkerAdmixedField()`. It serves the same FORMAT-like
fields that the existing VCF path expects:

- `ANC1..ANCk`: number of haplotypes from each ancestry per sample, independent of ALT genotype.
- `DS1..DSk`: ALT dosage carried on haplotypes from each ancestry.
- `DSALL`, `DS`, or `GT`: total ALT dosage `0/1/2`.

The reader supports two access modes. Random-like increasing
`global_variant_index` access is available for small indexed runs. For AoU-scale
step2 traversal, use the streaming iterator methods and pass dummy `0` indices
from R, VCF-style; this avoids loading the full `.mks` marker stream into R
memory. In admixed streaming mode the reader is intentionally lazy: it does not
advance after every field read, because SAIGE-TRACTOR queries several fields for
the same marker before moving to the next marker. Common marker bitsets and
ancestry blocks are buffered/cached for these repeated same-marker reads.

Current step2 reader optimizations:

- Payload files use explicit libc read buffers: `64 MiB` for `.bin` payloads and `8 MiB` for `.mks`/`.idx` sidecars.
- On Linux, payload and marker files are opened with `posix_fadvise(..., POSIX_FADV_SEQUENTIAL)` so the kernel can optimize readahead for chromosome traversal.
- The reader tracks current file offsets and skips redundant `fseeko()` calls when traversal is already sequential.
- Rare carriers and common bitsets are loaded in bulk reads, then decoded from memory.
- Same-marker results are cached across SAIGE-TRACTOR's repeated `ANC#`, `DS#`, and `DSALL` requests; ancestry-count matrices are also cached per ancestry block.

For `tractor_hybrid` input, the patched step2 scripts print a chunk-level timing
line:

```text
tractor_hybrid timing chunk 1: total=12.345s input_io=0.123s reader_decode=0.456s hybrid_get_marker=0.789s get_marker_overhead=0.210s marker_pvalue=9.876s reset_zero=0.012s condition_total=0.000s condition_cache_hit=0 condition_cache_miss=0 impute_qc=0.111s variance_ratio=0.222s joint_cct_spa=0.333s output_write=0.321s other_cpp_or_r=0.681s
```

`input_io` is time spent in tracked file seek/read calls inside the hybrid
reader. `reader_decode` is time spent constructing dosage and ancestry-count
vectors/matrices from packed data. `hybrid_get_marker` is the full C++ marker
fetch call, so it includes `input_io`, `reader_decode`, vector materialization,
and cache-copy overhead. `get_marker_overhead` is the part of marker fetch not
counted by file IO or packed decode. `marker_pvalue` wraps SAIGE's marker p-value
calculation, including SPA/Firth when triggered. `reset_zero` tracks repeated
vector/matrix zeroing before each marker, `impute_qc` tracks the impute/flip and
post-QC helper call, `variance_ratio` tracks variance-ratio assignment, and
`joint_cct_spa` tracks admixed joint/CCT/SPA-ER combination work outside the
per-ancestry p-value calls. `condition_total` wraps the conditioning-haplotype
setup and is diagnostic because it can contain nested marker reads or p-value
work when conditioning is active. For `tractor_hybrid`, conditioning setup is
cached across consecutive variants in the same simplified ancestry block;
`condition_cache_hit` and `condition_cache_miss` report that block-level reuse.
`output_write` wraps the admixed single-marker result writer. `other_cpp_or_r`
is the remaining chunk elapsed time outside the non-overlapping measured
sections.

Frozen version `.2` keeps the same statistical settings and adds a second
diagnostic line that splits `SAIGEClass::getMarkerPval()` internally:

```text
tractor_hybrid pvalue timing chunk 1: pvalue_score=1.234s score_fast_calls=100 score_slow_calls=0 scorefast_extract=0.111s scorefast_projection=0.222s scorefast_variance=0.333s scorefast_result=0.444s scorefast_gtilde=0.555s scorefast_gtilde_calls=100 pvalue_alloc=0.123s pvalue_getadjg=0.456s getadjg_calls=10 getadjg_accum=0.234s getadjg_projection=0.222s getadjg_spa=0.100s getadjg_spa_calls=3 getadjg_firth=0.050s getadjg_firth_calls=1 getadjg_condition=0.020s getadjg_condition_calls=1 getadjg_region=0.286s getadjg_region_calls=5 getadjg_other=0.000s getadjg_other_calls=0 pvalue_spa=0.789s spa_calls=10 pvalue_firth=2.345s firth_calls=3 pvalue_er=0.000s pvalue_condition_adjust=0.111s pvalue_region_finalize=0.222s
```

These fields are instrumentation only. They do not change score, SPA, Firth,
variance-ratio, filtering, or output calculations.
The current `.2` image also uses an optimized `scoreTestFast()` sparse path that
computes the same algebra without materializing the large `X1`, `A1`, `g1`, and
`res1` temporary copies for every call, plus a `getadjGFast()` projection path
that updates the adjusted genotype vector by column without constructing an
extra dense matrix-vector product temporary. In the region path, `.2` also lets
`scoreTestFast()` emit the full adjusted genotype vector (`scorefast_gtilde`) so
the later SPA/Firth/region logic can reuse it instead of calling
`getadjGFast()` again for the same marker.

Build dependency:

- htslib headers and library are required. On macOS with Homebrew, install with `brew install htslib`.

Static-style builds:

```bash
make static
make test-static
```

`make static` writes binaries to `bin-static/` and links `libhts.a` directly
when a static htslib archive is available. On Linux, you can request a fully
static executable with:

```bash
make static STATIC_FULLY=1
```

That requires static transitive dependencies for htslib, such as zlib,
libdeflate, bzip2, and xz/lzma. On macOS, the system runtime remains dynamically
linked by platform design, but the htslib dependency is still linked from
`libhts.a` rather than `libhts.dylib`.
- If htslib is installed in a non-standard location, build with:

```bash
make HTSLIB_CFLAGS="-I/path/to/htslib/include" HTSLIB_LIBS="-L/path/to/htslib/lib -lhts"
```

Arguments:

1. `genotype.phased.vcf.gz`: phased genotype VCF/BCF.
2. `flare.anc.vcf.gz`: FLARE local ancestry VCF/BCF with `FORMAT/AN1` and `FORMAT/AN2`.
3. `n_ancestries`: number of ancestry labels, max `32`.
4. `rare_threshold`: variants with `MAC <= rare_threshold` are stored sparse.
5. `out_prefix`: prefix for output files.

Outputs:

- `<prefix>.common.geno.bin`
- `<prefix>.common.variant.mks` (binary marker stream)
- `<prefix>.common.variant.idx` (binary offset index for the marker stream)
- `<prefix>.rare.carrier.bin`
- `<prefix>.rare.variant.mks` (binary marker stream)
- `<prefix>.rare.variant.idx` (binary offset index for the marker stream)
- `<prefix>.ancblock.bin`
- `<prefix>.ancblock.mks` (binary marker stream)
- `<prefix>.ancblock.idx` (binary offset index for the marker stream)
- `<prefix>.samples`
- `<prefix>.meta`

Input assumptions:

- Genotypes are diploid, phased, and non-missing. Missing, unphased, or non-diploid calls are rejected.
- FLARE `AN1` and `AN2` are scalar, non-missing integer FORMAT fields, encoded as `0..n_ancestries-1`.
- Only the FLARE hardcall ancestry fields are used. REF/ALT/GT in the LAI file are not interpreted as genotype data.
- Each LAI record represents the ancestry interval from the previous LAI position plus one through this LAI position, inclusive: `(previous_lai_pos, current_lai_pos]`.
- The first LAI record on a contig covers `1..current_lai_pos`.
- The last LAI hardcall on a contig is extended to the last genotype position on that contig, so the chromosome tail is not dropped when the LAI file lacks a terminal sentinel marker.
- Adjacent LAI intervals with identical haplotype ancestry masks are merged into one ancestry block.
- Duplicate LAI positions on the same contig are allowed; the last record at that position defines the ancestry hardcall for the interval ending there.
- Genotype and FLARE sample IDs must be identical and in the same order.
- Genotype and FLARE records are sorted by compatible contig order.
- FLARE can be a subset of genotype sites. Genotype contigs with no LAI records are skipped.

Implementation notes:

- Common variants are written as variant-major dense ALT haplotype bitsets.
- Rare variants are written as packed carriers: `uint32_t pos_index` plus `uint32_t anc_hap`.
- `anc_hap` stores ancestry in the high 5 bits and haplotype ID in the low 27 bits.
- Marker streams (`*.mks`) are binary, little-endian, and use length-prefixed strings, so allele/id lengths can vary naturally.
- Offset indexes (`*.idx`) are fixed-width binary tables pointing into both the matching `.mks` stream and the payload `.bin` file. Magic headers are `TRCMIDX2`, `TRRAIDX2`, and `TRANIDX2`.
- Common `.idx` records store `common_index`, `global_variant_index`, `mks_offset`, and `geno_offset`.
- Rare `.idx` records store `rare_index`, `global_variant_index`, `mks_offset`, `carrier_offset`, and `n_carriers`.
- Ancestry block `.idx` records store `block_id`, `mks_offset`, and `anc_offset`.
- Multi-allelic records are split logically by ALT allele, and each split ALT receives its own `global_variant_index`.
- Structural variants use the VCF start position (`POS`) as the marker position for ancestry-block assignment and chunk ownership. The converter does not interpret `INFO/END` or `SVLEN`.
- `<prefix>.samples` and `<prefix>.meta` are small sidecars used by readers and VCF roundtrip tooling.

Capacity limits:

- Offsets are `uint64_t`, so each `.mks` or payload `.bin` file can be addressed up to `2^64 - 1` bytes, about `16 EiB`.
- `global_variant_index` is `uint32_t`, supporting up to about `4.29B` split variants per prefix.
- `block_id` is `uint32_t`, supporting up to about `4.29B` ancestry blocks per prefix.
- `common_index` and `rare_index` are `uint64_t`.
- `hap_id` uses 27 bits inside `anc_hap`, supporting up to `134,217,728` haplotypes, or about `67M` diploid samples.
- Ancestry uses 5 bits inside `anc_hap`, supporting up to `32` ancestry labels.
- For AoU-scale data, a whole genome with about `1.6B` variants is under the `uint32_t` split-variant limit, and per-chromosome prefixes are substantially smaller.

Roundtrip scope:

- `tractor_hybrid_to_vcf` reconstructs a phased, split-biallelic BGZF-compressed VCF from the packed genotype bits/carriers and writes a tabix `.tbi` index.
- `tractor_hybrid_extract_region` writes a new packed prefix for variants with `start <= POS <= end`, renumbering marker indexes and clipping overlapping ancestry block metadata to the requested region.
- This verifies that split ALT haplotypes survive the packed representation.
- It does not restore original unsplit multi-allelic rows.
- Missing genotype and LAI calls are rejected by the forward converter, so they do not appear in roundtrip output.
- Genotype contigs skipped because they have no LAI records are not present in the packed files, so they cannot be emitted by the reverse converter.

Threshold estimation:

- `estimate_mac_threshold` reads phased diploid GT, splits multi-allelic sites by ALT, and builds a MAC distribution.
- `--sample-every N` gives a deterministic partial-site estimate by using the first VCF record and then every Nth record; the default `1` uses the whole input.
- Dense payload per split variant is `8 * ceil(2 * n_samples / 64)` bytes.
- Sparse payload per split variant is `8 * MAC` bytes.
- Storage-only break-even is therefore approximately `ceil(2 * n_samples / 64)`.
- The optional query model scores dense variants by scanned `uint64_t` words and sparse variants by scanned carriers.

The current converter intentionally stays close to the design document while fixing several prototype hazards:

- per-ALT `pos_index` is written correctly for multi-allelic sites;
- ancestry blocks are not reused across chromosomes;
- contig advancement uses the genotype header order instead of lexicographic chromosome strings;
- ancestry lookup for carriers is `O(1)` from the active FLARE block state.
