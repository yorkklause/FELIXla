# FELIXla / tractor-hybrid-tools

FELIXla is a haplotype-resolved, local-ancestry-aware genotype storage format.
It stores phased ALT alleles and inferred local ancestries under the same
explicit haplotype index, so ancestry-specific dosages can be returned by
direct query instead of being reconstructed from VCF and LAI files for every
phenotype run.

This repository currently provides the FELIXla v0 layout through the existing
`tractor_hybrid` packed files, plus standalone tools and SAIGE-TRACTOR adapter
code for reading that layout.

The main converter takes a phased genotype VCF/BCF plus FLARE local ancestry
VCF/BCF and writes a packed prefix that can be read back as VCF, subset by
region, or used by the patched SAIGE-TRACTOR step2 backend.

## Contents

- `src/`: standalone C++ tools.
- `saige_step2_adapter/`: SAIGE/SAIGE-TRACTOR step2 adapter source and patches.
- `prebuilt/linux-x86_64-static/`: dependency-free Linux x86_64 binaries.
- `resources/shapeit5_chunks/b38_4cM/`: GRCh38 autosome 4 cM chunk files.
- `scripts/felixla`: unified FELIXla command dispatcher.
- `scripts/shapeit_chunks_to_tractor_args.sh`: SHAPEIT5 chunk adapter for
  `flare_subset_to_tractor_hybrid`.
- `docker/`: tools-only and SAIGE-TRACTOR Docker builds.

## Quick Start

Build locally:

```bash
make
bin/felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --n-ancestries 5 \
  --make-felixla \
  --out hybrid/chr22
bin/felixla --felixla hybrid/chr22 --query chr22:160500 --ref A --alt G --nonzero-only
```

Use the prebuilt Linux x86_64 binaries:

```bash
export PATH="$PWD/prebuilt/linux-x86_64-static:$PATH"
flare_subset_to_tractor_hybrid genotype.phased.vcf.gz flare.anc.vcf.gz 5 512 hybrid/chr22
```

Use the tools-only Docker image:

```bash
docker pull kyuan1024/tractor-hybrid-tools:latest
docker run --rm -v "$PWD:/data" -w /data kyuan1024/tractor-hybrid-tools:latest \
  flare_subset_to_tractor_hybrid genotype.phased.vcf.gz flare.anc.vcf.gz 5 512 hybrid/chr22
```

The static binaries are intended for local VCF/BCF paths. Remote URL/S3/GCS
support is disabled in the bundled htslib so the binaries do not need conda or
shared runtime libraries.

## Common Workflows

### PLINK-Style FELIXla CLI

`felixla` is the preferred entrypoint and follows the PLINK convention: provide
input flags, action flags, and `--out`. Existing standalone commands remain
available for older scripts, and `felixla` also accepts old command names as
compatibility subcommands.

| PLINK-style command | Compatible standalone command |
| --- | --- |
| `felixla --phase-vcf ... --flare-vcf ... --make-felixla --out ...` | `flare_subset_to_tractor_hybrid ...` |
| `felixla --phase-vcf ... --rfmix-msp ... --make-felixla --out ...` | `rfmix_msp_to_tractor_hybrid ...` |
| `felixla --tractor-dosage-vcf ... --make-felixla --out ...` | `tractor_dosage_vcf_to_hybrid ...` |
| `felixla --felixla ... --export vcf --out ...` | `tractor_hybrid_to_vcf ...` |
| `felixla --felixla ... --region ... --make-felixla --out ...` | `tractor_hybrid_extract_region ...` |
| `felixla --felixla ... --query ...` | `felixla_query ...` |
| `felixla --vcf ... --admixture` | `calc_tractor_admixture ...` |
| `felixla --compare-vcfs A B` | `compare_vcfs A B` |
| `felixla --shapeit-args ...` | `shapeit_chunks_to_tractor_args ...` |
| `felixla --recommend-mac-threshold --n-samples N` | formula helper |

For compatibility, aliases such as `felixla flare_subset_to_tractor_hybrid ...`
and `felixla tractor_hybrid_to_vcf ...` dispatch to the same tools.

### FLARE VCF To Packed Hybrid

```bash
felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --n-ancestries 5 \
  --make-felixla \
  --out hybrid/chr22
```

Arguments are:

1. Phased genotype VCF/BCF.
2. FLARE local ancestry VCF/BCF with scalar integer `FORMAT/AN1` and
   `FORMAT/AN2`.
3. Number of ancestry labels, max `32`.
4. Output prefix.
5. Optional MAC threshold override: variants with `MAC <= threshold` are stored sparse.
   If `--mac-threshold` is omitted, `felixla` uses
   `auto`, and the htslib-backed writer computes `ceil(n_samples / 32)` from
   the input header. See [Choosing Mac Threshold](#choosing-mac-threshold).
6. Optional `chr:start-end` region.

For one region:

```bash
felixla \
  --phase-vcf phase/chr22.phased.vcf.gz \
  --flare-vcf flare/chr22.flare.vcf.gz \
  --n-ancestries 5 \
  --region chr22:1-50000000 \
  --make-felixla \
  --out hybrid/chr22.chunk0001
```

For PLINK-like sample and site filtering while writing from FLARE:

```bash
felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --n-ancestries 5 \
  --keep samples.keep \
  --extract sites.pvar \
  --make-felixla \
  --out hybrid/chr22.subset
```

`--keep` is one sample ID per line. Retained samples are written in genotype
VCF order. `--extract` accepts a PLINK2 `.pvar`-like file or VCF-like file and
uses the first variant columns `CHROM POS ID REF ALT`; VCF `QUAL/FILTER/INFO`,
`FORMAT`, and sample columns are ignored. Multi-allelic `ALT` values may be
comma-separated. REF must be known in the list, and a REF mismatch against the
genotype VCF is a fatal error.

### SHAPEIT5 4 cM Chunk Conversion

This repo vendors SHAPEIT5-style GRCh38 4 cM autosome chunks in
`resources/shapeit5_chunks/b38_4cM/`.

The chunk files have this layout:

```text
chunk_index    chrom    buffered_region    core_region    ...
```

Use `core_region` column 4 for conversion. Column 3 is a phasing buffer region
and overlaps neighboring chunks, so using it here would duplicate variants at
chunk boundaries.

Generate argument rows:

```bash
felixla \
  --shapeit-args \
  --chunks-dir resources/shapeit5_chunks/b38_4cM \
  --chrom-style chr \
  --phase-template 'phase/{chrom}.phased.vcf.gz' \
  --flare-template 'flare/{chrom}.flare.vcf.gz' \
  --n-ancestries 5 \
  --n-samples 100000 \
  --out-prefix-template 'hybrid/{chrom}.shapeit4cM.chunk{chunk0}' \
  > flare_subset.shapeit4cm.args.tsv
```

Rows match the converter argument order:

```text
source_phase_vcf    source_flare_vcf    n_ancestries    mac_threshold    out_prefix    region
```

Run in parallel:

```bash
xargs -a flare_subset.shapeit4cm.args.tsv -n 6 -P 8 flare_subset_to_tractor_hybrid
```

The bundled chunk files use chromosome names `1..22`. Use `--chrom-style chr`
for VCF contigs named `chr1..chr22`; use `--chrom-style keep` or omit it for
contigs named `1..22`.

Inside the tools Docker image:

```bash
docker run --rm -v "$PWD:/data" -w /data kyuan1024/tractor-hybrid-tools:latest \
  shapeit_chunks_to_tractor_args \
    --chunks-dir /opt/tractor-hybrid-tools/resources/shapeit5_chunks/b38_4cM \
    --chrom-style chr \
    --phase-template 'phase/{chrom}.phased.vcf.gz' \
    --flare-template 'flare/{chrom}.flare.vcf.gz' \
    --n-ancestries 5 \
    --mac-threshold 512 \
    --out-prefix-template 'hybrid/{chrom}.shapeit4cM.chunk{chunk0}' \
  > flare_subset.shapeit4cm.args.tsv
```

### RFMix MSP To Packed Hybrid

```bash
felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --rfmix-msp rfmix.msp.tsv.gz \
  --n-ancestries 5 \
  --make-felixla \
  --out hybrid/chr22
```

The MSP converter expects haplotype columns in VCF sample order, named like
`sample.0` and `sample.1`. It accepts exact chromosome names or simple
`chr`/non-`chr` equivalents, for example MSP `22` with VCF `chr22`.

Adjacent MSP rows often share an endpoint. The converter treats the first
interval as including `spos`, then assigns later shared `epos`/`spos`
boundaries to the previous interval.

### TRACTOR Dosage VCF To Packed Hybrid

For SAIGE-TRACTOR-style VCF/BCF files that already contain hardcall
`DS1..DSk` and `ANC1..ANCk` FORMAT fields:

```bash
felixla \
  --tractor-dosage-vcf tractor.step2.input.vcf.gz \
  --n-ancestries 5 \
  --make-felixla \
  --out hybrid/chr22
```

This reconstructs canonical haplotypes from ancestry counts and
ancestry-specific ALT dosages. It is equivalent for SAIGE-TRACTOR step2 fields
but does not recover the original phased haplotype order. Fractional imputed
dosages are rejected.

### Roundtrip And Region Extraction

Convert a packed prefix back to split-biallelic VCF:

```bash
felixla --felixla hybrid/chr22 --export vcf --out hybrid/chr22.roundtrip.vcf.gz
```

Extract a packed sub-region:

```bash
felixla \
  --felixla hybrid/chr22 \
  --region chr22:16000000-17000000 \
  --make-felixla \
  --out hybrid/chr22.16_17mb
```

The extractor also accepts split arguments:

```bash
felixla --felixla in_prefix --chr chr22 --from-bp 16000000 --to-bp 17000000 --make-felixla --out out_prefix
```

### Direct FELIXla Dosage Query

Query ancestry-specific ALT dosages from a packed prefix:

```bash
felixla --felixla hybrid/chr22 --query chr22:160500 --ref A --alt G
felixla --felixla hybrid/chr22 --query chr22:160500 --ref A --alt G --nonzero-only
felixla --felixla hybrid/chr22 --global-index 123456 --nonzero-only
```

Output is one row per matched split-biallelic variant and sample:

```text
global_variant_index    chr    pos    id    ref    alt    sample    DSALL    DS1 ... DSk
```

`DS1` corresponds to ancestry code `0`, `DS2` to ancestry code `1`, and so on.
For common variants, `felixla_query` computes `DSk` by intersecting the dense
ALT haplotype bitset with the ancestry `k` haplotype mask. For rare variants,
it uses the packed carrier list, where each carrier stores both haplotype ID
and ancestry code.

### Choosing Mac Threshold

For storage-optimal sparse/dense packing, the threshold does not require
scanning the VCF. It depends only on the sample count:

```text
mac_threshold = ceil(n_samples / 32)
```

For PLINK-style `felixla --make-felixla` commands, this is the default. The
wrapper passes `auto`, and the htslib-linked writer reads `n_samples` from the
input VCF/BCF header. Use `--mac-threshold` only to override it. Use
`--n-samples` for template-only jobs such as `--shapeit-args`, where there is
no concrete VCF/BCF header to read.

Equivalent integer shell formula:

```bash
n_samples=100000
felixla --recommend-mac-threshold --n-samples "$n_samples"
```

Reason: dense payload per split variant is
`8 * ceil(2 * n_samples / 64)` bytes, while sparse payload is `8 * MAC` bytes.
The break-even point is therefore `ceil(n_samples / 32)`. Scanning a VCF can
estimate how many variants fall below that threshold, but it is not needed to
choose the threshold itself.

### VCF Comparison

```bash
compare_vcfs expected.vcf.gz observed.vcf.gz --split-multiallelic
```

`--split-multiallelic` compares both inputs after logically splitting each ALT
allele, which is useful for checking a packed roundtrip against an original
multi-allelic VCF.

## SAIGE-TRACTOR Step2 Adapter

Adapter files:

- `saige_step2_adapter/TractorHybrid.hpp`
- `saige_step2_adapter/TractorHybrid.cpp`
- `saige_step2_adapter/Main.cpp.patch`
- `saige_step2_adapter/Geno.R.patch`
- `saige_step2_adapter/SAIGE_SPATest_Tractor.R.patch`
- `saige_step2_adapter/SAIGE_SPATest_Marker.R.patch`
- `docker/saigetractor-hybrid/Dockerfile`
- `docker/saigetractor-hybrid/patch_saigetractor_hybrid.py`

Build the local SAIGE-TRACTOR image:

```bash
docker build --platform linux/amd64 \
  -t saigetractor:1.4.9-tractor-hybrid \
  -f docker/saigetractor-hybrid/Dockerfile .
```

Run step2 with a packed prefix:

```bash
docker run --rm --platform linux/amd64 \
  -v /path/to/data:/data \
  saigetractor:1.4.9-tractor-hybrid \
  step2_SPAtests.R \
    --tractorHybridPrefix=/data/chr22 \
    --chrom=chr22 \
    --is_admixed=TRUE \
    --number_of_ancestry=5 \
    --markers_per_chunk=1000 \
    --GMMATmodelFile=/data/null.rda \
    --varianceRatioFile=/data/varianceRatio.txt \
    --SAIGEOutputFile=/data/chr22.tractor.out
```

`--markers_per_chunk=1000` is a practical setting for step2 chunk traversal.
SAIGE-TRACTOR 1.4.9 enforces `1000` as the minimum for single-variant tests.

The adapter adds a narrow `tractor_hybrid` genotype backend for step2:

- Ordinary SAIGE paths receive total ALT dosage `0/1/2`.
- SAIGE-TRACTOR admixed paths receive `ANC1..ANCk`, `DS1..DSk`, `DSALL`, `DS`,
  or `GT`-like fields.
- The reader supports streaming traversal to avoid loading a full marker stream
  into R memory.
- Repeated same-marker ancestry fields are cached, which matters because
  SAIGE-TRACTOR asks for several fields for each marker.

The local image also installs standalone tools in `/scripts/bin`.

Frozen Docker Hub image for the current stable `.2` line:

```text
kyuan1024/saigetractor:1.4.9-tractor-hybrid.2
sha256:9f872fbd2df8e2ee57a427fbbb28b529fb4165ef823e4894750b17f60910fb8e
```

## Build And Test

Dynamic local build:

```bash
make
make test
```

Extended keep/extract regression:

```bash
make test-intense
```

`test-intense` synthesizes multi-chromosome, multi-allelic phased VCF and
FLARE inputs, checks full roundtrips against the legacy positional converter,
verifies `--keep`, `--extract`, `--region`, and `felixla --query` against a
known truth table, and asserts that malformed keep/extract files fail with
specific errors.

Static-style build:

```bash
make static
make test-static
```

Build dependency:

- htslib headers and library are required.
- On macOS with Homebrew: `brew install htslib`.
- If htslib is installed in a non-standard location:

```bash
make HTSLIB_CFLAGS="-I/path/to/htslib/include" HTSLIB_LIBS="-L/path/to/htslib/lib -lhts"
```

`make static` writes binaries to `bin-static/` and links `libhts.a` directly
when a static htslib archive is available. On Linux, fully static executables
can be requested with:

```bash
make static STATIC_FULLY=1
```

## Input Assumptions

- Genotypes are diploid, phased, and non-missing.
- FLARE `AN1` and `AN2` are scalar, non-missing integer FORMAT fields encoded
  as `0..n_ancestries-1`.
- Only FLARE hardcall ancestry fields are used. REF/ALT/GT in the LAI file are
  not interpreted as genotype data.
- Genotype and LAI sample IDs must be identical and in the same order.
- Genotype and LAI records must be sorted by compatible contig order.
- FLARE can be a subset of genotype sites.
- Genotype contigs with no LAI records are skipped.
- Multi-allelic genotype records are split logically by ALT allele.
- Structural variants use VCF `POS` as the marker position; `INFO/END` and
  `SVLEN` are not interpreted.

FLARE interval convention:

- Each LAI record represents `(previous_lai_pos, current_lai_pos]`.
- The first LAI record on a contig covers `1..current_lai_pos`.
- The last LAI hardcall on a contig extends to the last genotype position on
  that contig.
- Adjacent LAI intervals with identical haplotype ancestry masks are merged.
- Duplicate LAI positions on the same contig are allowed; the last record at
  that position defines the interval ending there.

Region/chunk convention:

- Region strings are `chr:start-end`.
- Coordinates are 1-based inclusive.
- A variant belongs to a chunk when `start <= POS <= end`.
- Adjacent chunks should start at previous `end + 1`.

## FELIXla v0 Packed Output

FELIXla v0 is binary-compatible with the current `tractor_hybrid` packed
layout. The compatibility is intentional: existing converters and the
SAIGE-TRACTOR adapter can keep using the same files while FELIXla-facing tools
and documentation stabilize around the storage abstraction.

Each packed prefix writes:

- `<prefix>.common.geno.bin`
- `<prefix>.common.variant.mks`
- `<prefix>.common.variant.idx`
- `<prefix>.rare.carrier.bin`
- `<prefix>.rare.variant.mks`
- `<prefix>.rare.variant.idx`
- `<prefix>.ancblock.bin`
- `<prefix>.ancblock.mks`
- `<prefix>.ancblock.idx`
- `<prefix>.samples`
- `<prefix>.meta`

Storage notes:

- Haplotype IDs are explicit: sample `i` uses haplotypes `2*i` and `2*i+1`.
- Common variants are variant-major dense ALT haplotype bitsets.
- Rare variants are packed carriers: `uint32_t pos_index` plus
  `uint32_t anc_hap`.
- `anc_hap` stores ancestry in the high 5 bits and haplotype ID in the low
  27 bits.
- Local ancestry is stored as merged ancestry blocks. Each block contains one
  dense haplotype mask per ancestry, in the same haplotype coordinate system as
  the genotype bitsets.
- Marker streams (`*.mks`) are binary, little-endian, and use length-prefixed
  strings.
- Offset indexes (`*.idx`) are fixed-width binary tables pointing into the
  matching marker stream and payload file.
- Marker streams and offset indexes start with an 8-byte magic header used as a
  format identifier/version guard:

  | File kind | Magic | Meaning |
  | --- | --- | --- |
  | common variant marker stream | `TRCMMKS1` | TRACTOR common `*.mks`, version 1 |
  | rare variant marker stream | `TRRAMKS1` | TRACTOR rare `*.mks`, version 1 |
  | ancestry block marker stream | `TRANMKS1` | TRACTOR ancestry `*.mks`, version 1 |
  | common variant offset index | `TRCMIDX2` | TRACTOR common `*.idx`, version 2 |
  | rare variant offset index | `TRRAIDX2` | TRACTOR rare `*.idx`, version 2 |
  | ancestry block offset index | `TRANIDX2` | TRACTOR ancestry `*.idx`, version 2 |

- `<prefix>.ancblock.mks` stores ancestry block metadata after `TRANMKS1`.
  Each record is `uint32_t block_id`, length-prefixed `chr`,
  `int64_t start_pos`, `int64_t end_pos`, and `uint64_t anc_offset` into
  `<prefix>.ancblock.bin`.

Capacity limits:

- `global_variant_index` is `uint32_t`, supporting about `4.29B` split variants
  per prefix.
- `block_id` is `uint32_t`, supporting about `4.29B` ancestry blocks per
  prefix.
- Payload offsets are `uint64_t`.
- `hap_id` supports about `134M` haplotypes, or about `67M` diploid samples.
- Ancestry uses 5 bits, supporting up to `32` ancestry labels.

## Roundtrip Scope

`tractor_hybrid_to_vcf` reconstructs a phased, split-biallelic BGZF-compressed
VCF and writes a tabix index. It verifies that split ALT haplotypes survive the
packed representation, but it does not restore original unsplit multi-allelic
rows. Missing genotype and LAI calls are rejected by the forward converters, so
they do not appear in roundtrip output.

## Progress And Timing

Forward converters print progress to stderr. When an input `.csi` or `.tbi`
index exposes record statistics, progress includes `scanned records / total
records` and a percentage; otherwise it reports scanned records only.

The patched SAIGE-TRACTOR step2 scripts print chunk-level timing diagnostics
for `tractor_hybrid` input. These diagnostics are instrumentation only and do
not change score, SPA, Firth, variance-ratio, filtering, or output
calculations.
