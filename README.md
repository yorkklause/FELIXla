# FELIXla

**FELIX** (**F**ull-cohort **E**fficient **L**ocal ancestry-**I**ntegrated
mi**X**ed-model framework) is a framework for scalable local-ancestry-aware
association analysis. **FELIXla** is its C++ based command line tool and binary
storage format for haplotype-resolved, local-ancestry-aware genotype data. It
stores phased ALT alleles and inferred local ancestries under the same explicit
haplotype index, so ancestry-specific dosages can be returned by direct query
instead of being recomputed from phased genotype and local ancestry files for
every phenotype.

The ancestry-specific dosage for ancestry `k` at a variant is represented as a
bitwise intersection between the ALT haplotype vector and the ancestry `k`
haplotype mask. Local ancestry is piecewise constant along a chromosome, so
FELIXla stores ancestry transition blocks rather than one ancestry label per
variant. Common variants are stored as dense bit vectors, and rare variants are
stored as sparse carrier lists.

FELIXla v0 is binary-compatible with the `tractor_hybrid` packed layout used by
this repository. The compatibility is intentional: the same packed prefix can
be queried directly, converted back to split-biallelic VCF, or subset by
region.

## Citation

A manuscript describing FELIXla is in preparation. If you use this software in
a published analysis before a formal citation is available, please cite this
GitHub repository and record the commit hash used in the analysis.

## Getting Started

- Clone this repository using the following git command:

    `git clone https://github.com/yorkklause/FELIXla.git`

    Alternatively, download the source files from the GitHub website
    (https://github.com/yorkklause/FELIXla).

- **FELIXla** requires a C++17 compiler and **htslib** headers and library.

    On macOS with Homebrew:

    `brew install htslib`

    On Linux, install htslib through your system package manager or provide
    explicit compiler and linker flags as shown below.

- Compile FELIXla using the following command:

    `make`

- Once the compiling is done, the single executable `felixla` will be found in
  the `bin` folder. Type

    `bin/felixla --help` or `bin/felixla -h`

    to print a list of command-line options.

- FELIXla is a complete binary entry point. It does not invoke separate FELIXla
  helper executables at runtime.

- A Docker image is published through GitHub Container Registry:

    `docker pull ghcr.io/yorkklause/felixla:latest`

    The image entry point is `felixla`, so running

    `docker run --rm ghcr.io/yorkklause/felixla:latest --help`

    will print the FELIXla command-line help.

- A Linux x86_64 static binary is attached to GitHub Releases as
  `felixla-linux-x86_64-static`. The release binary is built as a single
  statically linked executable against a local-file htslib build without
  libcurl remote-URL support.

## Using FELIXla

The preferred entry point is the PLINK-style `felixla` command:

```
felixla \
  --phase-vcf PHASED_VCF \
  --flare-vcf FLARE_VCF \
  --n-ancestries N_ANCESTRIES \
  --make-felixla \
  --out OUT_PREFIX \
  [ --mac-threshold MAC_THRESHOLD ] \
  [ --region CHR:START-END ] \
  [ --keep SAMPLE_LIST ] \
  [ --extract SITE_LIST ]
```

- PHASED_VCF (required): Full path and filename of a phased diploid genotype
  VCF/BCF. Genotypes must be phased, diploid, and non-missing. Multi-allelic
  records are split logically by ALT allele.

- FLARE_VCF (required): Full path and filename of a FLARE local ancestry
  VCF/BCF. It must contain scalar integer `FORMAT/AN1` and `FORMAT/AN2` fields
  encoded as `0..N_ANCESTRIES-1`. Sample IDs must be identical to the genotype
  VCF and in the same order.

- N_ANCESTRIES (required): Number of local ancestry labels. The packed format
  stores ancestry codes in 5 bits, so at most 32 labels are supported.

- OUT_PREFIX (required): Prefix of the FELIXla output files.

- MAC_THRESHOLD (optional): Sparse/dense minor allele count threshold. Variants
  with `MAC <= MAC_THRESHOLD` are stored as sparse carrier lists; variants above
  the threshold are stored as dense bit vectors. Default is `auto`, which uses
  `ceil(n_samples / 32)` from the retained sample count.

- CHR:START-END (optional): 1-based inclusive region to convert. A variant
  belongs to the region when `START <= POS <= END`.

- SAMPLE_LIST (optional): A file containing one sample ID per line, analogous to
  PLINK `--keep`. Retained samples are written in genotype VCF order.

- SITE_LIST (optional): A PLINK2 `.pvar`-like file or VCF-like file used like
  PLINK `--extract`. FELIXla reads the first variant columns `CHROM POS ID REF
  ALT`; VCF `QUAL`, `FILTER`, `INFO`, `FORMAT`, and sample columns are ignored.
  The `ID` column is ignored: retained alleles are matched and de-duplicated
  only by `CHROM`, `POS`, `REF`, and `ALT`. Multi-allelic `ALT` values may be
  comma-separated. REF must be known, and a REF mismatch against the genotype
  VCF is a fatal error. When the genotype VCF/BCF has a tabix/CSI index,
  `--extract` first finds the first and last requested position on each
  chromosome, seeks to those bounded intervals, and then linearly scans inside
  them before exact allele-level filtering. Without an index, FELIXla falls
  back to a streaming scan and prints a warning.

The same binary also dispatches to compatibility subcommands:

```
felixla --phase-vcf PHASED_VCF --rfmix-msp MSP_FILE --n-ancestries N --make-felixla --out OUT_PREFIX
felixla --tractor-dosage-vcf DOSAGE_VCF --n-ancestries N --make-felixla --out OUT_PREFIX
felixla --felixla PREFIX --export vcf --out OUTPUT_VCF
felixla --felixla PREFIX --query CHR:POS --ref REF --alt ALT [ --nonzero-only ]
felixla --felixla PREFIX --region CHR:START-END --make-felixla --out OUT_PREFIX
felixla --vcf TRACTOR_VCF --admixture --out ADMIXTURE_TSV
felixla --compare-vcfs EXPECTED_VCF OBSERVED_VCF [ --split-multiallelic ]
felixla --recommend-mac-threshold --n-samples N_SAMPLES
```

For compatibility with older command lines, the original command names are
accepted as `felixla` subcommands, for example
`felixla flare_subset_to_tractor_hybrid ...` and
`felixla tractor_hybrid_to_vcf ...`. They are not installed as separate
executables.

## Output

A FELIXla packed prefix writes the following files:

- `<prefix>.common.geno.bin`: dense ALT haplotype bit vectors for common
  variants.

- `<prefix>.common.variant.mks`: marker records for common variants.

- `<prefix>.common.variant.idx`: fixed-width offsets into the common marker and
  genotype payload files.

- `<prefix>.rare.carrier.bin`: sparse carrier records for rare variants. Each
  carrier stores both haplotype ID and ancestry code.

- `<prefix>.rare.variant.mks`: marker records for rare variants.

- `<prefix>.rare.variant.idx`: fixed-width offsets into the rare marker and
  carrier payload files.

- `<prefix>.ancblock.bin`: ancestry-specific haplotype masks for merged local
  ancestry blocks.

- `<prefix>.ancblock.mks`: ancestry block coordinates and offsets.

- `<prefix>.ancblock.idx`: fixed-width offsets into the ancestry marker and
  payload files.

- `<prefix>.samples`: retained sample IDs, one sample per line.

- `<prefix>.meta`: format metadata, input provenance, sample count, ancestry
  count, word count, rare threshold, and optional region/filter provenance.

The split variant marker streams record `chr`, `pos`, split-biallelic `id`,
`ref`, `alt`, ALT allele index, global variant index, and MAC. Marker streams
and indexes are little-endian binary files with 8-byte magic headers.

A direct dosage query writes one row per matched split-biallelic variant and
sample:

```text
global_variant_index    chr    pos    id    ref    alt    sample    DSALL    DS1 ... DSk
```

`DS1` corresponds to ancestry code `0`, `DS2` to ancestry code `1`, and so on.
For common variants, `felixla` computes `DSk` by intersecting the dense ALT bit
vector with the ancestry `k` haplotype mask. For rare variants, it uses the
sparse carrier list directly.

## Example

The following commands build the binary, create a FELIXla prefix from the tiny
phased genotype and FLARE local ancestry fixtures, query one ancestry-specific
dosage vector, and export the packed data back to split-biallelic VCF:

```
make

bin/felixla \
  --phase-vcf testdata/tiny.genotypes.vcf \
  --flare-vcf testdata/tiny.flare.vcf \
  --n-ancestries 2 \
  --make-felixla \
  --out example/tiny

bin/felixla \
  --felixla example/tiny \
  --query chr1:160 \
  --ref A \
  --alt T

bin/felixla \
  --felixla example/tiny \
  --export vcf \
  --out example/tiny.roundtrip.vcf.gz
```

Sample and site filtering can be applied at conversion time:

```
bin/felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --n-ancestries 5 \
  --keep samples.keep \
  --extract sites.pvar \
  --make-felixla \
  --out hybrid/chr22.subset
```

The same command can be run through Docker by binding the working directory:

```
docker run --rm -v "$PWD":/data -w /data ghcr.io/yorkklause/felixla:latest \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --n-ancestries 5 \
  --keep samples.keep \
  --extract sites.pvar \
  --make-felixla \
  --out hybrid/chr22.subset
```

The SHAPEIT/GLIMPSE chunk mode can generate converter argument rows for
chromosome-scale jobs:

```
bin/felixla \
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

## Testing

Run the standard smoke test:

```
make test
```

Run the extended keep/extract regression test:

```
make test-intense
```

`test-intense` synthesizes multi-chromosome, multi-allelic phased genotype and
FLARE inputs, compares the PLINK-style command against the compatibility
subcommand path, verifies `--keep`, `--extract`, indexed bounded extract
scanning, `--region`, and `felixla --query` against a known truth table, and
checks that malformed keep/extract files fail with specific errors.

## Build Notes

If htslib is installed in a non-standard location, pass explicit flags:

```
make HTSLIB_CFLAGS="-I/path/to/htslib/include" HTSLIB_LIBS="-L/path/to/htslib/lib -lhts"
```

A static-style build can be requested with:

```
make static
make test-static
```

`make static` writes `bin-static/felixla` and links `libhts.a` directly when a
static htslib archive is available. On Linux, a fully static executable can be
requested with:

```
make static STATIC_FULLY=1
```

The Docker image is built from `docker/felixla/Dockerfile`. On pushes to
`main`, GitHub Actions publishes a multi-architecture image for `linux/amd64`
and `linux/arm64` at `ghcr.io/yorkklause/felixla:latest`.

On version tags, GitHub Actions also publishes
`felixla-linux-x86_64-static` and its SHA256 checksum to the corresponding
GitHub Release.

To build the image locally:

```
docker build -t felixla:local -f docker/felixla/Dockerfile .
```

## Input Assumptions

- Genotypes are diploid, phased, and non-missing.

- Genotype and local ancestry sample IDs are identical and in the same order.

- FLARE local ancestry uses non-missing integer `AN1` and `AN2` hard calls.

- FLARE can be a subset of genotype sites. The first LAI record on a contig
  covers `1..current_lai_pos`; later records represent
  `(previous_lai_pos, current_lai_pos]`; the final state extends to the last
  genotype position on that contig.

- Adjacent LAI intervals with identical haplotype ancestry masks are merged.

- Genotype contigs with no local ancestry records are skipped.

- Structural variants use VCF `POS` as the marker position; `INFO/END` and
  `SVLEN` are not interpreted.

## Support

Please direct questions or bug reports to Kai Yuan (kyuan@broadinstitute.org).
