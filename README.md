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

- Once compilation is done, the complete `felixla` executable and the separate
  `vcf_tbi_chunks` planning utility will be found in the `bin` folder. Type

    `bin/felixla --help` or `bin/felixla -h`

    to print a list of command-line options.

- FELIXla is a complete binary entry point. It does not invoke separate FELIXla
  helper executables at runtime.

- A Docker image is published through GitHub Container Registry:

    `docker pull ghcr.io/yorkklause/felixla:latest`

    The image entry point is `felixla`, so running

    `docker run --rm ghcr.io/yorkklause/felixla:latest --help`

    will print the FELIXla command-line help.

    The separate planner is also installed in the image:

    `docker run --rm --entrypoint vcf_tbi_chunks ghcr.io/yorkklause/felixla:latest --help`

- Linux x86_64 static binaries are attached to GitHub Releases as
  `felixla-linux-x86_64-static` and `vcf_tbi_chunks-linux-x86_64-static`.
  They are linked against a local-file htslib build without libcurl remote-URL
  support.

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
  encoded as `0..N_ANCESTRIES-1`. Sample IDs are matched by ID; FELIXla keeps
  the genotype/FLARE intersection in genotype VCF order.

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
  VCF is a fatal error. The site list may be unsorted: FELIXla sorts and groups
  it before conversion. Output follows genotype VCF record order and, within a
  multiallelic record, the genotype VCF ALT order. When the genotype VCF/BCF has a tabix/CSI index,
  `--extract` first finds the first and last requested position on each
  chromosome, seeks to those bounded intervals, and then linearly scans inside
  them before exact allele-level filtering. Without an index, FELIXla falls
  back to a streaming scan and prints a warning.

## Planning Parallel Chunks

`vcf_tbi_chunks` is a separate binary. It reads the phased VCF tabix index
directly through htslib and probes only the in-memory `.tbi`/`.csi` bins at
fixed-size boundaries. It does not open or decompress the phased VCF body and
does not scan variants. The VCF path is used to discover the sidecar index and
to render command templates; with an explicit `--tbi`, the VCF body does not
need to be locally readable. The utility does not invoke `tabix`, FELIXla, a
shell, or a job scheduler.

Generate 10 Mb chunks:

```bash
vcf_tbi_chunks \
  --phase-vcf genotype.phased.vcf.gz \
  --chunk-mb 10 \
  --out genotype.chunks.tsv
```

Chunk regions use the same 1-based inclusive coordinates as `felixla
--region`. Boundaries are aligned to the requested chunk length and adjacent
chunks begin at the previous end plus one:

```text
chr1:1-10000000
chr1:10000001-20000000
```

Thus there are no duplicated boundary variants. The first and last windows are
conservative index-derived bounds aligned to the requested chunk length. Tabix
bins are coarser than individual positions, so an edge window can be empty;
interior empty windows are also retained. This avoids any VCF data reads while
still guaranteeing that indexed records are covered. Use `--chunk-bp` when the
desired length is not an integer number of decimal megabases.

The manifest contains one row per task with global and per-contig chunk IDs,
`CHROM`, inclusive start/end, region text, conservative index-derived contig
bounds, and the indexed record count. `--chrom` may be repeated to select
contigs. The compatibility columns `contig_first_pos` and `contig_last_pos`
therefore contain aligned index bounds, not exact VCF record positions.

An optional command template writes a separate one-command-per-line file that
can be consumed by a scheduler or another parallel runner:

```bash
vcf_tbi_chunks \
  --phase-vcf genotype.phased.vcf.gz \
  --chunk-mb 10 \
  --out genotype.chunks.tsv \
  --command-template \
    'felixla --phase-vcf {phase_vcf_q} --flare-vcf flare.vcf.gz --n-ancestries 3 --region {region_q} --make-felixla --out out/{chrom}.chunk{chrom_chunk0}' \
  --commands-out genotype.commands.txt
```

Available template fields include `{phase_vcf}`, `{chrom}`, `{start}`, `{end}`,
`{region}`, `{global_chunk}`, `{global_chunk0}`, `{chrom_chunk}`, and
`{chrom_chunk0}`. `{phase_vcf_q}`, `{chrom_q}`, and `{region_q}` are
POSIX-shell-quoted forms.

### Remote object-store mounts

For region-parallel conversion, both the phased genotype VCF and the FLARE VCF
should be BGZF-compressed and have a readable `.tbi` or `.csi` index. FELIXla
seeks the genotype input to the requested/extracted interval. Its FLARE query
starts at `--region` START and reads forward only far enough to establish the
ancestry state; if no later FLARE record exists, it searches backward in
exponentially increasing windows for the preceding state. A warning containing
`scanning ... from the beginning` means indexed access was unavailable and the
job is taking the expensive compatibility path.

With `--extract` and `--region`, a plain PVAR/VCF site list is streamed but only
overlapping entries are retained in memory; unsorted lists remain supported.
Because a plain site list has no coordinate index, every worker still streams
that file once, so keep it on local storage when launching many regions.

Do not begin with dozens of region workers reading the same objects through a
Cloud Storage FUSE mount. Start with about 8 workers and increase concurrency
only while aggregate throughput improves. When possible, stage each
chromosome's two VCFs and their indexes once onto local SSD. If the mount is
under your control, Cloud Storage FUSE's
[file cache](https://docs.cloud.google.com/storage/docs/cloud-storage-fuse/file-caching)
and `--file-cache-cache-file-for-range-read=true` are intended for repeated
partial/random reads; size the cache to avoid
[cache thrashing](https://docs.cloud.google.com/storage/docs/cloud-storage-fuse/caching).

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

The large-sample packing microbenchmark is available separately from the
correctness suite:

```
make benchmark-pack
python3 tests/run_pack_benchmark.py --felixla bin/felixla --samples 400000 --records 100
python3 tests/run_pack_benchmark.py --felixla bin/felixla --samples 200000 --records 1000 --multiallelic-every 5 --drop-flare-every 100
python3 tests/run_pack_benchmark.py --felixla bin/felixla --samples 50000 --records 1000 --flare-every-record --flare-change-every 100
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

- Genotype and local ancestry sample IDs are matched by ID. If the two inputs
  differ, FELIXla keeps only their intersection and writes retained samples in
  genotype VCF order. `--keep` applies as an additional sample filter. Header
  intersections use sorted sample arrays, and retained sample columns are
  passed to htslib before record decoding.

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
