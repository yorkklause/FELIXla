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

FELIXla is distributed as a statically linked executable, so there is nothing
to install alongside it: no htslib, no shared libraries, no `LD_LIBRARY_PATH`.
Download it, mark it executable, and run it.

- Linux x86_64 static binaries are attached to every GitHub Release as
  `felixla-linux-x86_64-static` and `vcf_tbi_chunks-linux-x86_64-static`, with
  SHA256 checksums:

    ```bash
    curl -LO https://github.com/yorkklause/FELIXla/releases/latest/download/felixla-linux-x86_64-static
    chmod +x felixla-linux-x86_64-static
    ./felixla-linux-x86_64-static --help
    ```

    They carry their own htslib, built with
    [libdeflate](https://github.com/ebiggers/libdeflate) for fast BGZF
    decompression, and read local files only -- there is no libcurl remote-URL
    support.

- A Docker image is published through GitHub Container Registry, holding the
  same static executables:

    `docker pull ghcr.io/yorkklause/felixla:latest`

    The image entry point is `felixla`, so running

    `docker run --rm ghcr.io/yorkklause/felixla:latest --help`

    will print the FELIXla command-line help.

    The separate planner is also installed in the image:

    `docker run --rm --entrypoint vcf_tbi_chunks ghcr.io/yorkklause/felixla:latest --help`

- FELIXla is a complete binary entry point. It does not invoke separate FELIXla
  helper executables at runtime.

### Building from source

Building is only needed to develop FELIXla or to run it somewhere the released
binaries do not fit. It requires a C++17 compiler and **htslib** headers and
library:

- Clone this repository:

    `git clone https://github.com/yorkklause/FELIXla.git`

- Install htslib. On macOS with Homebrew, `brew install htslib`; on Linux, use
  the system package manager or pass explicit flags as shown under
  [Build Notes](#build-notes).

- Compile with `make`. The complete `felixla` executable and the separate
  `vcf_tbi_chunks` planning utility are written to `bin/`.

Prefer an htslib built with libdeflate; see
[Decompression Performance](#decompression-performance) for what it is worth.
To reproduce the self-contained release executables instead, see
[Build Notes](#build-notes).

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
  [ --extract SITE_LIST ] \
  [ --exclude SITE_LIST ] \
  [ --chr CONTIGS ] \
  [ --extract-bed BED_INTERVALS ] \
  [ --threads N_THREADS ]
```

- PHASED_VCF (required): Full path and filename of a phased diploid genotype
  VCF/BCF. Genotypes must be phased, diploid, and non-missing. The sex
  chromosomes are not supported: the packed format gives every sample exactly
  two haplotypes, so haploid genotypes cannot be represented, and a record on
  `X` or `Y` is refused with that explanation. Select the autosomes
  explicitly, for example `--chr 1-22`. Multi-allelic
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
  The `ID` column is ignored. Matching uses the source `CHROM` and `POS`, then
  compares each split `REF`/`ALT` after removing shared trailing and leading
  padding while retaining at least one base per allele. This permits equivalent
  padded multiallelic representations such as extract `ATT>AT` and genotype
  `AT>A`; output still uses the genotype VCF representation. Multi-allelic `ALT`
  values may be comma-separated. REF must be known. A raw REF difference is
  accepted only when at least one normalized split allele matches across all VCF
  records at that coordinate; otherwise it is a fatal REF mismatch. The site
  list may be unsorted: FELIXla sorts and groups it before conversion. Output
  follows genotype VCF record order and, within a multiallelic record, the
  genotype VCF ALT order. When the genotype VCF/BCF has a tabix/CSI index,
  `--extract` first finds the first and last requested position on each
  chromosome, seeks to those bounded intervals, and then linearly scans inside
  them before exact allele-level filtering. Without an index, FELIXla falls
  back to a streaming scan and prints a warning.

- `--exclude` SITE_LIST (optional): The same PVAR/VCF format and the same
  allele matching as `--extract`, naming alleles to drop rather than keep. It
  is applied after selection, so `--extract` and `--exclude` compose as they do
  in PLINK. Matching is allele-level: a multiallelic record keeps whatever ALTs
  the list does not name, and the record is dropped only when nothing remains.
  An entry matching no record in the input is silently inert. Unlike
  `--extract`, the list is not narrowed by `--region`, since a site outside the
  region is absent from the output regardless.

  Note that neither `--extract` nor `--exclude` uses the `ID` column: matching
  is on `CHROM`, `POS` and the normalized `REF`/`ALT`. A PLINK-style file of
  bare variant IDs is therefore rejected rather than silently misread.

- `--chr` CONTIGS (optional): Select whole contigs to convert.
  Comma-separated, with numeric ranges allowed, so `--chr 1-22` and
  `--chr chr1,chr2` both work. Each name is matched against the genotype VCF
  header as written, then with a `chr` prefix added, then removed, so the same
  command works on either naming convention. A name matching no contig in the
  header is a warning, not an error.

  Selecting whole contigs composes with the interval filters: with
  `--extract-bed`, intervals on unselected contigs are dropped; without it,
  each selected contig becomes one whole-contig interval. Combining a contig
  selection with a `--region` on an unselected contig is an error rather than
  a silently empty run.

- BED_INTERVALS (optional): A BED or gzip-compressed BED file. FELIXla reads
  the first three columns and ignores later columns, `track`/`browser` rows,
  blank lines, and comments. Coordinates follow the BED standard: 0-based,
  half-open `[START, END)`. Variant selection uses the VCF record position, so
  `chr1 99 100` selects a record at VCF position 100. Input intervals may be
  unsorted, duplicated, overlapping, or adjacent; FELIXla sorts and merges
  them. Every ALT of a selected multiallelic genotype record is retained unless
  `--extract` narrows the allele set. `--extract-bed`, `--extract`, and
  `--region` are intersected when combined. Ancestry blocks are clipped and
  split at BED gaps.

  With an indexed genotype VCF/BCF, FELIXla creates one bounded read span from
  the first through last retained BED interval on each chromosome, then applies
  the exact merged intervals with a linear cursor. It therefore does not issue
  one tabix request per BED row. This is intentional for large interval lists
  and remote object-store mounts. Without a usable index, it streams the input
  once and applies the same exact interval filter.

- N_THREADS (optional): Number of BGZF decompression threads shared by the
  genotype and FLARE readers. Default is `1`, which keeps one conversion in one
  core. Decompression is the largest single cost when both inputs are
  BGZF-compressed, so raising this is the most effective way to speed up a
  single conversion; it does not change the output. Parsing is serial and
  becomes the limit at roughly 8 threads. Threads and region workers compete
  for the same cores: use `--threads` when you run few conversions at a time,
  and leave it at `1` when you already saturate the machine with concurrent
  region jobs.

### Decompression Performance

Most of a conversion's work is BGZF decompression, and most of what is
decompressed is discarded: FELIXla reads `FORMAT/AN1` and `FORMAT/AN2` from the
FLARE VCF and nothing else, while a stock FLARE record also carries `GT`,
`ANP1` and `ANP2`. Dropping them is the single largest saving available and
does not change the output:

```bash
bcftools annotate -x '^FORMAT/AN1,FORMAT/AN2' -Oz -o flare.slim.vcf.gz FLARE_VCF
tabix -p vcf flare.slim.vcf.gz
```

FELIXla prints a one-line note when it sees a FLARE `FORMAT` carrying fields it
does not read.

Which htslib FELIXla is linked against matters as much as FELIXla's own code.
htslib uses
[libdeflate](https://github.com/ebiggers/libdeflate) when it is present at
build time, which is substantially faster than zlib for both inflate and the
per-block CRC. htslib's `configure` detects it automatically, so installing the
libdeflate development package before building htslib is enough; the published
static release binaries are built this way.

Verify a given build with `htsfile --version`, which lists libdeflate among its
features.

One conversion of a 50,000-sample, 2,000-variant BGZF fixture on one core of an
AMD EPYC 7742, each step leaving the output byte-identical:

| | wall |
|---|---|
| zlib htslib, stock FLARE | 6.1 s |
| libdeflate htslib, stock FLARE | 2.5 s |
| libdeflate htslib, stock FLARE, `--threads 8` | 1.0 s |
| libdeflate htslib, slim FLARE | 0.7 s |
| libdeflate htslib, slim FLARE, `--threads 4` | 0.4 s |

Slimming the FLARE input and linking against libdeflate each beat spending
eight cores on the stock input, and they compose, so prefer them before
`--threads`. Slim columns are also the shape the ancestry reader decodes eight
at a time, which is why the gain is larger than the file shrinks. After both,
roughly 70% of what is left is decompression, so `--threads` is what remains;
it is most useful when few conversions run at once.

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

## Concatenating FELIXla Chunks

After region jobs finish, concatenate complete FELIXla prefixes in genomic
order with the PLINK-style `--pmerge-list` interface:

```bash
felixla \
  --pmerge-list chr1.prefixes.txt \
  --make-felixla \
  --out chr1.merged
```

The list contains one prefix per line, without a component-file suffix:

```text
out/chr1.chunk0001
out/chr1.chunk0002
out/chr1.chunk0003
```

Blank lines and lines beginning with `#` are ignored. A trailing `.meta` is
also accepted. Paths are interpreted relative to the working directory. The
equivalent compatibility command is
`felixla concat chr1.prefixes.txt chr1.merged`.

To extract a different BED from each input and merge the retained records in
genomic coordinate order, use two tab-delimited columns on every data row:

```text
# FELIXla prefix<TAB>BED
out/source_B	beds/middle.bed
out/source_A	beds/outer.bed
```

List order does not determine record order in this mode. For example,
`outer.bed` may select `chr1:1-100` and `chr1:301-400`, while `middle.bed`
selects `chr1:101-300`; the merged output is ordered `1-400`. FELIXla reads
only the first three BED columns and interprets them as standard 0-based,
half-open intervals. Empty, malformed, duplicate, or overlapping intervals
within one BED are rejected; directly adjacent intervals are coalesced.

Prefix the BED path with `^` to select its complement, following the bcftools
targets convention:

```text
out/source_A	^beds/middle.bed
out/source_B	beds/middle.bed
```

For `^BED`, the selection universe is the ancestry-block coverage actually
present in that FELIXla prefix, not an inferred chromosome length. Plain BED
selections are likewise intersected with actual ancestry-block coverage.
Before writing output, FELIXla computes every effective selection and rejects
any overlap between inputs, reporting all affected prefixes. One-column and
two-column rows cannot be mixed in the same list.

Concatenation is a format-aware rewrite, not a bytewise `cat`. FELIXla remaps
global variant ordinals, common/rare local indexes, ancestry block IDs, marker
offsets, payload offsets, and sparse-carrier `pos_index` values. Prefixes must
have identical format version, sample order, sample/haplotype/word counts,
ancestry count, and rare threshold. Sample IDs and their order must match
exactly, including in BED-filtered mode. For a one-column list, regions and
actual records on the same chromosome must be non-overlapping and listed in
increasing coordinate order.

Concatenation runs in two strict phases. The first phase validates every listed
prefix without opening any output file. It does not stop after one bad prefix:
all problematic list entries and their detected errors are printed together,
and the merge is not started if any prefix fails. The checks
cover all 11 prefix files: required metadata, exact sample IDs and order, six
magic headers, fixed-width index records, marker/index agreement, contiguous
payload offsets, exact EOF, MAC versus dense popcount or sparse carrier count,
packed ancestry/haplotype ranges, ancestry-mask partitioning, and variant/block
ordering.

Only after the complete preflight passes does the second phase reread the
inputs, remap indexes, and write the merged data under a temporary prefix. In
two-column mode, ancestry blocks and common/rare variants are filtered and
coordinate-merged as bounded streams; the implementation neither materializes
all variants in memory nor invokes bcftools or another executable. The large
payloads are therefore read twice intentionally: once to establish that the
whole list is valid, then once to produce output. The temporary prefix is
published with `.meta` last, after the merge completes.

BED-filtered merging also simplifies ancestry blocks created by chunk or BED
boundaries. Two pieces are collapsed only when they are on the same chromosome,
their coordinates are directly adjacent, and every ancestry haplotype-mask byte
is identical. Pieces separated by even one unselected base remain distinct, so
an output block never bridges a BED gap. Metadata records the rule in
`concat_ancestry_block_simplification` and the number of removed boundaries in
`concat_ancestry_block_pieces_collapsed`.

Newly written FELIXla metadata includes `global_variants`, `common_variants`,
`rare_variants`, and `ancestry_blocks`. Concat cross-checks these declarations
when present. Older format-version-1 prefixes without the four count fields
remain supported; their counts are derived from and validated against the
binary marker/index streams.

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

A retained site costs about 41 bytes: contig names are interned, allele text is
held in one arena, and the per-site record is fixed-width. A 10-million-site
whole-genome list is therefore roughly 430 MB. That is per process, so budget
it against the worker count before raising concurrency, and narrow the list
with `--region` where possible -- sites outside the region are dropped as the
file is read and never held.

With `--extract-bed`, overlapping and adjacent intervals are merged first. The
genotype reader makes at most one indexed span query per retained chromosome
and performs exact BED membership checks while scanning that span. Sparse BED
files therefore avoid thousands of independent range requests; records in gaps
are read when they lie inside a chromosome span but are never written.

Do not begin with dozens of region workers reading the same objects through a
Cloud Storage FUSE mount. Start with about 8 workers and increase concurrency
only while aggregate throughput improves. When possible, stage each
chromosome's two VCFs and their indexes once onto local SSD. If the mount is
under your control, Cloud Storage FUSE's
[file cache](https://docs.cloud.google.com/storage/docs/cloud-storage-fuse/file-caching)
and `--file-cache-cache-file-for-range-read=true` are intended for repeated
partial/random reads; size the cache to avoid
[cache thrashing](https://docs.cloud.google.com/storage/docs/cloud-storage-fuse/caching).

Each region command is a separate process with its own sample map, ancestry
state, input buffers, and index handles. A 96-vCPU Workbench therefore must not
automatically be treated as 96 safe FELIXla workers. For Cloud Storage FUSE,
start with `MultiRun.sh cmds.sh 8`, then try 12 or 16 only if total throughput
increases and memory and swap remain stable. Shell messages such as `Killed`
with exit status 137 only establish that the process received `SIGKILL`; inspect
the cgroup events and the job runner logs to distinguish a memory limit from an
external cancellation. Partial prefixes from those commands must not be
concatenated.

`--export vcf` always writes BGZF with a tabix index beside it: a full-cohort
export is far too large to be worth keeping uncompressed. An `--out` ending in
`.vcf` is corrected to `.vcf.gz` with a note on stderr, and any other `--out`
gains the `.vcf.gz` suffix, so the two files written are always
`<out>.vcf.gz` and `<out>.vcf.gz.tbi`.

The same binary also dispatches to compatibility subcommands:

```
felixla --phase-vcf PHASED_VCF --rfmix-msp MSP_FILE --n-ancestries N --make-felixla --out OUT_PREFIX
felixla --tractor-dosage-vcf DOSAGE_VCF --n-ancestries N --make-felixla --out OUT_PREFIX
felixla --felixla PREFIX --export vcf --out OUTPUT_VCF_GZ
felixla --felixla PREFIX --query CHR:POS --ref REF --alt ALT [ --nonzero-only ]
felixla --felixla PREFIX --region CHR:START-END --make-felixla --out OUT_PREFIX
felixla --vcf TRACTOR_VCF --admixture --out ADMIXTURE_TSV
felixla --compare-vcfs EXPECTED_VCF OBSERVED_VCF [ --split-multiallelic ]
felixla --pmerge-list PREFIX_LIST --make-felixla --out OUT_PREFIX
felixla --recommend-mac-threshold --n-samples N_SAMPLES
```

For compatibility with older command lines, the original command names are
accepted as `felixla` subcommands, for example
`felixla flare_subset_to_tractor_hybrid ...` and
`felixla tractor_hybrid_to_vcf ...`. They are not installed as separate
executables.

## Output

`<prefix>.log` records the run and `<prefix>.meta` describes the data, and the
split is deliberate: anything that would differ between two runs producing the
same file belongs in the log, and anything a reader needs in order to
interpret the file belongs in the metadata.

The log holds the FELIXla version that ran, the command as invoked, the input
paths, every filter applied, the resulting counts, and a closing `completed`
line. A fatal error is appended there too, so a failed run leaves its own
explanation, and a log not ending in `completed` marks a prefix that did not
finish and must not be concatenated.

The metadata holds `format_version`, `n_samples`, `n_haps`, `n_words`,
`n_ancestries`, `rare_threshold`, the four variant and ancestry-block counts,
and `selected_region`. Every one of those is read by something:
`format_version`, `n_samples`, `n_haps`, `n_words` and `n_ancestries` are
required by FELIXassoc, and concat additionally reads `rare_threshold`, the
counts, and `selected_region`, which tells it the coordinate span the prefix
covers. Readers ignore keys they do not know, so the metadata can gain fields
without breaking them.


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
  count, word count, rare threshold, final common/rare/block counts, and
  optional region/filter/concatenation provenance.

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

Select many BED intervals at conversion time:

```bash
bin/felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --n-ancestries 5 \
  --extract-bed targets.bed \
  --make-felixla \
  --out hybrid/targets
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

Run the extended keep/extract/BED and concatenation regression tests:

```
make test-intense
```

`test-intense` synthesizes multi-chromosome, multi-allelic phased genotype and
FLARE inputs, compares the PLINK-style command against the compatibility
subcommand path, verifies `--keep`, `--extract`, `--extract-bed`, indexed
chromosome-span scanning, BED boundary and ancestry-gap behavior, all-filter
intersections, `--region`, concatenation, and `felixla --query` against known
truth, and checks every concat component plus malformed filter inputs for
fail-closed behavior. The BED suite also verifies 20,001 disjoint intervals and
byte-identical non-metadata output when BED covers the complete input.

`test-intense` also packs randomized fixtures twice and requires byte-identical
output both ways: once normally, and once with `FELIXLA_SCALAR_PATHS=1`, which
forces the column-at-a-time reader instead of the wide scanning paths the FLARE
reader and genotype packer select at run time. The fixtures deliberately
straddle the conditions those paths test, including ragged versus uniform FLARE
column widths, two-digit ancestry labels, the three `FORMAT` layouts, sample
counts either side of a 32-column word, sample subsets, and multiallelic
records. Setting `FELIXLA_SCALAR_PATHS=1` is also the way to confirm that a
suspected packing difference comes from the wide paths rather than the input.

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

A fully static build needs a static archive for everything htslib reports in
`pkg-config --libs --static htslib`, so the simplest htslib to build against is
one configured for exactly what FELIXla reads. `bz2` and `lzma` are only used
for CRAM, and remote-URL support is not needed to read local VCF/BCF:

```
./configure --prefix=$PREFIX --with-libdeflate \
    --disable-libcurl --disable-s3 --disable-gcs --disable-plugins \
    --disable-bz2 --disable-lzma
make libhts.a
make install-pkgconfig prefix=$PREFIX
install -d $PREFIX/lib $PREFIX/include/htslib
install -m644 libhts.a $PREFIX/lib/
install -m644 htslib/*.h $PREFIX/include/htslib/
```

Build the `libhts.a` target rather than the default one: `libhts.so` does not
link against a static `libdeflate.a` or `libz.a` unless those were compiled as
position-independent code, and FELIXla needs only the archive. Then point
FELIXla at it with `PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig make static
STATIC_FULLY=1`. The result is about 4 MB stripped, runs with no
`LD_LIBRARY_PATH`, and keeps both libdeflate decompression and `--threads`.

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
