# FELIXla

**FELIX** (**F**ull-cohort **E**fficient **L**ocal ancestry-**I**ntegrated
mi**X**ed-model framework) is a framework for scalable local-ancestry-aware
association analysis. **FELIXla** is its C++ based command line tool and binary
storage format for haplotype-resolved, local-ancestry-aware genotype data. It
stores phased ALT alleles and inferred local ancestries under the same explicit
haplotype index, so ancestry-specific dosages can be read straight out of the
packed files instead of being recomputed from phased genotype and local
ancestry files for every phenotype.

The ancestry-specific dosage for ancestry `k` at a variant is represented as a
bitwise intersection between the ALT haplotype vector and the ancestry `k`
haplotype mask. Local ancestry is piecewise constant along a chromosome, so
FELIXla stores ancestry transition blocks rather than one ancestry label per
variant. Common variants are stored as dense bit vectors, and rare variants are
stored as sparse carrier lists.

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

Every run names **one input**, **one output**, and optionally some
**parameters**:

```
felixla <input> <output> [parameters]
```

```text
Input       --phase-vcf --flare-vcf      phased genotypes plus FLARE ancestry
            --phase-vcf --rfmix-msp      phased genotypes plus an RFMix MSP file
            --tractor-dosage-vcf         a TRACTOR dosage VCF
            --felixla                    an existing packed prefix
            --merge-list                 a file of prefixes to concatenate

Output      --out                        the output path, always required
            --export-felixla             a packed prefix
            --export-vcf                 split-biallelic phased VCF.gz plus .tbi
            --export-lai                 local ancestry intervals, <out>.lai.gz
            --export-global-admixture    per-sample ancestry proportions
            --export-local-admixture     per-region ancestry proportions

Parameters  --keep --extract --exclude --extract-bed --exclude-bed
            --chr --region --mac --maf --anc-mac --anc-maf
            --n-ancestries --threads --version --help
```

A typical conversion:

```bash
felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --export-felixla \
  --out hybrid/chr22
```

Exactly one input and exactly one output are required, and combinations that
cannot mean anything are refused by name rather than silently ignored: an
export that reads a packed prefix says so, and a removed flag names what
replaced it.

Most filters apply while packing, which is the only point at which FELIXla
reads the source records, so they are accepted with `--phase-vcf` plus
`--flare-vcf` and refused elsewhere. `--region`, `--extract` and `--exclude`
are the exceptions: they also select out of a packed prefix, so

```bash
felixla --felixla hybrid/chr22 --extract sites.pvar \
        --export-felixla --out hybrid/chr22.subset
```

writes the same prefix that packing with `--extract sites.pvar` would have
written, without re-reading the source VCFs.

### Input

- **`--phase-vcf PATH`** -- a phased diploid genotype VCF/BCF, given together
  with `--flare-vcf` or `--rfmix-msp`. Genotypes must be phased, diploid, and
  non-missing. The sex chromosomes are not supported: the packed format gives
  every sample exactly two haplotypes, so haploid genotypes cannot be
  represented, and a record on `X` or `Y` is refused with that explanation.
  Select the autosomes explicitly, for example `--chr 1-22`. Multi-allelic
  records are split logically by ALT allele.

- **`--flare-vcf PATH`** -- a FLARE local ancestry VCF/BCF. It must carry
  scalar integer `FORMAT/AN1` and `FORMAT/AN2` fields encoded as
  `0..n_ancestries-1`. Sample IDs are matched by ID; FELIXla keeps the
  genotype/FLARE intersection in genotype VCF order.

- **`--rfmix-msp PATH`** -- an RFMix MSP file, given together with
  `--phase-vcf` in place of `--flare-vcf`. The ancestry count and the ancestry
  names both come from its `#Subpopulation order/codes:` header line. `--region`
  is not supported for this input; convert the whole file and extract
  afterwards.

- **`--tractor-dosage-vcf PATH`** -- a TRACTOR dosage VCF/BCF carrying
  `DS1`..`DSk` and `ANC1`..`ANCk` `FORMAT` fields. It is a complete input on
  its own: the
  dosages already carry both the genotypes and the ancestry assignment, so no
  `--phase-vcf` is given. The ancestry count comes from the `FORMAT`
  declarations in the header, which do not name the ancestries.

- **`--felixla PREFIX`** -- an existing FELIXla prefix, named without a
  component-file suffix. It is the input to every export, and, with
  `--export-felixla` plus `--region`, `--extract` or `--exclude`, to a new
  prefix holding part of it.

- **`--merge-list FILE`** -- a file of prefixes to concatenate in genomic
  order, one per line, or two tab-separated columns to take a different BED
  from each. See [Concatenating chunks](#concatenating-chunks).

#### What FELIXla assumes about its inputs

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

### Output

- **`--out PATH`** -- where the output goes. For `--export-felixla` it is the
  prefix the eleven component files are named from; for every other export it
  is the stem the suffix is added to.

- **`--export-felixla`** -- write a packed prefix. This is the conversion that
  the filters apply to.

  The sparse/dense storage threshold is not a command-line option. Variants
  whose ALT count is at or below `ceil(n_samples / 32)` are stored as sparse
  carrier lists and the rest as dense bit vectors, which is the crossover at
  which a carrier list stops being smaller than a bit vector.

- **`--export-vcf`** -- reconstruct a split-biallelic phased genotype VCF from
  a packed prefix.

  `--export-vcf` always writes BGZF with a tabix index beside it: a full-cohort
  export is far too large to be worth keeping uncompressed. An `--out` ending in
  `.vcf` is corrected to `.vcf.gz` with a note on stderr, and any other `--out`
  gains the `.vcf.gz` suffix, so the two files written are always
  `<out>.vcf.gz` and `<out>.vcf.gz.tbi`.

- **`--export-lai`**, **`--export-global-admixture`**,
  **`--export-local-admixture`** -- the three read-only views of a prefix's
  ancestry blocks, described under [Ancestry exports](#ancestry-exports).

#### The FELIXla prefix

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
`selected_region`, and one `ancestry_name_N` row per ancestry the input named,
numbered from one to match the ancestry exports. Every one of those is read by
something:
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

- `<prefix>.meta`: format metadata, sample count, ancestry count and names,
  word count, rare threshold, final common/rare/block counts, and the selected
  region where one was given.

The split variant marker streams record `chr`, `pos`, split-biallelic `id`,
`ref`, `alt`, ALT allele index, global variant index, and MAC. Marker streams
and indexes are little-endian binary files with 8-byte magic headers.

#### Ancestry exports

Three read-only views summarize a packed prefix's ancestry blocks. All three
take `--felixla PREFIX` and add their own suffix to `--out`, so `--out run`,
`--out run.lai` and `--out run.lai.gz` all write `run.lai.gz`.

**`--export-lai`** writes `<out>.lai.gz`, a BGZF-compressed table with one row
per ancestry block and one column per haplotype:

```text
##ANC1 = AFR
##ANC2 = EUR
#CHR	START	END	s0001_1	s0001_2	s0002_1	s0002_2
chr1	1	48210	1	2	1	1
chr1	48211	93004	1	2	2	1
```

The `##` lines name the ancestries, taking the names from the input header
where it gave them and falling back to `ANC1`, `ANC2`, ... where it did not.
Each sample contributes two columns, `_1` and `_2`, in the order the
haplotypes are packed. `START` and `END` are 1-based and inclusive, and the
blocks are sorted, non-overlapping, and cover every variant in the prefix.

**Ancestry codes in this file are 1-based**, so `1` is the ancestry the `##`
header calls `ANC1`. This differs from the codes inside the packed files and
from FLARE's own `AN1`/`AN2`, which both number from zero. The offset is
deliberate: `0` is reserved for a haplotype no caller labelled, so a future
format that admits missing ancestry can use it without any existing column
changing meaning.

**`--export-global-admixture`** writes `<out>.global.admixture.tsv`, one row
per sample:

```text
#ID	AFR	EUR
s0001	0.312500	0.687500
s0002	0.687500	0.312500
```

Proportions are weighted by base pairs rather than by blocks, so a region a
caller happened to split finely is not thereby counted more heavily, and each
row sums to one. The denominator is the span the prefix's ancestry blocks
actually cover, which for a region or contig subset is that subset, not the
genome: proportions from two prefixes are comparable only over the same span.

**`--export-local-admixture`** writes `<out>.local.admixture.tsv`, one row per
ancestry block, holding each ancestry's share of the cohort's haplotypes
there:

```text
#CHR	START	END	AFR	EUR
chr1	1	48210	0.750000	0.250000
chr1	48211	93004	0.500000	0.500000
```

The rows use the same blocks and the same coordinates as the LAI file, so the
two line up row for row.

### Parameters

#### Sample and variant filters

Most of these apply while packing, so they go with `--phase-vcf` plus
`--flare-vcf` and `--export-felixla`. `--region`, `--extract` and `--exclude`
also apply to `--felixla --export-felixla`, selecting out of a prefix that is
already packed.

- **`--keep FILE`** -- one sample ID per line, analogous to
  PLINK `--keep`. Retained samples are written in genotype VCF order.

- **`--extract FILE`** -- a PLINK2 `.pvar`-like file or VCF-like file used like
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

  `--extract` also selects out of an existing prefix, with
  `--felixla PREFIX --extract FILE --export-felixla`. The matching rules are
  the same ones, applied to the prefix's own split-biallelic markers instead of
  to VCF records, and the result is byte-identical to packing the source with
  the same list -- which the test suite asserts, component by component, so
  the two cannot drift apart. A listed coordinate the prefix holds under a
  different REF is a fatal mismatch there too.

- **`--exclude FILE`** -- the same PVAR/VCF format and the same
  allele matching as `--extract`, naming alleles to drop rather than keep. It
  is applied after selection, so `--extract` and `--exclude` compose as they do
  in PLINK. Matching is allele-level: a multiallelic record keeps whatever ALTs
  the list does not name, and the record is dropped only when nothing remains.
  An entry matching no record in the input is silently inert. Unlike
  `--extract`, the list is not narrowed by `--region`, since a site outside the
  region is absent from the output regardless. A REF that matches nothing is
  inert rather than fatal: an `--exclude` entry naming no record simply removes
  nothing.

  Like `--extract`, it also applies to `--felixla --export-felixla`, and the
  two compose there in the same order: `--extract` selects, then `--exclude`
  removes.

  Note that neither `--extract` nor `--exclude` uses the `ID` column: matching
  is on `CHROM`, `POS` and the normalized `REF`/`ALT`. A PLINK-style file of
  bare variant IDs is therefore rejected rather than silently misread.

- **`--extract-bed FILE`** -- a BED or gzip-compressed BED file. FELIXla reads
  the first three columns and ignores later columns, `track`/`browser` rows,
  blank lines, and comments. Coordinates follow the BED standard: 0-based,
  half-open `[START, END)`. Variant selection uses the VCF record position, so
  `chr1 99 100` selects a record at VCF position 100. Input intervals may be
  unsorted, duplicated, overlapping, or adjacent; FELIXla sorts and merges
  them. Every ALT of a selected multiallelic genotype record is retained unless
  `--extract` narrows the allele set. `--extract-bed`, `--extract`, and
  `--region` are intersected when combined. Ancestry blocks are clipped and
  split at BED gaps.

- **`--exclude-bed FILE`** -- BED intervals to drop, in the same
  0-based half-open convention as `--extract-bed`. Exclusion narrows the
  selection rather than filtering records separately, so the ancestry blocks
  are clipped at exclusion boundaries by the same rule that clips them at
  `--extract-bed` gaps, and a block never covers an excluded base. With no
  `--extract-bed`, the selection starts as the whole of every genotype contig,
  or as `--region` when one is given. Like `--exclude`, the list is not
  narrowed by `--region`.

  With an indexed genotype VCF/BCF, FELIXla creates one bounded read span from
  the first through last retained BED interval on each chromosome, then applies
  the exact merged intervals with a linear cursor. It therefore does not issue
  one tabix request per BED row. This is intentional for large interval lists
  and remote object-store mounts. Without a usable index, it streams the input
  once and applies the same exact interval filter.

- **`--chr LIST`** -- select whole contigs to convert.
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

- **`--region CHR:START-END`** -- 1-based inclusive region to convert. A variant
  belongs to the region when `START <= POS <= END`.

- **`--mac INT`**, **`--maf FLOAT`**, **`--anc-mac INT`**, **`--anc-maf F`**
  -- drop variants with
  too little variation to be worth carrying. All four are stated on the
  **minor** allele count or frequency, `min(alt, total - alt)`, so an allele
  carried by every haplotype is filtered as readily as one carried by none.
  Note this differs from the sparse/dense storage threshold, which is a
  layout decision made on the raw ALT count, not a filter: a variant below it
  is stored differently, not dropped.

  `--mac` and `--maf` count across all ancestries, as PLINK's do. `--anc-mac`
  and `--anc-maf` apply the threshold to each ancestry separately and keep the
  variant when **at least one** ancestry meets it on its own, each ancestry's
  frequency being taken over the haplotypes assigned to it. The two answer
  different questions: an allele spread thinly over several ancestries can
  clear `--mac` while no single ancestry carries enough of it for that
  ancestry's own test to be informative, and `--anc-mac` is what removes it.
  Combining them applies both. Dropped alleles consume no variant ordinal, and
  the count removed is reported in the log.

#### Other

- **`--n-ancestries INT`** -- number of local ancestry labels. The packed format
  stores ancestry codes in 5 bits, so at most 32 labels are supported. It is
  normally omitted: FELIXla takes the count from the header of whichever input
  defines it -- FLARE's `##ANCESTRY` lines, RFMix's `#Subpopulation
  order/codes:` line, or the `DS#`/`ANC#` `FORMAT` declarations of a TRACTOR
  dosage VCF -- and refuses the run if a value given here disagrees. Taking it
  from the data instead would be wrong: a region that happens to carry no
  haplotype of some ancestry would come out with a smaller count than its
  neighbours, and the two prefixes could then not be concatenated. When an
  input names its ancestries, the names are kept in `<prefix>.meta` and label
  the columns of the ancestry exports.

- **`--threads INT`** -- number of BGZF decompression threads shared by the
  genotype and FLARE readers. Default is `1`, which keeps one conversion in one
  core. Decompression is the largest single cost when both inputs are
  BGZF-compressed, so raising this is the most effective way to speed up a
  single conversion; it does not change the output. Parsing is serial and
  becomes the limit at roughly 8 threads. Threads and region workers compete
  for the same cores: use `--threads` when you run few conversions at a time,
  and leave it at `1` when you already saturate the machine with concurrent
  region jobs.

  Only two paths have anything to hand threads to: the FLARE conversion, where
  BGZF decompression is most of the work, and `--export-lai`, where BGZF
  compression is. The RFMix and dosage conversions, the region extract, the
  merge, and the other exports are single-threaded, and give `--threads` back
  as an error rather than accepting a tuning knob that would do nothing.

- **`--version`**, **`--help`** -- print the version or the flag summary and
  exit.

#### Decompression performance

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

## Region-parallel Workflows

### Planning chunks

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
    'felixla --phase-vcf {phase_vcf_q} --flare-vcf flare.vcf.gz --region {region_q} --export-felixla --out out/{chrom}.chunk{chrom_chunk0}' \
  --commands-out genotype.commands.txt
```

Available template fields include `{phase_vcf}`, `{chrom}`, `{start}`, `{end}`,
`{region}`, `{global_chunk}`, `{global_chunk0}`, `{chrom_chunk}`, and
`{chrom_chunk0}`. `{phase_vcf_q}`, `{chrom_q}`, and `{region_q}` are
POSIX-shell-quoted forms.

### Concatenating chunks

After region jobs finish, concatenate complete FELIXla prefixes in genomic
order with the PLINK-style `--merge-list` interface:

```bash
felixla \
  --merge-list chr1.prefixes.txt \
  --export-felixla \
  --out chr1.merged
```

The list contains one prefix per line, without a component-file suffix:

```text
out/chr1.chunk0001
out/chr1.chunk0002
out/chr1.chunk0003
```

Blank lines and lines beginning with `#` are ignored. A trailing `.meta` is
also accepted. Paths are interpreted relative to the working directory.

Ancestry names carry through the merge, and two chunks that name the same code
differently are refused before any merging starts: a prefix whose code `3`
means AFR in one half and EUR in the other would be silently wrong in every
export that reads it. A chunk that names none of its ancestries merges freely,
since its codes still line up.

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

## Example

The following commands build the binary, create a FELIXla prefix from the tiny
phased genotype and FLARE local ancestry fixtures, and export the packed data
back to split-biallelic VCF and as a local ancestry table:

```
make

bin/felixla \
  --phase-vcf testdata/tiny.genotypes.vcf \
  --flare-vcf testdata/tiny.flare.vcf \
  --export-felixla \
  --out example/tiny

bin/felixla \
  --felixla example/tiny \
  --export-vcf \
  --out example/tiny.roundtrip

bin/felixla \
  --felixla example/tiny \
  --export-lai \
  --out example/tiny
```

Sample and site filtering can be applied at conversion time:

```
bin/felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --keep samples.keep \
  --extract sites.pvar \
  --export-felixla \
  --out hybrid/chr22.subset
```

Select many BED intervals at conversion time:

```bash
bin/felixla \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --extract-bed targets.bed \
  --export-felixla \
  --out hybrid/targets
```

The same command can be run through Docker by binding the working directory:

```
docker run --rm -v "$PWD":/data -w /data ghcr.io/yorkklause/felixla:latest \
  --phase-vcf genotype.phased.vcf.gz \
  --flare-vcf flare.anc.vcf.gz \
  --keep samples.keep \
  --extract sites.pvar \
  --export-felixla \
  --out hybrid/chr22.subset
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
intersections, `--region`, concatenation, and packed dosages against known
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

`test-intense` also checks that selecting alleles out of a packed prefix lands
exactly where selecting them while packing lands: for each random case it packs
with `--extract`/`--exclude`, packs the whole input and then filters the prefix
with the same list, and requires the two results to be byte-identical component
by component. Two implementations of one filter can only stay in step if
something insists on it.

`test-intense` finally checks the three ancestry exports against three separate
oracles: the LAI rows against the ancestry masks decoded independently from
`.ancblock.bin`, the LAI rows at every genotype position against the FLARE
input that produced them, and both admixture tables against a fresh exact
aggregation of the LAI file. Nothing in that chain asks FELIXla to grade its
own output.

The readers the suites check against live in `tests/felixla_check.py`, which
parses the packed files, the LAI file and the admixture tables in Python with
no help from the binary. Keeping the oracle a separate implementation is the
point: a reader bug can no longer cancel out a writer bug.

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

## Support

Please direct questions or bug reports to Kai Yuan (kyuan@broadinstitute.org).
