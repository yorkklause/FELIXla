# Linux x86_64 Static Binaries

These binaries are prebuilt for Linux x86_64 systems such as All of Us
Workbench notebooks. They are fully static Linux ELF executables and do not
require conda, htslib, or system development headers.

Use them from a fresh checkout:

```bash
cd tractor-hybrid-tools
export PATH="$PWD/prebuilt/linux-x86_64-static:$PATH"
flare_subset_to_tractor_hybrid genotype.phased.vcf.gz flare.anc.vcf.gz 5 512 out/chr1
```

Or call a binary directly:

```bash
./prebuilt/linux-x86_64-static/flare_subset_to_tractor_hybrid \
  genotype.phased.vcf.gz \
  flare.anc.vcf.gz \
  3 \
  512 \
  out/chr1
```

The FLARE converter also accepts an optional final `chr:start-end` region:

```bash
./prebuilt/linux-x86_64-static/flare_subset_to_tractor_hybrid \
  genotype.phased.vcf.gz \
  flare.anc.vcf.gz \
  3 \
  512 \
  out/chr1.chunk0001 \
  chr1:1-50000000
```

For storage-optimal sparse/dense packing, choose `mac_threshold` from the
sample count:

```bash
n_samples=100000
mac_threshold=$(( (n_samples + 31) / 32 ))
echo "$mac_threshold"
```

For chunked conversion, this repo includes SHAPEIT5-style GRCh38 4 cM chunk
files under `resources/shapeit5_chunks/b38_4cM/` and a helper that converts
those chunks to `flare_subset_to_tractor_hybrid` argument rows:

```bash
scripts/shapeit_chunks_to_tractor_args.sh \
  --chunks-dir resources/shapeit5_chunks/b38_4cM \
  --chrom-style chr \
  --phase-template 'phase/{chrom}.phased.vcf.gz' \
  --flare-template 'flare/{chrom}.flare.vcf.gz' \
  --n-ancestries 5 \
  --mac-threshold 512 \
  --out-prefix-template 'hybrid/{chrom}.shapeit4cM.chunk{chunk0}' \
  > flare_subset.shapeit4cm.args.tsv

xargs -a flare_subset.shapeit4cm.args.tsv -n 6 -P 8 \
  ./prebuilt/linux-x86_64-static/flare_subset_to_tractor_hybrid
```

The RFMix MSP converter is also included:

```bash
./prebuilt/linux-x86_64-static/rfmix_msp_to_tractor_hybrid \
  genotype.phased.vcf.gz \
  rfmix.msp.tsv.gz \
  5 \
  512 \
  out/chr22
```

Extract a region from an existing packed prefix:

```bash
./prebuilt/linux-x86_64-static/tractor_hybrid_extract_region \
  out/chr22 \
  chr22:16000000-17000000 \
  out/chr22.16_17mb
```

The binaries were linked against htslib 1.23 with remote URL/S3/GCS support
disabled to keep them dependency-free. Use local VCF/BCF files; if your input is
in cloud storage, localize it into the Workbench environment first.

Verify checksums:

```bash
cd prebuilt/linux-x86_64-static
sha256sum -c SHA256SUMS
```
