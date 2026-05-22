# Linux x86_64 Static Binaries

These binaries are prebuilt for Linux x86_64 systems such as All of Us
Workbench notebooks. They are fully static musl ELF executables and do not
require conda, htslib, or system development headers.

Use them from a fresh checkout:

```bash
cd tractor-hybrid-tools
export PATH="$PWD/prebuilt/linux-x86_64-static:$PATH"
estimate_mac_threshold testdata/tiny.genotypes.vcf
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

The binaries were linked against htslib 1.23.1 with remote URL/S3/GCS support
disabled to keep them dependency-free. Use local VCF/BCF files; if your input is
in cloud storage, localize it into the Workbench environment first.

Verify checksums:

```bash
cd prebuilt/linux-x86_64-static
sha256sum -c SHA256SUMS
```
