#!/usr/bin/env python3
"""Selecting alleles out of a packed prefix must equal selecting while packing.

The criterion is byte equality of every payload component between

    felixla --phase-vcf G --flare-vcf F --extract S --export-felixla --out A
    felixla --phase-vcf G --flare-vcf F --export-felixla --out P
    felixla --felixla P --extract S --export-felixla --out B

so the two implementations of the same filter cannot drift: A and B must be
indistinguishable. The same is asserted for --exclude and for the two composed.

Adding --region to both routes compares only the variant components. The two
already bound the last ancestry block differently -- packing ends it at the
last genotype position inside the region, extracting clips it to the region
end -- and that predates the prefix-side site lists; it reproduces on v0.6.0
with no site list involved. Asserting either value here would freeze a
difference this suite is not about.
"""
import argparse, pathlib, random, subprocess, sys, tempfile

VARIANT_COMPONENTS = ["common.geno.bin", "common.variant.mks", "common.variant.idx",
                      "rare.carrier.bin", "rare.variant.mks", "rare.variant.idx",
                      "samples"]
ANCESTRY_COMPONENTS = ["ancblock.bin", "ancblock.mks", "ancblock.idx"]
BASES = "ACGT"


def build_inputs(d, rng, n_samples, n_records, n_anc):
    samples = ["s%03d" % i for i in range(n_samples)]
    variants = []
    pos = 100
    for r in range(n_records):
        pos += rng.randrange(3, 40)
        ref = rng.choice(BASES)
        if rng.random() < 0.25:
            # Padded representations, so the list and the VCF can disagree
            # about how many shared bases they carry.
            ref = ref + rng.choice(BASES)
            alts = [ref[0] + rng.choice(BASES)]
        else:
            alts = [b for b in rng.sample(BASES, 3) if b != ref][:rng.choice([1, 1, 2])]
        variants.append(("chr1" if r % 3 else "chr2", pos, ref, alts))
    variants.sort(key=lambda v: (v[0], v[1]))

    with (d / "g.vcf").open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
        for chrom, pos, ref, alts in variants:
            n = len(alts) + 1
            gts = "\t".join("%d|%d" % (rng.randrange(n), rng.randrange(n))
                            for _ in samples)
            out.write("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t%s\n"
                      % (chrom, pos, ref, ",".join(alts), gts))

    state = [(rng.randrange(n_anc), rng.randrange(n_anc)) for _ in samples]
    with (d / "f.vcf").open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n")
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="a">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="b">\n')
        for k in range(n_anc):
            out.write("##ANCESTRY=<ID=%d,Name=ANC%d>\n" % (k, k + 1))
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(samples) + "\n")
        for chrom, pos, _, _ in variants:
            if rng.random() < 0.25:
                state[rng.randrange(len(state))] = (rng.randrange(n_anc), rng.randrange(n_anc))
            out.write("%s\t%d\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t%s\n"
                      % (chrom, pos, "\t".join("%d:%d" % s for s in state)))
    return variants


def write_site_list(path, rng, variants, fraction):
    """A PVAR naming a random subset, sometimes padded, sometimes unsorted."""
    rows = []
    for chrom, pos, ref, alts in variants:
        chosen = [a for a in alts if rng.random() < fraction]
        if not chosen:
            continue
        if rng.random() < 0.3:
            # Same alleles, one extra shared base on each side.
            pad = rng.choice(BASES)
            rows.append((chrom, pos, ref + pad, [a + pad for a in chosen]))
        else:
            rows.append((chrom, pos, ref, chosen))
    if rng.random() < 0.3:
        rows.append(("chr1", 10 ** 7 + rng.randrange(1000), "A", ["T"]))  # inert
    rng.shuffle(rows)
    with path.open("w") as out:
        out.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for chrom, pos, ref, alts in rows:
            if len(alts) > 1 and rng.random() < 0.5:
                out.write("%s\t%d\t.\t%s\t%s\n" % (chrom, pos, ref, ",".join(alts)))
            else:
                for alt in alts:
                    out.write("%s\t%d\t.\t%s\t%s\n" % (chrom, pos, ref, alt))
    return sum(len(alts) for _, _, _, alts in rows)


def run(felixla, args):
    return subprocess.run([str(felixla)] + [str(a) for a in args],
                          capture_output=True, text=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin-dir", type=pathlib.Path, default=pathlib.Path("bin"))
    ap.add_argument("--cases", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20240928)
    a = ap.parse_args()
    a.felixla = a.bin_dir / "felixla"

    rng = random.Random(a.seed)
    failures = 0
    checked = 0
    for case in range(a.cases):
        n_samples = rng.choice([2, 9, 40])
        n_records = rng.choice([6, 20, 45])
        n_anc = rng.choice([2, 3])
        use_extract = rng.random() < 0.8
        use_exclude = rng.random() < 0.5 or not use_extract
        region = None
        if rng.random() < 0.3:
            region = "chr1:%d-%d" % (rng.randrange(100, 400), rng.randrange(500, 2000))

        with tempfile.TemporaryDirectory(prefix="felixla-prefix-extract.") as td:
            d = pathlib.Path(td)
            variants = build_inputs(d, rng, n_samples, n_records, n_anc)

            filters = []
            if use_extract:
                write_site_list(d / "keep.pvar", rng, variants, rng.choice([0.3, 0.6, 1.0]))
                filters += ["--extract", d / "keep.pvar"]
            if use_exclude:
                write_site_list(d / "drop.pvar", rng, variants, rng.choice([0.2, 0.5]))
                filters += ["--exclude", d / "drop.pvar"]
            region_args = ["--region", region] if region else []

            tag = ("case=%d n=%d r=%d anc=%d extract=%d exclude=%d region=%s"
                   % (case, n_samples, n_records, n_anc, use_extract, use_exclude, region))

            direct = run(a.felixla, ["--phase-vcf", d / "g.vcf", "--flare-vcf", d / "f.vcf"]
                         + region_args + filters + ["--export-felixla", "--out", d / "A"])
            full = run(a.felixla, ["--phase-vcf", d / "g.vcf", "--flare-vcf", d / "f.vcf",
                                   "--export-felixla", "--out", d / "P"])
            if full.returncode != 0:
                print("FAIL packing the whole input %s\n  %s" % (tag, full.stderr.strip()[-200:]))
                failures += 1
                continue
            via = run(a.felixla, ["--felixla", d / "P"] + region_args + filters
                      + ["--export-felixla", "--out", d / "B"])

            if (direct.returncode == 0) != (via.returncode == 0):
                print("FAIL one route refused %s\n  direct rc=%d %s\n  prefix rc=%d %s"
                      % (tag, direct.returncode, direct.stderr.strip()[-160:],
                         via.returncode, via.stderr.strip()[-160:]))
                failures += 1
                continue
            if direct.returncode != 0:
                continue

            compared = VARIANT_COMPONENTS + ([] if region else ANCESTRY_COMPONENTS)
            differing = [c for c in compared
                         if (d / ("A." + c)).read_bytes() != (d / ("B." + c)).read_bytes()]
            checked += 1
            if differing:
                print("FAIL %s\n  differing: %s" % (tag, differing))
                failures += 1

    print("%d/%d cases matched (%d compared byte-for-byte)"
          % (a.cases - failures, a.cases, checked))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
