#!/usr/bin/env python3
"""Pack the same inputs with and without the wide scanning paths.

The FLARE reader and the genotype packer each carry a wide path that is chosen
at run time: a 32-byte column scan, a fixed-stride read of uniform-width FLARE
columns, and a compare-and-gather that packs 32 genotype columns at once. Each
one must produce exactly what the column-at-a-time code produces, so this packs
every fixture twice -- once normally, once with FELIXLA_SCALAR_PATHS=1 -- and
requires the outputs to be byte-identical.

The fixtures deliberately straddle the conditions those paths test: ragged
versus uniform FLARE column widths, two-digit ancestry labels, the three
FORMAT layouts, sample counts either side of a 32-column word, sample subsets,
and multiallelic records.
"""

import argparse
import filecmp
import os
import pathlib
import random
import subprocess
import sys
import tempfile

COMPONENTS = (
    "common.geno.bin", "common.variant.mks", "common.variant.idx",
    "rare.carrier.bin", "rare.variant.mks", "rare.variant.idx",
    "ancblock.bin", "ancblock.mks", "ancblock.idx", "samples",
)

# "anc_only" is the slimmed FLARE layout, whose three-byte columns are the ones
# the ancestry reader decodes eight at a time.
LAYOUTS = ("gt_first", "anc_first", "anc_late", "anc_only")


def write_fixture(directory, rng, n_samples, n_records, n_ancestries, layout,
                  ragged, multiallelic_every, drop_flare_every):
    sample_ids = ["s%05d" % i for i in range(n_samples)]
    genotype = directory / "genotype.vcf"
    flare = directory / "flare.vcf"

    with genotype.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(sample_ids) + "\n")
        for record in range(n_records):
            position = 1000 + record * 7
            if multiallelic_every and record % multiallelic_every == 0:
                calls = ["%d|%d" % (rng.randrange(3), rng.randrange(3))
                         for _ in range(n_samples)]
                alts = "C,G"
            else:
                frequency = rng.choice([0.5, 0.2, 0.02, 0.001])
                calls = ["%d|%d" % (rng.random() < frequency,
                                    rng.random() < frequency)
                         for _ in range(n_samples)]
                alts = "C"
            out.write("chr1\t%d\tv%d\tA\t%s\t.\tPASS\t.\tGT\t%s\n"
                      % (position, record, alts, "\t".join(calls)))

    flare_samples = [i for i in range(n_samples)
                     if not drop_flare_every or (i + 1) % drop_flare_every != 0]
    ancestries = [(rng.randrange(n_ancestries), rng.randrange(n_ancestries))
                  for _ in flare_samples]
    format_tags = {
        "gt_first": "GT:AN1:AN2:ANP1:ANP2",
        "anc_first": "AN1:AN2:ANP1:ANP2",
        "anc_late": "GT:ANP1:AN1:AN2:ANP2",
        "anc_only": "AN1:AN2",
    }[layout]

    with flare.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="a1">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="a2">\n')
        out.write('##FORMAT=<ID=ANP1,Number=1,Type=Float,Description="p1">\n')
        out.write('##FORMAT=<ID=ANP2,Number=1,Type=Float,Description="p2">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                  + "\t".join(sample_ids[i] for i in flare_samples) + "\n")
        for record in range(n_records):
            if record and rng.random() < 0.3:
                for _ in range(rng.randrange(1, 4)):
                    j = rng.randrange(len(ancestries))
                    ancestries[j] = (rng.randrange(n_ancestries),
                                     rng.randrange(n_ancestries))
            columns = []
            for a1, a2 in ancestries:
                if ragged:
                    p1 = rng.choice(["1", "0.5", "0.99", "0.123"])
                    p2 = rng.choice(["0", "0.75", "0.9999"])
                else:
                    p1, p2 = "0.99", "0.98"
                call = "%d|%d" % (rng.randrange(2), rng.randrange(2))
                if layout == "gt_first":
                    columns.append("%s:%d:%d:%s:%s" % (call, a1, a2, p1, p2))
                elif layout == "anc_first":
                    columns.append("%d:%d:%s:%s" % (a1, a2, p1, p2))
                elif layout == "anc_only":
                    columns.append("%d:%d" % (a1, a2))
                else:
                    columns.append("%s:%s:%d:%d:%s" % (call, p1, a1, a2, p2))
            out.write("chr1\t%d\t.\tA\tC\t.\tPASS\t.\t%s\t%s\n"
                      % (1000 + record * 7, format_tags, "\t".join(columns)))

    return genotype, flare


def pack(felixla, genotype, flare, n_ancestries, prefix, keep, scalar_paths):
    command = [str(felixla),
               "--phase-vcf", str(genotype),
               "--flare-vcf", str(flare),
               "--n-ancestries", str(n_ancestries),
               "--export-felixla",
               "--out", str(prefix)]
    if keep:
        command += ["--keep", str(keep)]
    environment = dict(os.environ)
    if scalar_paths:
        environment["FELIXLA_SCALAR_PATHS"] = "1"
    else:
        environment.pop("FELIXLA_SCALAR_PATHS", None)
        command += ["--threads", "2"]
    return subprocess.run(command, env=environment, universal_newlines=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=pathlib.Path,
                        default=pathlib.Path("bin"))
    parser.add_argument("--cases", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()

    felixla = (args.bin_dir / "felixla").resolve()
    if not felixla.exists():
        print("missing %s" % felixla, file=sys.stderr)
        return 1

    rng = random.Random(args.seed)
    failures = 0
    for case in range(args.cases):
        n_samples = rng.choice([1, 2, 7, 31, 32, 33, 64, 65, 129])
        n_records = rng.choice([1, 3, 12])
        n_ancestries = rng.choice([2, 3, 5, 11, 17])
        layout = rng.choice(LAYOUTS)
        ragged = rng.random() < 0.5
        multiallelic_every = rng.choice([0, 0, 3])
        drop_flare_every = rng.choice([0, 0, 5])

        with tempfile.TemporaryDirectory(prefix="felixla-scalar-equiv.") as tmp:
            directory = pathlib.Path(tmp)
            genotype, flare = write_fixture(
                directory, rng, n_samples, n_records, n_ancestries, layout,
                ragged, multiallelic_every, drop_flare_every)

            keep = None
            if n_samples >= 4 and rng.random() < 0.4:
                kept = ["s%05d" % i for i in range(n_samples)
                        if rng.random() < 0.7]
                if len(kept) >= 2:
                    keep = directory / "keep.txt"
                    keep.write_text("\n".join(kept) + "\n")

            described = ("n_samples=%d n_records=%d n_ancestries=%d %s "
                         "ragged=%d multiallelic_every=%d drop_flare_every=%d"
                         % (n_samples, n_records, n_ancestries, layout,
                            ragged, multiallelic_every, drop_flare_every))

            wide = pack(felixla, genotype, flare, n_ancestries,
                        directory / "wide", keep, scalar_paths=False)
            scalar = pack(felixla, genotype, flare, n_ancestries,
                          directory / "scalar", keep, scalar_paths=True)

            if (wide.returncode == 0) != (scalar.returncode == 0):
                print("FAIL exit status differs: %s\n  wide=%d %s\n  scalar=%d %s"
                      % (described, wide.returncode, wide.stderr.strip()[-200:],
                         scalar.returncode, scalar.stderr.strip()[-200:]))
                failures += 1
                continue
            if wide.returncode != 0:
                continue

            for component in COMPONENTS:
                a = directory / ("wide." + component)
                b = directory / ("scalar." + component)
                if not a.exists() and not b.exists():
                    continue
                if not (a.exists() and b.exists()) or not filecmp.cmp(
                        a, b, shallow=False):
                    print("FAIL %s differs: %s" % (component, described))
                    failures += 1
                    break

    if failures:
        print("scalar equivalence FAILED in %d of %d cases"
              % (failures, args.cases))
        return 1
    print("scalar equivalence passed over %d cases" % args.cases)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
