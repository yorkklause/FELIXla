#!/usr/bin/env python3
"""Intense regression tests for --export-lai and the admixture exports.

Three checks per case, each with its own oracle, so that a mistake shared by
the writer and one reader still shows up in the others:

  1. The LAI rows must match the ancestry masks decoded here from
     .ancblock.bin -- the exporter's own decode, done independently.
  2. At every genotype position, the LAI row covering it must give each
     haplotype the ancestry the FLARE input states there -- the semantics,
     end to end, never consulting the packed files.
  3. The two admixture tables must equal a fresh aggregation of the LAI file,
     computed in exact rational arithmetic.
"""

from __future__ import annotations

import argparse
import fractions
import pathlib
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import felixla_check

TOLERANCE = 5e-7


def write_inputs(work, rng, samples, positions, n_ancestries, names):
    """Write a genotype and a FLARE VCF, and return the truth the FLARE states.

    A FLARE record sits at every genotype position, so the ancestry in force at
    a position is that record's, with no carrying forward to reason about.
    """
    with (work / "genotypes.vcf").open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        for pos in positions:
            calls = "\t".join(
                f"{rng.randrange(2)}|{rng.randrange(2)}" for _ in samples
            )
            out.write(f"chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t{calls}\n")

    state = [(rng.randrange(n_ancestries), rng.randrange(n_ancestries)) for _ in samples]
    truth = {}
    with (work / "flare.vcf").open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n")
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second">\n')
        for code, name in enumerate(names):
            out.write(f"##ANCESTRY=<ID={code},Name={name}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        for pos in positions:
            if rng.random() < 0.3:
                index = rng.randrange(len(state))
                state[index] = (rng.randrange(n_ancestries), rng.randrange(n_ancestries))
            labels = "\t".join(f"{a}:{b}" for a, b in state)
            out.write(f"chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t{labels}\n")
            # One-based, in haplotype order, as the LAI file writes them.
            truth[pos] = [code + 1 for pair in state for code in pair]
    return truth


def decode_masks(masks, n_haps):
    """The 1-based ancestry of every haplotype, 0 where no mask claims it."""
    codes = [0] * n_haps
    for index, words in enumerate(masks):
        for hap in range(n_haps):
            if (words[hap >> 6] >> (hap & 63)) & 1:
                codes[hap] = index + 1
    return codes


def check_case(felixla, work, rng, problems):
    n_samples = rng.choice([1, 2, 7, 33, 70])
    n_ancestries = rng.choice([1, 2, 3, 5])
    names = [f"ANCESTRY{index + 1}" for index in range(n_ancestries)]
    samples = [f"s{index:03d}" for index in range(n_samples)]
    positions = sorted(rng.sample(range(100, 9000), rng.choice([5, 20, 60])))
    truth = write_inputs(work, rng, samples, positions, n_ancestries, names)

    prefix = work / "packed"
    subprocess.run(
        [str(felixla), "--phase-vcf", str(work / "genotypes.vcf"),
         "--flare-vcf", str(work / "flare.vcf"), "--export-felixla",
         "--out", str(prefix)],
        check=True, capture_output=True, text=True,
    )
    for flag in ("--export-lai", "--export-global-admixture", "--export-local-admixture"):
        subprocess.run(
            [str(felixla), "--felixla", str(prefix), flag, "--out", str(prefix)],
            check=True, capture_output=True, text=True,
        )

    meta = felixla_check.read_meta(prefix)
    n_haps = int(meta["n_haps"])
    lai_names, columns, rows = felixla_check.read_lai(str(prefix) + ".lai.gz")

    expected_names = [(f"ANC{index + 1}", names[index]) for index in range(n_ancestries)]
    if lai_names != expected_names:
        problems.append(f"LAI header {lai_names} != {expected_names}")
    expected_columns = [f"{s}_{h}" for s in samples for h in (1, 2)]
    if columns != expected_columns:
        problems.append(f"LAI columns {columns[:6]} != {expected_columns[:6]}")

    # 1. Against the masks, decoded here straight from the packed payload.
    blocks = felixla_check.read_ancestry_blocks(
        prefix, int(meta["n_ancestries"]), int(meta["n_words"])
    )
    if len(blocks) != len(rows):
        problems.append(f"{len(rows)} LAI rows for {len(blocks)} ancestry blocks")
    else:
        for (chrom, start, end, masks), row in zip(blocks, rows):
            expected = (chrom, start, end, decode_masks(masks, n_haps))
            if expected != row:
                problems.append(
                    f"block {chrom}:{start}-{end} decodes to {expected[3][:6]}, "
                    f"LAI says {row[3][:6]}"
                )
                break

    # 2. Against FLARE, at every position the genotype VCF carries.
    for pos in positions:
        covering = [row for row in rows if row[1] <= pos <= row[2]]
        if len(covering) != 1:
            problems.append(f"{len(covering)} LAI rows cover position {pos}")
            break
        if covering[0][3] != truth[pos]:
            problems.append(
                f"at {pos} the LAI says {covering[0][3][:6]}, "
                f"FLARE says {truth[pos][:6]}"
            )
            break

    # 3. The admixture tables, re-aggregated from the LAI file.
    per_sample = [[0] * n_ancestries for _ in range(n_samples)]
    expected_local = []
    covered_bp = 0
    for chrom, start, end, codes in rows:
        length = end - start + 1
        covered_bp += length
        counts = [0] * n_ancestries
        for hap, code in enumerate(codes):
            counts[code - 1] += 1
            per_sample[hap // 2][code - 1] += length
        expected_local.append(
            (chrom, start, end, [fractions.Fraction(c, n_haps) for c in counts])
        )

    header, table = felixla_check.read_tsv(str(prefix) + ".local.admixture.tsv")
    if header != ["#CHR", "START", "END"] + names:
        problems.append(f"local admixture header {header}")
    elif len(table) != len(expected_local):
        problems.append(f"local admixture has {len(table)} rows, want {len(expected_local)}")
    else:
        for got, (chrom, start, end, want) in zip(table, expected_local):
            if got[:3] != [chrom, str(start), str(end)]:
                problems.append(f"local admixture coordinates {got[:3]}")
                break
            if any(abs(float(g) - float(w)) > TOLERANCE for g, w in zip(got[3:], want)):
                problems.append(
                    f"local admixture at {chrom}:{start} is {got[3:]}, "
                    f"want {[float(w) for w in want]}"
                )
                break

    header, table = felixla_check.read_tsv(str(prefix) + ".global.admixture.tsv")
    if header != ["#ID"] + names:
        problems.append(f"global admixture header {header}")
    elif len(table) != n_samples:
        problems.append(f"global admixture has {len(table)} rows, want {n_samples}")
    else:
        for index, got in enumerate(table):
            if got[0] != samples[index]:
                problems.append(f"global admixture row {index} is {got[0]}")
                break
            # Both haplotypes are covered by every block, so the row sums to one.
            want = [fractions.Fraction(bp, covered_bp * 2) for bp in per_sample[index]]
            if any(abs(float(g) - float(w)) > TOLERANCE for g, w in zip(got[1:], want)):
                problems.append(
                    f"global admixture for {samples[index]} is {got[1:]}, "
                    f"want {[float(w) for w in want]}"
                )
                break
            if abs(sum(float(value) for value in got[1:]) - 1.0) > 1e-5:
                problems.append(f"global admixture for {samples[index]} does not sum to one")
                break


def check_merge_keeps_names(felixla, work, problems):
    """A merged prefix must keep the ancestry names, and refuse a clash.

    Chromosome-scale prefixes are normally built by merging region chunks, so
    names that survived packing but not the merge would be lost on exactly the
    path that matters.
    """
    rng = random.Random(0)
    samples = [f"s{index:03d}" for index in range(6)]
    names = ["AFR", "EUR", "AMR"]
    positions = list(range(100, 900, 50))
    write_inputs(work, rng, samples, positions, len(names), names)

    chunks = []
    for index, region in enumerate(("chr1:1-400", "chr1:401-1000")):
        prefix = work / f"chunk{index}"
        subprocess.run(
            [str(felixla), "--phase-vcf", str(work / "genotypes.vcf"),
             "--flare-vcf", str(work / "flare.vcf"), "--region", region,
             "--export-felixla", "--out", str(prefix)],
            check=True, capture_output=True, text=True,
        )
        chunks.append(prefix)

    listing = work / "prefixes.txt"
    listing.write_text("".join(f"{prefix}\n" for prefix in chunks))
    merged = work / "merged"
    subprocess.run(
        [str(felixla), "--merge-list", str(listing), "--export-felixla",
         "--out", str(merged)],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        [str(felixla), "--felixla", str(merged), "--export-lai", "--out", str(merged)],
        check=True, capture_output=True, text=True,
    )

    meta = felixla_check.read_meta(merged)
    for index, name in enumerate(names):
        if meta.get(f"ancestry_name_{index + 1}") != name:
            problems.append(f"merged prefix lost ancestry_name_{index + 1}")
    lai_names, _, rows = felixla_check.read_lai(str(merged) + ".lai.gz")
    if lai_names != [(f"ANC{i + 1}", names[i]) for i in range(len(names))]:
        problems.append(f"merged LAI header {lai_names}")
    if not rows:
        problems.append("merged LAI has no rows")

    # Two chunks that disagree about what a code means must not merge at all.
    clashing = pathlib.Path(str(chunks[1]) + ".meta")
    clashing.write_text(
        clashing.read_text().replace("ancestry_name_2\tEUR", "ancestry_name_2\tSAS")
    )
    clash = subprocess.run(
        [str(felixla), "--merge-list", str(listing), "--export-felixla",
         "--out", str(work / "clash")],
        capture_output=True, text=True,
    )
    if clash.returncode == 0:
        problems.append("a merge of prefixes naming ancestry 2 differently was accepted")
    elif "named EUR" not in clash.stderr or "SAS" not in clash.stderr:
        problems.append(f"the clash was refused without naming it: {clash.stderr.strip()[-200:]}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=pathlib.Path, default=pathlib.Path("bin"))
    parser.add_argument("--cases", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20240607)
    args = parser.parse_args()

    felixla = args.bin_dir / "felixla"
    rng = random.Random(args.seed)
    failures = 0
    for case in range(args.cases):
        problems: list[str] = []
        with tempfile.TemporaryDirectory(prefix="felixla-ancestry-export.") as name:
            try:
                check_case(felixla, pathlib.Path(name), rng, problems)
            except subprocess.CalledProcessError as error:
                problems.append(f"{error.cmd[1:4]} failed: {error.stderr.strip()[-300:]}")
        if problems:
            failures += 1
            print(f"FAIL case {case}:")
            for problem in problems[:3]:
                print(f"  {problem}")

    problems = []
    with tempfile.TemporaryDirectory(prefix="felixla-ancestry-merge.") as name:
        try:
            check_merge_keeps_names(felixla, pathlib.Path(name), problems)
        except subprocess.CalledProcessError as error:
            problems.append(f"{error.cmd[1:4]} failed: {error.stderr.strip()[-300:]}")
    if problems:
        failures += 1
        print("FAIL merged prefix:")
        for problem in problems[:3]:
            print(f"  {problem}")

    print(f"{args.cases + 1 - failures}/{args.cases + 1} ancestry export cases matched")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
