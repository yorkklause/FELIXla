#!/usr/bin/env python3
"""Intense keep/extract/BED regression tests for FELIXla FLARE conversion."""

from __future__ import annotations

import argparse
import gzip
import math
import pathlib
import re
import shutil
import struct
import subprocess
import sys
import tempfile


N_ANCESTRIES = 3
ALLELES = ["A", "C", "G", "T"]


def fail(message: str) -> None:
    raise AssertionError(message)


def run(cmd: list[str], *, expect_fail: bool = False, contains: str | None = None) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if expect_fail:
        if result.returncode == 0:
            fail(f"expected failure but command succeeded: {' '.join(cmd)}")
        if contains and contains not in (result.stdout + result.stderr):
            fail(
                f"expected failed command output to contain {contains!r}\n"
                f"command: {' '.join(cmd)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        return result

    if result.returncode != 0:
        fail(
            f"command failed: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def opener(path: pathlib.Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def alt_list(ref: str, count: int, variant_index: int) -> list[str]:
    choices = [a for a in ALLELES if a != ref]
    offset = variant_index % len(choices)
    rotated = choices[offset:] + choices[:offset]
    return rotated[:count]


def build_inputs(work: pathlib.Path):
    samples = [f"s{i:03d}" for i in range(1, 38)]
    endpoints = {
        "chr1": [250, 500, 750, 1000],
        "chr2": [230, 520, 770, 1000],
    }
    chrom_offset = {"chr1": 0, "chr2": 1}

    def ancestry(chrom: str, pos: int, sample_index: int, hap: int) -> int:
        endpoint_index = next(i for i, end in enumerate(endpoints[chrom]) if pos <= end)
        if hap == 0:
            return (sample_index + endpoint_index + chrom_offset[chrom]) % N_ANCESTRIES
        return (2 * sample_index + endpoint_index + chrom_offset[chrom]) % N_ANCESTRIES

    records = []
    positions = {
        "chr1": [40, 85, 130, 175, 220, 265, 310, 355, 410, 465, 530, 585, 640, 695, 760, 825, 890, 955],
        "chr2": [35, 95, 155, 215, 275, 335, 395, 455, 535, 595, 655, 715, 775, 835, 895, 955],
    }
    variant_index = 0
    for chrom in ["chr1", "chr2"]:
        for pos in positions[chrom]:
            ref = ALLELES[(variant_index + (1 if chrom == "chr2" else 0)) % len(ALLELES)]
            alt_count = [1, 2, 1, 3, 1, 2, 1][variant_index % 7]
            alts = alt_list(ref, alt_count, variant_index)
            raw_id = "." if variant_index % 10 == 4 else f"v{variant_index:03d}"
            n_alleles = 1 + len(alts)
            gts: list[tuple[int, int]] = []
            for sample_index, _sample in enumerate(samples):
                if variant_index % 13 == 0:
                    a0 = 1 if sample_index == 0 else 0
                    a1 = 0
                elif variant_index % 17 == 0 and n_alleles > 2:
                    a0 = 2 if sample_index == len(samples) - 1 else 0
                    a1 = 0
                else:
                    a0 = (sample_index + 2 * variant_index) % n_alleles
                    a1 = (2 * sample_index + variant_index + 1) % n_alleles
                gts.append((a0, a1))
            records.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "id": raw_id,
                    "ref": ref,
                    "alts": alts,
                    "gts": gts,
                    "variant_index": variant_index,
                }
            )
            variant_index += 1

    genotype_path = work / "synthetic.phased.vcf"
    with genotype_path.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##contig=<ID=chr1>\n")
        out.write("##contig=<ID=chr2>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples))
        out.write("\n")
        for record in records:
            sample_gts = [f"{a0}|{a1}" for a0, a1 in record["gts"]]
            out.write(
                f"{record['chrom']}\t{record['pos']}\t{record['id']}\t"
                f"{record['ref']}\t{','.join(record['alts'])}\t.\tPASS\t.\tGT\t"
                + "\t".join(sample_gts)
                + "\n"
            )

    flare_path = work / "synthetic.flare.vcf"
    with flare_path.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##contig=<ID=chr1>\n")
        out.write("##contig=<ID=chr2>\n")
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First haplotype local ancestry">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second haplotype local ancestry">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples))
        out.write("\n")
        for chrom in ["chr1", "chr2"]:
            for pos in endpoints[chrom]:
                sample_ans = [
                    f"{ancestry(chrom, pos, sample_index, 0)}:{ancestry(chrom, pos, sample_index, 1)}"
                    for sample_index, _sample in enumerate(samples)
                ]
                out.write(f"{chrom}\t{pos}\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t")
                out.write("\t".join(sample_ans))
                out.write("\n")

    return samples, records, ancestry, genotype_path, flare_path


def split_id(record: dict, alt: str) -> str:
    if record["id"] == ".":
        return f"{record['chrom']}:{record['pos']}:{record['ref']}:{alt}"
    return f"{record['id']}_{record['ref']}_{alt}"


def split_records(
    records: list[dict],
    sample_indices: list[int],
    *,
    selected: set[tuple[str, int, str, str]] | None = None,
    region: tuple[str, int, int] | None = None,
    bed_intervals: list[tuple[str, int, int]] | None = None,
) -> list[dict]:
    result = []
    for record in records:
        if region:
            r_chrom, r_start, r_end = region
            if record["chrom"] != r_chrom or not (r_start <= record["pos"] <= r_end):
                continue
        if bed_intervals is not None and not any(
            record["chrom"] == chrom and start <= record["pos"] <= end
            for chrom, start, end in bed_intervals
        ):
            continue
        for alt_index, alt in enumerate(record["alts"], start=1):
            key = (record["chrom"], record["pos"], record["ref"], alt)
            if selected is not None and key not in selected:
                continue
            gts = []
            for sample_index in sample_indices:
                a0, a1 = record["gts"][sample_index]
                gts.append(f"{1 if a0 == alt_index else 0}|{1 if a1 == alt_index else 0}")
            result.append(
                {
                    "global_index": len(result),
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "id": split_id(record, alt),
                    "ref": record["ref"],
                    "alt": alt,
                    "alt_index": alt_index,
                    "gts": gts,
                    "record": record,
                }
            )
    return result


def read_roundtrip_vcf(path: pathlib.Path):
    with opener(path) as fh:
        lines = [line.rstrip("\n") for line in fh]
    header = next(line.split("\t") for line in lines if line.startswith("#CHROM"))
    rows = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows.append(
            {
                "chrom": fields[0],
                "pos": int(fields[1]),
                "id": fields[2],
                "ref": fields[3],
                "alt": fields[4],
                "gts": fields[9:],
            }
        )
    return header[9:], rows


def assert_vcf_matches(path: pathlib.Path, expected_samples: list[str], expected: list[dict]) -> None:
    actual_samples, actual_rows = read_roundtrip_vcf(path)
    if actual_samples != expected_samples:
        fail(f"samples mismatch for {path}: {actual_samples} != {expected_samples}")
    compact_actual = [
        (r["chrom"], r["pos"], r["id"], r["ref"], r["alt"], tuple(r["gts"]))
        for r in actual_rows
    ]
    compact_expected = [
        (r["chrom"], r["pos"], r["id"], r["ref"], r["alt"], tuple(r["gts"]))
        for r in expected
    ]
    if compact_actual != compact_expected:
        fail(f"VCF rows mismatch for {path}\nactual={compact_actual[:8]}\nexpected={compact_expected[:8]}")


def read_meta(prefix: pathlib.Path) -> dict[str, str]:
    meta = {}
    with pathlib.Path(str(prefix) + ".meta").open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2:
                meta[parts[0]] = parts[1]
    return meta


def read_ancestry_blocks(prefix: pathlib.Path) -> list[dict]:
    data = pathlib.Path(str(prefix) + ".ancblock.mks").read_bytes()
    if data[:8] != b"TRANMKS1":
        fail(f"bad ancestry marker magic for {prefix}: {data[:8]!r}")
    records = []
    offset = 8
    while offset < len(data):
        if offset + 8 > len(data):
            fail(f"truncated ancestry marker record in {prefix}")
        block_id, chrom_length = struct.unpack_from("<II", data, offset)
        offset += 8
        if offset + chrom_length + 24 > len(data):
            fail(f"truncated ancestry marker string in {prefix}")
        chrom = data[offset : offset + chrom_length].decode()
        offset += chrom_length
        start, end, anc_offset = struct.unpack_from("<qqQ", data, offset)
        offset += 24
        records.append(
            {
                "block_id": block_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "anc_offset": anc_offset,
            }
        )
    return records


def read_samples(prefix: pathlib.Path) -> list[str]:
    return pathlib.Path(str(prefix) + ".samples").read_text().splitlines()


def export_prefix(bin_dir: pathlib.Path, prefix: pathlib.Path) -> pathlib.Path:
    out_vcf = pathlib.Path(str(prefix) + ".roundtrip.vcf.gz")
    run([str(bin_dir / "felixla"), "to-vcf", str(prefix), str(out_vcf)])
    return out_vcf


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def build_spaced_indexed_vcf(bin_dir: pathlib.Path, work: pathlib.Path) -> pathlib.Path:
    genotype = work / "spaced.genotypes.vcf"
    flare = work / "spaced.flare.vcf"
    genotype.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr7>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr7\t16000000\tleft\tA\tT\t.\tPASS\t.\tGT\t0|1\n"
        "chr7\t34000000\tright\tG\tC\t.\tPASS\t.\tGT\t1|0\n"
    )
    flare.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr7>\n"
        '##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First ancestry">\n'
        '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second ancestry">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "chr7\t16000000\t.\tA\tT\t.\tPASS\t.\tAN1:AN2\t0:1\n"
    )
    prefix = work / "spaced"
    run(
        [
            str(bin_dir / "felixla"),
            "from-flare",
            str(genotype),
            str(flare),
            "2",
            "1",
            str(prefix),
        ]
    )
    return export_prefix(bin_dir, prefix)


def check_tbi_chunks(
    bin_dir: pathlib.Path,
    work: pathlib.Path,
    indexed_vcf: pathlib.Path,
    expected_variants: list[dict],
) -> None:
    chunker = bin_dir / "vcf_tbi_chunks"
    quoted_vcf = work / "quoted'phase.vcf.gz"
    shutil.copyfile(indexed_vcf, quoted_vcf)
    shutil.copyfile(pathlib.Path(str(indexed_vcf) + ".tbi"), pathlib.Path(str(quoted_vcf) + ".tbi"))

    manifest = work / "tbi.chunks.tsv"
    commands = work / "tbi.commands.txt"
    chunk_bp = 20_000
    run(
        [
            str(chunker),
            "--phase-vcf",
            str(quoted_vcf),
            "--chunk-bp",
            str(chunk_bp),
            "--out",
            str(manifest),
            "--command-template",
            "worker --vcf {phase_vcf_q} --region {region_q} --global {global_chunk0} --local {chrom_chunk0}",
            "--commands-out",
            str(commands),
        ]
    )

    header, *raw_rows = [line.split("\t") for line in manifest.read_text().splitlines()]
    expected_header = [
        "global_chunk",
        "chrom_chunk",
        "chrom",
        "start",
        "end",
        "region",
        "contig_first_pos",
        "contig_last_pos",
        "contig_records",
    ]
    if header != expected_header:
        fail(f"unexpected vcf_tbi_chunks header: {header}")

    by_chrom: dict[str, list[dict]] = {}
    for variant in expected_variants:
        by_chrom.setdefault(variant["chrom"], []).append(variant)

    expected_rows: list[list[str]] = []
    global_chunk = 0
    for chrom, variants in by_chrom.items():
        positions = [variant["pos"] for variant in variants]
        first = min(positions)
        last = max(positions)
        index_bin_bp = 16_384
        index_first = ((first - 1) // index_bin_bp) * index_bin_bp + 1
        index_last = ((last - 1) // index_bin_bp + 1) * index_bin_bp
        contig_start = ((index_first - 1) // chunk_bp) * chunk_bp + 1
        contig_end = ((index_last - 1) // chunk_bp + 1) * chunk_bp
        start = contig_start
        chrom_chunk = 0
        while start <= contig_end:
            end = start + chunk_bp - 1
            global_chunk += 1
            chrom_chunk += 1
            expected_rows.append(
                [
                    str(global_chunk),
                    str(chrom_chunk),
                    chrom,
                    str(start),
                    str(end),
                    f"{chrom}:{start}-{end}",
                    str(contig_start),
                    str(contig_end),
                    str(len(variants)),
                ]
            )
            start = end + 1

    if raw_rows != expected_rows:
        fail(f"vcf_tbi_chunks rows differ\nactual={raw_rows}\nexpected={expected_rows}")

    # Every split record belongs to exactly one closed chunk, including boundaries.
    for chrom, variants in by_chrom.items():
        chrom_rows = [row for row in raw_rows if row[2] == chrom]
        for variant in variants:
            owners = [row for row in chrom_rows if int(row[3]) <= variant["pos"] <= int(row[4])]
            if len(owners) != 1:
                fail(f"{chrom}:{variant['pos']} has {len(owners)} chunk owners: {owners}")

    expected_commands = []
    for row in expected_rows:
        expected_commands.append(
            f"worker --vcf {shell_quote(str(quoted_vcf))} --region {shell_quote(row[5])} "
            f"--global {int(row[0]):04d} --local {int(row[1]):04d}"
        )
    actual_commands = commands.read_text().splitlines()
    if actual_commands != expected_commands:
        fail(f"vcf_tbi_chunks commands differ\nactual={actual_commands}\nexpected={expected_commands}")

    chr2_manifest = work / "tbi.chr2.mb.tsv"
    run(
        [
            str(chunker),
            "--phase-vcf",
            str(quoted_vcf),
            "--tbi",
            str(quoted_vcf) + ".tbi",
            "--chrom",
            "chr2",
            "--chunk-mb",
            "1",
            "--out",
            str(chr2_manifest),
        ]
    )
    chr2_rows = [line.split("\t") for line in chr2_manifest.read_text().splitlines()]
    chr2_variants = by_chrom["chr2"]
    if chr2_rows[1:] != [[
        "1",
        "1",
        "chr2",
        "1",
        "1000000",
        "chr2:1-1000000",
        "1",
        "1000000",
        str(len(chr2_variants)),
    ]]:
        fail(f"unexpected --chrom/--chunk-mb output: {chr2_rows}")

    # Planning is index-only: an explicit index works even when the VCF body is absent.
    detached_index = work / "detached.tbi"
    shutil.copyfile(pathlib.Path(str(quoted_vcf) + ".tbi"), detached_index)
    absent_vcf = work / "absent.vcf.gz"
    detached_manifest = work / "detached.chunks.tsv"
    run(
        [
            str(chunker),
            "--phase-vcf",
            str(absent_vcf),
            "--tbi",
            str(detached_index),
            "--chunk-mb",
            "1",
            "--out",
            str(detached_manifest),
        ]
    )
    detached_rows = [line.split("\t") for line in detached_manifest.read_text().splitlines()]
    if detached_rows[0] != expected_header or len(detached_rows) != 3:
        fail(f"unexpected detached-index output: {detached_rows}")
    if [row[2:9] for row in detached_rows[1:]] != [
        ["chr1", "1", "1000000", "chr1:1-1000000", "1", "1000000", str(len(by_chrom["chr1"]))],
        ["chr2", "1", "1000000", "chr2:1-1000000", "1", "1000000", str(len(by_chrom["chr2"]))],
    ]:
        fail(f"detached index did not preserve chunk coverage: {detached_rows}")

    spaced_vcf = build_spaced_indexed_vcf(bin_dir, work)
    spaced_index = work / "spaced.detached.tbi"
    shutil.copyfile(pathlib.Path(str(spaced_vcf) + ".tbi"), spaced_index)
    spaced_manifest = work / "spaced.chunks.tsv"
    run(
        [
            str(chunker),
            "--phase-vcf",
            str(work / "spaced.absent.vcf.gz"),
            "--tbi",
            str(spaced_index),
            "--chunk-mb",
            "5",
            "--out",
            str(spaced_manifest),
        ]
    )
    spaced_rows = [line.split("\t") for line in spaced_manifest.read_text().splitlines()]
    expected_spaced = [expected_header]
    for chunk_index, start in enumerate(range(15_000_001, 35_000_001, 5_000_000), 1):
        end = start + 5_000_000 - 1
        expected_spaced.append(
            [
                str(chunk_index),
                str(chunk_index),
                "chr7",
                str(start),
                str(end),
                f"chr7:{start}-{end}",
                "15000001",
                "35000000",
                "2",
            ]
        )
    if spaced_rows != expected_spaced:
        fail(f"spaced index-only chunks differ\nactual={spaced_rows}\nexpected={expected_spaced}")

    run(
        [str(chunker), "--phase-vcf", str(quoted_vcf), "--chunk-bp", "0", "--out", str(work / "bad0.tsv")],
        expect_fail=True,
        contains="must be greater than zero",
    )
    run(
        [str(chunker), "--phase-vcf", str(quoted_vcf), "--chunk-bp", "100", "--chrom", "chr404", "--out", str(work / "badchrom.tsv")],
        expect_fail=True,
        contains="absent from the tabix index",
    )
    run(
        [str(chunker), "--phase-vcf", str(quoted_vcf), "--tbi", str(work / "missing.tbi"), "--chunk-bp", "100", "--out", str(work / "badindex.tsv")],
        expect_fail=True,
        contains="cannot load phased VCF index",
    )
    run(
        [str(chunker), "--phase-vcf", str(quoted_vcf), "--chunk-bp", "100", "--out", str(work / "badtemplate.tsv"), "--command-template", "echo {region}"],
        expect_fail=True,
        contains="must be supplied together",
    )


def make_selected_alleles(records: list[dict]) -> set[tuple[str, int, str, str]]:
    selected: set[tuple[str, int, str, str]] = set()
    multi_records = [record for record in records if len(record["alts"]) > 1]
    first_multi = multi_records[0]
    selected.add((first_multi["chrom"], first_multi["pos"], first_multi["ref"], first_multi["alts"][1]))

    second_multi = multi_records[3]
    for alt in second_multi["alts"][:2]:
        selected.add((second_multi["chrom"], second_multi["pos"], second_multi["ref"], alt))

    for record in records[5::6][:8]:
        selected.add((record["chrom"], record["pos"], record["ref"], record["alts"][0]))

    zero_after_keep = next(record for record in records if record["variant_index"] % 17 == 0 and len(record["alts"]) > 1)
    selected.add((zero_after_keep["chrom"], zero_after_keep["pos"], zero_after_keep["ref"], zero_after_keep["alts"][1]))
    return selected


def write_extract_files(
    work: pathlib.Path,
    selected: set[tuple[str, int, str, str]],
    stem: str = "extract.sites",
) -> tuple[pathlib.Path, pathlib.Path]:
    by_site: dict[tuple[str, int, str], list[str]] = {}
    for chrom, pos, ref, alt in sorted(selected):
        by_site.setdefault((chrom, pos, ref), []).append(alt)

    pvar = work / f"{stem}.pvar"
    with pvar.open("w") as out:
        out.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for (chrom, pos, ref), alts in sorted(by_site.items()):
            out.write(
                f"{chrom}\t{pos}\tpvar_id_must_be_ignored_{chrom}_{pos}"
                f"\t{ref}\t{','.join(alts)}\n"
            )

    vcf = work / f"{stem}.vcf"
    with vcf.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tignored\n")
        for (chrom, pos, ref), alts in sorted(by_site.items()):
            out.write(
                f"{chrom}\t{pos}\tvcf_id_must_be_ignored_{chrom}_{pos}"
                f"\t{ref}\t{','.join(alts)}\t.\tPASS\t.\tGT\t0|0\n"
            )

    vcfgz = work / f"{stem}.vcf.gz"
    with vcf.open("rb") as src, gzip.open(vcfgz, "wb") as dst:
        shutil.copyfileobj(src, dst)

    return pvar, vcfgz


def write_scrambled_extract_pvar(
    work: pathlib.Path,
    selected: set[tuple[str, int, str, str]],
) -> pathlib.Path:
    by_site: dict[tuple[str, int, str], list[str]] = {}
    for chrom, pos, ref, alt in sorted(selected):
        by_site.setdefault((chrom, pos, ref), []).append(alt)

    rows: list[tuple[str, int, str, str]] = []
    split_multiallelic = False
    for (chrom, pos, ref), alts in reversed(sorted(by_site.items())):
        reversed_alts = list(reversed(alts))
        if len(reversed_alts) > 1 and not split_multiallelic:
            for alt in reversed_alts:
                rows.append((chrom, pos, ref, alt))
            split_multiallelic = True
        else:
            rows.append((chrom, pos, ref, ",".join(reversed_alts)))

    if not split_multiallelic:
        fail("scrambled extract fixture lacks a multiallelic target site")

    pvar = work / "extract.sites.scrambled.pvar"
    with pvar.open("w") as out:
        out.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for row_i, (chrom, pos, ref, alts) in enumerate(rows):
            out.write(f"{chrom}\t{pos}\tignored_scrambled_{row_i}\t{ref}\t{alts}\n")
    return pvar


def write_keep_file(work: pathlib.Path, samples: list[str], keep_indices: list[int]) -> pathlib.Path:
    keep = work / "keep.samples"
    scrambled = [keep_indices[i] for i in [5, 0, 16, 2, 9, 20, 1, 12, 26, 3, 7, 10, 30, 4, 6, 8, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 24, 25, 27, 28, 29, 31, 32]]
    keep.write_text("\n".join(samples[i] for i in scrambled) + "\n")
    return keep


def write_reordered_subset_flare(
    work: pathlib.Path,
    flare_path: pathlib.Path,
    sample_order_indices: list[int],
) -> pathlib.Path:
    out_path = work / "synthetic.flare.reordered_subset.vcf"
    with flare_path.open() as src, out_path.open("w") as out:
        for line in src:
            line = line.rstrip("\n")
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                header_samples = fields[9:]
                reordered_samples = [header_samples[i] for i in sample_order_indices]
                out.write("\t".join(fields[:9] + reordered_samples) + "\n")
            elif line.startswith("#"):
                out.write(line + "\n")
            else:
                fields = line.split("\t")
                sample_values = fields[9:]
                reordered_values = [sample_values[i] for i in sample_order_indices]
                out.write("\t".join(fields[:9] + reordered_values) + "\n")
    return out_path


def write_extra_format_inputs(
    work: pathlib.Path,
    genotype_path: pathlib.Path,
    flare_path: pathlib.Path,
) -> tuple[pathlib.Path, pathlib.Path]:
    genotype_out = work / "synthetic.phased.extra_format.vcf"
    with genotype_path.open() as src, genotype_out.open("w") as out:
        for line in src:
            if line.startswith("#CHROM"):
                out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
                out.write(line)
            elif line.startswith("#"):
                out.write(line)
            else:
                fields = line.rstrip("\n").split("\t")
                fields[8] = "GT:DP"
                fields[9:] = [f"{value}:30" for value in fields[9:]]
                out.write("\t".join(fields) + "\n")

    flare_out = work / "synthetic.flare.extra_format.vcf"
    with flare_path.open() as src, flare_out.open("w") as out:
        for line in src:
            if line.startswith("#CHROM"):
                out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                out.write('##FORMAT=<ID=ANP1,Number=1,Type=Float,Description="First ancestry probability">\n')
                out.write('##FORMAT=<ID=ANP2,Number=1,Type=Float,Description="Second ancestry probability">\n')
                out.write(line)
            elif line.startswith("#"):
                out.write(line)
            else:
                fields = line.rstrip("\n").split("\t")
                fields[8] = "GT:AN1:AN2:ANP1:ANP2"
                fields[9:] = [f"0|0:{value}:0.99:0.98" for value in fields[9:]]
                out.write("\t".join(fields) + "\n")
    return genotype_out, flare_out


def append_extra_sample_to_first_record(
    source: pathlib.Path,
    destination: pathlib.Path,
    value: str,
) -> pathlib.Path:
    added = False
    with source.open() as src, destination.open("w") as out:
        for line in src:
            if not added and not line.startswith("#"):
                out.write(line.rstrip("\n") + "\t" + value + "\n")
                added = True
            else:
                out.write(line)
    if not added:
        fail(f"cannot add malformed sample column to empty VCF: {source}")
    return destination


def query_rows(
    bin_dir: pathlib.Path,
    prefix: pathlib.Path,
    split: dict,
    samples: list[str],
    sample_indices: list[int],
    ancestry,
    *,
    nonzero_only: bool = False,
) -> list[list[str]]:
    cmd = [
        str(bin_dir / "felixla"),
        "--felixla",
        str(prefix),
        "--query",
        f"{split['chrom']}:{split['pos']}",
        "--ref",
        split["ref"],
        "--alt",
        split["alt"],
    ]
    if nonzero_only:
        cmd.append("--nonzero-only")
    actual = run(cmd).stdout.strip().splitlines()
    expected_header = "global_variant_index\tchr\tpos\tid\tref\talt\tsample\tDSALL\tDS1\tDS2\tDS3"
    if not actual or actual[0] != expected_header:
        fail(f"unexpected query header for {prefix}: {actual[:1]}")
    actual_rows = [line.split("\t") for line in actual[1:]]

    expected_rows = []
    for out_i, sample_index in enumerate(sample_indices):
        gt = split["gts"][out_i]
        h0 = 1 if gt[0] == "1" else 0
        h1 = 1 if gt[2] == "1" else 0
        ds = [0, 0, 0]
        ds[ancestry(split["chrom"], split["pos"], sample_index, 0)] += h0
        ds[ancestry(split["chrom"], split["pos"], sample_index, 1)] += h1
        dsall = h0 + h1
        if nonzero_only and dsall == 0:
            continue
        expected_rows.append(
            [
                str(split["global_index"]),
                split["chrom"],
                str(split["pos"]),
                split["id"],
                split["ref"],
                split["alt"],
                samples[sample_index],
                str(dsall),
                *(str(x) for x in ds),
            ]
        )
    if actual_rows != expected_rows:
        fail(f"query mismatch for {prefix} {split['chrom']}:{split['pos']} {split['ref']}>{split['alt']}")
    return actual_rows


def check_queries(
    bin_dir: pathlib.Path,
    prefix: pathlib.Path,
    expected: list[dict],
    samples: list[str],
    sample_indices: list[int],
    ancestry,
) -> None:
    targets = expected[:3] + expected[len(expected) // 2: len(expected) // 2 + 2] + expected[-3:]
    seen = set()
    unique_targets = []
    for target in targets:
        key = (target["chrom"], target["pos"], target["ref"], target["alt"])
        if key not in seen:
            seen.add(key)
            unique_targets.append(target)
    for target in unique_targets:
        query_rows(bin_dir, prefix, target, samples, sample_indices, ancestry)

    nonzero_target = next(target for target in expected if any(gt != "0|0" for gt in target["gts"]))
    query_rows(bin_dir, prefix, nonzero_target, samples, sample_indices, ancestry, nonzero_only=True)


def build_with_cli(
    bin_dir: pathlib.Path,
    genotype_path: pathlib.Path,
    flare_path: pathlib.Path,
    prefix: pathlib.Path,
    *extra: str,
) -> pathlib.Path:
    run(
        [
            str(bin_dir / "felixla"),
            "--phase-vcf",
            str(genotype_path),
            "--flare-vcf",
            str(flare_path),
            "--n-ancestries",
            str(N_ANCESTRIES),
            *extra,
            "--make-felixla",
            "--out",
            str(prefix),
        ]
    )
    return export_prefix(bin_dir, prefix)


def check_int16_gt_encoding(bin_dir: pathlib.Path, work: pathlib.Path) -> None:
    samples = ["wide1", "wide2", "wide3"]
    alts = ["C", "G", "T"] + ["A" + "C" * length for length in range(1, 68)]
    target_index = 64
    target_alt = alts[target_index - 1]

    genotype = work / "int16_gt.phased.vcf"
    with genotype.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        out.write(
            f"chr1\t100\twide\tA\t{','.join(alts)}\t.\tPASS\t.\tGT\t"
            f"0|{target_index}\t{target_index}|0\t{target_index}|{target_index}\n"
        )

    flare = work / "int16_gt.flare.vcf"
    with flare.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First ancestry">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second ancestry">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        out.write("chr1\t100\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t0:1\t1:2\t2:0\n")

    extract = work / "int16_gt.extract.pvar"
    extract.write_text(f"#CHROM\tPOS\tID\tREF\tALT\nchr1\t100\tignored\tA\t{target_alt}\n")
    prefix = work / "int16_gt"
    roundtrip = build_with_cli(
        bin_dir,
        genotype,
        flare,
        prefix,
        "--extract",
        str(extract),
    )
    expected = [{
        "chrom": "chr1",
        "pos": 100,
        "id": f"wide_A_{target_alt}",
        "ref": "A",
        "alt": target_alt,
        "gts": ["0|1", "1|0", "1|1"],
    }]
    assert_vcf_matches(roundtrip, samples, expected)


def check_duplicate_flare_coordinate(bin_dir: pathlib.Path, work: pathlib.Path) -> None:
    samples = ["dup1", "dup2"]
    genotype = work / "duplicate_flare.phased.vcf"
    with genotype.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        out.write("chr1\t100\td1\tA\tC\t.\tPASS\t.\tGT\t0|1\t1|0\n")
        out.write("chr1\t200\td2\tG\tT\t.\tPASS\t.\tGT\t1|1\t0|0\n")

    flare = work / "duplicate_flare.flare.vcf"
    with flare.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First ancestry">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second ancestry">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        out.write("chr1\t100\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t0:1\t1:2\n")
        out.write("chr1\t200\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t2:1\t1:0\n")
        out.write("chr1\t200\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t0:1\t1:2\n")

    prefix = work / "duplicate_flare"
    roundtrip = build_with_cli(bin_dir, genotype, flare, prefix)
    expected = [
        {"chrom": "chr1", "pos": 100, "id": "d1_A_C", "ref": "A", "alt": "C", "gts": ["0|1", "1|0"]},
        {"chrom": "chr1", "pos": 200, "id": "d2_G_T", "ref": "G", "alt": "T", "gts": ["1|1", "0|0"]},
    ]
    assert_vcf_matches(roundtrip, samples, expected)
    expected_ancestry_bytes = N_ANCESTRIES * math.ceil(2 * len(samples) / 64) * 8
    ancestry_bytes = pathlib.Path(str(prefix) + ".ancblock.bin").stat().st_size
    if ancestry_bytes != expected_ancestry_bytes:
        fail(
            "duplicate FLARE coordinate created an unnecessary ancestry block: "
            f"{ancestry_bytes} bytes != {expected_ancestry_bytes}"
        )


def check_padded_multiallelic_extract(bin_dir: pathlib.Path, work: pathlib.Path) -> None:
    samples = ["pad1", "pad2", "pad3", "pad4"]
    genotype_alts = ["A", "AATT", "ATT", "ATTT", "ATTTT", "GT", "TT"]
    extract_alts = ["ATTT", "ATTTT", "AT", "ATTTTT", "AATTT", "TTT", "GTT", "A"]
    source_gts = [(1, 2), (3, 4), (5, 6), (7, 0)]

    genotype = work / "padded_multiallelic.phased.vcf"
    with genotype.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr19>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        # A different REF at the same coordinate must not cause an early mismatch.
        out.write(
            "chr19\t40176167\tunrelated\tC\tG\t.\tPASS\t.\tGT\t"
            "0|1\t0|0\t0|0\t0|0\n"
        )
        out.write(
            "chr19\t40176167\tpadded\tAT\t"
            + ",".join(genotype_alts)
            + "\t.\tPASS\t.\tGT\t"
            + "\t".join(f"{a0}|{a1}" for a0, a1 in source_gts)
            + "\n"
        )

    flare = work / "padded_multiallelic.flare.vcf"
    with flare.open("w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=chr19>\n")
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First ancestry">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second ancestry">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        out.write("chr19\t40176167\t.\tA\tC\t.\tPASS\t.\tAN1:AN2\t0:1\t1:2\t2:0\t0:2\n")

    extract = work / "padded_multiallelic.extract.pvar"
    extract.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\n"
        f"chr19\t40176167\t.\tATT\t{','.join(extract_alts)}\n"
    )

    prefix = work / "padded_multiallelic"
    roundtrip = build_with_cli(
        bin_dir,
        genotype,
        flare,
        prefix,
        "--extract",
        str(extract),
    )
    expected = []
    for alt_index, alt in enumerate(genotype_alts, start=1):
        expected.append({
            "chrom": "chr19",
            "pos": 40176167,
            "id": f"padded_AT_{alt}",
            "ref": "AT",
            "alt": alt,
            "gts": [
                f"{1 if a0 == alt_index else 0}|{1 if a1 == alt_index else 0}"
                for a0, a1 in source_gts
            ],
        })
    assert_vcf_matches(roundtrip, samples, expected)


def bounded_extract_expected_from_split_records(
    split_records_all: list[dict],
    selected: set[tuple[str, int, str, str]],
) -> list[dict]:
    expected = []
    for split in split_records_all:
        key = (split["chrom"], split["pos"], split["ref"], split["alt"])
        if key not in selected:
            continue
        copied = dict(split)
        copied["global_index"] = len(expected)
        copied["alt_index"] = 1
        copied["id"] = f"{split['id']}_{split['ref']}_{split['alt']}"
        expected.append(copied)
    return expected


def bounded_extract_scanned_records(
    split_records_all: list[dict],
    selected: set[tuple[str, int, str, str]],
) -> int:
    bounds: dict[str, list[int]] = {}
    for chrom, pos, _ref, _alt in selected:
        if chrom not in bounds:
            bounds[chrom] = [pos, pos]
        else:
            bounds[chrom][0] = min(bounds[chrom][0], pos)
            bounds[chrom][1] = max(bounds[chrom][1], pos)

    return sum(
        1 for split in split_records_all
        if split["chrom"] in bounds and
        bounds[split["chrom"]][0] <= split["pos"] <= bounds[split["chrom"]][1]
    )


def last_progress_record_count(stderr: str) -> int:
    matches = re.findall(r"Progress: records (\d+)", stderr)
    if not matches:
        fail(f"missing progress records line in stderr:\n{stderr}")
    return int(matches[-1])


def repacked_bed_expected(
    split_records_all: list[dict],
    bed_intervals: list[tuple[str, int, int]],
) -> list[dict]:
    expected = []
    for split in split_records_all:
        if not any(
            split["chrom"] == chrom and start <= split["pos"] <= end
            for chrom, start, end in bed_intervals
        ):
            continue
        copied = dict(split)
        copied["global_index"] = len(expected)
        copied["alt_index"] = 1
        copied["id"] = f"{split['id']}_{split['ref']}_{split['alt']}"
        expected.append(copied)
    return expected


def assert_ancestry_blocks_inside_bed(
    prefix: pathlib.Path,
    bed_intervals: list[tuple[str, int, int]],
) -> None:
    blocks = read_ancestry_blocks(prefix)
    if not blocks:
        fail(f"BED-selected output has no ancestry blocks: {prefix}")
    for block in blocks:
        if not any(
            block["chrom"] == chrom and
            start <= block["start"] <= block["end"] <= end
            for chrom, start, end in bed_intervals
        ):
            fail(f"ancestry block crosses an unselected BED gap: {block}")


def check_extract_bed_basic(
    bin_dir: pathlib.Path,
    work: pathlib.Path,
    samples: list[str],
    records: list[dict],
    ancestry,
    genotype_path: pathlib.Path,
    flare_path: pathlib.Path,
) -> tuple[pathlib.Path, list[tuple[str, int, int]]]:
    bed_path = work / "intervals.unsorted.bed"
    bed_path.write_text(
        "browser position chr1:1-1000\n"
        "chr2\t214\t335\tsecond-contig\n"
        "chr1\t175\t220\tadjacent\n"
        "chr1 39 40 one-base\n"
        "track name=felixla-test\n"
        "chr1\t84\t175\tleft\n"
        "# a comment between records\n"
        "chr2\t0\t35\tzero-start\n"
        "chr1\t500\t641\twide\n"
        "chr1\t500\t600\tduplicate-overlap\n"
    )
    merged_intervals = [
        ("chr1", 40, 40),
        ("chr1", 85, 220),
        ("chr1", 501, 641),
        ("chr2", 1, 35),
        ("chr2", 215, 335),
    ]
    expected = split_records(
        records,
        list(range(len(samples))),
        bed_intervals=merged_intervals,
    )
    prefix = work / "subset.extract_bed"
    roundtrip = build_with_cli(
        bin_dir,
        genotype_path,
        flare_path,
        prefix,
        "--extract-bed",
        str(bed_path),
    )
    assert_vcf_matches(roundtrip, samples, expected)
    check_queries(
        bin_dir,
        prefix,
        expected,
        samples,
        list(range(len(samples))),
        ancestry,
    )

    meta = read_meta(prefix)
    assert meta["extract_bed"] == str(bed_path), meta
    assert meta["extract_bed_coordinates"] == "0-based-half-open", meta
    assert meta["extract_bed_source_intervals"] == "7", meta
    assert meta["extract_bed_merged_intervals"] == "5", meta
    assert meta["extract_bed_selected_intervals"] == "5", meta
    assert_ancestry_blocks_inside_bed(prefix, merged_intervals)

    blocks = read_ancestry_blocks(prefix)
    block_at_40 = [b for b in blocks if b["chrom"] == "chr1" and b["start"] <= 40 <= b["end"]]
    block_at_85 = [b for b in blocks if b["chrom"] == "chr1" and b["start"] <= 85 <= b["end"]]
    if len(block_at_40) != 1 or len(block_at_85) != 1:
        fail(f"BED boundary ancestry blocks missing: pos40={block_at_40}, pos85={block_at_85}")
    if block_at_40[0]["block_id"] == block_at_85[0]["block_id"]:
        fail("one ancestry block incorrectly bridges the BED gap between chr1:40 and chr1:85")

    gz_path = work / "intervals.unsorted.bed.gz"
    with gzip.open(gz_path, "wt") as out:
        out.write(bed_path.read_text())
    gz_prefix = work / "subset.extract_bed_gz"
    gz_roundtrip = build_with_cli(
        bin_dir,
        genotype_path,
        flare_path,
        gz_prefix,
        "--extract-bed",
        str(gz_path),
    )
    assert_vcf_matches(gz_roundtrip, samples, expected)
    assert_ancestry_blocks_inside_bed(gz_prefix, merged_intervals)

    compatibility_prefix = work / "subset.extract_bed_compatibility"
    run(
        [
            str(bin_dir / "felixla"),
            "from-flare",
            str(genotype_path),
            str(flare_path),
            str(N_ANCESTRIES),
            "auto",
            str(compatibility_prefix),
            "--extract-bed",
            str(bed_path),
        ]
    )
    compatibility_vcf = export_prefix(bin_dir, compatibility_prefix)
    assert_vcf_matches(compatibility_vcf, samples, expected)
    assert_ancestry_blocks_inside_bed(compatibility_prefix, merged_intervals)
    return bed_path, merged_intervals


def check_indexed_and_large_extract_bed(
    bin_dir: pathlib.Path,
    work: pathlib.Path,
    samples: list[str],
    full_expected: list[dict],
    indexed_genotype: pathlib.Path,
    flare_path: pathlib.Path,
) -> None:
    sparse_bed = work / "indexed.sparse.bed"
    sparse_bed.write_text("chr1\t639\t640\nchr1\t129\t130\n")
    sparse_intervals = [("chr1", 130, 130), ("chr1", 640, 640)]
    sparse_expected = repacked_bed_expected(full_expected, sparse_intervals)
    sparse_prefix = work / "subset.indexed_extract_bed"
    sparse_result = run(
        [
            str(bin_dir / "felixla"),
            "--phase-vcf",
            str(indexed_genotype),
            "--flare-vcf",
            str(flare_path),
            "--n-ancestries",
            str(N_ANCESTRIES),
            "--extract-bed",
            str(sparse_bed),
            "--make-felixla",
            "--out",
            str(sparse_prefix),
        ]
    )
    if "Using indexed --extract-bed chromosome-span genotype reader" not in sparse_result.stderr:
        fail(f"indexed BED reader was not used:\n{sparse_result.stderr}")
    expected_scanned = sum(
        1 for split in full_expected
        if split["chrom"] == "chr1" and 130 <= split["pos"] <= 640
    )
    scanned = last_progress_record_count(sparse_result.stderr)
    if scanned != expected_scanned:
        fail(f"indexed BED span scanned {scanned} records; expected {expected_scanned}")
    sparse_vcf = export_prefix(bin_dir, sparse_prefix)
    assert_vcf_matches(sparse_vcf, samples, sparse_expected)
    assert_ancestry_blocks_inside_bed(sparse_prefix, sparse_intervals)

    large_bed = work / "large.unsorted.bed"
    large_rows = ["chr1\t129\t130\ttarget"]
    for i in reversed(range(20_000)):
        start0 = 2_000 + 2 * i
        large_rows.append(f"chr1\t{start0}\t{start0 + 1}\tinterval{i}")
    large_bed.write_text("\n".join(large_rows) + "\n")
    large_prefix = work / "subset.large_extract_bed"
    large_result = run(
        [
            str(bin_dir / "felixla"),
            "--phase-vcf",
            str(indexed_genotype),
            "--flare-vcf",
            str(flare_path),
            "--n-ancestries",
            str(N_ANCESTRIES),
            "--extract-bed",
            str(large_bed),
            "--make-felixla",
            "--out",
            str(large_prefix),
        ]
    )
    if "20001 exact interval(s) across 1 chromosome span(s)" not in large_result.stderr:
        fail(f"large BED was not collapsed to one indexed chromosome span:\n{large_result.stderr}")
    large_expected = repacked_bed_expected(full_expected, [("chr1", 130, 130)])
    large_vcf = export_prefix(bin_dir, large_prefix)
    assert_vcf_matches(large_vcf, samples, large_expected)
    large_meta = read_meta(large_prefix)
    assert large_meta["extract_bed_source_intervals"] == "20001", large_meta
    assert large_meta["extract_bed_merged_intervals"] == "20001", large_meta
    assert large_meta["extract_bed_selected_intervals"] == "20001", large_meta


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=pathlib.Path, default=pathlib.Path("bin"))
    parser.add_argument("--work-dir", type=pathlib.Path)
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()

    owned_tmp = None
    if args.work_dir:
        work = args.work_dir
        work.mkdir(parents=True, exist_ok=True)
    else:
        owned_tmp = tempfile.TemporaryDirectory(prefix="felixla-keep-extract-intense.")
        work = pathlib.Path(owned_tmp.name)

    try:
        bin_dir = args.bin_dir.resolve()
        if not (bin_dir / "felixla").exists():
            fail(f"missing tool: {bin_dir / 'felixla'}")
        if not (bin_dir / "vcf_tbi_chunks").exists():
            fail(f"missing tool: {bin_dir / 'vcf_tbi_chunks'}")
        unexpected = sorted(
            path.name for path in bin_dir.iterdir()
            if path.name not in {"felixla", "vcf_tbi_chunks"}
        )
        if unexpected:
            fail(f"unexpected files in {bin_dir}: {unexpected}")

        samples, records, ancestry, genotype_path, flare_path = build_inputs(work)
        check_int16_gt_encoding(bin_dir, work)
        check_duplicate_flare_coordinate(bin_dir, work)
        check_padded_multiallelic_extract(bin_dir, work)
        bed_path, merged_bed_intervals = check_extract_bed_basic(
            bin_dir,
            work,
            samples,
            records,
            ancestry,
            genotype_path,
            flare_path,
        )
        all_indices = list(range(len(samples)))
        keep_indices = [i for i in all_indices if i not in {5, 10, 26, 36}]
        keep_samples = [samples[i] for i in keep_indices]
        keep_path = write_keep_file(work, samples, keep_indices)
        selected = make_selected_alleles(records)
        pvar_path, vcfgz_path = write_extract_files(work, selected)
        scrambled_pvar_path = write_scrambled_extract_pvar(work, selected)

        legacy_prefix = work / "full.legacy"
        run(
            [
                str(bin_dir / "felixla"),
                "from-flare",
                str(genotype_path),
                str(flare_path),
                str(N_ANCESTRIES),
                "auto",
                str(legacy_prefix),
            ]
        )
        legacy_vcf = export_prefix(bin_dir, legacy_prefix)
        full_expected = split_records(records, all_indices)
        assert_vcf_matches(legacy_vcf, samples, full_expected)
        check_tbi_chunks(bin_dir, work, legacy_vcf, full_expected)
        check_queries(bin_dir, legacy_prefix, full_expected, samples, all_indices, ancestry)
        legacy_meta = read_meta(legacy_prefix)
        assert legacy_meta["n_samples"] == str(len(samples)), legacy_meta
        assert legacy_meta["n_words"] == "2", legacy_meta
        assert legacy_meta["rare_threshold"] == str(math.ceil(len(samples) / 32)), legacy_meta
        assert "keep_samples" not in legacy_meta and "extract_sites" not in legacy_meta, legacy_meta
        assert "extract_bed" not in legacy_meta, legacy_meta

        full_bed = work / "full_coverage.bed"
        full_bed.write_text("chr2\t0\t1000\nchr1\t0\t1000\n")
        full_bed_prefix = work / "full.extract_bed"
        full_bed_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            full_bed_prefix,
            "--extract-bed",
            str(full_bed),
        )
        assert_vcf_matches(full_bed_vcf, samples, full_expected)
        for suffix in [
            ".common.geno.bin",
            ".common.variant.mks",
            ".common.variant.idx",
            ".rare.carrier.bin",
            ".rare.variant.mks",
            ".rare.variant.idx",
            ".ancblock.bin",
            ".ancblock.mks",
            ".ancblock.idx",
            ".samples",
        ]:
            legacy_bytes = pathlib.Path(str(legacy_prefix) + suffix).read_bytes()
            bed_bytes = pathlib.Path(str(full_bed_prefix) + suffix).read_bytes()
            if bed_bytes != legacy_bytes:
                fail(f"full-coverage BED changed legacy output bytes for {suffix}")

        check_indexed_and_large_extract_bed(
            bin_dir,
            work,
            samples,
            full_expected,
            legacy_vcf,
            flare_path,
        )

        cli_prefix = work / "full.cli"
        cli_vcf = build_with_cli(bin_dir, genotype_path, flare_path, cli_prefix)
        assert_vcf_matches(cli_vcf, samples, full_expected)
        if read_roundtrip_vcf(cli_vcf) != read_roundtrip_vcf(legacy_vcf):
            fail("full CLI output differs from legacy positional output")

        extra_genotype, extra_flare = write_extra_format_inputs(
            work, genotype_path, flare_path
        )
        extra_format_prefix = work / "full.extra_format"
        extra_format_vcf = build_with_cli(
            bin_dir, extra_genotype, extra_flare, extra_format_prefix
        )
        if read_roundtrip_vcf(extra_format_vcf) != read_roundtrip_vcf(legacy_vcf):
            fail("extra genotype/FLARE FORMAT fields changed output")

        subset_expected = split_records(records, keep_indices, selected=selected)
        subset_prefix = work / "subset.vcfgz"
        subset_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            subset_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(vcfgz_path),
        )
        assert_vcf_matches(subset_vcf, keep_samples, subset_expected)
        check_queries(bin_dir, subset_prefix, subset_expected, samples, keep_indices, ancestry)
        subset_meta = read_meta(subset_prefix)
        assert subset_meta["n_samples"] == str(len(keep_indices)), subset_meta
        assert subset_meta["n_words"] == "2", subset_meta
        assert subset_meta["rare_threshold"] == str(math.ceil(len(keep_indices) / 32)), subset_meta
        assert subset_meta["keep_samples"] == str(keep_path), subset_meta
        assert subset_meta["extract_sites"] == str(vcfgz_path), subset_meta
        assert read_samples(subset_prefix) == keep_samples

        pvar_prefix = work / "subset.pvar"
        pvar_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            pvar_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(pvar_path),
        )
        assert_vcf_matches(pvar_vcf, keep_samples, subset_expected)
        if read_roundtrip_vcf(pvar_vcf) != read_roundtrip_vcf(subset_vcf):
            fail("PVAR extract output differs from gzipped VCF extract output")

        scrambled_prefix = work / "subset.scrambled_pvar"
        scrambled_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            scrambled_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(scrambled_pvar_path),
        )
        assert_vcf_matches(scrambled_vcf, keep_samples, subset_expected)
        if read_roundtrip_vcf(scrambled_vcf) != read_roundtrip_vcf(pvar_vcf):
            fail("scrambled multiallelic PVAR changed extracted variant order or values")

        missing_from_flare = {keep_indices[1], keep_indices[8], keep_indices[-2]}
        flare_subset_indices = [i for i in all_indices if i not in missing_from_flare]
        reordered_flare_indices = (
            flare_subset_indices[3::5]
            + flare_subset_indices[0::5]
            + flare_subset_indices[4::5]
            + flare_subset_indices[1::5]
            + flare_subset_indices[2::5]
        )
        reordered_flare_path = write_reordered_subset_flare(work, flare_path, reordered_flare_indices)
        intersection_indices = [i for i in keep_indices if i not in missing_from_flare]
        intersection_samples = [samples[i] for i in intersection_indices]
        intersection_expected = split_records(records, intersection_indices, selected=selected)
        intersection_prefix = work / "subset.reordered_flare_intersection"
        intersection_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            reordered_flare_path,
            intersection_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(pvar_path),
        )
        assert_vcf_matches(intersection_vcf, intersection_samples, intersection_expected)
        check_queries(bin_dir, intersection_prefix, intersection_expected, samples, intersection_indices, ancestry)
        intersection_meta = read_meta(intersection_prefix)
        assert intersection_meta["n_samples"] == str(len(intersection_indices)), intersection_meta
        assert intersection_meta["n_words"] == str(math.ceil(2 * len(intersection_indices) / 64)), intersection_meta
        assert intersection_meta["rare_threshold"] == str(math.ceil(len(intersection_indices) / 32)), intersection_meta
        assert read_samples(intersection_prefix) == intersection_samples

        indexed_selected = {
            (split["chrom"], split["pos"], split["ref"], split["alt"])
            for split in [full_expected[1], full_expected[7], full_expected[19], full_expected[-3]]
        }
        indexed_extract_pvar, _indexed_extract_vcfgz = write_extract_files(
            work,
            indexed_selected,
            "indexed.extract.sites",
        )
        indexed_expected = bounded_extract_expected_from_split_records(full_expected, indexed_selected)
        indexed_prefix = work / "subset.indexed_extract"
        indexed_result = run(
            [
                str(bin_dir / "felixla"),
                "--phase-vcf",
                str(legacy_vcf),
                "--flare-vcf",
                str(flare_path),
                "--n-ancestries",
                str(N_ANCESTRIES),
                "--extract",
                str(indexed_extract_pvar),
                "--make-felixla",
                "--out",
                str(indexed_prefix),
            ]
        )
        if "Using indexed --extract bounded genotype reader" not in indexed_result.stderr:
            fail(f"indexed bounded --extract reader was not used:\n{indexed_result.stderr}")
        scanned_records = last_progress_record_count(indexed_result.stderr)
        expected_scanned_records = bounded_extract_scanned_records(full_expected, indexed_selected)
        if scanned_records != expected_scanned_records:
            fail(
                f"indexed --extract scanned {scanned_records} records; "
                f"expected exactly {expected_scanned_records} records inside bounded extract intervals"
            )
        if scanned_records >= len(full_expected):
            fail(f"bounded indexed --extract scanned the full VCF: {scanned_records} of {len(full_expected)} records")
        indexed_vcf = export_prefix(bin_dir, indexed_prefix)
        assert_vcf_matches(indexed_vcf, samples, indexed_expected)

        indexed_keep_expected = bounded_extract_expected_from_split_records(
            split_records(records, keep_indices),
            indexed_selected,
        )
        indexed_keep_prefix = work / "subset.indexed_extract_keep"
        indexed_keep_vcf = build_with_cli(
            bin_dir,
            legacy_vcf,
            flare_path,
            indexed_keep_prefix,
            "--keep",
            str(keep_path),
            "--extract",
            str(indexed_extract_pvar),
        )
        assert_vcf_matches(indexed_keep_vcf, keep_samples, indexed_keep_expected)
        assert read_samples(indexed_keep_prefix) == keep_samples

        region = ("chr1", 200, 650)
        region_expected = split_records(records, keep_indices, selected=selected, region=region)
        region_prefix = work / "subset.region"
        region_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            region_prefix,
            "--region",
            f"{region[0]}:{region[1]}-{region[2]}",
            "--keep",
            str(keep_path),
            "--extract",
            str(vcfgz_path),
        )
        assert_vcf_matches(region_vcf, keep_samples, region_expected)
        check_queries(bin_dir, region_prefix, region_expected, samples, keep_indices, ancestry)

        bed_region = ("chr1", 100, 600)
        bed_region_intervals = [("chr1", 100, 220), ("chr1", 501, 600)]
        bed_intersection_expected = split_records(
            records,
            keep_indices,
            selected=selected,
            region=bed_region,
            bed_intervals=merged_bed_intervals,
        )
        bed_intersection_prefix = work / "subset.bed_region_extract_keep"
        bed_intersection_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            bed_intersection_prefix,
            "--region",
            f"{bed_region[0]}:{bed_region[1]}-{bed_region[2]}",
            "--keep",
            str(keep_path),
            "--extract",
            str(pvar_path),
            "--extract-bed",
            str(bed_path),
        )
        assert_vcf_matches(
            bed_intersection_vcf,
            keep_samples,
            bed_intersection_expected,
        )
        check_queries(
            bin_dir,
            bed_intersection_prefix,
            bed_intersection_expected,
            samples,
            keep_indices,
            ancestry,
        )
        assert_ancestry_blocks_inside_bed(
            bed_intersection_prefix,
            bed_region_intervals,
        )
        bed_intersection_meta = read_meta(bed_intersection_prefix)
        assert bed_intersection_meta["selected_region"] == "chr1:100-600", bed_intersection_meta
        assert bed_intersection_meta["extract_bed_source_intervals"] == "7", bed_intersection_meta
        assert bed_intersection_meta["extract_bed_merged_intervals"] == "2", bed_intersection_meta
        assert bed_intersection_meta["extract_bed_selected_intervals"] == "2", bed_intersection_meta

        bad_keep_dup = work / "bad.keep.dup"
        bad_keep_dup.write_text(f"{samples[0]}\n{samples[0]}\n")
        build_bad_base = [
            str(bin_dir / "felixla"),
            "--phase-vcf",
            str(genotype_path),
            "--flare-vcf",
            str(flare_path),
            "--n-ancestries",
            str(N_ANCESTRIES),
            "--make-felixla",
            "--out",
        ]
        run([*build_bad_base, str(work / "bad.keep.dup.out"), "--keep", str(bad_keep_dup)], expect_fail=True, contains="duplicate sample ID")

        bad_keep_extra = work / "bad.keep.extra"
        bad_keep_extra.write_text(f"{samples[0]} extra\n")
        run([*build_bad_base, str(work / "bad.keep.extra.out"), "--keep", str(bad_keep_extra)], expect_fail=True, contains="one sample ID per line")

        bad_keep_missing = work / "bad.keep.missing"
        bad_keep_missing.write_text("not_a_sample\n")
        run([*build_bad_base, str(work / "bad.keep.missing.out"), "--keep", str(bad_keep_missing)], expect_fail=True, contains="absent from genotype VCF")

        extra_genotype_sample = append_extra_sample_to_first_record(
            genotype_path,
            work / "bad.extra_genotype_sample.vcf",
            "0|0",
        )
        run(
            [
                str(bin_dir / "felixla"),
                "--phase-vcf",
                str(extra_genotype_sample),
                "--flare-vcf",
                str(flare_path),
                "--n-ancestries",
                str(N_ANCESTRIES),
                "--make-felixla",
                "--out",
                str(work / "bad.extra_genotype_sample.out"),
            ],
            expect_fail=True,
            contains="GT sample count mismatch",
        )

        extra_flare_sample = append_extra_sample_to_first_record(
            flare_path,
            work / "bad.extra_flare_sample.vcf",
            "0:0",
        )
        run(
            [
                str(bin_dir / "felixla"),
                "--phase-vcf",
                str(genotype_path),
                "--flare-vcf",
                str(extra_flare_sample),
                "--n-ancestries",
                str(N_ANCESTRIES),
                "--make-felixla",
                "--out",
                str(work / "bad.extra_flare_sample.out"),
            ],
            expect_fail=True,
            contains="FLARE sample count mismatch",
        )

        bad_ref_unknown = work / "bad.ref.unknown.pvar"
        first_key = next(iter(selected))
        bad_ref_unknown.write_text(f"#CHROM\tPOS\tID\tREF\tALT\n{first_key[0]}\t{first_key[1]}\t.\t.\t{first_key[3]}\n")
        run([*build_bad_base, str(work / "bad.ref.unknown.out"), "--extract", str(bad_ref_unknown)], expect_fail=True, contains="requires known REF")

        bad_ref_mismatch = work / "bad.ref.mismatch.pvar"
        wrong_ref = next(a for a in ALLELES if a != first_key[2])
        bad_ref_mismatch.write_text(f"#CHROM\tPOS\tID\tREF\tALT\n{first_key[0]}\t{first_key[1]}\t.\t{wrong_ref}\t{first_key[3]}\n")
        run([*build_bad_base, str(work / "bad.ref.mismatch.out"), "--extract", str(bad_ref_mismatch)], expect_fail=True, contains="REF mismatch")

        bad_ref_conflict = work / "bad.ref.conflict.pvar"
        bad_ref_conflict.write_text(
            f"#CHROM\tPOS\tID\tREF\tALT\n"
            f"{first_key[0]}\t{first_key[1]}\t.\t{first_key[2]}\t{first_key[3]}\n"
            f"{first_key[0]}\t{first_key[1]}\t.\t{wrong_ref}\t{first_key[3]}\n"
        )
        run([*build_bad_base, str(work / "bad.ref.conflict.out"), "--extract", str(bad_ref_conflict)], expect_fail=True, contains="conflicting REF")

        bad_dup_allele = work / "bad.dup.allele.pvar"
        bad_dup_allele.write_text(
            f"#CHROM\tPOS\tID\tREF\tALT\n"
            f"{first_key[0]}\t{first_key[1]}\tfirst_id\t{first_key[2]}\t{first_key[3]}\n"
            f"{first_key[0]}\t{first_key[1]}\tsecond_id\t{first_key[2]}\t{first_key[3]}\n"
        )
        run([*build_bad_base, str(work / "bad.dup.allele.out"), "--extract", str(bad_dup_allele)], expect_fail=True, contains="duplicate allele")

        bad_bed_empty = work / "bad.empty.bed"
        bad_bed_empty.write_text("# no intervals\ntrack name=empty\n")
        run(
            [*build_bad_base, str(work / "bad.empty.bed.out"), "--extract-bed", str(bad_bed_empty)],
            expect_fail=True,
            contains="interval list is empty",
        )

        bad_bed_short = work / "bad.short.bed"
        bad_bed_short.write_text("chr1\t10\n")
        run(
            [*build_bad_base, str(work / "bad.short.bed.out"), "--extract-bed", str(bad_bed_short)],
            expect_fail=True,
            contains="expects at least CHROM START END",
        )

        bad_bed_negative = work / "bad.negative.bed"
        bad_bed_negative.write_text("chr1\t-1\t10\n")
        run(
            [*build_bad_base, str(work / "bad.negative.bed.out"), "--extract-bed", str(bad_bed_negative)],
            expect_fail=True,
            contains="invalid --extract-bed START",
        )

        bad_bed_zero = work / "bad.zero_length.bed"
        bad_bed_zero.write_text("chr1\t10\t10\n")
        run(
            [*build_bad_base, str(work / "bad.zero_length.bed.out"), "--extract-bed", str(bad_bed_zero)],
            expect_fail=True,
            contains="must have positive length",
        )

        absent_bed = work / "absent_contig.bed"
        absent_bed.write_text("chrAbsent\t0\t100\n")
        absent_prefix = work / "subset.absent_bed"
        absent_vcf = build_with_cli(
            bin_dir,
            genotype_path,
            flare_path,
            absent_prefix,
            "--extract-bed",
            str(absent_bed),
        )
        assert_vcf_matches(absent_vcf, samples, [])
        absent_meta = read_meta(absent_prefix)
        assert absent_meta["global_variants"] == "0", absent_meta
        assert absent_meta["ancestry_blocks"] == "0", absent_meta

        run(
            [
                str(bin_dir / "felixla"),
                "--felixla",
                str(cli_prefix),
                "--export",
                "vcf",
                "--keep",
                str(keep_path),
                "--out",
                str(work / "bad.unsupported.vcf.gz"),
            ],
            expect_fail=True,
            contains="supported only for --phase-vcf + --flare-vcf --make-felixla",
        )

        run(
            [
                str(bin_dir / "felixla"),
                "--felixla",
                str(cli_prefix),
                "--export",
                "vcf",
                "--extract-bed",
                str(bed_path),
                "--out",
                str(work / "bad.unsupported_bed.vcf.gz"),
            ],
            expect_fail=True,
            contains="supported only for --phase-vcf + --flare-vcf --make-felixla",
        )

        print(f"intense keep/extract/BED regression passed in {work}")
    finally:
        if owned_tmp and not args.keep_work:
            owned_tmp.cleanup()

    return 0


if __name__ == "__main__":
    sys.exit(main())
