#!/usr/bin/env python3
"""Intense FELIXla concatenation and corruption-detection regression tests."""

from __future__ import annotations

import argparse
import pathlib
import resource
import shutil
import struct
import subprocess
import sys
import tempfile

from run_keep_extract_intense import build_inputs, read_meta


SUFFIXES = [
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
    ".meta",
]


def fail(message: str) -> None:
    raise AssertionError(message)


def run(cmd: list[str], *, expect_fail: bool = False) -> subprocess.CompletedProcess:
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if expect_fail:
        if result.returncode == 0:
            fail(f"expected failure but command succeeded: {' '.join(cmd)}")
        return result
    if result.returncode != 0:
        fail(
            f"command failed: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def build_prefix(
    felixla: pathlib.Path,
    genotype: pathlib.Path,
    flare: pathlib.Path,
    out: pathlib.Path,
    region: str | None = None,
) -> None:
    cmd = [
        str(felixla),
        "--phase-vcf",
        str(genotype),
        "--flare-vcf",
        str(flare),
        "--n-ancestries",
        "3",
        "--make-felixla",
        "--out",
        str(out),
    ]
    if region:
        cmd.extend(["--region", region])
    run(cmd)


def write_list(path: pathlib.Path, prefixes: list[pathlib.Path], *, meta_first: bool = False) -> None:
    rows = ["# FELIXla chunks", ""]
    for i, prefix in enumerate(prefixes):
        rows.append(str(prefix) + (".meta" if meta_first and i == 0 else ""))
    path.write_text("\n".join(rows) + "\n")


def write_bed_list(
    path: pathlib.Path,
    rows: list[tuple[pathlib.Path, pathlib.Path, bool]],
) -> None:
    lines = ["# FELIXla prefix-tab-BED merge", ""]
    for prefix, bed, complement in rows:
        lines.append(f"{prefix}\t{'^' if complement else ''}{bed}")
    path.write_text("\n".join(lines) + "\n")


def compare_binary_prefixes(expected: pathlib.Path, observed: pathlib.Path) -> None:
    for suffix in SUFFIXES:
        if suffix == ".meta":
            continue
        expected_bytes = pathlib.Path(str(expected) + suffix).read_bytes()
        observed_bytes = pathlib.Path(str(observed) + suffix).read_bytes()
        if expected_bytes != observed_bytes:
            fail(f"BED merge changed binary component {suffix}")


def clone_prefix(source: pathlib.Path, target: pathlib.Path) -> None:
    for suffix in SUFFIXES:
        shutil.copyfile(str(source) + suffix, str(target) + suffix)


def assert_no_output(prefix: pathlib.Path) -> None:
    present = [suffix for suffix in SUFFIXES if pathlib.Path(str(prefix) + suffix).exists()]
    if present:
        fail(f"failed concatenation published output files for {prefix}: {present}")
    temporary = list(prefix.parent.glob(prefix.name + ".concat.tmp.*"))
    if temporary:
        fail(f"failed concatenation left temporary files: {temporary}")


def concat_plink(felixla: pathlib.Path, list_path: pathlib.Path, out: pathlib.Path, *, fail_ok=False):
    result = run(
        [
            str(felixla),
            "--pmerge-list",
            str(list_path),
            "--make-felixla",
            "--out",
            str(out),
        ],
        expect_fail=fail_ok,
    )
    if fail_ok:
        assert_no_output(out)
    return result


def export_prefix(felixla: pathlib.Path, prefix: pathlib.Path, out: pathlib.Path) -> pathlib.Path:
    run([str(felixla), "--felixla", str(prefix), "--export", "vcf", "--out", str(out)])
    return pathlib.Path(str(out) + ".vcf.gz")


def compare_queries(
    felixla: pathlib.Path,
    direct: pathlib.Path,
    merged: pathlib.Path,
    global_count: int,
    work: pathlib.Path,
) -> None:
    indices = sorted({0, 1, 2, global_count // 3, global_count // 2, global_count - 1})
    for index in indices:
        direct_out = work / f"query.direct.{index}.tsv"
        merged_out = work / f"query.merged.{index}.tsv"
        run(
            [
                str(felixla),
                "--felixla",
                str(direct),
                "--global-index",
                str(index),
                "--out",
                str(direct_out),
            ]
        )
        run(
            [
                str(felixla),
                "--felixla",
                str(merged),
                "--global-index",
                str(index),
                "--out",
                str(merged_out),
            ]
        )
        if direct_out.read_bytes() != merged_out.read_bytes():
            fail(f"ancestry-specific query mismatch at global index {index}")


def corrupt_binary_tail(path: pathlib.Path) -> None:
    data = path.read_bytes()
    if not data:
        fail(f"test requires nonempty payload: {path}")
    path.write_bytes(data[:-1])


def expect_corrupt_failure(
    felixla: pathlib.Path,
    source: pathlib.Path,
    next_chunk: pathlib.Path,
    work: pathlib.Path,
    tag: str,
    mutate,
) -> None:
    corrupt = work / f"corrupt_{tag}_prefix"
    clone_prefix(source, corrupt)
    mutate(corrupt)
    list_path = work / f"corrupt.{tag}.list"
    write_list(list_path, [corrupt, next_chunk])
    concat_plink(felixla, list_path, work / f"corrupt.{tag}.out", fail_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", type=pathlib.Path, default=pathlib.Path("bin"))
    parser.add_argument("--work-dir", type=pathlib.Path)
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args()

    bin_dir = args.bin_dir.resolve()
    felixla = bin_dir / "felixla"
    if not felixla.exists():
        fail(f"missing FELIXla binary: {felixla}")

    owned_tmp = None
    if args.work_dir:
        work = args.work_dir.resolve()
        work.mkdir(parents=True, exist_ok=True)
    else:
        owned_tmp = tempfile.TemporaryDirectory(prefix="felixla-concat-intense.")
        work = pathlib.Path(owned_tmp.name)

    try:
        _samples, _records, _ancestry, genotype, flare = build_inputs(work)
        direct = work / "direct"
        build_prefix(felixla, genotype, flare, direct)

        source_a = work / "bed.source.a"
        source_b = work / "bed.source.b"
        clone_prefix(direct, source_a)
        clone_prefix(direct, source_b)

        include_a = work / "include.a.bed"
        include_a.write_text(
            "track name=source_a\n"
            "chr2\t0\t1000\tall_chr2\n"
            "chr1\t750\t1000\tright_2\n"
            "chr1\t0\t250\tleft\n"
            "chr1\t500\t750\tright_1\n"
        )
        include_b = work / "include.b.bed"
        include_b.write_text("chr1\t250\t500\tmiddle\n")
        interleaved_list = work / "bed.interleaved.list"
        # Deliberately list the middle source first; output must still be genomic.
        write_bed_list(
            interleaved_list,
            [(source_b, include_b, False), (source_a, include_a, False)],
        )
        interleaved = work / "bed.interleaved.merged"
        interleaved_result = concat_plink(felixla, interleaved_list, interleaved)
        if "exact sample IDs/order verified" not in interleaved_result.stderr:
            fail(f"BED concat did not report exact sample verification:\n{interleaved_result.stderr}")
        compare_binary_prefixes(direct, interleaved)
        interleaved_meta = read_meta(interleaved)
        if interleaved_meta.get("concat_bed_coordinates") != "0-based-half-open":
            fail(f"BED coordinate provenance is missing: {interleaved_meta}")
        if interleaved_meta.get("concat_bed_overlap_check") != "effective-selections-disjoint":
            fail(f"BED overlap provenance is missing: {interleaved_meta}")

        exclude_middle = work / "exclude.middle.bed"
        exclude_middle.write_text("chr1\t250\t500\n")
        complement_list = work / "bed.complement.list"
        write_bed_list(
            complement_list,
            [(source_a, exclude_middle, True), (source_b, include_b, False)],
        )
        complement_merged = work / "bed.complement.merged"
        concat_plink(felixla, complement_list, complement_merged)
        compare_binary_prefixes(direct, complement_merged)
        complement_meta = read_meta(complement_merged)
        if complement_meta.get("concat_bed_1") != f"^{exclude_middle}":
            fail(f"complement BED provenance is wrong: {complement_meta}")

        striped_a = work / "striped.a.bed"
        striped_b = work / "striped.b.bed"
        striped_rows = {striped_a: [], striped_b: []}
        for chrom in ["chr2", "chr1"]:
            for start0 in range(0, 1000, 50):
                owner = striped_a if (start0 // 50) % 2 == 0 else striped_b
                striped_rows[owner].append(f"{chrom}\t{start0}\t{start0 + 50}")
        for path, rows in striped_rows.items():
            path.write_text("\n".join(reversed(rows)) + "\n")
        striped_list = work / "bed.striped.list"
        write_bed_list(
            striped_list,
            [(source_b, striped_b, False), (source_a, striped_a, False)],
        )
        striped = work / "bed.striped.merged"
        concat_plink(felixla, striped_list, striped)
        striped_meta = read_meta(striped)
        direct_meta_for_bed = read_meta(direct)
        for key in ["global_variants", "common_variants", "rare_variants"]:
            if striped_meta[key] != direct_meta_for_bed[key]:
                fail(f"striped BED merge changed {key}: {striped_meta}")
        striped_vcf = export_prefix(felixla, striped, work / "bed.striped.export")
        direct_bed_vcf = export_prefix(felixla, direct, work / "bed.direct.export")
        striped_comparison = run(
            [str(felixla), "--compare-vcfs", str(direct_bed_vcf), str(striped_vcf)]
        )
        if "Differences:             0" not in striped_comparison.stdout:
            fail(f"striped BED merge changed exported genotypes:\n{striped_comparison.stdout}")
        compare_queries(
            felixla,
            direct,
            striped,
            int(direct_meta_for_bed["global_variants"]),
            work,
        )

        many_sources: list[pathlib.Path] = []
        many_beds: list[pathlib.Path] = []
        many_rows: list[list[str]] = [[] for _ in range(10)]
        for chrom in ["chr1", "chr2"]:
            for segment, start0 in enumerate(range(0, 1000, 50)):
                owner = segment % len(many_rows)
                many_rows[owner].append(f"{chrom}\t{start0}\t{start0 + 50}")
        for index, rows in enumerate(many_rows):
            source = work / f"bed.many.source.{index:02d}"
            bed = work / f"bed.many.{index:02d}.bed"
            clone_prefix(direct, source)
            bed.write_text("\n".join(reversed(rows)) + "\n")
            many_sources.append(source)
            many_beds.append(bed)
        many_list = work / "bed.many.list"
        write_bed_list(
            many_list,
            list(reversed([
                (source, bed, False)
                for source, bed in zip(many_sources, many_beds)
            ])),
        )
        original_nofile = resource.getrlimit(resource.RLIMIT_NOFILE)
        lowered_nofile = (
            64
            if original_nofile[0] == resource.RLIM_INFINITY
            else min(original_nofile[0], 64)
        )
        can_raise = original_nofile[1] == resource.RLIM_INFINITY or original_nofile[1] >= 76
        if can_raise:
            resource.setrlimit(resource.RLIMIT_NOFILE, (lowered_nofile, original_nofile[1]))
        try:
            many = work / "bed.many.merged"
            many_result = concat_plink(felixla, many_list, many)
        finally:
            if can_raise:
                resource.setrlimit(resource.RLIMIT_NOFILE, original_nofile)
        if can_raise and "raised open-file soft limit" not in many_result.stderr:
            fail(f"many-input BED merge did not raise its file limit:\n{many_result.stderr}")
        many_vcf = export_prefix(felixla, many, work / "bed.many.export")
        many_comparison = run(
            [str(felixla), "--compare-vcfs", str(direct_bed_vcf), str(many_vcf)]
        )
        if "Differences:             0" not in many_comparison.stdout:
            fail(f"10-way BED merge changed exported genotypes:\n{many_comparison.stdout}")

        overlap_a = work / "overlap.a.bed"
        overlap_b = work / "overlap.b.bed"
        overlap_a.write_text("chr1\t0\t500\n")
        overlap_b.write_text("chr1\t250\t750\n")
        bed_overlap_list = work / "bed.overlap.list"
        write_bed_list(
            bed_overlap_list,
            [(source_a, overlap_a, False), (source_b, overlap_b, False)],
        )
        bed_overlap_result = concat_plink(
            felixla, bed_overlap_list, work / "bed.overlap.out", fail_ok=True
        )
        overlap_report = bed_overlap_result.stdout + bed_overlap_result.stderr
        if "effective BED selection overlaps" not in overlap_report:
            fail(f"overlapping BED selections were not diagnosed:\n{overlap_report}")
        for source in [source_a, source_b]:
            if str(source) not in overlap_report:
                fail(f"BED overlap report omitted {source}:\n{overlap_report}")

        duplicate_bed = work / "duplicate.within.bed"
        duplicate_bed.write_text("chr1\t0\t300\nchr1\t200\t400\n")
        duplicate_bed_list = work / "bed.duplicate.list"
        write_bed_list(
            duplicate_bed_list,
            [(source_a, duplicate_bed, False), (source_b, include_b, False)],
        )
        duplicate_result = concat_plink(
            felixla, duplicate_bed_list, work / "bed.duplicate.out", fail_ok=True
        )
        if "overlaps or duplicates an earlier BED interval" not in (
            duplicate_result.stdout + duplicate_result.stderr
        ):
            fail("within-BED overlap was not rejected during preflight")

        unselected_corrupt = work / "bed.unselected.corrupt"
        clone_prefix(direct, unselected_corrupt)
        corrupt_binary_tail(pathlib.Path(str(unselected_corrupt) + ".ancblock.bin"))
        first_quarter = work / "first.quarter.bed"
        first_quarter.write_text("chr1\t0\t250\n")
        unselected_corrupt_list = work / "bed.unselected.corrupt.list"
        write_bed_list(
            unselected_corrupt_list,
            [
                (unselected_corrupt, first_quarter, False),
                (source_b, include_b, False),
            ],
        )
        unselected_corrupt_result = concat_plink(
            felixla,
            unselected_corrupt_list,
            work / "bed.unselected.corrupt.out",
            fail_ok=True,
        )
        if "unexpected EOF" not in (
            unselected_corrupt_result.stdout + unselected_corrupt_result.stderr
        ):
            fail("BED mode did not validate a corrupted payload outside its selected region")

        missing_bed = work / "does.not.exist.bed"
        malformed_bed = work / "malformed.bed"
        malformed_bed.write_text("chr1\t100\n")
        bad_beds_list = work / "bed.multiple_bad.list"
        write_bed_list(
            bad_beds_list,
            [(source_a, missing_bed, False), (source_b, malformed_bed, False)],
        )
        bad_beds_result = concat_plink(
            felixla, bad_beds_list, work / "bed.multiple_bad.out", fail_ok=True
        )
        bad_beds_report = bad_beds_result.stdout + bad_beds_result.stderr
        if "found 2 problematic prefix(es)" not in bad_beds_report:
            fail(f"multiple bad BEDs were not aggregated:\n{bad_beds_report}")
        for expected in ["cannot open BED", "must contain at least three BED columns"]:
            if expected not in bad_beds_report:
                fail(f"bad BED aggregate report omitted {expected!r}:\n{bad_beds_report}")

        sample_bad = work / "bed.sample.bad"
        clone_prefix(source_b, sample_bad)
        sample_bad_path = pathlib.Path(str(sample_bad) + ".samples")
        sample_bad_rows = sample_bad_path.read_text().splitlines()
        sample_bad_rows[0], sample_bad_rows[1] = sample_bad_rows[1], sample_bad_rows[0]
        sample_bad_path.write_text("\n".join(sample_bad_rows) + "\n")
        sample_bad_list = work / "bed.sample.bad.list"
        write_bed_list(
            sample_bad_list,
            [(source_a, include_a, False), (sample_bad, include_b, False)],
        )
        sample_bad_result = concat_plink(
            felixla, sample_bad_list, work / "bed.sample.bad.out", fail_ok=True
        )
        if "sample IDs/order differ" not in (sample_bad_result.stdout + sample_bad_result.stderr):
            fail("BED concat accepted mismatched sample order")

        mixed_list = work / "bed.mixed_columns.list"
        mixed_list.write_text(f"{source_a}\t{include_a}\n{source_b}\n")
        mixed_result = concat_plink(
            felixla, mixed_list, work / "bed.mixed_columns.out", fail_ok=True
        )
        if "mixes one-column and prefix-tab-BED rows" not in (
            mixed_result.stdout + mixed_result.stderr
        ):
            fail("mixed one-column/two-column merge list was not rejected")

        regions = [
            "chr1:1-250",
            "chr1:251-500",
            "chr1:501-750",
            "chr1:751-1000",
            "chr2:1-230",
            "chr2:231-520",
            "chr2:521-770",
            "chr2:771-1000",
            "chr2:1001-1100",  # A complete zero-variant/zero-block chunk is valid.
        ]
        chunks = []
        for index, region in enumerate(regions, start=1):
            prefix = work / f"chunk.{index:02d}"
            build_prefix(felixla, genotype, flare, prefix, region)
            chunks.append(prefix)

        list_path = work / "chunks.list"
        write_list(list_path, chunks, meta_first=True)
        merged = work / "merged"
        concat_plink(felixla, list_path, merged)

        direct_meta = read_meta(direct)
        merged_meta = read_meta(merged)
        for key in ["global_variants", "common_variants", "rare_variants"]:
            if merged_meta[key] != direct_meta[key]:
                fail(f"merged {key} mismatch: {merged_meta[key]} != {direct_meta[key]}")
        if merged_meta["concat_inputs"] != str(len(chunks)):
            fail(f"wrong concat_inputs: {merged_meta}")
        if "magic,index,offset,eof,mac,masks,ordering" not in merged_meta["integrity_checks"]:
            fail(f"missing integrity manifest: {merged_meta}")
        if merged_meta.get("concat_validation") != "full-preflight-then-merge":
            fail(f"missing two-phase validation manifest: {merged_meta}")

        direct_vcf = export_prefix(felixla, direct, work / "direct.export")
        merged_vcf = export_prefix(felixla, merged, work / "merged.export")
        comparison = run([str(felixla), "--compare-vcfs", str(direct_vcf), str(merged_vcf)])
        if "Differences:             0" not in comparison.stdout:
            fail(f"direct/merged VCF mismatch:\n{comparison.stdout}")
        compare_queries(
            felixla,
            direct,
            merged,
            int(direct_meta["global_variants"]),
            work,
        )

        one_list = work / "one.list"
        write_list(one_list, [merged])
        compat = work / "compat"
        run([str(felixla), "concat", str(one_list), str(compat)])
        if read_meta(compat)["global_variants"] != direct_meta["global_variants"]:
            fail("compatibility-mode concat changed global variant count")

        legacy = work / "legacy.no_counts"
        clone_prefix(chunks[0], legacy)
        legacy_meta_path = pathlib.Path(str(legacy) + ".meta")
        count_keys = {"global_variants", "common_variants", "rare_variants", "ancestry_blocks"}
        legacy_rows = [
            line
            for line in legacy_meta_path.read_text().splitlines()
            if line.split("\t", 1)[0] not in count_keys
        ]
        legacy_meta_path.write_text("\n".join(legacy_rows) + "\n")
        legacy_list = work / "legacy.list"
        write_list(legacy_list, [legacy, chunks[1]])
        concat_plink(felixla, legacy_list, work / "legacy.merged")

        extracted = work / "merged.extracted"
        run(
            [
                str(felixla),
                "--felixla",
                str(merged),
                "--region",
                "chr1:1-500",
                "--make-felixla",
                "--out",
                str(extracted),
            ]
        )
        extracted_meta = read_meta(extracted)
        if int(extracted_meta["global_variants"]) <= 0:
            fail(f"extracted record manifest was not refreshed: {extracted_meta}")
        extracted_list = work / "extracted.list"
        write_list(extracted_list, [extracted])
        concat_plink(felixla, extracted_list, work / "extracted.repacked")

        first = chunks[0]
        second = chunks[1]
        for suffix in SUFFIXES:
            def mutate(prefix: pathlib.Path, suffix=suffix) -> None:
                path = pathlib.Path(str(prefix) + suffix)
                if suffix == ".meta":
                    rows = [line for line in path.read_text().splitlines() if not line.startswith("n_words\t")]
                    path.write_text("\n".join(rows) + "\n")
                elif suffix == ".samples":
                    rows = path.read_text().splitlines()
                    rows[0], rows[1] = rows[1], rows[0]
                    path.write_text("\n".join(rows) + "\n")
                else:
                    corrupt_binary_tail(path)

            expect_corrupt_failure(
                felixla,
                first,
                second,
                work,
                suffix.lstrip(".").replace(".", "_"),
                mutate,
            )

        def corrupt_common_bits(prefix: pathlib.Path) -> None:
            path = pathlib.Path(str(prefix) + ".common.geno.bin")
            data = bytearray(path.read_bytes())
            data[0] ^= 1
            path.write_bytes(data)

        expect_corrupt_failure(felixla, first, second, work, "common_popcount", corrupt_common_bits)

        def corrupt_rare_global(prefix: pathlib.Path) -> None:
            path = pathlib.Path(str(prefix) + ".rare.carrier.bin")
            data = bytearray(path.read_bytes())
            struct.pack_into("<I", data, 0, 0xFFFFFFFF)
            path.write_bytes(data)

        expect_corrupt_failure(felixla, first, second, work, "rare_global", corrupt_rare_global)

        def corrupt_ancestry_partition(prefix: pathlib.Path) -> None:
            path = pathlib.Path(str(prefix) + ".ancblock.bin")
            data = bytearray(path.read_bytes())
            n_words = int(read_meta(prefix)["n_words"])
            first_mask = struct.unpack_from("<Q", data, 0)[0]
            bit = first_mask & -first_mask
            if bit == 0:
                fail("test fixture ancestry 0 mask unexpectedly empty")
            second_offset = n_words * 8
            second_mask = struct.unpack_from("<Q", data, second_offset)[0]
            struct.pack_into("<Q", data, second_offset, second_mask | bit)
            path.write_bytes(data)

        expect_corrupt_failure(
            felixla, first, second, work, "ancestry_partition", corrupt_ancestry_partition
        )

        def corrupt_index_offset(prefix: pathlib.Path) -> None:
            path = pathlib.Path(str(prefix) + ".common.variant.idx")
            data = bytearray(path.read_bytes())
            marker_offset = struct.unpack_from("<Q", data, 20)[0]
            struct.pack_into("<Q", data, 20, marker_offset + 1)
            path.write_bytes(data)

        expect_corrupt_failure(felixla, first, second, work, "index_offset", corrupt_index_offset)

        def corrupt_meta_count(prefix: pathlib.Path) -> None:
            path = pathlib.Path(str(prefix) + ".meta")
            rows = []
            for line in path.read_text().splitlines():
                if line.startswith("global_variants\t"):
                    key, value = line.split("\t")
                    line = f"{key}\t{int(value) + 1}"
                rows.append(line)
            path.write_text("\n".join(rows) + "\n")

        expect_corrupt_failure(felixla, first, second, work, "meta_count", corrupt_meta_count)

        overlap = work / "overlap.copy"
        clone_prefix(first, overlap)
        overlap_list = work / "overlap.list"
        write_list(overlap_list, [first, overlap])
        concat_plink(felixla, overlap_list, work / "overlap.out", fail_ok=True)

        missing = work / "missing.copy"
        clone_prefix(first, missing)
        pathlib.Path(str(missing) + ".ancblock.idx").unlink()
        missing_list = work / "missing.list"
        write_list(missing_list, [missing, second])
        concat_plink(felixla, missing_list, work / "missing.out", fail_ok=True)

        multi_bad_meta = work / "multi.badmeta"
        clone_prefix(chunks[1], multi_bad_meta)
        multi_meta_path = pathlib.Path(str(multi_bad_meta) + ".meta")
        multi_meta_path.write_text(
            "\n".join(
                line
                for line in multi_meta_path.read_text().splitlines()
                if not line.startswith("n_words\t")
            )
            + "\n"
        )

        multi_bad_samples = work / "multi.bad.samples"
        clone_prefix(chunks[2], multi_bad_samples)
        multi_samples_path = pathlib.Path(str(multi_bad_samples) + ".samples")
        multi_sample_rows = multi_samples_path.read_text().splitlines()
        multi_sample_rows[0], multi_sample_rows[1] = multi_sample_rows[1], multi_sample_rows[0]
        multi_samples_path.write_text("\n".join(multi_sample_rows) + "\n")

        multi_bad_binary = work / "multi.bad.binary"
        clone_prefix(chunks[3], multi_bad_binary)
        corrupt_binary_tail(pathlib.Path(str(multi_bad_binary) + ".ancblock.idx"))

        multi_bad_list = work / "multi.bad.list"
        write_list(
            multi_bad_list,
            [chunks[0], multi_bad_meta, multi_bad_samples, multi_bad_binary],
        )
        multi_result = concat_plink(
            felixla,
            multi_bad_list,
            work / "multi.bad.out",
            fail_ok=True,
        )
        multi_report = multi_result.stdout + multi_result.stderr
        for bad_prefix in [multi_bad_meta, multi_bad_samples, multi_bad_binary]:
            if str(bad_prefix) not in multi_report:
                fail(f"aggregate preflight report omitted {bad_prefix}:\n{multi_report}")
        if "found 3 problematic prefix(es)" not in multi_report or "merge was not started" not in multi_report:
            fail(f"aggregate preflight did not report all failures:\n{multi_report}")
        for expected_error in ["missing n_words", "sample IDs/order differ", "truncated ancestry index"]:
            if expected_error not in multi_report:
                fail(f"aggregate preflight omitted {expected_error!r}:\n{multi_report}")

        print(f"intense concat regression passed in {work}")
    finally:
        if owned_tmp and not args.keep_work:
            owned_tmp.cleanup()
    return 0


if __name__ == "__main__":
    sys.exit(main())
