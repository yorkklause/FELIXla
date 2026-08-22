#!/usr/bin/env python3
"""Intense FELIXla concatenation and corruption-detection regression tests."""

from __future__ import annotations

import argparse
import pathlib
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
