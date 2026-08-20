#!/usr/bin/env python3
"""Generate a large-sample FELIXla packing benchmark and report throughput."""

from __future__ import annotations

import argparse
import gzip
import json
import pathlib
import shutil
import subprocess
import tempfile
import time


def write_inputs(
    work: pathlib.Path,
    n_samples: int,
    n_records: int,
    flare_every_record: bool,
    flare_change_every: int,
    multiallelic_every: int,
    drop_flare_every: int,
) -> tuple[pathlib.Path, pathlib.Path, pathlib.Path]:
    samples = [f"s{i:07d}" for i in range(n_samples)]
    genotype = work / "benchmark.phased.vcf"
    flare = work / "benchmark.flare.vcf"
    extract = work / "benchmark.scrambled.pvar"

    common_gts = ["0|1" if i % 2 == 0 else "1|0" for i in range(n_samples)]
    rare_carriers = max(1, n_samples // 128)
    rare_gts = ["0|1" if i < rare_carriers else "0|0" for i in range(n_samples)]
    multiallelic_gts = [
        ("1|2", "2|0", "0|1", "0|0")[i % 4]
        for i in range(n_samples)
    ]
    flare_sample_indices = [
        i
        for i in range(n_samples)
        if drop_flare_every <= 0 or (i + 1) % drop_flare_every != 0
    ]

    with genotype.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples) + "\n")
        for record_i in range(n_records):
            pos = 100 + record_i
            is_multiallelic = multiallelic_every > 0 and record_i % multiallelic_every == 0
            if is_multiallelic:
                alts = "C,G"
                gts = multiallelic_gts
            else:
                alts = "C"
                gts = rare_gts if record_i % 4 == 0 else common_gts
            out.write(f"chr1\t{pos}\tv{record_i}\tA\t{alts}\t.\tPASS\t.\tGT\t")
            out.write("\t".join(gts) + "\n")

    with flare.open("w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##contig=<ID=chr1>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="First ancestry">\n')
        out.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Second ancestry">\n')
        out.write('##FORMAT=<ID=ANP1,Number=1,Type=Float,Description="First ancestry probability">\n')
        out.write('##FORMAT=<ID=ANP2,Number=1,Type=Float,Description="Second ancestry probability">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        out.write("\t".join(samples[i] for i in flare_sample_indices) + "\n")
        flare_record_indices = range(n_records) if flare_every_record else [n_records - 1]
        for record_i in flare_record_indices:
            ancestry_state = record_i // flare_change_every
            ancestry = [
                f"0|0:{(i + ancestry_state) % 3}:{(i + ancestry_state + 1) % 3}:0.99:0.98"
                for i in flare_sample_indices
            ]
            out.write(
                f"chr1\t{100 + record_i}\t.\tA\tC\t.\tPASS\t.\t"
                "GT:AN1:AN2:ANP1:ANP2\t"
            )
            out.write("\t".join(ancestry) + "\n")

    with extract.open("w") as out:
        out.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for record_i in reversed(range(n_records)):
            alts = (
                "C,G"
                if multiallelic_every > 0 and record_i % multiallelic_every == 0
                else "C"
            )
            out.write(f"chr1\t{100 + record_i}\tignored_{record_i}\tA\t{alts}\n")

    return genotype, flare, extract


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--felixla", type=pathlib.Path, default=pathlib.Path("bin/felixla"))
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--records", type=int, default=500)
    parser.add_argument("--flare-every-record", action="store_true")
    parser.add_argument(
        "--flare-change-every",
        type=int,
        default=1,
        metavar="N",
        help="change simulated ancestry labels every N FLARE records",
    )
    parser.add_argument(
        "--multiallelic-every",
        type=int,
        default=0,
        metavar="N",
        help="emit a three-allele genotype record every N records (0 disables)",
    )
    parser.add_argument(
        "--drop-flare-every",
        type=int,
        default=0,
        metavar="N",
        help="omit every Nth genotype sample from FLARE to benchmark intersection",
    )
    parser.add_argument("--gzip-inputs", action="store_true")
    parser.add_argument("--work-dir", type=pathlib.Path)
    args = parser.parse_args()

    if (
        args.samples <= 0
        or args.records <= 0
        or args.flare_change_every <= 0
        or args.multiallelic_every < 0
        or args.drop_flare_every < 0
    ):
        parser.error("sample/record counts and FLARE change interval must be positive")

    owned_tmp = None
    if args.work_dir:
        work = args.work_dir
        work.mkdir(parents=True, exist_ok=True)
    else:
        owned_tmp = tempfile.TemporaryDirectory(prefix="felixla-pack-benchmark.")
        work = pathlib.Path(owned_tmp.name)

    genotype, flare, extract = write_inputs(
        work,
        args.samples,
        args.records,
        args.flare_every_record,
        args.flare_change_every,
        args.multiallelic_every,
        args.drop_flare_every,
    )
    if args.gzip_inputs:
        compressed_paths = []
        for source in (genotype, flare):
            destination = source.with_suffix(source.suffix + ".gz")
            with source.open("rb") as src, gzip.open(destination, "wb") as dst:
                shutil.copyfileobj(src, dst)
            compressed_paths.append(destination)
        genotype, flare = compressed_paths
    prefix = work / "packed"
    command = [
        str(args.felixla.resolve()),
        "--phase-vcf",
        str(genotype),
        "--flare-vcf",
        str(flare),
        "--n-ancestries",
        "3",
        "--extract",
        str(extract),
        "--make-felixla",
        "--out",
        str(prefix),
    ]

    started = time.perf_counter()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed = time.perf_counter() - started
    if result.returncode != 0:
        raise RuntimeError(f"benchmark failed\n{result.stdout}\n{result.stderr}")

    sample_records = args.samples * args.records
    report = {
        "samples": args.samples,
        "records": args.records,
        "sample_records": sample_records,
        "flare_records": args.records if args.flare_every_record else 1,
        "flare_change_every": args.flare_change_every,
        "multiallelic_every": args.multiallelic_every,
        "drop_flare_every": args.drop_flare_every,
        "gzip_inputs": args.gzip_inputs,
        "elapsed_seconds": round(elapsed, 6),
        "records_per_second": round(args.records / elapsed, 3),
        "sample_records_per_second": round(sample_records / elapsed, 3),
    }
    print(json.dumps(report, sort_keys=True))
    if owned_tmp is not None:
        owned_tmp.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
