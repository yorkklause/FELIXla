#!/usr/bin/env python3
# designed by Kai, implemented by codex
#
# Generate 1-based inclusive chunk manifests for tractor_hybrid conversion
# wrapper scripts.

import argparse
import gzip
import re
import sys
from collections import OrderedDict
from pathlib import Path


GRCH38 = OrderedDict(
    [
        ("1", 248_956_422),
        ("2", 242_193_529),
        ("3", 198_295_559),
        ("4", 190_214_555),
        ("5", 181_538_259),
        ("6", 170_805_979),
        ("7", 159_345_973),
        ("8", 145_138_636),
        ("9", 138_394_717),
        ("10", 133_797_422),
        ("11", 135_086_622),
        ("12", 133_275_309),
        ("13", 114_364_328),
        ("14", 107_043_718),
        ("15", 101_991_189),
        ("16", 90_338_345),
        ("17", 83_257_441),
        ("18", 80_373_285),
        ("19", 58_617_616),
        ("20", 64_444_167),
        ("21", 46_709_983),
        ("22", 50_818_468),
    ]
)

GRCH37 = OrderedDict(
    [
        ("1", 249_250_621),
        ("2", 243_199_373),
        ("3", 198_022_430),
        ("4", 191_154_276),
        ("5", 180_915_260),
        ("6", 171_115_067),
        ("7", 159_138_663),
        ("8", 146_364_022),
        ("9", 141_213_431),
        ("10", 135_534_747),
        ("11", 135_006_516),
        ("12", 133_851_895),
        ("13", 115_169_878),
        ("14", 107_349_540),
        ("15", 102_531_392),
        ("16", 90_354_753),
        ("17", 81_195_210),
        ("18", 78_077_248),
        ("19", 59_128_983),
        ("20", 63_025_520),
        ("21", 48_129_895),
        ("22", 51_304_566),
    ]
)


CONTIG_RE = re.compile(r"##contig=<(.+)>")


def die(message):
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def open_text(path):
    path = str(path)
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_int(value, label):
    clean = value.replace(",", "").replace("_", "")
    try:
        parsed = int(clean)
    except ValueError:
        die(f"invalid integer for {label}: {value}")
    if parsed <= 0:
        die(f"{label} must be positive: {value}")
    return parsed


def split_contig_payload(payload):
    fields = []
    start = 0
    in_quote = False

    for i, char in enumerate(payload):
        if char == '"':
            in_quote = not in_quote
        elif char == "," and not in_quote:
            fields.append(payload[start:i])
            start = i + 1

    fields.append(payload[start:])
    return fields


def parse_vcf_contig_line(line):
    match = CONTIG_RE.match(line)
    if not match:
        return None

    values = {}
    for field in split_contig_payload(match.group(1)):
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        values[key.strip()] = value.strip().strip('"')

    contig = values.get("ID")
    length = values.get("length")
    if not contig or not length:
        return None

    return contig, parse_int(length, f"{contig} length")


def read_vcf_lengths(paths):
    lengths = OrderedDict()

    for path in paths:
        with open_text(path) as handle:
            for line in handle:
                if line.startswith("##contig=<"):
                    parsed = parse_vcf_contig_line(line.rstrip("\n\r"))
                    if not parsed:
                        continue
                    chrom, length = parsed
                    if chrom in lengths and lengths[chrom] != length:
                        die(f"conflicting contig length for {chrom}: {lengths[chrom]} vs {length}")
                    lengths[chrom] = length
                elif line.startswith("#CHROM"):
                    break

    if not lengths:
        die("no ##contig=<ID=...,length=...> lines found in VCF header(s)")

    return lengths


def read_chrom_sizes(path):
    lengths = OrderedDict()
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                die(f"chrom sizes line {line_no} needs at least 2 columns")
            chrom = fields[0]
            length = parse_int(fields[1], f"{chrom} length")
            if chrom in lengths and lengths[chrom] != length:
                die(f"conflicting chrom size for {chrom}: {lengths[chrom]} vs {length}")
            lengths[chrom] = length
    if not lengths:
        die(f"no chromosome sizes found in {path}")
    return lengths


def parse_chroms(text):
    if text.lower() in {"all", "*"}:
        return None

    chroms = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            left, right = item.split("-", 1)
            left_i = int(left.replace("chr", ""))
            right_i = int(right.replace("chr", ""))
            if left_i > right_i:
                die(f"invalid chromosome range: {item}")
            chroms.extend(str(i) for i in range(left_i, right_i + 1))
        else:
            chroms.append(item)
    if not chroms:
        die("--chroms did not contain any chromosomes")
    return chroms


def normalize_chrom(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def apply_chrom_style(chrom, style):
    base = normalize_chrom(chrom)
    if style == "chr":
        return "chr" + base
    if style == "nochr":
        return base
    return chrom


def select_chroms(lengths, requested, style):
    if requested is None:
        return OrderedDict((apply_chrom_style(chrom, style), length) for chrom, length in lengths.items())

    by_exact = dict(lengths)
    by_norm = {}
    for chrom, length in lengths.items():
        by_norm[normalize_chrom(chrom)] = (chrom, length)

    selected = OrderedDict()
    for wanted in requested:
        source = None
        if wanted in by_exact:
            source = (wanted, by_exact[wanted])
        else:
            norm = normalize_chrom(wanted)
            source = by_norm.get(norm)
        if source is None:
            die(f"requested chromosome {wanted} was not found in available contigs")
        chrom, length = source
        out_chrom = apply_chrom_style(chrom, style)
        selected[out_chrom] = length
    return selected


def built_in_lengths(build, style):
    source = GRCH38 if build == "GRCh38" else GRCH37
    if style == "auto":
        style = "chr"
    return OrderedDict((apply_chrom_style(chrom, style), length) for chrom, length in source.items())


def render_template(template, label, **values):
    try:
        return template.format(**values)
    except KeyError as exc:
        die(f"unknown placeholder in {label}: {exc}")
    except ValueError as exc:
        die(f"invalid {label} format: {exc}")


def template_prefix(template, chrom, chunk, start, end, region):
    return render_template(
        template,
        "--out-prefix-template",
        chrom=chrom,
        chunk=chunk,
        start=start,
        end=end,
        region=region,
    )


def emit_manifest(lengths, args):
    if args.header:
        if args.format == "tsv4":
            print("chrom\tstart\tend\tout_prefix")
        elif args.format == "tsv5":
            print("chrom\tstart\tend\tregion\tout_prefix")
        elif args.format == "regions":
            print("region")
        elif args.format == "flare":
            print(
                "chrom\tstart\tend\tregion\tsource_phase_vcf\tsource_flare_vcf\t"
                "chunk_phase_vcf\tchunk_flare_vcf\tn_ancestries\tmac_threshold\tout_prefix"
            )
        elif args.format == "flare-args":
            print("chunk_phase_vcf\tchunk_flare_vcf\tn_ancestries\tmac_threshold\tout_prefix")
        elif args.format == "flare-direct":
            print(
                "chrom\tstart\tend\tsource_phase_vcf\tsource_flare_vcf\t"
                "n_ancestries\tmac_threshold\tout_prefix\tregion"
            )
        elif args.format == "flare-direct-args":
            print("source_phase_vcf\tsource_flare_vcf\tn_ancestries\tmac_threshold\tout_prefix\tregion")

    chunk_bp = parse_int(args.chunk_bp, "--chunk-bp")
    chunk_index = 0

    for chrom, chrom_length in lengths.items():
        start = 1
        per_chrom_chunk = 0
        while start <= chrom_length:
            end = min(start + chunk_bp - 1, chrom_length)
            chunk_index += 1
            per_chrom_chunk += 1
            region = f"{chrom}:{start}-{end}"
            out_prefix = template_prefix(
                args.out_prefix_template,
                chrom=chrom,
                chunk=per_chrom_chunk if args.chunk_numbering == "per-chrom" else chunk_index,
                start=start,
                end=end,
                region=region,
            )
            template_values = {
                "chrom": chrom,
                "chunk": per_chrom_chunk if args.chunk_numbering == "per-chrom" else chunk_index,
                "start": start,
                "end": end,
                "region": region,
                "out_prefix": out_prefix,
            }

            if args.format == "tsv4":
                print(f"{chrom}\t{start}\t{end}\t{out_prefix}")
            elif args.format == "tsv5":
                print(f"{chrom}\t{start}\t{end}\t{region}\t{out_prefix}")
            elif args.format == "regions":
                print(region)
            elif args.format in {"flare", "flare-args"}:
                chunk_phase = render_template(
                    args.chunk_phase_template,
                    "--chunk-phase-template",
                    **template_values,
                )
                chunk_flare = render_template(
                    args.chunk_flare_template,
                    "--chunk-flare-template",
                    **template_values,
                )
                converter_args = (
                    f"{chunk_phase}\t{chunk_flare}\t"
                    f"{args.n_ancestries}\t{args.mac_threshold}\t{out_prefix}"
                )
                if args.format == "flare":
                    source_phase = render_template(
                        args.phase_template,
                        "--phase-template",
                        **template_values,
                    )
                    source_flare = render_template(
                        args.flare_template,
                        "--flare-template",
                        **template_values,
                    )
                    print(
                        f"{chrom}\t{start}\t{end}\t{region}\t"
                        f"{source_phase}\t{source_flare}\t{converter_args}"
                    )
                else:
                    print(converter_args)
            elif args.format in {"flare-direct", "flare-direct-args"}:
                source_phase = render_template(
                    args.phase_template,
                    "--phase-template",
                    **template_values,
                )
                source_flare = render_template(
                    args.flare_template,
                    "--flare-template",
                    **template_values,
                )
                converter_args = (
                    f"{source_phase}\t{source_flare}\t"
                    f"{args.n_ancestries}\t{args.mac_threshold}\t{out_prefix}\t{region}"
                )
                if args.format == "flare-direct":
                    print(
                        f"{chrom}\t{start}\t{end}\t"
                        f"{source_phase}\t{source_flare}\t"
                        f"{args.n_ancestries}\t{args.mac_threshold}\t{out_prefix}\t{region}"
                    )
                else:
                    print(converter_args)

            start = end + 1


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate chunk manifests with 1-based inclusive coordinates."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--vcf",
        action="append",
        help="VCF path whose ##contig length header lines define chromosome sizes. May be repeated.",
    )
    source.add_argument(
        "--chrom-sizes",
        help="Two-column chrom sizes file: chrom length.",
    )
    source.add_argument(
        "--build",
        choices=["GRCh38", "GRCh37"],
        help="Use built-in chromosome sizes for chromosomes 1-22.",
    )
    parser.add_argument(
        "--chroms",
        default="1-22",
        help="Chromosomes to include, for example 1-22, chr1-chr22, 1,2,22, or all. Default: 1-22.",
    )
    parser.add_argument(
        "--chrom-style",
        choices=["auto", "chr", "nochr"],
        default="auto",
        help="Output chromosome naming. auto keeps input names; built-in builds default to chr*. Default: auto.",
    )
    parser.add_argument(
        "--chunk-bp",
        default="50000000",
        help="Chunk size in bp. Default: 50000000.",
    )
    parser.add_argument(
        "--out-prefix-template",
        default="out/{chrom}.chunk{chunk:04d}",
        help="Python format template for output prefix. Placeholders: {chrom}, {chunk}, {start}, {end}, {region}.",
    )
    parser.add_argument(
        "--phase-template",
        default="phase/{chrom}.phased.vcf.gz",
        help="Source phased VCF path template for --format flare. Supports {chrom}, {chunk}, {start}, {end}, {region}, {out_prefix}.",
    )
    parser.add_argument(
        "--flare-template",
        default="flare/{chrom}.flare.vcf.gz",
        help="Source FLARE VCF path template for --format flare. Supports {chrom}, {chunk}, {start}, {end}, {region}, {out_prefix}.",
    )
    parser.add_argument(
        "--chunk-phase-template",
        default="{out_prefix}.phase.vcf.gz",
        help="Chunked phased VCF path template for --format flare and flare-args. Default: {out_prefix}.phase.vcf.gz.",
    )
    parser.add_argument(
        "--chunk-flare-template",
        default="{out_prefix}.flare.vcf.gz",
        help="Chunked FLARE VCF path template for --format flare and flare-args. Default: {out_prefix}.flare.vcf.gz.",
    )
    parser.add_argument(
        "--n-ancestries",
        help="n_ancestries argument for flare_subset_to_tractor_hybrid when using a flare* format.",
    )
    parser.add_argument(
        "--mac-threshold",
        help="mac_threshold argument for flare_subset_to_tractor_hybrid when using a flare* format.",
    )
    parser.add_argument(
        "--chunk-numbering",
        choices=["per-chrom", "global"],
        default="per-chrom",
        help="Whether {chunk} resets on each chromosome or increments globally. Default: per-chrom.",
    )
    parser.add_argument(
        "--format",
        choices=[
            "tsv4",
            "tsv5",
            "regions",
            "flare",
            "flare-args",
            "flare-direct",
            "flare-direct-args",
        ],
        default="tsv4",
        help=(
            "Output format. tsv4: chrom/start/end/out_prefix. tsv5 adds region. "
            "regions outputs chr:start-end only. flare adds split and converter columns. "
            "flare-args outputs only the 5 pre-split flare_subset_to_tractor_hybrid arguments. "
            "flare-direct and flare-direct-args target the converter's optional region argument."
        ),
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help="Print a header row.",
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.format in {"flare", "flare-args", "flare-direct", "flare-direct-args"}:
        if args.n_ancestries is None:
            die("--n-ancestries is required with flare output formats")
        if args.mac_threshold is None:
            die("--mac-threshold is required with flare output formats")
        args.n_ancestries = str(parse_int(args.n_ancestries, "--n-ancestries"))
        args.mac_threshold = str(parse_int(args.mac_threshold, "--mac-threshold"))

    requested = parse_chroms(args.chroms)
    if args.build:
        lengths = built_in_lengths(args.build, args.chrom_style)
        if requested is not None:
            lengths = select_chroms(lengths, requested, "auto")
    elif args.vcf:
        lengths = select_chroms(read_vcf_lengths(args.vcf), requested, args.chrom_style)
    else:
        lengths = select_chroms(read_chrom_sizes(args.chrom_sizes), requested, args.chrom_style)

    emit_manifest(lengths, args)


if __name__ == "__main__":
    main()
