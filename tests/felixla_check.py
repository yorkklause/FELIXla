"""Checks the suites run against FELIXla output, written independently of it.

These used to be done by asking FELIXla to check itself -- `--compare-vcfs` for
round trips and `--query` for dosages. Reading the packed files here instead
means a bug in the reader cannot hide a matching bug in the writer.
"""

import gzip
import pathlib
import struct


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_vcf_records(path):
    """(chrom, pos, id, ref, alt, [sample fields]) for each data line."""
    records = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            records.append((fields[0], int(fields[1]), fields[2], fields[3],
                            fields[4], fields[9:]))
    return records


def split_multiallelic(records):
    """One record per ALT, with each genotype recoded against that ALT alone.

    This is what a packed prefix holds, so comparing an original VCF with a
    round-tripped one means splitting the original the same way first.
    """
    out = []
    for chrom, pos, variant_id, ref, alt, samples in records:
        alts = alt.split(",")
        for index, one_alt in enumerate(alts, start=1):
            recoded = []
            for field in samples:
                call = field.split(":", 1)[0]
                separator = "|" if "|" in call else "/"
                alleles = call.split(separator)
                recoded.append(separator.join(
                    "." if a == "." else ("1" if a == str(index) else "0")
                    for a in alleles))
            out.append((chrom, pos, variant_id, ref, one_alt, recoded))
    return out


def compare_vcfs(lhs_path, rhs_path, max_diffs=20, split_lhs=False):
    """Differences between two VCFs, as human-readable strings."""
    lhs = read_vcf_records(lhs_path)
    if split_lhs:
        lhs = split_multiallelic(lhs)
    rhs = read_vcf_records(rhs_path)
    # The round trip keeps only GT, and rebuilds the ID from the split alleles.
    rhs = [(c, p, i, r, a, [f.split(":", 1)[0] for f in s])
           for (c, p, i, r, a, s) in rhs]
    lhs = [(c, p, i, r, a, [f.split(":", 1)[0] for f in s])
           for (c, p, i, r, a, s) in lhs]
    lhs = [(c, p, r, a, s) for (c, p, i, r, a, s) in lhs]
    rhs = [(c, p, r, a, s) for (c, p, i, r, a, s) in rhs]
    diffs = []
    if len(lhs) != len(rhs):
        diffs.append("record count %d vs %d" % (len(lhs), len(rhs)))
    for i, (a, b) in enumerate(zip(lhs, rhs)):
        if a == b:
            continue
        for label, x, y in zip(("CHROM", "POS", "REF", "ALT", "samples"), a, b):
            if x != y:
                diffs.append("record %d %s: %r vs %r" % (i, label, x, y))
                break
        if len(diffs) >= max_diffs:
            break
    return diffs


def _read_exact(handle, size, what):
    data = handle.read(size)
    if len(data) != size:
        raise AssertionError("short read of %s" % what)
    return data


def _read_string(handle):
    (length,) = struct.unpack("<I", _read_exact(handle, 4, "string length"))
    return _read_exact(handle, length, "string").decode()


def read_ancestry_blocks(prefix, n_ancestries, n_words):
    """(chrom, start, end, [mask words per ancestry]) for each block."""
    marks = pathlib.Path(str(prefix) + ".ancblock.mks").read_bytes()
    payload = pathlib.Path(str(prefix) + ".ancblock.bin").read_bytes()
    assert marks[:8] == b"TRANMKS1", marks[:8]

    blocks = []
    offset = 8
    while offset < len(marks):
        offset += 4                                   # block id
        (length,) = struct.unpack_from("<I", marks, offset)
        offset += 4
        chrom = marks[offset:offset + length].decode()
        offset += length
        start, end, payload_offset = struct.unpack_from("<qqQ", marks, offset)
        offset += 24

        masks = []
        for k in range(n_ancestries):
            base = payload_offset + k * n_words * 8
            masks.append(list(struct.unpack_from("<%dQ" % n_words, payload, base)))
        blocks.append((chrom, start, end, masks))
    return blocks


def read_meta(prefix):
    meta = {}
    with pathlib.Path(str(prefix) + ".meta").open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2:
                meta[parts[0]] = parts[1]
    return meta


def _read_variant_markers(prefix, kind):
    """One tuple per packed variant, in the marker file's own order.

    Field order follows write_common_mks_record / write_rare_mks_record:
      u64 local index, u32 global index, chrom, i64 pos, id, ref, alt,
      u32 alt index, [u32 block id (common only)], u64 payload offset,
      [u32 carrier count (rare only)], u32 mac
    """
    data = pathlib.Path(str(prefix) + ".%s.variant.mks" % kind).read_bytes()
    magic = b"TRCMMKS1" if kind == "common" else b"TRRAMKS1"
    assert data[:8] == magic, (kind, data[:8])

    out = []
    offset = 8

    def take_string(off):
        (length,) = struct.unpack_from("<I", data, off)
        off += 4
        return data[off:off + length].decode(), off + length

    while offset < len(data):
        offset += 8                                   # local index
        (global_index,) = struct.unpack_from("<I", data, offset)
        offset += 4
        chrom, offset = take_string(offset)
        (pos,) = struct.unpack_from("<q", data, offset)
        offset += 8
        variant_id, offset = take_string(offset)
        ref, offset = take_string(offset)
        alt, offset = take_string(offset)
        offset += 4                                   # alt index
        if kind == "common":
            offset += 4                               # block id
        (payload_offset,) = struct.unpack_from("<Q", data, offset)
        offset += 8
        n_carriers = None
        if kind == "rare":
            (n_carriers,) = struct.unpack_from("<I", data, offset)
            offset += 4
        offset += 4                                   # mac
        out.append({
            "global_index": global_index, "chrom": chrom, "pos": pos,
            "id": variant_id, "ref": ref, "alt": alt,
            "payload_offset": payload_offset, "n_carriers": n_carriers,
        })
    return out


def read_dosages(prefix):
    """Per-variant ancestry-specific dosages, computed from the packed files.

    Returns a dict keyed by (chrom, pos, ref, alt) holding, per sample index,
    (DSALL, [DS1..DSk]) -- the numbers FELIXassoc reads out of a prefix.
    """
    meta = read_meta(prefix)
    n_samples = int(meta["n_samples"])
    n_words = int(meta["n_words"])
    n_ancestries = int(meta["n_ancestries"])
    blocks = read_ancestry_blocks(prefix, n_ancestries, n_words)

    def ancestry_of(chrom, pos, hap):
        for block_chrom, start, end, masks in blocks:
            if block_chrom != chrom or not (start <= pos <= end):
                continue
            for k in range(n_ancestries):
                if masks[k][hap // 64] >> (hap % 64) & 1:
                    return k
            return None
        return None

    dosages = {}

    common_bits = pathlib.Path(str(prefix) + ".common.geno.bin").read_bytes()
    for marker in _read_variant_markers(prefix, "common"):
        words = struct.unpack_from("<%dQ" % n_words, common_bits,
                                   marker["payload_offset"])
        per_sample = []
        for i in range(n_samples):
            total = 0
            by_ancestry = [0] * n_ancestries
            for hap in (2 * i, 2 * i + 1):
                if not (words[hap // 64] >> (hap % 64) & 1):
                    continue
                total += 1
                k = ancestry_of(marker["chrom"], marker["pos"], hap)
                if k is not None:
                    by_ancestry[k] += 1
            per_sample.append((total, by_ancestry))
        dosages[(marker["chrom"], marker["pos"], marker["ref"],
                 marker["alt"])] = per_sample

    carriers = pathlib.Path(str(prefix) + ".rare.carrier.bin").read_bytes()
    for marker in _read_variant_markers(prefix, "rare"):
        per_sample = [(0, [0] * n_ancestries) for _ in range(n_samples)]
        base = marker["payload_offset"]
        for c in range(marker["n_carriers"]):
            _, anc_hap = struct.unpack_from("<II", carriers, base + c * 8)
            ancestry = (anc_hap >> 27) & 31
            hap = anc_hap & 0x07FFFFFF
            total, by_ancestry = per_sample[hap // 2]
            by_ancestry[ancestry] += 1
            per_sample[hap // 2] = (total + 1, by_ancestry)
        dosages[(marker["chrom"], marker["pos"], marker["ref"],
                 marker["alt"])] = per_sample

    return dosages
