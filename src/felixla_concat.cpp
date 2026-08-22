// felixla_concat.cpp
// designed by Kai, implemented by codex
//
// Concatenate complete FELIXla prefixes while validating and remapping every
// index, payload offset, variant ordinal, and ancestry block identifier.

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <unistd.h>

namespace {

constexpr uint32_t kPackedHapMask = 0x07ffffffu;
constexpr size_t kIoBufferBytes = 1u << 20;

constexpr std::array<uint8_t, 8> kCommonMksMagic = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
constexpr std::array<uint8_t, 8> kRareMksMagic = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
constexpr std::array<uint8_t, 8> kAncMksMagic = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
constexpr std::array<uint8_t, 8> kCommonIdxMagic = {'T', 'R', 'C', 'M', 'I', 'D', 'X', '2'};
constexpr std::array<uint8_t, 8> kRareIdxMagic = {'T', 'R', 'R', 'A', 'I', 'D', 'X', '2'};
constexpr std::array<uint8_t, 8> kAncIdxMagic = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};

const std::vector<std::string>& prefix_suffixes() {
    static const std::vector<std::string> suffixes = {
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
        ".meta"
    };
    return suffixes;
}

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

std::string trim(const std::string& value) {
    size_t begin = 0;
    while (begin < value.size() &&
           (value[begin] == ' ' || value[begin] == '\t' || value[begin] == '\r')) {
        ++begin;
    }
    size_t end = value.size();
    while (end > begin &&
           (value[end - 1] == ' ' || value[end - 1] == '\t' || value[end - 1] == '\r')) {
        --end;
    }
    return value.substr(begin, end - begin);
}

uint64_t parse_u64(const std::string& value, const std::string& field, const std::string& path) {
    if (value.empty() || value[0] == '-') {
        fail("invalid " + field + " in " + path + ": " + value);
    }
    size_t used = 0;
    unsigned long long parsed = 0;
    try {
        parsed = std::stoull(value, &used, 10);
    } catch (const std::exception&) {
        fail("invalid " + field + " in " + path + ": " + value);
    }
    if (used != value.size()) fail("invalid " + field + " in " + path + ": " + value);
    return static_cast<uint64_t>(parsed);
}

uint64_t checked_multiply(uint64_t a, uint64_t b, const std::string& what) {
    if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a) {
        fail("size overflow for " + what);
    }
    return a * b;
}

size_t checked_size_t(uint64_t value, const std::string& what) {
    if (value > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
        fail(what + " is too large for this platform");
    }
    return static_cast<size_t>(value);
}

struct Region {
    std::string chr;
    int64_t start = 0;
    int64_t end = 0;
};

std::optional<Region> parse_region(const std::string& value, const std::string& path) {
    size_t colon = value.find(':');
    size_t dash = value.find('-', colon == std::string::npos ? 0 : colon + 1);
    if (colon == std::string::npos || dash == std::string::npos || colon == 0 || dash <= colon + 1) {
        fail("invalid region in " + path + ": " + value);
    }

    Region region;
    region.chr = value.substr(0, colon);
    uint64_t start = parse_u64(value.substr(colon + 1, dash - colon - 1), "region start", path);
    uint64_t end = parse_u64(value.substr(dash + 1), "region end", path);
    if (start == 0 || start > end || end > static_cast<uint64_t>(std::numeric_limits<int64_t>::max())) {
        fail("invalid region in " + path + ": " + value);
    }
    region.start = static_cast<int64_t>(start);
    region.end = static_cast<int64_t>(end);
    return region;
}

struct DeclaredCounts {
    uint64_t global = 0;
    uint64_t common = 0;
    uint64_t rare = 0;
    uint64_t blocks = 0;
};

struct Meta {
    uint64_t format_version = 0;
    uint64_t n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t n_ancestries = 0;
    uint64_t rare_threshold = 0;
    std::optional<Region> region;
    std::optional<DeclaredCounts> counts;
};

Meta read_meta(const std::string& path) {
    std::ifstream in(path);
    if (!in) fail("cannot open " + path);

    std::unordered_map<std::string, std::string> values;
    std::string line;
    size_t line_no = 0;
    std::optional<Region> effective_region;
    while (std::getline(in, line)) {
        ++line_no;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        size_t tab = line.find('\t');
        if (tab == std::string::npos || tab == 0 || tab + 1 >= line.size()) {
            fail(path + ":" + std::to_string(line_no) + " is not a key-tab-value metadata row");
        }
        std::string key = line.substr(0, tab);
        std::string value = line.substr(tab + 1);
        auto previous = values.find(key);
        if (previous != values.end() && previous->second != value &&
            (key == "format_version" || key == "n_samples" || key == "n_haps" ||
             key == "n_words" || key == "n_ancestries" || key == "rare_threshold" ||
             key == "global_variants" || key == "common_variants" ||
             key == "rare_variants" || key == "ancestry_blocks")) {
            fail("conflicting " + key + " rows in " + path);
        }
        values[key] = value;
        if (key == "selected_region" || key == "extracted_region") {
            effective_region = parse_region(value, path);
        }
    }
    if (!in.eof()) fail("failed reading " + path);

    const char* required[] = {
        "format_version", "n_samples", "n_haps", "n_words", "n_ancestries", "rare_threshold"
    };
    for (const char* key : required) {
        if (!values.count(key)) fail("missing " + std::string(key) + " in " + path);
    }

    Meta meta;
    meta.format_version = parse_u64(values["format_version"], "format_version", path);
    meta.n_samples = parse_u64(values["n_samples"], "n_samples", path);
    meta.n_haps = parse_u64(values["n_haps"], "n_haps", path);
    meta.n_words = parse_u64(values["n_words"], "n_words", path);
    meta.n_ancestries = parse_u64(values["n_ancestries"], "n_ancestries", path);
    meta.rare_threshold = parse_u64(values["rare_threshold"], "rare_threshold", path);
    meta.region = effective_region;

    if (meta.format_version != 1) fail("unsupported format_version in " + path);
    if (meta.n_samples == 0) fail("n_samples must be positive in " + path);
    if (meta.n_haps != checked_multiply(meta.n_samples, 2, "n_haps")) {
        fail("n_haps must equal 2 * n_samples in " + path);
    }
    if (meta.n_haps > static_cast<uint64_t>(kPackedHapMask) + 1) {
        fail("n_haps exceeds the 27-bit packed haplotype limit in " + path);
    }
    if (meta.n_words != (meta.n_haps + 63) / 64) {
        fail("n_words does not match n_haps in " + path);
    }
    if (meta.n_ancestries == 0 || meta.n_ancestries > 32) {
        fail("n_ancestries must be in [1, 32] in " + path);
    }
    if (meta.rare_threshold > std::numeric_limits<uint32_t>::max()) {
        fail("rare_threshold exceeds uint32_t in " + path);
    }

    const char* count_keys[] = {
        "global_variants", "common_variants", "rare_variants", "ancestry_blocks"
    };
    size_t present = 0;
    for (const char* key : count_keys) present += values.count(key) ? 1 : 0;
    if (present != 0 && present != 4) fail("incomplete record-count manifest in " + path);
    if (present == 4) {
        DeclaredCounts counts;
        counts.global = parse_u64(values["global_variants"], "global_variants", path);
        counts.common = parse_u64(values["common_variants"], "common_variants", path);
        counts.rare = parse_u64(values["rare_variants"], "rare_variants", path);
        counts.blocks = parse_u64(values["ancestry_blocks"], "ancestry_blocks", path);
        if (counts.common > std::numeric_limits<uint64_t>::max() - counts.rare ||
            counts.global != counts.common + counts.rare) {
            fail("global_variants does not equal common_variants + rare_variants in " + path);
        }
        meta.counts = counts;
    }
    return meta;
}

uint64_t file_size(const std::string& path) {
    std::error_code ec;
    uintmax_t size = std::filesystem::file_size(path, ec);
    if (ec) fail("cannot stat " + path + ": " + ec.message());
    if (size > std::numeric_limits<uint64_t>::max()) fail("file is too large: " + path);
    return static_cast<uint64_t>(size);
}

class BinaryInput {
public:
    explicit BinaryInput(std::string path) : path_(std::move(path)), size_(file_size(path_)) {
        fp_ = std::fopen(path_.c_str(), "rb");
        if (!fp_) fail("cannot open " + path_ + ": " + std::strerror(errno));
        if (setvbuf(fp_, nullptr, _IOFBF, kIoBufferBytes) != 0) {
            int saved = errno;
            std::fclose(fp_);
            fp_ = nullptr;
            fail("setvbuf failed for " + path_ + ": " + std::strerror(saved));
        }
    }

    BinaryInput(const BinaryInput&) = delete;
    BinaryInput& operator=(const BinaryInput&) = delete;

    ~BinaryInput() {
        if (fp_) std::fclose(fp_);
    }

    uint64_t tell() const { return offset_; }
    uint64_t size() const { return size_; }
    uint64_t remaining() const { return size_ - offset_; }
    bool at_end() const { return offset_ == size_; }

    void read(void* data, size_t bytes) {
        if (static_cast<uint64_t>(bytes) > remaining()) {
            fail("unexpected EOF in " + path_ + " at byte " + std::to_string(offset_));
        }
        if (bytes != 0 && std::fread(data, 1, bytes, fp_) != bytes) {
            fail("failed reading " + path_ + " at byte " + std::to_string(offset_));
        }
        offset_ += bytes;
    }

    uint32_t read_u32() {
        uint8_t b[4];
        read(b, sizeof(b));
        return static_cast<uint32_t>(b[0]) |
               (static_cast<uint32_t>(b[1]) << 8) |
               (static_cast<uint32_t>(b[2]) << 16) |
               (static_cast<uint32_t>(b[3]) << 24);
    }

    uint64_t read_u64() {
        uint8_t b[8];
        read(b, sizeof(b));
        uint64_t value = 0;
        for (unsigned i = 0; i < 8; ++i) value |= static_cast<uint64_t>(b[i]) << (8 * i);
        return value;
    }

    int64_t read_i64() { return static_cast<int64_t>(read_u64()); }

    std::string read_string() {
        uint32_t length = read_u32();
        if (length > remaining()) fail("invalid string length in " + path_);
        std::string value(length, '\0');
        if (length != 0) read(value.data(), length);
        return value;
    }

    void check_magic(const std::array<uint8_t, 8>& expected) {
        uint8_t got[8];
        read(got, sizeof(got));
        if (!std::equal(expected.begin(), expected.end(), got)) fail("bad magic in " + path_);
    }

    void require_end() const {
        if (!at_end()) {
            fail("unexpected trailing bytes in " + path_ + " at byte " + std::to_string(offset_));
        }
    }

private:
    std::string path_;
    FILE* fp_ = nullptr;
    uint64_t size_ = 0;
    uint64_t offset_ = 0;
};

class BinaryOutput {
public:
    explicit BinaryOutput(std::string path) : path_(std::move(path)) {
        fp_ = std::fopen(path_.c_str(), "wb");
        if (!fp_) fail("cannot open " + path_ + ": " + std::strerror(errno));
        if (setvbuf(fp_, nullptr, _IOFBF, kIoBufferBytes) != 0) {
            int saved = errno;
            std::fclose(fp_);
            fp_ = nullptr;
            fail("setvbuf failed for " + path_ + ": " + std::strerror(saved));
        }
    }

    BinaryOutput(const BinaryOutput&) = delete;
    BinaryOutput& operator=(const BinaryOutput&) = delete;

    ~BinaryOutput() {
        if (fp_) std::fclose(fp_);
    }

    uint64_t tell() const { return offset_; }

    void write(const void* data, size_t bytes) {
        if (bytes != 0 && std::fwrite(data, 1, bytes, fp_) != bytes) {
            fail("failed writing " + path_ + " at byte " + std::to_string(offset_));
        }
        offset_ += bytes;
    }

    void write_u32(uint32_t value) {
        uint8_t b[4] = {
            static_cast<uint8_t>(value),
            static_cast<uint8_t>(value >> 8),
            static_cast<uint8_t>(value >> 16),
            static_cast<uint8_t>(value >> 24)
        };
        write(b, sizeof(b));
    }

    void write_u64(uint64_t value) {
        uint8_t b[8];
        for (unsigned i = 0; i < 8; ++i) b[i] = static_cast<uint8_t>(value >> (8 * i));
        write(b, sizeof(b));
    }

    void write_i64(int64_t value) { write_u64(static_cast<uint64_t>(value)); }

    void write_string(const std::string& value) {
        if (value.size() > std::numeric_limits<uint32_t>::max()) {
            fail("string is too long while writing " + path_);
        }
        write_u32(static_cast<uint32_t>(value.size()));
        write(value.data(), value.size());
    }

    void write_magic(const std::array<uint8_t, 8>& magic) { write(magic.data(), magic.size()); }

    void close() {
        if (!fp_) return;
        FILE* closing = fp_;
        fp_ = nullptr;
        if (std::fclose(closing) != 0) fail("failed closing " + path_ + ": " + std::strerror(errno));
        uint64_t observed = file_size(path_);
        if (observed != offset_) {
            fail("output size mismatch for " + path_ + ": expected " + std::to_string(offset_) +
                 ", observed " + std::to_string(observed));
        }
    }

private:
    std::string path_;
    FILE* fp_ = nullptr;
    uint64_t offset_ = 0;
};

uint32_t load_u32_le(const uint8_t* b) {
    return static_cast<uint32_t>(b[0]) |
           (static_cast<uint32_t>(b[1]) << 8) |
           (static_cast<uint32_t>(b[2]) << 16) |
           (static_cast<uint32_t>(b[3]) << 24);
}

uint64_t load_u64_le(const uint8_t* b) {
    uint64_t value = 0;
    for (unsigned i = 0; i < 8; ++i) value |= static_cast<uint64_t>(b[i]) << (8 * i);
    return value;
}

void store_u32_le(uint8_t* b, uint32_t value) {
    b[0] = static_cast<uint8_t>(value);
    b[1] = static_cast<uint8_t>(value >> 8);
    b[2] = static_cast<uint8_t>(value >> 16);
    b[3] = static_cast<uint8_t>(value >> 24);
}

struct InputInfo {
    std::string prefix;
    Meta meta;
};

std::vector<std::string> read_prefix_list(const std::string& path) {
    std::ifstream in(path);
    if (!in) fail("cannot open prefix list " + path);
    std::vector<std::string> prefixes;
    std::string line;
    size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        std::string prefix = trim(line);
        if (prefix.empty() || prefix[0] == '#') continue;
        if (prefix.find('\t') != std::string::npos || prefix.find('\n') != std::string::npos) {
            fail(path + ":" + std::to_string(line_no) + " must contain one prefix per line");
        }
        if (prefix.size() > 5 && prefix.compare(prefix.size() - 5, 5, ".meta") == 0) {
            prefix.resize(prefix.size() - 5);
        }
        prefixes.push_back(prefix);
    }
    if (!in.eof()) fail("failed reading prefix list " + path);
    if (prefixes.empty()) fail("prefix list is empty: " + path);
    return prefixes;
}

void preflight_prefix_files(const std::string& prefix) {
    for (const std::string& suffix : prefix_suffixes()) {
        std::string path = prefix + suffix;
        uint64_t size = file_size(path);
        if ((suffix.find(".mks") != std::string::npos ||
             suffix.find(".idx") != std::string::npos) && size < 8) {
            fail("truncated magic header in " + path);
        }
        if ((suffix == ".meta" || suffix == ".samples") && size == 0) {
            fail("empty required sidecar: " + path);
        }
    }
}

std::vector<std::string> read_samples(
    const std::string& path,
    uint64_t expected_count,
    const std::vector<std::string>* expected
) {
    std::ifstream in(path);
    if (!in) fail("cannot open " + path);
    std::vector<std::string> samples;
    if (!expected) samples.reserve(checked_size_t(expected_count, "sample count"));
    std::unordered_set<std::string> unique;
    if (!expected) unique.reserve(checked_size_t(expected_count, "sample count"));

    std::string line;
    uint64_t index = 0;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) fail("empty sample ID in " + path + " at line " + std::to_string(index + 1));
        if (index >= expected_count) fail("too many sample IDs in " + path);
        if (expected) {
            if (line != (*expected)[checked_size_t(index, "sample index")]) {
                fail("sample order mismatch in " + path + " at line " + std::to_string(index + 1));
            }
        } else {
            if (!unique.insert(line).second) fail("duplicate sample ID in " + path + ": " + line);
            samples.push_back(line);
        }
        ++index;
    }
    if (!in.eof()) fail("failed reading " + path);
    if (index != expected_count) {
        fail("sample count mismatch in " + path + ": meta declares " +
             std::to_string(expected_count) + ", file contains " + std::to_string(index));
    }
    return samples;
}

void require_same_meta(const InputInfo& first, const InputInfo& current) {
    const Meta& a = first.meta;
    const Meta& b = current.meta;
    if (a.format_version != b.format_version || a.n_samples != b.n_samples ||
        a.n_haps != b.n_haps || a.n_words != b.n_words ||
        a.n_ancestries != b.n_ancestries || a.rare_threshold != b.rare_threshold) {
        fail("incompatible FELIXla metadata between " + first.prefix + " and " + current.prefix +
             " (format/sample/haplotype/word/ancestry/MAC-threshold fields must match)");
    }
}

struct VariantRecord {
    bool common = false;
    uint64_t local_index = 0;
    uint32_t global_index = 0;
    std::string chr;
    int64_t pos = 0;
    std::string id;
    std::string ref;
    std::string alt;
    uint32_t alt_index = 0;
    uint32_t block_id = 0;
    uint64_t payload_offset = 0;
    uint32_t n_carriers = 0;
    uint32_t mac = 0;
};

struct AncRecord {
    uint32_t block_id = 0;
    std::string chr;
    int64_t start = 0;
    int64_t end = 0;
    uint64_t payload_offset = 0;
};

class OutputFiles {
public:
    explicit OutputFiles(const std::string& prefix)
        : common_bin(prefix + ".common.geno.bin"),
          common_mks(prefix + ".common.variant.mks"),
          common_idx(prefix + ".common.variant.idx"),
          rare_bin(prefix + ".rare.carrier.bin"),
          rare_mks(prefix + ".rare.variant.mks"),
          rare_idx(prefix + ".rare.variant.idx"),
          anc_bin(prefix + ".ancblock.bin"),
          anc_mks(prefix + ".ancblock.mks"),
          anc_idx(prefix + ".ancblock.idx") {
        common_mks.write_magic(kCommonMksMagic);
        common_idx.write_magic(kCommonIdxMagic);
        rare_mks.write_magic(kRareMksMagic);
        rare_idx.write_magic(kRareIdxMagic);
        anc_mks.write_magic(kAncMksMagic);
        anc_idx.write_magic(kAncIdxMagic);
    }

    void close() {
        common_bin.close();
        common_mks.close();
        common_idx.close();
        rare_bin.close();
        rare_mks.close();
        rare_idx.close();
        anc_bin.close();
        anc_mks.close();
        anc_idx.close();
    }

    BinaryOutput common_bin;
    BinaryOutput common_mks;
    BinaryOutput common_idx;
    BinaryOutput rare_bin;
    BinaryOutput rare_mks;
    BinaryOutput rare_idx;
    BinaryOutput anc_bin;
    BinaryOutput anc_mks;
    BinaryOutput anc_idx;
};

struct MergeState {
    uint64_t global_variants = 0;
    uint64_t common_variants = 0;
    uint64_t rare_variants = 0;
    uint64_t ancestry_blocks = 0;
    std::unordered_map<std::string, int64_t> last_variant_pos;
    std::unordered_map<std::string, int64_t> last_block_end;
};

struct CoordinateSpan {
    int64_t start = 0;
    int64_t end = 0;
};

struct PrefixSummary {
    DeclaredCounts counts;
    std::unordered_map<std::string, CoordinateSpan> variant_spans;
    std::unordered_map<std::string, CoordinateSpan> block_spans;
};

void extend_span(
    std::unordered_map<std::string, CoordinateSpan>& spans,
    const std::string& chr,
    int64_t start,
    int64_t end
) {
    auto inserted = spans.emplace(chr, CoordinateSpan{start, end});
    if (!inserted.second) {
        inserted.first->second.start = std::min(inserted.first->second.start, start);
        inserted.first->second.end = std::max(inserted.first->second.end, end);
    }
}

void validate_variant_fields(const VariantRecord& record, const Meta& meta, const std::string& path) {
    if (record.chr.empty() || record.pos <= 0 || record.ref.empty() || record.alt.empty()) {
        fail("invalid variant coordinate or allele in " + path);
    }
    if (record.ref == record.alt) fail("REF equals ALT in " + path);
    if (record.alt_index == 0) fail("ALT allele index is zero in " + path);
    if (record.mac > meta.n_haps) fail("variant MAC exceeds n_haps in " + path);
    if (record.common && record.mac <= meta.rare_threshold) {
        fail("common variant MAC is not above rare_threshold in " + path);
    }
    if (!record.common && record.mac > meta.rare_threshold) {
        fail("rare variant MAC exceeds rare_threshold in " + path);
    }
    if (!record.common && record.n_carriers != record.mac) {
        fail("rare marker n_carriers does not equal MAC in " + path);
    }
}

void write_common_marker(
    BinaryOutput& out,
    const VariantRecord& record,
    uint64_t local_index,
    uint32_t global_index,
    uint32_t block_id,
    uint64_t payload_offset
) {
    out.write_u64(local_index);
    out.write_u32(global_index);
    out.write_string(record.chr);
    out.write_i64(record.pos);
    out.write_string(record.id);
    out.write_string(record.ref);
    out.write_string(record.alt);
    out.write_u32(record.alt_index);
    out.write_u32(block_id);
    out.write_u64(payload_offset);
    out.write_u32(record.mac);
}

void write_rare_marker(
    BinaryOutput& out,
    const VariantRecord& record,
    uint64_t local_index,
    uint32_t global_index,
    uint64_t payload_offset
) {
    out.write_u64(local_index);
    out.write_u32(global_index);
    out.write_string(record.chr);
    out.write_i64(record.pos);
    out.write_string(record.id);
    out.write_string(record.ref);
    out.write_string(record.alt);
    out.write_u32(record.alt_index);
    out.write_u64(payload_offset);
    out.write_u32(record.n_carriers);
    out.write_u32(record.mac);
}

void write_anc_marker(
    BinaryOutput& out,
    const AncRecord& record,
    uint32_t block_id,
    uint64_t payload_offset
) {
    out.write_u32(block_id);
    out.write_string(record.chr);
    out.write_i64(record.start);
    out.write_i64(record.end);
    out.write_u64(payload_offset);
}

std::optional<VariantRecord> read_common_record(
    BinaryInput& mks,
    BinaryInput& idx,
    BinaryInput& payload,
    uint64_t expected_local,
    const std::string& prefix
) {
    if (mks.at_end()) return std::nullopt;
    uint64_t marker_offset = mks.tell();
    VariantRecord record;
    record.common = true;
    record.local_index = mks.read_u64();
    record.global_index = mks.read_u32();
    record.chr = mks.read_string();
    record.pos = mks.read_i64();
    record.id = mks.read_string();
    record.ref = mks.read_string();
    record.alt = mks.read_string();
    record.alt_index = mks.read_u32();
    record.block_id = mks.read_u32();
    record.payload_offset = mks.read_u64();
    record.mac = mks.read_u32();

    if (idx.remaining() < 28) fail("missing common index record for " + prefix);
    uint64_t idx_local = idx.read_u64();
    uint32_t idx_global = idx.read_u32();
    uint64_t idx_marker = idx.read_u64();
    uint64_t idx_payload = idx.read_u64();
    if (record.local_index != expected_local || idx_local != record.local_index ||
        idx_global != record.global_index || idx_marker != marker_offset ||
        idx_payload != record.payload_offset || record.payload_offset != payload.tell()) {
        fail("common marker/index/offset mismatch in " + prefix);
    }
    return record;
}

std::optional<VariantRecord> read_rare_record(
    BinaryInput& mks,
    BinaryInput& idx,
    BinaryInput& payload,
    uint64_t expected_local,
    const std::string& prefix
) {
    if (mks.at_end()) return std::nullopt;
    uint64_t marker_offset = mks.tell();
    VariantRecord record;
    record.common = false;
    record.local_index = mks.read_u64();
    record.global_index = mks.read_u32();
    record.chr = mks.read_string();
    record.pos = mks.read_i64();
    record.id = mks.read_string();
    record.ref = mks.read_string();
    record.alt = mks.read_string();
    record.alt_index = mks.read_u32();
    record.payload_offset = mks.read_u64();
    record.n_carriers = mks.read_u32();
    record.mac = mks.read_u32();

    if (idx.remaining() < 32) fail("missing rare index record for " + prefix);
    uint64_t idx_local = idx.read_u64();
    uint32_t idx_global = idx.read_u32();
    uint64_t idx_marker = idx.read_u64();
    uint64_t idx_payload = idx.read_u64();
    uint32_t idx_carriers = idx.read_u32();
    if (record.local_index != expected_local || idx_local != record.local_index ||
        idx_global != record.global_index || idx_marker != marker_offset ||
        idx_payload != record.payload_offset || idx_carriers != record.n_carriers ||
        record.payload_offset != payload.tell()) {
        fail("rare marker/index/offset mismatch in " + prefix);
    }
    return record;
}

bool position_is_covered(
    const std::unordered_map<std::string, std::vector<std::pair<int64_t, int64_t>>>& intervals,
    const std::string& chr,
    int64_t pos
) {
    auto found = intervals.find(chr);
    if (found == intervals.end()) return false;
    const auto& ranges = found->second;
    auto it = std::upper_bound(
        ranges.begin(), ranges.end(), pos,
        [](int64_t value, const std::pair<int64_t, int64_t>& range) { return value < range.first; }
    );
    if (it == ranges.begin()) return false;
    --it;
    return pos <= it->second;
}

uint64_t process_ancestry_blocks(
    const InputInfo& input,
    OutputFiles* output,
    MergeState& state,
    PrefixSummary& summary,
    std::vector<AncRecord>& blocks,
    std::unordered_map<std::string, std::vector<std::pair<int64_t, int64_t>>>& intervals
) {
    BinaryInput mks(input.prefix + ".ancblock.mks");
    BinaryInput idx(input.prefix + ".ancblock.idx");
    BinaryInput payload(input.prefix + ".ancblock.bin");
    mks.check_magic(kAncMksMagic);
    idx.check_magic(kAncIdxMagic);
    if (idx.remaining() % 20 != 0) fail("truncated ancestry index in " + input.prefix);

    uint64_t block_words = checked_multiply(input.meta.n_ancestries, input.meta.n_words, "ancestry block words");
    uint64_t block_bytes_u64 = checked_multiply(block_words, 8, "ancestry block bytes");
    size_t block_bytes = checked_size_t(block_bytes_u64, "ancestry block");
    std::vector<uint8_t> buffer(block_bytes);
    std::unordered_map<std::string, int64_t> local_last_end;
    std::unordered_set<std::string> seen_in_input;
    uint64_t input_blocks = 0;
    uint64_t output_block_base = state.ancestry_blocks;

    while (!mks.at_end()) {
        uint64_t marker_offset = mks.tell();
        AncRecord record;
        record.block_id = mks.read_u32();
        record.chr = mks.read_string();
        record.start = mks.read_i64();
        record.end = mks.read_i64();
        record.payload_offset = mks.read_u64();
        if (record.block_id != input_blocks || record.chr.empty() || record.start <= 0 ||
            record.start > record.end || record.payload_offset != payload.tell()) {
            fail("invalid ancestry marker sequence in " + input.prefix);
        }

        if (idx.remaining() < 20) fail("missing ancestry index record in " + input.prefix);
        uint32_t idx_block = idx.read_u32();
        uint64_t idx_marker = idx.read_u64();
        uint64_t idx_payload = idx.read_u64();
        if (idx_block != record.block_id || idx_marker != marker_offset ||
            idx_payload != record.payload_offset) {
            fail("ancestry marker/index mismatch in " + input.prefix);
        }

        auto local = local_last_end.find(record.chr);
        if (local != local_last_end.end() && record.start <= local->second) {
            fail("overlapping or unsorted ancestry blocks in " + input.prefix + " on " + record.chr);
        }
        if (seen_in_input.insert(record.chr).second) {
            auto previous = state.last_block_end.find(record.chr);
            if (previous != state.last_block_end.end() && record.start <= previous->second) {
                fail("ancestry block overlap between concatenated prefixes on " + record.chr +
                     " at " + std::to_string(record.start));
            }
        }
        local_last_end[record.chr] = record.end;
        state.last_block_end[record.chr] = record.end;

        payload.read(buffer.data(), buffer.size());
        uint64_t tail_bits = input.meta.n_haps % 64;
        uint64_t tail_mask = tail_bits == 0 ? std::numeric_limits<uint64_t>::max() :
                             ((uint64_t{1} << tail_bits) - 1);
        for (uint64_t word = 0; word < input.meta.n_words; ++word) {
            uint64_t combined = 0;
            for (uint64_t ancestry = 0; ancestry < input.meta.n_ancestries; ++ancestry) {
                uint64_t word_index = ancestry * input.meta.n_words + word;
                uint64_t bits = load_u64_le(buffer.data() + checked_size_t(word_index * 8, "ancestry word"));
                if ((combined & bits) != 0) {
                    fail("overlapping ancestry masks in block " + std::to_string(record.block_id) +
                         " of " + input.prefix);
                }
                combined |= bits;
            }
            uint64_t expected = word + 1 == input.meta.n_words ? tail_mask :
                                std::numeric_limits<uint64_t>::max();
            if (combined != expected) {
                fail("ancestry masks do not partition all haplotypes in block " +
                     std::to_string(record.block_id) + " of " + input.prefix);
            }
        }

        if (state.ancestry_blocks >= std::numeric_limits<uint32_t>::max()) {
            fail("concatenated ancestry block count exceeds uint32_t limit");
        }
        if (output) {
            uint32_t new_block = static_cast<uint32_t>(state.ancestry_blocks);
            uint64_t new_payload = output->anc_bin.tell();
            output->anc_bin.write(buffer.data(), buffer.size());
            uint64_t new_marker = output->anc_mks.tell();
            write_anc_marker(output->anc_mks, record, new_block, new_payload);
            output->anc_idx.write_u32(new_block);
            output->anc_idx.write_u64(new_marker);
            output->anc_idx.write_u64(new_payload);
        }

        blocks.push_back(record);
        intervals[record.chr].push_back({record.start, record.end});
        extend_span(summary.block_spans, record.chr, record.start, record.end);
        ++input_blocks;
        ++state.ancestry_blocks;
    }

    mks.require_end();
    idx.require_end();
    payload.require_end();
    if (state.ancestry_blocks != output_block_base + input_blocks) {
        fail("internal ancestry block count mismatch");
    }
    return input_blocks;
}

void validate_variant_order(
    const VariantRecord& record,
    const std::string& prefix,
    std::unordered_map<std::string, int64_t>& local_last,
    std::unordered_set<std::string>& seen_in_input,
    MergeState& state
) {
    auto local = local_last.find(record.chr);
    if (local != local_last.end() && record.pos < local->second) {
        fail("unsorted variants in " + prefix + " on " + record.chr);
    }
    if (seen_in_input.insert(record.chr).second) {
        auto previous = state.last_variant_pos.find(record.chr);
        if (previous != state.last_variant_pos.end() && record.pos <= previous->second) {
            fail("variant overlap or unsorted prefix list on " + record.chr +
                 " at " + std::to_string(record.pos));
        }
    }
    local_last[record.chr] = record.pos;
    state.last_variant_pos[record.chr] = record.pos;
}

std::pair<uint64_t, uint64_t> process_variants(
    const InputInfo& input,
    uint64_t output_block_base,
    const std::vector<AncRecord>& blocks,
    const std::unordered_map<std::string, std::vector<std::pair<int64_t, int64_t>>>& intervals,
    OutputFiles* output,
    MergeState& state,
    PrefixSummary& summary
) {
    BinaryInput common_mks(input.prefix + ".common.variant.mks");
    BinaryInput common_idx(input.prefix + ".common.variant.idx");
    BinaryInput common_bin(input.prefix + ".common.geno.bin");
    BinaryInput rare_mks(input.prefix + ".rare.variant.mks");
    BinaryInput rare_idx(input.prefix + ".rare.variant.idx");
    BinaryInput rare_bin(input.prefix + ".rare.carrier.bin");
    common_mks.check_magic(kCommonMksMagic);
    common_idx.check_magic(kCommonIdxMagic);
    rare_mks.check_magic(kRareMksMagic);
    rare_idx.check_magic(kRareIdxMagic);
    if (common_idx.remaining() % 28 != 0) fail("truncated common index in " + input.prefix);
    if (rare_idx.remaining() % 32 != 0) fail("truncated rare index in " + input.prefix);

    uint64_t common_bytes_u64 = checked_multiply(input.meta.n_words, 8, "common genotype bytes");
    size_t common_bytes = checked_size_t(common_bytes_u64, "common genotype");
    std::vector<uint8_t> common_buffer(common_bytes);
    std::vector<uint8_t> rare_buffer;
    uint64_t input_common = 0;
    uint64_t input_rare = 0;
    uint64_t input_global = 0;
    std::unordered_map<std::string, int64_t> local_last;
    std::unordered_set<std::string> seen_in_input;

    std::optional<VariantRecord> common = read_common_record(
        common_mks, common_idx, common_bin, input_common, input.prefix);
    std::optional<VariantRecord> rare = read_rare_record(
        rare_mks, rare_idx, rare_bin, input_rare, input.prefix);

    while (common || rare) {
        bool take_common = common && !rare;
        if (common && rare) {
            if (common->global_index == rare->global_index) {
                fail("global variant exists in both common and rare streams in " + input.prefix);
            }
            take_common = common->global_index < rare->global_index;
        }
        VariantRecord record = take_common ? *common : *rare;
        if (record.global_index != input_global) {
            fail("non-contiguous global variant index in " + input.prefix + ": expected " +
                 std::to_string(input_global) + ", observed " + std::to_string(record.global_index));
        }
        validate_variant_fields(record, input.meta, input.prefix);
        validate_variant_order(record, input.prefix, local_last, seen_in_input, state);
        extend_span(summary.variant_spans, record.chr, record.pos, record.pos);
        if (state.global_variants >= std::numeric_limits<uint32_t>::max()) {
            fail("concatenated global variant count exceeds uint32_t limit");
        }
        uint32_t new_global = static_cast<uint32_t>(state.global_variants);

        if (take_common) {
            if (record.block_id >= blocks.size()) {
                fail("common variant references missing ancestry block in " + input.prefix);
            }
            const AncRecord& block = blocks[record.block_id];
            if (record.chr != block.chr || record.pos < block.start || record.pos > block.end) {
                fail("common variant is outside its ancestry block in " + input.prefix);
            }
            common_bin.read(common_buffer.data(), common_buffer.size());
            uint64_t popcount = 0;
            for (uint64_t word = 0; word < input.meta.n_words; ++word) {
                uint64_t bits = load_u64_le(
                    common_buffer.data() + checked_size_t(word * 8, "common genotype word"));
                if (word + 1 == input.meta.n_words && input.meta.n_haps % 64 != 0) {
                    uint64_t valid = (uint64_t{1} << (input.meta.n_haps % 64)) - 1;
                    if ((bits & ~valid) != 0) fail("nonzero common genotype tail bits in " + input.prefix);
                }
                popcount += static_cast<uint64_t>(__builtin_popcountll(bits));
            }
            if (popcount != record.mac) fail("common genotype popcount does not equal MAC in " + input.prefix);

            uint64_t remapped_block_u64 = output_block_base + record.block_id;
            if (remapped_block_u64 > std::numeric_limits<uint32_t>::max()) {
                fail("remapped ancestry block ID exceeds uint32_t");
            }
            if (output) {
                uint64_t new_payload = output->common_bin.tell();
                output->common_bin.write(common_buffer.data(), common_buffer.size());
                uint64_t new_marker = output->common_mks.tell();
                uint32_t remapped_block = static_cast<uint32_t>(remapped_block_u64);
                write_common_marker(
                    output->common_mks,
                    record,
                    state.common_variants,
                    new_global,
                    remapped_block,
                    new_payload
                );
                output->common_idx.write_u64(state.common_variants);
                output->common_idx.write_u32(new_global);
                output->common_idx.write_u64(new_marker);
                output->common_idx.write_u64(new_payload);
            }
            ++state.common_variants;
            ++input_common;
            common = read_common_record(
                common_mks, common_idx, common_bin, input_common, input.prefix);
        } else {
            if (!position_is_covered(intervals, record.chr, record.pos)) {
                fail("rare variant is not covered by an ancestry block in " + input.prefix);
            }
            uint64_t rare_bytes_u64 = checked_multiply(record.n_carriers, 8, "rare carrier bytes");
            size_t rare_bytes = checked_size_t(rare_bytes_u64, "rare carrier payload");
            rare_buffer.resize(rare_bytes);
            rare_bin.read(rare_buffer.data(), rare_buffer.size());
            uint32_t previous_hap = 0;
            bool have_previous = false;
            for (uint32_t i = 0; i < record.n_carriers; ++i) {
                uint8_t* carrier = rare_buffer.data() + static_cast<size_t>(i) * 8;
                uint32_t old_global = load_u32_le(carrier);
                uint32_t anc_hap = load_u32_le(carrier + 4);
                uint32_t hap = anc_hap & kPackedHapMask;
                uint32_t ancestry = anc_hap >> 27;
                if (old_global != record.global_index) {
                    fail("rare carrier global index mismatch in " + input.prefix);
                }
                if (hap >= input.meta.n_haps || ancestry >= input.meta.n_ancestries) {
                    fail("rare carrier ancestry or haplotype is out of range in " + input.prefix);
                }
                if (have_previous && hap <= previous_hap) {
                    fail("rare carrier haplotypes are duplicated or unsorted in " + input.prefix);
                }
                previous_hap = hap;
                have_previous = true;
                if (output) store_u32_le(carrier, new_global);
            }

            if (output) {
                uint64_t new_payload = output->rare_bin.tell();
                output->rare_bin.write(rare_buffer.data(), rare_buffer.size());
                uint64_t new_marker = output->rare_mks.tell();
                write_rare_marker(
                    output->rare_mks,
                    record,
                    state.rare_variants,
                    new_global,
                    new_payload
                );
                output->rare_idx.write_u64(state.rare_variants);
                output->rare_idx.write_u32(new_global);
                output->rare_idx.write_u64(new_marker);
                output->rare_idx.write_u64(new_payload);
                output->rare_idx.write_u32(record.n_carriers);
            }
            ++state.rare_variants;
            ++input_rare;
            rare = read_rare_record(rare_mks, rare_idx, rare_bin, input_rare, input.prefix);
        }

        ++input_global;
        ++state.global_variants;
    }

    common_mks.require_end();
    common_idx.require_end();
    common_bin.require_end();
    rare_mks.require_end();
    rare_idx.require_end();
    rare_bin.require_end();
    return {input_common, input_rare};
}

PrefixSummary process_prefix(const InputInfo& input, OutputFiles* output, MergeState& state) {
    PrefixSummary summary;
    uint64_t block_base = state.ancestry_blocks;
    std::vector<AncRecord> blocks;
    std::unordered_map<std::string, std::vector<std::pair<int64_t, int64_t>>> intervals;
    uint64_t input_blocks = process_ancestry_blocks(
        input, output, state, summary, blocks, intervals);
    auto variant_counts = process_variants(
        input, block_base, blocks, intervals, output, state, summary);

    DeclaredCounts& actual = summary.counts;
    actual.common = variant_counts.first;
    actual.rare = variant_counts.second;
    actual.global = actual.common + actual.rare;
    actual.blocks = input_blocks;
    if (input.meta.counts) {
        const DeclaredCounts& declared = *input.meta.counts;
        if (declared.global != actual.global || declared.common != actual.common ||
            declared.rare != actual.rare || declared.blocks != actual.blocks) {
            fail("record counts in " + input.prefix + ".meta do not match its binary files");
        }
    }
    return summary;
}

class TempPrefixGuard {
public:
    explicit TempPrefixGuard(std::string prefix) : prefix_(std::move(prefix)) {}
    ~TempPrefixGuard() {
        if (committed_) return;
        for (const std::string& suffix : prefix_suffixes()) {
            std::remove((prefix_ + suffix).c_str());
        }
    }
    void commit() { committed_ = true; }

private:
    std::string prefix_;
    bool committed_ = false;
};

std::string make_temp_prefix(const std::string& out_prefix) {
    std::string base = out_prefix + ".concat.tmp." + std::to_string(static_cast<long long>(getpid()));
    for (unsigned attempt = 0; attempt < 1000; ++attempt) {
        std::string candidate = attempt == 0 ? base : base + "." + std::to_string(attempt);
        bool exists = false;
        for (const std::string& suffix : prefix_suffixes()) {
            std::error_code ec;
            if (std::filesystem::exists(candidate + suffix, ec)) {
                exists = true;
                break;
            }
            if (ec) fail("cannot inspect temporary output path: " + ec.message());
        }
        if (!exists) return candidate;
    }
    fail("cannot allocate a temporary output prefix for " + out_prefix);
}

void write_samples_file(const std::string& path, const std::vector<std::string>& samples) {
    std::ofstream out(path, std::ios::binary);
    if (!out) fail("cannot open " + path);
    for (const std::string& sample : samples) out << sample << '\n';
    out.close();
    if (!out) fail("failed writing " + path);
}

void write_meta_file(
    const std::string& path,
    const Meta& meta,
    const std::string& list_path,
    const std::vector<InputInfo>& inputs,
    const MergeState& state
) {
    std::ofstream out(path);
    if (!out) fail("cannot open " + path);
    out << "format_version\t" << meta.format_version << '\n';
    out << "n_samples\t" << meta.n_samples << '\n';
    out << "n_haps\t" << meta.n_haps << '\n';
    out << "n_words\t" << meta.n_words << '\n';
    out << "n_ancestries\t" << meta.n_ancestries << '\n';
    out << "rare_threshold\t" << meta.rare_threshold << '\n';
    out << "global_variants\t" << state.global_variants << '\n';
    out << "common_variants\t" << state.common_variants << '\n';
    out << "rare_variants\t" << state.rare_variants << '\n';
    out << "ancestry_blocks\t" << state.ancestry_blocks << '\n';
    out << "concat_list\t" << list_path << '\n';
    out << "concat_inputs\t" << inputs.size() << '\n';
    out << "concat_validation\tfull-preflight-then-merge\n";
    out << "integrity_checks\tmeta,samples,magic,index,offset,eof,mac,masks,ordering\n";
    for (size_t i = 0; i < inputs.size(); ++i) {
        out << "concat_source_" << (i + 1) << '\t' << inputs[i].prefix << '\n';
    }
    out.close();
    if (!out) fail("failed writing " + path);
}

void publish_temp_prefix(const std::string& temp_prefix, const std::string& out_prefix) {
    for (const std::string& suffix : prefix_suffixes()) {
        if (suffix == ".meta") continue;
        if (std::rename((temp_prefix + suffix).c_str(), (out_prefix + suffix).c_str()) != 0) {
            fail("cannot publish " + out_prefix + suffix + ": " + std::strerror(errno));
        }
    }
    if (std::rename((temp_prefix + ".meta").c_str(), (out_prefix + ".meta").c_str()) != 0) {
        fail("cannot publish " + out_prefix + ".meta: " + std::strerror(errno));
    }
}

void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s prefix_list.txt out_prefix\n\n"
        "The list contains one complete FELIXla prefix per line, in output order.\n"
        "Blank lines and lines beginning with # are ignored. A trailing .meta is accepted.\n",
        prog
    );
}

struct Candidate {
    std::string prefix;
    std::optional<InputInfo> input;
    std::optional<PrefixSummary> summary;
    std::vector<std::string> errors;
};

void add_candidate_error(Candidate& candidate, const std::string& message) {
    if (std::find(candidate.errors.begin(), candidate.errors.end(), message) == candidate.errors.end()) {
        candidate.errors.push_back(message);
    }
}

void validate_declared_region_order(std::vector<Candidate>& candidates) {
    struct Previous {
        int64_t end = 0;
        std::string prefix;
    };
    std::unordered_map<std::string, Previous> previous_by_chr;
    for (Candidate& candidate : candidates) {
        if (!candidate.input || !candidate.input->meta.region) continue;
        const Region& region = *candidate.input->meta.region;
        auto found = previous_by_chr.find(region.chr);
        if (found != previous_by_chr.end() && region.start <= found->second.end) {
            add_candidate_error(
                candidate,
                "declared region overlaps or precedes " + found->second.prefix + " on " +
                    region.chr + " at " + std::to_string(region.start)
            );
        }
        if (found == previous_by_chr.end() || region.end > found->second.end) {
            previous_by_chr[region.chr] = Previous{region.end, candidate.prefix};
        }
    }
}

void validate_actual_span_order(
    std::vector<Candidate>& candidates,
    bool ancestry_blocks
) {
    struct Previous {
        int64_t end = 0;
        std::string prefix;
    };
    std::unordered_map<std::string, Previous> previous_by_chr;
    for (Candidate& candidate : candidates) {
        if (!candidate.summary) continue;
        const auto& spans = ancestry_blocks ? candidate.summary->block_spans :
                                              candidate.summary->variant_spans;
        for (const auto& item : spans) {
            const std::string& chr = item.first;
            const CoordinateSpan& span = item.second;
            auto found = previous_by_chr.find(chr);
            if (found != previous_by_chr.end() && span.start <= found->second.end) {
                add_candidate_error(
                    candidate,
                    std::string(ancestry_blocks ? "ancestry-block" : "variant") +
                        " span overlaps or precedes " + found->second.prefix + " on " + chr +
                        " at " + std::to_string(span.start)
                );
            }
            if (found == previous_by_chr.end() || span.end > found->second.end) {
                previous_by_chr[chr] = Previous{span.end, candidate.prefix};
            }
        }
    }
}

int run_concat(const std::string& list_path, const std::string& out_prefix) {
    std::vector<std::string> prefixes = read_prefix_list(list_path);
    std::string normalized_out = std::filesystem::absolute(out_prefix).lexically_normal().string();
    std::vector<Candidate> candidates;
    candidates.reserve(prefixes.size());
    for (const std::string& prefix : prefixes) {
        Candidate candidate;
        candidate.prefix = prefix;
        candidates.push_back(std::move(candidate));
    }

    std::unordered_map<std::string, size_t> first_by_normalized_prefix;
    for (size_t i = 0; i < candidates.size(); ++i) {
        Candidate& candidate = candidates[i];
        std::cerr << "FELIXla concat preflight: [" << (i + 1) << '/' << candidates.size()
                  << "] " << candidate.prefix << '\n';
        std::string normalized_input =
            std::filesystem::absolute(candidate.prefix).lexically_normal().string();
        if (normalized_input == normalized_out) {
            add_candidate_error(candidate, "output prefix is also an input prefix");
        }
        auto duplicate = first_by_normalized_prefix.emplace(normalized_input, i);
        if (!duplicate.second) {
            add_candidate_error(
                candidate,
                "duplicate list entry; first occurrence is " +
                    candidates[duplicate.first->second].prefix
            );
        }
        try {
            preflight_prefix_files(candidate.prefix);
        } catch (const std::exception& error) {
            add_candidate_error(candidate, error.what());
        }
        try {
            candidate.input = InputInfo{
                candidate.prefix,
                read_meta(candidate.prefix + ".meta")
            };
        } catch (const std::exception& error) {
            add_candidate_error(candidate, error.what());
        }
    }

    const InputInfo* metadata_baseline = nullptr;
    for (Candidate& candidate : candidates) {
        if (!candidate.input) continue;
        if (!metadata_baseline) {
            metadata_baseline = &*candidate.input;
            continue;
        }
        try {
            require_same_meta(*metadata_baseline, *candidate.input);
        } catch (const std::exception& error) {
            add_candidate_error(candidate, error.what());
        }
    }

    std::vector<std::string> samples;
    std::string sample_baseline_prefix;
    for (Candidate& candidate : candidates) {
        if (!candidate.input) continue;
        try {
            std::vector<std::string> current = read_samples(
                candidate.prefix + ".samples",
                candidate.input->meta.n_samples,
                nullptr
            );
            if (samples.empty()) {
                samples = std::move(current);
                sample_baseline_prefix = candidate.prefix;
            } else if (current != samples) {
                add_candidate_error(
                    candidate,
                    "sample IDs/order differ from " + sample_baseline_prefix
                );
            }
        } catch (const std::exception& error) {
            add_candidate_error(candidate, error.what());
        }
    }

    for (Candidate& candidate : candidates) {
        if (!candidate.input) continue;
        try {
            MergeState validation_state;
            candidate.summary = process_prefix(*candidate.input, nullptr, validation_state);
            const DeclaredCounts& counts = candidate.summary->counts;
            std::cerr << "  valid binary structure: variants=" << counts.global
                      << " common=" << counts.common
                      << " rare=" << counts.rare
                      << " ancestry_blocks=" << counts.blocks << '\n';
        } catch (const std::exception& error) {
            add_candidate_error(candidate, error.what());
        }
    }

    validate_declared_region_order(candidates);
    validate_actual_span_order(candidates, true);
    validate_actual_span_order(candidates, false);

    size_t bad_prefixes = 0;
    for (const Candidate& candidate : candidates) bad_prefixes += candidate.errors.empty() ? 0 : 1;
    if (bad_prefixes != 0) {
        std::cerr << "ERROR: FELIXla concat preflight found " << bad_prefixes
                  << " problematic prefix(es) out of " << candidates.size()
                  << "; merge was not started.\n";
        for (size_t i = 0; i < candidates.size(); ++i) {
            const Candidate& candidate = candidates[i];
            if (candidate.errors.empty()) continue;
            std::cerr << "[" << (i + 1) << "] " << candidate.prefix << '\n';
            for (const std::string& error : candidate.errors) {
                std::cerr << "  - " << error << '\n';
            }
        }
        return 1;
    }

    std::vector<InputInfo> inputs;
    inputs.reserve(candidates.size());
    for (Candidate& candidate : candidates) inputs.push_back(std::move(*candidate.input));
    std::cerr << "FELIXla concat preflight passed for all " << inputs.size()
              << " prefix(es); starting merge.\n";

    std::string temp_prefix = make_temp_prefix(out_prefix);
    TempPrefixGuard guard(temp_prefix);
    MergeState state;
    {
        OutputFiles output(temp_prefix);
        for (size_t i = 0; i < inputs.size(); ++i) {
            std::cerr << "FELIXla concat merge: [" << (i + 1) << '/' << inputs.size()
                      << "] " << inputs[i].prefix << '\n';
            PrefixSummary added = process_prefix(inputs[i], &output, state);
            std::cerr << "  variants=" << added.counts.global
                      << " common=" << added.counts.common
                      << " rare=" << added.counts.rare
                      << " ancestry_blocks=" << added.counts.blocks << '\n';
        }
        output.close();
    }

    write_samples_file(temp_prefix + ".samples", samples);
    write_meta_file(temp_prefix + ".meta", inputs.front().meta, list_path, inputs, state);
    preflight_prefix_files(temp_prefix);
    publish_temp_prefix(temp_prefix, out_prefix);
    guard.commit();

    std::cerr << "Finished FELIXla concatenation.\n"
              << "Input prefixes:        " << inputs.size() << '\n'
              << "Global variants:       " << state.global_variants << '\n'
              << "Common variants:       " << state.common_variants << '\n'
              << "Rare variants:         " << state.rare_variants << '\n'
              << "Ancestry blocks:       " << state.ancestry_blocks << '\n';
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }
    try {
        return run_concat(argv[1], argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
