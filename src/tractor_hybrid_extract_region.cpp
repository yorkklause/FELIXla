// tractor_hybrid_extract_region.cpp
// designed by Kai, implemented by codex
//
// Select part of an existing FELIXla packed prefix into a new one: a 1-based
// inclusive genomic region, a PVAR/VCF list of alleles to keep or drop, or
// both. Ancestry blocks describe local ancestry along the coordinate axis, not
// per variant, so an allele-level selection passes them through untouched and
// only a region clips them.

#include <algorithm>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace {

struct Meta {
    uint64_t n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t n_ancestries = 0;
    std::vector<std::string> lines;
};

// --region is half-open [START, END): the first base is in, the last is not.
// `end` here is the last base actually covered, the interval closed, because
// the membership test and the block clip below are written that way; `label`
// is the half-open text, which is what goes into the metadata.
struct Region {
    std::string chr;
    int64_t start = 1;
    int64_t end = 0;
    std::string label;

    // An empty contig means the whole prefix, so a run that filters only by
    // site list needs no region and records none.
    bool whole() const { return chr.empty(); }
};

// One listed coordinate: its REF and the ALT alleles named there, sorted so a
// record's ALT can be found without scanning.
struct SitePosition {
    std::string ref;
    std::vector<std::string> alts;
    uint64_t line_no = 0;
};

struct NormalizedAllele {
    int64_t pos = 0;
    std::string_view ref;
    std::string_view alt;
};

struct VariantRecord {
    bool common = false;
    uint64_t local_index = 0;
    uint32_t global_variant_index = 0;
    std::string chr;
    int64_t pos = 0;
    std::string id;
    std::string ref;
    std::string alt;
    uint32_t alt_index = 0;
    uint32_t block_id = UINT32_MAX;
    uint64_t payload_offset = 0;
    uint32_t n_carriers = 0;
    uint32_t mac = 0;
};

struct AncBlockRecord {
    uint32_t old_block_id = 0;
    uint32_t new_block_id = 0;
    std::string chr;
    int64_t start = 0;
    int64_t end = 0;
    uint64_t anc_offset = 0;
};

static constexpr uint32_t kPackedHapMask = 0x07ffffffu;

[[noreturn]] static void die(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::fputs("ERROR: ", stderr);
    std::vfprintf(stderr, fmt, args);
    std::fputc('\n', stderr);
    va_end(args);
    std::exit(1);
}

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;

    while (true) {
        size_t tab = line.find('\t', start);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }

        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }

    return fields;
}

static std::string strip_commas(std::string value) {
    value.erase(std::remove(value.begin(), value.end(), ','), value.end());
    return value;
}

static int64_t parse_i64_string(const std::string& value, const char* field) {
    std::string clean = strip_commas(value);
    char* end = nullptr;
    long long parsed = std::strtoll(clean.c_str(), &end, 10);
    if (!end || *end != '\0') {
        die("invalid integer for %s: %s", field, value.c_str());
    }
    return static_cast<int64_t>(parsed);
}

static uint64_t parse_u64(const std::string& value, const char* field) {
    char* end = nullptr;
    unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
    if (!end || *end != '\0') {
        die("invalid unsigned integer for %s: %s", field, value.c_str());
    }
    return static_cast<uint64_t>(parsed);
}

static FILE* open_file_or_die(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) die("cannot open %s", path.c_str());

    if (setvbuf(fp, nullptr, _IOFBF, 1 << 24) != 0) {
        die("setvbuf failed for %s", path.c_str());
    }

    return fp;
}

static uint64_t tell_or_die(FILE* fp, const char* path) {
    off_t offset = ftello(fp);
    if (offset < 0) die("ftello failed for %s", path);
    return static_cast<uint64_t>(offset);
}

static void seek_or_die(FILE* fp, uint64_t offset, const char* path) {
    if (offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max())) {
        die("offset too large for this platform in %s", path);
    }

    if (fseeko(fp, static_cast<off_t>(offset), SEEK_SET) != 0) {
        die("seek failed in %s", path);
    }
}

static void read_exact(FILE* fp, void* data, size_t bytes, const char* path) {
    if (bytes == 0) return;
    if (std::fread(data, 1, bytes, fp) != bytes) {
        die("unexpected EOF while reading %s", path);
    }
}

static void write_exact(FILE* fp, const void* data, size_t bytes, const char* what) {
    if (bytes == 0) return;
    if (std::fwrite(data, 1, bytes, fp) != bytes) {
        die("failed writing %s", what);
    }
}

static uint32_t read_u32_le(FILE* fp, const char* path) {
    uint8_t b[4];
    read_exact(fp, b, sizeof(b), path);
    return static_cast<uint32_t>(b[0]) |
           (static_cast<uint32_t>(b[1]) << 8) |
           (static_cast<uint32_t>(b[2]) << 16) |
           (static_cast<uint32_t>(b[3]) << 24);
}

static uint64_t read_u64_le(FILE* fp, const char* path) {
    uint8_t b[8];
    read_exact(fp, b, sizeof(b), path);
    return static_cast<uint64_t>(b[0]) |
           (static_cast<uint64_t>(b[1]) << 8) |
           (static_cast<uint64_t>(b[2]) << 16) |
           (static_cast<uint64_t>(b[3]) << 24) |
           (static_cast<uint64_t>(b[4]) << 32) |
           (static_cast<uint64_t>(b[5]) << 40) |
           (static_cast<uint64_t>(b[6]) << 48) |
           (static_cast<uint64_t>(b[7]) << 56);
}

static int64_t read_i64_le(FILE* fp, const char* path) {
    return static_cast<int64_t>(read_u64_le(fp, path));
}

static std::string read_string(FILE* fp, const char* path) {
    uint32_t len = read_u32_le(fp, path);
    std::string value(len, '\0');
    if (len > 0) read_exact(fp, value.data(), len, path);
    return value;
}

static void write_u32_le(FILE* fp, uint32_t value, const char* what) {
    uint8_t b[4] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24)
    };
    write_exact(fp, b, sizeof(b), what);
}

static void write_u64_le(FILE* fp, uint64_t value, const char* what) {
    uint8_t b[8] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24),
        static_cast<uint8_t>(value >> 32),
        static_cast<uint8_t>(value >> 40),
        static_cast<uint8_t>(value >> 48),
        static_cast<uint8_t>(value >> 56)
    };
    write_exact(fp, b, sizeof(b), what);
}

static void write_i64_le(FILE* fp, int64_t value, const char* what) {
    write_u64_le(fp, static_cast<uint64_t>(value), what);
}

static void write_string(FILE* fp, const std::string& value, const char* what) {
    if (value.size() > std::numeric_limits<uint32_t>::max()) {
        die("string too long while writing %s", what);
    }

    write_u32_le(fp, static_cast<uint32_t>(value.size()), what);
    write_exact(fp, value.data(), value.size(), what);
}

static void read_magic(FILE* fp, const std::string& path, const char magic[8]) {
    char got[8];
    read_exact(fp, got, sizeof(got), path.c_str());
    if (std::memcmp(got, magic, sizeof(got)) != 0) {
        die("bad magic in %s", path.c_str());
    }
}

static void write_magic(FILE* fp, const char magic[8]) {
    write_exact(fp, magic, 8, "magic header");
}

static bool try_record_start(FILE* fp, const std::string& path) {
    int first = std::fgetc(fp);
    if (first == EOF) return false;
    if (std::ungetc(first, fp) == EOF) die("ungetc failed for %s", path.c_str());
    return true;
}

static bool chrom_matches(const std::string& observed, const std::string& requested) {
    if (observed == requested) return true;

    if (observed.size() > 3 && observed.compare(0, 3, "chr") == 0) {
        return observed.substr(3) == requested;
    }

    if (requested.size() > 3 && requested.compare(0, 3, "chr") == 0) {
        return requested.substr(3) == observed;
    }

    return false;
}

static bool in_region(const std::string& chr, int64_t pos, const Region& region) {
    if (region.whole()) return true;
    return chrom_matches(chr, region.chr) && pos >= region.start && pos <= region.end;
}

// A list may spell a contig with or without the chr prefix. There is no VCF
// header here to resolve against, so both the list and the prefix's own
// markers are keyed on the name with the prefix removed, which collapses the
// two spellings onto one key without preferring either.
static std::string contig_key(const std::string& chr) {
    if (chr.size() > 3 && chr.compare(0, 3, "chr") == 0) return chr.substr(3);
    return chr;
}

// A one-column line is a gnomAD-style variant ID: CHROM, POS, REF and ALT
// joined by ':' or '-'. Split from the right, three separators back, because
// REF and ALT never contain one but a contig name can -- HLA-A*01:01-100-A-T
// has to come apart as HLA-A*01:01 / 100 / A / T.
static bool split_variant_id(const std::string& line, std::vector<std::string>& fields) {
    size_t cut[3];
    size_t found = 0;
    for (size_t i = line.size(); i-- > 0 && found < 3;) {
        if (line[i] == ':' || line[i] == '-') cut[found++] = i;
    }
    if (found < 3) return false;
    if (cut[0] == line.size() - 1 || cut[1] + 1 == cut[0] || cut[2] + 1 == cut[1] ||
        cut[2] == 0) {
        return false;
    }

    fields.assign(5, std::string());
    fields[0] = line.substr(0, cut[2]);
    fields[1] = line.substr(cut[2] + 1, cut[1] - cut[2] - 1);
    fields[3] = line.substr(cut[1] + 1, cut[0] - cut[1] - 1);
    fields[4] = line.substr(cut[0] + 1);
    return true;
}

static bool is_symbolic_allele(const std::string& value) {
    return value.empty() || value[0] == '<' || value == "." ||
           value.find_first_of("[]") != std::string::npos;
}

// The same normalization the packer applies when a site list writes an allele
// with different padding from the genotype VCF: strip shared trailing bases,
// then shared leading ones, keeping at least one base on each side. Stripping
// from the front moves the coordinate, so the position travels with the pair.
static NormalizedAllele normalize_allele(int64_t pos, const std::string& ref,
                                         const std::string& alt) {
    std::string_view r(ref);
    std::string_view a(alt);
    if (is_symbolic_allele(ref) || is_symbolic_allele(alt)) {
        return NormalizedAllele{pos, r, a};
    }
    while (r.size() > 1 && a.size() > 1 && r.back() == a.back()) {
        r.remove_suffix(1);
        a.remove_suffix(1);
    }
    while (r.size() > 1 && a.size() > 1 && r.front() == a.front()) {
        r.remove_prefix(1);
        a.remove_prefix(1);
        ++pos;
    }
    return NormalizedAllele{pos, r, a};
}

// A PVAR/VCF-style list of alleles, matched against the markers of a packed
// prefix exactly as the packer matches it against genotype records: by CHROM
// and POS, then by ALT where the REF agrees, and otherwise by comparing both
// sides normalized.
class SiteList {
public:
    void load(const std::string& path, const char* flag, bool ref_mismatch_fatal) {
        flag_ = flag;
        path_ = path;
        ref_mismatch_fatal_ = ref_mismatch_fatal;
        active_ = true;

        std::ifstream in(path);
        if (!in) die("cannot open %s list %s", flag, path.c_str());

        std::string line;
        uint64_t line_no = 0;
        while (std::getline(in, line)) {
            ++line_no;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#') continue;

            std::vector<std::string> fields;
            if (line.find('\t') == std::string::npos) {
                if (!split_variant_id(line, fields)) {
                    die("%s:%llu is neither PVAR/VCF columns CHROM POS ID REF ALT "
                        "nor a CHROM:POS:REF:ALT variant ID",
                        path.c_str(), static_cast<unsigned long long>(line_no));
                }
            } else {
                fields = split_tabs(line);
                if (fields.size() < 5) {
                    die("%s:%llu has %zu columns; CHROM POS ID REF ALT are required",
                        path.c_str(), static_cast<unsigned long long>(line_no),
                        fields.size());
                }
            }
            const std::string& chr = fields[0];
            int64_t pos = parse_i64_string(fields[1], "site list POS");
            const std::string& ref = fields[3];
            if (chr.empty()) die("%s:%llu has an empty CHROM", path.c_str(),
                                 static_cast<unsigned long long>(line_no));
            if (pos <= 0) die("%s:%llu has a non-positive POS", path.c_str(),
                              static_cast<unsigned long long>(line_no));
            if (ref.empty() || ref == ".") {
                die("%s:%llu has an unknown REF", path.c_str(),
                    static_cast<unsigned long long>(line_no));
            }

            SitePosition& position = positions_[Key{contig_key(chr), pos}];
            if (position.ref.empty()) {
                position.ref = ref;
                position.line_no = line_no;
            } else if (position.ref != ref) {
                die("conflicting REF values in %s list %s at %s:%lld: %s vs %s",
                    flag, path.c_str(), chr.c_str(), static_cast<long long>(pos),
                    position.ref.c_str(), ref.c_str());
            }

            size_t start = 0;
            while (start <= fields[4].size()) {
                size_t comma = fields[4].find(',', start);
                std::string alt = fields[4].substr(
                    start, comma == std::string::npos ? std::string::npos : comma - start);
                if (alt.empty() || alt == ".") {
                    die("%s:%llu has an empty ALT", path.c_str(),
                        static_cast<unsigned long long>(line_no));
                }
                auto at = std::lower_bound(position.alts.begin(), position.alts.end(), alt);
                if (at != position.alts.end() && *at == alt) {
                    die("duplicate allele in %s list %s at %s:%lld: %s>%s",
                        flag, path.c_str(), chr.c_str(), static_cast<long long>(pos),
                        ref.c_str(), alt.c_str());
                }
                position.alts.insert(at, alt);
                ++allele_count_;
                if (comma == std::string::npos) break;
                start = comma + 1;
            }
        }
        if (!in.eof()) die("failed reading %s list %s", flag, path.c_str());
        // A list naming nothing is a mistake worth stopping for rather than a
        // filter that keeps everything or drops nothing, and the packer
        // refuses it the same way.
        if (positions_.empty()) {
            die("%s site list is empty: %s", flag, path.c_str());
        }
    }

    bool active() const { return active_; }
    const std::string& path() const { return path_; }
    size_t position_count() const { return positions_.size(); }
    uint64_t allele_count() const { return allele_count_; }

    // True when the list names this marker's allele. Also accumulates the
    // per-coordinate evidence that --extract needs in order to refuse a list
    // whose REF disagrees with the prefix.
    bool matches(const std::string& chr, int64_t pos, const std::string& ref,
                 const std::string& alt) {
        auto it = positions_.find(Key{contig_key(chr), pos});
        if (it == positions_.end()) return false;
        const SitePosition& position = it->second;

        note_coordinate(chr, pos, ref, position);

        bool raw_ref = position.ref == ref;
        if (raw_ref) {
            auto at = std::lower_bound(position.alts.begin(), position.alts.end(), alt);
            bool exact = at != position.alts.end() && *at == alt;
            normalized_matched_ = normalized_matched_ || exact;
            return exact;
        }

        NormalizedAllele query = normalize_allele(pos, ref, alt);
        for (const std::string& candidate_alt : position.alts) {
            NormalizedAllele candidate = normalize_allele(pos, position.ref, candidate_alt);
            if (candidate.pos == query.pos && candidate.ref == query.ref &&
                candidate.alt == query.alt) {
                normalized_matched_ = true;
                return true;
            }
        }
        return false;
    }

    // Called once the marker stream has left a coordinate, and once at the end.
    void finish_coordinate() {
        if (!open_) return;
        if (ref_mismatch_fatal_ && !raw_ref_matched_ && !normalized_matched_) {
            die("REF mismatch for %s site %s:%lld: list has %s, prefix has %s",
                flag_, open_chr_.c_str(), static_cast<long long>(open_pos_),
                open_list_ref_.c_str(), first_prefix_ref_.c_str());
        }
        open_ = false;
        raw_ref_matched_ = false;
        normalized_matched_ = false;
        first_prefix_ref_.clear();
    }

private:
    struct Key {
        std::string chr;
        int64_t pos;
        bool operator==(const Key& other) const {
            return pos == other.pos && chr == other.chr;
        }
    };
    struct KeyHash {
        size_t operator()(const Key& key) const {
            return std::hash<std::string>()(key.chr) ^
                   (std::hash<int64_t>()(key.pos) * 1099511628211ULL);
        }
    };

    void note_coordinate(const std::string& chr, int64_t pos, const std::string& ref,
                         const SitePosition& position) {
        if (open_ && (open_pos_ != pos || open_chr_ != chr)) finish_coordinate();
        if (!open_) {
            open_ = true;
            open_chr_ = chr;
            open_pos_ = pos;
            open_list_ref_ = position.ref;
        }
        bool raw_ref = position.ref == ref;
        if (!raw_ref && first_prefix_ref_.empty()) first_prefix_ref_ = ref;
        raw_ref_matched_ = raw_ref_matched_ || raw_ref;
    }

    bool active_ = false;
    bool ref_mismatch_fatal_ = false;
    const char* flag_ = "--extract";
    std::string path_;
    std::unordered_map<Key, SitePosition, KeyHash> positions_;
    uint64_t allele_count_ = 0;

    bool open_ = false;
    std::string open_chr_;
    int64_t open_pos_ = 0;
    std::string open_list_ref_;
    std::string first_prefix_ref_;
    bool raw_ref_matched_ = false;
    bool normalized_matched_ = false;
};

static bool overlaps_region(const AncBlockRecord& block, const Region& region) {
    if (region.whole()) return true;
    return chrom_matches(block.chr, region.chr) &&
           block.end >= region.start &&
           block.start <= region.end;
}

static bool ancestry_blocks_cover_position(
    const std::vector<AncBlockRecord>& blocks,
    const std::string& chr,
    int64_t pos
) {
    for (const AncBlockRecord& block : blocks) {
        if (chrom_matches(block.chr, chr) && pos >= block.start && pos <= block.end) {
            return true;
        }
    }
    return false;
}

static Region make_region(const std::string& chr, int64_t start, int64_t end_exclusive,
                          const std::string& source) {
    if (chr.empty()) die("region chromosome is empty");
    if (start <= 0) die("region START must be positive: %s", source.c_str());
    if (end_exclusive <= start) {
        die("region END must be greater than START: regions are half-open "
            "[START, END), so %s covers nothing", source.c_str());
    }
    Region region;
    region.chr = chr;
    region.start = start;
    region.end = end_exclusive - 1;
    region.label = chr + ":" + std::to_string(start) + "-" + std::to_string(end_exclusive);
    return region;
}

static Region parse_region_string(const std::string& region_text) {
    size_t colon = region_text.find(':');
    size_t dash = region_text.find('-', colon == std::string::npos ? 0 : colon + 1);
    if (colon == std::string::npos || dash == std::string::npos || dash <= colon + 1) {
        die("region must look like chr:start-end: %s", region_text.c_str());
    }

    return make_region(
        region_text.substr(0, colon),
        parse_i64_string(region_text.substr(colon + 1, dash - colon - 1), "region start"),
        parse_i64_string(region_text.substr(dash + 1), "region end"),
        region_text
    );
}

static Meta read_meta(const std::string& path) {
    std::ifstream in(path);
    if (!in) die("cannot open %s", path.c_str());

    std::unordered_map<std::string, std::string> kv;
    Meta meta;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        meta.lines.push_back(line);
        std::vector<std::string> fields = split_tabs(line);
        if (fields.size() >= 2) kv[fields[0]] = fields[1];
    }

    if (!kv.count("n_samples")) die("missing n_samples in %s", path.c_str());
    if (!kv.count("n_haps")) die("missing n_haps in %s", path.c_str());
    if (!kv.count("n_words")) die("missing n_words in %s", path.c_str());
    if (!kv.count("n_ancestries")) die("missing n_ancestries in %s", path.c_str());

    meta.n_samples = parse_u64(kv["n_samples"], "n_samples");
    meta.n_haps = parse_u64(kv["n_haps"], "n_haps");
    meta.n_words = parse_u64(kv["n_words"], "n_words");
    meta.n_ancestries = parse_u64(kv["n_ancestries"], "n_ancestries");

    if (meta.n_samples == 0) die("n_samples must be positive");
    if (meta.n_haps != meta.n_samples * 2) die("n_haps must equal 2 * n_samples");
    if (meta.n_words != (meta.n_haps + 63) / 64) die("n_words does not match n_haps");
    if (meta.n_ancestries == 0 || meta.n_ancestries > 32) die("n_ancestries must be in [1, 32]");

    return meta;
}

static void copy_text_file(const std::string& in_path, const std::string& out_path) {
    std::ifstream in(in_path, std::ios::binary);
    if (!in) die("cannot open %s", in_path.c_str());
    std::ofstream out(out_path, std::ios::binary);
    if (!out) die("cannot open %s", out_path.c_str());
    out << in.rdbuf();
    if (!out) die("failed writing %s", out_path.c_str());
}

static void write_meta(
    const Meta& meta,
    const std::string& out_path,
    const std::string& in_prefix,
    const Region& region,
    uint64_t global_variants,
    uint64_t common_variants,
    uint64_t rare_variants,
    uint64_t ancestry_blocks
) {
    std::ofstream out(out_path);
    if (!out) die("cannot open %s", out_path.c_str());

    for (const std::string& line : meta.lines) {
        size_t tab = line.find('\t');
        std::string key = line.substr(0, tab);
        if (key == "global_variants" || key == "common_variants" ||
            key == "rare_variants" || key == "ancestry_blocks") {
            continue;
        }
        out << line << '\n';
    }
    out << "source_hybrid_prefix\t" << in_prefix << '\n';
    if (!region.whole()) {
        out << "extracted_region\t" << region.label << '\n';
        out << "region_convention\thalf-open\n";
    }
    out << "global_variants\t" << global_variants << '\n';
    out << "common_variants\t" << common_variants << '\n';
    out << "rare_variants\t" << rare_variants << '\n';
    out << "ancestry_blocks\t" << ancestry_blocks << '\n';
    if (!out) die("failed writing %s", out_path.c_str());
}

static void read_common_record(
    FILE* fp,
    const std::string& path,
    VariantRecord& record
) {
    record.common = true;
    record.local_index = read_u64_le(fp, path.c_str());
    record.global_variant_index = read_u32_le(fp, path.c_str());
    record.chr = read_string(fp, path.c_str());
    record.pos = read_i64_le(fp, path.c_str());
    record.id = read_string(fp, path.c_str());
    record.ref = read_string(fp, path.c_str());
    record.alt = read_string(fp, path.c_str());
    record.alt_index = read_u32_le(fp, path.c_str());
    record.block_id = read_u32_le(fp, path.c_str());
    record.payload_offset = read_u64_le(fp, path.c_str());
    record.mac = read_u32_le(fp, path.c_str());
}

static void read_rare_record(
    FILE* fp,
    const std::string& path,
    VariantRecord& record
) {
    record.common = false;
    record.local_index = read_u64_le(fp, path.c_str());
    record.global_variant_index = read_u32_le(fp, path.c_str());
    record.chr = read_string(fp, path.c_str());
    record.pos = read_i64_le(fp, path.c_str());
    record.id = read_string(fp, path.c_str());
    record.ref = read_string(fp, path.c_str());
    record.alt = read_string(fp, path.c_str());
    record.alt_index = read_u32_le(fp, path.c_str());
    record.payload_offset = read_u64_le(fp, path.c_str());
    record.n_carriers = read_u32_le(fp, path.c_str());
    record.mac = read_u32_le(fp, path.c_str());
}

static AncBlockRecord read_anc_record(FILE* fp, const std::string& path) {
    AncBlockRecord record;
    record.old_block_id = read_u32_le(fp, path.c_str());
    record.chr = read_string(fp, path.c_str());
    record.start = read_i64_le(fp, path.c_str());
    record.end = read_i64_le(fp, path.c_str());
    record.anc_offset = read_u64_le(fp, path.c_str());
    return record;
}

static void write_common_record(
    FILE* fp,
    uint64_t local_index,
    uint32_t global_variant_index,
    const VariantRecord& record,
    uint32_t new_block_id,
    uint64_t geno_offset
) {
    write_u64_le(fp, local_index, "common mks");
    write_u32_le(fp, global_variant_index, "common mks");
    write_string(fp, record.chr, "common mks chr");
    write_i64_le(fp, record.pos, "common mks pos");
    write_string(fp, record.id, "common mks id");
    write_string(fp, record.ref, "common mks ref");
    write_string(fp, record.alt, "common mks alt");
    write_u32_le(fp, record.alt_index, "common mks");
    write_u32_le(fp, new_block_id, "common mks");
    write_u64_le(fp, geno_offset, "common mks");
    write_u32_le(fp, record.mac, "common mks");
}

static void write_rare_record(
    FILE* fp,
    uint64_t local_index,
    uint32_t global_variant_index,
    const VariantRecord& record,
    uint64_t carrier_offset
) {
    write_u64_le(fp, local_index, "rare mks");
    write_u32_le(fp, global_variant_index, "rare mks");
    write_string(fp, record.chr, "rare mks chr");
    write_i64_le(fp, record.pos, "rare mks pos");
    write_string(fp, record.id, "rare mks id");
    write_string(fp, record.ref, "rare mks ref");
    write_string(fp, record.alt, "rare mks alt");
    write_u32_le(fp, record.alt_index, "rare mks");
    write_u64_le(fp, carrier_offset, "rare mks");
    write_u32_le(fp, record.n_carriers, "rare mks");
    write_u32_le(fp, record.mac, "rare mks");
}

static void write_anc_record(FILE* fp, const AncBlockRecord& record, uint64_t anc_offset) {
    write_u32_le(fp, record.new_block_id, "ancestry mks");
    write_string(fp, record.chr, "ancestry mks chr");
    write_i64_le(fp, record.start, "ancestry mks start");
    write_i64_le(fp, record.end, "ancestry mks end");
    write_u64_le(fp, anc_offset, "ancestry mks offset");
}

static void copy_fixed_bytes(
    FILE* in_fp,
    const std::string& in_path,
    uint64_t offset,
    FILE* out_fp,
    const char* out_what,
    uint64_t bytes
) {
    seek_or_die(in_fp, offset, in_path.c_str());
    std::vector<uint8_t> buffer(static_cast<size_t>(bytes));
    read_exact(in_fp, buffer.data(), buffer.size(), in_path.c_str());
    write_exact(out_fp, buffer.data(), buffer.size(), out_what);
}

static void read_selected_anc_blocks(
    const std::string& in_prefix,
    const Meta& meta,
    const Region& region,
    std::vector<AncBlockRecord>& selected_blocks,
    std::unordered_map<uint32_t, uint32_t>& block_id_map
) {
    const std::string anc_mks_path = in_prefix + ".ancblock.mks";
    FILE* fp = open_file_or_die(anc_mks_path, "rb");
    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    read_magic(fp, anc_mks_path, anc_mks_magic);

    while (try_record_start(fp, anc_mks_path)) {
        AncBlockRecord block = read_anc_record(fp, anc_mks_path);
        if (!overlaps_region(block, region)) continue;

        // Only a region clips a block. An allele-level selection does not
        // change where a haplotype's ancestry switches, so the blocks it
        // leaves behind must describe the same spans they always did.
        if (!region.whole()) {
            block.start = std::max(block.start, region.start);
            block.end = std::min(block.end, region.end);
        }
        block.new_block_id = static_cast<uint32_t>(selected_blocks.size());

        if (block_id_map.count(block.old_block_id)) {
            die("duplicate ancestry block id in %s: %u", anc_mks_path.c_str(), block.old_block_id);
        }
        if (selected_blocks.size() > std::numeric_limits<uint32_t>::max()) {
            die("too many selected ancestry blocks");
        }
        block_id_map[block.old_block_id] = block.new_block_id;
        selected_blocks.push_back(block);
    }

    std::fclose(fp);

    if (meta.n_ancestries > 0 && meta.n_words > 0 && selected_blocks.empty()) {
        std::fprintf(
            stderr,
            "WARNING: no ancestry blocks overlap %s; output variants may be empty or unreadable if selected variants need ancestry blocks.\n",
            region.label.c_str()
        );
    }
}

static void write_selected_anc_blocks(
    const std::string& in_prefix,
    const std::string& out_prefix,
    const Meta& meta,
    const std::vector<AncBlockRecord>& selected_blocks
) {
    const std::string in_anc_bin_path = in_prefix + ".ancblock.bin";
    const std::string out_anc_bin_path = out_prefix + ".ancblock.bin";
    const std::string out_anc_mks_path = out_prefix + ".ancblock.mks";
    const std::string out_anc_idx_path = out_prefix + ".ancblock.idx";

    FILE* in_bin = open_file_or_die(in_anc_bin_path, "rb");
    FILE* out_bin = open_file_or_die(out_anc_bin_path, "wb");
    FILE* out_mks = open_file_or_die(out_anc_mks_path, "wb");
    FILE* out_idx = open_file_or_die(out_anc_idx_path, "wb");

    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    const char anc_idx_magic[8] = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};
    write_magic(out_mks, anc_mks_magic);
    write_magic(out_idx, anc_idx_magic);

    uint64_t block_bytes = meta.n_ancestries * meta.n_words * sizeof(uint64_t);
    for (const AncBlockRecord& block : selected_blocks) {
        uint64_t anc_offset = tell_or_die(out_bin, out_anc_bin_path.c_str());
        copy_fixed_bytes(in_bin, in_anc_bin_path, block.anc_offset, out_bin, "ancestry block", block_bytes);

        uint64_t mks_offset = tell_or_die(out_mks, out_anc_mks_path.c_str());
        write_anc_record(out_mks, block, anc_offset);

        write_u32_le(out_idx, block.new_block_id, "ancestry idx");
        write_u64_le(out_idx, mks_offset, "ancestry idx");
        write_u64_le(out_idx, anc_offset, "ancestry idx");
    }

    std::fclose(in_bin);
    std::fclose(out_bin);
    std::fclose(out_mks);
    std::fclose(out_idx);
}

static const VariantRecord* next_variant(
    FILE* common_mks_fp,
    const std::string& common_mks_path,
    FILE* rare_mks_fp,
    const std::string& rare_mks_path,
    bool& have_common,
    bool& have_rare,
    VariantRecord& common_record,
    VariantRecord& rare_record
) {
    if (!have_common && try_record_start(common_mks_fp, common_mks_path)) {
        read_common_record(common_mks_fp, common_mks_path, common_record);
        have_common = true;
    }

    if (!have_rare && try_record_start(rare_mks_fp, rare_mks_path)) {
        read_rare_record(rare_mks_fp, rare_mks_path, rare_record);
        have_rare = true;
    }

    if (!have_common && !have_rare) return nullptr;

    bool take_common = have_common && !have_rare;
    if (have_common && have_rare) {
        if (common_record.global_variant_index == rare_record.global_variant_index) {
            die("global variant %u exists in both common and rare streams", common_record.global_variant_index);
        }
        take_common = common_record.global_variant_index < rare_record.global_variant_index;
    }

    if (take_common) {
        have_common = false;
        return &common_record;
    }
    have_rare = false;
    return &rare_record;
}

static void write_selected_variants(
    const std::string& in_prefix,
    const std::string& out_prefix,
    const Meta& meta,
    const Region& region,
    SiteList& extract_sites,
    SiteList& exclude_sites,
    const std::unordered_map<uint32_t, uint32_t>& block_id_map,
    const std::vector<AncBlockRecord>& selected_blocks,
    uint64_t& common_written,
    uint64_t& rare_written,
    uint64_t& total_written
) {
    const std::string in_common_mks_path = in_prefix + ".common.variant.mks";
    const std::string in_rare_mks_path = in_prefix + ".rare.variant.mks";
    const std::string in_common_bin_path = in_prefix + ".common.geno.bin";
    const std::string in_rare_bin_path = in_prefix + ".rare.carrier.bin";
    const std::string out_common_mks_path = out_prefix + ".common.variant.mks";
    const std::string out_rare_mks_path = out_prefix + ".rare.variant.mks";
    const std::string out_common_idx_path = out_prefix + ".common.variant.idx";
    const std::string out_rare_idx_path = out_prefix + ".rare.variant.idx";
    const std::string out_common_bin_path = out_prefix + ".common.geno.bin";
    const std::string out_rare_bin_path = out_prefix + ".rare.carrier.bin";

    FILE* in_common_mks = open_file_or_die(in_common_mks_path, "rb");
    FILE* in_rare_mks = open_file_or_die(in_rare_mks_path, "rb");
    FILE* in_common_bin = open_file_or_die(in_common_bin_path, "rb");
    FILE* in_rare_bin = open_file_or_die(in_rare_bin_path, "rb");
    FILE* out_common_mks = open_file_or_die(out_common_mks_path, "wb");
    FILE* out_rare_mks = open_file_or_die(out_rare_mks_path, "wb");
    FILE* out_common_idx = open_file_or_die(out_common_idx_path, "wb");
    FILE* out_rare_idx = open_file_or_die(out_rare_idx_path, "wb");
    FILE* out_common_bin = open_file_or_die(out_common_bin_path, "wb");
    FILE* out_rare_bin = open_file_or_die(out_rare_bin_path, "wb");

    const char common_mks_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
    const char rare_mks_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
    const char common_idx_magic[8] = {'T', 'R', 'C', 'M', 'I', 'D', 'X', '2'};
    const char rare_idx_magic[8] = {'T', 'R', 'R', 'A', 'I', 'D', 'X', '2'};

    read_magic(in_common_mks, in_common_mks_path, common_mks_magic);
    read_magic(in_rare_mks, in_rare_mks_path, rare_mks_magic);
    write_magic(out_common_mks, common_mks_magic);
    write_magic(out_rare_mks, rare_mks_magic);
    write_magic(out_common_idx, common_idx_magic);
    write_magic(out_rare_idx, rare_idx_magic);

    bool have_common = false;
    bool have_rare = false;
    VariantRecord common_record;
    VariantRecord rare_record;

    uint64_t common_bytes = meta.n_words * sizeof(uint64_t);
    std::vector<uint8_t> common_payload(static_cast<size_t>(common_bytes));

    while (true) {
        const VariantRecord* next = next_variant(
            in_common_mks,
            in_common_mks_path,
            in_rare_mks,
            in_rare_mks_path,
            have_common,
            have_rare,
            common_record,
            rare_record
        );
        if (!next) break;
        const VariantRecord& record = *next;
        if (!in_region(record.chr, record.pos, region)) continue;

        // --extract selects, --exclude then removes, so the two compose the
        // way they do in PLINK and the way they do when packing.
        if (extract_sites.active() &&
            !extract_sites.matches(record.chr, record.pos, record.ref, record.alt)) {
            continue;
        }
        if (exclude_sites.active() &&
            exclude_sites.matches(record.chr, record.pos, record.ref, record.alt)) {
            continue;
        }

        if (total_written > std::numeric_limits<uint32_t>::max()) {
            die("selected variant count exceeds uint32_t limit");
        }
        uint32_t new_global_index = static_cast<uint32_t>(total_written);

        if (record.common) {
            auto block_it = block_id_map.find(record.block_id);
            if (block_it == block_id_map.end()) {
                die(
                    "common variant %s:%lld references ancestry block %u outside selected region",
                    record.chr.c_str(),
                    static_cast<long long>(record.pos),
                    record.block_id
                );
            }

            seek_or_die(in_common_bin, record.payload_offset, in_common_bin_path.c_str());
            read_exact(in_common_bin, common_payload.data(), common_payload.size(), in_common_bin_path.c_str());

            uint64_t geno_offset = tell_or_die(out_common_bin, out_common_bin_path.c_str());
            write_exact(out_common_bin, common_payload.data(), common_payload.size(), "common genotype payload");

            uint64_t mks_offset = tell_or_die(out_common_mks, out_common_mks_path.c_str());
            write_common_record(
                out_common_mks,
                common_written,
                new_global_index,
                record,
                block_it->second,
                geno_offset
            );
            write_u64_le(out_common_idx, common_written, "common idx");
            write_u32_le(out_common_idx, new_global_index, "common idx");
            write_u64_le(out_common_idx, mks_offset, "common idx");
            write_u64_le(out_common_idx, geno_offset, "common idx");

            ++common_written;
        } else {
            if (!ancestry_blocks_cover_position(selected_blocks, record.chr, record.pos)) {
                die(
                    "rare variant %s:%lld is not covered by any selected ancestry block",
                    record.chr.c_str(),
                    static_cast<long long>(record.pos)
                );
            }

            uint64_t carrier_offset = tell_or_die(out_rare_bin, out_rare_bin_path.c_str());
            seek_or_die(in_rare_bin, record.payload_offset, in_rare_bin_path.c_str());

            for (uint32_t i = 0; i < record.n_carriers; ++i) {
                uint32_t old_pos_index = read_u32_le(in_rare_bin, in_rare_bin_path.c_str());
                uint32_t anc_hap = read_u32_le(in_rare_bin, in_rare_bin_path.c_str());
                if (old_pos_index != record.global_variant_index) {
                    die(
                        "rare carrier pos_index mismatch: carrier has %u, marker has %u",
                        old_pos_index,
                        record.global_variant_index
                    );
                }
                uint32_t hap_id = anc_hap & kPackedHapMask;
                if (hap_id >= meta.n_haps) {
                    die("rare carrier hap_id %u exceeds n_haps", hap_id);
                }

                write_u32_le(out_rare_bin, new_global_index, "rare carrier");
                write_u32_le(out_rare_bin, anc_hap, "rare carrier");
            }

            uint64_t mks_offset = tell_or_die(out_rare_mks, out_rare_mks_path.c_str());
            write_rare_record(out_rare_mks, rare_written, new_global_index, record, carrier_offset);
            write_u64_le(out_rare_idx, rare_written, "rare idx");
            write_u32_le(out_rare_idx, new_global_index, "rare idx");
            write_u64_le(out_rare_idx, mks_offset, "rare idx");
            write_u64_le(out_rare_idx, carrier_offset, "rare idx");
            write_u32_le(out_rare_idx, record.n_carriers, "rare idx");

            ++rare_written;
        }

        ++total_written;
    }

    extract_sites.finish_coordinate();
    exclude_sites.finish_coordinate();

    std::fclose(in_common_mks);
    std::fclose(in_rare_mks);
    std::fclose(in_common_bin);
    std::fclose(in_rare_bin);
    std::fclose(out_common_mks);
    std::fclose(out_rare_mks);
    std::fclose(out_common_idx);
    std::fclose(out_rare_idx);
    std::fclose(out_common_bin);
    std::fclose(out_rare_bin);
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s in_prefix chr:start-end out_prefix [options]\n"
        "  %s in_prefix chr start end out_prefix [options]\n"
        "  %s in_prefix out_prefix --extract FILE | --exclude FILE\n\n"
        "Coordinates are 1-based and inclusive.\n\n"
        "Options:\n"
        "  --extract FILE   PVAR/VCF alleles to keep\n"
        "  --exclude FILE   PVAR/VCF alleles to drop, applied after --extract\n\n"
        "Example:\n"
        "  %s chr22 chr22:16000000-17000000 chr22.region\n",
        prog,
        prog,
        prog,
        prog
    );
}

} // namespace

int main(int argc, char** argv) {
    // The site-list flags may appear anywhere after the program name; what is
    // left is the historical positional form, with the region now optional.
    std::vector<std::string> positional;
    std::string extract_path;
    std::string exclude_path;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--extract" || arg == "--exclude") {
            if (i + 1 >= argc) die("%s requires a file", arg.c_str());
            std::string& target = arg == "--extract" ? extract_path : exclude_path;
            if (!target.empty()) die("%s was given twice", arg.c_str());
            target = argv[++i];
        } else if (arg.size() > 2 && arg[0] == '-' && arg[1] == '-') {
            print_usage(argv[0]);
            die("unknown option: %s", arg.c_str());
        } else {
            positional.push_back(arg);
        }
    }

    bool have_sites = !extract_path.empty() || !exclude_path.empty();
    if (positional.size() != 2 && positional.size() != 3 && positional.size() != 5) {
        print_usage(argv[0]);
        return 1;
    }
    if (positional.size() == 2 && !have_sites) {
        print_usage(argv[0]);
        die("give a region, a site list, or both");
    }

    std::string in_prefix = positional.front();
    std::string out_prefix = positional.back();
    Region region;
    if (positional.size() == 3) {
        region = parse_region_string(positional[1]);
    } else if (positional.size() == 5) {
        region = make_region(
            positional[1],
            parse_i64_string(positional[2], "region start"),
            parse_i64_string(positional[3], "region end"),
            positional[1] + " " + positional[2] + " " + positional[3]
        );
    } else {
        region.label = "whole prefix";
    }

    SiteList extract_sites;
    SiteList exclude_sites;
    // A listed coordinate the prefix holds under a different REF is a mistake
    // worth stopping for when selecting, and inert when dropping: an --exclude
    // entry that matches nothing simply removes nothing.
    if (!extract_path.empty()) extract_sites.load(extract_path, "--extract", true);
    if (!exclude_path.empty()) exclude_sites.load(exclude_path, "--exclude", false);

    Meta meta = read_meta(in_prefix + ".meta");
    copy_text_file(in_prefix + ".samples", out_prefix + ".samples");

    std::vector<AncBlockRecord> selected_blocks;
    std::unordered_map<uint32_t, uint32_t> block_id_map;
    read_selected_anc_blocks(in_prefix, meta, region, selected_blocks, block_id_map);
    write_selected_anc_blocks(in_prefix, out_prefix, meta, selected_blocks);

    uint64_t common_written = 0;
    uint64_t rare_written = 0;
    uint64_t total_written = 0;
    write_selected_variants(
        in_prefix,
        out_prefix,
        meta,
        region,
        extract_sites,
        exclude_sites,
        block_id_map,
        selected_blocks,
        common_written,
        rare_written,
        total_written
    );
    write_meta(
        meta,
        out_prefix + ".meta",
        in_prefix,
        region,
        total_written,
        common_written,
        rare_written,
        selected_blocks.size()
    );

    std::fprintf(stderr, "Finished.\n");
    std::fprintf(stderr, "Region:                %s\n", region.label.c_str());
    if (extract_sites.active()) {
        std::fprintf(stderr, "Extract list:          %s (%llu alleles)\n",
            extract_sites.path().c_str(),
            static_cast<unsigned long long>(extract_sites.allele_count()));
    }
    if (exclude_sites.active()) {
        std::fprintf(stderr, "Exclude list:          %s (%llu alleles)\n",
            exclude_sites.path().c_str(),
            static_cast<unsigned long long>(exclude_sites.allele_count()));
    }
    std::fprintf(stderr, "Global variants:       %llu\n", static_cast<unsigned long long>(total_written));
    std::fprintf(stderr, "Common variants:       %llu\n", static_cast<unsigned long long>(common_written));
    std::fprintf(stderr, "Rare variants:         %llu\n", static_cast<unsigned long long>(rare_written));
    std::fprintf(stderr, "Ancestry blocks:       %llu\n", static_cast<unsigned long long>(selected_blocks.size()));

    return 0;
}
