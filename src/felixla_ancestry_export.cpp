// felixla_ancestry_export.cpp
// designed by Kai, implemented by claude
//
// The three read-only views of a packed prefix's ancestry blocks:
//
//   lai               <out>.lai.gz, one row per block, one column per
//                     haplotype, holding that haplotype's ancestry code.
//   global-admixture  <out>.global.admixture.tsv, one row per sample, holding
//                     the base-pair weighted share of each ancestry.
//   local-admixture   <out>.local.admixture.tsv, one row per block, holding
//                     each ancestry's share of the cohort's haplotypes there.
//
// Ancestry codes are written 1-based: the packed files number ancestries from
// zero, but zero is reserved in the LAI file for a haplotype no caller labelled.
// Block coordinates are 1-based inclusive, as they are stored.

#include <htslib/bgzf.h>

#include <algorithm>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

enum class Mode { Lai, GlobalAdmixture, LocalAdmixture };

struct Meta {
    uint64_t n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t n_ancestries = 0;
    // Empty where the input did not name that ancestry; index is the 0-based
    // code, so ancestry_names[0] is the name meta calls ancestry_name_1.
    std::vector<std::string> ancestry_names;
};

struct Block {
    uint32_t id = 0;
    std::string chr;
    int64_t start = 0;
    int64_t end = 0;
};

[[noreturn]] void die(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::fputs("ERROR: ", stderr);
    std::vfprintf(stderr, fmt, args);
    std::fputc('\n', stderr);
    va_end(args);
    std::exit(1);
}

uint64_t parse_u64(const std::string& value, const char* field) {
    char* end = nullptr;
    unsigned long long parsed = std::strtoull(value.c_str(), &end, 10);
    if (!end || end == value.c_str() || *end != '\0') {
        die("invalid unsigned integer for %s: %s", field, value.c_str());
    }
    return static_cast<uint64_t>(parsed);
}

uint64_t load_u64_le(const uint8_t* bytes) {
    return static_cast<uint64_t>(bytes[0]) |
           (static_cast<uint64_t>(bytes[1]) << 8) |
           (static_cast<uint64_t>(bytes[2]) << 16) |
           (static_cast<uint64_t>(bytes[3]) << 24) |
           (static_cast<uint64_t>(bytes[4]) << 32) |
           (static_cast<uint64_t>(bytes[5]) << 40) |
           (static_cast<uint64_t>(bytes[6]) << 48) |
           (static_cast<uint64_t>(bytes[7]) << 56);
}

Meta read_meta(const std::string& path) {
    std::ifstream in(path);
    if (!in) die("cannot open %s", path.c_str());

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        kv[line.substr(0, tab)] = line.substr(tab + 1);
    }

    Meta meta;
    for (const char* key : {"n_samples", "n_haps", "n_words", "n_ancestries"}) {
        if (!kv.count(key)) die("missing %s in %s", key, path.c_str());
    }
    meta.n_samples = parse_u64(kv["n_samples"], "n_samples");
    meta.n_haps = parse_u64(kv["n_haps"], "n_haps");
    meta.n_words = parse_u64(kv["n_words"], "n_words");
    meta.n_ancestries = parse_u64(kv["n_ancestries"], "n_ancestries");

    if (meta.n_samples == 0) die("n_samples must be positive in %s", path.c_str());
    if (meta.n_haps != meta.n_samples * 2) die("n_haps must equal 2 * n_samples");
    if (meta.n_words != (meta.n_haps + 63) / 64) die("n_words does not match n_haps");
    if (meta.n_ancestries == 0 || meta.n_ancestries > 32) {
        die("n_ancestries must be in [1, 32]");
    }

    meta.ancestry_names.resize(static_cast<size_t>(meta.n_ancestries));
    for (uint64_t i = 0; i < meta.n_ancestries; ++i) {
        auto it = kv.find("ancestry_name_" + std::to_string(i + 1));
        if (it != kv.end()) meta.ancestry_names[static_cast<size_t>(i)] = it->second;
    }
    return meta;
}

// The column heading for an ancestry: its name where the caller supplied one,
// and the positional label otherwise, so the header always has n_ancestries
// distinct columns whether or not the input named them.
std::string ancestry_label(const Meta& meta, uint64_t index) {
    const std::string& name = meta.ancestry_names[static_cast<size_t>(index)];
    if (!name.empty()) return name;
    return "ANC" + std::to_string(index + 1);
}

std::vector<std::string> read_samples(const std::string& path, uint64_t expected_n) {
    std::ifstream in(path);
    if (!in) die("cannot open %s", path.c_str());

    std::vector<std::string> samples;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) samples.push_back(line);
    }
    if (samples.size() != expected_n) {
        die("sample count mismatch in %s: got %llu, expected %llu", path.c_str(),
            static_cast<unsigned long long>(samples.size()),
            static_cast<unsigned long long>(expected_n));
    }
    return samples;
}

FILE* open_or_die(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) die("cannot open %s", path.c_str());
    if (setvbuf(fp, nullptr, _IOFBF, 1 << 22) != 0) die("setvbuf failed for %s", path.c_str());
    return fp;
}

void read_exact(FILE* fp, void* data, size_t bytes, const char* path) {
    if (bytes == 0) return;
    if (std::fread(data, 1, bytes, fp) != bytes) die("unexpected EOF while reading %s", path);
}

uint32_t read_u32_le(FILE* fp, const char* path) {
    uint8_t bytes[4];
    read_exact(fp, bytes, sizeof(bytes), path);
    return static_cast<uint32_t>(bytes[0]) | (static_cast<uint32_t>(bytes[1]) << 8) |
           (static_cast<uint32_t>(bytes[2]) << 16) | (static_cast<uint32_t>(bytes[3]) << 24);
}

uint64_t read_u64_le(FILE* fp, const char* path) {
    uint8_t bytes[8];
    read_exact(fp, bytes, sizeof(bytes), path);
    return load_u64_le(bytes);
}

std::string read_string(FILE* fp, const char* path) {
    uint32_t len = read_u32_le(fp, path);
    std::string value(len, '\0');
    if (len > 0) read_exact(fp, value.data(), len, path);
    return value;
}

// Walks .ancblock.mks and .ancblock.bin together. The markers are written in
// payload order, so the payload is read sequentially and its recorded offset
// checked rather than seeked to -- a mismatch means the pair disagree.
class BlockReader {
public:
    BlockReader(const std::string& prefix, const Meta& meta)
        : mks_path_(prefix + ".ancblock.mks"),
          bin_path_(prefix + ".ancblock.bin"),
          meta_(meta),
          mks_(open_or_die(mks_path_, "rb")),
          bin_(open_or_die(bin_path_, "rb")) {
        char magic[8];
        read_exact(mks_, magic, sizeof(magic), mks_path_.c_str());
        if (std::memcmp(magic, "TRANMKS1", sizeof(magic)) != 0) {
            die("bad ancestry marker magic in %s", mks_path_.c_str());
        }
        block_words_ = meta_.n_ancestries * meta_.n_words;
        if (block_words_ > (uint64_t{1} << 40)) die("ancestry block is implausibly large");
        payload_.resize(static_cast<size_t>(block_words_ * 8));
        labels_.resize(static_cast<size_t>(meta_.n_haps));
    }

    ~BlockReader() {
        if (mks_) std::fclose(mks_);
        if (bin_) std::fclose(bin_);
    }

    // Reads the next block into current() and labels(); false at end of file.
    bool next() {
        int first = std::fgetc(mks_);
        if (first == EOF) {
            // The payload must end with the markers, or one of the pair is truncated.
            if (std::fgetc(bin_) != EOF) die("unread ancestry payload left in %s", bin_path_.c_str());
            return false;
        }
        if (std::ungetc(first, mks_) == EOF) die("ungetc failed for %s", mks_path_.c_str());

        current_ = Block{};
        current_.id = read_u32_le(mks_, mks_path_.c_str());
        current_.chr = read_string(mks_, mks_path_.c_str());
        current_.start = static_cast<int64_t>(read_u64_le(mks_, mks_path_.c_str()));
        current_.end = static_cast<int64_t>(read_u64_le(mks_, mks_path_.c_str()));
        uint64_t payload_offset = read_u64_le(mks_, mks_path_.c_str());

        if (current_.id != seen_) die("ancestry block ids are not consecutive in %s", mks_path_.c_str());
        if (current_.chr.empty() || current_.start <= 0 || current_.start > current_.end) {
            die("invalid ancestry block %u in %s", current_.id, mks_path_.c_str());
        }
        if (payload_offset != seen_ * payload_.size()) {
            die("ancestry block %u points at offset %llu in %s, expected %llu",
                current_.id, static_cast<unsigned long long>(payload_offset), bin_path_.c_str(),
                static_cast<unsigned long long>(seen_ * payload_.size()));
        }
        read_exact(bin_, payload_.data(), payload_.size(), bin_path_.c_str());
        decode();
        ++seen_;
        return true;
    }

    const Block& current() const { return current_; }
    // labels()[hap_id] is the 1-based ancestry code, or 0 where unlabelled.
    const std::vector<uint8_t>& labels() const { return labels_; }
    uint64_t seen() const { return seen_; }

private:
    void decode() {
        std::fill(labels_.begin(), labels_.end(), uint8_t{0});
        for (uint64_t ancestry = 0; ancestry < meta_.n_ancestries; ++ancestry) {
            const uint8_t* base = payload_.data() + static_cast<size_t>(ancestry * meta_.n_words * 8);
            for (uint64_t word = 0; word < meta_.n_words; ++word) {
                uint64_t bits = load_u64_le(base + static_cast<size_t>(word * 8));
                while (bits) {
                    uint64_t hap = word * 64 + static_cast<uint64_t>(__builtin_ctzll(bits));
                    bits &= bits - 1;
                    if (hap >= meta_.n_haps) {
                        die("ancestry block %u labels haplotype %llu past the last one",
                            current_.id, static_cast<unsigned long long>(hap));
                    }
                    if (labels_[static_cast<size_t>(hap)] != 0) {
                        die("ancestry block %u gives haplotype %llu two ancestries",
                            current_.id, static_cast<unsigned long long>(hap));
                    }
                    labels_[static_cast<size_t>(hap)] = static_cast<uint8_t>(ancestry + 1);
                }
            }
        }
    }

    std::string mks_path_;
    std::string bin_path_;
    Meta meta_;
    FILE* mks_ = nullptr;
    FILE* bin_ = nullptr;
    uint64_t block_words_ = 0;
    uint64_t seen_ = 0;
    Block current_;
    std::vector<uint8_t> payload_;
    std::vector<uint8_t> labels_;
};

// Appends a non-negative integer without going through printf, which shows up
// in a profile once there is a column per haplotype.
void append_u64(std::string& out, uint64_t value) {
    char digits[20];
    int n = 0;
    do {
        digits[n++] = static_cast<char>('0' + value % 10);
        value /= 10;
    } while (value);
    while (n > 0) out.push_back(digits[--n]);
}

void append_i64(std::string& out, int64_t value) {
    if (value < 0) {
        out.push_back('-');
        append_u64(out, static_cast<uint64_t>(-(value + 1)) + 1);
        return;
    }
    append_u64(out, static_cast<uint64_t>(value));
}

void append_proportion(std::string& out, uint64_t numerator, uint64_t denominator) {
    char buffer[32];
    double value = denominator ? static_cast<double>(numerator) / static_cast<double>(denominator) : 0.0;
    int n = std::snprintf(buffer, sizeof(buffer), "%.6f", value);
    if (n < 0 || static_cast<size_t>(n) >= sizeof(buffer)) die("failed formatting a proportion");
    out.append(buffer, static_cast<size_t>(n));
}

constexpr size_t kFlushBytes = 4u << 20;

class BgzfWriter {
public:
    BgzfWriter(const std::string& path, int threads) : path_(path) {
        fp_ = bgzf_open(path.c_str(), "w");
        if (!fp_) die("cannot write %s", path.c_str());
        if (threads > 1 && bgzf_mt(fp_, threads, 256) != 0) {
            die("cannot start %d compression threads for %s", threads, path.c_str());
        }
        buffer_.reserve(kFlushBytes + (1u << 16));
    }

    std::string& buffer() { return buffer_; }

    void maybe_flush() {
        if (buffer_.size() >= kFlushBytes) flush();
    }

    void flush() {
        if (buffer_.empty()) return;
        ssize_t wrote = bgzf_write(fp_, buffer_.data(), buffer_.size());
        if (wrote < 0 || static_cast<size_t>(wrote) != buffer_.size()) {
            die("failed writing %s", path_.c_str());
        }
        buffer_.clear();
    }

    void close() {
        flush();
        if (fp_ && bgzf_close(fp_) != 0) die("failed closing %s", path_.c_str());
        fp_ = nullptr;
    }

    ~BgzfWriter() {
        if (fp_) bgzf_close(fp_);
    }

private:
    std::string path_;
    BGZF* fp_ = nullptr;
    std::string buffer_;
};

void export_lai(const std::string& prefix, const std::string& out_path, const Meta& meta,
                const std::vector<std::string>& samples, int threads) {
    BgzfWriter writer(out_path, threads);
    std::string& out = writer.buffer();

    for (uint64_t i = 0; i < meta.n_ancestries; ++i) {
        out += "##ANC";
        append_u64(out, i + 1);
        out += " = ";
        out += ancestry_label(meta, i);
        out.push_back('\n');
    }
    out += "#CHR\tSTART\tEND";
    for (const std::string& sample : samples) {
        out.push_back('\t');
        out += sample;
        out += "_1\t";
        out += sample;
        out += "_2";
    }
    out.push_back('\n');

    BlockReader reader(prefix, meta);
    while (reader.next()) {
        const Block& block = reader.current();
        out += block.chr;
        out.push_back('\t');
        append_i64(out, block.start);
        out.push_back('\t');
        append_i64(out, block.end);
        const std::vector<uint8_t>& labels = reader.labels();
        for (uint64_t hap = 0; hap < meta.n_haps; ++hap) {
            out.push_back('\t');
            append_u64(out, labels[static_cast<size_t>(hap)]);
        }
        out.push_back('\n');
        writer.maybe_flush();
    }
    writer.close();

    std::fprintf(stderr, "Wrote %s: %llu blocks x %llu haplotypes.\n", out_path.c_str(),
                 static_cast<unsigned long long>(reader.seen()),
                 static_cast<unsigned long long>(meta.n_haps));
}

void export_global_admixture(const std::string& prefix, const std::string& out_path,
                             const Meta& meta, const std::vector<std::string>& samples) {
    // Base pairs, not blocks: a caller that splits a chromosome finely in one
    // place must not thereby weight that place more heavily.
    std::vector<uint64_t> bp(static_cast<size_t>(meta.n_samples * meta.n_ancestries), 0);
    uint64_t covered_bp = 0;

    BlockReader reader(prefix, meta);
    while (reader.next()) {
        const Block& block = reader.current();
        uint64_t length = static_cast<uint64_t>(block.end - block.start) + 1;
        covered_bp += length;
        const std::vector<uint8_t>& labels = reader.labels();
        for (uint64_t hap = 0; hap < meta.n_haps; ++hap) {
            uint8_t label = labels[static_cast<size_t>(hap)];
            if (label == 0) continue;
            size_t row = static_cast<size_t>((hap >> 1) * meta.n_ancestries);
            bp[row + label - 1] += length;
        }
    }
    if (covered_bp == 0) die("the prefix has no ancestry blocks to summarise");

    FILE* fp = open_or_die(out_path, "w");
    std::string out = "#ID";
    for (uint64_t i = 0; i < meta.n_ancestries; ++i) {
        out.push_back('\t');
        out += ancestry_label(meta, i);
    }
    out.push_back('\n');
    // Both haplotypes of a sample are covered by every block, so the row sums
    // to one exactly where no haplotype went unlabelled.
    uint64_t per_sample_bp = covered_bp * 2;
    for (uint64_t s = 0; s < meta.n_samples; ++s) {
        out += samples[static_cast<size_t>(s)];
        for (uint64_t i = 0; i < meta.n_ancestries; ++i) {
            out.push_back('\t');
            append_proportion(out, bp[static_cast<size_t>(s * meta.n_ancestries + i)], per_sample_bp);
        }
        out.push_back('\n');
        if (out.size() >= kFlushBytes) {
            if (std::fwrite(out.data(), 1, out.size(), fp) != out.size()) {
                die("failed writing %s", out_path.c_str());
            }
            out.clear();
        }
    }
    if (!out.empty() && std::fwrite(out.data(), 1, out.size(), fp) != out.size()) {
        die("failed writing %s", out_path.c_str());
    }
    if (std::fclose(fp) != 0) die("failed closing %s", out_path.c_str());

    std::fprintf(stderr, "Wrote %s: %llu samples over %llu bp of ancestry blocks.\n",
                 out_path.c_str(), static_cast<unsigned long long>(meta.n_samples),
                 static_cast<unsigned long long>(covered_bp));
}

void export_local_admixture(const std::string& prefix, const std::string& out_path,
                            const Meta& meta) {
    FILE* fp = open_or_die(out_path, "w");
    std::string out = "#CHR\tSTART\tEND";
    for (uint64_t i = 0; i < meta.n_ancestries; ++i) {
        out.push_back('\t');
        out += ancestry_label(meta, i);
    }
    out.push_back('\n');

    std::vector<uint64_t> counts(static_cast<size_t>(meta.n_ancestries), 0);
    BlockReader reader(prefix, meta);
    while (reader.next()) {
        const Block& block = reader.current();
        std::fill(counts.begin(), counts.end(), uint64_t{0});
        const std::vector<uint8_t>& labels = reader.labels();
        for (uint64_t hap = 0; hap < meta.n_haps; ++hap) {
            uint8_t label = labels[static_cast<size_t>(hap)];
            if (label != 0) ++counts[label - 1];
        }
        out += block.chr;
        out.push_back('\t');
        append_i64(out, block.start);
        out.push_back('\t');
        append_i64(out, block.end);
        for (uint64_t i = 0; i < meta.n_ancestries; ++i) {
            out.push_back('\t');
            append_proportion(out, counts[static_cast<size_t>(i)], meta.n_haps);
        }
        out.push_back('\n');
        if (out.size() >= kFlushBytes) {
            if (std::fwrite(out.data(), 1, out.size(), fp) != out.size()) {
                die("failed writing %s", out_path.c_str());
            }
            out.clear();
        }
    }
    if (!out.empty() && std::fwrite(out.data(), 1, out.size(), fp) != out.size()) {
        die("failed writing %s", out_path.c_str());
    }
    if (std::fclose(fp) != 0) die("failed closing %s", out_path.c_str());

    std::fprintf(stderr, "Wrote %s: %llu blocks.\n", out_path.c_str(),
                 static_cast<unsigned long long>(reader.seen()));
}

}  // namespace

// argv: felixla_ancestry_export <mode> <in_prefix> <out_path> [threads]
int felixla_ancestry_export_main(int argc, char** argv) {
    if (argc < 4 || argc > 5) {
        std::fprintf(stderr,
            "Usage: %s <lai|global-admixture|local-admixture> <in_prefix> <out_path> [threads]\n",
            argc > 0 ? argv[0] : "felixla_ancestry_export");
        return 2;
    }

    std::string mode_name = argv[1];
    Mode mode;
    if (mode_name == "lai") mode = Mode::Lai;
    else if (mode_name == "global-admixture") mode = Mode::GlobalAdmixture;
    else if (mode_name == "local-admixture") mode = Mode::LocalAdmixture;
    else die("unknown ancestry export mode: %s", mode_name.c_str());

    std::string prefix = argv[2];
    std::string out_path = argv[3];
    int threads = 1;
    if (argc == 5) {
        threads = static_cast<int>(parse_u64(argv[4], "threads"));
        if (threads < 1) threads = 1;
    }

    Meta meta = read_meta(prefix + ".meta");
    std::vector<std::string> samples = read_samples(prefix + ".samples", meta.n_samples);

    switch (mode) {
        case Mode::Lai:
            export_lai(prefix, out_path, meta, samples, threads);
            break;
        case Mode::GlobalAdmixture:
            export_global_admixture(prefix, out_path, meta, samples);
            break;
        case Mode::LocalAdmixture:
            export_local_admixture(prefix, out_path, meta);
            break;
    }
    return 0;
}
