#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr uint64_t kMaxVcfPosition =
    static_cast<uint64_t>(std::numeric_limits<int32_t>::max());

struct Options {
    std::string phase_vcf;
    std::string tbi_path;
    std::string out_path;
    std::string command_template;
    std::string commands_out_path;
    uint64_t chunk_bp = 0;
    int chunk_width = 4;
    std::vector<std::string> chroms;
};

struct ContigExtent {
    std::string chrom;
    uint64_t first_pos = 0;
    uint64_t last_pos = 0;
    std::optional<uint64_t> record_count;
};

struct Chunk {
    uint64_t global_chunk = 0;
    uint64_t chrom_chunk = 0;
    std::string chrom;
    uint64_t start = 0;
    uint64_t end = 0;
    uint64_t contig_first_pos = 0;
    uint64_t contig_last_pos = 0;
    std::optional<uint64_t> contig_records;
};

[[noreturn]] void die(const std::string& message) {
    std::cerr << "ERROR: " << message << '\n';
    std::exit(1);
}

void usage(std::ostream& out) {
    out <<
R"(Usage:
  vcf_tbi_chunks --phase-vcf FILE.vcf.gz --chunk-mb INT --out MANIFEST.tsv [options]
  vcf_tbi_chunks --phase-vcf FILE.vcf.gz --chunk-bp INT --out MANIFEST.tsv [options]

Read a phased VCF tabix index, find each contig's first and last VCF position,
and generate non-overlapping 1-based inclusive chunks aligned to the requested
chunk length.

Required:
  --phase-vcf PATH          BGZF-compressed phased VCF.
  --chunk-mb INT            Chunk length in decimal megabases (1 MB = 1,000,000 bp).
  --chunk-bp INT            Chunk length in bp; mutually exclusive with --chunk-mb.
  --out PATH                Output TSV manifest. Use - for stdout.

Index and filtering:
  --tbi PATH                Explicit .tbi/.csi path. Default: discover beside the VCF.
  --chrom CHROM             Include one contig. May be repeated.
  --chunk-width INT         Zero-padding width for command placeholders. Default: 4.

Optional command generation:
  --command-template STR    Emit one rendered command per chunk.
  --commands-out PATH       Command output file; required with --command-template.

Command placeholders:
  {phase_vcf} {phase_vcf_q} {chrom} {chrom_q} {start} {end}
  {region} {region_q} {global_chunk} {global_chunk0}
  {chrom_chunk} {chrom_chunk0} {contig_first_pos} {contig_last_pos}

The *_q placeholders are POSIX-shell quoted. The tool only writes the manifest
and commands; it does not invoke tabix, FELIXla, a shell, or a scheduler.
)";
}

uint64_t parse_positive_u64(const std::string& text, const std::string& flag) {
    if (text.empty()) die(flag + " requires a positive integer");
    uint64_t value = 0;
    for (char c : text) {
        if (c < '0' || c > '9') die(flag + " requires a positive integer: " + text);
        uint64_t digit = static_cast<uint64_t>(c - '0');
        if (value > (std::numeric_limits<uint64_t>::max() - digit) / 10) {
            die(flag + " is too large: " + text);
        }
        value = value * 10 + digit;
    }
    if (value == 0) die(flag + " must be greater than zero");
    return value;
}

int parse_chunk_width(const std::string& text) {
    uint64_t value = parse_positive_u64(text, "--chunk-width");
    if (value > 20) die("--chunk-width must be at most 20");
    return static_cast<int>(value);
}

std::string require_value(int& i, int argc, char** argv, const std::string& flag) {
    if (i + 1 >= argc) die(flag + " requires a value");
    return argv[++i];
}

Options parse_args(int argc, char** argv) {
    Options options;
    bool saw_chunk_bp = false;
    bool saw_chunk_mb = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage(std::cout);
            std::exit(0);
        }
        if (arg == "--version") {
            std::cout << "vcf_tbi_chunks v0\n";
            std::exit(0);
        }
        if (arg == "--phase-vcf" || arg == "--vcf") {
            options.phase_vcf = require_value(i, argc, argv, arg);
        } else if (arg == "--tbi" || arg == "--index") {
            options.tbi_path = require_value(i, argc, argv, arg);
        } else if (arg == "--chunk-bp") {
            if (saw_chunk_mb) die("use only one of --chunk-bp and --chunk-mb");
            options.chunk_bp = parse_positive_u64(require_value(i, argc, argv, arg), arg);
            saw_chunk_bp = true;
        } else if (arg == "--chunk-mb") {
            if (saw_chunk_bp) die("use only one of --chunk-bp and --chunk-mb");
            uint64_t chunk_mb = parse_positive_u64(require_value(i, argc, argv, arg), arg);
            if (chunk_mb > std::numeric_limits<uint64_t>::max() / 1000000ULL) {
                die("--chunk-mb is too large");
            }
            options.chunk_bp = chunk_mb * 1000000ULL;
            saw_chunk_mb = true;
        } else if (arg == "--out") {
            options.out_path = require_value(i, argc, argv, arg);
        } else if (arg == "--chrom" || arg == "--chr") {
            options.chroms.push_back(require_value(i, argc, argv, arg));
        } else if (arg == "--chunk-width") {
            options.chunk_width = parse_chunk_width(require_value(i, argc, argv, arg));
        } else if (arg == "--command-template") {
            options.command_template = require_value(i, argc, argv, arg);
        } else if (arg == "--commands-out") {
            options.commands_out_path = require_value(i, argc, argv, arg);
        } else {
            die("unknown argument: " + arg);
        }
    }

    if (options.phase_vcf.empty()) die("--phase-vcf is required");
    if (!saw_chunk_bp && !saw_chunk_mb) die("--chunk-bp or --chunk-mb is required");
    if (options.chunk_bp > kMaxVcfPosition) {
        die("chunk length exceeds the maximum standard VCF coordinate");
    }
    if (options.out_path.empty()) die("--out is required");
    if (options.command_template.empty() != options.commands_out_path.empty()) {
        die("--command-template and --commands-out must be supplied together");
    }
    if (options.out_path == "-" && options.commands_out_path == "-") {
        die("manifest and commands cannot both be written to stdout");
    }
    if (!options.commands_out_path.empty() && options.out_path == options.commands_out_path) {
        die("--out and --commands-out must be different paths");
    }
    if (options.command_template.find('\n') != std::string::npos ||
        options.command_template.find('\r') != std::string::npos) {
        die("--command-template must be a single line");
    }

    std::unordered_set<std::string> unique_chroms;
    for (const std::string& chrom : options.chroms) {
        if (chrom.empty()) die("--chrom cannot be empty");
        if (!unique_chroms.insert(chrom).second) die("duplicate --chrom: " + chrom);
    }
    return options;
}

uint64_t parse_vcf_position(const kstring_t& line) {
    const char* begin = line.s;
    const char* end = line.s + line.l;
    const char* first_tab = nullptr;
    const char* second_tab = nullptr;
    for (const char* p = begin; p < end; ++p) {
        if (*p != '\t') continue;
        if (!first_tab) {
            first_tab = p;
        } else {
            second_tab = p;
            break;
        }
    }
    if (!first_tab || !second_tab || first_tab + 1 == second_tab) {
        die("indexed VCF record does not contain CHROM and POS columns");
    }

    uint64_t pos = 0;
    for (const char* p = first_tab + 1; p < second_tab; ++p) {
        if (*p < '0' || *p > '9') die("indexed VCF record has a non-integer POS");
        uint64_t digit = static_cast<uint64_t>(*p - '0');
        if (pos > (kMaxVcfPosition - digit) / 10) {
            die("indexed VCF POS exceeds the standard VCF coordinate range");
        }
        pos = pos * 10 + digit;
    }
    if (pos == 0) die("indexed VCF POS must be positive");
    return pos;
}

class TabixVcfReader {
public:
    TabixVcfReader(const std::string& vcf_path, const std::string& index_path) {
        fp_ = hts_open(vcf_path.c_str(), "r");
        if (!fp_) die("cannot open phased VCF: " + vcf_path);
        if (!hts_get_bgzfp(fp_)) {
            hts_close(fp_);
            fp_ = nullptr;
            die("--phase-vcf must be a BGZF-compressed VCF: " + vcf_path);
        }

        tbx_ = tbx_index_load2(vcf_path.c_str(), index_path.empty() ? nullptr : index_path.c_str());
        if (!tbx_) {
            hts_close(fp_);
            fp_ = nullptr;
            if (index_path.empty()) {
                die("cannot load .tbi/.csi index beside phased VCF: " + vcf_path);
            }
            die("cannot load phased VCF index: " + index_path);
        }
        if ((tbx_->conf.preset & 0xffff) != TBX_VCF) {
            tbx_destroy(tbx_);
            tbx_ = nullptr;
            hts_close(fp_);
            fp_ = nullptr;
            die("tabix index is not configured for VCF records");
        }
    }

    TabixVcfReader(const TabixVcfReader&) = delete;
    TabixVcfReader& operator=(const TabixVcfReader&) = delete;

    ~TabixVcfReader() {
        std::free(line_.s);
        if (tbx_) tbx_destroy(tbx_);
        if (fp_) hts_close(fp_);
    }

    std::vector<std::string> sequence_names() const {
        int count = 0;
        const char** names = tbx_seqnames(tbx_, &count);
        if (count < 0 || (count > 0 && !names)) die("cannot read contig names from tabix index");
        std::vector<std::string> result;
        result.reserve(static_cast<size_t>(count));
        for (int i = 0; i < count; ++i) result.emplace_back(names[i]);
        std::free(names);
        return result;
    }

    std::optional<uint64_t> record_count(int tid) const {
        uint64_t mapped = 0;
        uint64_t unmapped = 0;
        if (hts_idx_get_stat(tbx_->idx, tid, &mapped, &unmapped) != 0) return std::nullopt;
        return mapped;
    }

    std::optional<uint64_t> first_position(int tid) {
        hts_itr_t* itr = tbx_itr_queryi(
            tbx_, tid, 0, static_cast<hts_pos_t>(kMaxVcfPosition)
        );
        if (!itr) return std::nullopt;
        int rc = tbx_itr_next(fp_, tbx_, itr, &line_);
        tbx_itr_destroy(itr);
        if (rc < 0) return std::nullopt;
        return parse_vcf_position(line_);
    }

    bool has_position_at_or_after(int tid, uint64_t target) {
        if (target == 0 || target > kMaxVcfPosition) return false;
        hts_pos_t begin = static_cast<hts_pos_t>(target - 1);
        hts_itr_t* itr = tbx_itr_queryi(
            tbx_, tid, begin, static_cast<hts_pos_t>(kMaxVcfPosition)
        );
        if (!itr) return false;

        bool found = false;
        while (tbx_itr_next(fp_, tbx_, itr, &line_) >= 0) {
            if (parse_vcf_position(line_) >= target) {
                found = true;
                break;
            }
        }
        tbx_itr_destroy(itr);
        return found;
    }

    uint64_t last_position(int tid, uint64_t first_pos) {
        if (first_pos == kMaxVcfPosition) return first_pos;

        uint64_t low = first_pos;
        uint64_t high = first_pos > kMaxVcfPosition / 2
            ? kMaxVcfPosition
            : first_pos * 2;
        if (high == low) ++high;

        while (has_position_at_or_after(tid, high)) {
            low = high;
            if (high == kMaxVcfPosition) return high;
            high = high > kMaxVcfPosition / 2
                ? kMaxVcfPosition
                : high * 2;
        }

        while (low + 1 < high) {
            uint64_t middle = low + (high - low) / 2;
            if (has_position_at_or_after(tid, middle)) {
                low = middle;
            } else {
                high = middle;
            }
        }
        return low;
    }

private:
    htsFile* fp_ = nullptr;
    tbx_t* tbx_ = nullptr;
    kstring_t line_{0, 0, nullptr};
};

std::string zero_pad(uint64_t value, int width) {
    std::ostringstream out;
    out << std::setw(width) << std::setfill('0') << value;
    return out.str();
}

void replace_all(std::string& value, const std::string& from, const std::string& to) {
    size_t pos = 0;
    while ((pos = value.find(from, pos)) != std::string::npos) {
        value.replace(pos, from.size(), to);
        pos += to.size();
    }
}

std::string shell_quote(const std::string& value) {
    std::string result = "'";
    for (char c : value) {
        if (c == '\'') {
            result += "'\"'\"'";
        } else {
            result.push_back(c);
        }
    }
    result.push_back('\'');
    return result;
}

std::string region_for(const Chunk& chunk) {
    return chunk.chrom + ":" + std::to_string(chunk.start) + "-" + std::to_string(chunk.end);
}

std::string render_command(
    std::string command,
    const Options& options,
    const Chunk& chunk
) {
    const std::string region = region_for(chunk);
    const std::vector<std::pair<std::string, std::string>> replacements = {
        {"{phase_vcf_q}", shell_quote(options.phase_vcf)},
        {"{phase_vcf}", options.phase_vcf},
        {"{chrom_q}", shell_quote(chunk.chrom)},
        {"{chrom}", chunk.chrom},
        {"{start}", std::to_string(chunk.start)},
        {"{end}", std::to_string(chunk.end)},
        {"{region_q}", shell_quote(region)},
        {"{region}", region},
        {"{global_chunk0}", zero_pad(chunk.global_chunk, options.chunk_width)},
        {"{global_chunk}", std::to_string(chunk.global_chunk)},
        {"{chrom_chunk0}", zero_pad(chunk.chrom_chunk, options.chunk_width)},
        {"{chrom_chunk}", std::to_string(chunk.chrom_chunk)},
        {"{contig_first_pos}", std::to_string(chunk.contig_first_pos)},
        {"{contig_last_pos}", std::to_string(chunk.contig_last_pos)},
    };
    for (const auto& replacement : replacements) {
        replace_all(command, replacement.first, replacement.second);
    }
    return command;
}

std::vector<ContigExtent> read_extents(
    TabixVcfReader& reader,
    const std::vector<std::string>& selected_chroms
) {
    std::unordered_set<std::string> selected(selected_chroms.begin(), selected_chroms.end());
    std::unordered_set<std::string> seen;
    std::vector<std::string> names = reader.sequence_names();
    std::vector<ContigExtent> extents;

    for (size_t tid = 0; tid < names.size(); ++tid) {
        const std::string& chrom = names[tid];
        if (!selected.empty() && selected.count(chrom) == 0) continue;
        seen.insert(chrom);

        std::optional<uint64_t> first = reader.first_position(static_cast<int>(tid));
        if (!first) {
            std::cerr << "WARNING: index contig has no VCF records; skipping " << chrom << '\n';
            continue;
        }
        ContigExtent extent;
        extent.chrom = chrom;
        extent.first_pos = *first;
        extent.last_pos = reader.last_position(static_cast<int>(tid), *first);
        extent.record_count = reader.record_count(static_cast<int>(tid));
        extents.push_back(std::move(extent));
    }

    for (const std::string& chrom : selected_chroms) {
        if (seen.count(chrom) == 0) die("--chrom is absent from the tabix index: " + chrom);
    }
    if (extents.empty()) die("no indexed VCF records matched the requested contigs");
    return extents;
}

std::vector<Chunk> make_chunks(
    const std::vector<ContigExtent>& extents,
    uint64_t chunk_bp
) {
    std::vector<Chunk> chunks;
    uint64_t global_chunk = 0;

    for (const ContigExtent& extent : extents) {
        uint64_t aligned_start = ((extent.first_pos - 1) / chunk_bp) * chunk_bp + 1;
        uint64_t chrom_chunk = 0;
        for (uint64_t start = aligned_start; start <= extent.last_pos;) {
            if (start > std::numeric_limits<uint64_t>::max() - (chunk_bp - 1)) {
                die("chunk coordinate overflow on " + extent.chrom);
            }
            uint64_t end = start + chunk_bp - 1;
            ++global_chunk;
            ++chrom_chunk;
            chunks.push_back(Chunk{
                global_chunk,
                chrom_chunk,
                extent.chrom,
                start,
                end,
                extent.first_pos,
                extent.last_pos,
                extent.record_count,
            });
            if (end == std::numeric_limits<uint64_t>::max()) break;
            start = end + 1;
        }

        const Chunk& first_chunk = chunks[chunks.size() - static_cast<size_t>(chrom_chunk)];
        const Chunk& last_chunk = chunks.back();
        std::cerr << extent.chrom
                  << ": first=" << extent.first_pos
                  << ", last=" << extent.last_pos
                  << ", aligned=" << first_chunk.start << '-' << last_chunk.end
                  << ", chunks=" << chrom_chunk;
        if (extent.record_count) std::cerr << ", records=" << *extent.record_count;
        std::cerr << '\n';
    }
    return chunks;
}

void write_outputs(
    const Options& options,
    const std::vector<Chunk>& chunks
) {
    std::ofstream manifest_file;
    std::ostream* manifest = &std::cout;
    if (options.out_path != "-") {
        manifest_file.open(options.out_path);
        if (!manifest_file) die("cannot open manifest output: " + options.out_path);
        manifest = &manifest_file;
    }

    std::ofstream commands_file;
    std::ostream* commands = nullptr;
    if (!options.command_template.empty()) {
        if (options.commands_out_path == "-") {
            commands = &std::cout;
        } else {
            commands_file.open(options.commands_out_path);
            if (!commands_file) die("cannot open command output: " + options.commands_out_path);
            commands = &commands_file;
        }
    }

    *manifest
        << "global_chunk\tchrom_chunk\tchrom\tstart\tend\tregion\t"
        << "contig_first_pos\tcontig_last_pos\tcontig_records\n";
    for (const Chunk& chunk : chunks) {
        *manifest
            << chunk.global_chunk << '\t'
            << chunk.chrom_chunk << '\t'
            << chunk.chrom << '\t'
            << chunk.start << '\t'
            << chunk.end << '\t'
            << region_for(chunk) << '\t'
            << chunk.contig_first_pos << '\t'
            << chunk.contig_last_pos << '\t';
        if (chunk.contig_records) {
            *manifest << *chunk.contig_records;
        } else {
            *manifest << '.';
        }
        *manifest << '\n';

        if (commands) *commands << render_command(options.command_template, options, chunk) << '\n';
    }

    if (!*manifest) die("failed while writing manifest: " + options.out_path);
    if (commands && !*commands) die("failed while writing commands: " + options.commands_out_path);
}

}  // namespace

int main(int argc, char** argv) {
    Options options = parse_args(argc, argv);
    TabixVcfReader reader(options.phase_vcf, options.tbi_path);
    std::vector<ContigExtent> extents = read_extents(reader, options.chroms);
    std::vector<Chunk> chunks = make_chunks(extents, options.chunk_bp);
    write_outputs(options, chunks);
    std::cerr << "Wrote " << chunks.size() << " chunk(s) to " << options.out_path;
    if (!options.commands_out_path.empty()) {
        std::cerr << " and commands to " << options.commands_out_path;
    }
    std::cerr << ".\n";
    return 0;
}
