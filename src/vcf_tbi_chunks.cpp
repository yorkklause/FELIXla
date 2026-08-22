#include <htslib/hts.h>
#include <htslib/tbx.h>

#include <algorithm>
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
    uint64_t first_chunk_start = 0;
    uint64_t last_chunk_end = 0;
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

Read a phased VCF tabix index, find the first and last fixed-size chunk that may
contain indexed records, and generate non-overlapping 1-based inclusive regions.
The VCF body is not opened or decompressed. Bounds are conservative at tabix/CSI
bin granularity, so an edge chunk can be empty.

Required:
  --phase-vcf PATH          Phased VCF path used to discover its sidecar index.
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

class TabixIndexReader {
public:
    TabixIndexReader(const std::string& vcf_path, const std::string& index_path) {
        tbx_ = tbx_index_load2(vcf_path.c_str(), index_path.empty() ? nullptr : index_path.c_str());
        if (!tbx_) {
            if (index_path.empty()) {
                die("cannot load .tbi/.csi index beside phased VCF: " + vcf_path);
            }
            die("cannot load phased VCF index: " + index_path);
        }
        if ((tbx_->conf.preset & 0xffff) != TBX_VCF) {
            tbx_destroy(tbx_);
            tbx_ = nullptr;
            die("tabix index is not configured for VCF records");
        }
    }

    TabixIndexReader(const TabixIndexReader&) = delete;
    TabixIndexReader& operator=(const TabixIndexReader&) = delete;

    ~TabixIndexReader() {
        if (tbx_) tbx_destroy(tbx_);
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

    bool has_index_chunks(int tid, uint64_t begin, uint64_t end) const {
        if (begin >= end || end > kMaxVcfPosition) return false;
        hts_itr_t* itr = tbx_itr_queryi(
            tbx_, tid, static_cast<hts_pos_t>(begin), static_cast<hts_pos_t>(end)
        );
        if (!itr) die("cannot query phased VCF index");
        // Do not call tbx_itr_next: candidate offsets are already in the loaded index.
        bool found = !itr->finished && itr->n_off > 0;
        tbx_itr_destroy(itr);
        return found;
    }

    std::optional<std::pair<uint64_t, uint64_t>> chunk_bounds(
        int tid,
        uint64_t chunk_bp
    ) const {
        if (!has_index_chunks(tid, 0, kMaxVcfPosition)) return std::nullopt;

        const uint64_t chunk_count = (kMaxVcfPosition - 1) / chunk_bp + 1;
        uint64_t low = 0;
        uint64_t high = chunk_count - 1;
        while (low < high) {
            uint64_t middle = low + (high - low) / 2;
            uint64_t prefix_end = std::min((middle + 1) * chunk_bp, kMaxVcfPosition);
            if (has_index_chunks(tid, 0, prefix_end)) {
                high = middle;
            } else {
                low = middle + 1;
            }
        }
        const uint64_t first_chunk = low;

        low = first_chunk;
        high = chunk_count - 1;
        while (low < high) {
            uint64_t middle = low + (high - low + 1) / 2;
            uint64_t suffix_begin = middle * chunk_bp;
            if (has_index_chunks(tid, suffix_begin, kMaxVcfPosition)) {
                low = middle;
            } else {
                high = middle - 1;
            }
        }
        const uint64_t last_chunk = low;

        const uint64_t first = first_chunk * chunk_bp + 1;
        const uint64_t last = std::min((last_chunk + 1) * chunk_bp, kMaxVcfPosition);
        return std::make_pair(first, last);
    }

private:
    tbx_t* tbx_ = nullptr;
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
    TabixIndexReader& reader,
    const std::vector<std::string>& selected_chroms,
    uint64_t chunk_bp
) {
    std::unordered_set<std::string> selected(selected_chroms.begin(), selected_chroms.end());
    std::unordered_set<std::string> seen;
    std::vector<std::string> names = reader.sequence_names();
    std::vector<ContigExtent> extents;

    for (size_t tid = 0; tid < names.size(); ++tid) {
        const std::string& chrom = names[tid];
        if (!selected.empty() && selected.count(chrom) == 0) continue;
        seen.insert(chrom);

        std::optional<std::pair<uint64_t, uint64_t>> bounds =
            reader.chunk_bounds(static_cast<int>(tid), chunk_bp);
        if (!bounds) {
            std::cerr << "WARNING: index contig has no VCF records; skipping " << chrom << '\n';
            continue;
        }
        ContigExtent extent;
        extent.chrom = chrom;
        extent.first_chunk_start = bounds->first;
        extent.last_chunk_end = bounds->second;
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
        uint64_t aligned_start =
            ((extent.first_chunk_start - 1) / chunk_bp) * chunk_bp + 1;
        uint64_t chrom_chunk = 0;
        for (uint64_t start = aligned_start; start <= extent.last_chunk_end;) {
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
                extent.first_chunk_start,
                extent.last_chunk_end,
                extent.record_count,
            });
            if (end == std::numeric_limits<uint64_t>::max()) break;
            start = end + 1;
        }

        std::cerr << extent.chrom
                  << ": index_chunk_bounds=" << extent.first_chunk_start
                  << '-' << extent.last_chunk_end
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
    TabixIndexReader reader(options.phase_vcf, options.tbi_path);
    std::vector<ContigExtent> extents = read_extents(reader, options.chroms, options.chunk_bp);
    std::vector<Chunk> chunks = make_chunks(extents, options.chunk_bp);
    write_outputs(options, chunks);
    std::cerr << "Wrote " << chunks.size() << " chunk(s) to " << options.out_path;
    if (!options.commands_out_path.empty()) {
        std::cerr << " and commands to " << options.commands_out_path;
    }
    std::cerr << ".\n";
    return 0;
}
