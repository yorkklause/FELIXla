// flare_subset_to_tractor_hybrid.cpp
// designed by Kai, implemented by codex
//
// Convert phased genotype VCF/BCF plus FLARE local ancestry VCF/BCF into
// ancestry-aware packed files for a SAIGE-TRACTOR genotype backend.

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <unistd.h>

struct RareCarrierPacked {
    uint32_t pos_index; // split biallelic variant ordinal
    uint32_t anc_hap;   // high 5 bits ancestry, low 27 bits hap_id
};

struct AltCarrierBuilder {
    uint32_t mac = 0;
    std::vector<uint32_t> sparse_haps;
    std::vector<uint64_t> dense_bits;
};

struct AncestryState {
    std::vector<std::vector<uint64_t>> masks;
    std::vector<int8_t> hap_ancestry;
};

struct OpenAncestryBlock {
    bool active = false;
    uint32_t block_id = 0;
    int geno_rid = -1;
    std::string chr;
    int64_t start_pos = -1;
    int64_t end_pos = -1;
    AncestryState state;
};

struct AncestryChange {
    uint32_t hap_id = 0;
    int8_t ancestry = -1;
};

struct LaiRecord {
    bool valid = false;
    bool reset_state = false;
    bool merged_duplicate = false;
    int geno_rid = -1;
    std::string chr;
    int64_t pos = 0;
    std::vector<AncestryChange> changes;
};

struct FlareDeltaDecoder {
    int geno_rid = -1;
    std::vector<int8_t> hap_ancestry;
};

struct Region {
    bool active = false;
    std::string label;
    std::string chr;
    int64_t start = 0;
    int64_t end = 0;
    int geno_rid = -1;
};

struct KeepSamples {
    bool active = false;
    std::string path;
    std::unordered_set<std::string> ids;
};

struct SampleSelection {
    std::vector<int> genotype_raw_indices;
    std::vector<int> flare_raw_indices;
    std::vector<int> genotype_text_indices;
    std::vector<int> flare_text_indices;
    std::vector<std::string> sample_ids;
    int genotype_raw_sample_count = 0;
    int flare_raw_sample_count = 0;
    bool genotype_subset_active = false;
    bool flare_subset_active = false;
    std::string genotype_decode_samples;
    std::string flare_decode_samples;
    bool keep_active = false;
    std::string keep_path;
};

struct FastTextVcfReader {
    htsFile* fp = nullptr;
    bool owns_fp = false;
    tbx_t* tbx = nullptr;
    hts_itr_t* itr = nullptr;
    std::vector<std::string> queries;
    size_t query_index = 0;
    kstring_t line = {0, 0, nullptr};
};

struct FastGenotypeRecord {
    int rid = -1;
    int64_t pos = 0;
    char* chr = nullptr;
    char* id = nullptr;
    char* ref = nullptr;
    std::vector<char*> alleles;
    char* samples = nullptr;
    size_t samples_length = 0;
    int gt_field_index = -1;
    bool gt_only = false;
};

struct ExtractPosition {
    std::string chr;
    int64_t pos = 0;
    std::string ref;
    std::vector<std::string> alts;
    int geno_rid = -1;
    uint64_t line_no = 0;
};

struct ExtractSites {
    bool active = false;
    std::string path;
    std::unordered_set<std::string> contigs;
    std::vector<ExtractPosition> positions;
    uint64_t allele_count = 0;
};

static constexpr uint32_t kMaxPackedHapId = (1u << 27) - 1u;
static constexpr uint64_t kProgressRecordInterval = 10000;
static constexpr int kProgressSecondsInterval = 2;
static constexpr int64_t kMaxVcfCoordinate = 2147483647LL;
static constexpr int kInputBlockSize = 8 * 1024 * 1024;

[[noreturn]] static void die(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    std::fputs("ERROR: ", stderr);
    std::vfprintf(stderr, fmt, args);
    std::fputc('\n', stderr);
    va_end(args);
    std::exit(1);
}

struct IndexRecordCount {
    bool available = false;
    uint64_t total = 0;
};

static hts_idx_t* load_any_vcf_index(const char* path) {
    hts_idx_t* idx = bcf_index_load3(path, nullptr, HTS_IDX_SILENT_FAIL);
    if (idx) return idx;

    idx = hts_idx_load3(path, nullptr, HTS_FMT_CSI, HTS_IDX_SILENT_FAIL);
    if (idx) return idx;
    return hts_idx_load3(path, nullptr, HTS_FMT_TBI, HTS_IDX_SILENT_FAIL);
}

static int64_t parse_i64_string(const std::string& text, const char* label) {
    if (text.empty()) die("%s is empty", label);

    char* end = nullptr;
    long long value = std::strtoll(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0') {
        die("invalid integer for %s: %s", label, text.c_str());
    }
    if (value <= 0) {
        die("%s must be positive: %s", label, text.c_str());
    }
    return static_cast<int64_t>(value);
}

static Region parse_region_string(const std::string& region_text) {
    size_t colon = region_text.find(':');
    size_t dash = region_text.find('-', colon == std::string::npos ? 0 : colon + 1);
    if (colon == std::string::npos || dash == std::string::npos || dash <= colon + 1) {
        die("region must look like chr:start-end: %s", region_text.c_str());
    }

    Region region;
    region.active = true;
    region.chr = region_text.substr(0, colon);
    region.start = parse_i64_string(region_text.substr(colon + 1, dash - colon - 1), "region start");
    region.end = parse_i64_string(region_text.substr(dash + 1), "region end");
    if (region.chr.empty()) die("region chromosome is empty");
    if (region.start > region.end) {
        die("invalid region coordinates: %s", region_text.c_str());
    }
    region.label = region.chr + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
    return region;
}

class TextLineReader {
public:
    TextLineReader(const std::string& path, const char* label)
        : path_(path) {
        fp_ = hts_open(path.c_str(), "r");
        if (!fp_) die("cannot open %s: %s", label, path.c_str());
    }

    TextLineReader(const TextLineReader&) = delete;
    TextLineReader& operator=(const TextLineReader&) = delete;

    ~TextLineReader() {
        std::free(line_.s);
        if (fp_) hts_close(fp_);
    }

    bool next(std::string& out) {
        int ret = hts_getline(fp_, '\n', &line_);
        if (ret >= 0) {
            out.assign(line_.s, static_cast<size_t>(line_.l));
            return true;
        }
        if (ret == -1) return false;
        die("error reading text file: %s", path_.c_str());
    }

private:
    std::string path_;
    htsFile* fp_ = nullptr;
    kstring_t line_{0, 0, nullptr};
};

static std::string strip_trailing_cr(std::string value) {
    if (!value.empty() && value.back() == '\r') {
        value.pop_back();
    }
    return value;
}

static std::vector<std::string> split_commas(const std::string& value) {
    std::vector<std::string> parts;
    size_t begin = 0;
    while (begin <= value.size()) {
        size_t comma = value.find(',', begin);
        size_t end = comma == std::string::npos ? value.size() : comma;
        std::string part = value.substr(begin, end - begin);
        if (part.empty()) die("empty ALT allele in --extract list");
        parts.push_back(part);
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    return parts;
}

static KeepSamples load_keep_samples(const std::string& path) {
    KeepSamples keep;
    if (path.empty()) return keep;

    keep.active = true;
    keep.path = path;

    TextLineReader reader(path, "--keep sample list");

    std::string line;
    uint64_t line_no = 0;
    while (reader.next(line)) {
        ++line_no;
        line = strip_trailing_cr(line);
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string sample_id;
        std::string extra;
        iss >> sample_id;
        if (sample_id.empty()) continue;
        if (iss >> extra) {
            die("--keep expects one sample ID per line at %s:%llu", path.c_str(),
                static_cast<unsigned long long>(line_no));
        }
        if (!keep.ids.insert(sample_id).second) {
            die("duplicate sample ID in --keep list at %s:%llu: %s", path.c_str(),
                static_cast<unsigned long long>(line_no), sample_id.c_str());
        }
    }
    if (keep.ids.empty()) {
        die("--keep sample list is empty: %s", path.c_str());
    }

    return keep;
}

static bool sample_headers_identical(bcf_hdr_t* ghdr, bcf_hdr_t* ahdr) {
    int n_genotype = bcf_hdr_nsamples(ghdr);
    int n_flare = bcf_hdr_nsamples(ahdr);
    if (n_genotype != n_flare) return false;

    for (int i = 0; i < n_genotype; ++i) {
        if (std::strcmp(ghdr->samples[i], ahdr->samples[i]) != 0) {
            return false;
        }
    }
    return true;
}

struct HeaderSample {
    std::string id;
    int raw_index = -1;
};

static std::vector<HeaderSample> sorted_header_samples(
    bcf_hdr_t* hdr,
    const char* label
) {
    int n_samples = bcf_hdr_nsamples(hdr);
    std::vector<HeaderSample> samples;
    samples.reserve(static_cast<size_t>(n_samples));

    for (int i = 0; i < n_samples; ++i) {
        samples.push_back(HeaderSample{hdr->samples[i], i});
    }
    std::sort(
        samples.begin(),
        samples.end(),
        [](const HeaderSample& a, const HeaderSample& b) {
            return a.id < b.id;
        }
    );

    for (size_t i = 1; i < samples.size(); ++i) {
        if (samples[i - 1].id == samples[i].id) {
            die("duplicate sample ID in %s VCF header: %s", label, samples[i].id.c_str());
        }
    }
    return samples;
}

static SampleSelection build_sample_selection(
    bcf_hdr_t* ghdr,
    bcf_hdr_t* ahdr,
    const std::string& keep_path
) {
    SampleSelection selection;
    selection.genotype_raw_sample_count = bcf_hdr_nsamples(ghdr);
    selection.flare_raw_sample_count = bcf_hdr_nsamples(ahdr);
    if (selection.genotype_raw_sample_count <= 0) die("genotype VCF has no samples");
    if (selection.flare_raw_sample_count <= 0) die("FLARE VCF has no samples");

    KeepSamples keep = load_keep_samples(keep_path);
    selection.keep_active = keep.active;
    selection.keep_path = keep.path;
    std::unordered_set<std::string> keep_missing_from_genotype = keep.ids;

    bool identical = sample_headers_identical(ghdr, ahdr);

    selection.genotype_raw_indices.reserve(static_cast<size_t>(selection.genotype_raw_sample_count));
    selection.flare_raw_indices.reserve(static_cast<size_t>(selection.genotype_raw_sample_count));
    selection.sample_ids.reserve(static_cast<size_t>(selection.genotype_raw_sample_count));

    if (identical) {
        (void)sorted_header_samples(ghdr, "genotype");
        for (int i = 0; i < selection.genotype_raw_sample_count; ++i) {
            std::string sample(ghdr->samples[i]);
            if (keep.active) {
                if (keep.ids.count(sample) == 0) continue;
                keep_missing_from_genotype.erase(sample);
            }
            selection.genotype_raw_indices.push_back(i);
            selection.flare_raw_indices.push_back(i);
            selection.sample_ids.push_back(std::move(sample));
        }
    } else {
        std::vector<HeaderSample> genotype_samples =
            sorted_header_samples(ghdr, "genotype");
        std::vector<HeaderSample> flare_samples =
            sorted_header_samples(ahdr, "FLARE");
        if (keep.active) {
            for (const HeaderSample& sample : genotype_samples) {
                keep_missing_from_genotype.erase(sample.id);
            }
        }

        struct MatchedSample {
            std::string id;
            int genotype_raw_index = -1;
            int flare_raw_index = -1;
        };
        std::vector<MatchedSample> matches;
        matches.reserve(std::min(genotype_samples.size(), flare_samples.size()));

        size_t genotype_i = 0;
        size_t flare_i = 0;
        while (genotype_i < genotype_samples.size() && flare_i < flare_samples.size()) {
            const HeaderSample& genotype_sample = genotype_samples[genotype_i];
            const HeaderSample& flare_sample = flare_samples[flare_i];
            if (genotype_sample.id < flare_sample.id) {
                ++genotype_i;
            } else if (flare_sample.id < genotype_sample.id) {
                ++flare_i;
            } else {
                if (!keep.active || keep.ids.count(genotype_sample.id) != 0) {
                    matches.push_back(MatchedSample{
                        genotype_sample.id,
                        genotype_sample.raw_index,
                        flare_sample.raw_index
                    });
                }
                ++genotype_i;
                ++flare_i;
            }
        }

        std::sort(
            matches.begin(),
            matches.end(),
            [](const MatchedSample& a, const MatchedSample& b) {
                return a.genotype_raw_index < b.genotype_raw_index;
            }
        );
        for (MatchedSample& match : matches) {
            selection.genotype_raw_indices.push_back(match.genotype_raw_index);
            selection.flare_raw_indices.push_back(match.flare_raw_index);
            selection.sample_ids.push_back(std::move(match.id));
        }
    }

    if (keep.active && !keep_missing_from_genotype.empty()) {
        die("sample ID in --keep list is absent from genotype VCF: %s",
            keep_missing_from_genotype.begin()->c_str());
    }

    if (selection.sample_ids.empty()) {
        if (keep.active) {
            die("--keep retained zero samples after genotype/FLARE intersection: %s",
                keep.path.c_str());
        }
        die("genotype and FLARE sample intersection is empty");
    }

    if (keep.active) {
        std::fprintf(
            stderr,
            "Loaded --keep sample list: retaining %llu of %llu requested sample(s) after genotype/FLARE intersection.\n",
            static_cast<unsigned long long>(selection.sample_ids.size()),
            static_cast<unsigned long long>(keep.ids.size())
        );
    }

    if (!identical || static_cast<int>(selection.sample_ids.size()) != selection.genotype_raw_sample_count ||
        static_cast<int>(selection.sample_ids.size()) != selection.flare_raw_sample_count) {
        std::fprintf(
            stderr,
            "Using genotype/FLARE sample intersection: retaining %llu sample(s) in genotype VCF order (genotype %d, FLARE %d).\n",
            static_cast<unsigned long long>(selection.sample_ids.size()),
            selection.genotype_raw_sample_count,
            selection.flare_raw_sample_count
        );
    }

    return selection;
}

static std::string comma_join_samples(const std::vector<std::string>& samples) {
    std::string joined;
    size_t total_bytes = samples.empty() ? 0 : samples.size() - 1;
    for (const std::string& sample : samples) total_bytes += sample.size();
    joined.reserve(total_bytes);

    for (size_t i = 0; i < samples.size(); ++i) {
        if (i != 0) joined.push_back(',');
        joined += samples[i];
    }
    return joined;
}

static void set_header_samples_or_die(
    bcf_hdr_t* hdr,
    const std::string& samples,
    const char* label
) {
    int ret = bcf_hdr_set_samples(hdr, samples.c_str(), 0);
    if (ret < 0) {
        die("failed to configure %s sample subset", label);
    }
    if (ret > 0) {
        die("failed to configure %s sample subset: selected sample is absent", label);
    }
}

static void apply_decode_sample_subsets(
    bcf_hdr_t* ghdr,
    bcf_hdr_t* ahdr,
    SampleSelection& selection
) {
    int n_selected = static_cast<int>(selection.sample_ids.size());

    if (n_selected < selection.genotype_raw_sample_count) {
        selection.genotype_subset_active = true;
        selection.genotype_decode_samples = comma_join_samples(selection.sample_ids);
        set_header_samples_or_die(
            ghdr,
            selection.genotype_decode_samples,
            "genotype VCF"
        );
        selection.genotype_raw_indices.resize(static_cast<size_t>(n_selected));
        for (int i = 0; i < n_selected; ++i) {
            selection.genotype_raw_indices[static_cast<size_t>(i)] = i;
        }
        selection.genotype_raw_sample_count = n_selected;
    }

    if (n_selected < selection.flare_raw_sample_count) {
        selection.flare_subset_active = true;

        std::vector<std::pair<int, int>> raw_and_output;
        raw_and_output.reserve(static_cast<size_t>(n_selected));
        for (int output_i = 0; output_i < n_selected; ++output_i) {
            raw_and_output.emplace_back(
                selection.flare_raw_indices[static_cast<size_t>(output_i)],
                output_i
            );
        }
        std::sort(raw_and_output.begin(), raw_and_output.end());

        std::vector<std::string> flare_source_order_ids;
        flare_source_order_ids.reserve(static_cast<size_t>(n_selected));
        std::vector<int> compact_flare_indices(static_cast<size_t>(n_selected), -1);
        for (int compact_i = 0; compact_i < n_selected; ++compact_i) {
            int output_i = raw_and_output[static_cast<size_t>(compact_i)].second;
            flare_source_order_ids.push_back(
                selection.sample_ids[static_cast<size_t>(output_i)]
            );
            compact_flare_indices[static_cast<size_t>(output_i)] = compact_i;
        }

        selection.flare_decode_samples = comma_join_samples(flare_source_order_ids);
        set_header_samples_or_die(
            ahdr,
            selection.flare_decode_samples,
            "FLARE VCF"
        );
        selection.flare_raw_indices = std::move(compact_flare_indices);
        selection.flare_raw_sample_count = n_selected;
    }
}

static ExtractSites load_extract_sites(const std::string& path) {
    ExtractSites sites;
    if (path.empty()) return sites;

    sites.active = true;
    sites.path = path;

    TextLineReader reader(path, "--extract site list");

    std::string line;
    uint64_t line_no = 0;
    while (reader.next(line)) {
        ++line_no;
        line = strip_trailing_cr(line);
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string chr;
        std::string pos_text;
        std::string ignored_id;
        std::string ref;
        std::string alts;
        // Consume the ID column only to preserve PVAR/VCF column alignment.
        // Site filtering is keyed exclusively by CHROM, POS, REF, and ALT.
        if (!(iss >> chr >> pos_text >> ignored_id >> ref >> alts)) {
            die("--extract expects PVAR/VCF columns CHROM POS ID REF ALT at %s:%llu",
                path.c_str(), static_cast<unsigned long long>(line_no));
        }
        if (chr == "CHROM" || chr == "#CHROM") continue;

        int64_t pos = parse_i64_string(pos_text, "--extract POS");
        if (ref.empty() || ref == ".") {
            die("--extract requires known REF at %s:%llu", path.c_str(),
                static_cast<unsigned long long>(line_no));
        }
        if (alts.empty() || alts == ".") {
            die("--extract requires known ALT at %s:%llu", path.c_str(),
                static_cast<unsigned long long>(line_no));
        }

        ExtractPosition position;
        position.chr = chr;
        position.pos = pos;
        position.ref = ref;
        position.alts = split_commas(alts);
        position.line_no = line_no;
        std::sort(position.alts.begin(), position.alts.end());
        for (size_t i = 1; i < position.alts.size(); ++i) {
            if (position.alts[i - 1] == position.alts[i]) {
                die(
                    "duplicate allele in --extract list at %s:%llu: %s:%lld %s>%s",
                    path.c_str(),
                    static_cast<unsigned long long>(line_no),
                    chr.c_str(),
                    static_cast<long long>(pos),
                    ref.c_str(),
                    position.alts[i].c_str()
                );
            }
        }
        sites.allele_count += position.alts.size();
        sites.contigs.insert(position.chr);
        sites.positions.push_back(std::move(position));
    }

    if (sites.positions.empty()) {
        die("--extract site list is empty: %s", path.c_str());
    }

    return sites;
}

static IndexRecordCount get_index_record_count(const char* path) {
    IndexRecordCount result;
    hts_idx_t* idx = load_any_vcf_index(path);
    if (!idx) return result;

    int nseq = hts_idx_nseq(idx);
    if (nseq <= 0) {
        hts_idx_destroy(idx);
        return result;
    }

    uint64_t total = 0;
    for (int tid = 0; tid < nseq; ++tid) {
        uint64_t mapped = 0;
        uint64_t unmapped = 0;
        if (hts_idx_get_stat(idx, tid, &mapped, &unmapped) != 0) {
            hts_idx_destroy(idx);
            return result;
        }

        if (mapped > std::numeric_limits<uint64_t>::max() - total) {
            total = std::numeric_limits<uint64_t>::max();
        } else {
            total += mapped;
        }
    }

    hts_idx_destroy(idx);
    result.available = true;
    result.total = total;
    return result;
}

class ProgressReporter {
public:
    ProgressReporter(
        const char* input_label,
        const char* input_path,
        const char* converted_label,
        const char* region_label = nullptr
    )
        : input_label_(input_label),
          converted_label_(converted_label),
          total_(region_label ? IndexRecordCount{} : get_index_record_count(input_path)),
          last_report_(std::chrono::steady_clock::now()),
          stderr_is_tty_(isatty(fileno(stderr)) != 0) {
        if (region_label) {
            std::fprintf(
                stderr,
                "Progress: %s region %s; reporting scanned region records only.\n",
                input_label_,
                region_label
            );
        } else if (total_.available) {
            std::fprintf(
                stderr,
                "Progress: %s index reports %llu records.\n",
                input_label_,
                static_cast<unsigned long long>(total_.total)
            );
        } else {
            std::fprintf(
                stderr,
                "Progress: %s index record count unavailable; reporting scanned records only.\n",
                input_label_
            );
        }
    }

    void record_scanned() {
        ++records_scanned_;
    }

    void maybe_report(uint64_t converted, uint64_t common, uint64_t rare, bool force = false) {
        auto now = std::chrono::steady_clock::now();
        bool count_due = records_scanned_ >= next_record_report_;
        bool time_due = now - last_report_ >= std::chrono::seconds(kProgressSecondsInterval);

        if (!force && !count_due && !time_due) {
            return;
        }

        print(converted, common, rare);
        last_report_ = now;
        printed_ = true;

        if (records_scanned_ >= next_record_report_) {
            if (records_scanned_ > std::numeric_limits<uint64_t>::max() - kProgressRecordInterval) {
                next_record_report_ = std::numeric_limits<uint64_t>::max();
            } else {
                next_record_report_ = records_scanned_ + kProgressRecordInterval;
            }
        }
    }

    void finish(uint64_t converted, uint64_t common, uint64_t rare) {
        maybe_report(converted, common, rare, true);
        if (stderr_is_tty_ && printed_) {
            std::fputc('\n', stderr);
        }
    }

private:
    void print(uint64_t converted, uint64_t common, uint64_t rare) const {
        const char* prefix = stderr_is_tty_ ? "\r" : "";
        const char* suffix = stderr_is_tty_ ? "" : "\n";

        if (total_.available && total_.total > 0) {
            uint64_t capped_records = std::min(records_scanned_, total_.total);
            double pct = 100.0 * static_cast<double>(capped_records) /
                         static_cast<double>(total_.total);
            std::fprintf(
                stderr,
                "%sProgress: records %llu/%llu (%.1f%%), converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(total_.total),
                pct,
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        } else if (total_.available) {
            std::fprintf(
                stderr,
                "%sProgress: records %llu/0, converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        } else {
            std::fprintf(
                stderr,
                "%sProgress: records %llu, converted %llu %s, common %llu, rare %llu%s",
                prefix,
                static_cast<unsigned long long>(records_scanned_),
                static_cast<unsigned long long>(converted),
                converted_label_,
                static_cast<unsigned long long>(common),
                static_cast<unsigned long long>(rare),
                suffix
            );
        }

        std::fflush(stderr);
    }

    const char* input_label_;
    const char* converted_label_;
    IndexRecordCount total_;
    uint64_t records_scanned_ = 0;
    uint64_t next_record_report_ = kProgressRecordInterval;
    std::chrono::steady_clock::time_point last_report_;
    bool stderr_is_tty_ = false;
    bool printed_ = false;
};

static inline uint32_t pack_anc_hap(uint32_t ancestry, uint32_t hap_id) {
    return ((ancestry & 31u) << 27) | (hap_id & 0x07ffffffu);
}

static inline void set_bit(std::vector<uint64_t>& bits, uint32_t h) {
    bits[h >> 6] |= (1ULL << (h & 63));
}

static inline void clear_bit(std::vector<uint64_t>& bits, uint32_t h) {
    bits[h >> 6] &= ~(1ULL << (h & 63));
}

static int integer_format_width(const bcf_fmt_t* fmt, const char* label) {
    if (fmt->type == BCF_BT_INT8) return 1;
    if (fmt->type == BCF_BT_INT16) return 2;
    if (fmt->type == BCF_BT_INT32) return 4;
    if (fmt->type == BCF_BT_INT64) return 8;
    die("%s FORMAT field is not integer encoded", label);
}

static int64_t integer_format_value(
    const bcf_fmt_t* fmt,
    int sample_index,
    int value_index,
    const char* label
) {
    if (sample_index < 0 || value_index < 0 || value_index >= fmt->n) {
        die("internal %s FORMAT index out of range", label);
    }
    int width = integer_format_width(fmt, label);
    const uint8_t* value_ptr = fmt->p +
        static_cast<size_t>(sample_index) * static_cast<size_t>(fmt->size) +
        static_cast<size_t>(value_index) * static_cast<size_t>(width);
    uint8_t* next = nullptr;
    int64_t value = bcf_dec_int1(value_ptr, fmt->type, &next);
    if (!next) die("failed to decode %s FORMAT integer", label);
    return value;
}

static bool integer_format_missing(const bcf_fmt_t* fmt, int64_t value) {
    if (fmt->type == BCF_BT_INT8) return value == bcf_int8_missing;
    if (fmt->type == BCF_BT_INT16) return value == bcf_int16_missing;
    if (fmt->type == BCF_BT_INT32) return value == bcf_int32_missing;
    if (fmt->type == BCF_BT_INT64) return value == bcf_int64_missing;
    return false;
}

static bool integer_format_vector_end(const bcf_fmt_t* fmt, int64_t value) {
    if (fmt->type == BCF_BT_INT8) return value == bcf_int8_vector_end;
    if (fmt->type == BCF_BT_INT16) return value == bcf_int16_vector_end;
    if (fmt->type == BCF_BT_INT32) return value == bcf_int32_vector_end;
    if (fmt->type == BCF_BT_INT64) return value == bcf_int64_vector_end;
    return false;
}

static bcf_fmt_t* require_integer_format(
    bcf_hdr_t* hdr,
    bcf1_t* rec,
    const char* tag
) {
    bcf_unpack(rec, BCF_UN_FMT);
    bcf_fmt_t* fmt = bcf_get_fmt(hdr, rec, tag);
    if (!fmt) die("missing FORMAT/%s", tag);
    int width = integer_format_width(fmt, tag);
    if (fmt->n <= 0 || fmt->size < fmt->n * width || !fmt->p) {
        die("invalid FORMAT/%s storage", tag);
    }
    return fmt;
}

static uint64_t tell_or_die(FILE* fp, const char* path) {
    off_t offset = ftello(fp);
    if (offset < 0) die("ftello failed for %s", path);
    return static_cast<uint64_t>(offset);
}

static FILE* open_output_or_die(const std::string& path, const char* mode) {
    FILE* fp = std::fopen(path.c_str(), mode);
    if (!fp) die("cannot open output file: %s", path.c_str());

    if (setvbuf(fp, nullptr, _IOFBF, 1 << 24) != 0) {
        die("setvbuf failed for output file: %s", path.c_str());
    }

    return fp;
}

static void write_exact(FILE* fp, const void* data, size_t bytes, const char* what) {
    if (bytes == 0) return;
    if (std::fwrite(data, 1, bytes, fp) != bytes) {
        die("failed writing %s", what);
    }
}

static void write_u32_le(FILE* fp, uint32_t value, const char* what) {
    uint8_t bytes[4] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24)
    };
    write_exact(fp, bytes, sizeof(bytes), what);
}

static void write_u64_le(FILE* fp, uint64_t value, const char* what) {
    uint8_t bytes[8] = {
        static_cast<uint8_t>(value),
        static_cast<uint8_t>(value >> 8),
        static_cast<uint8_t>(value >> 16),
        static_cast<uint8_t>(value >> 24),
        static_cast<uint8_t>(value >> 32),
        static_cast<uint8_t>(value >> 40),
        static_cast<uint8_t>(value >> 48),
        static_cast<uint8_t>(value >> 56)
    };
    write_exact(fp, bytes, sizeof(bytes), what);
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

static void write_magic_header(FILE* fp, const char magic[8]) {
    write_exact(fp, magic, 8, "index header");
}

static void write_common_offset_idx_record(
    FILE* fp,
    uint64_t common_index,
    uint32_t global_variant_index,
    uint64_t mks_offset,
    uint64_t geno_offset
) {
    write_u64_le(fp, common_index, "common offset index");
    write_u32_le(fp, global_variant_index, "common offset index");
    write_u64_le(fp, mks_offset, "common offset index");
    write_u64_le(fp, geno_offset, "common offset index");
}

static void write_rare_offset_idx_record(
    FILE* fp,
    uint64_t rare_index,
    uint32_t global_variant_index,
    uint64_t mks_offset,
    uint64_t carrier_offset,
    uint32_t n_carriers
) {
    write_u64_le(fp, rare_index, "rare offset index");
    write_u32_le(fp, global_variant_index, "rare offset index");
    write_u64_le(fp, mks_offset, "rare offset index");
    write_u64_le(fp, carrier_offset, "rare offset index");
    write_u32_le(fp, n_carriers, "rare offset index");
}

static void write_anc_offset_idx_record(
    FILE* fp,
    uint32_t block_id,
    uint64_t mks_offset,
    uint64_t anc_offset
) {
    write_u32_le(fp, block_id, "ancestry offset index");
    write_u64_le(fp, mks_offset, "ancestry offset index");
    write_u64_le(fp, anc_offset, "ancestry offset index");
}

static void write_common_mks_record(
    FILE* fp,
    uint64_t common_index,
    uint32_t global_variant_index,
    const char* chr,
    int64_t pos,
    const std::string& id,
    const char* ref,
    const char* alt,
    uint32_t alt_index,
    uint32_t block_id,
    uint64_t geno_offset,
    uint32_t mac
) {
    write_u64_le(fp, common_index, "common mks");
    write_u32_le(fp, global_variant_index, "common mks");
    write_string(fp, chr, "common mks chr");
    write_i64_le(fp, pos, "common mks pos");
    write_string(fp, id, "common mks id");
    write_string(fp, ref, "common mks ref");
    write_string(fp, alt, "common mks alt");
    write_u32_le(fp, alt_index, "common mks");
    write_u32_le(fp, block_id, "common mks");
    write_u64_le(fp, geno_offset, "common mks");
    write_u32_le(fp, mac, "common mks");
}

static void write_rare_mks_record(
    FILE* fp,
    uint64_t rare_index,
    uint32_t global_variant_index,
    const char* chr,
    int64_t pos,
    const std::string& id,
    const char* ref,
    const char* alt,
    uint32_t alt_index,
    uint64_t carrier_offset,
    uint32_t n_carriers,
    uint32_t mac
) {
    write_u64_le(fp, rare_index, "rare mks");
    write_u32_le(fp, global_variant_index, "rare mks");
    write_string(fp, chr, "rare mks chr");
    write_i64_le(fp, pos, "rare mks pos");
    write_string(fp, id, "rare mks id");
    write_string(fp, ref, "rare mks ref");
    write_string(fp, alt, "rare mks alt");
    write_u32_le(fp, alt_index, "rare mks");
    write_u64_le(fp, carrier_offset, "rare mks");
    write_u32_le(fp, n_carriers, "rare mks");
    write_u32_le(fp, mac, "rare mks");
}

static void write_anc_mks_record(
    FILE* fp,
    uint32_t block_id,
    const char* chr,
    int64_t start_pos,
    int64_t end_pos,
    uint64_t anc_offset
) {
    write_u32_le(fp, block_id, "ancestry mks");
    write_string(fp, chr, "ancestry mks chr");
    write_i64_le(fp, start_pos, "ancestry mks start");
    write_i64_le(fp, end_pos, "ancestry mks end");
    write_u64_le(fp, anc_offset, "ancestry mks offset");
}

static void write_sidecars(
    const std::string& samples_path,
    const std::string& meta_path,
    const char* geno_vcf,
    const char* flare_vcf,
    const std::vector<std::string>& sample_ids,
    uint64_t n_haps,
    int n_words,
    int n_ancestries,
    int rare_threshold,
    const char* selected_region,
    const char* keep_path,
    const char* extract_path
) {
    FILE* samples_fp = open_output_or_die(samples_path, "w");
    for (const std::string& sample_id : sample_ids) {
        std::fprintf(samples_fp, "%s\n", sample_id.c_str());
    }
    std::fclose(samples_fp);

    FILE* meta_fp = open_output_or_die(meta_path, "w");
    std::fprintf(meta_fp, "format_version\t1\n");
    std::fprintf(meta_fp, "n_samples\t%llu\n",
        static_cast<unsigned long long>(sample_ids.size()));
    std::fprintf(meta_fp, "n_haps\t%llu\n", static_cast<unsigned long long>(n_haps));
    std::fprintf(meta_fp, "n_words\t%d\n", n_words);
    std::fprintf(meta_fp, "n_ancestries\t%d\n", n_ancestries);
    std::fprintf(meta_fp, "rare_threshold\t%d\n", rare_threshold);
    std::fprintf(meta_fp, "source_genotype\t%s\n", geno_vcf);
    std::fprintf(meta_fp, "source_flare\t%s\n", flare_vcf);
    if (selected_region) {
        std::fprintf(meta_fp, "selected_region\t%s\n", selected_region);
    }
    if (keep_path) {
        std::fprintf(meta_fp, "keep_samples\t%s\n", keep_path);
    }
    if (extract_path) {
        std::fprintf(meta_fp, "extract_sites\t%s\n", extract_path);
    }
    std::fclose(meta_fp);
}

static void write_anc_block(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    uint32_t block_id,
    const char* chr,
    int64_t start_pos,
    int64_t end_pos,
    const AncestryState& state,
    int n_ancestries,
    int n_words
) {
    uint64_t offset = tell_or_die(anc_fp, anc_bin_path);

    for (int k = 0; k < n_ancestries; ++k) {
        size_t wrote = std::fwrite(
            state.masks[k].data(),
            sizeof(uint64_t),
            static_cast<size_t>(n_words),
            anc_fp
        );

        if (wrote != static_cast<size_t>(n_words)) {
            die("failed writing ancestry block %u", block_id);
        }
    }

    uint64_t mks_offset = tell_or_die(anc_mks_fp, anc_mks_path);
    write_anc_mks_record(anc_mks_fp, block_id, chr, start_pos, end_pos, offset);
    write_anc_offset_idx_record(anc_idx_fp, block_id, mks_offset, offset);
}

static void decode_flare_changes(
    bcf_hdr_t* ahdr,
    bcf1_t* arec,
    int flare_raw_n_samples,
    const std::vector<int>& flare_raw_indices,
    int n_ancestries,
    int flare_geno_rid,
    FlareDeltaDecoder& decoder,
    LaiRecord& out
) {
    bcf_fmt_t* an1 = require_integer_format(ahdr, arec, "AN1");
    bcf_fmt_t* an2 = require_integer_format(ahdr, arec, "AN2");
    if (arec->n_sample != flare_raw_n_samples || an1->n != 1 || an2->n != 1) {
        die("FLARE VCF must contain scalar FORMAT/AN1 and FORMAT/AN2");
    }

    int n_output_samples = static_cast<int>(flare_raw_indices.size());
    size_t n_haps = static_cast<size_t>(2 * n_output_samples);
    out.reset_state = decoder.geno_rid != flare_geno_rid ||
                      decoder.hap_ancestry.size() != n_haps;
    out.changes.clear();
    if (out.reset_state) {
        decoder.geno_rid = flare_geno_rid;
        decoder.hap_ancestry.assign(n_haps, -1);
        out.changes.reserve(n_haps);
    }

    for (int out_i = 0; out_i < n_output_samples; ++out_i) {
        int raw_i = flare_raw_indices[out_i];
        int64_t raw_a1 = integer_format_value(an1, raw_i, 0, "AN1");
        int64_t raw_a2 = integer_format_value(an2, raw_i, 0, "AN2");

        if (integer_format_missing(an1, raw_a1) || integer_format_vector_end(an1, raw_a1)) {
            die("missing FORMAT/AN1 in FLARE record at sample index %d", raw_i);
        }

        if (integer_format_missing(an2, raw_a2) || integer_format_vector_end(an2, raw_a2)) {
            die("missing FORMAT/AN2 in FLARE record at sample index %d", raw_i);
        }

        int a1 = static_cast<int>(raw_a1);
        int a2 = static_cast<int>(raw_a2);

        if (a1 < 0 || a1 >= n_ancestries) {
            die("FORMAT/AN1 value out of range at sample index %d: %d", raw_i, a1);
        }

        if (a2 < 0 || a2 >= n_ancestries) {
            die("FORMAT/AN2 value out of range at sample index %d: %d", raw_i, a2);
        }

        uint32_t hap0 = static_cast<uint32_t>(2 * out_i);
        uint32_t hap1 = static_cast<uint32_t>(2 * out_i + 1);
        int8_t next_a1 = static_cast<int8_t>(a1);
        int8_t next_a2 = static_cast<int8_t>(a2);
        if (decoder.hap_ancestry[hap0] != next_a1) {
            decoder.hap_ancestry[hap0] = next_a1;
            out.changes.push_back(AncestryChange{hap0, next_a1});
        }
        if (decoder.hap_ancestry[hap1] != next_a2) {
            decoder.hap_ancestry[hap1] = next_a2;
            out.changes.push_back(AncestryChange{hap1, next_a2});
        }
    }
}

static void reset_ancestry_state(
    AncestryState& state,
    int n_ancestries,
    int n_words,
    size_t n_haps
) {
    state.masks.resize(static_cast<size_t>(n_ancestries));
    for (std::vector<uint64_t>& mask : state.masks) {
        mask.assign(static_cast<size_t>(n_words), 0);
    }
    state.hap_ancestry.assign(n_haps, -1);
}

static void apply_ancestry_changes(
    AncestryState& state,
    const LaiRecord& record,
    int n_ancestries,
    int n_words,
    size_t n_haps
) {
    if (record.reset_state || state.hap_ancestry.size() != n_haps) {
        reset_ancestry_state(state, n_ancestries, n_words, n_haps);
    }

    for (const AncestryChange& change : record.changes) {
        if (change.hap_id >= n_haps || change.ancestry < 0 ||
            change.ancestry >= n_ancestries) {
            die("internal FLARE ancestry delta out of range");
        }
        int8_t previous = state.hap_ancestry[change.hap_id];
        if (previous >= 0) clear_bit(state.masks[previous], change.hap_id);
        set_bit(state.masks[change.ancestry], change.hap_id);
        state.hap_ancestry[change.hap_id] = change.ancestry;
    }
}

static bool ancestry_record_changes_state(
    const AncestryState& state,
    const LaiRecord& record,
    int current_rid
) {
    if (record.reset_state || current_rid != record.geno_rid) return true;
    if (record.changes.empty()) return false;
    if (!record.merged_duplicate) return true;

    std::vector<int8_t> final_ancestry = state.hap_ancestry;
    for (const AncestryChange& change : record.changes) {
        if (change.hap_id >= final_ancestry.size()) {
            die("internal duplicate-coordinate FLARE delta out of range");
        }
        final_ancestry[change.hap_id] = change.ancestry;
    }
    return final_ancestry != state.hap_ancestry;
}

static void close_open_block(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words
) {
    if (!block.active) return;

    write_anc_block(
        anc_fp,
        anc_bin_path,
        anc_mks_fp,
        anc_mks_path,
        anc_idx_fp,
        block.block_id,
        block.chr.c_str(),
        block.start_pos,
        block.end_pos,
        block.state,
        n_ancestries,
        n_words
    );

    ++n_blocks_written;
    block.active = false;
}

static void update_open_block_from_flare(
    FILE* anc_fp,
    const char* anc_bin_path,
    FILE* anc_mks_fp,
    const char* anc_mks_path,
    FILE* anc_idx_fp,
    OpenAncestryBlock& block,
    uint32_t& n_blocks_written,
    int n_ancestries,
    int n_words,
    int flare_geno_rid,
    const char* flare_chr,
    int64_t interval_start,
    int64_t interval_end,
    bool state_changed
) {
    if (interval_start > interval_end) {
        die(
            "invalid LAI interval for %s: %lld > %lld",
            flare_chr,
            static_cast<long long>(interval_start),
            static_cast<long long>(interval_end)
        );
    }

    bool same_chromosome = block.geno_rid == flare_geno_rid;
    if (block.active && same_chromosome && !state_changed) {
        block.end_pos = std::max(block.end_pos, interval_end);
        return;
    }

    if (block.active) {
        close_open_block(
            anc_fp,
            anc_bin_path,
            anc_mks_fp,
            anc_mks_path,
            anc_idx_fp,
            block,
            n_blocks_written,
            n_ancestries,
            n_words
        );
    }

    block.active = true;
    block.block_id = n_blocks_written;
    block.geno_rid = flare_geno_rid;
    block.chr = flare_chr;
    block.start_pos = interval_start;
    block.end_pos = interval_end;
}

static std::string make_split_id(
    const char* raw_id,
    const char* chr,
    int64_t pos,
    const char* ref,
    const char* alt
) {
    if (raw_id == nullptr || std::strcmp(raw_id, ".") == 0) {
        return std::string(chr) + ":" + std::to_string(pos) + ":" + ref + ":" + alt;
    }

    return std::string(raw_id) + "_" + ref + "_" + alt;
}

static void add_alt_haplotype(
    std::vector<AltCarrierBuilder>& builders_by_alt,
    const std::vector<uint8_t>& selected_alts,
    const AncestryState& state,
    int allele,
    uint32_t hap_id,
    int rare_threshold,
    int n_words
) {
    if (allele <= 0 || allele >= static_cast<int>(builders_by_alt.size()) ||
        !selected_alts[static_cast<size_t>(allele)]) {
        return;
    }

    int8_t ancestry = state.hap_ancestry[hap_id];
    if (ancestry < 0) return;

    AltCarrierBuilder& builder = builders_by_alt[static_cast<size_t>(allele)];
    ++builder.mac;

    if (!builder.dense_bits.empty()) {
        set_bit(builder.dense_bits, hap_id);
        return;
    }

    if (builder.sparse_haps.size() < static_cast<size_t>(rare_threshold)) {
        builder.sparse_haps.push_back(hap_id);
        return;
    }

    builder.dense_bits.assign(static_cast<size_t>(n_words), 0);
    for (uint32_t previous_hap : builder.sparse_haps) {
        set_bit(builder.dense_bits, previous_hap);
    }
    builder.sparse_haps.clear();
    set_bit(builder.dense_bits, hap_id);
}

static inline void add_biallelic_haplotype(
    AltCarrierBuilder& builder,
    const AncestryState& state,
    uint32_t hap_id,
    int rare_threshold,
    int n_words
) {
    if (state.hap_ancestry[hap_id] < 0) return;
    ++builder.mac;

    if (!builder.dense_bits.empty()) {
        set_bit(builder.dense_bits, hap_id);
    } else if (builder.sparse_haps.size() < static_cast<size_t>(rare_threshold)) {
        builder.sparse_haps.push_back(hap_id);
    } else {
        builder.dense_bits.assign(static_cast<size_t>(n_words), 0);
        for (uint32_t previous_hap : builder.sparse_haps) {
            set_bit(builder.dense_bits, previous_hap);
        }
        builder.sparse_haps.clear();
        set_bit(builder.dense_bits, hap_id);
    }
}

static void validate_extra_ploidy(
    const bcf_fmt_t* gt,
    int sample_index,
    const char* chr,
    int64_t pos
) {
    for (int p = 2; p < gt->n; ++p) {
        int64_t value = integer_format_value(gt, sample_index, p, "GT");
        if (integer_format_vector_end(gt, value)) break;
        die(
            "non-diploid genotype at %s:%lld sample index %d",
            chr,
            static_cast<long long>(pos),
            sample_index
        );
    }
}

static void process_genotypes(
    bcf_hdr_t* ghdr,
    bcf1_t* grec,
    const OpenAncestryBlock& block,
    int genotype_raw_n_samples,
    const std::vector<int>& genotype_raw_indices,
    const std::vector<uint8_t>& selected_alts,
    int rare_threshold,
    int n_words,
    std::vector<AltCarrierBuilder>& builders_by_alt
) {
    const char* chr = bcf_hdr_id2name(ghdr, grec->rid);
    int64_t pos = static_cast<int64_t>(grec->pos) + 1;

    bcf_fmt_t* gt = require_integer_format(ghdr, grec, "GT");
    if (grec->n_sample != genotype_raw_n_samples) {
        die("GT sample count mismatch at %s:%lld", chr, static_cast<long long>(pos));
    }
    if (gt->n < 2) {
        die("expected diploid GT at %s:%lld", chr, static_cast<long long>(pos));
    }

    if (selected_alts.size() != static_cast<size_t>(grec->n_allele)) {
        die("internal selected ALT mask mismatch at %s:%lld", chr, static_cast<long long>(pos));
    }

    builders_by_alt.resize(static_cast<size_t>(grec->n_allele));
    for (AltCarrierBuilder& builder : builders_by_alt) {
        builder.mac = 0;
        builder.sparse_haps.clear();
        builder.dense_bits.clear();
    }

    int n_output_samples = static_cast<int>(genotype_raw_indices.size());
    bool biallelic_fast_path = grec->n_allele == 2 && selected_alts[1] != 0;
    for (int out_i = 0; out_i < n_output_samples; ++out_i) {
        int raw_i = genotype_raw_indices[out_i];
        int64_t raw_g0 = integer_format_value(gt, raw_i, 0, "GT");
        int64_t raw_g1 = integer_format_value(gt, raw_i, 1, "GT");

        validate_extra_ploidy(gt, raw_i, chr, pos);

        if (integer_format_vector_end(gt, raw_g0) ||
            integer_format_vector_end(gt, raw_g1) ||
            integer_format_missing(gt, raw_g0) ||
            integer_format_missing(gt, raw_g1) ||
            raw_g0 < 0 || raw_g0 > std::numeric_limits<int32_t>::max() ||
            raw_g1 < 0 || raw_g1 > std::numeric_limits<int32_t>::max()) {
            die(
                "missing genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                raw_i
            );
        }

        int32_t g0 = static_cast<int32_t>(raw_g0);
        int32_t g1 = static_cast<int32_t>(raw_g1);
        if (bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1)) {
            die(
                "missing genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                raw_i
            );
        }

        if (bcf_gt_allele(g0) >= grec->n_allele ||
            bcf_gt_allele(g1) >= grec->n_allele) {
            die(
                "GT allele out of range at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                raw_i
            );
        }

        int allele0 = bcf_gt_allele(g0);
        int allele1 = bcf_gt_allele(g1);

        if (!bcf_gt_is_phased(g1)) {
            die(
                "unphased genotype at %s:%lld sample index %d",
                chr,
                static_cast<long long>(pos),
                raw_i
            );
        }

        if (biallelic_fast_path) {
            AltCarrierBuilder& builder = builders_by_alt[1];
            if (allele0 == 1) {
                add_biallelic_haplotype(
                    builder,
                    block.state,
                    static_cast<uint32_t>(2 * out_i),
                    rare_threshold,
                    n_words
                );
            }
            if (allele1 == 1) {
                add_biallelic_haplotype(
                    builder,
                    block.state,
                    static_cast<uint32_t>(2 * out_i + 1),
                    rare_threshold,
                    n_words
                );
            }
        } else {
            add_alt_haplotype(
                builders_by_alt,
                selected_alts,
                block.state,
                allele0,
                static_cast<uint32_t>(2 * out_i),
                rare_threshold,
                n_words
            );
            add_alt_haplotype(
                builders_by_alt,
                selected_alts,
                block.state,
                allele1,
                static_cast<uint32_t>(2 * out_i + 1),
                rare_threshold,
                n_words
            );
        }
    }
}

static std::vector<std::string> split_region_queries(const std::string& query) {
    return query.empty() ? std::vector<std::string>{} : split_commas(query);
}

static void init_fast_text_reader(
    FastTextVcfReader& reader,
    htsFile* sequential_fp,
    const char* path,
    const std::string& indexed_query
) {
    if (indexed_query.empty()) {
        reader.fp = sequential_fp;
        return;
    }

    reader.fp = hts_open(path, "r");
    reader.tbx = tbx_index_load3(path, nullptr, HTS_IDX_SILENT_FAIL);
    if (!reader.fp || !reader.tbx) {
        if (reader.fp) hts_close(reader.fp);
        if (reader.tbx) tbx_destroy(reader.tbx);
        reader = FastTextVcfReader{};
        reader.fp = sequential_fp;
        return;
    }
    reader.owns_fp = true;
    reader.queries = split_region_queries(indexed_query);
    hts_set_opt(reader.fp, HTS_OPT_BLOCK_SIZE, kInputBlockSize);
}

static int fast_text_next_line(FastTextVcfReader& reader) {
    if (!reader.tbx) {
        return hts_getline(reader.fp, '\n', &reader.line) >= 0 ? 0 : -1;
    }

    while (true) {
        if (reader.itr && tbx_itr_next(reader.fp, reader.tbx, reader.itr, &reader.line) >= 0) {
            return 0;
        }
        if (reader.itr) {
            tbx_itr_destroy(reader.itr);
            reader.itr = nullptr;
        }
        if (reader.query_index >= reader.queries.size()) return -1;
        reader.itr = tbx_itr_querys(
            reader.tbx,
            reader.queries[reader.query_index++].c_str()
        );
        if (!reader.itr) continue;
    }
}

static void destroy_fast_text_reader(FastTextVcfReader& reader) {
    if (reader.itr) tbx_itr_destroy(reader.itr);
    if (reader.tbx) tbx_destroy(reader.tbx);
    if (reader.owns_fp && reader.fp) hts_close(reader.fp);
    std::free(reader.line.s);
    reader = FastTextVcfReader{};
}

static char* find_and_terminate_tab(char* cursor, const char* label) {
    char* tab = std::strchr(cursor, '\t');
    if (!tab) die("malformed VCF record: missing %s column separator", label);
    *tab = '\0';
    return tab + 1;
}

static int find_format_field_index(const char* format, const char* tag) {
    int index = 0;
    const char* start = format;
    size_t tag_len = std::strlen(tag);
    while (*start) {
        const char* end = std::strchr(start, ':');
        size_t len = end ? static_cast<size_t>(end - start) : std::strlen(start);
        if (len == tag_len && std::memcmp(start, tag, len) == 0) return index;
        if (!end) break;
        start = end + 1;
        ++index;
    }
    return -1;
}

static void split_fast_vcf_fixed_fields(char* line, char* fields[9], char*& samples) {
    char* cursor = line;
    for (int i = 0; i < 9; ++i) {
        fields[i] = cursor;
        if (i < 8) cursor = find_and_terminate_tab(cursor, "fixed VCF");
    }
    samples = find_and_terminate_tab(fields[8], "FORMAT");
}

static int64_t parse_fast_positive_pos(const char* text) {
    char* end = nullptr;
    long long value = std::strtoll(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0) {
        die("invalid VCF POS: %s", text);
    }
    return static_cast<int64_t>(value);
}

static void parse_fast_genotype_record(
    FastTextVcfReader& reader,
    bcf_hdr_t* ghdr,
    FastGenotypeRecord& record
) {
    char* fields[9];
    split_fast_vcf_fixed_fields(reader.line.s, fields, record.samples);
    record.samples_length = reader.line.l -
        static_cast<size_t>(record.samples - reader.line.s);
    record.chr = fields[0];
    record.rid = bcf_hdr_name2id(ghdr, record.chr);
    record.pos = parse_fast_positive_pos(fields[1]);
    record.id = fields[2];
    record.ref = fields[3];
    record.gt_field_index = find_format_field_index(fields[8], "GT");
    record.gt_only = std::strcmp(fields[8], "GT") == 0;
    if (record.gt_field_index < 0) {
        die("missing FORMAT/GT at %s:%lld", record.chr, static_cast<long long>(record.pos));
    }

    record.alleles.clear();
    record.alleles.push_back(record.ref);
    char* alt = fields[4];
    while (true) {
        record.alleles.push_back(alt);
        char* comma = std::strchr(alt, ',');
        if (!comma) break;
        *comma = '\0';
        alt = comma + 1;
    }
}

static const char* locate_sample_subfield(
    const char* sample_start,
    const char* sample_end,
    int field_index,
    const char*& field_end
) {
    const char* start = sample_start;
    for (int field = 0; field < field_index; ++field) {
        const char* colon = static_cast<const char*>(
            std::memchr(start, ':', static_cast<size_t>(sample_end - start))
        );
        if (!colon) return nullptr;
        start = colon + 1;
    }
    const char* colon = static_cast<const char*>(
        std::memchr(start, ':', static_cast<size_t>(sample_end - start))
    );
    field_end = colon ? colon : sample_end;
    return start;
}

static int parse_fast_nonnegative_int(const char* start, const char* end, const char* label) {
    if (start >= end || *start == '.') die("missing %s", label);
    int value = 0;
    for (const char* p = start; p < end; ++p) {
        if (*p < '0' || *p > '9') die("invalid %s value", label);
        if (value > (std::numeric_limits<int>::max() - (*p - '0')) / 10) {
            die("%s value out of range", label);
        }
        value = value * 10 + (*p - '0');
    }
    return value;
}

static void parse_fast_phased_gt(
    const char* start,
    const char* end,
    int& allele0,
    int& allele1,
    const char* chr,
    int64_t pos,
    int sample_index
) {
    const char* separator = static_cast<const char*>(
        std::memchr(start, '|', static_cast<size_t>(end - start))
    );
    if (!separator) {
        if (std::memchr(start, '/', static_cast<size_t>(end - start))) {
            die("unphased genotype at %s:%lld sample index %d", chr,
                static_cast<long long>(pos), sample_index);
        }
        die("expected diploid GT at %s:%lld", chr, static_cast<long long>(pos));
    }
    if (std::memchr(separator + 1, '|', static_cast<size_t>(end - separator - 1)) ||
        std::memchr(separator + 1, '/', static_cast<size_t>(end - separator - 1))) {
        die("non-diploid genotype at %s:%lld sample index %d", chr,
            static_cast<long long>(pos), sample_index);
    }
    allele0 = parse_fast_nonnegative_int(start, separator, "genotype");
    allele1 = parse_fast_nonnegative_int(separator + 1, end, "genotype");
}

static void reset_carrier_builders(
    std::vector<AltCarrierBuilder>& builders_by_alt,
    size_t n_alleles
) {
    builders_by_alt.resize(n_alleles);
    for (AltCarrierBuilder& builder : builders_by_alt) {
        builder.mac = 0;
        builder.sparse_haps.clear();
        builder.dense_bits.clear();
    }
}

static inline void append_carrier_word(
    AltCarrierBuilder& builder,
    uint64_t carrier_word,
    size_t word_index,
    int rare_threshold,
    int n_words
) {
    if (carrier_word == 0) return;

    uint32_t word_mac = static_cast<uint32_t>(__builtin_popcountll(carrier_word));
    builder.mac += word_mac;
    if (!builder.dense_bits.empty()) {
        builder.dense_bits[word_index] = carrier_word;
        return;
    }

    if (builder.mac <= static_cast<uint32_t>(rare_threshold)) {
        uint64_t remaining = carrier_word;
        uint32_t word_hap = static_cast<uint32_t>(word_index * 64);
        while (remaining != 0) {
            uint32_t bit = static_cast<uint32_t>(__builtin_ctzll(remaining));
            builder.sparse_haps.push_back(word_hap + bit);
            remaining &= remaining - 1;
        }
        return;
    }

    builder.dense_bits.assign(static_cast<size_t>(n_words), 0);
    for (uint32_t previous_hap : builder.sparse_haps) {
        set_bit(builder.dense_bits, previous_hap);
    }
    builder.sparse_haps.clear();
    builder.dense_bits[word_index] = carrier_word;
}

static bool process_compact_gt_words(
    const FastGenotypeRecord& record,
    int genotype_text_n_samples,
    const std::vector<int>& genotype_text_indices,
    const std::vector<uint8_t>& selected_alts,
    int rare_threshold,
    int n_words,
    std::vector<AltCarrierBuilder>& builders_by_alt
) {
    if (!record.gt_only) return false;
    size_t expected_length = genotype_text_n_samples > 0
        ? static_cast<size_t>(genotype_text_n_samples) * 4 - 1
        : 0;
    if (record.samples_length != expected_length) return false;

    bool identity_sample_order =
        genotype_text_indices.size() == static_cast<size_t>(genotype_text_n_samples) &&
        !genotype_text_indices.empty() &&
        genotype_text_indices.front() == 0 &&
        genotype_text_indices.back() == genotype_text_n_samples - 1;
    bool biallelic = record.alleles.size() == 2 && selected_alts[1] != 0;
    if (identity_sample_order) {
        const char* cursor = record.samples;
        int out_i = 0;
        if (biallelic) {
            AltCarrierBuilder& builder = builders_by_alt[1];
            for (size_t word_index = 0; word_index < static_cast<size_t>(n_words);
                 ++word_index) {
                int samples_in_word = std::min(32, genotype_text_n_samples - out_i);
                uint64_t carrier_word = 0;
                for (int sample_in_word = 0; sample_in_word < samples_in_word;
                     ++sample_in_word, ++out_i) {
                    char expected_delimiter =
                        out_i + 1 == genotype_text_n_samples ? '\0' : '\t';
                    if (cursor[0] < '0' || cursor[0] > '9' || cursor[1] != '|' ||
                        cursor[2] < '0' || cursor[2] > '9' ||
                        cursor[3] != expected_delimiter) {
                        return false;
                    }
                    int allele0 = cursor[0] - '0';
                    int allele1 = cursor[2] - '0';
                    if (allele0 >= 2 || allele1 >= 2) {
                        die("GT allele out of range at %s:%lld sample index %d", record.chr,
                            static_cast<long long>(record.pos), out_i);
                    }
                    carrier_word |= static_cast<uint64_t>(allele0 == 1)
                                    << (2 * sample_in_word);
                    carrier_word |= static_cast<uint64_t>(allele1 == 1)
                                    << (2 * sample_in_word + 1);
                    cursor += expected_delimiter == '\t' ? 4 : 3;
                }
                append_carrier_word(
                    builder, carrier_word, word_index, rare_threshold, n_words);
            }
            return true;
        }

        std::vector<uint64_t> words_by_alt(record.alleles.size(), 0);
        for (size_t word_index = 0; word_index < static_cast<size_t>(n_words);
             ++word_index) {
            std::fill(words_by_alt.begin(), words_by_alt.end(), 0);
            int samples_in_word = std::min(32, genotype_text_n_samples - out_i);
            for (int sample_in_word = 0; sample_in_word < samples_in_word;
                 ++sample_in_word, ++out_i) {
                char expected_delimiter =
                    out_i + 1 == genotype_text_n_samples ? '\0' : '\t';
                if (cursor[0] < '0' || cursor[0] > '9' || cursor[1] != '|' ||
                    cursor[2] < '0' || cursor[2] > '9' ||
                    cursor[3] != expected_delimiter) {
                    return false;
                }
                int allele0 = cursor[0] - '0';
                int allele1 = cursor[2] - '0';
                if (allele0 >= static_cast<int>(record.alleles.size()) ||
                    allele1 >= static_cast<int>(record.alleles.size())) {
                    die("GT allele out of range at %s:%lld sample index %d", record.chr,
                        static_cast<long long>(record.pos), out_i);
                }
                if (allele0 > 0 && selected_alts[static_cast<size_t>(allele0)] != 0) {
                    words_by_alt[static_cast<size_t>(allele0)] |=
                        uint64_t{1} << (2 * sample_in_word);
                }
                if (allele1 > 0 && selected_alts[static_cast<size_t>(allele1)] != 0) {
                    words_by_alt[static_cast<size_t>(allele1)] |=
                        uint64_t{1} << (2 * sample_in_word + 1);
                }
                cursor += expected_delimiter == '\t' ? 4 : 3;
            }
            for (size_t alt_index = 1; alt_index < words_by_alt.size(); ++alt_index) {
                if (selected_alts[alt_index] == 0) continue;
                append_carrier_word(
                    builders_by_alt[alt_index],
                    words_by_alt[alt_index],
                    word_index,
                    rare_threshold,
                    n_words
                );
            }
        }
        return true;
    }

    int out_i = 0;
    int n_output_samples = static_cast<int>(genotype_text_indices.size());
    if (biallelic) {
        AltCarrierBuilder& builder = builders_by_alt[1];
        for (size_t word_index = 0; word_index < static_cast<size_t>(n_words);
             ++word_index) {
            int samples_in_word = std::min(32, n_output_samples - out_i);
            uint64_t carrier_word = 0;
            for (int sample_in_word = 0; sample_in_word < samples_in_word;
                 ++sample_in_word, ++out_i) {
                int raw_i = genotype_text_indices[static_cast<size_t>(out_i)];
                const char* sample = record.samples + static_cast<size_t>(raw_i) * 4;
                char expected_delimiter =
                    raw_i + 1 == genotype_text_n_samples ? '\0' : '\t';
                if (sample[0] < '0' || sample[0] > '9' || sample[1] != '|' ||
                    sample[2] < '0' || sample[2] > '9' ||
                    sample[3] != expected_delimiter) {
                    return false;
                }
                int allele0 = sample[0] - '0';
                int allele1 = sample[2] - '0';
                if (allele0 >= 2 || allele1 >= 2) {
                    die("GT allele out of range at %s:%lld sample index %d", record.chr,
                        static_cast<long long>(record.pos), raw_i);
                }
                carrier_word |= static_cast<uint64_t>(allele0 == 1)
                                << (2 * sample_in_word);
                carrier_word |= static_cast<uint64_t>(allele1 == 1)
                                << (2 * sample_in_word + 1);
            }
            append_carrier_word(
                builder, carrier_word, word_index, rare_threshold, n_words);
        }
        return true;
    }

    std::vector<uint64_t> words_by_alt(record.alleles.size(), 0);
    for (size_t word_index = 0; word_index < static_cast<size_t>(n_words);
         ++word_index) {
        std::fill(words_by_alt.begin(), words_by_alt.end(), 0);
        int samples_in_word = std::min(32, n_output_samples - out_i);
        for (int sample_in_word = 0; sample_in_word < samples_in_word;
             ++sample_in_word, ++out_i) {
            int raw_i = genotype_text_indices[static_cast<size_t>(out_i)];
            const char* sample = record.samples + static_cast<size_t>(raw_i) * 4;
            char expected_delimiter = raw_i + 1 == genotype_text_n_samples ? '\0' : '\t';
            if (sample[0] < '0' || sample[0] > '9' || sample[1] != '|' ||
                sample[2] < '0' || sample[2] > '9' || sample[3] != expected_delimiter) {
                return false;
            }
            int allele0 = sample[0] - '0';
            int allele1 = sample[2] - '0';
            if (allele0 >= static_cast<int>(record.alleles.size()) ||
                allele1 >= static_cast<int>(record.alleles.size())) {
                die("GT allele out of range at %s:%lld sample index %d", record.chr,
                    static_cast<long long>(record.pos), raw_i);
            }
            if (allele0 > 0 && selected_alts[static_cast<size_t>(allele0)] != 0) {
                words_by_alt[static_cast<size_t>(allele0)] |=
                    uint64_t{1} << (2 * sample_in_word);
            }
            if (allele1 > 0 && selected_alts[static_cast<size_t>(allele1)] != 0) {
                words_by_alt[static_cast<size_t>(allele1)] |=
                    uint64_t{1} << (2 * sample_in_word + 1);
            }
        }
        for (size_t alt_index = 1; alt_index < words_by_alt.size(); ++alt_index) {
            if (selected_alts[alt_index] == 0) continue;
            append_carrier_word(
                builders_by_alt[alt_index],
                words_by_alt[alt_index],
                word_index,
                rare_threshold,
                n_words
            );
        }
    }
    return true;
}

static void process_fast_text_genotypes(
    const FastGenotypeRecord& record,
    const OpenAncestryBlock& block,
    int genotype_text_n_samples,
    const std::vector<int>& genotype_text_indices,
    const std::vector<uint8_t>& selected_alts,
    int rare_threshold,
    int n_words,
    std::vector<AltCarrierBuilder>& builders_by_alt
) {
    if (selected_alts.size() != record.alleles.size()) {
        die("internal selected ALT mask mismatch at %s:%lld", record.chr,
            static_cast<long long>(record.pos));
    }
    reset_carrier_builders(builders_by_alt, record.alleles.size());
    bool biallelic_fast_path = record.alleles.size() == 2 && selected_alts[1] != 0;

    bool identity_sample_order =
        genotype_text_indices.size() == static_cast<size_t>(genotype_text_n_samples) &&
        !genotype_text_indices.empty() &&
        genotype_text_indices.front() == 0 &&
        genotype_text_indices.back() == genotype_text_n_samples - 1;
    if (process_compact_gt_words(
            record,
            genotype_text_n_samples,
            genotype_text_indices,
            selected_alts,
            rare_threshold,
            n_words,
            builders_by_alt)) {
        return;
    }
    reset_carrier_builders(builders_by_alt, record.alleles.size());
    if (identity_sample_order && record.gt_field_index == 0) {
        const char* cursor = record.samples;
        for (int out_i = 0; out_i < genotype_text_n_samples; ++out_i) {
            int allele0 = 0;
            int allele1 = 0;
            const char* sample_end = nullptr;
            char delimiter = '\0';
            bool compact_gt =
                cursor[0] >= '0' && cursor[0] <= '9' &&
                cursor[1] == '|' &&
                cursor[2] >= '0' && cursor[2] <= '9' &&
                (cursor[3] == '\t' || cursor[3] == ':' || cursor[3] == '\0');
            if (compact_gt) {
                delimiter = cursor[3];
                allele0 = cursor[0] - '0';
                allele1 = cursor[2] - '0';
            } else {
                sample_end = std::strchr(cursor, '\t');
                if (!sample_end) sample_end = cursor + std::strlen(cursor);
                const char* gt_end = nullptr;
                const char* gt_start = locate_sample_subfield(cursor, sample_end, 0, gt_end);
                parse_fast_phased_gt(
                    gt_start,
                    gt_end,
                    allele0,
                    allele1,
                    record.chr,
                    record.pos,
                    out_i
                );
                delimiter = *gt_end;
            }

            if (allele0 >= static_cast<int>(record.alleles.size()) ||
                allele1 >= static_cast<int>(record.alleles.size())) {
                die("GT allele out of range at %s:%lld sample index %d", record.chr,
                    static_cast<long long>(record.pos), out_i);
            }
            if (biallelic_fast_path) {
                AltCarrierBuilder& builder = builders_by_alt[1];
                if (allele0 == 1) {
                    add_biallelic_haplotype(builder, block.state,
                        static_cast<uint32_t>(2 * out_i), rare_threshold, n_words);
                }
                if (allele1 == 1) {
                    add_biallelic_haplotype(builder, block.state,
                        static_cast<uint32_t>(2 * out_i + 1), rare_threshold, n_words);
                }
            } else {
                add_alt_haplotype(builders_by_alt, selected_alts, block.state, allele0,
                    static_cast<uint32_t>(2 * out_i), rare_threshold, n_words);
                add_alt_haplotype(builders_by_alt, selected_alts, block.state, allele1,
                    static_cast<uint32_t>(2 * out_i + 1), rare_threshold, n_words);
            }

            if (compact_gt && delimiter == ':') {
                sample_end = std::strchr(cursor + 4, '\t');
                if (!sample_end) sample_end = cursor + std::strlen(cursor);
            }
            bool has_next_sample = compact_gt
                ? (delimiter == '\t' || (delimiter == ':' && *sample_end == '\t'))
                : *sample_end == '\t';
            if ((out_i + 1 == genotype_text_n_samples) == has_next_sample) {
                die("GT sample count mismatch at %s:%lld", record.chr,
                    static_cast<long long>(record.pos));
            }
            if (!compact_gt) {
                cursor = *sample_end ? sample_end + 1 : sample_end;
            } else if (delimiter == '\t') {
                cursor += 4;
            } else if (delimiter == '\0') {
                if (out_i + 1 != genotype_text_n_samples) {
                    die("GT sample count mismatch at %s:%lld", record.chr,
                        static_cast<long long>(record.pos));
                }
                cursor += 3;
            } else {
                cursor = *sample_end ? sample_end + 1 : sample_end;
            }
        }
        if (*cursor != '\0') {
            die("GT sample count mismatch at %s:%lld", record.chr,
                static_cast<long long>(record.pos));
        }
        return;
    }

    const char* sample_cursor = record.samples;
    int raw_cursor = 0;
    for (size_t out_i = 0; out_i < genotype_text_indices.size(); ++out_i) {
        int target_raw = genotype_text_indices[out_i];
        while (raw_cursor < target_raw) {
            const char* tab = std::strchr(sample_cursor, '\t');
            if (!tab) die("GT sample count mismatch at %s:%lld", record.chr,
                static_cast<long long>(record.pos));
            sample_cursor = tab + 1;
            ++raw_cursor;
        }
        const char* sample_end = std::strchr(sample_cursor, '\t');
        if (!sample_end) sample_end = sample_cursor + std::strlen(sample_cursor);
        const char* gt_end = nullptr;
        const char* gt_start = locate_sample_subfield(
            sample_cursor,
            sample_end,
            record.gt_field_index,
            gt_end
        );
        if (!gt_start) {
            die("missing FORMAT/GT at %s:%lld sample index %d", record.chr,
                static_cast<long long>(record.pos), target_raw);
        }
        int allele0 = 0;
        int allele1 = 0;
        if (gt_end - gt_start == 3 && gt_start[1] == '|' &&
            gt_start[0] >= '0' && gt_start[0] <= '9' &&
            gt_start[2] >= '0' && gt_start[2] <= '9') {
            allele0 = gt_start[0] - '0';
            allele1 = gt_start[2] - '0';
        } else {
            parse_fast_phased_gt(
                gt_start,
                gt_end,
                allele0,
                allele1,
                record.chr,
                record.pos,
                target_raw
            );
        }
        if (allele0 >= static_cast<int>(record.alleles.size()) ||
            allele1 >= static_cast<int>(record.alleles.size())) {
            die("GT allele out of range at %s:%lld sample index %d", record.chr,
                static_cast<long long>(record.pos), target_raw);
        }

        if (biallelic_fast_path) {
            AltCarrierBuilder& builder = builders_by_alt[1];
            if (allele0 == 1) {
                add_biallelic_haplotype(builder, block.state,
                    static_cast<uint32_t>(2 * out_i), rare_threshold, n_words);
            }
            if (allele1 == 1) {
                add_biallelic_haplotype(builder, block.state,
                    static_cast<uint32_t>(2 * out_i + 1), rare_threshold, n_words);
            }
        } else {
            add_alt_haplotype(builders_by_alt, selected_alts, block.state, allele0,
                static_cast<uint32_t>(2 * out_i), rare_threshold, n_words);
            add_alt_haplotype(builders_by_alt, selected_alts, block.state, allele1,
                static_cast<uint32_t>(2 * out_i + 1), rare_threshold, n_words);
        }

        bool has_next_sample = *sample_end == '\t';
        if (target_raw + 1 == genotype_text_n_samples && has_next_sample) {
            die("GT sample count mismatch at %s:%lld", record.chr,
                static_cast<long long>(record.pos));
        }
        sample_cursor = *sample_end ? sample_end + 1 : sample_end;
        raw_cursor = target_raw + 1;
    }

    while (raw_cursor < genotype_text_n_samples && *sample_cursor != '\0') {
        const char* tab = std::strchr(sample_cursor, '\t');
        if (tab && raw_cursor + 1 == genotype_text_n_samples) {
            die("GT sample count mismatch at %s:%lld", record.chr,
                static_cast<long long>(record.pos));
        }
        ++raw_cursor;
        if (!tab) {
            sample_cursor += std::strlen(sample_cursor);
            break;
        }
        sample_cursor = tab + 1;
    }
    if (raw_cursor != genotype_text_n_samples || *sample_cursor != '\0') {
        die("GT sample count mismatch at %s:%lld", record.chr,
            static_cast<long long>(record.pos));
    }
}

static void decode_fast_flare_samples(
    const char* samples,
    int an1_index,
    int an2_index,
    const std::vector<int>& raw_to_output,
    int n_output_samples,
    int n_ancestries,
    int flare_geno_rid,
    FlareDeltaDecoder& decoder,
    LaiRecord& out
) {
    size_t n_haps = static_cast<size_t>(2 * n_output_samples);
    out.reset_state = decoder.geno_rid != flare_geno_rid ||
                      decoder.hap_ancestry.size() != n_haps;
    out.changes.clear();
    if (out.reset_state) {
        decoder.geno_rid = flare_geno_rid;
        decoder.hap_ancestry.assign(n_haps, -1);
        out.changes.reserve(n_haps);
    }

    const char* cursor = samples;
    for (size_t raw_i = 0; raw_i < raw_to_output.size(); ++raw_i) {
        int out_i = raw_to_output[raw_i];
        const char* sample_end = nullptr;
        const char* next_cursor = nullptr;
        bool has_next_sample = false;
        if (out_i >= 0) {
            int a1 = 0;
            int a2 = 0;
            bool compact_ancestry_first =
                an1_index == 0 && an2_index == 1 &&
                cursor[0] >= '0' && cursor[0] <= '9' &&
                cursor[1] == ':' &&
                cursor[2] >= '0' && cursor[2] <= '9' &&
                (cursor[3] == '\t' || cursor[3] == ':' || cursor[3] == '\0');
            bool compact_flare_layout =
                an1_index == 1 && an2_index == 2 &&
                cursor[0] >= '0' && cursor[0] <= '9' &&
                cursor[1] == '|' &&
                cursor[2] >= '0' && cursor[2] <= '9' &&
                cursor[3] == ':' &&
                cursor[4] >= '0' && cursor[4] <= '9' &&
                cursor[5] == ':' &&
                cursor[6] >= '0' && cursor[6] <= '9' &&
                (cursor[7] == '\t' || cursor[7] == ':' || cursor[7] == '\0');
            if (compact_ancestry_first) {
                a1 = cursor[0] - '0';
                a2 = cursor[2] - '0';
                if (cursor[3] == '\t') {
                    has_next_sample = true;
                    next_cursor = cursor + 4;
                } else if (cursor[3] == '\0') {
                    next_cursor = cursor + 3;
                } else {
                    sample_end = std::strchr(cursor + 4, '\t');
                    has_next_sample = sample_end != nullptr;
                    next_cursor = sample_end ? sample_end + 1
                                             : cursor + std::strlen(cursor);
                }
            } else if (compact_flare_layout) {
                a1 = cursor[4] - '0';
                a2 = cursor[6] - '0';
                if (cursor[7] == '\t') {
                    has_next_sample = true;
                    next_cursor = cursor + 8;
                } else if (cursor[7] == '\0') {
                    next_cursor = cursor + 7;
                } else {
                    sample_end = std::strchr(cursor + 8, '\t');
                    has_next_sample = sample_end != nullptr;
                    next_cursor = sample_end ? sample_end + 1
                                             : cursor + std::strlen(cursor);
                }
            } else {
                sample_end = std::strchr(cursor, '\t');
                if (!sample_end) sample_end = cursor + std::strlen(cursor);
                const char* an1_end = nullptr;
                const char* an2_end = nullptr;
                const char* an1_start = locate_sample_subfield(
                    cursor, sample_end, an1_index, an1_end);
                const char* an2_start = locate_sample_subfield(
                    cursor, sample_end, an2_index, an2_end);
                if (!an1_start || !an2_start) {
                    die("missing FORMAT/AN1 or FORMAT/AN2 in FLARE record at sample index %llu",
                        static_cast<unsigned long long>(raw_i));
                }
                a1 = parse_fast_nonnegative_int(an1_start, an1_end, "FORMAT/AN1");
                a2 = parse_fast_nonnegative_int(an2_start, an2_end, "FORMAT/AN2");
                has_next_sample = *sample_end == '\t';
                next_cursor = *sample_end ? sample_end + 1 : sample_end;
            }
            if (a1 >= n_ancestries) {
                die("FORMAT/AN1 value out of range at sample index %llu: %d",
                    static_cast<unsigned long long>(raw_i), a1);
            }
            if (a2 >= n_ancestries) {
                die("FORMAT/AN2 value out of range at sample index %llu: %d",
                    static_cast<unsigned long long>(raw_i), a2);
            }
            uint32_t hap0 = static_cast<uint32_t>(2 * out_i);
            uint32_t hap1 = hap0 + 1;
            int8_t next_a1 = static_cast<int8_t>(a1);
            int8_t next_a2 = static_cast<int8_t>(a2);
            if (decoder.hap_ancestry[hap0] != next_a1) {
                decoder.hap_ancestry[hap0] = next_a1;
                out.changes.push_back(AncestryChange{hap0, next_a1});
            }
            if (decoder.hap_ancestry[hap1] != next_a2) {
                decoder.hap_ancestry[hap1] = next_a2;
                out.changes.push_back(AncestryChange{hap1, next_a2});
            }
        } else {
            sample_end = std::strchr(cursor, '\t');
            if (!sample_end) sample_end = cursor + std::strlen(cursor);
            has_next_sample = *sample_end == '\t';
            next_cursor = *sample_end ? sample_end + 1 : sample_end;
        }

        if ((raw_i + 1 == raw_to_output.size()) == has_next_sample) {
            die("FLARE sample count mismatch");
        }
        cursor = next_cursor;
    }
}

static bool read_next_fast_lai_record(
    FastTextVcfReader& reader,
    bcf_hdr_t* ghdr,
    const std::vector<int>& raw_to_output,
    int n_output_samples,
    int n_ancestries,
    FlareDeltaDecoder& decoder,
    std::unordered_set<std::string>& warned_missing_flare_contigs,
    LaiRecord& out
) {
    while (fast_text_next_line(reader) == 0) {
        char* fields[9];
        char* samples = nullptr;
        split_fast_vcf_fixed_fields(reader.line.s, fields, samples);
        int geno_rid = bcf_hdr_name2id(ghdr, fields[0]);
        if (geno_rid < 0) {
            if (warned_missing_flare_contigs.insert(fields[0]).second) {
                std::fprintf(stderr,
                    "WARNING: skipping FLARE contig absent from genotype header: %s\n",
                    fields[0]);
            }
            continue;
        }
        int an1_index = find_format_field_index(fields[8], "AN1");
        int an2_index = find_format_field_index(fields[8], "AN2");
        if (an1_index < 0 || an2_index < 0) {
            die("FLARE VCF must contain scalar FORMAT/AN1 and FORMAT/AN2");
        }
        out.valid = true;
        out.merged_duplicate = false;
        out.geno_rid = geno_rid;
        out.chr = fields[0];
        out.pos = parse_fast_positive_pos(fields[1]);
        decode_fast_flare_samples(
            samples,
            an1_index,
            an2_index,
            raw_to_output,
            n_output_samples,
            n_ancestries,
            geno_rid,
            decoder,
            out
        );
        return true;
    }
    out = LaiRecord{};
    return false;
}

static int read_next_fast_genotype_record(
    FastTextVcfReader& reader,
    bcf_hdr_t* ghdr,
    FastGenotypeRecord& record,
    const Region& region
) {
    while (fast_text_next_line(reader) == 0) {
        parse_fast_genotype_record(reader, ghdr, record);
        if (record.rid < 0) {
            if (region.active) continue;
            die("VCF contig absent from genotype header: %s", record.chr);
        }
        if (!region.active ||
            (record.rid == region.geno_rid &&
             record.pos >= region.start && record.pos <= region.end)) {
            return 0;
        }
    }
    return -1;
}

static int read_next_record(htsFile* fp, bcf_hdr_t* hdr, bcf1_t* rec, hts_itr_t* itr = nullptr) {
    rec->max_unpack = BCF_UN_STR | BCF_UN_FMT;
    int ret = itr ? bcf_itr_next(fp, itr, rec) : bcf_read(fp, hdr, rec);
    if (ret == 0) bcf_unpack(rec, BCF_UN_STR);
    return ret;
}

static bool add_contig_to_header_if_missing(
    bcf_hdr_t* hdr,
    const std::string& contig,
    const char* header_label
) {
    if (contig.empty()) return false;
    if (bcf_hdr_name2id(hdr, contig.c_str()) >= 0) return false;

    std::string line = "##contig=<ID=" + contig + ">";
    if (bcf_hdr_append(hdr, line.c_str()) != 0 || bcf_hdr_sync(hdr) != 0) {
        die(
            "failed to add inferred contig to %s header: %s",
            header_label,
            contig.c_str()
        );
    }
    return true;
}

static int add_contigs_from_header_if_missing(
    bcf_hdr_t* dst,
    bcf_hdr_t* src,
    const char* dst_label,
    const char* src_label
) {
    int added = 0;
    int n_contigs = src ? src->n[BCF_DT_CTG] : 0;
    for (int i = 0; i < n_contigs; ++i) {
        const char* contig = bcf_hdr_id2name(src, i);
        if (!contig || !*contig) continue;
        if (add_contig_to_header_if_missing(dst, contig, dst_label)) {
            ++added;
        }
    }
    if (added > 0) {
        std::fprintf(
            stderr,
            "Added %d %s contig(s) to %s header.\n",
            added,
            src_label,
            dst_label
        );
    }
    return added;
}

static int add_contigs_from_extract_if_missing(
    bcf_hdr_t* dst,
    const ExtractSites& extract_sites,
    const char* dst_label
) {
    if (!extract_sites.active) return 0;

    std::unordered_set<std::string> seen;
    int added = 0;
    for (const ExtractPosition& position : extract_sites.positions) {
        if (!seen.insert(position.chr).second) continue;
        if (add_contig_to_header_if_missing(dst, position.chr, dst_label)) {
            ++added;
        }
    }

    if (added > 0) {
        std::fprintf(
            stderr,
            "Added %d --extract contig(s) to %s header.\n",
            added,
            dst_label
        );
    }
    return added;
}

static bool extract_position_in_region(
    const ExtractPosition& position,
    const Region& region
) {
    if (!region.active) return true;
    return position.chr == region.chr &&
           position.pos >= region.start &&
           position.pos <= region.end;
}

static int extract_position_sort_rid(bcf_hdr_t* hdr, const std::string& chr) {
    int rid = bcf_hdr_name2id(hdr, chr.c_str());
    if (rid >= 0) return rid;
    return std::numeric_limits<int>::max();
}

static bool extract_position_less(const ExtractPosition& a, const ExtractPosition& b) {
    if (a.geno_rid != b.geno_rid) return a.geno_rid < b.geno_rid;
    if (a.chr != b.chr) return a.chr < b.chr;
    if (a.pos != b.pos) return a.pos < b.pos;
    return a.ref < b.ref;
}

static void merge_extract_alts(
    ExtractPosition& destination,
    const ExtractPosition& source,
    const std::string& path
) {
    std::vector<std::string> merged;
    merged.reserve(destination.alts.size() + source.alts.size());

    size_t destination_i = 0;
    size_t source_i = 0;
    while (destination_i < destination.alts.size() && source_i < source.alts.size()) {
        const std::string& destination_alt = destination.alts[destination_i];
        const std::string& source_alt = source.alts[source_i];
        if (destination_alt < source_alt) {
            merged.push_back(destination_alt);
            ++destination_i;
        } else if (source_alt < destination_alt) {
            merged.push_back(source_alt);
            ++source_i;
        } else {
            die(
                "duplicate allele in --extract list at %s:%llu: %s:%lld %s>%s",
                path.c_str(),
                static_cast<unsigned long long>(source.line_no),
                source.chr.c_str(),
                static_cast<long long>(source.pos),
                source.ref.c_str(),
                source_alt.c_str()
            );
        }
    }
    merged.insert(
        merged.end(),
        destination.alts.begin() + static_cast<std::ptrdiff_t>(destination_i),
        destination.alts.end()
    );
    merged.insert(
        merged.end(),
        source.alts.begin() + static_cast<std::ptrdiff_t>(source_i),
        source.alts.end()
    );
    destination.alts = std::move(merged);
}

static void prepare_extract_positions(ExtractSites& extract_sites, bcf_hdr_t* ghdr) {
    if (!extract_sites.active) return;

    for (ExtractPosition& position : extract_sites.positions) {
        position.geno_rid = bcf_hdr_name2id(ghdr, position.chr.c_str());
        if (position.geno_rid < 0) {
            position.geno_rid = extract_position_sort_rid(ghdr, position.chr);
        }
    }

    bool already_sorted = std::is_sorted(
        extract_sites.positions.begin(),
        extract_sites.positions.end(),
        extract_position_less
    );
    if (!already_sorted) {
        std::sort(
            extract_sites.positions.begin(),
            extract_sites.positions.end(),
            extract_position_less
        );
    }

    size_t write_i = 0;
    for (size_t read_i = 0; read_i < extract_sites.positions.size(); ++read_i) {
        ExtractPosition& source = extract_sites.positions[read_i];
        if (write_i == 0 ||
            extract_sites.positions[write_i - 1].geno_rid != source.geno_rid ||
            extract_sites.positions[write_i - 1].chr != source.chr ||
            extract_sites.positions[write_i - 1].pos != source.pos) {
            if (write_i != read_i) {
                extract_sites.positions[write_i] = std::move(source);
            }
            ++write_i;
            continue;
        }

        ExtractPosition& destination = extract_sites.positions[write_i - 1];
        if (destination.ref != source.ref) {
            die(
                "conflicting REF values in --extract list at %s:%llu for %s:%lld: %s vs %s",
                extract_sites.path.c_str(),
                static_cast<unsigned long long>(source.line_no),
                source.chr.c_str(),
                static_cast<long long>(source.pos),
                destination.ref.c_str(),
                source.ref.c_str()
            );
        }
        merge_extract_alts(destination, source, extract_sites.path);
    }
    extract_sites.positions.resize(write_i);

    std::fprintf(
        stderr,
        "Prepared --extract site list: %llu position(s), %llu split allele(s), %s input order; using linear merge.\n",
        static_cast<unsigned long long>(extract_sites.positions.size()),
        static_cast<unsigned long long>(extract_sites.allele_count),
        already_sorted ? "preserved" : "sorted"
    );
}

static int compare_extract_position_to_record(
    const ExtractPosition& position,
    int record_rid,
    int64_t record_pos
) {
    if (position.geno_rid != record_rid) {
        return position.geno_rid < record_rid ? -1 : 1;
    }
    if (position.pos != record_pos) {
        return position.pos < record_pos ? -1 : 1;
    }
    return 0;
}

static const ExtractPosition* match_next_extract_position(
    const std::vector<ExtractPosition>& positions,
    size_t& cursor,
    int record_rid,
    const char* record_chr,
    int64_t record_pos,
    const char* record_ref
) {
    while (cursor < positions.size() &&
           compare_extract_position_to_record(positions[cursor], record_rid, record_pos) < 0) {
        ++cursor;
    }

    if (cursor >= positions.size() ||
        compare_extract_position_to_record(positions[cursor], record_rid, record_pos) != 0) {
        return nullptr;
    }

    // Keep the cursor on an equal position until the VCF advances. This allows
    // multiple VCF records with the same CHROM/POS/REF to contribute different
    // ALT alleles from one multiallelic extract site.
    const ExtractPosition& position = positions[cursor];
    if (!record_ref || !*record_ref || std::strcmp(record_ref, ".") == 0) {
        die(
            "genotype VCF has unknown REF at --extract position %s:%lld",
            record_chr,
            static_cast<long long>(record_pos)
        );
    }
    if (position.ref != record_ref) {
        die(
            "REF mismatch for --extract site %s:%lld: list has %s, genotype VCF has %s",
            record_chr,
            static_cast<long long>(record_pos),
            position.ref.c_str(),
            record_ref
        );
    }
    return &position;
}

static bool fill_selected_alt_mask(
    bcf1_t* rec,
    const ExtractPosition* extract_position,
    std::vector<uint8_t>& selected
) {
    selected.assign(static_cast<size_t>(rec->n_allele), 0);
    if (!extract_position) return false;

    bool any_selected = false;
    for (int alt_idx = 1; alt_idx < rec->n_allele; ++alt_idx) {
        const char* alt = rec->d.allele[alt_idx];
        auto it = std::lower_bound(
            extract_position->alts.begin(),
            extract_position->alts.end(),
            alt,
            [](const std::string& value, const char* query) {
                return value.compare(query) < 0;
            }
        );
        bool matched = it != extract_position->alts.end() && *it == alt;
        selected[static_cast<size_t>(alt_idx)] = static_cast<uint8_t>(matched);
        any_selected = any_selected || matched;
    }
    return any_selected;
}

static bool fill_selected_alt_mask_fast(
    const FastGenotypeRecord& rec,
    const ExtractPosition* extract_position,
    std::vector<uint8_t>& selected
) {
    selected.assign(rec.alleles.size(), 0);
    if (!extract_position) return false;
    bool any_selected = false;
    for (size_t alt_idx = 1; alt_idx < rec.alleles.size(); ++alt_idx) {
        const char* alt = rec.alleles[alt_idx];
        auto it = std::lower_bound(
            extract_position->alts.begin(),
            extract_position->alts.end(),
            alt,
            [](const std::string& value, const char* query) {
                return value.compare(query) < 0;
            }
        );
        bool matched = it != extract_position->alts.end() && *it == alt;
        selected[alt_idx] = static_cast<uint8_t>(matched);
        any_selected = any_selected || matched;
    }
    return any_selected;
}

static std::string region_contig_token(const std::string& chr) {
    if (chr.find_first_of(":-") == std::string::npos) {
        return chr;
    }
    return "{" + chr + "}";
}

static void append_region_query(
    std::string& query,
    const std::string& chr,
    int64_t start,
    int64_t end
) {
    if (!query.empty()) query.push_back(',');
    query += region_contig_token(chr);
    query.push_back(':');
    query += std::to_string(start);
    query.push_back('-');
    query += std::to_string(end);
}

static std::string build_extract_genotype_region_query(
    const ExtractSites& extract_sites,
    const Region& region,
    uint64_t& n_positions,
    uint64_t& n_intervals
) {
    n_positions = 0;
    n_intervals = 0;

    std::string query;
    int64_t last_pos = -1;
    std::string last_chr;
    std::string interval_chr;
    int64_t interval_start = 0;
    int64_t interval_end = 0;
    bool have_interval = false;

    auto flush_interval = [&]() {
        if (!have_interval) return;
        append_region_query(query, interval_chr, interval_start, interval_end);
        ++n_intervals;
    };

    for (const ExtractPosition& position : extract_sites.positions) {
        if (!extract_position_in_region(position, region)) continue;
        if (position.chr == last_chr && position.pos == last_pos) continue;

        if (!have_interval) {
            interval_chr = position.chr;
            interval_start = position.pos;
            interval_end = position.pos;
            have_interval = true;
        } else if (position.chr == interval_chr) {
            interval_end = position.pos;
        } else {
            flush_interval();
            interval_chr = position.chr;
            interval_start = position.pos;
            interval_end = position.pos;
            have_interval = true;
        }

        last_chr = position.chr;
        last_pos = position.pos;
        ++n_positions;
    }
    flush_interval();

    return query;
}

static std::string build_extract_flare_region_query(
    const ExtractSites& extract_sites,
    const Region& region,
    uint64_t& n_contigs
) {
    n_contigs = 0;

    std::string query;
    std::string last_chr;
    for (const ExtractPosition& position : extract_sites.positions) {
        if (!extract_position_in_region(position, region)) continue;
        if (position.chr == last_chr) continue;
        append_region_query(query, position.chr, 1, kMaxVcfCoordinate);
        last_chr = position.chr;
        ++n_contigs;
    }

    return query;
}

static int infer_contigs_from_records(
    const char* vcf_path,
    bcf_hdr_t* dst,
    const char* dst_label
) {
    htsFile* fp = bcf_open(vcf_path, "r");
    if (!fp) die("cannot reopen VCF for contig inference: %s", vcf_path);

    bcf_hdr_t* scan_hdr = bcf_hdr_read(fp);
    if (!scan_hdr) die("cannot read VCF header for contig inference");

    bcf1_t* rec = bcf_init();
    if (!rec) die("failed to allocate VCF record for contig inference");

    std::unordered_set<std::string> seen;
    int added = 0;
    int ret = 0;
    while ((ret = read_next_record(fp, scan_hdr, rec)) == 0) {
        const char* chr = bcf_hdr_id2name(scan_hdr, rec->rid);
        if (!chr || !*chr) {
            die("cannot resolve chromosome during contig inference");
        }
        std::string contig(chr);
        if (seen.insert(contig).second &&
            add_contig_to_header_if_missing(dst, contig, dst_label)) {
            ++added;
        }
    }

    if (ret < -1) {
        die("error while scanning VCF for contig inference");
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(scan_hdr);
    bcf_close(fp);

    std::fprintf(
        stderr,
        "Inferred %d contig(s) for %s header by scanning %s.\n",
        added,
        dst_label,
        vcf_path
    );
    return added;
}

static bool ends_with(const std::string& value, const char* suffix) {
    size_t n = std::strlen(suffix);
    return value.size() >= n &&
           value.compare(value.size() - n, n, suffix) == 0;
}

static bool looks_indexable_variant_path(const char* path) {
    std::string value(path);
    bool indexed_suffix = ends_with(value, ".bcf") ||
                          ends_with(value, ".vcf.gz") ||
                          ends_with(value, ".vcf.bgz") ||
                          ends_with(value, ".bcf.gz");
    if (!indexed_suffix) return false;

    htsFile* probe = hts_open(path, "r");
    if (!probe) return false;
    const htsFormat* format = hts_get_format(probe);
    bool supported = format &&
                     (format->format == bcf ||
                      (format->format == vcf && format->compression == bgzf));
    hts_close(probe);
    return supported;
}

static bool init_synced_region_reader(
    bcf_srs_t** out,
    const char* path,
    const std::string& query,
    const char* label,
    const char* query_label = nullptr,
    const std::string* decode_samples = nullptr
) {
    const char* shown_query = query_label ? query_label : query.c_str();

    if (!looks_indexable_variant_path(path)) {
        std::fprintf(
            stderr,
            "WARNING: %s is not bgzip-compressed VCF/BCF; falling back to scan/filter.\n",
            label
        );
        return false;
    }

    bcf_srs_t* reader = bcf_sr_init();
    if (!reader) die("failed to initialize synced reader for %s", label);
    reader->max_unpack = BCF_UN_STR | BCF_UN_FMT;

    if (bcf_sr_set_regions(reader, query.c_str(), 0) != 0) {
        std::fprintf(
            stderr,
            "WARNING: cannot set %s region %s; falling back to scan/filter.\n",
            label,
            shown_query
        );
        bcf_sr_destroy(reader);
        return false;
    }

    if (decode_samples &&
        !bcf_sr_set_samples(reader, decode_samples->c_str(), 0)) {
        bcf_sr_destroy(reader);
        die("failed to configure indexed %s sample subset", label);
    }

    if (!bcf_sr_add_reader(reader, path)) {
        std::fprintf(
            stderr,
            "WARNING: cannot use indexed %s reader for %s (%s); falling back to scan/filter.\n",
            label,
            shown_query,
            bcf_sr_strerror(reader->errnum)
        );
        bcf_sr_destroy(reader);
        return false;
    }

    if (reader->nreaders > 0 && reader->readers[0].file) {
        hts_set_opt(reader->readers[0].file, HTS_OPT_BLOCK_SIZE, kInputBlockSize);
    }

    *out = reader;
    return true;
}

static int read_next_synced_record(bcf_srs_t* reader, bcf_hdr_t* hdr, bcf1_t*& rec) {
    while (bcf_sr_next_line(reader)) {
        bcf1_t* line = bcf_sr_get_line(reader, 0);
        if (!line) continue;
        if (hdr->keep_samples && line->n_sample != bcf_hdr_nsamples(hdr) &&
            bcf_subset_format(hdr, line) != 0) {
            die("failed to subset indexed VCF record samples");
        }
        bcf_unpack(line, BCF_UN_STR);
        rec = line;
        return 0;
    }
    return -1;
}

static bool record_in_region(bcf1_t* rec, const Region& region) {
    if (!region.active) return true;
    if (rec->rid != region.geno_rid) return false;
    int64_t pos = static_cast<int64_t>(rec->pos) + 1;
    return pos >= region.start && pos <= region.end;
}

static int read_next_genotype_record(
    htsFile* fp,
    bcf_hdr_t* hdr,
    bcf1_t*& rec,
    bcf1_t* storage,
    bcf_srs_t* synced_reader,
    const Region& region
) {
    if (synced_reader) {
        int ret = 0;
        while ((ret = read_next_synced_record(synced_reader, hdr, rec)) == 0) {
            if (record_in_region(rec, region)) return 0;
        }
        return ret;
    }

    rec = storage;
    int ret = 0;
    while ((ret = read_next_record(fp, hdr, rec)) == 0) {
        if (record_in_region(rec, region)) return 0;
    }
    return ret;
}

static bool read_next_lai_record(
    htsFile* fp,
    bcf_hdr_t* ahdr,
    bcf_hdr_t* ghdr,
    bcf1_t*& rec,
    bcf1_t* storage,
    bcf_srs_t* synced_reader,
    int flare_raw_n_samples,
    const std::vector<int>& flare_raw_indices,
    int n_ancestries,
    FlareDeltaDecoder& decoder,
    std::unordered_set<std::string>& warned_missing_flare_contigs,
    LaiRecord& out
) {
    if (!synced_reader) rec = storage;
    while ((synced_reader ? read_next_synced_record(synced_reader, ahdr, rec)
                          : read_next_record(fp, ahdr, rec)) == 0) {
        const char* a_chr = bcf_hdr_id2name(ahdr, rec->rid);
        int a_geno_rid = bcf_hdr_name2id(ghdr, a_chr);

        if (a_geno_rid < 0) {
            std::string contig(a_chr);
            if (warned_missing_flare_contigs.insert(contig).second) {
                std::fprintf(
                    stderr,
                    "WARNING: skipping FLARE contig absent from genotype header: %s\n",
                    a_chr
                );
            }
            continue;
        }

        out.valid = true;
        out.merged_duplicate = false;
        out.geno_rid = a_geno_rid;
        out.chr = a_chr;
        out.pos = static_cast<int64_t>(rec->pos) + 1;

        decode_flare_changes(
            ahdr,
            rec,
            flare_raw_n_samples,
            flare_raw_indices,
            n_ancestries,
            a_geno_rid,
            decoder,
            out
        );

        return true;
    }

    out = LaiRecord{};
    return false;
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s genotype.phased.vcf.gz flare.anc.vcf.gz n_ancestries rare_threshold|auto out_prefix [chr:start-end] [--keep samples.txt] [--extract sites.pvar|sites.vcf]\n\n"
        "Example:\n"
        "  %s chr1.phased.vcf.gz chr1.flare.anc.vcf.gz 3 auto chr1\n"
        "  %s chr22.phased.vcf.gz chr22.flare.anc.vcf.gz 5 512 chr22.1 chr22:1-50000000\n",
        prog,
        prog,
        prog
    );
}

static int parse_rare_threshold_arg(const char* text) {
    if (std::strcmp(text, "auto") == 0 || std::strcmp(text, "default") == 0) {
        return -1;
    }

    char* end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0') {
        die("rare_threshold must be a non-negative integer or auto: %s", text);
    }
    if (value < 0 || value > std::numeric_limits<int>::max()) {
        die("rare_threshold out of range: %s", text);
    }
    return static_cast<int>(value);
}

static int default_rare_threshold_from_samples(int n_samples) {
    return (n_samples + 31) / 32;
}

int main(int argc, char** argv) {
    if (argc < 6) {
        print_usage(argv[0]);
        return 1;
    }

    const char* geno_vcf = argv[1];
    const char* flare_vcf = argv[2];
    int n_ancestries = std::atoi(argv[3]);
    int rare_threshold = parse_rare_threshold_arg(argv[4]);
    const char* out_prefix = argv[5];
    Region region;
    std::string keep_path;
    std::string extract_path;
    for (int argi = 6; argi < argc; ++argi) {
        std::string arg = argv[argi];
        if (arg == "--keep") {
            if (argi + 1 >= argc) die("--keep requires a value");
            keep_path = argv[++argi];
        } else if (arg == "--extract") {
            if (argi + 1 >= argc) die("--extract requires a value");
            extract_path = argv[++argi];
        } else if (!arg.empty() && arg[0] == '-') {
            die("unknown option: %s", arg.c_str());
        } else {
            if (region.active) die("multiple region arguments supplied");
            region = parse_region_string(arg);
        }
    }

    if (n_ancestries <= 0 || n_ancestries > 32) {
        die("n_ancestries must be in [1, 32]");
    }

    htsFile* gfp = bcf_open(geno_vcf, "r");
    htsFile* afp = bcf_open(flare_vcf, "r");

    if (!gfp || !afp) {
        die("cannot open input VCF/BCF");
    }

    hts_set_opt(gfp, HTS_OPT_BLOCK_SIZE, kInputBlockSize);
    hts_set_opt(afp, HTS_OPT_BLOCK_SIZE, kInputBlockSize);

    bcf_hdr_t* ghdr = bcf_hdr_read(gfp);
    bcf_hdr_t* ahdr = bcf_hdr_read(afp);

    if (!ghdr || !ahdr) {
        die("cannot read input headers");
    }

    const htsFormat* genotype_format = hts_get_format(gfp);
    const htsFormat* flare_format = hts_get_format(afp);
    bool fast_genotype_text = genotype_format && genotype_format->format == vcf;
    bool fast_flare_text = flare_format && flare_format->format == vcf;

    SampleSelection sample_selection = build_sample_selection(ghdr, ahdr, keep_path);
    int genotype_text_n_samples = sample_selection.genotype_raw_sample_count;
    int flare_text_n_samples = sample_selection.flare_raw_sample_count;
    sample_selection.genotype_text_indices = sample_selection.genotype_raw_indices;
    sample_selection.flare_text_indices = sample_selection.flare_raw_indices;
    apply_decode_sample_subsets(ghdr, ahdr, sample_selection);
    int genotype_raw_n_samples = sample_selection.genotype_raw_sample_count;
    int flare_raw_n_samples = sample_selection.flare_raw_sample_count;
    int n_samples = static_cast<int>(sample_selection.sample_ids.size());
    if (n_samples <= 0) {
        die("retained zero samples");
    }

    ExtractSites extract_sites = load_extract_sites(extract_path);
    if (extract_sites.active) {
        std::fprintf(
            stderr,
            "Loaded --extract site list: retaining up to %llu split allele(s).\n",
            static_cast<unsigned long long>(extract_sites.allele_count)
        );
    }

    if (ghdr->n[BCF_DT_CTG] == 0) {
        std::fprintf(
            stderr,
            "WARNING: genotype VCF header has no ##contig lines; adding FLARE/header region contigs to the in-memory genotype header.\n"
        );
        add_contigs_from_header_if_missing(ghdr, ahdr, "genotype", "FLARE");
        add_contigs_from_extract_if_missing(ghdr, extract_sites, "genotype");
        if (!region.active && ghdr->n[BCF_DT_CTG] == 0) {
            infer_contigs_from_records(geno_vcf, ghdr, "genotype");
        }
    }
    add_contigs_from_extract_if_missing(ghdr, extract_sites, "genotype");

    bcf_srs_t* genotype_region_reader = nullptr;
    bcf_srs_t* flare_region_reader = nullptr;
    bool skip_genotype_loop = false;
    bool bounded_extract_reader = false;
    std::string progress_scope;
    std::string genotype_indexed_query;
    std::string flare_indexed_query;

    if (region.active) {
        add_contig_to_header_if_missing(ghdr, region.chr, "genotype");
        region.geno_rid = bcf_hdr_name2id(ghdr, region.chr.c_str());
        if (region.geno_rid < 0) {
            die("failed to add region chromosome to genotype header: %s", region.chr.c_str());
        }
    }
    prepare_extract_positions(extract_sites, ghdr);

    if (extract_sites.active) {
        uint64_t n_extract_positions = 0;
        uint64_t n_extract_intervals = 0;
        std::string genotype_query = build_extract_genotype_region_query(
            extract_sites,
            region,
            n_extract_positions,
            n_extract_intervals
        );

        if (n_extract_positions == 0) {
            skip_genotype_loop = true;
            progress_scope = "--extract bounded intervals";
            std::fprintf(
                stderr,
                "No --extract positions overlap the selected conversion scope; skipping genotype scan.\n"
            );
        } else {
            bool have_genotype_extract_reader = init_synced_region_reader(
                &genotype_region_reader,
                geno_vcf,
                genotype_query,
                "genotype VCF",
                "--extract bounded intervals",
                sample_selection.genotype_subset_active
                    ? &sample_selection.genotype_decode_samples
                    : nullptr
            );

            if (have_genotype_extract_reader) {
                genotype_indexed_query = genotype_query;
                bounded_extract_reader = true;
                progress_scope = "--extract bounded intervals";
                std::fprintf(
                    stderr,
                    "Using indexed --extract bounded genotype reader: %llu target position(s) across %llu interval(s), %llu split allele(s).\n",
                    static_cast<unsigned long long>(n_extract_positions),
                    static_cast<unsigned long long>(n_extract_intervals),
                    static_cast<unsigned long long>(extract_sites.allele_count)
                );
            } else {
                std::fprintf(
                    stderr,
                    "WARNING: --extract could not use genotype VCF bounded random access; scanning records and filtering by CHROM/POS/REF/ALT.\n"
                );
            }
        }
    }

    if (!extract_sites.active && region.active) {
        std::string flare_query =
            region.chr + ":1-" + std::to_string(kMaxVcfCoordinate);
        bool have_genotype_region_reader = init_synced_region_reader(
            &genotype_region_reader,
            geno_vcf,
            region.label,
            "genotype VCF",
            nullptr,
            sample_selection.genotype_subset_active
                ? &sample_selection.genotype_decode_samples
                : nullptr
        );
        if (have_genotype_region_reader) {
            progress_scope = region.label;
            genotype_indexed_query = region.label;
        }
        bool have_flare_region_reader = init_synced_region_reader(
            &flare_region_reader,
            flare_vcf,
            flare_query,
            "FLARE VCF",
            nullptr,
            sample_selection.flare_subset_active
                ? &sample_selection.flare_decode_samples
                : nullptr
        );

        if (!have_genotype_region_reader || !have_flare_region_reader) {
            if (genotype_region_reader) {
                bcf_sr_destroy(genotype_region_reader);
                genotype_region_reader = nullptr;
            }
            if (flare_region_reader) {
                bcf_sr_destroy(flare_region_reader);
                flare_region_reader = nullptr;
            }
            genotype_indexed_query.clear();
            flare_indexed_query.clear();
            std::fprintf(
                stderr,
                "WARNING: indexed region reader unavailable; scanning inputs and filtering %s.\n",
                region.label.c_str()
            );
        } else {
            flare_indexed_query = flare_query;
        }
    } else if (region.active) {
        std::string flare_query =
            region.chr + ":1-" + std::to_string(kMaxVcfCoordinate);
        bool have_flare_region_reader = init_synced_region_reader(
            &flare_region_reader,
            flare_vcf,
            flare_query,
            "FLARE VCF",
            nullptr,
            sample_selection.flare_subset_active
                ? &sample_selection.flare_decode_samples
                : nullptr
        );
        if (!have_flare_region_reader) {
            std::fprintf(
                stderr,
                "WARNING: indexed FLARE reader unavailable; scanning FLARE records while using --extract filter.\n"
            );
        } else {
            flare_indexed_query = flare_query;
        }
    } else if (extract_sites.active && !skip_genotype_loop) {
        uint64_t n_extract_contigs = 0;
        std::string flare_query = build_extract_flare_region_query(
            extract_sites,
            region,
            n_extract_contigs
        );
        bool have_flare_extract_reader = init_synced_region_reader(
            &flare_region_reader,
            flare_vcf,
            flare_query,
            "FLARE VCF",
            "--extract contig list",
            sample_selection.flare_subset_active
                ? &sample_selection.flare_decode_samples
                : nullptr
        );
        if (have_flare_extract_reader) {
            flare_indexed_query = flare_query;
            std::fprintf(
                stderr,
                "Using indexed FLARE reader for --extract: %llu contig(s).\n",
                static_cast<unsigned long long>(n_extract_contigs)
            );
        } else {
            std::fprintf(
                stderr,
                "WARNING: indexed FLARE reader unavailable for --extract contigs; scanning FLARE records.\n"
            );
        }
    }

    FastTextVcfReader fast_genotype_reader;
    FastTextVcfReader fast_flare_reader;
    std::vector<int> flare_text_raw_to_output;
    if (fast_genotype_text && !skip_genotype_loop) {
        int previous_index = -1;
        for (int raw_index : sample_selection.genotype_text_indices) {
            if (raw_index <= previous_index || raw_index >= genotype_text_n_samples) {
                die("internal genotype text sample indices are not sorted and in range");
            }
            previous_index = raw_index;
        }
        if (genotype_region_reader) {
            bcf_sr_destroy(genotype_region_reader);
            genotype_region_reader = nullptr;
        }
        init_fast_text_reader(
            fast_genotype_reader,
            gfp,
            geno_vcf,
            genotype_indexed_query
        );
    }
    if (fast_flare_text && !skip_genotype_loop) {
        if (flare_region_reader) {
            bcf_sr_destroy(flare_region_reader);
            flare_region_reader = nullptr;
        }
        init_fast_text_reader(
            fast_flare_reader,
            afp,
            flare_vcf,
            flare_indexed_query
        );
        flare_text_raw_to_output.assign(
            static_cast<size_t>(flare_text_n_samples),
            -1
        );
        for (size_t out_i = 0; out_i < sample_selection.flare_text_indices.size(); ++out_i) {
            int raw_i = sample_selection.flare_text_indices[out_i];
            if (raw_i < 0 || raw_i >= flare_text_n_samples) {
                die("internal FLARE text sample index out of range");
            }
            flare_text_raw_to_output[static_cast<size_t>(raw_i)] = static_cast<int>(out_i);
        }
    }
    if (fast_genotype_text || fast_flare_text) {
        std::fprintf(
            stderr,
            "Using FELIXla direct text VCF parser for%s%s input.\n",
            fast_genotype_text ? " genotype" : "",
            fast_flare_text ? " FLARE" : ""
        );
    }

    ProgressReporter progress(
        "genotype VCF",
        geno_vcf,
        "split variants",
        progress_scope.empty() ? (region.active ? region.label.c_str() : nullptr) : progress_scope.c_str()
    );

    if (rare_threshold < 0) {
        rare_threshold = default_rare_threshold_from_samples(n_samples);
        std::fprintf(
            stderr,
            "Using rare_threshold=%d = ceil(%d / 32) from retained sample count.\n",
            rare_threshold,
            n_samples
        );
    }

    uint64_t n_haps = static_cast<uint64_t>(n_samples) * 2ULL;
    if (n_haps > static_cast<uint64_t>(kMaxPackedHapId) + 1ULL) {
        die("hap_id exceeds 27-bit packed limit");
    }

    int n_words = static_cast<int>((n_haps + 63ULL) / 64ULL);

    std::string common_bin = std::string(out_prefix) + ".common.geno.bin";
    std::string common_mks = std::string(out_prefix) + ".common.variant.mks";
    std::string common_idx = std::string(out_prefix) + ".common.variant.idx";
    std::string rare_bin = std::string(out_prefix) + ".rare.carrier.bin";
    std::string rare_mks = std::string(out_prefix) + ".rare.variant.mks";
    std::string rare_idx = std::string(out_prefix) + ".rare.variant.idx";
    std::string anc_bin = std::string(out_prefix) + ".ancblock.bin";
    std::string anc_mks = std::string(out_prefix) + ".ancblock.mks";
    std::string anc_idx = std::string(out_prefix) + ".ancblock.idx";
    std::string samples_path = std::string(out_prefix) + ".samples";
    std::string meta_path = std::string(out_prefix) + ".meta";

    FILE* common_fp = open_output_or_die(common_bin, "wb");
    FILE* common_mks_fp = open_output_or_die(common_mks, "wb");
    FILE* common_idx_fp = open_output_or_die(common_idx, "wb");
    FILE* rare_fp = open_output_or_die(rare_bin, "wb");
    FILE* rare_mks_fp = open_output_or_die(rare_mks, "wb");
    FILE* rare_idx_fp = open_output_or_die(rare_idx, "wb");
    FILE* anc_fp = open_output_or_die(anc_bin, "wb");
    FILE* anc_mks_fp = open_output_or_die(anc_mks, "wb");
    FILE* anc_idx_fp = open_output_or_die(anc_idx, "wb");

    const char common_mks_magic[8] = {'T', 'R', 'C', 'M', 'M', 'K', 'S', '1'};
    const char rare_mks_magic[8] = {'T', 'R', 'R', 'A', 'M', 'K', 'S', '1'};
    const char anc_mks_magic[8] = {'T', 'R', 'A', 'N', 'M', 'K', 'S', '1'};
    const char common_idx_magic[8] = {'T', 'R', 'C', 'M', 'I', 'D', 'X', '2'};
    const char rare_idx_magic[8] = {'T', 'R', 'R', 'A', 'I', 'D', 'X', '2'};
    const char anc_idx_magic[8] = {'T', 'R', 'A', 'N', 'I', 'D', 'X', '2'};
    write_magic_header(common_mks_fp, common_mks_magic);
    write_magic_header(rare_mks_fp, rare_mks_magic);
    write_magic_header(anc_mks_fp, anc_mks_magic);
    write_magic_header(common_idx_fp, common_idx_magic);
    write_magic_header(rare_idx_fp, rare_idx_magic);
    write_magic_header(anc_idx_fp, anc_idx_magic);

    write_sidecars(
        samples_path,
        meta_path,
        geno_vcf,
        flare_vcf,
        sample_selection.sample_ids,
        n_haps,
        n_words,
        n_ancestries,
        rare_threshold,
        region.active ? region.label.c_str() : nullptr,
        keep_path.empty() ? nullptr : keep_path.c_str(),
        extract_path.empty() ? nullptr : extract_path.c_str()
    );

    bcf1_t* grec_storage = bcf_init();
    bcf1_t* arec_storage = bcf_init();
    bcf1_t* grec = grec_storage;
    bcf1_t* arec = arec_storage;
    FastGenotypeRecord fast_grec;

    OpenAncestryBlock block;
    std::vector<AltCarrierBuilder> builders_by_alt;
    std::vector<uint8_t> selected_alts;
    std::vector<RareCarrierPacked> rare_carrier_batch;
    std::unordered_set<std::string> warned_missing_flare_contigs;
    const std::vector<ExtractPosition>& ordered_extract_positions = extract_sites.positions;
    size_t extract_cursor = 0;

    uint32_t n_blocks_written = 0;
    uint64_t common_index = 0;
    uint64_t rare_index = 0;
    uint32_t global_variant_index = 0;

    FlareDeltaDecoder flare_decoder;
    LaiRecord lai_record;
    auto next_lai_record = [&]() -> bool {
        if (fast_flare_text) {
            return read_next_fast_lai_record(
                fast_flare_reader,
                ghdr,
                flare_text_raw_to_output,
                n_samples,
                n_ancestries,
                flare_decoder,
                warned_missing_flare_contigs,
                lai_record
            );
        }
        return read_next_lai_record(
            afp,
            ahdr,
            ghdr,
            arec,
            arec_storage,
            flare_region_reader,
            flare_raw_n_samples,
            sample_selection.flare_raw_indices,
            n_ancestries,
            flare_decoder,
            warned_missing_flare_contigs,
            lai_record
        );
    };
    auto next_genotype_record = [&]() -> int {
        if (fast_genotype_text) {
            return read_next_fast_genotype_record(
                fast_genotype_reader,
                ghdr,
                fast_grec,
                region
            );
        }
        return read_next_genotype_record(
            gfp,
            ghdr,
            grec,
            grec_storage,
            genotype_region_reader,
            region
        );
    };
    bool has_lai_record = false;
    if (!skip_genotype_loop) {
        has_lai_record = next_lai_record();
    }
    int last_flare_geno_rid = -1;
    int64_t last_flare_pos = 0;

    int ancestry_state_rid = -1;
    std::string ancestry_state_chr;
    bool ancestry_state_ready = false;

    while (!skip_genotype_loop && next_genotype_record() == 0) {
        progress.record_scanned();

        int g_rid = fast_genotype_text ? fast_grec.rid : grec->rid;
        const char* g_chr = fast_genotype_text
            ? fast_grec.chr
            : bcf_hdr_id2name(ghdr, grec->rid);
        int64_t g_pos = fast_genotype_text
            ? fast_grec.pos
            : static_cast<int64_t>(grec->pos) + 1;
        int g_n_allele = fast_genotype_text
            ? static_cast<int>(fast_grec.alleles.size())
            : grec->n_allele;
        const char* g_ref = fast_genotype_text
            ? fast_grec.ref
            : (grec->n_allele > 0 ? grec->d.allele[0] : nullptr);
        const char* g_raw_id = fast_genotype_text ? fast_grec.id : grec->d.id;
        auto g_allele = [&](int allele_index) -> const char* {
            return fast_genotype_text
                ? fast_grec.alleles[static_cast<size_t>(allele_index)]
                : grec->d.allele[allele_index];
        };

        while (has_lai_record) {
            bool current_block_covers =
                block.active &&
                block.geno_rid == g_rid &&
                g_pos <= block.end_pos;

            if (lai_record.geno_rid > g_rid ||
                (lai_record.geno_rid == g_rid && current_block_covers)) {
                break;
            }

            LaiRecord interval_record = std::move(lai_record);
            has_lai_record = next_lai_record();

            while (has_lai_record &&
                   lai_record.geno_rid == interval_record.geno_rid &&
                   lai_record.pos == interval_record.pos) {
                interval_record.reset_state =
                    interval_record.reset_state || lai_record.reset_state;
                interval_record.merged_duplicate = true;
                interval_record.changes.insert(
                    interval_record.changes.end(),
                    lai_record.changes.begin(),
                    lai_record.changes.end()
                );
                has_lai_record = next_lai_record();
            }

            if (last_flare_geno_rid > interval_record.geno_rid) {
                die("FLARE records are not sorted by genotype header contig order");
            }

            int64_t interval_start = 1;
            if (last_flare_geno_rid == interval_record.geno_rid) {
                if (interval_record.pos <= last_flare_pos) {
                    die(
                        "FLARE records are not strictly increasing at %s:%lld",
                        interval_record.chr.c_str(),
                        static_cast<long long>(interval_record.pos)
                    );
                }
                interval_start = last_flare_pos + 1;
            }

            int64_t interval_end = interval_record.pos;
            bool state_changed = ancestry_record_changes_state(
                block.state,
                interval_record,
                ancestry_state_rid
            );
            if (block.active &&
                (state_changed || block.geno_rid != interval_record.geno_rid)) {
                close_open_block(
                    anc_fp,
                    anc_bin.c_str(),
                    anc_mks_fp,
                    anc_mks.c_str(),
                    anc_idx_fp,
                    block,
                    n_blocks_written,
                    n_ancestries,
                    n_words
                );
            }
            apply_ancestry_changes(
                block.state,
                interval_record,
                n_ancestries,
                n_words,
                static_cast<size_t>(n_haps)
            );
            ancestry_state_rid = interval_record.geno_rid;
            ancestry_state_chr = interval_record.chr;
            ancestry_state_ready = true;

            int64_t output_start = interval_start;
            int64_t output_end = interval_end;
            bool emit_interval = true;

            if (region.active) {
                if (interval_record.geno_rid != region.geno_rid) {
                    emit_interval = false;
                } else if (interval_end < region.start) {
                    emit_interval = false;
                } else if (interval_start > region.end) {
                    emit_interval = false;
                } else {
                    output_start = std::max(interval_start, region.start);
                    output_end = std::min(interval_end, region.end);
                    emit_interval = output_start <= output_end;
                }
            }

            if (!region.active && extract_sites.active && bounded_extract_reader &&
                extract_sites.contigs.count(interval_record.chr) == 0) {
                emit_interval = false;
            }

            if (emit_interval) {
                update_open_block_from_flare(
                    anc_fp,
                    anc_bin.c_str(),
                    anc_mks_fp,
                    anc_mks.c_str(),
                    anc_idx_fp,
                    block,
                    n_blocks_written,
                    n_ancestries,
                    n_words,
                    interval_record.geno_rid,
                    interval_record.chr.c_str(),
                    output_start,
                    output_end,
                    state_changed
                );
            }

            last_flare_geno_rid = interval_record.geno_rid;
            last_flare_pos = interval_record.pos;
        }

        if (region.active && !block.active && !has_lai_record &&
            ancestry_state_ready && ancestry_state_rid == g_rid) {
            update_open_block_from_flare(
                anc_fp,
                anc_bin.c_str(),
                anc_mks_fp,
                anc_mks.c_str(),
                anc_idx_fp,
                block,
                n_blocks_written,
                n_ancestries,
                n_words,
                ancestry_state_rid,
                ancestry_state_chr.c_str(),
                region.start,
                region.end,
                false
            );
        }

        if (block.active && block.geno_rid == g_rid && g_pos > block.end_pos &&
            (!has_lai_record || lai_record.geno_rid > g_rid)) {
            block.end_pos = g_pos;
        }

        const char* ref = g_ref;
        const ExtractPosition* extract_position = nullptr;
        if (extract_sites.active) {
            extract_position = match_next_extract_position(
                ordered_extract_positions,
                extract_cursor,
                g_rid,
                g_chr,
                g_pos,
                ref
            );
            if (!extract_position) {
                progress.maybe_report(global_variant_index, common_index, rare_index);
                continue;
            }
        }

        if (!block.active || block.geno_rid != g_rid || g_pos > block.end_pos) {
            progress.maybe_report(global_variant_index, common_index, rare_index);
            continue;
        }

        if (g_n_allele < 2) {
            progress.maybe_report(global_variant_index, common_index, rare_index);
            continue;
        }

        const char* raw_id = g_raw_id;
        if (extract_sites.active) {
            bool selected_any = fast_genotype_text
                ? fill_selected_alt_mask_fast(fast_grec, extract_position, selected_alts)
                : fill_selected_alt_mask(grec, extract_position, selected_alts);
            if (!selected_any) {
                progress.maybe_report(global_variant_index, common_index, rare_index);
                continue;
            }
        } else {
            selected_alts.assign(static_cast<size_t>(g_n_allele), 1);
            selected_alts[0] = 0;
        }

        if (fast_genotype_text) {
            process_fast_text_genotypes(
                fast_grec,
                block,
                genotype_text_n_samples,
                sample_selection.genotype_text_indices,
                selected_alts,
                rare_threshold,
                n_words,
                builders_by_alt
            );
        } else {
            process_genotypes(
                ghdr,
                grec,
                block,
                genotype_raw_n_samples,
                sample_selection.genotype_raw_indices,
                selected_alts,
                rare_threshold,
                n_words,
                builders_by_alt
            );
        }

        if (builders_by_alt.empty()) {
            progress.maybe_report(global_variant_index, common_index, rare_index);
            continue;
        }

        for (int alt_idx = 1; alt_idx < g_n_allele; ++alt_idx) {
            if (!selected_alts[static_cast<size_t>(alt_idx)]) continue;
            const char* alt = g_allele(alt_idx);

            if (global_variant_index == std::numeric_limits<uint32_t>::max()) {
                die("global split variant index exceeds uint32_t limit");
            }

            const AltCarrierBuilder& builder = builders_by_alt[static_cast<size_t>(alt_idx)];
            uint32_t mac = builder.mac;
            std::string split_id = make_split_id(raw_id, g_chr, g_pos, ref, alt);

            if (mac <= static_cast<uint32_t>(rare_threshold)) {
                uint64_t carrier_offset = tell_or_die(rare_fp, rare_bin.c_str());

                rare_carrier_batch.clear();
                rare_carrier_batch.reserve(builder.sparse_haps.size());
                for (uint32_t hap_id : builder.sparse_haps) {
                    int8_t ancestry = block.state.hap_ancestry[hap_id];
                    rare_carrier_batch.push_back(RareCarrierPacked{
                        global_variant_index,
                        pack_anc_hap(static_cast<uint32_t>(ancestry), hap_id)
                    });
                }
                write_exact(
                    rare_fp,
                    rare_carrier_batch.data(),
                    rare_carrier_batch.size() * sizeof(RareCarrierPacked),
                    "rare carriers"
                );

                uint64_t rare_mks_offset = tell_or_die(rare_mks_fp, rare_mks.c_str());
                write_rare_mks_record(
                    rare_mks_fp,
                    rare_index,
                    global_variant_index,
                    g_chr,
                    g_pos,
                    split_id,
                    ref,
                    alt,
                    static_cast<uint32_t>(alt_idx),
                    carrier_offset,
                    mac,
                    mac
                );
                write_rare_offset_idx_record(
                    rare_idx_fp,
                    rare_index,
                    global_variant_index,
                    rare_mks_offset,
                    carrier_offset,
                    mac
                );

                ++rare_index;
            } else {
                uint64_t geno_offset = tell_or_die(common_fp, common_bin.c_str());
                size_t wrote = std::fwrite(
                    builder.dense_bits.data(),
                    sizeof(uint64_t),
                    static_cast<size_t>(n_words),
                    common_fp
                );

                if (wrote != static_cast<size_t>(n_words)) {
                    die("failed writing common genotype for split variant %u", global_variant_index);
                }

                uint64_t common_mks_offset = tell_or_die(common_mks_fp, common_mks.c_str());
                write_common_mks_record(
                    common_mks_fp,
                    common_index,
                    global_variant_index,
                    g_chr,
                    g_pos,
                    split_id,
                    ref,
                    alt,
                    static_cast<uint32_t>(alt_idx),
                    block.block_id,
                    geno_offset,
                    mac
                );
                write_common_offset_idx_record(
                    common_idx_fp,
                    common_index,
                    global_variant_index,
                    common_mks_offset,
                    geno_offset
                );

                ++common_index;
            }

            ++global_variant_index;
        }

        progress.maybe_report(global_variant_index, common_index, rare_index);
    }

    close_open_block(
        anc_fp,
        anc_bin.c_str(),
        anc_mks_fp,
        anc_mks.c_str(),
        anc_idx_fp,
        block,
        n_blocks_written,
        n_ancestries,
        n_words
    );

    bcf_destroy(grec_storage);
    bcf_destroy(arec_storage);
    if (fast_genotype_text) destroy_fast_text_reader(fast_genotype_reader);
    if (fast_flare_text) destroy_fast_text_reader(fast_flare_reader);
    if (genotype_region_reader) bcf_sr_destroy(genotype_region_reader);
    if (flare_region_reader) bcf_sr_destroy(flare_region_reader);
    bcf_hdr_destroy(ghdr);
    bcf_hdr_destroy(ahdr);
    bcf_close(gfp);
    bcf_close(afp);

    std::fclose(common_fp);
    std::fclose(common_mks_fp);
    std::fclose(common_idx_fp);
    std::fclose(rare_fp);
    std::fclose(rare_mks_fp);
    std::fclose(rare_idx_fp);
    std::fclose(anc_fp);
    std::fclose(anc_mks_fp);
    std::fclose(anc_idx_fp);

    progress.finish(global_variant_index, common_index, rare_index);

    std::fprintf(stderr, "Finished.\n");
    std::fprintf(stderr, "Global split variants: %u\n", global_variant_index);
    std::fprintf(stderr, "Common variants:       %llu\n", static_cast<unsigned long long>(common_index));
    std::fprintf(stderr, "Rare variants:         %llu\n", static_cast<unsigned long long>(rare_index));
    std::fprintf(stderr, "Ancestry blocks:       %u\n", n_blocks_written);

    return 0;
}
