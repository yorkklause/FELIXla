// estimate_mac_threshold.cpp
// designed by Kai, implemented by codex
//
// Estimate a sparse/dense MAC threshold for the ancestry-aware packed backend.

#include <htslib/hts.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <chrono>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <unistd.h>

struct Options {
    const char* genotype_path = nullptr;
    double storage_weight = 1.0;
    double query_weight = 0.0;
    double dense_word_cost = 1.0;
    double sparse_carrier_cost = 1.0;
    uint64_t max_records = 0;
};

struct ScanStats {
    int n_samples = 0;
    uint64_t n_haps = 0;
    uint64_t n_words = 0;
    uint64_t records_seen = 0;
    uint64_t split_variants = 0;
    uint64_t total_alt_carriers = 0;
    uint32_t max_mac = 0;
    std::map<uint32_t, uint64_t> mac_counts;
};

struct ThresholdScore {
    uint32_t threshold = 0;
    uint64_t rare_variants = 0;
    uint64_t common_variants = 0;
    uint64_t rare_carriers = 0;
    uint64_t common_dense_words = 0;
    long double payload_bytes = 0.0L;
    long double query_units = 0.0L;
    long double score = 0.0L;
};

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

static constexpr uint64_t kProgressRecordInterval = 10000;
static constexpr int kProgressSecondsInterval = 2;

static hts_idx_t* load_any_vcf_index(const char* path) {
    hts_idx_t* idx = hts_idx_load3(path, nullptr, HTS_FMT_CSI, HTS_IDX_SILENT_FAIL);
    if (idx) return idx;
    return hts_idx_load3(path, nullptr, HTS_FMT_TBI, HTS_IDX_SILENT_FAIL);
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
    ProgressReporter(const char* input_path, uint64_t max_records)
        : total_(get_index_record_count(input_path)),
          start_time_(std::chrono::steady_clock::now()),
          last_report_(start_time_),
          stderr_is_tty_(isatty(fileno(stderr)) != 0) {
        if (total_.available && max_records != 0 && max_records < total_.total) {
            total_.total = max_records;
        }

        if (total_.available) {
            std::fprintf(
                stderr,
                "Progress: genotype VCF index reports %llu records%s.\n",
                static_cast<unsigned long long>(total_.total),
                max_records != 0 ? " after --max-records cap" : ""
            );
        } else {
            std::fprintf(
                stderr,
                "Progress: genotype VCF index record count unavailable; reporting scanned records only.\n"
            );
        }
    }

    void maybe_report(const ScanStats& stats, bool force = false) {
        auto now = std::chrono::steady_clock::now();
        bool count_due = stats.records_seen >= next_record_report_;
        bool time_due = now - last_report_ >= std::chrono::seconds(kProgressSecondsInterval);

        if (!force && !count_due && !time_due) {
            return;
        }

        print(stats);
        last_report_ = now;
        printed_ = true;

        if (stats.records_seen >= next_record_report_) {
            if (stats.records_seen > std::numeric_limits<uint64_t>::max() - kProgressRecordInterval) {
                next_record_report_ = std::numeric_limits<uint64_t>::max();
            } else {
                next_record_report_ = stats.records_seen + kProgressRecordInterval;
            }
        }
    }

    void finish(const ScanStats& stats) {
        maybe_report(stats, true);
        if (stderr_is_tty_ && printed_) {
            std::fputc('\n', stderr);
        }
    }

private:
    static std::string format_duration(double seconds, bool round_up) {
        if (seconds < 0.0) seconds = 0.0;

        uint64_t rounded = static_cast<uint64_t>(seconds);
        if (round_up && seconds > static_cast<double>(rounded)) {
            ++rounded;
        }

        uint64_t hours = rounded / 3600ULL;
        uint64_t minutes = (rounded % 3600ULL) / 60ULL;
        uint64_t secs = rounded % 60ULL;

        char buf[64];
        if (hours > 0) {
            std::snprintf(
                buf,
                sizeof(buf),
                "%lluh%02llum%02llus",
                static_cast<unsigned long long>(hours),
                static_cast<unsigned long long>(minutes),
                static_cast<unsigned long long>(secs)
            );
        } else if (minutes > 0) {
            std::snprintf(
                buf,
                sizeof(buf),
                "%llum%02llus",
                static_cast<unsigned long long>(minutes),
                static_cast<unsigned long long>(secs)
            );
        } else {
            std::snprintf(
                buf,
                sizeof(buf),
                "%llus",
                static_cast<unsigned long long>(secs)
            );
        }

        return buf;
    }

    static std::string format_rate(double records_per_second) {
        char buf[64];
        if (records_per_second >= 1000000.0) {
            std::snprintf(buf, sizeof(buf), "%.2fM records/s", records_per_second / 1000000.0);
        } else if (records_per_second >= 1000.0) {
            std::snprintf(buf, sizeof(buf), "%.2fk records/s", records_per_second / 1000.0);
        } else {
            std::snprintf(buf, sizeof(buf), "%.1f records/s", records_per_second);
        }

        return buf;
    }

    void print(const ScanStats& stats) const {
        const char* prefix = stderr_is_tty_ ? "\r" : "";
        const char* suffix = stderr_is_tty_ ? "" : "\n";
        auto now = std::chrono::steady_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(now - start_time_).count();
        double records_per_second =
            elapsed_seconds > 0.0
                ? static_cast<double>(stats.records_seen) / elapsed_seconds
                : 0.0;
        std::string elapsed = format_duration(elapsed_seconds, false);
        std::string rate = format_rate(records_per_second);
        std::string eta = "unknown";

        if (total_.available && total_.total > 0) {
            uint64_t capped_records = std::min(stats.records_seen, total_.total);
            double pct = 100.0 * static_cast<double>(capped_records) /
                         static_cast<double>(total_.total);
            if (records_per_second > 0.0 && capped_records < total_.total) {
                double remaining_records = static_cast<double>(total_.total - capped_records);
                eta = format_duration(remaining_records / records_per_second, true);
            } else {
                eta = "0s";
            }
            std::fprintf(
                stderr,
                "%sProgress: records %llu/%llu (%.1f%%), split ALT variants %llu, max MAC %u, elapsed %s, ETA %s, rate %s%s",
                prefix,
                static_cast<unsigned long long>(stats.records_seen),
                static_cast<unsigned long long>(total_.total),
                pct,
                static_cast<unsigned long long>(stats.split_variants),
                stats.max_mac,
                elapsed.c_str(),
                eta.c_str(),
                rate.c_str(),
                suffix
            );
        } else if (total_.available) {
            std::fprintf(
                stderr,
                "%sProgress: records %llu/0, split ALT variants %llu, max MAC %u, elapsed %s, ETA unknown, rate %s%s",
                prefix,
                static_cast<unsigned long long>(stats.records_seen),
                static_cast<unsigned long long>(stats.split_variants),
                stats.max_mac,
                elapsed.c_str(),
                rate.c_str(),
                suffix
            );
        } else {
            std::fprintf(
                stderr,
                "%sProgress: records %llu, split ALT variants %llu, max MAC %u, elapsed %s, ETA unknown, rate %s%s",
                prefix,
                static_cast<unsigned long long>(stats.records_seen),
                static_cast<unsigned long long>(stats.split_variants),
                stats.max_mac,
                elapsed.c_str(),
                rate.c_str(),
                suffix
            );
        }

        std::fflush(stderr);
    }

    IndexRecordCount total_;
    uint64_t next_record_report_ = kProgressRecordInterval;
    std::chrono::steady_clock::time_point start_time_;
    std::chrono::steady_clock::time_point last_report_;
    bool stderr_is_tty_ = false;
    bool printed_ = false;
};

static uint64_t parse_u64_arg(const char* value, const char* name) {
    char* end = nullptr;
    unsigned long long parsed = std::strtoull(value, &end, 10);
    if (!end || *end != '\0') die("invalid %s: %s", name, value);
    return static_cast<uint64_t>(parsed);
}

static double parse_double_arg(const char* value, const char* name) {
    char* end = nullptr;
    double parsed = std::strtod(value, &end);
    if (!end || *end != '\0' || parsed < 0.0) die("invalid %s: %s", name, value);
    return parsed;
}

static void print_usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage:\n"
        "  %s genotype.phased.vcf.gz [options]\n\n"
        "Options:\n"
        "  --storage-weight FLOAT       Weight for payload bytes (default: 1.0)\n"
        "  --query-weight FLOAT         Weight for per-query work units (default: 0.0)\n"
        "  --dense-word-cost FLOAT      Cost per dense uint64 word scanned (default: 1.0)\n"
        "  --sparse-carrier-cost FLOAT  Cost per sparse carrier scanned (default: 1.0)\n"
        "  --max-records N              Scan only the first N VCF records (default: all)\n",
        prog
    );
}

static Options parse_options(int argc, char** argv) {
    Options opt;
    if (argc < 2) {
        print_usage(argv[0]);
        std::exit(1);
    }

    opt.genotype_path = argv[1];

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        auto require_value = [&](const char* name) -> const char* {
            if (i + 1 >= argc) die("%s requires a value", name);
            return argv[++i];
        };

        if (arg == "--storage-weight") {
            opt.storage_weight = parse_double_arg(require_value("--storage-weight"), "--storage-weight");
        } else if (arg == "--query-weight") {
            opt.query_weight = parse_double_arg(require_value("--query-weight"), "--query-weight");
        } else if (arg == "--dense-word-cost") {
            opt.dense_word_cost = parse_double_arg(require_value("--dense-word-cost"), "--dense-word-cost");
        } else if (arg == "--sparse-carrier-cost") {
            opt.sparse_carrier_cost = parse_double_arg(require_value("--sparse-carrier-cost"), "--sparse-carrier-cost");
        } else if (arg == "--max-records") {
            opt.max_records = parse_u64_arg(require_value("--max-records"), "--max-records");
        } else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else {
            die("unknown option: %s", arg.c_str());
        }
    }

    if (opt.storage_weight == 0.0 && opt.query_weight == 0.0) {
        die("at least one of --storage-weight or --query-weight must be positive");
    }

    return opt;
}

static void validate_extra_ploidy(const int32_t* sample_gt, int ploidy, const char* chr, int64_t pos, int sample_index) {
    for (int p = 2; p < ploidy; ++p) {
        int32_t gt = sample_gt[p];
        if (gt == bcf_int32_vector_end) break;
        die("non-diploid genotype at %s:%lld sample index %d", chr, static_cast<long long>(pos), sample_index);
    }
}

static ScanStats scan_mac_distribution(const char* path, uint64_t max_records) {
    htsFile* fp = bcf_open(path, "r");
    if (!fp) die("cannot open genotype VCF/BCF: %s", path);

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) die("cannot read genotype header: %s", path);

    ScanStats stats;
    stats.n_samples = bcf_hdr_nsamples(hdr);
    if (stats.n_samples <= 0) die("genotype VCF has no samples");
    stats.n_haps = static_cast<uint64_t>(stats.n_samples) * 2ULL;
    stats.n_words = (stats.n_haps + 63ULL) / 64ULL;

    bcf1_t* rec = bcf_init();
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;
    ProgressReporter progress(path, max_records);

    while (bcf_read(fp, hdr, rec) == 0) {
        if (max_records != 0 && stats.records_seen >= max_records) break;
        ++stats.records_seen;

        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele < 2) {
            progress.maybe_report(stats);
            continue;
        }

        const char* chr = bcf_hdr_id2name(hdr, rec->rid);
        int64_t pos = static_cast<int64_t>(rec->pos) + 1;

        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) die("missing FORMAT/GT at %s:%lld", chr, static_cast<long long>(pos));
        if (ngt % stats.n_samples != 0) {
            die("GT field length is not divisible by sample count at %s:%lld", chr, static_cast<long long>(pos));
        }

        int ploidy = ngt / stats.n_samples;
        if (ploidy < 2) die("expected diploid GT at %s:%lld", chr, static_cast<long long>(pos));

        std::vector<uint32_t> mac_by_alt(static_cast<size_t>(rec->n_allele), 0);

        for (int i = 0; i < stats.n_samples; ++i) {
            const int32_t* sample_gt = gt_arr + static_cast<size_t>(i) * ploidy;
            validate_extra_ploidy(sample_gt, ploidy, chr, pos, i);

            int32_t g0 = sample_gt[0];
            int32_t g1 = sample_gt[1];
            if (g0 == bcf_int32_vector_end || g1 == bcf_int32_vector_end ||
                bcf_gt_is_missing(g0) || bcf_gt_is_missing(g1)) {
                die("missing genotype at %s:%lld sample index %d", chr, static_cast<long long>(pos), i);
            }
            if (!bcf_gt_is_phased(g1)) {
                die("unphased genotype at %s:%lld sample index %d", chr, static_cast<long long>(pos), i);
            }

            int allele0 = bcf_gt_allele(g0);
            int allele1 = bcf_gt_allele(g1);
            if (allele0 > 0 && allele0 < rec->n_allele) ++mac_by_alt[allele0];
            if (allele1 > 0 && allele1 < rec->n_allele) ++mac_by_alt[allele1];
        }

        for (int alt = 1; alt < rec->n_allele; ++alt) {
            uint32_t mac = mac_by_alt[alt];
            ++stats.mac_counts[mac];
            ++stats.split_variants;
            stats.total_alt_carriers += mac;
            stats.max_mac = std::max(stats.max_mac, mac);
        }

        progress.maybe_report(stats);
    }

    if (gt_arr) std::free(gt_arr);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    progress.finish(stats);
    return stats;
}

static ThresholdScore score_threshold(const ScanStats& stats, const Options& opt, uint32_t threshold) {
    ThresholdScore result;
    result.threshold = threshold;

    for (const auto& [mac, count] : stats.mac_counts) {
        if (mac <= threshold) {
            result.rare_variants += count;
            result.rare_carriers += static_cast<uint64_t>(mac) * count;
        } else {
            result.common_variants += count;
        }
    }

    result.common_dense_words = result.common_variants * stats.n_words;
    result.payload_bytes =
        static_cast<long double>(result.rare_carriers) * 8.0L +
        static_cast<long double>(result.common_dense_words) * 8.0L;
    result.query_units =
        static_cast<long double>(result.rare_carriers) * opt.sparse_carrier_cost +
        static_cast<long double>(result.common_dense_words) * opt.dense_word_cost;
    result.score =
        static_cast<long double>(opt.storage_weight) * result.payload_bytes +
        static_cast<long double>(opt.query_weight) * result.query_units;

    return result;
}

static ThresholdScore find_best_threshold(const ScanStats& stats, const Options& opt) {
    ThresholdScore best = score_threshold(stats, opt, 0);

    for (const auto& [mac, _count] : stats.mac_counts) {
        ThresholdScore candidate = score_threshold(stats, opt, mac);
        if (candidate.score < best.score ||
            (candidate.score == best.score && candidate.threshold > best.threshold)) {
            best = candidate;
        }
    }

    return best;
}

static std::string format_bytes(long double bytes) {
    const char* units[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
    int unit = 0;
    while (bytes >= 1024.0L && unit < 5) {
        bytes /= 1024.0L;
        ++unit;
    }

    char buf[128];
    std::snprintf(buf, sizeof(buf), "%.2Lf %s", bytes, units[unit]);
    return buf;
}

static uint32_t model_break_even_threshold(const ScanStats& stats, const Options& opt) {
    long double dense_cost =
        static_cast<long double>(opt.storage_weight) * 8.0L * stats.n_words +
        static_cast<long double>(opt.query_weight) * opt.dense_word_cost * stats.n_words;
    long double sparse_per_mac =
        static_cast<long double>(opt.storage_weight) * 8.0L +
        static_cast<long double>(opt.query_weight) * opt.sparse_carrier_cost;

    if (sparse_per_mac <= 0.0L) return 0;
    long double threshold = dense_cost / sparse_per_mac;
    if (threshold > static_cast<long double>(UINT32_MAX)) return UINT32_MAX;
    return static_cast<uint32_t>(threshold);
}

static void print_score_row(const char* label, const ThresholdScore& score) {
    std::cout << std::left << std::setw(14) << label
              << std::right << std::setw(10) << score.threshold
              << std::setw(14) << score.rare_variants
              << std::setw(14) << score.common_variants
              << std::setw(18) << score.rare_carriers
              << std::setw(16) << format_bytes(score.payload_bytes)
              << std::setw(18) << std::fixed << std::setprecision(0) << score.query_units
              << "\n";
}

int main(int argc, char** argv) {
    Options opt = parse_options(argc, argv);
    ScanStats stats = scan_mac_distribution(opt.genotype_path, opt.max_records);

    if (stats.split_variants == 0) die("no split ALT variants found");

    uint32_t break_even = model_break_even_threshold(stats, opt);
    ThresholdScore best = find_best_threshold(stats, opt);

    std::cout << "MAC threshold estimator\n";
    std::cout << "Samples:              " << stats.n_samples << "\n";
    std::cout << "Haplotypes:           " << stats.n_haps << "\n";
    std::cout << "Dense words/variant:  " << stats.n_words << "\n";
    std::cout << "Dense bytes/variant:  " << stats.n_words * 8ULL << "\n";
    std::cout << "VCF records scanned:  " << stats.records_seen << "\n";
    std::cout << "Split ALT variants:   " << stats.split_variants << "\n";
    std::cout << "Max MAC observed:     " << stats.max_mac << "\n";
    std::cout << "Mean MAC observed:    "
              << std::fixed << std::setprecision(2)
              << static_cast<double>(stats.total_alt_carriers) / static_cast<double>(stats.split_variants)
              << "\n";
    std::cout << "Model break-even MAC: " << break_even << "\n";
    std::cout << "Recommended threshold:" << best.threshold << "\n\n";

    std::cout << std::left << std::setw(14) << "candidate"
              << std::right << std::setw(10) << "threshold"
              << std::setw(14) << "rare_vars"
              << std::setw(14) << "common_vars"
              << std::setw(18) << "rare_carriers"
              << std::setw(16) << "payload"
              << std::setw(18) << "query_units"
              << "\n";

    print_score_row("best", best);
    print_score_row("break_even", score_threshold(stats, opt, break_even));

    std::vector<uint32_t> candidates = {64, 128, 256, 512, 1024, 2048, 4096, 8192};
    for (uint32_t candidate : candidates) {
        if (candidate <= stats.max_mac || candidate <= break_even) {
            print_score_row(("T=" + std::to_string(candidate)).c_str(), score_threshold(stats, opt, candidate));
        }
    }

    return 0;
}
