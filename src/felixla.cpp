#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

int felixla_pack_main(int argc, char** argv);
int felixla_rfmix_main(int argc, char** argv);
int felixla_dosage_main(int argc, char** argv);
int felixla_to_vcf_main(int argc, char** argv);
int felixla_extract_main(int argc, char** argv);
int felixla_concat_main(int argc, char** argv);
int felixla_ancestry_export_main(int argc, char** argv);

namespace {

using ToolMain = int (*)(int, char**);

[[noreturn]] void die(const std::string& message) {
    std::cerr << "ERROR: " << message << '\n';
    std::exit(1);
}

void usage(std::ostream& out) {
    out <<
R"(Usage:
  felixla [input] [output] [filters]

Examples:
  felixla --phase-vcf genotype.phased.vcf.gz \
          --flare-vcf flare.anc.vcf.gz \
          --export-felixla \
          --out hybrid/chr22

  felixla --felixla hybrid/chr22 --export-vcf --out hybrid/chr22.roundtrip
  felixla --felixla hybrid/chr22 --export-lai --out hybrid/chr22
  felixla --merge-list chunks.txt --export-felixla --out hybrid/chr1

Input:
  --phase-vcf PATH              Phased diploid genotype VCF/BCF.
  --flare-vcf PATH              FLARE local ancestry VCF/BCF with AN1/AN2.
  --rfmix-msp PATH              RFMix MSP file.
  --tractor-dosage-vcf PATH     TRACTOR dosage VCF/BCF with DS#/ANC# fields.
  --felixla PREFIX              Existing FELIXla prefix.
  --merge-list FILE             One prefix per line, or prefix<TAB>[^]BED.

Output:
  --out PATH                    Output prefix.
  --export-felixla              Write a FELIXla prefix.
  --export-vcf                  Write split-biallelic phased genotype VCF.gz + .tbi.
  --export-lai                  Write local ancestry intervals as <out>.lai.gz.
  --export-global-admixture     Write per-sample ancestry proportions.
  --export-local-admixture      Write per-region ancestry proportions.

Filters:
  --keep FILE                   Sample IDs to retain, one per line.
  --extract FILE                PVAR/VCF alleles to keep; ID column ignored.
  --exclude FILE                PVAR/VCF alleles to drop; applied after --extract.
  --extract-bed FILE            BED intervals to keep; 0-based half-open.
  --exclude-bed FILE            BED intervals to drop; same convention.
  --chr LIST                    Contigs to convert; ranges like 1-22 allowed.
  --region CHR:START-END        Single region, 1-based half-open [START, END).
  --mac INT / --maf FLOAT       Minimum minor allele count/frequency overall.
  --anc-mac INT / --anc-maf F   Same, required of at least one ancestry alone.

Other:
  --n-ancestries INT            Optional; otherwise taken from the input header.
  --threads INT                 BGZF threads for --flare-vcf packing and
                                --export-lai. Default: 1.
  --version, --help
)";
}

void print_version() {
    std::cout
        << "FELIXla CLI v0 for FELIX "
        << "(Full-cohort Efficient Local ancestry-Integrated miXed-model framework), "
        << "tractor_hybrid-compatible layout format_version=1\n";
}

bool starts_with_dash(const std::string& value) {
    return !value.empty() && value[0] == '-';
}

std::string require_value(int& i, int argc, char** argv, const std::string& flag) {
    if (i + 1 >= argc) die(flag + " requires a value");
    return argv[++i];
}

std::vector<char*> make_argv(std::vector<std::string>& args) {
    std::vector<char*> argv;
    argv.reserve(args.size() + 1);
    for (std::string& arg : args) argv.push_back(arg.data());
    argv.push_back(nullptr);
    return argv;
}

int run_tool(ToolMain tool, std::vector<std::string> args) {
    std::vector<char*> argv = make_argv(args);
    return tool(static_cast<int>(args.size()), argv.data());
}

// Export is always BGZF with a tabix index beside it, so a plain .vcf name
// would describe the file wrongly. A full-cohort export is far too large to be
// worth writing uncompressed, so the name is corrected rather than honoured.
std::string vcf_export_path(const std::string& out) {
    if (out.size() >= 7 && out.compare(out.size() - 7, 7, ".vcf.gz") == 0) return out;
    if (out.size() >= 8 && out.compare(out.size() - 8, 8, ".vcf.bgz") == 0) return out;
    if (out.size() >= 4 && out.compare(out.size() - 4, 4, ".vcf") == 0) {
        std::string corrected = out + ".gz";
        std::cerr << "FELIXla: --export vcf writes BGZF with a tabix index; "
                  << "writing " << corrected << " rather than " << out << ".\n";
        return corrected;
    }
    return out + ".vcf.gz";
}

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// The ancestry views each have exactly one representation, so --out names the
// prefix and the suffix is ours to add: .lai.gz is always BGZF, the admixture
// tables are always plain TSV. A name that already reaches partway there --
// run.lai for a file that will be compressed -- is completed rather than
// appended to, so --out run.lai and --out run land on the same file.
std::string suffixed_export_path(const std::string& out, const std::string& suffix,
                                 const char* what) {
    if (ends_with(out, suffix)) return out;

    std::string stem = out;
    for (size_t dot = suffix.find('.', 1); dot != std::string::npos;
         dot = suffix.find('.', dot + 1)) {
        std::string partial = suffix.substr(0, dot);
        if (ends_with(out, partial)) stem = out.substr(0, out.size() - partial.size());
    }
    if (stem == out && (ends_with(out, ".gz") || ends_with(out, ".tsv"))) {
        stem = out.substr(0, out.size() - (ends_with(out, ".gz") ? 3 : 4));
    }

    std::string corrected = stem + suffix;
    if (corrected != out + suffix) {
        std::cerr << "FELIXla: " << what << " writes " << suffix << "; writing "
                  << corrected << " rather than " << out << ".\n";
    }
    return corrected;
}

struct PlinkArgs {
    // Input
    std::string phase_vcf;
    std::string flare_vcf;
    std::string rfmix_msp;
    std::string tractor_dosage_vcf;
    std::string felixla_prefix;
    std::string merge_list;

    // Output
    std::string out_path;
    bool export_felixla = false;
    bool export_vcf = false;
    bool export_lai = false;
    bool export_global_admixture = false;
    bool export_local_admixture = false;

    // Filters
    std::string keep_path;
    std::string extract_path;
    std::string exclude_path;
    std::string extract_bed_path;
    std::string exclude_bed_path;
    std::string chr;
    std::string region;
    std::string mac, maf, anc_mac, anc_maf;

    // Other
    std::string n_ancestries;
    std::string threads;
};

PlinkArgs parse_plink_args(int argc, char** argv) {
    PlinkArgs args;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage(std::cout);
            std::exit(0);
        } else if (arg == "--version") {
            print_version();
            std::exit(0);

        // Input
        } else if (arg == "--phase-vcf") {
            args.phase_vcf = require_value(i, argc, argv, arg);
        } else if (arg == "--flare-vcf") {
            args.flare_vcf = require_value(i, argc, argv, arg);
        } else if (arg == "--rfmix-msp") {
            args.rfmix_msp = require_value(i, argc, argv, arg);
        } else if (arg == "--tractor-dosage-vcf") {
            args.tractor_dosage_vcf = require_value(i, argc, argv, arg);
        } else if (arg == "--felixla") {
            args.felixla_prefix = require_value(i, argc, argv, arg);
        } else if (arg == "--merge-list") {
            args.merge_list = require_value(i, argc, argv, arg);

        // Output
        } else if (arg == "--out") {
            args.out_path = require_value(i, argc, argv, arg);
        } else if (arg == "--export-felixla") {
            args.export_felixla = true;
        } else if (arg == "--export-vcf") {
            args.export_vcf = true;
        } else if (arg == "--export-lai") {
            args.export_lai = true;
        } else if (arg == "--export-global-admixture") {
            args.export_global_admixture = true;
        } else if (arg == "--export-local-admixture") {
            args.export_local_admixture = true;

        // Filters
        } else if (arg == "--keep") {
            args.keep_path = require_value(i, argc, argv, arg);
        } else if (arg == "--extract") {
            args.extract_path = require_value(i, argc, argv, arg);
        } else if (arg == "--exclude") {
            args.exclude_path = require_value(i, argc, argv, arg);
        } else if (arg == "--extract-bed") {
            args.extract_bed_path = require_value(i, argc, argv, arg);
        } else if (arg == "--exclude-bed") {
            args.exclude_bed_path = require_value(i, argc, argv, arg);
        } else if (arg == "--chr") {
            args.chr = require_value(i, argc, argv, arg);
        } else if (arg == "--region") {
            args.region = require_value(i, argc, argv, arg);
        } else if (arg == "--mac") {
            args.mac = require_value(i, argc, argv, arg);
        } else if (arg == "--maf") {
            args.maf = require_value(i, argc, argv, arg);
        } else if (arg == "--anc-mac") {
            args.anc_mac = require_value(i, argc, argv, arg);
        } else if (arg == "--anc-maf") {
            args.anc_maf = require_value(i, argc, argv, arg);

        // Other
        } else if (arg == "--n-ancestries") {
            args.n_ancestries = require_value(i, argc, argv, arg);
        } else if (arg == "--threads") {
            args.threads = require_value(i, argc, argv, arg);

        // Removed, with an explanation rather than "unknown flag"
        } else if (arg == "--make-felixla" || arg == "--make-tractor-hybrid") {
            die(arg + " is now --export-felixla");
        } else if (arg == "--export") {
            die("--export vcf is now --export-vcf");
        } else if (arg == "--pmerge-list" || arg == "--concat-list" ||
                   arg == "--concat" || arg == "--concate" ||
                   arg == "--concatenate") {
            die(arg + " is now --merge-list");
        } else if (arg == "--genotype-vcf") {
            die("--genotype-vcf is now --phase-vcf");
        } else if (arg == "--lai-vcf") {
            die("--lai-vcf is now --flare-vcf");
        } else if (arg == "--tractor-hybrid") {
            die("--tractor-hybrid is now --felixla");
        } else if (arg == "--n-ancestry" || arg == "--ancestries") {
            die(arg + " is now --n-ancestries");
        } else if (arg == "--mac-threshold" || arg == "--recommend-mac-threshold" ||
                   arg == "--choose-mac-threshold" || arg == "--n-samples") {
            die(arg + " has been removed: the sparse/dense threshold is always "
                      "ceil(n_samples / 32), which is where a carrier list "
                      "stops being smaller than a bit vector");
        } else if (arg == "--query" || arg == "--admixture" ||
                   arg == "--compare-vcfs" || arg == "--shapeit-args") {
            die(arg + " has been removed");
        } else {
            die("unknown FELIXla flag: " + arg);
        }
    }
    return args;
}

int run_plink_style(int argc, char** argv) {
    PlinkArgs args = parse_plink_args(argc, argv);

    int n_outputs = args.export_felixla + args.export_vcf + args.export_lai +
                    args.export_global_admixture + args.export_local_admixture;
    if (n_outputs == 0) die("no output requested; see --help for --export-*");
    if (n_outputs > 1) die("choose one --export-* at a time");
    if (args.out_path.empty()) die("--out is required");

    bool from_flare = !args.phase_vcf.empty() && !args.flare_vcf.empty();
    bool from_rfmix = !args.phase_vcf.empty() && !args.rfmix_msp.empty();
    bool from_dosage = !args.tractor_dosage_vcf.empty();
    bool from_prefix = !args.felixla_prefix.empty();
    bool from_merge = !args.merge_list.empty();

    int n_inputs = from_flare + from_rfmix + from_dosage + from_prefix + from_merge;
    if (n_inputs == 0) die("no input given; see --help");
    if (n_inputs > 1) die("choose one input at a time");

    // --extract and --exclude also select alleles out of a packed prefix, so
    // they are the two filters that do not need the source records.
    bool has_site_list = !args.extract_path.empty() || !args.exclude_path.empty();
    bool has_pack_only_filter =
        !args.keep_path.empty() || !args.extract_bed_path.empty() ||
        !args.exclude_bed_path.empty() || !args.chr.empty() ||
        !args.mac.empty() || !args.maf.empty() ||
        !args.anc_mac.empty() || !args.anc_maf.empty();
    bool has_filter = has_site_list || has_pack_only_filter;
    if (has_pack_only_filter && !from_flare) {
        die("--keep, the BED filters and the frequency filters are applied "
            "while packing; they are supported only for --phase-vcf + "
            "--flare-vcf --export-felixla");
    }
    if (has_site_list && !from_flare && !(from_prefix && args.export_felixla)) {
        die("--extract and --exclude are supported for --phase-vcf + "
            "--flare-vcf --export-felixla and for --felixla --export-felixla");
    }

    // Only two paths have anything to hand threads to: the FLARE conversion,
    // where BGZF decompression is most of the work, and the LAI export, where
    // BGZF compression is. Everywhere else the flag would be accepted and do
    // nothing, which reads as a tuning knob that does not work.
    bool threads_apply = (from_flare && args.export_felixla) || args.export_lai;
    if (!args.threads.empty() && !threads_apply) {
        die("--threads applies to --phase-vcf + --flare-vcf --export-felixla "
            "and to --export-lai; the other conversions and exports are "
            "single-threaded");
    }

    if (args.export_felixla) {
        if (from_merge) {
            if (has_filter || !args.region.empty()) {
                die("--merge-list cannot be combined with a filter");
            }
            return run_tool(felixla_concat_main,
                            {"felixla_concat", args.merge_list, args.out_path});
        }
        if (from_flare) {
            std::vector<std::string> pack_args = {
                "flare_subset_to_tractor_hybrid",
                args.phase_vcf,
                args.flare_vcf,
                args.n_ancestries.empty() ? "auto" : args.n_ancestries,
                "auto",
                args.out_path
            };
            if (!args.region.empty()) pack_args.push_back(args.region);
            const std::pair<const char*, const std::string&> options[] = {
                {"--keep", args.keep_path},
                {"--extract", args.extract_path},
                {"--exclude", args.exclude_path},
                {"--extract-bed", args.extract_bed_path},
                {"--exclude-bed", args.exclude_bed_path},
                {"--chr", args.chr},
                {"--mac", args.mac},
                {"--maf", args.maf},
                {"--anc-mac", args.anc_mac},
                {"--anc-maf", args.anc_maf},
                {"--threads", args.threads},
            };
            for (const auto& option : options) {
                if (option.second.empty()) continue;
                pack_args.push_back(option.first);
                pack_args.push_back(option.second);
            }
            return run_tool(felixla_pack_main, pack_args);
        }
        if (from_rfmix) {
            if (!args.region.empty()) {
                die("--region is not supported for --rfmix-msp conversion");
            }
            return run_tool(felixla_rfmix_main, {
                "rfmix_msp_to_tractor_hybrid", args.phase_vcf, args.rfmix_msp,
                args.n_ancestries.empty() ? "auto" : args.n_ancestries,
                "auto", args.out_path});
        }
        if (from_dosage) {
            return run_tool(felixla_dosage_main, {
                "tractor_dosage_vcf_to_hybrid", args.tractor_dosage_vcf,
                args.n_ancestries.empty() ? "auto" : args.n_ancestries,
                "auto", args.out_path});
        }
        // From an existing prefix: select a region, a set of alleles, or both.
        if (args.region.empty() && !has_site_list) {
            die("--felixla --export-felixla needs --region, --extract or "
                "--exclude to select with");
        }
        std::vector<std::string> extract_args = {
            "tractor_hybrid_extract_region", args.felixla_prefix};
        if (!args.region.empty()) extract_args.push_back(args.region);
        extract_args.push_back(args.out_path);
        if (!args.extract_path.empty()) {
            extract_args.push_back("--extract");
            extract_args.push_back(args.extract_path);
        }
        if (!args.exclude_path.empty()) {
            extract_args.push_back("--exclude");
            extract_args.push_back(args.exclude_path);
        }
        return run_tool(felixla_extract_main, extract_args);
    }

    if (!from_prefix) {
        die("--export-vcf, --export-lai and the admixture exports read a "
            "packed prefix; give --felixla");
    }
    // Subsetting happens while packing, so a filter here would have to be
    // silently ignored; extract the region into its own prefix first.
    if (!args.region.empty()) {
        die("--region is not supported for this export; extract the region "
            "with --felixla --region --export-felixla first");
    }
    if (args.export_vcf) {
        return run_tool(felixla_to_vcf_main, {
            "tractor_hybrid_to_vcf", args.felixla_prefix,
            vcf_export_path(args.out_path)});
    }
    if (args.export_lai) {
        return run_tool(felixla_ancestry_export_main, {
            "felixla_ancestry_export", "lai", args.felixla_prefix,
            suffixed_export_path(args.out_path, ".lai.gz", "--export-lai"),
            args.threads.empty() ? "1" : args.threads});
    }
    if (args.export_global_admixture) {
        return run_tool(felixla_ancestry_export_main, {
            "felixla_ancestry_export", "global-admixture", args.felixla_prefix,
            suffixed_export_path(args.out_path, ".global.admixture.tsv",
                                 "--export-global-admixture")});
    }
    return run_tool(felixla_ancestry_export_main, {
        "felixla_ancestry_export", "local-admixture", args.felixla_prefix,
        suffixed_export_path(args.out_path, ".local.admixture.tsv",
                             "--export-local-admixture")});
}

std::vector<std::string> command_args(const std::string& command, int argc, char** argv) {
    std::vector<std::string> args;
    args.push_back(command);
    for (int i = 2; i < argc; ++i) args.push_back(argv[i]);
    return args;
}

int run_compat_command(int argc, char** argv) {
    std::string cmd = argv[1];
    if (cmd == "-h" || cmd == "--help" || cmd == "help") {
        usage(std::cout);
        return 0;
    }
    if (cmd == "--version" || cmd == "version") {
        print_version();
        return 0;
    }
    if (cmd == "from-flare" || cmd == "build-flare" || cmd == "pack-flare" ||
        cmd == "flare" || cmd == "flare_subset_to_tractor_hybrid") {
        return run_tool(felixla_pack_main, command_args("flare_subset_to_tractor_hybrid", argc, argv));
    }
    if (cmd == "from-rfmix" || cmd == "rfmix" || cmd == "rfmix-msp" ||
        cmd == "rfmix_msp_to_tractor_hybrid") {
        return run_tool(felixla_rfmix_main, command_args("rfmix_msp_to_tractor_hybrid", argc, argv));
    }
    if (cmd == "from-tractor-dosage" || cmd == "from-dosage" ||
        cmd == "tractor-dosage" || cmd == "dosage" ||
        cmd == "tractor_dosage_vcf_to_hybrid") {
        return run_tool(felixla_dosage_main, command_args("tractor_dosage_vcf_to_hybrid", argc, argv));
    }
    if (cmd == "to-vcf" || cmd == "export-vcf" || cmd == "vcf" ||
        cmd == "tractor_hybrid_to_vcf") {
        return run_tool(felixla_to_vcf_main, command_args("tractor_hybrid_to_vcf", argc, argv));
    }
    if (cmd == "extract" || cmd == "extract-region" || cmd == "region" ||
        cmd == "tractor_hybrid_extract_region") {
        return run_tool(felixla_extract_main, command_args("tractor_hybrid_extract_region", argc, argv));
    }
    if (cmd == "concat" || cmd == "concate" || cmd == "concatenate" || cmd == "pmerge") {
        return run_tool(felixla_concat_main, command_args("felixla_concat", argc, argv));
    }

    usage(std::cerr);
    die("unknown FELIXla command: " + cmd);
}

}  // namespace

// The PLINK-style front end rewrites its flags into each tool's positional
// form, so a tool that wants to record how it was invoked needs the original.
std::string g_felixla_invocation;

int main(int argc, char** argv) {
    for (int i = 0; i < argc; ++i) {
        if (i > 0) g_felixla_invocation.push_back(' ');
        g_felixla_invocation += argv[i];
    }
    if (argc == 1) {
        usage(std::cerr);
        return 2;
    }
    std::string first = argv[1];
    if (starts_with_dash(first)) return run_plink_style(argc, argv);
    return run_compat_command(argc, argv);
}
