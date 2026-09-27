// The run record every converter writes beside its prefix.
//
// A prefix is described by two files with different jobs: <prefix>.log is about
// the run -- the build, the invocation, the inputs, the filters, what came out
// -- while <prefix>.meta is about the packed data, holding only what a reader
// needs in order to interpret it. Anything that would differ between two runs
// producing the same file belongs here rather than there.
#ifndef FELIXLA_LOG_H
#define FELIXLA_LOG_H

#include <cerrno>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <string>

#include "felixla_version.h"

// Defined by the front end in felixla.cpp: the flags rewrite into each tool's
// positional form, so recording argv here would report the translation rather
// than what was typed.
extern std::string g_felixla_invocation;

namespace felixla {

inline FILE*& log_sink() {
    static FILE* sink = nullptr;
    return sink;
}

inline std::string utc_now() {
    std::time_t now = std::time(nullptr);
    std::tm utc{};
    if (!gmtime_r(&now, &utc)) return std::string();
    char buffer[32];
    if (!std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &utc)) {
        return std::string();
    }
    return std::string(buffer);
}

inline std::string join_argv(int argc, char** argv) {
    std::string joined;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) joined.push_back(' ');
        joined += argv[i];
    }
    return joined;
}

inline void log_line(const char* fmt, ...) {
    FILE* sink = log_sink();
    if (!sink) return;
    va_list args;
    va_start(args, fmt);
    std::vfprintf(sink, fmt, args);
    va_end(args);
    std::fputc('\n', sink);
    std::fflush(sink);
}

// Called from each tool's die(), so a failed run leaves its own explanation.
inline void log_fatal(const char* fmt, va_list args) {
    FILE* sink = log_sink();
    if (!sink) return;
    std::fputs("ERROR: ", sink);
    std::vfprintf(sink, fmt, args);
    std::fputc('\n', sink);
    std::fflush(sink);
}

// Opens <out_prefix>.log and writes the lines every converter shares. Returns
// false rather than exiting, since each tool has its own die().
inline bool log_open(const char* out_prefix, int argc, char** argv) {
    std::string path = std::string(out_prefix) + ".log";
    FILE* sink = std::fopen(path.c_str(), "w");
    if (!sink) return false;
    log_sink() = sink;

    log_line("FELIXla %s", FELIXLA_VERSION);
    std::string started = utc_now();
    if (!started.empty()) log_line("started %s", started.c_str());
    log_line("command %s",
        g_felixla_invocation.empty() ? join_argv(argc, argv).c_str()
                                     : g_felixla_invocation.c_str());
    log_line("out %s", out_prefix);
    return true;
}

// A log not ending in `completed` marks a run that did not finish, whose
// prefix must not be concatenated.
inline void log_finish() {
    std::string finished = utc_now();
    log_line("completed %s", finished.empty() ? "" : finished.c_str());
    if (log_sink()) {
        std::fclose(log_sink());
        log_sink() = nullptr;
    }
}

}  // namespace felixla

#endif
