#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include <jpeglib.h>

namespace jpegultrascan {

enum class DcScans {
    Upstream,     // jpegultrascan's default: one DC scan per component competes with one for all
    Interleaved,  // only one DC scan for all components, as jpegultrascan -i -1
};

// the host's loop: runs body(i) once for every i < n, on any of its threads (possibly all on
// the calling one), in any order, and returns when all are done or rethrows the first exception
using ParallelFor = std::function<void(std::size_t n, const std::function<void(std::size_t)>& body)>;

// the most successive-approximation bits libjpeg allows, for 12-bit images; 10 for 8-bit ones
inline constexpr int kMostBits = 13;

struct Options {
    int max_bits = kMostBits;  // most successive-approximation bits tried per band (jpegultrascan's
                               // default is 3); capped as libjpeg does. more never makes the
                               // script bigger, as the search sees every candidate of fewer
    unsigned threads = 0;   // without parallel_for, the threads to start (0: half the hardware's);
                            // with it, the threads behind it, which only sizes the pieces of work
                            // (0: unknown). the result never depends on it
    DcScans dc = DcScans::Upstream;
    std::size_t memory_limit = 0;  // most bytes a call may take, estimated before it allocates:
                                   // beyond it, Unsupported. 0: no limit
    ParallelFor parallel_for;              // empty: threads of our own
    std::function<bool()> stop_requested;  // polled between pieces of work, from any thread
                                           // running them; true stops with Interrupted
};

enum class Status {
    Optimized,    // optimize: the output is smaller; find_best_script: a script was found
    NotSmaller,   // optimize: re-encoding didn't make the image smaller
    Unsupported,  // libjpeg can't re-encode the image progressively, it would take more than
                  // memory_limit, or the linked libjpeg is mozjpeg and no script of at most 64
                  // scans was found (it can't encode more)
    Interrupted,  // stop_requested returned true
    Error,
};

struct Script {
    std::vector<jpeg_scan_info> scans;
    std::size_t bytes = 0;  // predicted size of the scans with their tables
};

struct ScriptResult {
    Status status = Status::Error;
    Script script;      // when Optimized, and when Unsupported for mozjpeg's 64 scans
    std::string error;  // why, when Unsupported or Error
};

// the scan script with the smallest total size for an image read with jpeg_read_coefficients;
// source and coefficients must stay alive and unchanged meanwhile
ScriptResult find_best_script(j_decompress_ptr source, jvirt_barray_ptr* coefficients, const Options& options);

enum class Markers {
    Keep,      // every APPn and COM segment, in order and as it is
    Jpegtran,  // what `jpegtran -copy all` writes, for comparing with jpegultrascan.pl
};

struct Result {
    Status status = Status::Error;
    std::vector<std::uint8_t> data;  // the re-encoded image, when Optimized
    Script script;
    std::string error;               // why, when Unsupported or Error
};

// reads one JPEG image, searches, encodes, and checks that the output decodes to the input's
// coefficients. bytes after the end of the image are not kept.
Result optimize(std::span<const std::uint8_t> jpeg, const Options& options, Markers markers = Markers::Keep);

}  // namespace jpegultrascan
