#include <charconv>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "jpeg_io.hpp"
#include "pipeline.hpp"
#include "scan_script.hpp"

namespace {

using namespace jpegultrascan;
using namespace jpegultrascan::detail;

void usage() {
    std::fputs("jpegultrascancpp " JPEGULTRASCAN_VERSION "\n"
               "usage: jpegultrascancpp [options] input.jpg output.jpg\n"
               "  -b N               most successive-approximation bits per band (default: as many as\n"
               "                     libjpeg allows, 10 for 8-bit images; jpegultrascan.pl's is 3)\n"
               "  -t N               threads, 0 for half the cores (default 0)\n"
               "  --dc MODE          upstream (default) or interleaved (like -i -1)\n"
               "  --markers MODE     keep (default): every APPn/COM segment as it is\n"
               "                     jpegtran: what jpegtran -copy all writes\n"
               "  --scans FILE       skip the search, encode with this script (jpegtran -scans syntax)\n"
               "  --log FILE         write what jpegultrascan.pl -v prints\n"
               "  --sizes            print the size of each scan of the output\n"
               "  -q                 quiet\n"
               "The output is the input again when re-encoding doesn't make it smaller.\n",
               stderr);
}

// a whole decimal number up to max
bool parse_number(const char* text, unsigned max, unsigned& out) {
    const char* end = text + std::strlen(text);
    unsigned value = 0;
    const auto [stop, ec] = std::from_chars(text, end, value);
    if (ec != std::errc() || stop != end || stop == text || value > max) return false;
    out = value;
    return true;
}

bool read_file(const char* path, std::vector<std::uint8_t>& data) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    data.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    return !in.bad();
}

bool write_file(const char* path, const std::vector<std::uint8_t>& data) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) return false;
    out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size()));
    out.close();
    return !out.fail();
}

bool write_text(const char* path, const std::vector<std::string>& lines) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    for (const std::string& l : lines) out << l << '\n';
    out.close();
    return !out.fail();
}

int run(int argc, const char* const* argv) {
    Options options;
    const char* scans_path = nullptr;
    const char* log_path = nullptr;
    detail::Markers markers = detail::Markers::Keep;
    bool print_sizes = false;
    bool quiet = false;
    std::vector<const char*> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const bool has_value = i + 1 < argc;
        if (arg == "-b" && has_value) {
            unsigned bits = 0;
            if (!parse_number(argv[++i], 1000, bits)) {
                usage();
                return 2;
            }
            options.max_bits = static_cast<int>(bits);
        } else if (arg == "-t" && has_value) {
            if (!parse_number(argv[++i], 1024, options.threads)) {
                usage();
                return 2;
            }
        } else if (arg == "--dc" && has_value) {
            const std::string m = argv[++i];
            if (m == "upstream") {
                options.dc = DcScans::Upstream;
            } else if (m == "interleaved") {
                options.dc = DcScans::Interleaved;
            } else {
                usage();
                return 2;
            }
        } else if (arg == "--markers" && has_value) {
            const std::string m = argv[++i];
            if (m == "keep") {
                markers = detail::Markers::Keep;
            } else if (m == "jpegtran") {
                markers = detail::Markers::Jpegtran;
            } else {
                usage();
                return 2;
            }
        } else if (arg == "--scans" && has_value) {
            scans_path = argv[++i];
        } else if (arg == "--log" && has_value) {
            log_path = argv[++i];
        } else if (arg == "--sizes") {
            print_sizes = true;
        } else if (arg == "-q") {
            quiet = true;
        } else if (arg == "-h" || arg == "--help") {
            usage();
            return 0;
        } else if (!arg.empty() && arg[0] == '-' && arg != "-") {
            usage();
            return 2;
        } else {
            files.push_back(argv[i]);
        }
    }
    if (files.size() != 2) {
        usage();
        return 2;
    }

    std::vector<std::uint8_t> input;
    if (!read_file(files[0], input)) {
        std::fprintf(stderr, "can't read %s\n", files[0]);
        return 1;
    }
    Decoder source;
    if (!source.read(input, true)) {
        std::fprintf(stderr, "%s: %s\n", files[0], source.error().c_str());
        return 1;
    }
    if (source.warnings() && !quiet) {
        std::fprintf(stderr, "%s: libjpeg warns: %s\n", files[0], source.last_message());
    }

    const auto start = std::chrono::steady_clock::now();
    std::vector<jpeg_scan_info> scans;
    VerboseLog log;
    log.width = source.info()->num_components * 2 + 10;
    log.width2 = static_cast<int>(std::to_string(input.size()).size());
    std::string error;
    if (scans_path) {
        std::vector<std::uint8_t> text;
        if (!read_file(scans_path, text)) {
            std::fprintf(stderr, "can't read %s\n", scans_path);
            return 1;
        }
        const auto parsed = parse_scan_script(std::string_view(reinterpret_cast<const char*>(text.data()), text.size()));
        if (!parsed) {
            std::fprintf(stderr, "bad scan script %s\n", scans_path);
            return 1;
        }
        scans = *parsed;
    } else {
        ScriptResult found = detail::find_best_script(source.source(), options, log_path ? &log : nullptr);
        if (found.status != Status::Optimized) {
            std::fprintf(stderr, "%s: %s\n", files[0], found.error.c_str());
            return 1;
        }
        scans = std::move(found.script.scans);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    // beyond mozjpeg's limit the input stays as it is; a script given by hand fails instead
    const bool encodable = scans.size() <= kMaxScans || scans_path != nullptr;
    std::vector<std::uint8_t> output;
    if (!encodable) {
        output = input;
    } else if (!encode(source.source(), scans, markers, output, error)) {
        std::fprintf(stderr, "%s: %s\n", files[0], error.c_str());
        return 1;
    }
    Decoder check;
    std::string why;
    if (encodable &&
        (!check.read(output, true) || !same_image(source, check, markers == detail::Markers::Keep, why))) {
        std::fprintf(stderr, "%s: output check failed: %s\n", files[0],
                     why.empty() ? check.error().c_str() : why.c_str());
        return 1;
    }
    const std::size_t encoded = output.size();
    const bool improved = encodable && encoded <= input.size();
    if (!write_file(files[1], improved ? output : input)) {
        std::fprintf(stderr, "can't write %s\n", files[1]);
        return 1;
    }
    if (log_path) {
        std::vector<std::string> lines = log.sizes;
        lines.emplace_back("--");
        lines.insert(lines.end(), log.summary.begin(), log.summary.end());
        lines.emplace_back("Best scans:");
        for (const auto& s : scans) lines.push_back(format_scan(s));
        lines.push_back("Size: " + std::to_string(encoded) + " bytes");
        char change[64];
        std::snprintf(change, sizeof change, "Change: %+.6f%%",
                      (static_cast<double>(encoded) / static_cast<double>(input.size()) - 1) * 100);
        lines.emplace_back(change);
        if (!write_text(log_path, lines)) {
            std::fprintf(stderr, "can't write %s\n", log_path);
            return 1;
        }
    }
    if (print_sizes) {
        const auto sizes = scan_sizes(output);
        for (std::size_t i = 0; i < sizes.size(); ++i) {
            std::printf("%s %zu\n", i < scans.size() ? format_scan(scans[i]).c_str() : "?", sizes[i]);
        }
    }
    if (!quiet) {
        std::fprintf(stderr, "%s: %zu -> %zu bytes, %zu scans, %.2f s search%s\n", files[0], input.size(),
                     improved ? encoded : input.size(), scans.size(), seconds,
                     !encodable ? " (kept the input: mozjpeg can't encode that many scans)"
                     : improved ? "" : " (kept the input)");
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {  // NOLINT(misc-const-correctness): the signature main must have
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "jpegultrascancpp: %s\n", e.what());
        return 1;
    }
}
