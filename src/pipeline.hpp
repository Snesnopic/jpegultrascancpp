#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "jpeg_io.hpp"
#include "jpegultrascan/jpegultrascan.hpp"
#include "search.hpp"

namespace jpegultrascan::detail {

int max_bits_for(const jpeg_decompress_struct& info, int requested);
SearchOptions search_options_for(const jpeg_decompress_struct& info, const Options& options, VerboseLog* log);

// passes on what the engine measures and keeps it, so that the search for a script of fewer
// scans runs on the same sizes without measuring anything again
class Kept : public SizeSource {
public:
    explicit Kept(SizeSource& source) : source_(source) {}

    bool measure_ac(int component, AcSizes& out) override;
    bool measure_dc(const std::vector<int>& components, std::vector<std::uint64_t>& out) override;
    bool measure_sequential(const std::vector<std::vector<int>>& parts, std::vector<std::uint64_t>& out) override;
    const char* error() const override { return source_.error(); }

private:
    SizeSource& source_;
    std::map<int, AcSizes> ac_;
    std::map<std::vector<int>, std::vector<std::uint64_t>> dc_;
    std::map<std::vector<std::vector<int>>, std::vector<std::uint64_t>> sequential_;
};

// the search and, for a script beyond max_scans, the smallest script within max_scans on the
// same sizes; the first script when not even that fits
std::optional<SearchResult> search_within(SizeSource& source, int num_components, const SearchOptions& options,
                                          std::size_t max_scans, std::string& error);

// Optimized with the script, Unsupported, Interrupted or Error
ScriptResult find_best_script(const Source& source, const Options& options, VerboseLog* log,
                              std::size_t max_scans = kMaxScans);

}  // namespace jpegultrascan::detail
