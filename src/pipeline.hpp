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

// passes on what the engine measures and keeps it, so that the search can run again with fewer
// bits without measuring anything: a scan's size doesn't depend on max_bits
class Kept : public SizeSource {
public:
    Kept(SizeSource& source, int max_bits) : source_(source), max_bits_(max_bits), bits_(max_bits) {}

    // from now on, answer from what was kept as a search with at most `bits` asks
    void limit(int bits) { bits_ = bits; }

    bool measure_ac(int component, AcSizes& out) override;
    bool measure_dc(const std::vector<int>& components, std::vector<std::uint64_t>& out) override;
    bool measure_sequential(const std::vector<std::vector<int>>& parts, std::vector<std::uint64_t>& out) override;
    const char* error() const override { return source_.error(); }

private:
    SizeSource& source_;
    int max_bits_;
    int bits_;
    std::map<int, AcSizes> ac_;
    std::map<std::vector<int>, std::vector<std::uint64_t>> dc_;
    std::map<std::vector<std::vector<int>>, std::vector<std::uint64_t>> sequential_;
};

// the search and, for a script beyond max_scans, the search again with fewer bits on the same
// sizes until one fits; it may still not fit
std::optional<SearchResult> search_within(SizeSource& source, int num_components, const SearchOptions& options,
                                          std::size_t max_scans, std::string& error);

// Optimized with the script, Unsupported, Interrupted or Error
ScriptResult find_best_script(const Source& source, const Options& options, VerboseLog* log,
                              std::size_t max_scans = kMaxScans);

}  // namespace jpegultrascan::detail
