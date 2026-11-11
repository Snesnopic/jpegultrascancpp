#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

#include <jpeglib.h>

namespace jpegultrascan::detail {

constexpr int kMaxCoef = DCTSIZE2 - 1;

// the passes of a band: first passes Al = 0..max_bits, then refinements Ah = 1..max_bits
// (Al = Ah - 1)
constexpr int pass_count(int max_bits) { return 2 * max_bits + 1; }
constexpr int pass_index(int max_bits, int ah, int al) { return ah == 0 ? al : max_bits + ah; }

// sizes of every AC scan of one component, 0 where not measured
class AcSizes {
public:
    AcSizes() = default;
    explicit AcSizes(int max_bits);
    std::uint64_t get(int ss, int se, int ah, int al) const { return sizes_[index(ss, se, ah, al)]; }
    std::uint64_t& at(int ss, int se, int ah, int al) { return sizes_[index(ss, se, ah, al)]; }
    int max_bits() const { return max_bits_; }

private:
    std::size_t index(int ss, int se, int ah, int al) const {
        return (static_cast<std::size_t>(ss) * DCTSIZE2 + static_cast<std::size_t>(se)) *
                   static_cast<std::size_t>(pass_count(max_bits_)) +
               static_cast<std::size_t>(pass_index(max_bits_, ah, al));
    }
    int max_bits_ = 0;
    std::vector<std::uint64_t> sizes_;
};

// what `jpegultrascan.pl -v` prints: the size of each scan the first time it is measured,
// then its summaries. the perl writes the two through separate buffers, so only the order
// within each of them is meaningful.
struct VerboseLog {
    int width = 0;   // ($planes + 1) * 2 + 10
    int width2 = 0;  // digits of the input size
    std::vector<std::string> sizes;
    std::vector<std::string> summary;
    std::unordered_set<std::string> seen;
    void scan(const std::string& text, std::uint64_t size);
};

// where the search gets its sizes from: the bytes each scan takes in an encoded file, its
// DHT segments and SOS header included
class SizeSource {
public:
    SizeSource() = default;
    SizeSource(const SizeSource&) = delete;
    SizeSource& operator=(const SizeSource&) = delete;
    SizeSource(SizeSource&&) = delete;
    SizeSource& operator=(SizeSource&&) = delete;
    virtual ~SizeSource() = default;
    // every band 1 <= Ss <= Se <= 63 of one component, each pass
    virtual bool measure_ac(int component, AcSizes& out) = 0;
    // DC scans of a set of components, indexed by pass_index
    virtual bool measure_dc(const std::vector<int>& components, std::vector<std::uint64_t>& out) = 0;
    // one sequential scan per part
    virtual bool measure_sequential(const std::vector<std::vector<int>>& parts, std::vector<std::uint64_t>& out) = 0;
    virtual const char* error() const { return ""; }
};

enum class DcMode {
    Upstream,     // one DC scan per component competes with one interleaved DC scan
    Interleaved,  // only the interleaved DC scan, as jpegultrascan.pl -i -1
};

struct SearchOptions {
    int max_bits = 3;
    DcMode dc = DcMode::Upstream;
    VerboseLog* log = nullptr;
    std::function<void()> check_stop;  // called between the steps of search_fitting; throws to stop it
};

struct SearchResult {
    std::vector<jpeg_scan_info> scans;  // in the order they go in the file
    std::uint64_t bytes = 0;            // sum of the scan sizes
};

std::optional<SearchResult> search(SizeSource& source, int num_components, const SearchOptions& options,
                                   std::string& error);

// the smallest script of at most max_scans scans among those search looks at, on the same
// sizes; nullopt when not even the fewest scans fit. its time grows with max_scans squared,
// so it is meant for small limits such as mozjpeg's 64
std::optional<SearchResult> search_fitting(SizeSource& source, int num_components, const SearchOptions& options,
                                           std::size_t max_scans, std::string& error);

// the most search_fitting allocates
std::size_t fitting_memory(int num_components, int max_bits, std::size_t max_scans);

}  // namespace jpegultrascan::detail
