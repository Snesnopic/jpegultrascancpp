#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "jpeg_io.hpp"
#include "parallel.hpp"
#include "search.hpp"

namespace jpegultrascan::detail {

// past this many tokens (16 MB) a piece of work stops keeping what the first pass of an AC scan
// would emit, and the second pass walks the blocks again
inline constexpr std::size_t kMaxTokens = std::size_t{1} << 22;

// engine B: works out the bytes each scan would take the way libjpeg writes it (optimal
// Huffman tables, DHT and SOS segments, 0xFF stuffing) without encoding anything
class DirectEngine : public SizeSource {
public:
    // fails when libjpeg would refuse to encode the source; error() says why. the loop's
    // stop throws Stopped out of the measure functions
    DirectEngine(const Source& source, int max_bits, Loop loop, VerboseLog* log, std::size_t max_tokens = kMaxTokens);

    bool ok() const { return error_.empty(); }
    // the most the engine allocates for the source, with `pieces` pieces of work at once
    static std::size_t memory_for(const Source& source, int max_bits, unsigned pieces);
    bool measure_ac(int component, AcSizes& out) override;
    bool measure_dc(const std::vector<int>& components, std::vector<std::uint64_t>& out) override;
    bool measure_sequential(const std::vector<std::vector<int>>& parts, std::vector<std::uint64_t>& out) override;
    const char* error() const override { return error_.c_str(); }

    struct Job {
        jpeg_scan_info scan{};
        bool sequential = false;  // part of a sequential script: jchuff.c, not jcphuff.c
        std::uint64_t size = 0;
        bool failed = false;
    };
    // sizes any scans; AC scans of one component at a time, with Al <= max_bits
    bool measure(std::vector<Job>& jobs);

private:
    // one component's blocks in raster order, and per Al up to max_bits a mask per block: bit
    // k set when |coefficient k in zigzag order| >> Al != 0
    struct AcPlane {
        int component = -1;
        std::vector<const JCOEF*> blocks;
        std::vector<std::uint64_t> masks;  // masks[al * blocks.size() + block]

        const std::uint64_t* of(int al) const { return masks.data() + static_cast<std::size_t>(al) * blocks.size(); }
    };

    // what a piece of work reuses from one AC scan to the next
    struct Scratch {
        std::vector<std::uint32_t> busy;
        std::vector<std::uint32_t> tokens;
    };

    void size_of(Job& job) const;
    // an AC scan of the plane's component, given its blocks with something in the band
    bool size_ac(const jpeg_scan_info& scan, std::span<const std::uint32_t> busy, std::vector<std::uint32_t>& tokens,
                 std::uint64_t& size) const;
    bool measure_group(int ss, int al, AcSizes& out, Scratch& scratch) const;
    void build_plane(int component);

    const Source& source_;
    int max_bits_;
    Loop loop_;
    VerboseLog* log_;
    std::size_t max_tokens_;
    int max_coef_bits_;
    std::vector<int> dc_tables_;
    std::vector<int> ac_tables_;
    AcPlane plane_;
    std::string error_;
};

}  // namespace jpegultrascan::detail
