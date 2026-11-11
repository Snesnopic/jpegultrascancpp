// Follows the progressive and sequential Huffman encoders of libjpeg-turbo 3.2.0 (jcphuff.c,
// jchuff.c including jpeg_gen_optimal_table and jpeg_make_c_derived_tbl, jcmarker.c and the
// MCU walk of jctrans.c), which is based in part on the work of the Independent JPEG Group.
// See licenses/libjpeg-turbo.

#include "engine_direct.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <numeric>

#include "scan_script.hpp"

namespace jpegultrascan::detail {

namespace {

constexpr int kNatural[DCTSIZE2] = {
    0,  1,  8,  16, 9,  2,  3,  10, 17, 24, 32, 25, 18, 11, 4,  5,  12, 19, 26, 33, 40, 48,
    41, 34, 27, 20, 13, 6,  7,  14, 21, 28, 35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23,
    30, 37, 44, 51, 58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63,
};

constexpr unsigned kMaxCorrBits = 1000;
constexpr unsigned kMaxEobrun = 0x7FFF;
constexpr int kDummyAl = 9;

std::uint64_t band_mask(int ss, int se) {
    const std::uint64_t upto = se == kMaxCoef ? ~std::uint64_t{0} : (std::uint64_t{1} << (se + 1)) - 1;
    return upto & ~((std::uint64_t{1} << ss) - 1);
}

unsigned magnitude(int v) { return static_cast<unsigned>(v < 0 ? -v : v); }

// the bits of a value as emit_bits gets them: the complement of the magnitude when negative
std::uint32_t value_bits(int v, unsigned a) { return v < 0 ? ~a : a; }

std::uint64_t dht_bytes(int symbols) { return 2 + 2 + 1 + 16 + static_cast<std::uint64_t>(symbols); }
std::uint64_t sos_bytes(int components) { return 2 + 2 + 1 + 2 * static_cast<std::uint64_t>(components) + 3; }

// how many bytes of a word are 0xFF, each of which gets a stuffed zero after it
int ff_bytes(std::uint32_t word) {
    const std::uint32_t t = ~word;
    // the top bit of each byte of t that is zero
    return std::popcount(~(((t & 0x7F7F7F7FU) + 0x7F7F7F7FU) | t | 0x7F7F7F7FU));
}

// emit_bits and flush_bits, counting bytes instead of storing them, four at a time
class ByteCounter {
public:
    // size <= 32
    void put(std::uint64_t code, int size) {
        buffer_ = (buffer_ << size) | (code & ((std::uint64_t{1} << size) - 1));
        bits_ += size;
        if (bits_ >= 32) {
            bits_ -= 32;
            bytes_ += 4 + static_cast<std::uint64_t>(ff_bytes(static_cast<std::uint32_t>(buffer_ >> bits_)));
        }
    }
    void flush() {
        put(0x7F, 7);
        for (; bits_ >= 8; bits_ -= 8) bytes_ += ((buffer_ >> (bits_ - 8)) & 0xFF) == 0xFF ? 2 : 1;
        buffer_ = 0;
        bits_ = 0;
    }
    std::uint64_t bytes() const { return bytes_; }

private:
    std::uint64_t buffer_ = 0;
    int bits_ = 0;
    std::uint64_t bytes_ = 0;
};

struct HuffCode {
    std::array<std::uint32_t, 256> code{};
    std::array<std::uint8_t, 256> size{};
    int symbols = 0;
};

using Counts = std::array<std::int64_t, 257>;

// jpeg_gen_optimal_table, then the codes jpeg_make_c_derived_tbl assigns
bool optimal_code(Counts freq, HuffCode& out) {
    constexpr int kMaxClen = 32;
    int bits[kMaxClen + 2] = {};
    int bit_pos[kMaxClen + 1] = {};
    int codesize[257] = {};
    int nz_index[257] = {};
    int others[257];
    std::fill(std::begin(others), std::end(others), -1);
    int huffval[256] = {};

    freq[256] = 1;
    int num_nz = 0;
    for (int i = 0; i < 257; ++i) {
        if (freq[i]) {
            nz_index[num_nz] = i;
            freq[num_nz] = freq[i];
            ++num_nz;
        }
    }
    for (;;) {
        // the two smallest counts, the larger symbol on ties
        int c1 = -1, c2 = -1;
        std::int64_t v = 1000000000, v2 = 1000000000;
        for (int i = 0; i < num_nz; ++i) {
            if (freq[i] <= v2) {
                if (freq[i] <= v) {
                    c2 = c1;
                    v2 = v;
                    v = freq[i];
                    c1 = i;
                } else {
                    v2 = freq[i];
                    c2 = i;
                }
            }
        }
        if (c2 < 0) break;
        freq[c1] += freq[c2];
        freq[c2] = 1000000001;
        ++codesize[c1];
        while (others[c1] >= 0) {
            c1 = others[c1];
            ++codesize[c1];
        }
        others[c1] = c2;
        ++codesize[c2];
        while (others[c2] >= 0) {
            c2 = others[c2];
            ++codesize[c2];
        }
    }
    for (int i = 0; i < num_nz; ++i) {
        if (codesize[i] > kMaxClen) return false;
        ++bits[codesize[i]];
    }
    int p = 0;
    for (int i = 1; i <= kMaxClen; ++i) {
        bit_pos[i] = p;
        p += bits[i];
    }
    int i = kMaxClen;
    for (; i > 16; --i) {
        while (bits[i] > 0) {
            int j = i - 2;
            while (bits[j] == 0) --j;
            bits[i] -= 2;
            ++bits[i - 1];
            bits[j + 1] += 2;
            --bits[j];
        }
    }
    // drop the pseudo-symbol 256 from the longest length
    while (bits[i] == 0) --i;
    --bits[i];
    for (int k = 0; k < num_nz - 1; ++k) huffval[bit_pos[codesize[k]]++] = nz_index[k];

    std::uint32_t code = 0;
    p = 0;
    for (int len = 1; len <= 16; ++len) {
        for (int k = 0; k < bits[len]; ++k) {
            const int s = huffval[p++];
            out.code[static_cast<std::size_t>(s)] = code++;
            out.size[static_cast<std::size_t>(s)] = static_cast<std::uint8_t>(len);
        }
        code <<= 1;
    }
    out.symbols = p;
    return true;
}

// first pass of a scan: symbol counts per table
struct Gather {
    static constexpr bool kBits = false;
    std::array<Counts, NUM_HUFF_TBLS> dc{};
    std::array<Counts, NUM_HUFF_TBLS> ac{};
    void dc_symbol(int tbl, int s, std::uint32_t = 0, int = 0) {
        ++dc[static_cast<std::size_t>(tbl)][static_cast<std::size_t>(s)];
    }
    void ac_symbol(int tbl, int s, std::uint32_t = 0, int = 0) {
        ++ac[static_cast<std::size_t>(tbl)][static_cast<std::size_t>(s)];
    }
    static void bits(std::uint32_t, int) {}
};

// second pass: the bytes. a symbol's extra bits go out with its code
struct Emit {
    static constexpr bool kBits = true;
    const std::array<HuffCode, NUM_HUFF_TBLS>* dc = nullptr;
    const std::array<HuffCode, NUM_HUFF_TBLS>* ac = nullptr;
    ByteCounter out;
    void dc_symbol(int tbl, int s, std::uint32_t v = 0, int n = 0) { put((*dc)[static_cast<std::size_t>(tbl)], s, v, n); }
    void ac_symbol(int tbl, int s, std::uint32_t v = 0, int n = 0) { put((*ac)[static_cast<std::size_t>(tbl)], s, v, n); }
    void bits(std::uint32_t v, int n) { out.put(v, n); }
    void put(const HuffCode& t, int s, std::uint32_t v, int n) {
        const auto i = static_cast<std::size_t>(s);
        out.put((std::uint64_t{t.code[i]} << n) | (v & ((std::uint64_t{1} << n) - 1)), t.size[i] + n);
    }
};

// the first pass of an AC scan, keeping what the second would emit: a symbol with its extra
// bits (symbol, how many and the bits from bit 0, 8 and 13), or bits alone (bit 31 set)
struct Record {
    static constexpr bool kBits = true;
    Counts counts{};
    std::uint32_t* next = nullptr;
    std::uint32_t* end = nullptr;
    bool kept = true;
    void ac_symbol(int, int s, std::uint32_t v = 0, int n = 0) {
        ++counts[static_cast<std::size_t>(s)];
        keep(static_cast<std::uint32_t>(s), v, n);
    }
    void bits(std::uint32_t v, int n) { keep(std::uint32_t{1} << 31, v, n); }
    void keep(std::uint32_t head, std::uint32_t v, int n) {
        if (next == end) {
            kept = false;
            return;
        }
        *next++ = head | static_cast<std::uint32_t>(n) << 8 | (v & ((std::uint32_t{1} << n) - 1)) << 13;
    }
};

// the second pass from what the first kept
std::uint64_t replay(std::span<const std::uint32_t> tokens, const HuffCode& code) {
    ByteCounter out;
    for (const std::uint32_t t : tokens) {
        const auto n = static_cast<int>((t >> 8) & 31);
        const std::uint64_t v = (t >> 13) & 0xFFFF;
        if (t >> 31) {
            out.put(v, n);
        } else {
            const std::size_t s = t & 0xFF;
            out.put((std::uint64_t{code.code[s]} << n) | v, code.size[s] + n);
        }
    }
    out.flush();
    return out.bytes();
}

struct Context {
    const jpeg_decompress_struct* info;
    const std::vector<Plane>* planes;
    const std::vector<int>* dc_tables;
    const std::vector<int>* ac_tables;
    int max_coef_bits;
};

unsigned ceil_div(unsigned a, unsigned b) { return (a + b - 1) / b; }

// the blocks of a scan in the order jctrans.c hands them to the entropy coder: a single
// component's blocks row by row, or MCU by MCU with dummy blocks at the right and bottom
// edges (zero AC, the DC of the block before)
template <class F>
void walk(const Context& ctx, const jpeg_scan_info& scan, const F& f) {
    const std::vector<Plane>& planes = *ctx.planes;
    if (scan.comps_in_scan == 1) {
        const Plane& p = planes[static_cast<std::size_t>(scan.component_index[0])];
        for (unsigned r = 0; r < p.height_in_blocks; ++r) {
            for (unsigned c = 0; c < p.width_in_blocks; ++c) f(0, static_cast<const JCOEF*>(p.rows[r][c]));
        }
        return;
    }
    const jpeg_decompress_struct& info = *ctx.info;
    const unsigned mcus_per_row = ceil_div(info.image_width, static_cast<unsigned>(info.max_h_samp_factor * DCTSIZE));
    const unsigned mcu_rows = ceil_div(info.image_height, static_cast<unsigned>(info.max_v_samp_factor * DCTSIZE));
    JBLOCK dummy[C_MAX_BLOCKS_IN_MCU] = {};
    for (unsigned my = 0; my < mcu_rows; ++my) {
        for (unsigned mx = 0; mx < mcus_per_row; ++mx) {
            int blkn = 0;
            const JCOEF* prev = nullptr;
            for (int ci = 0; ci < scan.comps_in_scan; ++ci) {
                const int comp = scan.component_index[ci];
                const jpeg_component_info& cinfo = info.comp_info[comp];
                const Plane& p = planes[static_cast<std::size_t>(comp)];
                const auto h = static_cast<unsigned>(cinfo.h_samp_factor);
                const auto v = static_cast<unsigned>(cinfo.v_samp_factor);
                const unsigned last_col_width = (p.width_in_blocks % h) ? p.width_in_blocks % h : h;
                const unsigned last_row_height = (p.height_in_blocks % v) ? p.height_in_blocks % v : v;
                const unsigned blockcnt = mx + 1 < mcus_per_row ? h : last_col_width;
                for (unsigned y = 0; y < v; ++y) {
                    unsigned x = 0;
                    if (my + 1 < mcu_rows || y < last_row_height) {
                        const JBLOCK* row = p.rows[my * v + y] + static_cast<std::size_t>(mx) * h;
                        for (; x < blockcnt; ++x, ++blkn) {
                            prev = row[x];
                            f(ci, prev);
                        }
                    }
                    for (; x < h; ++x, ++blkn) {
                        dummy[blkn][0] = prev[0];
                        prev = dummy[blkn];
                        f(ci, prev);
                    }
                }
            }
        }
    }
}

// encode_mcu_DC_first
template <class Sink>
bool dc_first(const Context& ctx, const jpeg_scan_info& scan, Sink& sink) {
    int last[MAX_COMPS_IN_SCAN] = {};
    int tables[MAX_COMPS_IN_SCAN] = {};
    for (int ci = 0; ci < scan.comps_in_scan; ++ci) {
        tables[ci] = (*ctx.dc_tables)[static_cast<std::size_t>(scan.component_index[ci])];
    }
    bool ok = true;
    walk(ctx, scan, [&](int ci, const JCOEF* block) {
        const int value = static_cast<int>(block[0]) >> scan.Al;
        const int diff = value - last[ci];
        last[ci] = value;
        const unsigned a = magnitude(diff);
        const int nbits = std::bit_width(a);
        if (nbits > ctx.max_coef_bits + 1) ok = false;
        sink.dc_symbol(tables[ci], nbits, value_bits(diff, a), nbits);
    });
    return ok;
}

// encode_mcu_DC_refine
template <class Sink>
void dc_refine(const Context& ctx, const jpeg_scan_info& scan, Sink& sink) {
    walk(ctx, scan, [&](int, const JCOEF* block) {
        sink.bits(static_cast<std::uint32_t>(static_cast<int>(block[0]) >> scan.Al), 1);
    });
}

// an AC scan of one component, walked through the blocks with something in the band (by
// index, in raster order): all the others do is lengthen the EOB run
struct AcWalk {
    const JCOEF* const* blocks = nullptr;
    std::size_t count = 0;                      // the component's blocks
    std::span<const std::uint32_t> busy;
    const std::uint64_t* masks = nullptr;       // for the scan's Al
    const std::uint64_t* next_masks = nullptr;  // for Al + 1, in a refinement: what the earlier passes coded
    int table = 0;
    int max_coef_bits = 0;
};

// n blocks in a row with nothing in the band, each one more in the run
template <class EmitRun>
void skip_blocks(std::size_t n, unsigned& eobrun, const EmitRun& emit_eobrun) {
    while (n) {
        const auto take = static_cast<unsigned>(std::min<std::size_t>(n, kMaxEobrun - eobrun));
        eobrun += take;
        n -= take;
        if (eobrun == kMaxEobrun) emit_eobrun();
    }
}

// encode_mcu_AC_first, the nonzero coefficients found through the masks
template <class Sink>
bool ac_first(const AcWalk& w, const jpeg_scan_info& scan, Sink& sink) {
    const std::uint64_t band = band_mask(scan.Ss, scan.Se);
    const int al = scan.Al;
    unsigned eobrun = 0;
    const auto emit_eobrun = [&] {
        const int nbits = std::bit_width(eobrun) - 1;
        sink.ac_symbol(w.table, nbits << 4, eobrun, nbits);
        eobrun = 0;
    };
    bool ok = true;
    std::size_t next = 0;
    for (const std::uint32_t b : w.busy) {
        skip_blocks(b - next, eobrun, emit_eobrun);
        next = std::size_t{b} + 1;
        const JCOEF* block = w.blocks[b];
        if (eobrun) emit_eobrun();
        int k = scan.Ss;
        for (std::uint64_t z = w.masks[b] & band; z; z &= z - 1) {
            const int pos = std::countr_zero(z);
            int run = pos - k;
            for (; run > 15; run -= 16) sink.ac_symbol(w.table, 0xF0);
            const int v = block[kNatural[pos]];
            const unsigned a = magnitude(v) >> al;
            const int nbits = std::bit_width(a);
            if (nbits > w.max_coef_bits) ok = false;
            sink.ac_symbol(w.table, (run << 4) + nbits, value_bits(v, a), nbits);
            k = pos + 1;
        }
        if (k <= scan.Se && ++eobrun == kMaxEobrun) emit_eobrun();
    }
    skip_blocks(w.count - next, eobrun, emit_eobrun);
    if (eobrun) emit_eobrun();
    return ok;
}

// encode_mcu_AC_refine: correction bits wait in a buffer until the next symbol
template <class Sink>
void ac_refine(const AcWalk& w, const jpeg_scan_info& scan, Sink& sink) {
    const std::uint64_t band = band_mask(scan.Ss, scan.Se);
    const int al = scan.Al;
    std::array<std::uint8_t, kMaxCorrBits> buffer{};
    unsigned eobrun = 0;
    unsigned be = 0;
    const auto emit_buffered = [&](unsigned from, unsigned n) {
        if constexpr (Sink::kBits) {
            while (n) {
                const unsigned m = std::min(n, 16U);
                std::uint32_t v = 0;
                for (unsigned i = 0; i < m; ++i) v = (v << 1) | buffer[from + i];
                sink.bits(v, static_cast<int>(m));
                from += m;
                n -= m;
            }
        }
    };
    const auto emit_eobrun = [&] {
        if (!eobrun) return;
        const int nbits = std::bit_width(eobrun) - 1;
        sink.ac_symbol(w.table, nbits << 4, eobrun, nbits);
        eobrun = 0;
        emit_buffered(0, be);
        be = 0;
    };
    std::size_t next = 0;
    for (const std::uint32_t b : w.busy) {
        skip_blocks(b - next, eobrun, emit_eobrun);
        next = std::size_t{b} + 1;
        const JCOEF* block = w.blocks[b];
        const std::uint64_t now = w.masks[b] & band;
        const std::uint64_t fresh = now & ~w.next_masks[b];
        // the last coefficient that becomes nonzero; zero runs past it fold into the EOB
        const int eob = fresh ? 63 - std::countl_zero(fresh) : -1;
        unsigned br_start = be;
        unsigned br = 0;
        int run = 0;
        int k = scan.Ss;
        for (std::uint64_t z = now; z; z &= z - 1) {
            const int pos = std::countr_zero(z);
            run += pos - k;
            k = pos + 1;
            while (run > 15 && pos <= eob) {
                emit_eobrun();
                sink.ac_symbol(w.table, 0xF0);
                run -= 16;
                emit_buffered(br_start, br);
                br_start = 0;
                br = 0;
            }
            // pos is below 64: it comes from a nonzero mask
            const int v = block[kNatural[pos]];  // NOLINT(clang-analyzer-security.ArrayBound)
            const unsigned a = magnitude(v) >> al;
            if (a > 1) {
                buffer[br_start + br++] = static_cast<std::uint8_t>(a & 1);
                continue;
            }
            emit_eobrun();
            sink.ac_symbol(w.table, (run << 4) + 1, v < 0 ? 0 : 1, 1);
            emit_buffered(br_start, br);
            br_start = 0;
            br = 0;
            run = 0;
        }
        run += scan.Se + 1 - k;
        if (run > 0 || br > 0) {
            ++eobrun;
            be += br;
            if (eobrun == kMaxEobrun || be > kMaxCorrBits - DCTSIZE2 + 1) emit_eobrun();
        }
    }
    skip_blocks(w.count - next, eobrun, emit_eobrun);
    emit_eobrun();
}

// htest_one_block and encode_one_block of jchuff.c
template <class Sink>
bool sequential(const Context& ctx, const jpeg_scan_info& scan, Sink& sink) {
    int last[MAX_COMPS_IN_SCAN] = {};
    int dc_tbl[MAX_COMPS_IN_SCAN] = {};
    int ac_tbl[MAX_COMPS_IN_SCAN] = {};
    for (int ci = 0; ci < scan.comps_in_scan; ++ci) {
        dc_tbl[ci] = (*ctx.dc_tables)[static_cast<std::size_t>(scan.component_index[ci])];
        ac_tbl[ci] = (*ctx.ac_tables)[static_cast<std::size_t>(scan.component_index[ci])];
    }
    bool ok = true;
    walk(ctx, scan, [&](int ci, const JCOEF* block) {
        const int diff = block[0] - last[ci];
        last[ci] = block[0];
        unsigned a = magnitude(diff);
        int nbits = std::bit_width(a);
        if (nbits > ctx.max_coef_bits + 1) ok = false;
        sink.dc_symbol(dc_tbl[ci], nbits, value_bits(diff, a), nbits);
        int run = 0;
        for (int k = 1; k < DCTSIZE2; ++k) {
            const int v = block[kNatural[k]];
            if (!v) {
                ++run;
                continue;
            }
            for (; run > 15; run -= 16) sink.ac_symbol(ac_tbl[ci], 0xF0);
            a = magnitude(v);
            nbits = std::bit_width(a);
            if (nbits > ctx.max_coef_bits) ok = false;
            sink.ac_symbol(ac_tbl[ci], (run << 4) + nbits, value_bits(v, a), nbits);
            run = 0;
        }
        if (run > 0) sink.ac_symbol(ac_tbl[ci], 0);
    });
    return ok;
}

}  // namespace

DirectEngine::DirectEngine(const Source& source, int max_bits, Loop loop, VerboseLog* log, std::size_t max_tokens)
    : source_(source), max_bits_(max_bits), loop_(std::move(loop)), log_(log), max_tokens_(max_tokens),
      max_coef_bits_(source.info->data_precision + 2) {
    huffman_tables(source, dc_tables_, ac_tables_, error_);
}

std::size_t DirectEngine::memory_for(const Source& source, int max_bits, unsigned pieces) {
    std::size_t blocks = 0;
    for (const Plane& p : source.planes) blocks = std::max(blocks, std::size_t{p.width_in_blocks} * p.height_in_blocks);
    const auto bits = static_cast<std::size_t>(max_bits) + 1;
    // block pointers and masks of one component
    const std::size_t plane = blocks * (sizeof(const JCOEF*) + bits * sizeof(std::uint64_t));
    // a piece's list of blocks, and its tokens: no more than kMaxTokens, nor 128 per block
    const std::size_t piece = blocks * sizeof(std::uint32_t) + std::min(kMaxTokens, (blocks + 1) * 128) * sizeof(std::uint32_t);
    // the sizes of every AC scan and the search's memo, for each component and once more
    const std::size_t sizes = (source.planes.size() + 1) * DCTSIZE2 * DCTSIZE2 * (2 * bits - 1) * 48;
    // the search for a script within mozjpeg's limit, when the first one is beyond it
    std::size_t fitting = 0;
    if constexpr (kMaxScans != SIZE_MAX) fitting = fitting_memory(static_cast<int>(source.planes.size()), max_bits, kMaxScans);
    return plane + std::size_t{std::max(pieces, 1U)} * piece + sizes + fitting;
}

void DirectEngine::build_plane(int component) {
    if (plane_.component == component) return;
    const Plane& p = source_.planes[static_cast<std::size_t>(component)];
    const std::size_t count = std::size_t{p.width_in_blocks} * p.height_in_blocks;
    plane_.component = component;
    plane_.blocks.clear();
    plane_.blocks.reserve(count);
    for (unsigned r = 0; r < p.height_in_blocks; ++r) {
        for (unsigned c = 0; c < p.width_in_blocks; ++c) plane_.blocks.push_back(p.rows[r][c]);
    }
    plane_.masks.assign(count * (static_cast<std::size_t>(max_bits_) + 1), 0);
    for (std::size_t b = 0; b < count; ++b) {
        const JCOEF* block = plane_.blocks[b];
        for (int k = 1; k < DCTSIZE2; ++k) {
            const int width = std::bit_width(magnitude(block[kNatural[k]]));
            for (int al = 0; al < width && al <= max_bits_; ++al) {
                plane_.masks[static_cast<std::size_t>(al) * count + b] |= std::uint64_t{1} << k;
            }
        }
    }
}

bool DirectEngine::size_ac(const jpeg_scan_info& scan, std::span<const std::uint32_t> busy,
                           std::vector<std::uint32_t>& tokens, std::uint64_t& size) const {
    const AcWalk walk{plane_.blocks.data(),
                      plane_.blocks.size(),
                      busy,
                      plane_.of(scan.Al),
                      scan.Ah ? plane_.of(scan.Al + 1) : nullptr,
                      ac_tables_[static_cast<std::size_t>(scan.component_index[0])],
                      max_coef_bits_};
    const std::size_t room = std::min(max_tokens_, (busy.size() + 1) * 128);
    if (tokens.size() < room) tokens.resize(room);
    Record record;
    record.next = tokens.data();
    record.end = tokens.data() + room;
    if (scan.Ah) ac_refine(walk, scan, record);
    else if (!ac_first(walk, scan, record)) return false;
    HuffCode code;
    if (!optimal_code(record.counts, code)) return false;
    std::uint64_t bytes = 0;
    if (record.kept) {
        bytes = replay({tokens.data(), record.next}, code);
    } else {
        std::array<HuffCode, NUM_HUFF_TBLS> ac{};
        ac[static_cast<std::size_t>(walk.table)] = code;
        Emit emit;
        emit.ac = &ac;
        if (scan.Ah) ac_refine(walk, scan, emit);
        else ac_first(walk, scan, emit);
        emit.out.flush();
        bytes = emit.out.bytes();
    }
    size = dht_bytes(code.symbols) + sos_bytes(1) + bytes;
    return true;
}

void DirectEngine::size_of(Job& job) const {
    const jpeg_scan_info& scan = job.scan;
    if (scan.Ss > 0 && !job.sequential) {
        const std::uint64_t* masks = plane_.of(scan.Al);
        const std::uint64_t band = band_mask(scan.Ss, scan.Se);
        Scratch scratch;
        for (std::size_t b = 0; b < plane_.blocks.size(); ++b) {
            if (masks[b] & band) scratch.busy.push_back(static_cast<std::uint32_t>(b));
        }
        job.failed = !size_ac(scan, scratch.busy, scratch.tokens, job.size);
        return;
    }
    const Context ctx{source_.info, &source_.planes, &dc_tables_, &ac_tables_, max_coef_bits_};
    // which tables the scan's DHT segments carry
    bool dc_used[NUM_HUFF_TBLS] = {};
    bool ac_used[NUM_HUFF_TBLS] = {};
    for (int ci = 0; ci < scan.comps_in_scan; ++ci) {
        const auto comp = static_cast<std::size_t>(scan.component_index[ci]);
        if (scan.Ss == 0 && scan.Ah == 0) dc_used[dc_tables_[comp]] = true;
        if (scan.Se) ac_used[ac_tables_[comp]] = true;
    }
    const auto code = [&](auto& sink) {
        if (job.sequential) return sequential(ctx, scan, sink);
        if (scan.Ah == 0) return dc_first(ctx, scan, sink);
        dc_refine(ctx, scan, sink);
        return true;
    };
    Gather gather;
    if (!code(gather)) {
        job.failed = true;
        return;
    }
    std::array<HuffCode, NUM_HUFF_TBLS> dc{}, ac{};
    std::uint64_t dht = 0;
    for (int t = 0; t < NUM_HUFF_TBLS; ++t) {
        const auto i = static_cast<std::size_t>(t);
        if (dc_used[t]) {
            if (!optimal_code(gather.dc[i], dc[i])) job.failed = true;
            dht += dht_bytes(dc[i].symbols);
        }
        if (ac_used[t]) {
            if (!optimal_code(gather.ac[i], ac[i])) job.failed = true;
            dht += dht_bytes(ac[i].symbols);
        }
    }
    if (job.failed) return;
    Emit emit;
    emit.dc = &dc;
    emit.ac = &ac;
    code(emit);
    emit.out.flush();
    job.size = dht + sos_bytes(scan.comps_in_scan) + emit.out.bytes();
}

bool DirectEngine::measure(std::vector<Job>& jobs) {
    if (!error_.empty()) return false;
    loop_.check_stop();
    int ac_component = -1;
    for (const Job& job : jobs) {
        const jpeg_scan_info& s = job.scan;
        if (s.comps_in_scan > MAX_COMPS_IN_SCAN) {
            error_ = "Too many color components: " + std::to_string(s.comps_in_scan) + ", max " +
                     std::to_string(MAX_COMPS_IN_SCAN);
            return false;
        }
        int blocks = 0;
        for (int ci = 0; ci < s.comps_in_scan; ++ci) {
            const jpeg_component_info& c = source_.info->comp_info[s.component_index[ci]];
            blocks += c.h_samp_factor * c.v_samp_factor;
        }
        if (s.comps_in_scan > 1 && blocks > C_MAX_BLOCKS_IN_MCU) {
            error_ = "Sampling factors too large for interleaved scan";
            return false;
        }
        if (s.Ss > 0 && !job.sequential) {
            if (ac_component >= 0 && ac_component != s.component_index[0]) {
                error_ = "AC scans of more than one component";
                return false;
            }
            if (s.Al + (s.Ah ? 1 : 0) > max_bits_) {
                error_ = "Al beyond max_bits";
                return false;
            }
            ac_component = s.component_index[0];
        }
    }
    if (ac_component >= 0) build_plane(ac_component);
    parallel_for(loop_, jobs.size(), [&](std::size_t i) {
        loop_.check_stop();
        size_of(jobs[i]);
    });
    if (std::any_of(jobs.begin(), jobs.end(), [](const Job& job) { return job.failed; })) {
        error_ = "DCT coefficient out of range";
        return false;
    }
    if (log_) {
        for (const Job& job : jobs) log_->scan(format_scan(job.scan), job.size);
    }
    return true;
}

// every AC scan of the component that starts at ss and codes bit al: the first pass and,
// below max_bits, the refinement Ah = al + 1, for each Se from 63 down. the blocks with
// something in the band are listed once, and the list shrinks with the band
bool DirectEngine::measure_group(int ss, int al, AcSizes& out, Scratch& scratch) const {
    std::vector<std::uint32_t>& busy = scratch.busy;
    const std::uint64_t* masks = plane_.of(al);
    const std::uint64_t from = band_mask(ss, kMaxCoef);
    busy.clear();
    for (std::size_t b = 0; b < plane_.blocks.size(); ++b) {
        if (masks[b] & from) busy.push_back(static_cast<std::uint32_t>(b));
    }
    jpeg_scan_info scan{};
    scan.comps_in_scan = 1;
    scan.component_index[0] = plane_.component;
    scan.Ss = ss;
    scan.Al = al;
    for (int se = kMaxCoef; se >= ss; --se) {
        loop_.check_stop();
        const std::uint64_t band = band_mask(ss, se);
        if (se < kMaxCoef) std::erase_if(busy, [&](std::uint32_t b) { return !(masks[b] & band); });
        scan.Se = se;
        scan.Ah = 0;
        if (!size_ac(scan, busy, scratch.tokens, out.at(ss, se, 0, al))) return false;
        if (al < max_bits_) {
            scan.Ah = al + 1;
            if (!size_ac(scan, busy, scratch.tokens, out.at(ss, se, al + 1, al))) return false;
        }
    }
    return true;
}

namespace {

DirectEngine::Job make_job(const std::vector<int>& comps, int ss, int se, int ah, int al, bool sequential = false) {
    DirectEngine::Job job;
    job.scan.comps_in_scan = static_cast<int>(comps.size());
    for (std::size_t i = 0; i < comps.size() && i < MAX_COMPS_IN_SCAN; ++i) job.scan.component_index[i] = comps[i];
    job.scan.Ss = ss;
    job.scan.Se = se;
    job.scan.Ah = ah;
    job.scan.Al = al;
    job.sequential = sequential;
    return job;
}

// the "0 1 2: 0 0 0 9;" scans jpegultrascan.pl puts in front of the scans it measures: only for the
// log, which has their sizes too
void add_dummies(std::vector<DirectEngine::Job>& jobs, const std::vector<int>& comps) {
    for (std::size_t i = 0; i < comps.size(); i += MAX_COMPS_IN_SCAN) {
        const std::vector<int> group(comps.begin() + static_cast<std::ptrdiff_t>(i),
                                     comps.begin() + static_cast<std::ptrdiff_t>(std::min(comps.size(), i + MAX_COMPS_IN_SCAN)));
        jobs.push_back(make_job(group, 0, 0, 0, kDummyAl));
    }
}

}  // namespace

bool DirectEngine::measure_ac(int component, AcSizes& out) {
    if (!error_.empty()) return false;
    loop_.check_stop();
    if (log_) {
        std::vector<int> all(static_cast<std::size_t>(source_.info->num_components));
        std::iota(all.begin(), all.end(), 0);
        std::vector<Job> dummies;
        add_dummies(dummies, all);
        if (!measure(dummies)) return false;
    }
    build_plane(component);
    const auto passes = static_cast<std::size_t>(max_bits_) + 1;
    const std::size_t groups = static_cast<std::size_t>(kMaxCoef) * passes;
    // a few pieces per thread even out the groups' work, which falls as ss grows; every group
    // is a piece when the loop's threads aren't known
    const std::size_t pieces = loop_.threads ? std::min<std::size_t>(groups, std::size_t{loop_.threads} * 4) : groups;
    std::vector<char> failed(pieces, 0);
    parallel_for(loop_, pieces, [&](std::size_t piece) {
        Scratch scratch;
        for (std::size_t g = piece; g < groups; g += pieces) {
            const int ss = 1 + static_cast<int>(g / passes);
            const int al = static_cast<int>(g % passes);
            if (!measure_group(ss, al, out, scratch)) {
                failed[piece] = 1;
                return;
            }
        }
    });
    if (std::find(failed.begin(), failed.end(), 1) != failed.end()) {
        error_ = "DCT coefficient out of range";
        return false;
    }
    if (log_) {
        const std::vector<int> plane{component};
        for (int ss = 1; ss <= kMaxCoef; ++ss) {
            for (int se = ss; se <= kMaxCoef; ++se) {
                for (int al = 0; al <= max_bits_; ++al) {
                    log_->scan(format_scan(make_job(plane, ss, se, 0, al).scan), out.get(ss, se, 0, al));
                }
                for (int ah = 1; ah <= max_bits_; ++ah) {
                    log_->scan(format_scan(make_job(plane, ss, se, ah, ah - 1).scan), out.get(ss, se, ah, ah - 1));
                }
            }
        }
    }
    return true;
}

bool DirectEngine::measure_dc(const std::vector<int>& components, std::vector<std::uint64_t>& out) {
    std::vector<Job> jobs;
    for (int al = 0; al <= max_bits_; ++al) jobs.push_back(make_job(components, 0, 0, 0, al));
    for (int ah = 1; ah <= max_bits_; ++ah) jobs.push_back(make_job(components, 0, 0, ah, ah - 1));
    const std::size_t count = jobs.size();
    if (log_) {
        std::vector<int> others;
        for (int c = 0; c < source_.info->num_components; ++c) {
            if (std::find(components.begin(), components.end(), c) == components.end()) others.push_back(c);
        }
        add_dummies(jobs, others);
    }
    if (!measure(jobs)) return false;
    for (std::size_t i = 0; i < count; ++i) {
        const jpeg_scan_info& s = jobs[i].scan;
        out[static_cast<std::size_t>(pass_index(max_bits_, s.Ah, s.Al))] = jobs[i].size;
    }
    return true;
}

bool DirectEngine::measure_sequential(const std::vector<std::vector<int>>& parts, std::vector<std::uint64_t>& out) {
    std::vector<Job> jobs;
    jobs.reserve(parts.size());
    for (const auto& part : parts) jobs.push_back(make_job(part, 0, kMaxCoef, 0, 0, true));
    if (!measure(jobs)) return false;
    out.clear();
    out.reserve(jobs.size());
    for (const Job& job : jobs) out.push_back(job.size);
    return true;
}

}  // namespace jpegultrascan::detail
