// Modified from jpegultrascan 1.4.0 (choosebest, trybitsplits, trysplits, partitionscans and
// the final scan ordering), Copyright 2015-2024 Aaron Kaluszka, Apache-2.0.
// Changed: ported from Perl to C++; sizes come from a SizeSource instead of jpegtran runs.

#include "search.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <span>
#include <stdexcept>

namespace jpegultrascan::detail {

namespace {

// perl's ~0, the starting minimum of trysplits and of the final choice
constexpr std::uint64_t kInfinity = std::numeric_limits<std::uint64_t>::max();

std::uint64_t add(std::uint64_t a, std::uint64_t b) { return a > kInfinity - b ? kInfinity : a + b; }

struct Scan {
    std::string plane;  // component indices separated by spaces, as the perl writes them
    int ss = 0;
    int se = 0;
    int ah = 0;
    int al = 0;
    bool sequential = false;

    std::string text() const {
        if (sequential) return plane + ";";
        return plane + ": " + std::to_string(ss) + ' ' + std::to_string(se) + ' ' + std::to_string(ah) + ' ' +
               std::to_string(al) + ';';
    }
};

// a candidate's scans are the perl's concatenated string: here a tree of concatenations,
// with a node made only when the candidate wins
struct Node {
    int left = -1;
    int right = -1;
    int scan = -1;
};

struct Candidate {
    std::uint64_t size = 0;  // 0 plays perl's undef
    std::size_t count = 0;   // the newlines choosebest counts
    int node = -1;           // -1 is the empty string
};

class Arena {
public:
    int leaf(Scan s) {
        scans_.push_back(std::move(s));
        nodes_.push_back({-1, -1, static_cast<int>(scans_.size() - 1)});
        return static_cast<int>(nodes_.size() - 1);
    }
    int concat(int a, int b) {
        if (a < 0) return b;
        if (b < 0) return a;
        nodes_.push_back({a, b, -1});
        return static_cast<int>(nodes_.size() - 1);
    }
    void collect(int node, std::vector<const Scan*>& out) const {
        // iterative in-order walk: the trees can be deep
        std::vector<int> stack;
        if (node >= 0) stack.push_back(node);
        while (!stack.empty()) {
            const Node& n = nodes_[static_cast<std::size_t>(stack.back())];
            stack.pop_back();
            if (n.scan >= 0) {
                out.push_back(&scans_[static_cast<std::size_t>(n.scan)]);
            } else {
                stack.push_back(n.right);
                stack.push_back(n.left);
            }
        }
    }
    std::string text(int node) const {
        std::vector<const Scan*> scans;
        collect(node, scans);
        std::string s;
        for (const Scan* scan : scans) s += scan->text() + "\n";
        return s;
    }

private:
    std::vector<Scan> scans_;
    std::vector<Node> nodes_;
};

enum class Pick { Left, Right, Neither };

// perl: a tie goes to the right side only when the left one has more scans, so between
// equal sizes and equal counts the candidate looked at last wins
Pick choosebest(std::uint64_t leftsize, std::size_t leftcount, std::uint64_t rightsize, std::size_t rightcount) {
    if (!leftsize && !rightsize) return Pick::Neither;
    if (!rightsize || (leftsize && leftsize < rightsize)) return Pick::Left;
    if (!leftsize || leftsize > rightsize || leftcount > rightcount) return Pick::Right;
    return Pick::Left;
}

template <class MakeNode>
void consider(Candidate& best, std::uint64_t size, std::size_t count, MakeNode make_node) {
    switch (choosebest(size, count, best.size, best.count)) {
    case Pick::Left:
        best = {size, count, make_node()};
        break;
    case Pick::Right:
        break;
    case Pick::Neither:
        best = {};
        break;
    }
}

// trysplits and trybitsplits for one plane, memoized as in the perl
class BandSearch {
public:
    using SizeFn = std::function<std::uint64_t(int ss, int se, int ah, int al)>;

    BandSearch(Arena& arena, std::string plane, int max_bits, SizeFn size)
        : arena_(arena), plane_(std::move(plane)), max_bits_(max_bits), size_(std::move(size)),
          splits_(DCTSIZE2 * DCTSIZE2),
          bits_(static_cast<std::size_t>(DCTSIZE2 * DCTSIZE2 * pass_count(max_bits))) {}

    Candidate trybitsplits(int a, int b, int c, int d) {
        auto& slot = bits_[(static_cast<std::size_t>(a) * DCTSIZE2 + static_cast<std::size_t>(b)) *
                               static_cast<std::size_t>(pass_count(max_bits_)) +
                           static_cast<std::size_t>(pass_index(max_bits_, c, d))];
        if (slot) return *slot;
        Candidate min{size_(a, b, c, d), 1, arena_.leaf({plane_, a, b, c, d, false})};
        for (int s = a; s < b; ++s) {
            const Candidate x = trybitsplits(a, s, c, d);
            const Candidate y = trybitsplits(s + 1, b, c, d);
            consider(min, add(x.size, y.size), x.count + y.count, [&] { return arena_.concat(x.node, y.node); });
        }
        slot = min;
        return min;
    }

    Candidate trysplits(int a, int b) {
        auto& slot = splits_[static_cast<std::size_t>(a) * DCTSIZE2 + static_cast<std::size_t>(b)];
        if (slot) return *slot;
        Candidate min{kInfinity, 0, -1};
        std::vector<Candidate> parts;
        for (int n = 0; n <= max_bits_; ++n) {
            parts.clear();
            parts.push_back(trybitsplits(a, b, 0, n));
            std::uint64_t size = parts.back().size;
            std::size_t count = parts.back().count;
            for (int i = n; i > 0; --i) {
                parts.push_back(trybitsplits(a, b, i, i - 1));
                size = add(size, parts.back().size);
                count += parts.back().count;
            }
            consider(min, size, count, [&] {
                int node = -1;
                for (const Candidate& p : parts) node = arena_.concat(node, p.node);
                return node;
            });
        }
        for (int s = a; s < b; ++s) {
            const Candidate x = trysplits(a, s);
            const Candidate y = trysplits(s + 1, b);
            consider(min, add(x.size, y.size), x.count + y.count, [&] { return arena_.concat(x.node, y.node); });
        }
        slot = min;
        return min;
    }

private:
    Arena& arena_;
    std::string plane_;
    int max_bits_;
    SizeFn size_;
    std::vector<std::optional<Candidate>> splits_;
    std::vector<std::optional<Candidate>> bits_;
};

std::vector<int> components_of(const std::string& plane) {
    std::vector<int> out;
    int v = -1;
    for (const char c : plane) {
        if (c >= '0' && c <= '9') {
            v = (v < 0 ? 0 : v * 10) + (c - '0');
        } else if (v >= 0) {
            out.push_back(v);
            v = -1;
        }
    }
    if (v >= 0) out.push_back(v);
    return out;
}

// the fields the perl's final ordering reads from a line: /^([^:]*)(?:: (\d+) (\d+) \d+ (\d+))?/,
// with the missing ones of a sequential line taken as 0
struct Line {
    std::string prefix;
    bool has_numbers = false;
    int ss = 0;
    int se = 0;
    int al = 0;
    const Scan* scan = nullptr;
};

Line parse_line(const Scan& s) {
    Line l;
    l.scan = &s;
    if (s.sequential) {
        l.prefix = s.plane + ";";
    } else {
        l.prefix = s.plane;
        l.has_numbers = true;
        l.ss = s.ss;
        l.se = s.se;
        l.al = s.al;
    }
    return l;
}

// the scans in the order the perl writes them
std::optional<SearchResult> in_file_order(const std::vector<const Scan*>& chosen, std::uint64_t bytes,
                                          std::string& error) {
    // stable sort by Se, then by plane; sequential lines by their text
    std::vector<Line> lines;
    lines.reserve(chosen.size());
    for (const Scan* s : chosen) lines.push_back(parse_line(*s));
    std::stable_sort(lines.begin(), lines.end(), [](const Line& a, const Line& b) {
        if (!a.has_numbers || !b.has_numbers) return a.prefix < b.prefix;
        if (a.se != b.se) return a.se < b.se;
        return a.prefix < b.prefix;
    });

    // i goes before j when it codes higher bits of an overlapping band of the same plane, or
    // when j's plane sorts after i's and j reaches at least as far as i
    const std::size_t n = lines.size();
    std::vector<std::vector<std::size_t>> prev(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            const Line& a = lines[i];
            const Line& b = lines[j];
            const bool overlap = (a.ss <= b.ss && b.ss <= a.se) || (b.ss <= a.ss && a.ss <= b.se);
            if ((b.prefix == a.prefix && overlap && a.al > b.al) || (b.prefix > a.prefix && b.se >= a.se)) {
                prev[j].push_back(i);
            }
        }
    }
    // repeatedly take the first line, in sorted order, whose predecessors are all out
    std::vector<std::size_t> order;
    std::vector<bool> done(n, false);
    while (order.size() < n) {
        bool progress = false;
        for (std::size_t i = 0; i < n; ++i) {
            if (done[i]) continue;
            if (!std::all_of(prev[i].begin(), prev[i].end(), [&](std::size_t k) { return done[k]; })) continue;
            order.push_back(i);
            done[i] = true;
            progress = true;
            break;
        }
        if (!progress) {
            error = "scan order has a cycle";
            return std::nullopt;
        }
    }

    SearchResult result;
    result.bytes = bytes;
    for (const std::size_t i : order) {
        const Scan& s = *lines[i].scan;
        jpeg_scan_info info{};
        const std::vector<int> comps = components_of(s.plane);
        if (comps.empty() || comps.size() > MAX_COMPS_IN_SCAN) {
            error = "bad scan " + s.text();
            return std::nullopt;
        }
        info.comps_in_scan = static_cast<int>(comps.size());
        for (std::size_t k = 0; k < comps.size(); ++k) info.component_index[k] = comps[k];
        if (s.sequential) {
            info.Ss = 0;
            info.Se = kMaxCoef;
        } else {
            info.Ss = s.ss;
            info.Se = s.se;
            info.Ah = s.ah;
            info.Al = s.al;
        }
        result.scans.push_back(info);
    }
    return result;
}

// search_fitting makes the same choices as BandSearch and the partitions, keeping the best size
// for each count of scans instead of only the best size

// costs[k]: the smallest size of something coded in k scans, kInfinity when it can't be
using Costs = std::vector<std::uint64_t>;
using CostSpan = std::span<const std::uint64_t>;

// 0 is perl's undef
std::uint64_t cost(std::uint64_t size) { return size ? size : kInfinity; }

// the best of x and y together, for each count up to cap
Costs combine(CostSpan x, CostSpan y, std::size_t cap) {
    Costs out;
    if (x.empty() || y.empty()) return out;
    out.assign(std::min(cap, x.size() + y.size() - 2) + 1, kInfinity);
    for (std::size_t i = 0; i < x.size() && i <= cap; ++i) {
        if (x[i] == kInfinity) continue;
        const std::size_t last = std::min(y.size() - 1, cap - i);
        for (std::size_t j = 0; j <= last; ++j) {
            const std::uint64_t v = add(x[i], y[j]);
            out[i + j] = std::min(out[i + j], v);
        }
    }
    return out;
}

void keep_lower(Costs& into, CostSpan from) {
    if (into.size() < from.size()) into.resize(from.size(), kInfinity);
    for (std::size_t k = 0; k < from.size(); ++k) into[k] = std::min(into[k], from[k]);
}

// how many of the k scans combine took from x to make total
std::size_t from_left(CostSpan x, CostSpan y, std::size_t k, std::uint64_t total) {
    if (total != kInfinity) {
        for (std::size_t i = 0; i <= k && i < x.size(); ++i) {
            if (k - i < y.size() && add(x[i], y[k - i]) == total) return i;
        }
    }
    throw std::logic_error("search_fitting: a count of scans has no way back");
}

// the most scans a script can have: every coefficient in a band of its own in every pass, and
// every component's DC with every bit
std::size_t most_scans(int num_components, int max_bits) {
    const auto passes = static_cast<std::size_t>(pass_count(max_bits));
    return static_cast<std::size_t>(num_components) * (kMaxCoef * passes + static_cast<std::size_t>(max_bits) + 1);
}

// one component's AC scans by count: each pass of a band split into bands of its own
// (trybitsplits), a band coded with n bits as a first pass and n refinements (trysplits' first
// loop), and 1..63 split into such bands (trysplits' second loop)
class FittingBands {
public:
    FittingBands(const AcSizes& sizes, std::size_t cap, const std::function<void()>& check_stop)
        : sizes_(sizes), bits_(sizes.max_bits()), cap_(cap) {
        const Layout l = layout(bits_, cap_, &split_at_, &whole_at_, &top_at_);
        per_pass_ = l.split;
        split_.assign(passes() * l.split, kInfinity);
        whole_.assign(l.whole, kInfinity);
        top_.assign(l.top, kInfinity);
        for (std::size_t pass = 0; pass < passes(); ++pass) {
            if (check_stop) check_stop();
            fill_split(pass);
        }
        for (int a = 1; a <= kMaxCoef; ++a) {
            if (check_stop) check_stop();
            for (int b = a; b <= kMaxCoef; ++b) fill_whole(a, b);
        }
        fill_top();
    }

    CostSpan by_count() const { return top(kMaxCoef); }

    // appends the scans of the best way with `count` scans
    void collect(std::size_t count, const std::string& plane, std::vector<Scan>& out) const {
        collect_top(kMaxCoef, count, plane, out);
    }

    static std::size_t memory(int max_bits, std::size_t cap) {
        const Layout l = layout(max_bits, cap, nullptr, nullptr, nullptr);
        const auto passes = static_cast<std::size_t>(pass_count(max_bits));
        // the tables, their offsets, and the chains of refinements of one band
        return (passes * l.split + l.whole + l.top + (static_cast<std::size_t>(max_bits) + 3) * (cap + 1)) *
                   sizeof(std::uint64_t) +
               (2 * DCTSIZE2 * DCTSIZE2 + DCTSIZE2) * sizeof(std::size_t);
    }

private:
    struct Layout {
        std::size_t split = 0;  // one pass
        std::size_t whole = 0;
        std::size_t top = 0;
    };

    static std::size_t id(int a, int b) { return static_cast<std::size_t>(a) * DCTSIZE2 + static_cast<std::size_t>(b); }
    static std::size_t length(int a, int b) { return static_cast<std::size_t>(b) - static_cast<std::size_t>(a) + 1; }

    // a..b of one pass in k bands: k <= the coefficients; a..b coded whole or 1..b: each
    // coefficient takes at most max_bits + 1 scans
    static std::size_t split_size(int a, int b, std::size_t cap) { return std::min(length(a, b), cap) + 1; }
    static std::size_t whole_size(int a, int b, int bits, std::size_t cap) {
        return std::min((static_cast<std::size_t>(bits) + 1) * length(a, b), cap) + 1;
    }

    static Layout layout(int bits, std::size_t cap, std::vector<std::size_t>* split_at, std::vector<std::size_t>* whole_at,
                         std::vector<std::size_t>* top_at) {
        Layout l;
        if (split_at) split_at->assign(DCTSIZE2 * DCTSIZE2, 0);
        if (whole_at) whole_at->assign(DCTSIZE2 * DCTSIZE2, 0);
        if (top_at) top_at->assign(DCTSIZE2, 0);
        for (int a = 1; a <= kMaxCoef; ++a) {
            for (int b = a; b <= kMaxCoef; ++b) {
                if (split_at) (*split_at)[id(a, b)] = l.split;
                if (whole_at) (*whole_at)[id(a, b)] = l.whole;
                l.split += split_size(a, b, cap);
                l.whole += whole_size(a, b, bits, cap);
            }
        }
        for (int b = 0; b <= kMaxCoef; ++b) {
            if (top_at) (*top_at)[static_cast<std::size_t>(b)] = l.top;
            l.top += whole_size(1, b, bits, cap);
        }
        return l;
    }

    std::size_t passes() const { return static_cast<std::size_t>(pass_count(bits_)); }
    std::size_t pass_at(int ah, int al) const { return static_cast<std::size_t>(pass_index(bits_, ah, al)); }
    int ah_of(std::size_t pass) const { return pass <= static_cast<std::size_t>(bits_) ? 0 : static_cast<int>(pass) - bits_; }
    int al_of(std::size_t pass) const {
        return pass <= static_cast<std::size_t>(bits_) ? static_cast<int>(pass) : static_cast<int>(pass) - bits_ - 1;
    }
    std::uint64_t size(std::size_t pass, int a, int b) const { return cost(sizes_.get(a, b, ah_of(pass), al_of(pass))); }

    CostSpan split(std::size_t pass, int a, int b) const {
        return {split_.data() + pass * per_pass_ + split_at_[id(a, b)], split_size(a, b, cap_)};
    }
    CostSpan whole(int a, int b) const { return {whole_.data() + whole_at_[id(a, b)], whole_size(a, b, bits_, cap_)}; }
    CostSpan top(int b) const {
        return {top_.data() + top_at_[static_cast<std::size_t>(b)], whole_size(1, b, bits_, cap_)};
    }

    void fill_split(std::size_t pass) {
        for (int a = 1; a <= kMaxCoef; ++a) {
            for (int b = a; b <= kMaxCoef; ++b) {
                std::uint64_t* out = split_.data() + pass * per_pass_ + split_at_[id(a, b)];
                out[1] = size(pass, a, b);
                const std::size_t most = split_size(a, b, cap_) - 1;
                for (std::size_t k = 2; k <= most; ++k) {
                    std::uint64_t m = kInfinity;
                    for (int s = a + static_cast<int>(k) - 2; s < b; ++s) {
                        m = std::min(m, add(split(pass, a, s)[k - 1], size(pass, s + 1, b)));
                    }
                    out[k] = m;
                }
            }
        }
    }

    // chain[n]: refinements n..1 of a..b, by count
    std::vector<Costs> chains(int a, int b) const {
        std::vector<Costs> chain{Costs{0}};
        for (int n = 1; n <= bits_; ++n) chain.push_back(combine(split(pass_at(n, n - 1), a, b), chain.back(), cap_));
        return chain;
    }

    void fill_whole(int a, int b) {
        const CostSpan none = split(pass_at(0, 0), a, b);
        Costs best(none.begin(), none.end());
        Costs chain{0};
        for (int n = 1; n <= bits_; ++n) {
            chain = combine(split(pass_at(n, n - 1), a, b), chain, cap_);
            keep_lower(best, combine(split(pass_at(0, n), a, b), chain, cap_));
        }
        std::copy(best.begin(), best.end(), whole_.begin() + static_cast<std::ptrdiff_t>(whole_at_[id(a, b)]));
    }

    void fill_top() {
        top_[top_at_[0]] = 0;
        for (int b = 1; b <= kMaxCoef; ++b) {
            Costs best;
            for (int s = 0; s < b; ++s) keep_lower(best, combine(top(s), whole(s + 1, b), cap_));
            std::copy(best.begin(), best.end(), top_.begin() + static_cast<std::ptrdiff_t>(top_at_[static_cast<std::size_t>(b)]));
        }
    }

    void collect_top(int b, std::size_t k, const std::string& plane, std::vector<Scan>& out) const {
        if (b == 0) return;
        const std::uint64_t total = top(b)[k];
        for (int s = 0; s < b && total != kInfinity; ++s) {
            const CostSpan left = top(s);
            const CostSpan right = whole(s + 1, b);
            for (std::size_t i = 0; i <= k && i < left.size(); ++i) {
                if (k - i < right.size() && add(left[i], right[k - i]) == total) {
                    collect_top(s, i, plane, out);
                    collect_whole(s + 1, b, k - i, plane, out);
                    return;
                }
            }
        }
        throw std::logic_error("search_fitting: a band split has no way back");
    }

    void collect_whole(int a, int b, std::size_t k, const std::string& plane, std::vector<Scan>& out) const {
        const std::uint64_t total = whole(a, b)[k];
        const std::vector<Costs> chain = chains(a, b);
        for (int n = 0; n <= bits_ && total != kInfinity; ++n) {
            const CostSpan first = split(pass_at(0, n), a, b);
            for (std::size_t i = 0; i <= k && i < first.size(); ++i) {
                const Costs& refinements = chain[static_cast<std::size_t>(n)];
                if (k - i >= refinements.size() || add(first[i], refinements[k - i]) != total) continue;
                collect_split(pass_at(0, n), a, b, i, plane, out);
                std::size_t rest = k - i;
                for (int r = n; r >= 1; --r) {
                    const auto below = static_cast<std::size_t>(r - 1);
                    const std::size_t j = from_left(split(pass_at(r, r - 1), a, b), chain[below], rest,
                                                    chain[static_cast<std::size_t>(r)][rest]);
                    collect_split(pass_at(r, r - 1), a, b, j, plane, out);
                    rest -= j;
                }
                return;
            }
        }
        throw std::logic_error("search_fitting: a band has no way back");
    }

    void collect_split(std::size_t pass, int a, int b, std::size_t k, const std::string& plane,
                       std::vector<Scan>& out) const {
        if (k == 1) {
            out.push_back({plane, a, b, ah_of(pass), al_of(pass), false});
            return;
        }
        const std::uint64_t total = split(pass, a, b)[k];
        for (int s = a + static_cast<int>(k) - 2; s < b && total != kInfinity; ++s) {
            if (add(split(pass, a, s)[k - 1], size(pass, s + 1, b)) == total) {
                collect_split(pass, a, s, k - 1, plane, out);
                out.push_back({plane, s + 1, b, ah_of(pass), al_of(pass), false});
                return;
            }
        }
        throw std::logic_error("search_fitting: a pass split has no way back");
    }

    const AcSizes& sizes_;
    int bits_;
    std::size_t cap_;
    std::vector<std::size_t> split_at_;
    std::vector<std::size_t> whole_at_;
    std::vector<std::size_t> top_at_;
    std::size_t per_pass_ = 0;
    Costs split_;
    Costs whole_;
    Costs top_;
};

}  // namespace

AcSizes::AcSizes(int max_bits)
    : max_bits_(max_bits), sizes_(static_cast<std::size_t>(DCTSIZE2 * DCTSIZE2 * pass_count(max_bits)), 0) {}

void VerboseLog::scan(const std::string& text, std::uint64_t size) {
    if (!seen.insert(text).second) return;
    std::string line = text;
    if (static_cast<int>(line.size()) < width) line.append(static_cast<std::size_t>(width) - line.size(), ' ');
    std::string number = std::to_string(size);
    if (static_cast<int>(number.size()) < width2) number.insert(0, static_cast<std::size_t>(width2) - number.size(), ' ');
    sizes.push_back(line + ' ' + number);
}

std::optional<SearchResult> search(SizeSource& source, int num_components, const SearchOptions& options,
                                   std::string& error) {
    const int planes = num_components - 1;
    const int max_bits = options.max_bits;
    VerboseLog* log = options.log;
    Arena arena;

    // ac scans, one component at a time
    Candidate acall{0, 0, -1};
    for (int p = 0; p <= planes; ++p) {
        if (log) log->summary.push_back("Calculating sizes of AC scans, plane " + std::to_string(p));
        const auto sizes = std::make_shared<AcSizes>(max_bits);
        if (!source.measure_ac(p, *sizes)) {
            error = std::string("measuring AC scans failed: ") + source.error();
            return std::nullopt;
        }
        if (log) log->summary.push_back("Searching for best AC scan, plane " + std::to_string(p));
        BandSearch bands(arena, std::to_string(p), max_bits,
                         [sizes](int ss, int se, int ah, int al) { return sizes->get(ss, se, ah, al); });
        const Candidate best = bands.trysplits(1, kMaxCoef);
        acall = {add(acall.size, best.size), acall.count + best.count, arena.concat(acall.node, best.node)};
        if (log) {
            log->summary.push_back("Best AC scan, plane " + std::to_string(p) + ":");
            std::vector<const Scan*> scans;
            arena.collect(best.node, scans);
            for (const Scan* s : scans) log->summary.push_back(s->text());
            log->summary.push_back("Size: " + std::to_string(best.size));
        }
    }

    // dc scans and the sequential scripts, partition by partition
    if (log) log->summary.emplace_back("Calculating sizes and searching for best DC scan");
    Candidate best{kInfinity, 0, -1};
    std::map<std::string, Candidate> dc_best;
    std::map<std::string, std::uint64_t> sequential;
    const auto partition = [&](const std::vector<std::string>& parts) -> bool {
        Candidate dcall{0, 0, -1};
        for (const std::string& part : parts) {
            auto it = dc_best.find(part);
            if (it == dc_best.end()) {
                const auto dc = std::make_shared<std::vector<std::uint64_t>>(static_cast<std::size_t>(pass_count(max_bits)), 0);
                if (!source.measure_dc(components_of(part), *dc)) {
                    error = std::string("measuring DC scans failed: ") + source.error();
                    return false;
                }
                BandSearch dcsearch(arena, part, max_bits, [dc, max_bits](int, int, int ah, int al) {
                    return (*dc)[static_cast<std::size_t>(pass_index(max_bits, ah, al))];
                });
                it = dc_best.emplace(part, dcsearch.trysplits(0, 0)).first;
            }
            const Candidate& c = it->second;
            dcall = {add(dcall.size, c.size), dcall.count + c.count, arena.concat(dcall.node, c.node)};
        }
        consider(best, add(dcall.size, acall.size), dcall.count + acall.count,
                 [&] { return arena.concat(dcall.node, acall.node); });

        bool measured = true;
        for (const std::string& part : parts) measured = measured && sequential.contains(part + ";");
        if (!measured) {
            std::vector<std::vector<int>> comps;
            comps.reserve(parts.size());
            for (const std::string& part : parts) comps.push_back(components_of(part));
            std::vector<std::uint64_t> sizes;
            if (!source.measure_sequential(comps, sizes) || sizes.size() != parts.size()) {
                error = std::string("measuring sequential scans failed: ") + source.error();
                return false;
            }
            for (std::size_t i = 0; i < parts.size(); ++i) sequential.emplace(parts[i] + ";", sizes[i]);
        }
        std::uint64_t allsize = 0;
        for (const std::string& part : parts) allsize = add(allsize, sequential[part + ";"]);
        consider(best, allsize, parts.size(), [&] {
            int node = -1;
            for (const std::string& part : parts) node = arena.concat(node, arena.leaf({part, 0, 0, 0, 0, true}));
            return node;
        });
        return true;
    };
    std::string all;
    std::vector<std::string> separate;
    for (int p = 0; p <= planes; ++p) {
        all += (p ? " " : "") + std::to_string(p);
        separate.push_back(std::to_string(p));
    }
    if (planes > 0 && options.dc == DcMode::Upstream && !partition(separate)) return std::nullopt;
    if (!partition({all})) return std::nullopt;
    if (best.node < 0 || !best.size || best.size == kInfinity) {
        error = "no scan script";
        return std::nullopt;
    }

    std::vector<const Scan*> chosen;
    arena.collect(best.node, chosen);
    if (log) {
        // perl: the best script with the AC part cut out, when it has one
        const std::string text = arena.text(best.node);
        const std::string ac = arena.text(acall.node);
        const std::size_t at = text.find(ac);
        if (at != std::string::npos && !ac.empty()) {
            const std::string dc = text.substr(0, at) + text.substr(at + ac.size());
            log->summary.emplace_back("Best DC scan:");
            std::size_t start = 0;
            while (start < dc.size()) {
                const std::size_t end = dc.find('\n', start);
                log->summary.push_back(dc.substr(start, end - start));
                start = end + 1;
            }
            log->summary.push_back("Size: " + std::to_string(best.size - acall.size));
        }
    }
    return in_file_order(chosen, best.size, error);
}

std::optional<SearchResult> search_fitting(SizeSource& source, int num_components, const SearchOptions& options,
                                           std::size_t max_scans, std::string& error) {
    const int planes = num_components - 1;
    const int max_bits = options.max_bits;
    const std::size_t cap = std::min(max_scans, most_scans(num_components, max_bits));
    const auto check_stop = [&options] {
        if (options.check_stop) options.check_stop();
    };
    if (cap == 0) {
        error = "no scan script has 0 scans";
        return std::nullopt;
    }

    // each component's AC scans, the best size for each count
    std::vector<AcSizes> ac;
    std::vector<Costs> ac_costs;
    for (int p = 0; p <= planes; ++p) {
        check_stop();
        ac.emplace_back(max_bits);
        if (!source.measure_ac(p, ac.back())) {
            error = std::string("measuring AC scans failed: ") + source.error();
            return std::nullopt;
        }
        const FittingBands bands(ac.back(), cap, options.check_stop);
        ac_costs.emplace_back(bands.by_count().begin(), bands.by_count().end());
    }
    std::vector<Costs> ac_prefix{Costs{0}};
    for (const Costs& c : ac_costs) ac_prefix.push_back(combine(ac_prefix.back(), c, cap));

    // the same partitions as search, each with its DC scans or as sequential scans
    std::vector<std::vector<std::string>> partitions;
    std::string all;
    std::vector<std::string> separate;
    for (int p = 0; p <= planes; ++p) {
        all += (p ? " " : "") + std::to_string(p);
        separate.push_back(std::to_string(p));
    }
    if (planes > 0 && options.dc == DcMode::Upstream) partitions.push_back(separate);
    partitions.push_back({all});

    struct Choice {
        std::uint64_t size = kInfinity;
        std::size_t count = 0;
        std::size_t partition = 0;
        bool sequential = false;
    };
    Choice best;
    const auto consider_choice = [&best](const Choice& c) {
        if (c.size != kInfinity && (c.size < best.size || (c.size == best.size && c.count < best.count))) best = c;
    };
    std::vector<std::vector<Costs>> dc_costs(partitions.size());
    std::vector<std::vector<Costs>> dc_prefix(partitions.size());
    for (std::size_t t = 0; t < partitions.size(); ++t) {
        check_stop();
        dc_prefix[t].push_back(Costs{0});
        std::vector<std::vector<int>> comps;
        for (const std::string& part : partitions[t]) {
            comps.push_back(components_of(part));
            std::vector<std::uint64_t> dc(static_cast<std::size_t>(pass_count(max_bits)), 0);
            if (!source.measure_dc(comps.back(), dc)) {
                error = std::string("measuring DC scans failed: ") + source.error();
                return std::nullopt;
            }
            // trysplits(0, 0): n bits are a first pass and n refinements
            Costs d(std::min(cap, static_cast<std::size_t>(max_bits) + 1) + 1, kInfinity);
            for (int n = 0; n <= max_bits && static_cast<std::size_t>(n) + 1 <= cap; ++n) {
                std::uint64_t size = cost(dc[static_cast<std::size_t>(pass_index(max_bits, 0, n))]);
                for (int r = 1; r <= n; ++r) size = add(size, cost(dc[static_cast<std::size_t>(pass_index(max_bits, r, r - 1))]));
                d[static_cast<std::size_t>(n) + 1] = size;
            }
            dc_prefix[t].push_back(combine(dc_prefix[t].back(), d, cap));
            dc_costs[t].push_back(std::move(d));
        }
        const Costs total = combine(ac_prefix.back(), dc_prefix[t].back(), cap);
        for (std::size_t k = 0; k < total.size(); ++k) consider_choice({total[k], k, t, false});

        std::vector<std::uint64_t> sizes;
        if (!source.measure_sequential(comps, sizes) || sizes.size() != comps.size()) {
            error = std::string("measuring sequential scans failed: ") + source.error();
            return std::nullopt;
        }
        std::uint64_t size = 0;
        for (const std::uint64_t s : sizes) size = add(size, cost(s));
        if (comps.size() <= cap) consider_choice({size, comps.size(), t, true});
    }
    if (best.size == kInfinity) {
        error = "no scan script of at most " + std::to_string(max_scans) + " scans";
        return std::nullopt;
    }

    // back from the best count to the scans: dc first, then each component's AC, as search
    // puts them before ordering
    std::vector<Scan> scans;
    const std::vector<std::string>& parts = partitions[best.partition];
    if (best.sequential) {
        for (const std::string& part : parts) scans.push_back({part, 0, 0, 0, 0, true});
    } else {
        const std::vector<Costs>& dprefix = dc_prefix[best.partition];
        const std::size_t ac_count = from_left(ac_prefix.back(), dprefix.back(), best.count, best.size);
        std::vector<std::size_t> dc_counts(parts.size());
        std::size_t rest = best.count - ac_count;
        for (std::size_t i = parts.size(); i-- > 0;) {
            const std::size_t before = from_left(dprefix[i], dc_costs[best.partition][i], rest, dprefix[i + 1][rest]);
            dc_counts[i] = rest - before;
            rest = before;
        }
        std::vector<std::size_t> ac_counts(ac_costs.size());
        rest = ac_count;
        for (std::size_t c = ac_costs.size(); c-- > 0;) {
            const std::size_t before = from_left(ac_prefix[c], ac_costs[c], rest, ac_prefix[c + 1][rest]);
            ac_counts[c] = rest - before;
            rest = before;
        }
        for (std::size_t i = 0; i < parts.size(); ++i) {
            const int n = static_cast<int>(dc_counts[i]) - 1;
            scans.push_back({parts[i], 0, 0, 0, n, false});
            for (int r = n; r >= 1; --r) scans.push_back({parts[i], 0, 0, r, r - 1, false});
        }
        for (std::size_t c = 0; c < ac.size(); ++c) {
            check_stop();
            const FittingBands bands(ac[c], cap, options.check_stop);
            bands.collect(ac_counts[c], std::to_string(c), scans);
        }
    }
    std::vector<const Scan*> chosen;
    chosen.reserve(scans.size());
    for (const Scan& s : scans) chosen.push_back(&s);
    return in_file_order(chosen, best.size, error);
}

std::size_t fitting_memory(int num_components, int max_bits, std::size_t max_scans) {
    const std::size_t cap = std::min(max_scans, most_scans(num_components, max_bits));
    const auto components = static_cast<std::size_t>(num_components);
    // each component's AC sizes and costs, the prefixes, and one FittingBands at a time
    const auto sizes = static_cast<std::size_t>(DCTSIZE2 * DCTSIZE2 * pass_count(max_bits));
    return (components * sizes + 4 * (components + 2) * (cap + 1)) * sizeof(std::uint64_t) +
           FittingBands::memory(max_bits, cap);
}

}  // namespace jpegultrascan::detail
