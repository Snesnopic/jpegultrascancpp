#include "jpegultrascan/jpegultrascan.hpp"

#include <algorithm>

#include "engine_direct.hpp"
#include "parallel.hpp"
#include "pipeline.hpp"

namespace jpegultrascan {

namespace detail {

int max_bits_for(const jpeg_decompress_struct& info, int requested) {
    // libjpeg's limit on Ah and Al
    return std::clamp(requested, 0, info.data_precision == 12 ? kMostBits : 10);
}

SearchOptions search_options_for(const jpeg_decompress_struct& info, const Options& options, VerboseLog* log) {
    SearchOptions search_options;
    search_options.max_bits = max_bits_for(info, options.max_bits);
    search_options.dc = options.dc == DcScans::Interleaved ? DcMode::Interleaved : DcMode::Upstream;
    search_options.log = log;
    return search_options;
}

bool Kept::measure_ac(int component, AcSizes& out) {
    if (const auto it = ac_.find(component); it != ac_.end()) {
        out = it->second;
        return true;
    }
    if (!source_.measure_ac(component, out)) return false;
    ac_.emplace(component, out);
    return true;
}

bool Kept::measure_dc(const std::vector<int>& components, std::vector<std::uint64_t>& out) {
    if (const auto it = dc_.find(components); it != dc_.end()) {
        out = it->second;
        return true;
    }
    if (!source_.measure_dc(components, out)) return false;
    dc_.emplace(components, out);
    return true;
}

bool Kept::measure_sequential(const std::vector<std::vector<int>>& parts, std::vector<std::uint64_t>& out) {
    if (const auto it = sequential_.find(parts); it != sequential_.end()) {
        out = it->second;
        return true;
    }
    if (!source_.measure_sequential(parts, out)) return false;
    sequential_.emplace(parts, out);
    return true;
}

std::optional<SearchResult> search_within(SizeSource& source, int num_components, const SearchOptions& options,
                                          std::size_t max_scans, std::string& error) {
    Kept kept(source);
    auto found = search(kept, num_components, options, error);
    if (!found || found->scans.size() <= max_scans) return found;
    SearchOptions fitting = options;
    fitting.log = nullptr;
    std::string why;
    const auto fits = search_fitting(kept, num_components, fitting, max_scans, why);
    return fits ? fits : found;
}

ScriptResult find_best_script(const Source& source, const Options& options, VerboseLog* log, std::size_t max_scans) {
    ScriptResult result;
    try {
        const Loop loop = loop_for(options);
        SearchOptions search_options = search_options_for(*source.info, options, log);
        search_options.check_stop = [&loop] { loop.check_stop(); };
        DirectEngine engine(source, search_options.max_bits, loop, log);
        if (!engine.ok()) {
            result.status = Status::Unsupported;
            result.error = engine.error();
            return result;
        }
        auto found = search_within(engine, source.info->num_components, search_options, max_scans, result.error);
        if (!found) {
            // the engine fails on what libjpeg wouldn't encode
            result.status = engine.ok() ? Status::Error : Status::Unsupported;
            return result;
        }
        result.script = {std::move(found->scans), static_cast<std::size_t>(found->bytes)};
        result.status = Status::Optimized;
    } catch (const Stopped&) {
        result = {};
        result.status = Status::Interrupted;
    } catch (const std::exception& e) {
        result = {};
        result.error = e.what();
    }
    return result;
}

namespace {

// past memory_limit with `extra` bytes besides what the search takes
bool over_memory(const Source& source, const Options& options, std::size_t extra, ScriptResult& result) {
    if (!options.memory_limit) return false;
    const unsigned pieces = options.parallel_for ? options.threads : thread_count(options.threads);
    const std::size_t engine = DirectEngine::memory_for(source, max_bits_for(*source.info, options.max_bits), pieces);
    if (engine + extra <= options.memory_limit) return false;
    result.status = Status::Unsupported;
    result.error = "image too large for the memory limit";
    return true;
}

}  // namespace

}  // namespace detail

ScriptResult find_best_script(j_decompress_ptr source, jvirt_barray_ptr* coefficients, const Options& options) {
    ScriptResult result;
    try {
        detail::Source src{source, coefficients, {}};
        if (!detail::collect_planes(src, result.error)) return result;
        if (detail::over_memory(src, options, 0, result)) return result;
        result = detail::find_best_script(src, options, nullptr);
    } catch (const std::exception& e) {
        result = {};
        result.error = e.what();
        return result;
    }
    if (result.status == Status::Optimized && result.script.scans.size() > detail::kMaxScans) {
        result.status = Status::Unsupported;
        result.error = "mozjpeg can't encode more than " + std::to_string(detail::kMaxScans) + " scans";
    }
    return result;
}

Result optimize(std::span<const std::uint8_t> jpeg, const Options& options, Markers markers) {
    Result result;
    try {
        detail::Decoder source;
        if (!source.read(jpeg, true, options.memory_limit)) {
            if (source.over_limit()) result.status = Status::Unsupported;
            result.error = source.error();
            return result;
        }
        // besides the search: the input's coefficients, the output (which may grow to twice its
        // size), and the coefficients again to check it
        const std::size_t besides = 2 * detail::coefficient_bytes(*source.info()) + 2 * jpeg.size();
        ScriptResult found;
        if (!detail::over_memory(source.source(), options, besides, found)) {
            found = detail::find_best_script(source.source(), options, nullptr);
        }
        result.status = found.status;
        result.script = std::move(found.script);
        result.error = std::move(found.error);
        if (result.status != Status::Optimized) return result;
        if (result.script.scans.size() > detail::kMaxScans) {
            result.status = Status::Unsupported;
            result.error = "mozjpeg can't encode more than " + std::to_string(detail::kMaxScans) + " scans";
            return result;
        }
        // encoding and checking take seconds on a very large image
        const auto stopped = [&] {
            if (!options.stop_requested || !options.stop_requested()) return false;
            result = {};
            result.status = Status::Interrupted;
            return true;
        };
        if (stopped()) return result;
        result.status = Status::Error;
        const detail::Markers mode = markers == Markers::Keep ? detail::Markers::Keep : detail::Markers::Jpegtran;
        if (!detail::encode(source.source(), result.script.scans, mode, result.data, result.error)) {
            result.data.clear();
            return result;
        }
        if (stopped()) return result;
        detail::Decoder check;
        if (!check.read(result.data, true)) {
            result.error = "output doesn't decode: " + check.error();
            result.data.clear();
            return result;
        }
        if (!detail::same_image(source, check, mode == detail::Markers::Keep, result.error)) {
            result.error = "output differs from input: " + result.error;
            result.data.clear();
            return result;
        }
        result.status = result.data.size() < jpeg.size() ? Status::Optimized : Status::NotSmaller;
        if (result.status == Status::NotSmaller) result.data.clear();
    } catch (const std::exception& e) {
        result = {};
        result.error = e.what();
    }
    return result;
}

}  // namespace jpegultrascan
