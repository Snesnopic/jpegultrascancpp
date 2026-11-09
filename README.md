# jpegultrascancpp

A C++20 port of [jpegultrascan](https://github.com/MegaByte/jpegultrascan) 1.4.0 by Aaron
Kaluszka: it finds the progressive scan script that makes a JPEG smallest, trying every way
to split the coefficients into bands and bits of successive approximation. The coefficients
never change, so the image stays exactly the same.

Meant to become a submodule of [chisel](https://github.com/Snesnopic/chisel).

## Status

Work in progress. The search is a line-by-line port of the Perl script's, but instead of
encoding thousands of candidate scripts it works out each scan's size the way libjpeg-turbo's
Huffman encoder would, without encoding. The result is the Perl script's output byte for byte,
every scan size included, on all the images tested so far.

## Building

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build

Standalone, the build downloads and compiles libjpeg-turbo 3.2.0 as a static library. When a
parent project already defines `JPEG::JPEG` (chisel does, with mozjpeg), that one is used
instead: a binary can link only one libjpeg. mozjpeg is detected from its `jpeglib.h`, and its
own scan search is switched off (`JCP_FASTEST`) before any encoding. mozjpeg 4.1.5 can't
encode more than 64 scans: `select_scan_parameters` in `jcmaster.c` records each scan's Al in
`actual_Al[64]` (`jcmaster.h`) even with `JCP_FASTEST`, and overruns it past that. The search
then tries fewer bits, and if no script fits, the input is kept. With a fixed mozjpeg, for
example in a fork, it is enough to raise `kMaxScans` in `src/jpeg_io.hpp`.

Options: `JPEGULTRASCAN_BUILD_CLI`, `JPEGULTRASCAN_MOZJPEG` (`AUTO`, `ON`, `OFF`).

Binaries for Linux, macOS, Windows and FreeBSD are built by the GitHub workflow; pushing a
`v*` tag publishes them as a release.

## Command line

    jpegultrascancpp [options] input.jpg output.jpg

`-b N` sets the most successive-approximation bits per band (by default as many as libjpeg
allows, which never gives a bigger file than fewer; upstream's default is 3), `-t N` the threads
(by default half the cores; the result never depends on it), `--dc interleaved` allows only one DC scan for
all components (upstream's `-i -1`), `--markers jpegtran` writes what `jpegtran -copy all`
writes instead of keeping every segment as it is. The output is checked by decoding it and
comparing every coefficient with the input's; if it isn't smaller, the input is copied.

## Library

```cpp
#include <jpegultrascan/jpegultrascan.hpp>

jpegultrascan::Options options;
jpegultrascan::ScriptResult found = jpegultrascan::find_best_script(srcinfo, coef_arrays, options);
if (found.status == jpegultrascan::Status::Optimized) {
    // set cinfo.scan_info / num_scans from found.script.scans, then jpeg_write_coefficients
}
```

`optimize()` does the whole job on a JPEG in memory: search, encode, check. Both report
`Optimized`, `NotSmaller` (optimize only), `Unsupported`, `Interrupted` or `Error`.

The library is reentrant, never prints and never ends the process. It starts no threads of its
own when the host passes its loop as `Options::parallel_for`, which then runs independent pieces
of the search on its threads, in any order or all on the calling thread, with the same result.
`Options::stop_requested` is polled between pieces; exceptions the host's loop throws go back to
the host. `Options::memory_limit` makes an image that would need more `Unsupported` before
anything big is allocated. With mozjpeg, a script beyond its 64 scans is searched again with
fewer bits, on the sizes already measured, until one fits.

Link `jpegultrascan::core`.

## License

Apache-2.0, like jpegultrascan: see `LICENSE` and `NOTICE`. Files derived from the Perl
script say so at the top.
