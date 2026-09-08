/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  NNLS-Chroma / Chordino

  This file copyright 2026 Keith Adler.  GPL-2.0-or-later, see COPYING.

  WebAssembly entry point: analyse a mono float buffer and return the
  chords as JSON.  Built with Emscripten; see wasm/chordino-wasm.js.
*/

#include "cli/analyze.h"

#include <cstdlib>
#include <cstring>

extern "C" {

// flags: bit 0 no NNLS, bit 1 local tuning, bit 2 Harte syntax,
//        bit 3 include harmonic change curve.
// boostN < 0 keeps the default (0.1).
const char *chordino_analyze_json(const float *mono, int frames, float sampleRate,
                                  int flags, float boostN)
{
    chordino::Options opts;
    opts.useNNLS = !(flags & 1);
    opts.tuneLocal = (flags & 2) != 0;
    opts.harteSyntax = (flags & 4) != 0;
    if (boostN >= 0.f) opts.boostN = boostN;
    chordino::Result r;
    std::string err;
    std::string out;
    if (frames <= 0 || sampleRate <= 0.f ||
        !chordino::analyzeChords(mono, (size_t)frames, sampleRate, opts, r, &err)) {
        out = "{\"error\": \"" + (err.empty() ? std::string("bad input") : err) + "\"}";
    } else {
        out = chordino::toJson(r, sampleRate, (flags & 8) != 0);
    }
    char *buf = (char *)malloc(out.size() + 1);
    memcpy(buf, out.c_str(), out.size() + 1);
    return buf;
}

void chordino_free(const char *p)
{
    free((void *)p);
}

}
