# Changelog

## 1.2.0 (2026-09-08)

First release since 1.1 (2020). The analysis code is unchanged; the plugin
version number reported to hosts stays at 5.

- **Prebuilt binaries for every platform** from GitHub Actions: universal
  macOS (Apple silicon + Intel), Windows x64, Linux x86_64 and WebAssembly.
  Previously only an Intel macOS 1.1 build was offered.
- **CMake build**, replacing the per-platform Makefiles. The Vamp SDK is
  fetched automatically; Boost is no longer required (the chord dictionary
  parser uses a plain tokenizer).
- **Visual Studio support** (issue #1): builds with MSVC, with a `.def`
  file for the plugin export.
- **`chordino` command-line tool**: chords with timings from WAV, AIFF,
  FLAC, MP3 or Ogg Vorbis, in text, CSV, MIREX `.lab` or JSON, plus chord
  notes, chromagram and tuning output. No Vamp host or external libraries
  needed.
- **WebAssembly module** with a small JavaScript wrapper and a demo page.
- Chord dictionary loading is more robust: an empty or malformed
  `chord.dict` falls back to the built-in dictionary instead of yielding
  only "N", and `;` as well as `,` separates fields.
- Tests: synthesized chord progression checked through the in-process
  analysis, the plugin library via the Vamp host SDK, the CLI in every
  output format, and the WebAssembly build. Continuous integration on
  Linux, macOS, Windows and Emscripten.
- Deprecated `sprintf` calls replaced with `snprintf`.

## 1.1 (2020)

See `releasenotes-0.3.txt` and earlier for the history of the original
project.
