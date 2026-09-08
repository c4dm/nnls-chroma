# NNLS Chroma and Chordino

Chromagram and chord estimation plugins for audio, from the Centre for
Digital Music at Queen Mary University of London, by Matthias Mauch.

- **NNLS Chroma** turns audio into treble and bass chromagrams using an
  approximate transcription (non-negative least squares over a dictionary
  of harmonic note profiles).
- **Chordino** estimates a chord sequence from those chromagrams, with
  HMM/Viterbi smoothing. It was the reference chord estimator in MIREX 2010.
- **Tuning** estimates the concert pitch of a recording.

All three are [Vamp](https://vamp-plugins.org) plugins for hosts such as
Sonic Visualiser, Sonic Annotator and Audacity. Since 1.2 the same code
also ships as a command-line tool and as a WebAssembly module, so you can
get chords out of a file without installing a host at all.

```
$ chordino song.mp3
    0.000  N
    0.046  Ab
    1.625  Bm6
    2.926  Cm
    ...
```

## Downloads

Prebuilt packages for every release are on the
[releases page](https://github.com/keithadler/nnls-chroma/releases):

| Package | Contents |
| --- | --- |
| `nnls-chroma-<ver>-macos-universal.zip` | plugin (`nnls-chroma.dylib`, Apple silicon + Intel), `chordino` tool |
| `nnls-chroma-<ver>-windows-x64.zip` | plugin (`nnls-chroma.dll`), `chordino.exe` |
| `nnls-chroma-<ver>-linux-x86_64.tar.gz` | plugin (`nnls-chroma.so`), `chordino` tool |
| `chordino-<ver>-wasm.zip` | `chordino.wasm` + JavaScript wrapper and a demo page |

To install the plugin, copy `nnls-chroma.dylib`/`.dll`/`.so` together with
`nnls-chroma.n3` and `nnls-chroma.cat` into your Vamp plugin folder:

- macOS: `~/Library/Audio/Plug-Ins/Vamp`
- Windows: `C:\Program Files\Vamp Plugins`
- Linux: `~/vamp`, `~/.vamp` or `/usr/local/lib/vamp`

## The chordino tool

```
chordino [options] AUDIOFILE
```

Reads WAV, AIFF, FLAC, MP3 and Ogg Vorbis (mixed to mono) and prints the
chord sequence with timings. Output formats:

| Option | Output |
| --- | --- |
| *(default)* | `time  chord`, one line per chord change |
| `--format csv` | `time,duration,chord` |
| `--format lab` | `start<TAB>end<TAB>chord`, the MIREX `.lab` layout |
| `--format json` | chords, chord notes (MIDI numbers) and duration; add `--change` for the harmonic change curve |
| `--notes` | the chord notes instead of labels |
| `--chroma` | the treble and bass chromagrams as CSV, one row per frame |
| `--tuning` | the estimated concert pitch, e.g. `439.9 Hz` |

Analysis options mirror the plugin parameters: `--no-nnls`, `--rollon`,
`--local-tuning`, `--whitening`, `--shape`, `--boost-n`, `--harte`
(labels such as `C:min` instead of `Cm`) and `--normalise` for the chroma
output. `chordino --help` lists them all. The defaults are the ones used
for the MIREX 2010 submission, the same as the plugin's.

The tool has no dependencies beyond the C++ runtime; the decoders are
compiled in.

## In the browser

`chordino-<ver>-wasm.zip` contains the plugin compiled with Emscripten.
Decoding is left to the browser:

```js
import { loadChordino } from './chordino-wasm.js';
const chordino = await loadChordino();
const audio = await ctx.decodeAudioData(arrayBuffer);   // Web Audio
const mono = audio.getChannelData(0);                    // mix down if stereo
const { chords } = chordino.analyze(mono, audio.sampleRate, { harte: false });
// chords: [{ time, duration, label }, ...]
```

Open `index.html` from the package (served over HTTP) for a drop-a-file
demo. Analysis of a four-minute track takes a few seconds.

## Building

Requires CMake 3.16+ and a C++11 compiler. The Vamp SDK is fetched
automatically (or point `NNLS_VAMP_SDK_DIR` at an unpacked copy). There
are no other dependencies; Boost is no longer needed.

```
cmake -B build
cmake --build build
ctest --test-dir build
```

This produces `build/plugin/nnls-chroma.<dylib|dll|so>` and
`build/chordino`. `cmake --install build` copies the plugin into the
per-user Vamp folder on macOS and Windows, or `lib/vamp` under the
install prefix on Linux.

For a universal macOS binary add
`-DCMAKE_OSX_ARCHITECTURES="arm64;x86_64"`. On Windows, use the Visual
Studio developer prompt (or `ilammy/msvc-dev-cmd` in CI) with
`-G Ninja`; MinGW also works. For WebAssembly, configure with
`emcmake cmake` and run `node test.mjs` in `build-wasm/wasm`.

The tests synthesize a C, F, G, A minor progression and check it comes
back from the in-process analysis, from the plugin library loaded through
the Vamp host SDK, and from the command-line tool in every format.

The original Makefiles (`Makefile.linux`, `Makefile.osx`, `Makefile.mingw`)
are kept for reference; note they still expect Boost headers that the
sources no longer use.

## Plugin documentation

Parameters and outputs of the three plugins are described in
[README](README) (plain text) and in the RDF description
`nnls-chroma.n3`.

## Citation

If you use these plugins in research, please cite:

> Matthias Mauch and Simon Dixon, "Approximate Note Transcription for
> the Improved Identification of Difficult Chords", in Proceedings of the
> 11th International Society for Music Information Retrieval Conference
> (ISMIR 2010), 2010.

See [CITATION](CITATION).

## Licence

GNU General Public License, version 2 or later. See [COPYING](COPYING).
Copyright 2008-2020 Matthias Mauch and QMUL; command-line tool, tests,
build system and WebAssembly port copyright 2026 Keith Adler.
Third-party decoders in `third_party/` (dr_libs, stb_vorbis) are public
domain / MIT-0.
