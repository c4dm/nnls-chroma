/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  NNLS-Chroma / Chordino

  Audio feature extraction plugins for chromagram and chord
  estimation.

  Centre for Digital Music, Queen Mary University of London.
  This file copyright 2026 Keith Adler.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of the
  License, or (at your option) any later version.  See the file
  COPYING included with this distribution for more information.
*/

/*
  In-process use of the Chordino, NNLS Chroma and Tuning plugins: run
  them over a mono float buffer without a Vamp host, using the Vamp
  host SDK adapters to do the framing and FFT.  Used by the chordino
  command-line tool, the tests and the WebAssembly module.
*/

#ifndef CHORDINO_ANALYZE_H
#define CHORDINO_ANALYZE_H

#include <string>
#include <vector>

namespace chordino {

struct Options {
    bool useNNLS;        // approximate transcription (NNLS) on/off
    float rollon;        // spectral roll-on, 0..5 (percent)
    bool tuneLocal;      // local rather than global tuning
    float whitening;     // spectral whitening, 0..1
    float spectralShape; // 0.5..0.9
    float boostN;        // boost likelihood of "N" (no chord), 0..1
    bool harteSyntax;    // chord labels in Harte syntax (C:min) rather than plain (Cm)
    int chromaNormalise; // 0 none, 1 max, 2 L1, 3 L2 (chroma output only)

    Options() : useNNLS(true), rollon(0.f), tuneLocal(false),
                whitening(1.f), spectralShape(0.7f), boostN(0.1f),
                harteSyntax(false), chromaNormalise(0) {}
};

struct Chord {
    double time;       // seconds
    double duration;   // seconds, until the next chord or end of audio
    std::string label; // e.g. "C", "Am", "G7", "N" (no chord)
};

struct ChordNote {
    double time;
    double duration;
    int midi;          // MIDI note number
};

struct ChromaFrame {
    double time;
    float chroma[12];  // A, Bb, B, C, ... Ab (treble range)
    float bass[12];    // same order, bass range
};

struct Result {
    double duration;   // seconds of audio analysed
    std::vector<Chord> chords;
    std::vector<ChordNote> notes;
    std::vector<double> harmonicChange; // one per analysis frame
    std::vector<double> frameTimes;     // matching harmonicChange
};

// Estimated concert pitch in Hz (the Tuning plugin's global estimate).
float estimateTuning(const float *mono, size_t frames, float sampleRate,
                     const Options &opts, std::string *error);

// Chords, chord notes and harmonic change from the Chordino plugin.
bool analyzeChords(const float *mono, size_t frames, float sampleRate,
                   const Options &opts, Result &out, std::string *error);

// Treble and bass chromagrams from the NNLS Chroma plugin.
bool analyzeChroma(const float *mono, size_t frames, float sampleRate,
                   const Options &opts, std::vector<ChromaFrame> &out,
                   std::string *error);

// JSON rendering of a Result (chords and notes; harmonic change omitted
// unless includeChange is set).
std::string toJson(const Result &r, float sampleRate, bool includeChange);

// Names of the twelve chroma bins in output order.
extern const char *const kChromaNames[12];

} // namespace chordino

#endif
