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
  chordino: print the chords of an audio file, with timings, using the
  Chordino plugin compiled in directly (no Vamp host required).
*/

#include "cli/analyze.h"
#include "cli/audiofile.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#ifndef CHORDINO_VERSION
#define CHORDINO_VERSION "dev"
#endif

using namespace chordino;
using std::string;

static void usage(FILE *to)
{
    fprintf(to,
"usage: chordino [options] AUDIOFILE\n"
"\n"
"Estimate the chords in an audio file (WAV, AIFF, FLAC, MP3, Ogg Vorbis)\n"
"and print them with timings.\n"
"\n"
"Output:\n"
"  -f, --format FMT     text (default), csv, json or lab\n"
"                       text: one \"time  chord\" line per chord change\n"
"                       csv:  time,duration,chord\n"
"                       json: chords, chord notes and duration\n"
"                       lab:  start<TAB>end<TAB>chord (MIREX style)\n"
"  -o, --output FILE    write to FILE instead of standard output\n"
"      --notes          print the chord notes (MIDI numbers) instead of labels\n"
"      --chroma         print the treble and bass chromagrams as CSV instead\n"
"      --tuning         print the estimated concert pitch and exit\n"
"      --change         include the harmonic change curve in JSON output\n"
"      --harte          label chords in Harte syntax (C:min) instead of Cm\n"
"\n"
"Analysis (defaults are Chordino's, as used for MIREX 2010):\n"
"      --no-nnls        linear spectral mapping instead of NNLS transcription\n"
"      --rollon PCT     spectral roll-on, 0..5 (default 0)\n"
"      --local-tuning   use local rather than global tuning\n"
"      --whitening X    spectral whitening, 0..1 (default 1.0)\n"
"      --shape X        spectral shape, 0.5..0.9 (default 0.7)\n"
"      --boost-n X      boost the no-chord (N) label, 0..1 (default 0.1)\n"
"      --normalise MODE chroma normalisation: none, max, l1, l2 (chroma output)\n"
"\n"
"  -h, --help           this text\n"
"  -V, --version        print the version\n");
}

static bool parseFloat(const char *s, float &out)
{
    char *end = 0;
    out = (float)strtod(s, &end);
    return end && *end == '\0' && end != s;
}

static string chordText(const Result &r)
{
    string s;
    char buf[128];
    for (size_t i = 0; i < r.chords.size(); ++i) {
        snprintf(buf, sizeof buf, "%9.3f  %s\n", r.chords[i].time, r.chords[i].label.c_str());
        s += buf;
    }
    return s;
}

static string chordCsv(const Result &r)
{
    string s = "time,duration,chord\n";
    char buf[160];
    for (size_t i = 0; i < r.chords.size(); ++i) {
        snprintf(buf, sizeof buf, "%.6f,%.6f,\"%s\"\n", r.chords[i].time, r.chords[i].duration, r.chords[i].label.c_str());
        s += buf;
    }
    return s;
}

static string chordLab(const Result &r)
{
    string s;
    char buf[160];
    for (size_t i = 0; i < r.chords.size(); ++i) {
        if (r.chords[i].duration <= 0.0) continue;
        snprintf(buf, sizeof buf, "%.6f\t%.6f\t%s\n", r.chords[i].time,
                 r.chords[i].time + r.chords[i].duration, r.chords[i].label.c_str());
        s += buf;
    }
    return s;
}

static string notesText(const Result &r, const string &format)
{
    string s;
    char buf[160];
    if (format == "csv") s = "time,duration,midi\n";
    for (size_t i = 0; i < r.notes.size(); ++i) {
        if (format == "csv") {
            snprintf(buf, sizeof buf, "%.6f,%.6f,%d\n", r.notes[i].time, r.notes[i].duration, r.notes[i].midi);
        } else if (format == "lab") {
            snprintf(buf, sizeof buf, "%.6f\t%.6f\t%d\n", r.notes[i].time, r.notes[i].time + r.notes[i].duration, r.notes[i].midi);
        } else {
            snprintf(buf, sizeof buf, "%9.3f  %9.3f  %d\n", r.notes[i].time, r.notes[i].duration, r.notes[i].midi);
        }
        s += buf;
    }
    return s;
}

static string chromaCsv(const std::vector<ChromaFrame> &frames)
{
    string s = "time";
    for (int b = 0; b < 12; ++b) { s += ","; s += kChromaNames[b]; }
    for (int b = 0; b < 12; ++b) { s += ",bass_"; s += kChromaNames[b]; }
    s += "\n";
    char buf[64];
    for (size_t i = 0; i < frames.size(); ++i) {
        snprintf(buf, sizeof buf, "%.6f", frames[i].time);
        s += buf;
        for (int b = 0; b < 12; ++b) { snprintf(buf, sizeof buf, ",%.5f", frames[i].chroma[b]); s += buf; }
        for (int b = 0; b < 12; ++b) { snprintf(buf, sizeof buf, ",%.5f", frames[i].bass[b]); s += buf; }
        s += "\n";
    }
    return s;
}

int main(int argc, char **argv)
{
    Options opts;
    string format = "text";
    string outputPath;
    string inputPath;
    bool wantNotes = false, wantChroma = false, wantTuning = false, wantChange = false;

    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        bool hasNext = i + 1 < argc;
        if (a == "-h" || a == "--help") { usage(stdout); return 0; }
        if (a == "-V" || a == "--version") { printf("chordino %s (Chordino plugin by Matthias Mauch, QMUL)\n", CHORDINO_VERSION); return 0; }
        if ((a == "-f" || a == "--format") && hasNext) { format = argv[++i]; continue; }
        if (a.compare(0, 9, "--format=") == 0) { format = a.substr(9); continue; }
        if ((a == "-o" || a == "--output") && hasNext) { outputPath = argv[++i]; continue; }
        if (a == "--notes") { wantNotes = true; continue; }
        if (a == "--chroma") { wantChroma = true; continue; }
        if (a == "--tuning") { wantTuning = true; continue; }
        if (a == "--change") { wantChange = true; continue; }
        if (a == "--harte") { opts.harteSyntax = true; continue; }
        if (a == "--no-nnls") { opts.useNNLS = false; continue; }
        if (a == "--local-tuning") { opts.tuneLocal = true; continue; }
        if ((a == "--rollon" || a == "--whitening" || a == "--shape" || a == "--boost-n") && hasNext) {
            float v;
            if (!parseFloat(argv[i + 1], v)) { fprintf(stderr, "chordino: %s needs a number\n", a.c_str()); return 2; }
            ++i;
            if (a == "--rollon") opts.rollon = v;
            else if (a == "--whitening") opts.whitening = v;
            else if (a == "--shape") opts.spectralShape = v;
            else opts.boostN = v;
            continue;
        }
        if (a == "--normalise" || a == "--normalize") {
            if (!hasNext) { fprintf(stderr, "chordino: %s needs a mode\n", a.c_str()); return 2; }
            string m = argv[++i];
            if (m == "none") opts.chromaNormalise = 0;
            else if (m == "max") opts.chromaNormalise = 1;
            else if (m == "l1") opts.chromaNormalise = 2;
            else if (m == "l2") opts.chromaNormalise = 3;
            else { fprintf(stderr, "chordino: unknown normalisation %s\n", m.c_str()); return 2; }
            continue;
        }
        if (a.size() > 1 && a[0] == '-') {
            fprintf(stderr, "chordino: unknown option %s\n", a.c_str());
            usage(stderr);
            return 2;
        }
        if (!inputPath.empty()) { fprintf(stderr, "chordino: only one input file is supported\n"); return 2; }
        inputPath = a;
    }

    if (inputPath.empty()) { usage(stderr); return 2; }
    if (format != "text" && format != "csv" && format != "json" && format != "lab") {
        fprintf(stderr, "chordino: unknown format %s\n", format.c_str());
        return 2;
    }

    AudioData audio;
    string error;
    if (!readAudioFile(inputPath, audio, &error)) {
        fprintf(stderr, "chordino: %s\n", error.c_str());
        return 1;
    }

    string text;
    if (wantTuning) {
        float hz = estimateTuning(&audio.mono[0], audio.mono.size(), audio.sampleRate, opts, &error);
        if (hz <= 0.f) { fprintf(stderr, "chordino: %s\n", error.c_str()); return 1; }
        char buf[64];
        snprintf(buf, sizeof buf, "%.1f Hz\n", hz);
        text = buf;
    } else if (wantChroma) {
        std::vector<ChromaFrame> frames;
        if (!analyzeChroma(&audio.mono[0], audio.mono.size(), audio.sampleRate, opts, frames, &error)) {
            fprintf(stderr, "chordino: %s\n", error.c_str());
            return 1;
        }
        text = chromaCsv(frames);
    } else {
        Result result;
        if (!analyzeChords(&audio.mono[0], audio.mono.size(), audio.sampleRate, opts, result, &error)) {
            fprintf(stderr, "chordino: %s\n", error.c_str());
            return 1;
        }
        if (format == "json") text = toJson(result, audio.sampleRate, wantChange);
        else if (wantNotes) text = notesText(result, format);
        else if (format == "csv") text = chordCsv(result);
        else if (format == "lab") text = chordLab(result);
        else text = chordText(result);
    }

    if (outputPath.empty()) {
        fputs(text.c_str(), stdout);
    } else {
        std::ofstream out(outputPath.c_str(), std::ios::binary);
        if (!out) { fprintf(stderr, "chordino: cannot write %s\n", outputPath.c_str()); return 1; }
        out << text;
    }
    return 0;
}
