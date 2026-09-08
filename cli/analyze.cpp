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

#include "analyze.h"

#include <vamp-hostsdk/PluginInputDomainAdapter.h>
#include <vamp-hostsdk/PluginBufferingAdapter.h>

#include "Chordino.h"
#include "NNLSChroma.h"
#include "Tuning.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <iomanip>

using namespace Vamp;
using namespace Vamp::HostExt;

namespace chordino {

const char *const kChromaNames[12] = {
    "A", "Bb", "B", "C", "C#", "D", "Eb", "E", "F", "F#", "G", "Ab"
};

namespace {

void applyOptions(Plugin *p, const Options &o)
{
    // Identifiers are those declared by the plugins; unknown ones are
    // ignored by setParameter, so it is safe to set the union here.
    p->setParameter("useNNLS", o.useNNLS ? 1.f : 0.f);
    p->setParameter("rollon", o.rollon);
    p->setParameter("tuningmode", o.tuneLocal ? 1.f : 0.f);
    p->setParameter("whitening", o.whitening);
    p->setParameter("s", o.spectralShape);
    p->setParameter("boostn", o.boostN);
    p->setParameter("usehartesyntax", o.harteSyntax ? 1.f : 0.f);
    p->setParameter("chromanormalize", (float)o.chromaNormalise);
}

int outputIndex(Plugin *p, const std::string &id)
{
    Plugin::OutputList outputs = p->getOutputDescriptors();
    for (int i = 0; i < (int)outputs.size(); ++i) {
        if (outputs[i].identifier == id) return i;
    }
    return -1;
}

// Wrap a plugin in the input-domain and buffering adapters, feed it the
// whole buffer, and return every feature it produced keyed by output.
// The adapter takes ownership of the plugin.
bool runPlugin(Plugin *plugin, const float *mono, size_t frames,
               float sampleRate, Plugin::FeatureSet &all, std::string *error)
{
    PluginInputDomainAdapter *ia = new PluginInputDomainAdapter(plugin);
    ia->setProcessTimestampMethod(PluginInputDomainAdapter::ShiftData);
    PluginBufferingAdapter *adapter = new PluginBufferingAdapter(ia);

    int blocksize = (int)adapter->getPreferredBlockSize();
    if (blocksize <= 0) blocksize = 16384;

    if (!adapter->initialise(1, blocksize, blocksize)) {
        if (error) *error = "failed to initialise " + plugin->getIdentifier();
        delete adapter;
        return false;
    }

    std::vector<float> block(blocksize);
    size_t pos = 0;
    while (pos < frames) {
        size_t count = frames - pos;
        if (count > (size_t)blocksize) count = blocksize;
        for (size_t i = 0; i < (size_t)blocksize; ++i) {
            block[i] = (i < count) ? mono[pos + i] : 0.f;
        }
        const float *ptr = &block[0];
        RealTime ts = RealTime::frame2RealTime((long)pos, (unsigned int)sampleRate);
        Plugin::FeatureSet fs = adapter->process(&ptr, ts);
        for (Plugin::FeatureSet::iterator i = fs.begin(); i != fs.end(); ++i) {
            all[i->first].insert(all[i->first].end(), i->second.begin(), i->second.end());
        }
        pos += count;
    }
    Plugin::FeatureSet fs = adapter->getRemainingFeatures();
    for (Plugin::FeatureSet::iterator i = fs.begin(); i != fs.end(); ++i) {
        all[i->first].insert(all[i->first].end(), i->second.begin(), i->second.end());
    }
    delete adapter;
    return true;
}

double seconds(const RealTime &t)
{
    return t.sec + t.nsec / 1e9;
}

std::string jsonEscape(const std::string &s)
{
    std::string out;
    for (size_t i = 0; i < s.size(); ++i) {
        char c = s[i];
        if (c == '"' || c == '\\') { out += '\\'; out += c; }
        else if ((unsigned char)c < 0x20) {
            char buf[8];
            snprintf(buf, sizeof buf, "\\u%04x", (unsigned)(unsigned char)c);
            out += buf;
        } else out += c;
    }
    return out;
}

std::string num(double v)
{
    char buf[64];
    snprintf(buf, sizeof buf, "%.6f", v);
    // trim trailing zeros for compactness
    std::string s(buf);
    size_t dot = s.find('.');
    if (dot != std::string::npos) {
        size_t end = s.size();
        while (end > dot + 2 && s[end - 1] == '0') --end;
        s.erase(end);
    }
    return s;
}

} // namespace

float estimateTuning(const float *mono, size_t frames, float sampleRate,
                     const Options &opts, std::string *error)
{
    Tuning *tuning = new Tuning(sampleRate);
    applyOptions(tuning, opts);
    int idx = outputIndex(tuning, "tuning");
    Plugin::FeatureSet all;
    if (!runPlugin(tuning, mono, frames, sampleRate, all, error)) return 0.f;
    if (idx < 0 || all[idx].empty()) {
        if (error) *error = "no tuning estimate produced";
        return 0.f;
    }
    // Label is "%0.1f Hz"
    return (float)atof(all[idx][0].label.c_str());
}

bool analyzeChords(const float *mono, size_t frames, float sampleRate,
                   const Options &opts, Result &out, std::string *error)
{
    Chordino *plugin = new Chordino(sampleRate);
    applyOptions(plugin, opts);
    int chordIdx = outputIndex(plugin, "simplechord");
    int notesIdx = outputIndex(plugin, "chordnotes");
    int changeIdx = outputIndex(plugin, "harmonicchange");
    if (chordIdx < 0 || notesIdx < 0 || changeIdx < 0) {
        if (error) *error = "Chordino outputs not found";
        delete plugin;
        return false;
    }

    Plugin::FeatureSet all;
    if (!runPlugin(plugin, mono, frames, sampleRate, all, error)) return false;

    out = Result();
    out.duration = (double)frames / sampleRate;

    const Plugin::FeatureList &chords = all[chordIdx];
    for (size_t i = 0; i < chords.size(); ++i) {
        Chord c;
        c.time = seconds(chords[i].timestamp);
        c.label = chords[i].label;
        c.duration = 0.0;
        out.chords.push_back(c);
    }
    // Chordino emits a final "N" at the last frame; keep it (it marks the
    // end) but give every chord a duration up to the next one.
    for (size_t i = 0; i < out.chords.size(); ++i) {
        double next = (i + 1 < out.chords.size()) ? out.chords[i + 1].time : out.duration;
        if (next < out.chords[i].time) next = out.chords[i].time;
        out.chords[i].duration = next - out.chords[i].time;
    }

    const Plugin::FeatureList &notes = all[notesIdx];
    for (size_t i = 0; i < notes.size(); ++i) {
        if (notes[i].values.empty()) continue;
        ChordNote n;
        n.time = seconds(notes[i].timestamp);
        n.duration = notes[i].hasDuration ? seconds(notes[i].duration) : 0.0;
        n.midi = (int)floor(notes[i].values[0] + 0.5f);
        out.notes.push_back(n);
    }

    const Plugin::FeatureList &change = all[changeIdx];
    for (size_t i = 0; i < change.size(); ++i) {
        out.frameTimes.push_back(seconds(change[i].timestamp));
        out.harmonicChange.push_back(change[i].values.empty() ? 0.0 : change[i].values[0]);
    }
    return true;
}

bool analyzeChroma(const float *mono, size_t frames, float sampleRate,
                   const Options &opts, std::vector<ChromaFrame> &out,
                   std::string *error)
{
    NNLSChroma *plugin = new NNLSChroma(sampleRate);
    applyOptions(plugin, opts);
    int chromaIdx = outputIndex(plugin, "chroma");
    int bassIdx = outputIndex(plugin, "basschroma");
    if (chromaIdx < 0 || bassIdx < 0) {
        if (error) *error = "NNLS Chroma outputs not found";
        delete plugin;
        return false;
    }
    Plugin::FeatureSet all;
    if (!runPlugin(plugin, mono, frames, sampleRate, all, error)) return false;

    const Plugin::FeatureList &treble = all[chromaIdx];
    const Plugin::FeatureList &bass = all[bassIdx];
    out.clear();
    for (size_t i = 0; i < treble.size(); ++i) {
        ChromaFrame f;
        f.time = seconds(treble[i].timestamp);
        for (int b = 0; b < 12; ++b) {
            f.chroma[b] = b < (int)treble[i].values.size() ? treble[i].values[b] : 0.f;
            f.bass[b] = (i < bass.size() && b < (int)bass[i].values.size()) ? bass[i].values[b] : 0.f;
        }
        out.push_back(f);
    }
    return true;
}

std::string toJson(const Result &r, float sampleRate, bool includeChange)
{
    std::ostringstream s;
    s << "{\n  \"sampleRate\": " << num(sampleRate)
      << ",\n  \"duration\": " << num(r.duration)
      << ",\n  \"chords\": [";
    for (size_t i = 0; i < r.chords.size(); ++i) {
        s << (i ? ",\n    " : "\n    ")
          << "{\"time\": " << num(r.chords[i].time)
          << ", \"duration\": " << num(r.chords[i].duration)
          << ", \"label\": \"" << jsonEscape(r.chords[i].label) << "\"}";
    }
    s << "\n  ],\n  \"notes\": [";
    for (size_t i = 0; i < r.notes.size(); ++i) {
        s << (i ? ",\n    " : "\n    ")
          << "{\"time\": " << num(r.notes[i].time)
          << ", \"duration\": " << num(r.notes[i].duration)
          << ", \"midi\": " << r.notes[i].midi << "}";
    }
    s << "\n  ]";
    if (includeChange) {
        s << ",\n  \"harmonicChange\": [";
        for (size_t i = 0; i < r.harmonicChange.size(); ++i) {
            s << (i ? ", " : "") << "[" << num(r.frameTimes[i]) << ", "
              << num(r.harmonicChange[i]) << "]";
        }
        s << "]";
    }
    s << "\n}\n";
    return s.str();
}

} // namespace chordino
