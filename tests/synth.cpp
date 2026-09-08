// Test-only chord progression synthesizer.  GPL-2.0-or-later.
#include "synth.h"
#include <cmath>

std::vector<TestChord> testProgression()
{
    std::vector<TestChord> p;
    int c[]  = { 36, 48, 52, 55, 60 };      // C2 C3 E3 G3 C4
    int f[]  = { 41, 53, 57, 60, 65 };      // F2 F3 A3 C4 F4
    int g[]  = { 43, 55, 59, 62, 67 };      // G2 G3 B3 D4 G4
    int am[] = { 45, 57, 60, 64, 69 };      // A2 A3 C4 E4 A4
    TestChord t;
    t.seconds = 2.0;
    t.label = "C";  t.midi.assign(c, c + 5);   p.push_back(t);
    t.label = "F";  t.midi.assign(f, f + 5);   p.push_back(t);
    t.label = "G";  t.midi.assign(g, g + 5);   p.push_back(t);
    t.label = "Am"; t.midi.assign(am, am + 5); p.push_back(t);
    return p;
}

std::vector<float> renderProgression(const std::vector<TestChord> &prog, float sampleRate)
{
    std::vector<float> out;
    const double pi = 3.14159265358979323846;
    for (size_t k = 0; k < prog.size(); ++k) {
        size_t n = (size_t)(prog[k].seconds * sampleRate);
        size_t start = out.size();
        out.resize(start + n, 0.f);
        for (size_t j = 0; j < prog[k].midi.size(); ++j) {
            double f0 = 440.0 * pow(2.0, (prog[k].midi[j] - 69) / 12.0);
            for (size_t i = 0; i < n; ++i) {
                double t = i / (double)sampleRate;
                // a few decaying harmonics so it looks like a plucked note
                double env = exp(-t * 0.8) * (1.0 - exp(-t * 200.0));
                double s = 0.0;
                for (int h = 1; h <= 6; ++h) {
                    if (f0 * h > sampleRate * 0.45) break;
                    s += sin(2.0 * pi * f0 * h * t) * pow(0.6, h - 1);
                }
                out[start + i] += (float)(s * env * 0.12);
            }
        }
    }
    return out;
}
