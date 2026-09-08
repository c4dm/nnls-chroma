// Test-only chord progression synthesizer.  GPL-2.0-or-later.
#ifndef CHORDINO_TEST_SYNTH_H
#define CHORDINO_TEST_SYNTH_H
#include <string>
#include <vector>

struct TestChord {
    std::string label;      // expected Chordino label
    std::vector<int> midi;  // notes to play
    double seconds;
};

// The progression every test uses: C, F, G, Am, two seconds each.
std::vector<TestChord> testProgression();

// Render the progression to mono float samples.
std::vector<float> renderProgression(const std::vector<TestChord> &prog, float sampleRate);

#endif
