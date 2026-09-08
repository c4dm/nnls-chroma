// Checks that the in-process Chordino analysis finds the synthesized
// progression, then writes it as a WAV for the CLI tests.  GPL-2.0-or-later.
#include "cli/analyze.h"
#include "cli/audiofile.h"
#include "synth.h"

#include <cstdio>
#include <cstring>
#include <map>
#include <string>

using namespace chordino;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s\n", msg); ++failures; } else printf("ok: %s\n", msg); } while (0)

// The chord that covers most of [start, end) according to the result.
static std::string dominantChord(const Result &r, double start, double end)
{
    std::map<std::string, double> cover;
    for (size_t i = 0; i < r.chords.size(); ++i) {
        double a = r.chords[i].time, b = a + r.chords[i].duration;
        double lo = a > start ? a : start, hi = b < end ? b : end;
        if (hi > lo) cover[r.chords[i].label] += hi - lo;
    }
    std::string best;
    double bestCover = 0;
    for (std::map<std::string, double>::iterator i = cover.begin(); i != cover.end(); ++i) {
        if (i->second > bestCover) { bestCover = i->second; best = i->first; }
    }
    return best;
}

static void runAtRate(float rate, const char *wavPath)
{
    std::vector<TestChord> prog = testProgression();
    std::vector<float> audio = renderProgression(prog, rate);
    printf("-- %.0f Hz, %zu samples\n", rate, audio.size());

    Options opts;
    Result r;
    std::string err;
    CHECK(analyzeChords(&audio[0], audio.size(), rate, opts, r, &err), "analyzeChords succeeds");
    if (!err.empty()) printf("   error: %s\n", err.c_str());
    CHECK(r.chords.size() >= prog.size(), "at least one chord per segment");
    CHECK(r.duration > 7.9 && r.duration < 8.1, "duration is 8 s");

    double t = 0;
    for (size_t k = 0; k < prog.size(); ++k) {
        // skip the first 0.4 s of each segment: attack plus analysis latency
        std::string got = dominantChord(r, t + 0.4, t + prog[k].seconds);
        char msg[128];
        snprintf(msg, sizeof msg, "segment %zu is %s (got %s)", k, prog[k].label.c_str(), got.c_str());
        CHECK(got == prog[k].label, msg);
        t += prog[k].seconds;
    }
    CHECK(!r.notes.empty(), "chord notes produced");
    CHECK(r.harmonicChange.size() == r.frameTimes.size() && !r.frameTimes.empty(), "harmonic change per frame");

    // Harte syntax changes the minor label only
    Options harte;
    harte.harteSyntax = true;
    Result rh;
    CHECK(analyzeChords(&audio[0], audio.size(), rate, harte, rh, &err), "harte analysis succeeds");
    CHECK(dominantChord(rh, 6.4, 8.0) == "A:min", "harte label for A minor is A:min");

    // JSON contains every chord label
    std::string json = toJson(r, rate, true);
    CHECK(json.find("\"label\": \"Am\"") != std::string::npos, "json has Am");
    CHECK(json.find("\"harmonicChange\"") != std::string::npos, "json has harmonic change when asked");

    // Chroma: the C segment should have its strongest treble bin at C (index 3)
    std::vector<ChromaFrame> chroma;
    CHECK(analyzeChroma(&audio[0], audio.size(), rate, opts, chroma, &err), "analyzeChroma succeeds");
    int cWins = 0, cFrames = 0;
    for (size_t i = 0; i < chroma.size(); ++i) {
        if (chroma[i].time < 0.4 || chroma[i].time > 1.9) continue;
        ++cFrames;
        int best = 0;
        for (int b = 1; b < 12; ++b) if (chroma[i].chroma[b] > chroma[i].chroma[best]) best = b;
        if (best == 3) ++cWins;
    }
    CHECK(cFrames > 0 && cWins * 2 > cFrames, "C is the strongest chroma bin in the C segment");

    float hz = estimateTuning(&audio[0], audio.size(), rate, opts, &err);
    printf("   tuning estimate %.1f Hz\n", hz);
    CHECK(hz > 437.f && hz < 443.f, "tuning estimate near 440 Hz");

    if (wavPath) {
        CHECK(writeWavMono(wavPath, &audio[0], audio.size(), rate, &err), "wrote progression wav");
        AudioData back;
        CHECK(readAudioFile(wavPath, back, &err), "read progression wav back");
        CHECK(back.mono.size() == audio.size() && back.sampleRate == rate, "wav round trip preserves length and rate");
    }
}

int main(int argc, char **argv)
{
    const char *wavPath = argc > 1 ? argv[1] : NULL;
    runAtRate(44100.f, wavPath);
    runAtRate(48000.f, NULL);
    runAtRate(22050.f, NULL);
    printf("%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
