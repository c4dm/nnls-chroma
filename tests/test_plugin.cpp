// Loads the built nnls-chroma plugin library through the Vamp host SDK
// (VAMP_PATH must point at the build's plugin directory), checks the three
// plugins are listed, and compares Chordino's chords with the in-process
// analysis.  GPL-2.0-or-later.
#include "cli/analyze.h"
#include "synth.h"

#include <vamp-hostsdk/PluginLoader.h>
#include <vamp-hostsdk/PluginInputDomainAdapter.h>
#include <vamp-hostsdk/PluginWrapper.h>

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace Vamp;
using namespace Vamp::HostExt;
using namespace chordino;

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { printf("FAIL: %s\n", msg); ++failures; } else printf("ok: %s\n", msg); } while (0)

int main()
{
    const char *vp = getenv("VAMP_PATH");
    printf("VAMP_PATH=%s\n", vp ? vp : "(unset)");

    PluginLoader *loader = PluginLoader::getInstance();
    std::vector<PluginLoader::PluginKey> keys = loader->listPlugins();
    bool haveChroma = false, haveChordino = false, haveTuning = false;
    for (size_t i = 0; i < keys.size(); ++i) {
        printf("  found %s\n", keys[i].c_str());
        if (keys[i] == "nnls-chroma:nnls-chroma") haveChroma = true;
        if (keys[i] == "nnls-chroma:chordino") haveChordino = true;
        if (keys[i] == "nnls-chroma:tuning") haveTuning = true;
    }
    CHECK(haveChroma, "nnls-chroma plugin listed");
    CHECK(haveChordino, "chordino plugin listed");
    CHECK(haveTuning, "tuning plugin listed");

    const float rate = 44100.f;
    std::vector<TestChord> prog = testProgression();
    std::vector<float> audio = renderProgression(prog, rate);

    Plugin *plugin = loader->loadPlugin("nnls-chroma:chordino", rate,
                                        PluginLoader::ADAPT_INPUT_DOMAIN | PluginLoader::ADAPT_BUFFER_SIZE);
    CHECK(plugin != NULL, "chordino loads from the library");
    if (!plugin) return 1;

    // Use the same framing as the in-process analysis (audio shifted by
    // half a block so the first frame is centred on time zero).
    PluginWrapper *wrapper = dynamic_cast<PluginWrapper *>(plugin);
    PluginInputDomainAdapter *ida = wrapper ? wrapper->getWrapper<PluginInputDomainAdapter>() : NULL;
    CHECK(ida != NULL, "input domain adapter present");
    if (ida) ida->setProcessTimestampMethod(PluginInputDomainAdapter::ShiftData);

    CHECK(plugin->getName() == "Chordino", "plugin name is Chordino");
    CHECK(plugin->getPluginVersion() >= 5, "plugin version >= 5");
    Plugin::ParameterList params = plugin->getParameterDescriptors();
    CHECK(params.size() == 7, "chordino has 7 parameters");

    int blocksize = (int)plugin->getPreferredBlockSize();
    CHECK(plugin->initialise(1, blocksize, blocksize), "plugin initialises");

    int chordIdx = -1;
    Plugin::OutputList outputs = plugin->getOutputDescriptors();
    for (size_t i = 0; i < outputs.size(); ++i) if (outputs[i].identifier == "simplechord") chordIdx = (int)i;
    CHECK(chordIdx >= 0, "simplechord output present");

    std::vector<float> block(blocksize);
    Plugin::FeatureList chords;
    size_t pos = 0;
    while (pos < audio.size()) {
        size_t count = audio.size() - pos;
        if (count > (size_t)blocksize) count = blocksize;
        for (int i = 0; i < blocksize; ++i) block[i] = (size_t)i < count ? audio[pos + i] : 0.f;
        const float *ptr = &block[0];
        Plugin::FeatureSet fs = plugin->process(&ptr, RealTime::frame2RealTime((long)pos, (unsigned)rate));
        chords.insert(chords.end(), fs[chordIdx].begin(), fs[chordIdx].end());
        pos += count;
    }
    Plugin::FeatureSet fs = plugin->getRemainingFeatures();
    chords.insert(chords.end(), fs[chordIdx].begin(), fs[chordIdx].end());
    delete plugin;

    Result inProcess;
    std::string err;
    CHECK(analyzeChords(&audio[0], audio.size(), rate, Options(), inProcess, &err), "in-process analysis succeeds");
    CHECK(chords.size() == inProcess.chords.size(), "same number of chords from library and in-process");
    bool same = chords.size() == inProcess.chords.size();
    for (size_t i = 0; same && i < chords.size(); ++i) {
        if (chords[i].label != inProcess.chords[i].label) same = false;
        double t = chords[i].timestamp.sec + chords[i].timestamp.nsec / 1e9;
        if (fabs(t - inProcess.chords[i].time) > 1e-6) same = false;
    }
    CHECK(same, "library and in-process chords match");
    for (size_t i = 0; i < chords.size(); ++i) {
        printf("  %s %s\n", chords[i].timestamp.toString().c_str(), chords[i].label.c_str());
    }

    printf("%s (%d failures)\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
