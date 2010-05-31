
// Remember to use a different guard symbol in each header!
#ifndef _NNLS_CHROMA_
#define _NNLS_CHROMA_

#include <vamp-sdk/Plugin.h>
#include <list>
// #include "FFT.h"
using namespace std;

using std::string;

class ChordTranscriberData;



class NNLSChroma : public Vamp::Plugin
{
public:
    NNLSChroma(float inputSampleRate);
    virtual ~NNLSChroma();

    string getIdentifier() const;
    string getName() const;
    string getDescription() const;
    string getMaker() const;
    int getPluginVersion() const;
    string getCopyright() const;

    InputDomain getInputDomain() const;
    size_t getPreferredBlockSize() const;
    size_t getPreferredStepSize() const;
    size_t getMinChannelCount() const;
    size_t getMaxChannelCount() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(string identifier) const;
    void setParameter(string identifier, float value);

    ProgramList getPrograms() const;
    string getCurrentProgram() const;
    void selectProgram(string name);

    OutputList getOutputDescriptors() const;

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    // plugin-specific data and methods go here
	int frameCount;
    FeatureList m_fl;
    size_t m_blockSize;
    size_t m_stepSize;
    int m_lengthOfNoteIndex;
	float m_meanTuning0;
	float m_meanTuning1;
	float m_meanTuning2;
    float m_localTuning0;
    float m_localTuning1;
    float m_localTuning2;
    float m_paling;
	float m_preset;
    vector<float> m_localTuning;
	vector<float> m_kernelValue;
	vector<int> m_kernelFftIndex;
	vector<int> m_kernelNoteIndex;
	float *m_dict;
    bool m_tuneLocal;
    int m_dictID;
    // list< vector< float > > *logfreqSpecList;
};



#endif
