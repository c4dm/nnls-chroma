
#include "NNLSChroma.h"
#include <cmath>
#include <list>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdio>
// #include "cblas.h"
#include "chorddict.cpp"
using namespace std;

const float sinvalue = 0.866025404;
const float cosvalue = -0.5;
const float hammingwind[19] = {0.0082, 0.0110, 0.0191, 0.0316, 0.0470, 0.0633, 0.0786, 0.0911, 0.0992, 0.1020, 0.0992, 0.0911, 0.0786, 0.0633, 0.0470, 0.0316, 0.0191, 0.0110, 0.0082};
const float basswindow[] = {0.001769, 0.015848, 0.043608, 0.084265, 0.136670, 0.199341, 0.270509, 0.348162, 0.430105, 0.514023, 0.597545, 0.678311, 0.754038, 0.822586, 0.882019, 0.930656, 0.967124, 0.990393, 0.999803, 0.995091, 0.976388, 0.944223, 0.899505, 0.843498, 0.777785, 0.704222, 0.624888, 0.542025, 0.457975, 0.375112, 0.295778, 0.222215, 0.156502, 0.100495, 0.055777, 0.023612, 0.004909, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};
const float treblewindow[] = {0.000350, 0.003144, 0.008717, 0.017037, 0.028058, 0.041719, 0.057942, 0.076638, 0.097701, 0.121014, 0.146447, 0.173856, 0.203090, 0.233984, 0.266366, 0.300054, 0.334860, 0.370590, 0.407044, 0.444018, 0.481304, 0.518696, 0.555982, 0.592956, 0.629410, 0.665140, 0.699946, 0.733634, 0.766016, 0.796910, 0.826144, 0.853553, 0.878986, 0.902299, 0.923362, 0.942058, 0.958281, 0.971942, 0.982963, 0.991283, 0.996856, 0.999650, 0.999650, 0.996856, 0.991283, 0.982963, 0.971942, 0.958281, 0.942058, 0.923362, 0.902299, 0.878986, 0.853553, 0.826144, 0.796910, 0.766016, 0.733634, 0.699946, 0.665140, 0.629410, 0.592956, 0.555982, 0.518696, 0.481304, 0.444018, 0.407044, 0.370590, 0.334860, 0.300054, 0.266366, 0.233984, 0.203090, 0.173856, 0.146447, 0.121014, 0.097701, 0.076638, 0.057942, 0.041719, 0.028058, 0.017037, 0.008717, 0.003144, 0.000350};
const char* notenames[24] = {"A  (bass)","Bb (bass)","B  (bass)","C  (bass)","C# (bass)","D  (bass)","Eb (bass)","E  (bass)","F  (bass)","F# (bass)","G  (bass)","Ab (bass)",
"A","Bb","B","C","C#","D","Eb","E","F","F#","G","Ab"};
const vector<float> hw(hammingwind, hammingwind+19);
const int nNote = 256;

/** Special Convolution
special convolution is as long as the convolvee, i.e. the first argument. in the valid core part of the 
convolution it contains the usual convolution values, but the pads at the beginning (ending) have the same values
as the first (last) valid convolution bin.
**/

const bool debug_on = false;

vector<float> SpecialConvolution(vector<float> convolvee, vector<float> kernel)
{
    float s;
    int m, n;
    int lenConvolvee = convolvee.size();
    int lenKernel = kernel.size();

    vector<float> Z(256,0);
    assert(lenKernel % 2 != 0); // no exception handling !!!
    
    for (n = lenKernel - 1; n < lenConvolvee; n++) {
    	s=0.0;
    	for (m = 0; m < lenKernel; m++) {
            // cerr << "m = " << m << ", n = " << n << ", n-m = " << (n-m) << '\n';
            s += convolvee[n-m] * kernel[m];
            // if (debug_on) cerr << "--> s = " << s << '\n';
    	}
        // cerr << n - lenKernel/2 << endl;
        Z[n -lenKernel/2] = s;
    }
    
    // fill upper and lower pads
    for (n = 0; n < lenKernel/2; n++) Z[n] = Z[lenKernel/2];    
    for (n = lenConvolvee; n < lenConvolvee +lenKernel/2; n++) Z[n - lenKernel/2] = 
        Z[lenConvolvee - lenKernel/2 -  1];
    return Z;
}

// vector<float> FftBin2Frequency(vector<float> binnumbers, int fs, int blocksize)
// {
// 	vector<float> freq(binnumbers.size, 0.0);
// 	for (unsigned i = 0; i < binnumbers.size; ++i) {
// 		freq[i] = (binnumbers[i]-1.0) * fs * 1.0 / blocksize;	
// 	}
// 	return freq;
// }

float cospuls(float x, float centre, float width) 
{
	float recipwidth = 1.0/width;
	if (abs(x - centre) <= 0.5 * width) {
		return cos((x-centre)*2*M_PI*recipwidth)*.5+.5;
	}
	return 0.0;
}

float pitchCospuls(float x, float centre, int binsperoctave) 
{
	float warpedf = -binsperoctave * (log2(centre) - log2(x));
	float out = cospuls(warpedf, 0.0, 2.0);
	// now scale to correct for note density
	float c = log(2.0)/binsperoctave;
	if (x > 0) {
		out = out / (c * x);
	} else {
		out = 0;
	}
	return out;
}

bool logFreqMatrix(int fs, int blocksize, float *outmatrix) {
	
	int binspersemitone = 3; // this must be 3
	int minoctave = 0; // this must be 0
	int maxoctave = 7; // this must be 7
	int oversampling = 20;
	
	// linear frequency vector
	vector<float> fft_f;
	for (int i = 0; i < blocksize/2; ++i) {
		fft_f.push_back(i * (fs * 1.0 / blocksize));
	}
	float fft_width = fs * 2.0 / blocksize;
	
	// linear oversampled frequency vector
	vector<float> oversampled_f;
	for (unsigned int i = 0; i < oversampling * blocksize/2; ++i) {
		oversampled_f.push_back(i * ((fs * 1.0 / blocksize) / oversampling));
	}
	
	// pitch-spaced frequency vector
	int minMIDI = 21 + minoctave * 12 - 1; // this includes one additional semitone!
	int maxMIDI = 21 + maxoctave * 12; // this includes one additional semitone!
	vector<float> cq_f;
	float oob = 1.0/binspersemitone; // one over binspersemitone
	cq_f.push_back(440 * pow(2.0,0.083333 * (minMIDI-69))); // 0.083333 is approx 1/12
	cq_f.push_back(440 * pow(2.0,0.083333 * (minMIDI+oob-69)));
	for (int i = minMIDI + 1; i < maxMIDI; ++i) {
		for (int k = -1; k < 2; ++k)	 {
			cq_f.push_back(440 * pow(2.0,0.083333333333 * (i+oob*k-69)));
		}
	}
	cq_f.push_back(440 * pow(2.0,0.083333 * (minMIDI-oob-69)));
	cq_f.push_back(440 * pow(2.0,0.083333 * (maxMIDI-69)));

	// limit the linear vectors to the frequencies used
	float maxfreq = *cq_f.end() * pow(2.0,0.083333);
	while (*oversampled_f.end() > maxfreq) {
		oversampled_f.pop_back();
	}
	while (*fft_f.end() > maxfreq) {
		fft_f.pop_back();
	}
	
	int nFFT = fft_f.size();
	
	// for (int i = 0; i < 100; i++) {
	// 	cerr << oversampled_f[i] << " " << cospuls(oversampled_f[i],fft_f[1],fft_width) << endl;
	// }
	
	vector<float> fft_activation;
	for (int iOS = 0; iOS < 2 * oversampling; ++iOS) {
		float cosp = cospuls(oversampled_f[iOS],fft_f[1],fft_width);
		fft_activation.push_back(cosp);
		// cerr << cosp << endl;
	}
	
	float cq_activation;
	for (int iFFT = 1; iFFT < nFFT; ++iFFT) {
		// find frequency stretch where the oversampled vector can be non-zero (i.e. in a window of width fft_width around the current frequency)
		int curr_start = oversampling * iFFT - oversampling;
		int curr_end = oversampling * iFFT + oversampling; // don't know if I should add "+1" here
		// cerr << oversampled_f[curr_start] << " " << fft_f[iFFT] << " " << oversampled_f[curr_end] << endl;
		for (unsigned iCQ = 0; iCQ < cq_f.size(); ++iCQ) {
			outmatrix[iFFT + nFFT * iCQ] = 0;
			if (cq_f[iCQ] * pow(2.0, 0.084) + fft_width/2 > fft_f[iFFT] && cq_f[iCQ] * pow(2.0, -0.084 * 2) - fft_width/2 < fft_f[iFFT]) { // within a generous neighbourhood
				for (int iOS = curr_start; iOS < curr_end; ++iOS) {
					cq_activation = pitchCospuls(oversampled_f[iOS],cq_f[iCQ],binspersemitone*12);
					// cerr << oversampled_f[iOS] << " " << cq_f[iCQ] << " " << cq_activation << endl;
					outmatrix[iFFT + nFFT * iCQ] += cq_activation * fft_activation[iOS-curr_start];
				}				
				// if (iCQ == 1 || iCQ == 2) {
				// 	cerr << " " << outmatrix[iFFT + nFFT * iCQ] << endl;
				// }
			}
		}
	}
	return true;	
}


NNLSChroma::NNLSChroma(float inputSampleRate) :
  Plugin(inputSampleRate),
  m_fl(0),
  m_blockSize(0),
  m_stepSize(0),
  m_lengthOfNoteIndex(0),
  m_meanTuning0(0),
  m_meanTuning1(0),
  m_meanTuning2(0),
  m_localTuning0(0),
  m_localTuning1(0),
  m_localTuning2(0),
  m_paling(0),
  m_localTuning(0),
  m_kernelValue(0),
  m_kernelFftIndex(0),
  m_kernelNoteIndex(0),
  m_tuneLocal(false),
  m_dictID(0)
{
	if (debug_on) cerr << "--> NNLSChroma" << endl;
}


NNLSChroma::~NNLSChroma()
{
		if (debug_on) cerr << "--> ~NNLSChroma" << endl;
}

string
NNLSChroma::getIdentifier() const
{
	if (debug_on) cerr << "--> getIdentifier" << endl;
    return "nnls_chroma";
}

string
NNLSChroma::getName() const
{
		if (debug_on) cerr << "--> getName" << endl;
    return "NNLS Chroma";
}

string
NNLSChroma::getDescription() const
{
    // Return something helpful here!
	if (debug_on) cerr << "--> getDescription" << endl;
    return "";
}

string
NNLSChroma::getMaker() const
{
		if (debug_on) cerr << "--> getMaker" << endl;
    // Your name here
    return "Matthias Mauch";
}

int
NNLSChroma::getPluginVersion() const
{
		if (debug_on) cerr << "--> getPluginVersion" << endl;
    // Increment this each time you release a version that behaves
    // differently from the previous one
    return 1;
}

string
NNLSChroma::getCopyright() const
{
		if (debug_on) cerr << "--> getCopyright" << endl;
    // This function is not ideally named.  It does not necessarily
    // need to say who made the plugin -- getMaker does that -- but it
    // should indicate the terms under which it is distributed.  For
    // example, "Copyright (year). All Rights Reserved", or "GPL"
    return "Copyright (2010). All rights reserved.";
}

NNLSChroma::InputDomain
NNLSChroma::getInputDomain() const
{
		if (debug_on) cerr << "--> getInputDomain" << endl;
    return FrequencyDomain;
}

size_t
NNLSChroma::getPreferredBlockSize() const
{
		if (debug_on) cerr << "--> getPreferredBlockSize" << endl;
    return 16384; // 0 means "I can handle any block size"
}

size_t 
NNLSChroma::getPreferredStepSize() const
{
		if (debug_on) cerr << "--> getPreferredStepSize" << endl;
    return 2048; // 0 means "anything sensible"; in practice this
              // means the same as the block size for TimeDomain
              // plugins, or half of it for FrequencyDomain plugins
}

size_t
NNLSChroma::getMinChannelCount() const
{
	if (debug_on) cerr << "--> getMinChannelCount" << endl;
    return 1;
}

size_t
NNLSChroma::getMaxChannelCount() const
{
		if (debug_on) cerr << "--> getMaxChannelCount" << endl;
    return 1;
}

NNLSChroma::ParameterList
NNLSChroma::getParameterDescriptors() const
{
		if (debug_on) cerr << "--> getParameterDescriptors" << endl;
    ParameterList list;

    ParameterDescriptor d0;
    d0.identifier = "notedict";
    d0.name = "note dictionary";
    d0.description = "Notes in different note dictionaries differ by their spectral shapes.";
    d0.unit = "";
    d0.minValue = 0;
    d0.maxValue = 2;
    d0.defaultValue = 0;
    d0.isQuantized = true;
    d0.valueNames.push_back("s = 0.6");
    d0.valueNames.push_back("s = 0.9");
    d0.valueNames.push_back("s linearly spaced");
    d0.valueNames.push_back("no NNLS");
    d0.quantizeStep = 1.0;
    list.push_back(d0);

    ParameterDescriptor d1;
    d1.identifier = "tuningmode";
    d1.name = "tuning mode";
    d1.description = "Tuning can be performed locally or on the whole extraction area.";
    d1.unit = "";
    d1.minValue = 0;
    d1.maxValue = 1;
    d1.defaultValue = 1;
    d1.isQuantized = true;
    d1.valueNames.push_back("global tuning");
    d1.valueNames.push_back("local tuning");
    d1.quantizeStep = 1.0;
    list.push_back(d1);

    ParameterDescriptor d2;
    d2.identifier = "paling";
    d2.name = "spectral paling";
    d2.description = "Spectral paling: no paling - 0; whitening - 1.";
    d2.unit = "";
    d2.minValue = 0;
    d2.maxValue = 1;
    d2.defaultValue = 0.5;
    d2.isQuantized = false;
    // d1.valueNames.push_back("global tuning");
    // d1.valueNames.push_back("local tuning");
    // d1.quantizeStep = 0.1;
    list.push_back(d2);

    return list;
}

float
NNLSChroma::getParameter(string identifier) const
{
		if (debug_on) cerr << "--> getParameter" << endl;
    if (identifier == "notedict") {
        return m_dictID; 
    }
    
    if (identifier == "paling") {
        return m_paling; 
    }
    
    if (identifier == "tuningmode") {
        if (m_tuneLocal) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    return 0;
    
}

void
NNLSChroma::setParameter(string identifier, float value) 
{
		if (debug_on) cerr << "--> setParameter" << endl;
    if (identifier == "notedict") {
        m_dictID = (int) value;
    }
    
    if (identifier == "paling") {
        m_paling = value;
    }
    
    if (identifier == "tuningmode") {
        m_tuneLocal = (value > 0) ? true : false;
        // cerr << "m_tuneLocal :" << m_tuneLocal << endl;
    }
}

NNLSChroma::ProgramList
NNLSChroma::getPrograms() const
{
		if (debug_on) cerr << "--> getPrograms" << endl;
    ProgramList list;

    // If you have no programs, return an empty list (or simply don't
    // implement this function or getCurrentProgram/selectProgram)

    return list;
}

string
NNLSChroma::getCurrentProgram() const
{
		if (debug_on) cerr << "--> getCurrentProgram" << endl;
    return ""; // no programs
}

void
NNLSChroma::selectProgram(string name)
{
		if (debug_on) cerr << "--> selectProgram" << endl;
}


NNLSChroma::OutputList
NNLSChroma::getOutputDescriptors() const
{
		if (debug_on) cerr << "--> getOutputDescriptors" << endl;
    OutputList list;
    
    // Make chroma names for the binNames property
    vector<string> chromanames;
    vector<string> bothchromanames;
    for (int iNote = 0; iNote < 24; iNote++) {
        bothchromanames.push_back(notenames[iNote]);
        if (iNote < 12) {
            chromanames.push_back(notenames[iNote]);
        }
    }
    
	int nNote = 84;

    // See OutputDescriptor documentation for the possibilities here.
    // Every plugin must have at least one output.

    OutputDescriptor d0;
    d0.identifier = "tuning";
    d0.name = "Tuning";
    d0.description = "The concert pitch.";
    d0.unit = "Hz";
    d0.hasFixedBinCount = true;
    d0.binCount = 0;
    d0.hasKnownExtents = true;
	d0.minValue = 427.47;
	d0.maxValue = 452.89;
    d0.isQuantized = false;
    d0.sampleType = OutputDescriptor::VariableSampleRate;
    d0.hasDuration = false;
    list.push_back(d0);
	
	OutputDescriptor d1;
    d1.identifier = "logfreqspec";
    d1.name = "Log-Frequency Spectrum";
    d1.description = "A Log-Frequency Spectrum (constant Q) that is obtained by cosine filter mapping.";
    d1.unit = "";
    d1.hasFixedBinCount = true;
    d1.binCount = nNote;
    d1.hasKnownExtents = false;
    d1.isQuantized = false;
    d1.sampleType = OutputDescriptor::FixedSampleRate;
    d1.hasDuration = false;
    d1.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d1);

	OutputDescriptor d2;
    d2.identifier = "tunedlogfreqspec";
    d2.name = "Tuned Log-Frequency Spectrum";
    d2.description = "A Log-Frequency Spectrum (constant Q) that is obtained by cosine filter mapping, then its tuned using the estimated tuning frequency.";
    d2.unit = "";
    d2.hasFixedBinCount = true;
    d2.binCount = 256;
    d2.hasKnownExtents = false;
    d2.isQuantized = false;
    d2.sampleType = OutputDescriptor::FixedSampleRate;
    d2.hasDuration = false;
    d2.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d2);
    
    OutputDescriptor d3;
    d3.identifier = "semitonespectrum";
    d3.name = "Semitone Spectrum";
    d3.description = "A semitone-spaced log-frequency spectrum derived from the third-of-a-semitone-spaced tuned log-frequency spectrum.";
    d3.unit = "";
    d3.hasFixedBinCount = true;
    d3.binCount = 84;
    d3.hasKnownExtents = false;
    d3.isQuantized = false;
    d3.sampleType = OutputDescriptor::FixedSampleRate;
    d3.hasDuration = false;
    d3.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d3);
    
    OutputDescriptor d4;
    d4.identifier = "chroma";
    d4.name = "Chromagram";
    d4.description = "Tuning-adjusted chromagram from NNLS soft transcription, with an emphasis on the medium note range.";
    d4.unit = "";
    d4.hasFixedBinCount = true;
    d4.binCount = 12;
    d4.binNames = chromanames;
    d4.hasKnownExtents = false;
    d4.isQuantized = false;
    d4.sampleType = OutputDescriptor::FixedSampleRate;
    d4.hasDuration = false;
    d4.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d4);
    
    OutputDescriptor d5;
    d5.identifier = "basschroma";
    d5.name = "Bass Chromagram";
    d5.description = "Tuning-adjusted bass chromagram from NNLS soft transcription, with an emphasis on the bass note range.";
    d5.unit = "";
    d5.hasFixedBinCount = true;
    d5.binCount = 12;
    d5.binNames = chromanames;
    d5.hasKnownExtents = false;
    d5.isQuantized = false;
    d5.sampleType = OutputDescriptor::FixedSampleRate;
    d5.hasDuration = false;
    d5.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d5);
    
    OutputDescriptor d6;
    d6.identifier = "bothchroma";
    d6.name = "Chromagram and Bass Chromagram";
    d6.description = "Tuning-adjusted chromagram and bass chromagram (stacked on top of each other) from NNLS soft transcription.";
    d6.unit = "";
    d6.hasFixedBinCount = true;
    d6.binCount = 24;
    d6.binNames = bothchromanames;
    d6.hasKnownExtents = false;
    d6.isQuantized = false;
    d6.sampleType = OutputDescriptor::FixedSampleRate;
    d6.hasDuration = false;
    d6.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d6);
    
    OutputDescriptor d7;
    d7.identifier = "simplechord";
    d7.name = "Simple Chord Estimate";
    d7.description = "A simple chord estimate based on the inner product of chord templates with the smoothed chroma.";
    d7.unit = "";
    d7.hasFixedBinCount = true;
    d7.binCount = 0;
    d7.hasKnownExtents = false;
    d7.isQuantized = false;
    d7.sampleType = OutputDescriptor::VariableSampleRate;
    d7.hasDuration = false;
    d7.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d7);
    
    OutputDescriptor d8;
    d8.identifier = "inconsistency";
    d8.name = "Harmonic inconsistency value";
    d8.description = "Harmonic inconsistency. Indicates music if low, non-music or speech when high.";
    d8.unit = "";
    d8.hasFixedBinCount = true;
    d8.binCount = 1;
    d8.hasKnownExtents = false;
    d8.isQuantized = false;
    d8.sampleType = OutputDescriptor::FixedSampleRate;
    d8.hasDuration = false;
    d8.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d8);
    
    OutputDescriptor d9;
    d9.identifier = "inconsistencysegment";
    d9.name = "Harmonic inconsistency segmenter";
    d9.description = "Segments the audio based on the harmonic inconsistency value into speech and music.";
    d9.unit = "";
    d9.hasFixedBinCount = true;
    d9.binCount = 0;
    d9.hasKnownExtents = true;
    d9.minValue = 0.1;
	d9.maxValue = 0.9;
    d9.isQuantized = false;
    d9.sampleType = OutputDescriptor::VariableSampleRate;
    d9.hasDuration = false;
    d9.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d9);

    OutputDescriptor d10;
    d10.identifier = "localtuning";
    d10.name = "Local tuning";
    d10.description = "";
    d10.unit = "Hz";
    d10.hasFixedBinCount = true;
    d10.binCount = 1;
    d10.hasKnownExtents = true;
	d10.minValue = 427.47;
	d10.maxValue = 452.89;
    d10.isQuantized = false;
    d10.sampleType = OutputDescriptor::OneSamplePerStep;
    d10.hasDuration = false;
    d10.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d10);

    return list;
}


bool
NNLSChroma::initialise(size_t channels, size_t stepSize, size_t blockSize)
{	
		if (debug_on) cerr << "--> initialise";
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    frameCount = 0;
	int tempn = 256 * m_blockSize/2;
	cerr << tempn;
	float *tempkernel = new float[tempn];
	
	// vector<float> m_kernelValue;
	// vector<int> m_kernelFftIndex;
	// vector<int> m_kernelNoteIndex;
	
	
	logFreqMatrix(m_inputSampleRate, m_blockSize, tempkernel);

	for (unsigned iNote = 0; iNote < nNote; ++iNote) { // I don't know if this is wise: manually making a sparse matrix
			for (unsigned iFFT = 0; iFFT < blockSize/2; ++iFFT) {
				if (tempkernel[iFFT + blockSize/2 * iNote] > 0) {
					m_kernelValue.push_back(tempkernel[iFFT + blockSize/2 * iNote]);
					m_kernelFftIndex.push_back(iFFT);
					m_kernelNoteIndex.push_back(iNote);				
				}
			}
		}
	delete tempkernel;
	// int count = 0;
	// for (vector<float>::iterator it = m_kernelValue.begin(); it != m_kernelValue.end(); ++it) {
	// 	cerr << m_kernelFftIndex[count] << " -- " << m_kernelNoteIndex[count] << " -- " << m_kernelValue[count] << endl;
	// 	count++;
	// }
    return true;
}

void
NNLSChroma::reset()
{
		if (debug_on) cerr << "--> reset";
    // Clear buffers, reset stored values, etc
	    frameCount = 0;
        m_dictID = 0;
}

NNLSChroma::FeatureSet
NNLSChroma::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{   
		if (debug_on) cerr << "--> process" << endl;
	// int nNote = 84; // TODO: this should be globally set and/or depend on the kernel matrix
    
	frameCount++;   
	float *magnitude = new float[m_blockSize/2];
	
	Feature f10; // local tuning
	
	const float *fbuf = inputBuffers[0];	
	
	// make magnitude
	for (size_t iBin = 0; iBin < m_blockSize/2; iBin++) {
		magnitude[iBin] = sqrt(fbuf[2 * iBin] * fbuf[2 * iBin] + 
			fbuf[2 * iBin + 1] * fbuf[2 * iBin + 1]);
	}
	
	
	// note magnitude mapping using pre-calculated matrix
	float *nm  = new float[nNote]; // note magnitude
	for (size_t iNote = 0; iNote < nNote; iNote++) {
		nm[iNote] = 0; // initialise as 0
	}
	int binCount = 0;
	for (vector<float>::iterator it = m_kernelValue.begin(); it != m_kernelValue.end(); ++it) {
		// cerr << ".";
		nm[m_kernelNoteIndex[binCount]] += magnitude[m_kernelFftIndex[binCount]] * m_kernelValue[binCount];	
		binCount++;	
	}
	cerr << nm[20];
	cerr << endl;
	
	
    float one_over_N = 1.0/frameCount;
    // update means of complex tuning variables
    m_meanTuning0 *= float(frameCount-1)*one_over_N;
    m_meanTuning1 *= float(frameCount-1)*one_over_N;
    m_meanTuning2 *= float(frameCount-1)*one_over_N;
	
    for (int iTone = 0; iTone < 160; iTone = iTone + 3) {
        m_meanTuning0 += nm[iTone + 0]*one_over_N;
    	m_meanTuning1 += nm[iTone + 1]*one_over_N;
    	m_meanTuning2 += nm[iTone + 2]*one_over_N;
        m_localTuning0 *= 0.99994; m_localTuning0 += nm[iTone + 0];
        m_localTuning1 *= 0.99994; m_localTuning1 += nm[iTone + 1];
        m_localTuning2 *= 0.99994; m_localTuning2 += nm[iTone + 2];
    }
	
    // if (m_tuneLocal) {
    	// local tuning
        float localTuningImag = sinvalue * m_localTuning1 - sinvalue * m_localTuning2;
        float localTuningReal = m_localTuning0 + cosvalue * m_localTuning1 + cosvalue * m_localTuning2;
        float normalisedtuning = atan2(localTuningImag, localTuningReal)/(2*M_PI);
        m_localTuning.push_back(normalisedtuning);
        float tuning440 = 440 * pow(2,normalisedtuning/12);
        f10.values.push_back(tuning440);
    // }
    
	Feature f1; // logfreqspec
	f1.hasTimestamp = true;
    f1.timestamp = timestamp;
	for (size_t iNote = 0; iNote < nNote; iNote++) {
		f1.values.push_back(nm[iNote]);
	}
	
	FeatureSet fs;
	fs[1].push_back(f1);
    fs[10].push_back(f10);

    // deletes
    delete[] magnitude;
    delete[] nm;

    m_fl.push_back(f1); // remember note magnitude for getRemainingFeatures
	return fs;	
}

NNLSChroma::FeatureSet
NNLSChroma::getRemainingFeatures()
{
		if (debug_on) cerr << "--> getRemainingFeatures" << endl;
    FeatureSet fsOut;
	// 
	// /**  Calculate Tuning
	// calculate tuning from (using the angle of the complex number defined by the 
	// cumulative mean real and imag values)
	// **/
	// float meanTuningImag = sinvalue * m_meanTuning1 - sinvalue * m_meanTuning2;
	//     float meanTuningReal = m_meanTuning0 + cosvalue * m_meanTuning1 + cosvalue * m_meanTuning2;
	//     float cumulativetuning = 440 * pow(2,atan2(meanTuningImag, meanTuningReal)/(24*M_PI));
	//     float normalisedtuning = atan2(meanTuningImag, meanTuningReal)/(2*M_PI);
	//     int intShift = floor(normalisedtuning * 3);
	//     float intFactor = normalisedtuning * 3 - intShift; // intFactor is a really bad name for this
	//     
	//     char buffer0 [50];
	// 
	//     sprintf(buffer0, "estimated tuning: %0.1f Hz", cumulativetuning);
	//     
	//     // cerr << "normalisedtuning: " << normalisedtuning << '\n';
	//     
	//     // push tuning to FeatureSet fsOut
	// Feature f0; // tuning
	// f0.hasTimestamp = true;
	//     f0.timestamp = Vamp::RealTime::frame2RealTime(0, lrintf(m_inputSampleRate));;
	//     f0.label = buffer0;
	//     fsOut[0].push_back(f0);  
	//     
	//     /** Tune Log-Frequency Spectrogram
	//     calculate a tuned log-frequency spectrogram (f2): use the tuning estimated above (kinda f0) to 
	//     perform linear interpolation on the existing log-frequency spectrogram (kinda f1).
	//     **/    
	// 
	//     float tempValue = 0;
	//     float dbThreshold = 0; // relative to the background spectrum
	//     float thresh = pow(10,dbThreshold/20);
	//     // cerr << "tune local ? " << m_tuneLocal << endl;
	//     int count = 0;
	// 
	//     for (FeatureList::iterator i = m_fl.begin(); i != m_fl.end(); ++i) {
	//         Feature f1 = *i;
	//         Feature f2; // tuned log-frequency spectrum
	//         f2.hasTimestamp = true;
	//         f2.timestamp = f1.timestamp;
	//         f2.values.push_back(0.0); f2.values.push_back(0.0); // set lower edge to zero
	// 
	//         if (m_tuneLocal) {
	//             intShift = floor(m_localTuning[count] * 3);
	//             intFactor = m_localTuning[count] * 3 - intShift; // intFactor is a really bad name for this
	//         }
	//         
	//         // cerr << intShift << " " << intFactor << endl;
	//         
	//         for (int k = 2; k < f1.values.size() - 3; ++k) { // interpolate all inner bins
	//             tempValue = f1.values[k + intShift] * (1-intFactor) + f1.values[k+intShift+1] * intFactor;
	//             f2.values.push_back(tempValue);
	//         }
	//         
	//         f2.values.push_back(0.0); f2.values.push_back(0.0); f2.values.push_back(0.0); // upper edge
	//         vector<float> runningmean = SpecialConvolution(f2.values,hw);
	//         vector<float> runningstd;
	//         for (int i = 0; i < 256; i++) { // first step: squared values into vector (variance)
	//             runningstd.push_back((f2.values[i] - runningmean[i]) * (f2.values[i] - runningmean[i]));
	//         }
	//         runningstd = SpecialConvolution(runningstd,hw); // second step convolve
	//         for (int i = 0; i < 256; i++) { 
	//             runningstd[i] = sqrt(runningstd[i]); // square root to finally have running std
	//             if (runningstd[i] > 0) {
	//                 f2.values[i] = (f2.values[i] / runningmean[i]) > thresh ? 
	//                     (f2.values[i] - runningmean[i]) / pow(runningstd[i],m_paling) : 0;
	//             }
	//             if (f2.values[i] < 0) {
	//                 cerr << "ERROR: negative value in logfreq spectrum" << endl;
	//             }
	//         }
	//         fsOut[2].push_back(f2);
	//         count++;
	//     }
	//     
	//     /** Semitone spectrum and chromagrams
	//     Semitone-spaced log-frequency spectrum derived from the tuned log-freq spectrum above. the spectrum
	//     is inferred using a non-negative least squares algorithm.
	//     Three different kinds of chromagram are calculated, "treble", "bass", and "both" (which means 
	//     bass and treble stacked onto each other).
	//     **/
	//     // taucs_ccs_matrix* A_original_ordering = taucs_construct_sorted_ccs_matrix(nnlsdict06, nnls_m, nnls_n);
	//     
	//     vector<vector<float> > chordogram;
	//     vector<float> oldchroma = vector<float>(12,0);
	//     vector<float> oldbasschroma = vector<float>(12,0);
	//     count = 0;
	// 
	//     for (FeatureList::iterator it = fsOut[2].begin(); it != fsOut[2].end(); ++it) {
	//         Feature f2 = *it; // logfreq spectrum
	//         Feature f3; // semitone spectrum
	//         Feature f4; // treble chromagram
	//         Feature f5; // bass chromagram
	//         Feature f6; // treble and bass chromagram
	// 
	//         f3.hasTimestamp = true;
	//         f3.timestamp = f2.timestamp;
	//         
	//         f4.hasTimestamp = true;
	//         f4.timestamp = f2.timestamp;
	//         
	//         f5.hasTimestamp = true;
	//         f5.timestamp = f2.timestamp;
	//         
	//         f6.hasTimestamp = true;
	//         f6.timestamp = f2.timestamp;
	//         
	// 	float b[256];
	// 
	//         bool some_b_greater_zero = false;
	//         for (int i = 0; i < 256; i++) {
	//             b[i] = f2.values[i];
	//             if (b[i] > 0) {
	//                 some_b_greater_zero = true;
	//             }            
	//         }
	//     
	//         // here's where the non-negative least squares algorithm calculates the note activation x
	// 
	//         vector<float> chroma = vector<float>(12, 0);
	//         vector<float> basschroma = vector<float>(12, 0);
	//         if (some_b_greater_zero) {
	//         }
	// 
	//         f4.values = chroma;
	//         f5.values = basschroma;
	//         chroma.insert(chroma.begin(), basschroma.begin(), basschroma.end()); // just stack the both chromas 
	//         f6.values = chroma; 
	//         
	//         // local chord estimation
	//         vector<float> currentChordSalience;
	//         float tempchordvalue = 0;
	//         float sumchordvalue = 0;
	//         int nChord = nChorddict / 24;
	//         for (int iChord = 0; iChord < nChord; iChord++) {
	//             tempchordvalue = 0;
	//             for (int iBin = 0; iBin < 12; iBin++) {
	//                 tempchordvalue += chorddict[24 * iChord + iBin] * chroma[iBin];
	//             }
	//             for (int iBin = 12; iBin < 24; iBin++) {
	//                 tempchordvalue += chorddict[24 * iChord + iBin] * chroma[iBin];
	//             }
	//             sumchordvalue+=tempchordvalue;
	//             currentChordSalience.push_back(tempchordvalue);
	//         }
	//         for (int iChord = 0; iChord < nChord; iChord++) {
	//             currentChordSalience[iChord] /= sumchordvalue;
	//         }
	//         chordogram.push_back(currentChordSalience);
	//         
	//         fsOut[3].push_back(f3);
	//         fsOut[4].push_back(f4);
	//         fsOut[5].push_back(f5);
	//         fsOut[6].push_back(f6);
	// // if (x) free(x);
	//         // delete[] b;
	//         count++;
	//     }
	//     // cerr << m_stepSize << endl<< endl;
	//     count = 0;
	//     int kernelwidth = (49 * 2048) / m_stepSize;
	//     int nChord = nChorddict / 24;
	//     int musicitykernelwidth = (50 * 2048) / m_stepSize;
	//     
	//     /* Simple chord estimation
	//     I just take the local chord estimates ("currentChordSalience") and average them over time, then
	//     take the maximum. Very simple, don't do this at home...
	//     */
	//     vector<int> chordSequence;
	//     for (FeatureList::iterator it = fsOut[6].begin(); it != fsOut[6].end(); ++it) {
	//     
	//         int startIndex = max(count - kernelwidth/2 + 1,0);
	//         int endIndex = min(int(chordogram.size()), startIndex + kernelwidth - 1 + 1);
	//         vector<float> temp = vector<float>(nChord,0);
	//         for (int iChord = 0; iChord < nChord; iChord++) {
	//             float val = 0;
	//             for (int i = startIndex; i < endIndex; i++) {
	//                 val += chordogram[i][iChord] * 
	//                     (kernelwidth - abs(i - startIndex - kernelwidth * 0.5)); // weigthed sum (triangular window)
	//             }
	//             temp[iChord] = val; // sum
	//         }
	//         
	//         // get maximum for "chord estimate"
	//         
	//         float bestChordValue = 0;
	//         int bestChordIndex = nChord-1; // "no chord" is default
	//         for (int iChord = 0; iChord < nChord; iChord++) {
	//             if (temp[iChord] > bestChordValue) {
	//                 bestChordValue = temp[iChord];
	//                 bestChordIndex = iChord;
	//             }
	//         }
	//         // cerr << bestChordIndex << endl;
	//         chordSequence.push_back(bestChordIndex);
	//         count++;
	//     }
	//     // mode filter on chordSequence
	//     count = 0;
	//     int oldChordIndex = -1;
	//     for (FeatureList::iterator it = fsOut[6].begin(); it != fsOut[6].end(); ++it) {
	//         Feature f6 = *it;
	//         Feature f7; // chord estimate
	//         
	//         f7.hasTimestamp = true;
	//         f7.timestamp = f6.timestamp;
	//         
	//         vector<int> chordCount = vector<int>(121,0);
	//         
	//         int maxChordCount = 0;
	//         int maxChordIndex = 120;
	//         int startIndex = max(count - kernelwidth/2,0);
	//         int endIndex = min(int(chordogram.size()), startIndex + kernelwidth - 1);
	//         for (int i = startIndex; i < endIndex; i++) {
	//             chordCount[chordSequence[i]]++;
	//             if (chordCount[chordSequence[i]] > maxChordCount) {
	//                 maxChordCount++;
	//                 maxChordIndex = chordSequence[i];
	//             }
	//         }
	//         if (oldChordIndex != maxChordIndex) {
	//             oldChordIndex = maxChordIndex;
	// 
	//             char buffer1 [50];
	//             if (maxChordIndex < nChord - 1) {
	//                 sprintf(buffer1, "%s%s", notenames[maxChordIndex % 12 + 12], chordtypes[maxChordIndex]);
	//             } else {
	//                 sprintf(buffer1, "N");
	//             }
	//             f7.label = buffer1;
	//             fsOut[7].push_back(f7);
	//         }
	//         count++;
	//     }
	//     // musicity
	//     count = 0;
	//     int oldlabeltype = 0; // start value is 0, music is 1, speech is 2
	//     vector<float> musicityValue; 
	//     for (FeatureList::iterator it = fsOut[4].begin(); it != fsOut[4].end(); ++it) {
	//         Feature f4 = *it;
	//         
	//         int startIndex = max(count - musicitykernelwidth/2,0);
	//         int endIndex = min(int(chordogram.size()), startIndex + musicitykernelwidth - 1);
	//         float chromasum = 0;
	//         float diffsum = 0;
	//         for (int k = 0; k < 12; k++) {
	//             for (int i = startIndex + 1; i < endIndex; i++) {
	//                 chromasum += pow(fsOut[4][i].values[k],2);
	//                 diffsum += abs(fsOut[4][i-1].values[k] - fsOut[4][i].values[k]);
	//             }
	//         }
	//         diffsum /= chromasum;
	//         musicityValue.push_back(diffsum);        
	//         count++;
	//     }
	//     
	//     float musicityThreshold = 0.44;
	//     if (m_stepSize == 4096) {
	//         musicityThreshold = 0.74;
	//     }
	//     if (m_stepSize == 4410) {
	//         musicityThreshold = 0.77;
	//     }
	//     
	//     count = 0;
	//     for (FeatureList::iterator it = fsOut[4].begin(); it != fsOut[4].end(); ++it) {
	//         Feature f4 = *it;
	//         Feature f8; // musicity
	//         Feature f9; // musicity segmenter
	//         
	//         f8.hasTimestamp = true;
	//         f8.timestamp = f4.timestamp;
	//         f9.hasTimestamp = true;
	//         f9.timestamp = f4.timestamp;    
	//         
	//         int startIndex = max(count - musicitykernelwidth/2,0);
	//         int endIndex = min(int(chordogram.size()), startIndex + musicitykernelwidth - 1);
	//         int musicityCount = 0;
	//         for (int i = startIndex; i <= endIndex; i++) {
	//             if (musicityValue[i] > musicityThreshold) musicityCount++;
	//         }
	//         bool isSpeech = (2 * musicityCount > endIndex - startIndex + 1); 
	//         
	//         if (isSpeech) {
	//             if (oldlabeltype != 2) {
	//                 f9.label = "Speech";
	//                 fsOut[9].push_back(f9);
	//                 oldlabeltype = 2;
	//             }
	//         } else {
	//             if (oldlabeltype != 1) {
	//                 f9.label = "Music";
	//                 fsOut[9].push_back(f9);
	//                 oldlabeltype = 1;
	//             }
	//         }
	//         f8.values.push_back(musicityValue[count]);
	//         fsOut[8].push_back(f8);
	//         count++;
	//      }
     return fsOut;     

}

