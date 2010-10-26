/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  NNLS-Chroma / Chordino

  Audio feature extraction plugins for chromagram and chord
  estimation.

  Centre for Digital Music, Queen Mary University of London.
  This file copyright 2008-2010 Matthias Mauch and QMUL.
    
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of the
  License, or (at your option) any later version.  See the file
  COPYING included with this distribution for more information.
*/

#include "NNLSBase.h"

#include "chromamethods.h"

#include <cstdlib>
#include <fstream>
#include <cmath>

#include <algorithm>

const bool debug_on = false;

const vector<float> hw(hammingwind, hammingwind+19);

NNLSBase::NNLSBase(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_logSpectrum(0),
    m_blockSize(0),
    m_stepSize(0),
    m_lengthOfNoteIndex(0),
    m_meanTuning0(0),
    m_meanTuning1(0),
    m_meanTuning2(0),
    m_localTuning0(0),
    m_localTuning1(0),
    m_localTuning2(0),
    m_whitening(1.0),
    m_preset(0.0),
    m_localTuning(0),
    m_kernelValue(0),
    m_kernelFftIndex(0),
    m_kernelNoteIndex(0),
    m_dict(0),
    m_tuneLocal(0),
    m_chorddict(0),
    m_chordnames(0),
    m_doNormalizeChroma(0),
    m_rollon(0),
	m_s(0.7),
	m_useNNLS(1),
	m_useHMM(1)
{
    if (debug_on) cerr << "--> NNLSBase" << endl;

    // make the *note* dictionary matrix
    m_dict = new float[nNote * 84];
    for (unsigned i = 0; i < nNote * 84; ++i) m_dict[i] = 0.0;
    dictionaryMatrix(m_dict, 0.7);
	
    // get the *chord* dictionary from file (if the file exists)
    m_chordnames = chordDictionary(&m_chorddict);
}


NNLSBase::~NNLSBase()
{
    if (debug_on) cerr << "--> ~NNLSBase" << endl;
    delete [] m_dict;
}

string
NNLSBase::getMaker() const
{
    if (debug_on) cerr << "--> getMaker" << endl;
    // Your name here
    return "Matthias Mauch";
}

int
NNLSBase::getPluginVersion() const
{
    if (debug_on) cerr << "--> getPluginVersion" << endl;
    // Increment this each time you release a version that behaves
    // differently from the previous one
    return 1;
}

string
NNLSBase::getCopyright() const
{
    if (debug_on) cerr << "--> getCopyright" << endl;
    // This function is not ideally named.  It does not necessarily
    // need to say who made the plugin -- getMaker does that -- but it
    // should indicate the terms under which it is distributed.  For
    // example, "Copyright (year). All Rights Reserved", or "GPL"
    return "GPL";
}

NNLSBase::InputDomain
NNLSBase::getInputDomain() const
{
    if (debug_on) cerr << "--> getInputDomain" << endl;
    return FrequencyDomain;
}

size_t
NNLSBase::getPreferredBlockSize() const
{
    if (debug_on) cerr << "--> getPreferredBlockSize" << endl;
    return 16384; // 0 means "I can handle any block size"
}

size_t 
NNLSBase::getPreferredStepSize() const
{
    if (debug_on) cerr << "--> getPreferredStepSize" << endl;
    return 2048; // 0 means "anything sensible"; in practice this
    // means the same as the block size for TimeDomain
    // plugins, or half of it for FrequencyDomain plugins
}

size_t
NNLSBase::getMinChannelCount() const
{
    if (debug_on) cerr << "--> getMinChannelCount" << endl;
    return 1;
}

size_t
NNLSBase::getMaxChannelCount() const
{
    if (debug_on) cerr << "--> getMaxChannelCount" << endl;
    return 1;
}

NNLSBase::ParameterList
NNLSBase::getParameterDescriptors() const
{
    if (debug_on) cerr << "--> getParameterDescriptors" << endl;
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "useNNLS";
    d.name = "use approximate transcription (NNLS)";
    d.description = "Toggles approximate transcription (NNLS).";
    d.unit = "";
    d.minValue = 0.0;
    d.maxValue = 1.0;
    d.defaultValue = 1.0;
    d.isQuantized = true;
	d.quantizeStep = 1.0;
    list.push_back(d);

    ParameterDescriptor d0;
    d0.identifier = "rollon";
    d0.name = "spectral roll-on";
    d0.description = "Consider the cumulative energy spectrum (from low to high frequencies). All bins below the first bin whose cumulative energy exceeds the quantile [spectral roll on] x [total energy] will be set to 0. A value of 0 means that no bins will be changed.";
    d0.unit = "%";
    d0.minValue = 0;
    d0.maxValue = 5;
    d0.defaultValue = 0;
    d0.isQuantized = true;
	d0.quantizeStep = 0.5;
    list.push_back(d0);

    ParameterDescriptor d1;
    d1.identifier = "tuningmode";
    d1.name = "tuning mode";
    d1.description = "Tuning can be performed locally or on the whole extraction segment. Local tuning is only advisable when the tuning is likely to change over the audio, for example in podcasts, or in a cappella singing.";
    d1.unit = "";
    d1.minValue = 0;
    d1.maxValue = 1;
    d1.defaultValue = 0;
    d1.isQuantized = true;
    d1.valueNames.push_back("global tuning");
    d1.valueNames.push_back("local tuning");
    d1.quantizeStep = 1.0;
    list.push_back(d1);

    ParameterDescriptor d2;
    d2.identifier = "whitening";
    d2.name = "spectral whitening";
    d2.description = "Spectral whitening: no whitening - 0; whitening - 1.";
    d2.unit = "";
    d2.isQuantized = true;
    d2.minValue = 0.0;
    d2.maxValue = 1.0;
    d2.defaultValue = 1.0;
    d2.isQuantized = false;
    list.push_back(d2);

    ParameterDescriptor d3;
    d3.identifier = "s";
    d3.name = "spectral shape";
    d3.description = "Determines how individual notes in the note dictionary look: higher values mean more dominant higher harmonics.";
    d3.unit = "";
    d3.minValue = 0.5;
    d3.maxValue = 0.9;
    d3.defaultValue = 0.7;
    d3.isQuantized = false;
    list.push_back(d3);

    ParameterDescriptor d4;
    d4.identifier = "chromanormalize";
    d4.name = "chroma normalization";
    d4.description = "How shall the chroma vector be normalized?";
    d4.unit = "";
    d4.minValue = 0;
    d4.maxValue = 3;
    d4.defaultValue = 0;
    d4.isQuantized = true;
    d4.valueNames.push_back("none");
    d4.valueNames.push_back("maximum norm");
    d4.valueNames.push_back("L1 norm");
    d4.valueNames.push_back("L2 norm");
    d4.quantizeStep = 1.0;
    list.push_back(d4);

    return list;
}

float
NNLSBase::getParameter(string identifier) const
{
    if (debug_on) cerr << "--> getParameter" << endl;
    if (identifier == "useNNLS") {
        return m_useNNLS; 
    }
    
    if (identifier == "whitening") {
        return m_whitening; 
    }

    if (identifier == "s") {
        return m_s; 
    }

    if (identifier == "rollon") {
        return m_rollon; 
    }
    
    if (identifier == "tuningmode") {
        if (m_tuneLocal) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    if (identifier == "preset") {
        return m_preset;
    }
    if (identifier == "chromanormalize") {
        return m_doNormalizeChroma;
    }
    
    if (identifier == "useHMM") {
        return m_useHMM;
    }
    
    return 0;
    
}

void
NNLSBase::setParameter(string identifier, float value) 
{
    if (debug_on) cerr << "--> setParameter" << endl;
    if (identifier == "useNNLS") {
        m_useNNLS = (int) value;
    }
    
    if (identifier == "whitening") {
        m_whitening = value;
    }
    
    if (identifier == "s") {
        m_s = value;
    }

    if (identifier == "useHMM") {
        m_useHMM = value;
    }
    
    if (identifier == "tuningmode") {
        // m_tuneLocal = (value > 0) ? true : false;
        m_tuneLocal = value;
        // cerr << "m_tuneLocal :" << m_tuneLocal << endl;
    }
    // if (identifier == "preset") {
    //     m_preset = value;
    //     if (m_preset == 0.0) {
    //         m_tuneLocal = false;
    //         m_whitening = 1.0;
    //         m_dictID = 0.0;
    //     }
    //     if (m_preset == 1.0) {
    //         m_tuneLocal = false;
    //         m_whitening = 1.0;
    //         m_dictID = 1.0;
    //     }
    //     if (m_preset == 2.0) {
    //         m_tuneLocal = false;
    //         m_whitening = 0.7;
    //         m_dictID = 0.0;
    //     }
    // }
    if (identifier == "chromanormalize") {
        m_doNormalizeChroma = value;
    }
	 
    if (identifier == "rollon") {
        m_rollon = value;
    }
}

NNLSBase::ProgramList
NNLSBase::getPrograms() const
{
    if (debug_on) cerr << "--> getPrograms" << endl;
    ProgramList list;

    // If you have no programs, return an empty list (or simply don't
    // implement this function or getCurrentProgram/selectProgram)

    return list;
}

string
NNLSBase::getCurrentProgram() const
{
    if (debug_on) cerr << "--> getCurrentProgram" << endl;
    return ""; // no programs
}

void
NNLSBase::selectProgram(string name)
{
    if (debug_on) cerr << "--> selectProgram" << endl;
}


bool
NNLSBase::initialise(size_t channels, size_t stepSize, size_t blockSize)
{	
    if (debug_on) {
        cerr << "--> initialise";
    }
	
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;
    m_blockSize = blockSize;
    m_stepSize = stepSize;
    m_frameCount = 0;
    int tempn = 256 * m_blockSize/2;
    // cerr << "length of tempkernel : " <<  tempn << endl;
    float *tempkernel;

    tempkernel = new float[tempn];

    logFreqMatrix(m_inputSampleRate, m_blockSize, tempkernel);
    m_kernelValue.clear();
    m_kernelFftIndex.clear();
    m_kernelNoteIndex.clear();
    int countNonzero = 0;
    for (unsigned iNote = 0; iNote < nNote; ++iNote) { // I don't know if this is wise: manually making a sparse matrix
        for (unsigned iFFT = 0; iFFT < blockSize/2; ++iFFT) {
            if (tempkernel[iFFT + blockSize/2 * iNote] > 0) {
                m_kernelValue.push_back(tempkernel[iFFT + blockSize/2 * iNote]);
                if (tempkernel[iFFT + blockSize/2 * iNote] > 0) {
                    countNonzero++;
                }
                m_kernelFftIndex.push_back(iFFT);
                m_kernelNoteIndex.push_back(iNote);				
            }
        }
    }
    // cerr << "nonzero count : " << countNonzero << endl;
    delete [] tempkernel;
/*
    ofstream myfile;
    myfile.open ("matrix.txt");
    // myfile << "Writing this to a file.\n";	
    for (int i = 0; i < nNote * 84; ++i) {
        myfile << m_dict[i] << endl;		
    }
    myfile.close();
*/
    return true;
}

void
NNLSBase::reset()
{
    if (debug_on) cerr << "--> reset";
	
    // Clear buffers, reset stored values, etc
    m_frameCount = 0;
    // m_dictID = 0;
    m_logSpectrum.clear();
    m_meanTuning0 = 0;
    m_meanTuning1 = 0;
    m_meanTuning2 = 0;
    m_localTuning0 = 0;
    m_localTuning1 = 0;
    m_localTuning2 = 0;
    m_localTuning.clear();
}

void
NNLSBase::baseProcess(const float *const *inputBuffers, Vamp::RealTime timestamp)
{   
    m_frameCount++;   
    float *magnitude = new float[m_blockSize/2];
	
    const float *fbuf = inputBuffers[0];	
    float energysum = 0;
    // make magnitude
    float maxmag = -10000;
    for (size_t iBin = 0; iBin < m_blockSize/2; iBin++) {
        magnitude[iBin] = sqrt(fbuf[2 * iBin] * fbuf[2 * iBin] + 
                               fbuf[2 * iBin + 1] * fbuf[2 * iBin + 1]);		
        if (maxmag < magnitude[iBin]) maxmag = magnitude[iBin];
        if (m_rollon > 0) {
            energysum += pow(magnitude[iBin],2);
        }
    }
	
    float cumenergy = 0;
    if (m_rollon > 0) {
        for (size_t iBin = 2; iBin < m_blockSize/2; iBin++) {
            cumenergy +=  pow(magnitude[iBin],2);
            if (cumenergy < energysum * m_rollon / 100) magnitude[iBin-2] = 0;
            else break;
        }
    }
	
    if (maxmag < 2) {
        // cerr << "timestamp " << timestamp << ": very low magnitude, setting magnitude to all zeros" << endl;
        for (size_t iBin = 0; iBin < m_blockSize/2; iBin++) {
            magnitude[iBin] = 0;
        }
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
        // cerr << m_kernelFftIndex[binCount] << " -- " << magnitude[m_kernelFftIndex[binCount]] << " -- "<< m_kernelValue[binCount] << endl;
        binCount++;	
    }
    // cerr << nm[20];
    // cerr << endl;
	
	
    float one_over_N = 1.0/m_frameCount;
    // update means of complex tuning variables
    m_meanTuning0 *= float(m_frameCount-1)*one_over_N;
    m_meanTuning1 *= float(m_frameCount-1)*one_over_N;
    m_meanTuning2 *= float(m_frameCount-1)*one_over_N;
	
    for (int iTone = 0; iTone < 160; iTone = iTone + 3) {
        m_meanTuning0 += nm[iTone + 0]*one_over_N;
    	m_meanTuning1 += nm[iTone + 1]*one_over_N;
    	m_meanTuning2 += nm[iTone + 2]*one_over_N;
        float ratioOld = 0.997;
        m_localTuning0 *= ratioOld; m_localTuning0 += nm[iTone + 0] * (1 - ratioOld);
        m_localTuning1 *= ratioOld; m_localTuning1 += nm[iTone + 1] * (1 - ratioOld);
        m_localTuning2 *= ratioOld; m_localTuning2 += nm[iTone + 2] * (1 - ratioOld);
    }
	
    // if (m_tuneLocal) {
    // local tuning
    float localTuningImag = sinvalue * m_localTuning1 - sinvalue * m_localTuning2;
    float localTuningReal = m_localTuning0 + cosvalue * m_localTuning1 + cosvalue * m_localTuning2;
    float normalisedtuning = atan2(localTuningImag, localTuningReal)/(2*M_PI);
    m_localTuning.push_back(normalisedtuning);
    
    Feature f1; // logfreqspec
    f1.hasTimestamp = true;
    f1.timestamp = timestamp;
    for (size_t iNote = 0; iNote < nNote; iNote++) {
        f1.values.push_back(nm[iNote]);
    }

    // deletes
    delete[] magnitude;
    delete[] nm;

    m_logSpectrum.push_back(f1); // remember note magnitude
}


#ifdef NOT_DEFINED

NNLSBase::FeatureSet
NNLSBase::getRemainingFeatures()
{
    if (debug_on) cerr << "--> getRemainingFeatures" << endl;
    FeatureSet fsOut;
    if (m_logSpectrum.size() == 0) return fsOut;
    int nChord = m_chordnames.size();
    // 
    /**  Calculate Tuning
         calculate tuning from (using the angle of the complex number defined by the 
         cumulative mean real and imag values)
    **/
    float meanTuningImag = sinvalue * m_meanTuning1 - sinvalue * m_meanTuning2;
    float meanTuningReal = m_meanTuning0 + cosvalue * m_meanTuning1 + cosvalue * m_meanTuning2;
    float cumulativetuning = 440 * pow(2,atan2(meanTuningImag, meanTuningReal)/(24*M_PI));
    float normalisedtuning = atan2(meanTuningImag, meanTuningReal)/(2*M_PI);
    int intShift = floor(normalisedtuning * 3);
    float intFactor = normalisedtuning * 3 - intShift; // intFactor is a really bad name for this
		    
    char buffer0 [50];
		
    sprintf(buffer0, "estimated tuning: %0.1f Hz", cumulativetuning);
		    
    // cerr << "normalisedtuning: " << normalisedtuning << '\n';
		    
    // push tuning to FeatureSet fsOut
    Feature f0; // tuning
    f0.hasTimestamp = true;
    f0.timestamp = Vamp::RealTime::frame2RealTime(0, lrintf(m_inputSampleRate));;
    f0.label = buffer0;
    fsOut[0].push_back(f0);  
		    
    /** Tune Log-Frequency Spectrogram
        calculate a tuned log-frequency spectrogram (f2): use the tuning estimated above (kinda f0) to 
        perform linear interpolation on the existing log-frequency spectrogram (kinda f1).
    **/
    cerr << endl << "[NNLS Chroma Plugin] Tuning Log-Frequency Spectrogram ... ";
					
    float tempValue = 0;
    float dbThreshold = 0; // relative to the background spectrum
    float thresh = pow(10,dbThreshold/20);
    // cerr << "tune local ? " << m_tuneLocal << endl;
    int count = 0;
		
    for (FeatureList::iterator i = m_logSpectrum.begin(); i != m_logSpectrum.end(); ++i) {
        Feature f1 = *i;
        Feature f2; // tuned log-frequency spectrum
        f2.hasTimestamp = true;
        f2.timestamp = f1.timestamp;
        f2.values.push_back(0.0); f2.values.push_back(0.0); // set lower edge to zero
		
        if (m_tuneLocal == 1.0) {
            intShift = floor(m_localTuning[count] * 3);
            intFactor = m_localTuning[count] * 3 - intShift; // intFactor is a really bad name for this
        }
		        
        // cerr << intShift << " " << intFactor << endl;
		        
        for (unsigned k = 2; k < f1.values.size() - 3; ++k) { // interpolate all inner bins
            tempValue = f1.values[k + intShift] * (1-intFactor) + f1.values[k+intShift+1] * intFactor;
            f2.values.push_back(tempValue);
        }
		        
        f2.values.push_back(0.0); f2.values.push_back(0.0); f2.values.push_back(0.0); // upper edge
        vector<float> runningmean = SpecialConvolution(f2.values,hw);
        vector<float> runningstd;
        for (int i = 0; i < 256; i++) { // first step: squared values into vector (variance)
            runningstd.push_back((f2.values[i] - runningmean[i]) * (f2.values[i] - runningmean[i]));
        }
        runningstd = SpecialConvolution(runningstd,hw); // second step convolve
        for (int i = 0; i < 256; i++) { 
            runningstd[i] = sqrt(runningstd[i]); // square root to finally have running std
            if (runningstd[i] > 0) {
                // f2.values[i] = (f2.values[i] / runningmean[i]) > thresh ? 
                // 		                    (f2.values[i] - runningmean[i]) / pow(runningstd[i],m_whitening) : 0;
                f2.values[i] = (f2.values[i] - runningmean[i]) > 0 ?
                    (f2.values[i] - runningmean[i]) / pow(runningstd[i],m_whitening) : 0;
            }
            if (f2.values[i] < 0) {
                cerr << "ERROR: negative value in logfreq spectrum" << endl;
            }
        }
        fsOut[2].push_back(f2);
        count++;
    }
    cerr << "done." << endl;
	    
    /** Semitone spectrum and chromagrams
        Semitone-spaced log-frequency spectrum derived from the tuned log-freq spectrum above. the spectrum
        is inferred using a non-negative least squares algorithm.
        Three different kinds of chromagram are calculated, "treble", "bass", and "both" (which means 
        bass and treble stacked onto each other).
    **/
    if (m_useNNLS == 0) {
        cerr << "[NNLS Chroma Plugin] Mapping to semitone spectrum and chroma ... ";
    } else {
        cerr << "[NNLS Chroma Plugin] Performing NNLS and mapping to chroma ... ";
    }

	    
    vector<vector<float> > chordogram;
    vector<vector<int> > scoreChordogram;
    vector<float> chordchange = vector<float>(fsOut[2].size(),0);
    vector<float> oldchroma = vector<float>(12,0);
    vector<float> oldbasschroma = vector<float>(12,0);
    count = 0;

    for (FeatureList::iterator it = fsOut[2].begin(); it != fsOut[2].end(); ++it) {
        Feature f2 = *it; // logfreq spectrum
        Feature f3; // semitone spectrum
        Feature f4; // treble chromagram
        Feature f5; // bass chromagram
        Feature f6; // treble and bass chromagram
	
        f3.hasTimestamp = true;
        f3.timestamp = f2.timestamp;
	        
        f4.hasTimestamp = true;
        f4.timestamp = f2.timestamp;
	        
        f5.hasTimestamp = true;
        f5.timestamp = f2.timestamp;
	        
        f6.hasTimestamp = true;
        f6.timestamp = f2.timestamp;
	        
			float b[256];
	
        bool some_b_greater_zero = false;
        float sumb = 0;
        for (int i = 0; i < 256; i++) {
            // b[i] = m_dict[(256 * count + i) % (256 * 84)];
            b[i] = f2.values[i];
            sumb += b[i];
            if (b[i] > 0) {
                some_b_greater_zero = true;
            }            
        }
	    
        // here's where the non-negative least squares algorithm calculates the note activation x
	
        vector<float> chroma = vector<float>(12, 0);
        vector<float> basschroma = vector<float>(12, 0);
        float currval;
        unsigned iSemitone = 0;
			
        if (some_b_greater_zero) {
            if (m_useNNLS == 0) {
                for (unsigned iNote = 2; iNote < nNote - 2; iNote += 3) {
                    currval = 0;
                    currval += b[iNote + 1 + -1] * 0.5;						
                    currval += b[iNote + 1 +  0] * 1.0;						
                    currval += b[iNote + 1 +  1] * 0.5;						
                    f3.values.push_back(currval);
                    chroma[iSemitone % 12] += currval * treblewindow[iSemitone];
                    basschroma[iSemitone % 12] += currval * basswindow[iSemitone];
                    iSemitone++;
                }
		        
            } else {
					float x[84+1000];
                for (int i = 1; i < 1084; ++i) x[i] = 1.0;
                vector<int> signifIndex;
                int index=0;
                sumb /= 84.0;
                for (unsigned iNote = 2; iNote < nNote - 2; iNote += 3) {
                    float currval = 0;
                    currval += b[iNote + 1 + -1];						
                    currval += b[iNote + 1 +  0];						
                    currval += b[iNote + 1 +  1];
                    if (currval > 0) signifIndex.push_back(index);
                    f3.values.push_back(0); // fill the values, change later
                    index++;
                }
				    float rnorm;
				    float w[84+1000];
				    float zz[84+1000];
                int indx[84+1000];
                int mode;
                int dictsize = 256*signifIndex.size();
                // cerr << "dictsize is " << dictsize << "and values size" << f3.values.size()<< endl;
					float *curr_dict = new float[dictsize];
                for (unsigned iNote = 0; iNote < signifIndex.size(); ++iNote) {
                    for (unsigned iBin = 0; iBin < 256; iBin++) {
                        curr_dict[iNote * 256 + iBin] = 1.0 * m_dict[signifIndex[iNote] * 256 + iBin];
                    }
                }
					nnls(curr_dict, nNote, nNote, signifIndex.size(), b, x, &rnorm, w, zz, indx, &mode);
                delete [] curr_dict;
                for (unsigned iNote = 0; iNote < signifIndex.size(); ++iNote) {
                    f3.values[signifIndex[iNote]] = x[iNote];
                    // cerr << mode << endl;
                    chroma[signifIndex[iNote] % 12] += x[iNote] * treblewindow[signifIndex[iNote]];
                    basschroma[signifIndex[iNote] % 12] += x[iNote] * basswindow[signifIndex[iNote]];
                }
            }	
        }
			
            
	        
			
        f4.values = chroma; 
        f5.values = basschroma;
        chroma.insert(chroma.begin(), basschroma.begin(), basschroma.end()); // just stack the both chromas 
        f6.values = chroma; 
	        
        if (m_doNormalizeChroma > 0) {
            vector<float> chromanorm = vector<float>(3,0);			
            switch (int(m_doNormalizeChroma)) {
            case 0: // should never end up here
                break;
            case 1:
                chromanorm[0] = *max_element(f4.values.begin(), f4.values.end());
                chromanorm[1] = *max_element(f5.values.begin(), f5.values.end());
                chromanorm[2] = max(chromanorm[0], chromanorm[1]);
                break;
            case 2:
                for (vector<float>::iterator it = f4.values.begin(); it != f4.values.end(); ++it) {
                    chromanorm[0] += *it; 						
                }
                for (vector<float>::iterator it = f5.values.begin(); it != f5.values.end(); ++it) {
                    chromanorm[1] += *it; 						
                }
                for (vector<float>::iterator it = f6.values.begin(); it != f6.values.end(); ++it) {
                    chromanorm[2] += *it; 						
                }
                break;
            case 3:
                for (vector<float>::iterator it = f4.values.begin(); it != f4.values.end(); ++it) {
                    chromanorm[0] += pow(*it,2); 						
                }
                chromanorm[0] = sqrt(chromanorm[0]);
                for (vector<float>::iterator it = f5.values.begin(); it != f5.values.end(); ++it) {
                    chromanorm[1] += pow(*it,2); 						
                }
                chromanorm[1] = sqrt(chromanorm[1]);
                for (vector<float>::iterator it = f6.values.begin(); it != f6.values.end(); ++it) {
                    chromanorm[2] += pow(*it,2); 						
                }
                chromanorm[2] = sqrt(chromanorm[2]);
                break;
            }
            if (chromanorm[0] > 0) {
                for (int i = 0; i < f4.values.size(); i++) {
                    f4.values[i] /= chromanorm[0];
                }
            }
            if (chromanorm[1] > 0) {
                for (int i = 0; i < f5.values.size(); i++) {
                    f5.values[i] /= chromanorm[1];
                }
            }
            if (chromanorm[2] > 0) {
                for (int i = 0; i < f6.values.size(); i++) {
                    f6.values[i] /= chromanorm[2];
                }
            }
				
        }
	
        // local chord estimation
        vector<float> currentChordSalience;
        float tempchordvalue = 0;
        float sumchordvalue = 0;
	        
        for (int iChord = 0; iChord < nChord; iChord++) {
            tempchordvalue = 0;
            for (int iBin = 0; iBin < 12; iBin++) {
                tempchordvalue += m_chorddict[24 * iChord + iBin] * chroma[iBin];
            }
            for (int iBin = 12; iBin < 24; iBin++) {
                tempchordvalue += m_chorddict[24 * iChord + iBin] * chroma[iBin];
            }
            sumchordvalue+=tempchordvalue;
            currentChordSalience.push_back(tempchordvalue);
        }
        if (sumchordvalue > 0) {
            for (int iChord = 0; iChord < nChord; iChord++) {
                currentChordSalience[iChord] /= sumchordvalue;
            }
        } else {
            currentChordSalience[nChord-1] = 1.0;
        }
        chordogram.push_back(currentChordSalience);
	        
        fsOut[3].push_back(f3);
        fsOut[4].push_back(f4);
        fsOut[5].push_back(f5);
        fsOut[6].push_back(f6);
        count++;
    }
    cerr << "done." << endl;
		

    /* Simple chord estimation
       I just take the local chord estimates ("currentChordSalience") and average them over time, then
       take the maximum. Very simple, don't do this at home...
    */
    cerr << "[NNLS Chroma Plugin] Chord Estimation ... ";
    count = 0; 
    int halfwindowlength = m_inputSampleRate / m_stepSize;
    vector<int> chordSequence;
    for (FeatureList::iterator it = fsOut[6].begin(); it != fsOut[6].end(); ++it) { // initialise the score chordogram
        vector<int> temp = vector<int>(nChord,0);
        scoreChordogram.push_back(temp);
    }
    for (FeatureList::iterator it = fsOut[6].begin(); it < fsOut[6].end()-2*halfwindowlength-1; ++it) {		
        int startIndex = count + 1;
        int endIndex = count + 2 * halfwindowlength;
			
        float chordThreshold = 2.5/nChord;//*(2*halfwindowlength+1);
            
        vector<int> chordCandidates;
        for (unsigned iChord = 0; iChord < nChord-1; iChord++) {
            // float currsum = 0;
            // for (unsigned iFrame = startIndex; iFrame < endIndex; ++iFrame) {
            //  currsum += chordogram[iFrame][iChord];
            // }
            //                 if (currsum > chordThreshold) chordCandidates.push_back(iChord);
            for (unsigned iFrame = startIndex; iFrame < endIndex; ++iFrame) {
                if (chordogram[iFrame][iChord] > chordThreshold) {
                    chordCandidates.push_back(iChord);
                    break;
                }                    
            }
        }
        chordCandidates.push_back(nChord-1);
        // cerr << chordCandidates.size() << endl;          
	        
        float maxval = 0; // will be the value of the most salient *chord change* in this frame
        float maxindex = 0; //... and the index thereof
        unsigned bestchordL = nChord-1; // index of the best "left" chord
        unsigned bestchordR = nChord-1; // index of the best "right" chord
 	 		
        for (int iWF = 1; iWF < 2*halfwindowlength; ++iWF) {
            // now find the max values on both sides of iWF
            // left side:
            float maxL = 0;
            unsigned maxindL = nChord-1;
            for (unsigned kChord = 0; kChord < chordCandidates.size(); kChord++) {
                unsigned iChord = chordCandidates[kChord];
                float currsum = 0;
                for (unsigned iFrame = 0; iFrame < iWF-1; ++iFrame) {
                    currsum += chordogram[count+iFrame][iChord];
                }
                if (iChord == nChord-1) currsum *= 0.8;
                if (currsum > maxL) {
                    maxL = currsum;
                    maxindL = iChord;
                }
            }				
            // right side:
            float maxR = 0;
            unsigned maxindR = nChord-1;
            for (unsigned kChord = 0; kChord < chordCandidates.size(); kChord++) {
                unsigned iChord = chordCandidates[kChord];
                float currsum = 0;
                for (unsigned iFrame = iWF-1; iFrame < 2*halfwindowlength; ++iFrame) {
                    currsum += chordogram[count+iFrame][iChord];
                }
                if (iChord == nChord-1) currsum *= 0.8;
                if (currsum > maxR) {
                    maxR = currsum;
                    maxindR = iChord;
                }
            }
            if (maxL+maxR > maxval) {					
                maxval = maxL+maxR;
                maxindex = iWF;
                bestchordL = maxindL;
                bestchordR = maxindR;
            }
				
        }
        // cerr << "maxindex: " << maxindex << ", bestchordR is " << bestchordR << ", of frame " << count << endl;
        // add a score to every chord-frame-point that was part of a maximum 
        for (unsigned iFrame = 0; iFrame < maxindex-1; ++iFrame) {
            scoreChordogram[iFrame+count][bestchordL]++;
        }
        for (unsigned iFrame = maxindex-1; iFrame < 2*halfwindowlength; ++iFrame) {
            scoreChordogram[iFrame+count][bestchordR]++;
        }
        if (bestchordL != bestchordR) chordchange[maxindex+count] += (halfwindowlength - abs(maxindex-halfwindowlength)) * 2.0 / halfwindowlength;
        count++;	
    }
    // cerr << "*******  agent finished   *******" << endl;
    count = 0;
    for (FeatureList::iterator it = fsOut[6].begin(); it != fsOut[6].end(); ++it) { 
        float maxval = 0; // will be the value of the most salient chord in this frame
        float maxindex = 0; //... and the index thereof
        for (unsigned iChord = 0; iChord < nChord; iChord++) {
            if (scoreChordogram[count][iChord] > maxval) {
                maxval = scoreChordogram[count][iChord];
                maxindex = iChord;
                // cerr << iChord << endl;
            }
        }
        chordSequence.push_back(maxindex);
        // cerr << "before modefilter, maxindex: " << maxindex << endl;
        count++;
    }
    // cerr << "*******  mode filter done *******" << endl;

	
    // mode filter on chordSequence
    count = 0;
    string oldChord = "";
    for (FeatureList::iterator it = fsOut[6].begin(); it != fsOut[6].end(); ++it) {
        Feature f6 = *it;
        Feature f7; // chord estimate
        f7.hasTimestamp = true;
        f7.timestamp = f6.timestamp;
        Feature f8; // chord estimate
        f8.hasTimestamp = true;
        f8.timestamp = f6.timestamp;
			
        vector<int> chordCount = vector<int>(nChord,0);
        int maxChordCount = 0;
        int maxChordIndex = nChord-1;
        string maxChord;
        int startIndex = max(count - halfwindowlength/2,0);
        int endIndex = min(int(chordogram.size()), count + halfwindowlength/2);
        for (int i = startIndex; i < endIndex; i++) {				
            chordCount[chordSequence[i]]++;
            if (chordCount[chordSequence[i]] > maxChordCount) {
                // cerr << "start index " << startIndex << endl;
                maxChordCount++;
                maxChordIndex = chordSequence[i];
                maxChord = m_chordnames[maxChordIndex];
            }
        }
        // chordSequence[count] = maxChordIndex;
        // cerr << maxChordIndex << endl;
        f8.values.push_back(chordchange[count]/(halfwindowlength*2));
        // cerr << chordchange[count] << endl;
        fsOut[9].push_back(f8);
        if (oldChord != maxChord) {
            oldChord = maxChord;
	
            // char buffer1 [50];
            // if (maxChordIndex < nChord - 1) {
            //     sprintf(buffer1, "%s%s", notenames[maxChordIndex % 12 + 12], chordtypes[maxChordIndex]);
            // } else {
            //     sprintf(buffer1, "N");
            // }
            // f7.label = buffer1;
            f7.label = m_chordnames[maxChordIndex];
            fsOut[7].push_back(f7);
        }
        count++;
    }
    Feature f7; // last chord estimate
    f7.hasTimestamp = true;
    f7.timestamp = fsOut[6][fsOut[6].size()-1].timestamp;
    f7.label = "N";
    fsOut[7].push_back(f7);
    cerr << "done." << endl;
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

#endif
