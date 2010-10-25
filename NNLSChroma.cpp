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

#include "NNLSChroma.h"

#include "chromamethods.h"

#include <cstdlib>
#include <fstream>
#include <cmath>

#include <algorithm>

const bool debug_on = false;

const vector<float> hw(hammingwind, hammingwind+19);

NNLSChroma::NNLSChroma(float inputSampleRate) :
    NNLSBase(inputSampleRate)
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
    return "nnls-chroma"; 
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
    if (debug_on) cerr << "--> getDescription" << endl;
    return "This plugin provides a number of features derived from a DFT-based log-frequency amplitude spectrum: some variants of the log-frequency spectrum, including a semitone spectrum derived from approximate transcription using the NNLS algorithm; and based on this semitone spectrum, different chroma features.";
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
            chromanames.push_back(notenames[iNote+12]);
        }
    }
    
    int index = 0;

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
    m_outputLogSpec = index++;

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
    m_outputTunedSpec = index++;
    
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
    m_outputSemiSpec = index++;
    
    OutputDescriptor d4;
    d4.identifier = "chroma";
    d4.name = "Chromagram";
    d4.description = "Tuning-adjusted chromagram from NNLS approximate transcription, with an emphasis on the medium note range.";
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
    m_outputChroma = index++;
    
    OutputDescriptor d5;
    d5.identifier = "basschroma";
    d5.name = "Bass Chromagram";
    d5.description = "Tuning-adjusted bass chromagram from NNLS approximate transcription, with an emphasis on the bass note range.";
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
    m_outputBassChroma = index++;
    
    OutputDescriptor d6;
    d6.identifier = "bothchroma";
    d6.name = "Chromagram and Bass Chromagram";
    d6.description = "Tuning-adjusted chromagram and bass chromagram (stacked on top of each other) from NNLS approximate transcription.";
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
    m_outputBothChroma = index++;
  
    return list;
}


bool
NNLSChroma::initialise(size_t channels, size_t stepSize, size_t blockSize)
{	
    if (debug_on) {
        cerr << "--> initialise";
    }

    if (!NNLSBase::initialise(channels, stepSize, blockSize)) {
        return false;
    }

    return true;
}

void
NNLSChroma::reset()
{
    if (debug_on) cerr << "--> reset";
    NNLSBase::reset();
}

NNLSChroma::FeatureSet
NNLSChroma::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{   
    if (debug_on) cerr << "--> process" << endl;

    NNLSBase::baseProcess(inputBuffers, timestamp);
	
    FeatureSet fs;
    fs[m_outputLogSpec].push_back(m_logSpectrum[m_logSpectrum.size()-1]);
    return fs;	
}

NNLSChroma::FeatureSet
NNLSChroma::getRemainingFeatures()
{
    if (debug_on) cerr << "--> getRemainingFeatures" << endl;
    FeatureSet fsOut;
    if (m_logSpectrum.size() == 0) return fsOut;
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
		
        if (m_tuneLocal) {
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
        fsOut[m_outputTunedSpec].push_back(f2);
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

	    
    vector<float> oldchroma = vector<float>(12,0);
    vector<float> oldbasschroma = vector<float>(12,0);
    count = 0;

    for (FeatureList::iterator it = fsOut[m_outputTunedSpec].begin(); it != fsOut[m_outputTunedSpec].end(); ++it) {
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
	
        fsOut[m_outputSemiSpec].push_back(f3);
        fsOut[m_outputChroma].push_back(f4);
        fsOut[m_outputBassChroma].push_back(f5);
        fsOut[m_outputBothChroma].push_back(f6);
        count++;
    }
    cerr << "done." << endl;

    return fsOut;     

}

