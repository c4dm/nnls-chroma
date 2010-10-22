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

#include "Chordino.h"

#include "chromamethods.h"

#include <cstdlib>
#include <fstream>
#include <cmath>

#include <algorithm>

const bool debug_on = false;

const vector<float> hw(hammingwind, hammingwind+19);

Chordino::Chordino(float inputSampleRate) :
    NNLSBase(inputSampleRate)
{
    if (debug_on) cerr << "--> Chordino" << endl;
}

Chordino::~Chordino()
{
    if (debug_on) cerr << "--> ~Chordino" << endl;
}

string
Chordino::getIdentifier() const
{
    if (debug_on) cerr << "--> getIdentifier" << endl;
    return "chordino";
}

string
Chordino::getName() const
{
    if (debug_on) cerr << "--> getName" << endl;
    return "Chordino";
}

string
Chordino::getDescription() const
{
    if (debug_on) cerr << "--> getDescription" << endl;
    return "This plugin provides a number of features derived from a log-frequency amplitude spectrum of the DFT: some variants of the log-frequency spectrum, including a semitone spectrum derived from approximate transcription using the NNLS algorithm; based on this semitone spectrum, chroma features and a simple chord estimate.";
}

Chordino::OutputList
Chordino::getOutputDescriptors() const
{
    if (debug_on) cerr << "--> getOutputDescriptors" << endl;
    OutputList list;
    
    int index = 0;

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
    m_outputChords = index++;
    
    OutputDescriptor d8;
    d8.identifier = "harmonicchange";
    d8.name = "Harmonic change value";
    d8.description = "Harmonic change.";
    d8.unit = "";
    d8.hasFixedBinCount = true;
    d8.binCount = 1;
    d8.hasKnownExtents = true;
    d8.minValue = 0.0;
    d8.maxValue = 0.999;
    d8.isQuantized = false;
    d8.sampleType = OutputDescriptor::FixedSampleRate;
    d8.hasDuration = false;
    // d8.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d8);
    m_outputHarmonicChange = index++;
  
    return list;
}

bool
Chordino::initialise(size_t channels, size_t stepSize, size_t blockSize)
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
Chordino::reset()
{
    if (debug_on) cerr << "--> reset";
    NNLSBase::reset();
}

Chordino::FeatureSet
Chordino::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{   
    if (debug_on) cerr << "--> process" << endl;

    NNLSBase::baseProcess(inputBuffers, timestamp);

    return FeatureSet();
}

Chordino::FeatureSet
Chordino::getRemainingFeatures()
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
		    
		    
    /** Tune Log-Frequency Spectrogram
        calculate a tuned log-frequency spectrogram (f2): use the tuning estimated above (kinda f0) to 
        perform linear interpolation on the existing log-frequency spectrogram (kinda f1).
    **/
    cerr << endl << "[Chordino Plugin] Tuning Log-Frequency Spectrogram ... ";
					
    float tempValue = 0;
    float dbThreshold = 0; // relative to the background spectrum
    float thresh = pow(10,dbThreshold/20);
    // cerr << "tune local ? " << m_tuneLocal << endl;
    int count = 0;
		
    FeatureList tunedSpec;

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
                // 		                    (f2.values[i] - runningmean[i]) / pow(runningstd[i],m_paling) : 0;
                f2.values[i] = (f2.values[i] - runningmean[i]) > 0 ?
                    (f2.values[i] - runningmean[i]) / pow(runningstd[i],m_paling) : 0;
            }
            if (f2.values[i] < 0) {
                cerr << "ERROR: negative value in logfreq spectrum" << endl;
            }
        }
        tunedSpec.push_back(f2);
        count++;
    }
    cerr << "done." << endl;
	    
    /** Semitone spectrum and chromagrams
        Semitone-spaced log-frequency spectrum derived from the tuned log-freq spectrum above. the spectrum
        is inferred using a non-negative least squares algorithm.
        Three different kinds of chromagram are calculated, "treble", "bass", and "both" (which means 
        bass and treble stacked onto each other).
    **/
    if (m_dictID == 1) {
        cerr << "[Chordino Plugin] Mapping to semitone spectrum and chroma ... ";
    } else {
        cerr << "[Chordino Plugin] Performing NNLS and mapping to chroma ... ";
    }

	    
    vector<vector<float> > chordogram;
    vector<vector<int> > scoreChordogram;
    vector<float> chordchange = vector<float>(tunedSpec.size(),0);
    count = 0;

    FeatureList chromaList;

    for (FeatureList::iterator it = tunedSpec.begin(); it != tunedSpec.end(); ++it) {
        Feature f2 = *it; // logfreq spectrum
        Feature f6; // treble and bass chromagram

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
            if (m_dictID == 1) {
                for (unsigned iNote = 2; iNote < nNote - 2; iNote += 3) {
                    currval = 0;
                    currval += b[iNote + 1 + -1] * 0.5;
                    currval += b[iNote + 1 +  0] * 1.0; 
                    currval += b[iNote + 1 +  1] * 0.5; 
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
                    index++;
                }
                float rnorm;
                float w[84+1000];
                float zz[84+1000];
                int indx[84+1000];
                int mode;
                int dictsize = 256*signifIndex.size();
                float *curr_dict = new float[dictsize];
                for (unsigned iNote = 0; iNote < signifIndex.size(); ++iNote) {
                    for (unsigned iBin = 0; iBin < 256; iBin++) {
                        curr_dict[iNote * 256 + iBin] = 1.0 * m_dict[signifIndex[iNote] * 256 + iBin];
                    }
                }
                nnls(curr_dict, nNote, nNote, signifIndex.size(), b, x, &rnorm, w, zz, indx, &mode);
                delete [] curr_dict;
                for (unsigned iNote = 0; iNote < signifIndex.size(); ++iNote) {
                    // cerr << mode << endl;
                    chroma[signifIndex[iNote] % 12] += x[iNote] * treblewindow[signifIndex[iNote]];
                    basschroma[signifIndex[iNote] % 12] += x[iNote] * basswindow[signifIndex[iNote]];
                }
            }	
        }

        vector<float> origchroma = chroma;
        chroma.insert(chroma.begin(), basschroma.begin(), basschroma.end()); // just stack the both chromas 
        f6.values = chroma;
 
        if (m_doNormalizeChroma > 0) {
            vector<float> chromanorm = vector<float>(3,0);			
            switch (int(m_doNormalizeChroma)) {
            case 0: // should never end up here
                break;
            case 1:
                chromanorm[0] = *max_element(origchroma.begin(), origchroma.end());
                chromanorm[1] = *max_element(basschroma.begin(), basschroma.end());
                chromanorm[2] = max(chromanorm[0], chromanorm[1]);
                break;
            case 2:
                for (vector<float>::iterator it = chroma.begin(); it != chroma.end(); ++it) {
                    chromanorm[2] += *it; 						
                }
                break;
            case 3:
                for (vector<float>::iterator it = chroma.begin(); it != chroma.end(); ++it) {
                    chromanorm[2] += pow(*it,2); 						
                }
                chromanorm[2] = sqrt(chromanorm[2]);
                break;
            }
            if (chromanorm[2] > 0) {
                for (int i = 0; i < chroma.size(); i++) {
                    f6.values[i] /= chromanorm[2];
                }
            }
        }

        chromaList.push_back(f6);

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
	        
        count++;
    }
    cerr << "done." << endl;
		

    /* Simple chord estimation
       I just take the local chord estimates ("currentChordSalience") and average them over time, then
       take the maximum. Very simple, don't do this at home...
    */
    cerr << "[Chordino Plugin] Chord Estimation ... ";
    count = 0; 
    int halfwindowlength = m_inputSampleRate / m_stepSize;
    vector<int> chordSequence;

    for (FeatureList::iterator it = chromaList.begin(); it != chromaList.end(); ++it) { // initialise the score chordogram
        vector<int> temp = vector<int>(nChord,0);
        scoreChordogram.push_back(temp);
    }

    for (FeatureList::iterator it = chromaList.begin(); it < chromaList.end()-2*halfwindowlength-1; ++it) {		
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
//        cerr << chordCandidates.size() << endl;          
	        
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
//        cerr << "maxindex: " << maxindex << ", bestchordR is " << bestchordR << ", of frame " << count << endl;
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
//    cerr << "*******  agent finished   *******" << endl;
    count = 0;
    for (FeatureList::iterator it = chromaList.begin(); it != chromaList.end(); ++it) { 
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
//    cerr << "*******  mode filter done *******" << endl;

	
    // mode filter on chordSequence
    count = 0;
    string oldChord = "";
    for (FeatureList::iterator it = chromaList.begin(); it != chromaList.end(); ++it) {
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
        fsOut[m_outputHarmonicChange].push_back(f8);
        if (oldChord != maxChord) {
            oldChord = maxChord;
            f7.label = m_chordnames[maxChordIndex];
            fsOut[m_outputChords].push_back(f7);
        }
        count++;
    }
    Feature f7; // last chord estimate
    f7.hasTimestamp = true;
    f7.timestamp = chromaList[chromaList.size()-1].timestamp;
    f7.label = "N";
    fsOut[m_outputChords].push_back(f7);
    cerr << "done." << endl;

    return fsOut;     

}

