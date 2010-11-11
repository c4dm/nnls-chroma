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
#include "viterbi.h"

#include <cstdlib>
#include <fstream>
#include <cmath>

#include <algorithm>

const bool debug_on = false;

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
    return "Chordino provides a simple chord transcription based on NNLS Chroma (as in the NNLS Chroma plugin). Chord profiles given by the user in the file chord.dict are used to calculate frame-wise chord similarities. Two simple (non-state-of-the-art!) algorithms are available that smooth these to provide a chord transcription: a simple chord change method, and a standard HMM/Viterbi approach.";
}

Chordino::ParameterList
Chordino::getParameterDescriptors() const
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

    ParameterDescriptor d4;
    d4.identifier = "useHMM";
    d4.name = "HMM (Viterbi decoding)";
    d4.description = "Turns on Viterbi decoding (when off, the simple chord estimator is used).";
    d4.unit = "";
    d4.minValue = 0.0;
    d4.maxValue = 1.0;
    d4.defaultValue = 1.0;
    d4.isQuantized = true;
	d4.quantizeStep = 1.0;
    list.push_back(d4);

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

    // ParameterDescriptor d4;
    // d4.identifier = "chromanormalize";
    // d4.name = "chroma normalization";
    // d4.description = "How shall the chroma vector be normalized?";
    // d4.unit = "";
    // d4.minValue = 0;
    // d4.maxValue = 3;
    // d4.defaultValue = 0;
    // d4.isQuantized = true;
    // d4.valueNames.push_back("none");
    // d4.valueNames.push_back("maximum norm");
    // d4.valueNames.push_back("L1 norm");
    // d4.valueNames.push_back("L2 norm");
    // d4.quantizeStep = 1.0;
    // list.push_back(d4);

    return list;
}

Chordino::OutputList
Chordino::getOutputDescriptors() const
{
    if (debug_on) cerr << "--> getOutputDescriptors" << endl;
    OutputList list;
    
    int index = 0;

    OutputDescriptor d7;
    d7.identifier = "simplechord";
    d7.name = "Chord Estimate";
    d7.description = "Estimated chord times and labels. Two simple (non-state-of-the-art!) algorithms are available that smooth these to provide a chord transcription: a simple chord change method, and a standard HMM/Viterbi approach.";
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
    d8.name = "Harmonic Change Value";
    d8.description = "An indication of the likelihood of harmonic change. Depends on the chord dictionary. Calculation is different depending on whether the Viterbi algorithm is used for chord estimation, or the simple chord estimate.";
    d8.unit = "";
    d8.hasFixedBinCount = true;
    d8.binCount = 1;
    d8.hasKnownExtents = false;
    // d8.minValue = 0.0;
    // d8.maxValue = 0.999;
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
    cerr << hw[0] << hw[1] << endl;
    if (debug_on) cerr << "--> getRemainingFeatures" << endl;
    FeatureSet fsOut;
    if (m_logSpectrum.size() == 0) return fsOut;
    int nChord = m_chordnames.size();
    // 
    /**  Calculate Tuning
         calculate tuning from (using the angle of the complex number defined by the 
         cumulative mean real and imag values)
    **/
    float meanTuningImag = 0;
    float meanTuningReal = 0;
    for (int iBPS = 0; iBPS < nBPS; ++iBPS) {
        meanTuningReal += m_meanTunings[iBPS] * cosvalues[iBPS];
        meanTuningImag += m_meanTunings[iBPS] * sinvalues[iBPS];
    }
    float cumulativetuning = 440 * pow(2,atan2(meanTuningImag, meanTuningReal)/(24*M_PI));
    float normalisedtuning = atan2(meanTuningImag, meanTuningReal)/(2*M_PI);
    int intShift = floor(normalisedtuning * 3);
    float floatShift = normalisedtuning * 3 - intShift; // floatShift is a really bad name for this
		    
    char buffer0 [50];
		
    sprintf(buffer0, "estimated tuning: %0.1f Hz", cumulativetuning);
		    
		    
    /** Tune Log-Frequency Spectrogram
        calculate a tuned log-frequency spectrogram (currentTunedSpec): use the tuning estimated above (kinda f0) to 
        perform linear interpolation on the existing log-frequency spectrogram (kinda currentLogSpectum).
    **/
    cerr << endl << "[Chordino Plugin] Tuning Log-Frequency Spectrogram ... ";
					
    float tempValue = 0;
    float dbThreshold = 0; // relative to the background spectrum
    float thresh = pow(10,dbThreshold/20);
    // cerr << "tune local ? " << m_tuneLocal << endl;
    int count = 0;
		
    FeatureList tunedSpec;
    int nFrame = m_logSpectrum.size();
    
    vector<Vamp::RealTime> timestamps;

    for (FeatureList::iterator i = m_logSpectrum.begin(); i != m_logSpectrum.end(); ++i) {
        Feature currentLogSpectum = *i;
        Feature currentTunedSpec; // tuned log-frequency spectrum
        currentTunedSpec.hasTimestamp = true;
        currentTunedSpec.timestamp = currentLogSpectum.timestamp;
        timestamps.push_back(currentLogSpectum.timestamp);
        currentTunedSpec.values.push_back(0.0); currentTunedSpec.values.push_back(0.0); // set lower edge to zero
		
        if (m_tuneLocal) {
            intShift = floor(m_localTuning[count] * 3);
            floatShift = m_localTuning[count] * 3 - intShift; // floatShift is a really bad name for this
        }
		        
        // cerr << intShift << " " << floatShift << endl;
		        
        for (unsigned k = 2; k < currentLogSpectum.values.size() - 3; ++k) { // interpolate all inner bins
            tempValue = currentLogSpectum.values[k + intShift] * (1-floatShift) + currentLogSpectum.values[k+intShift+1] * floatShift;
            currentTunedSpec.values.push_back(tempValue);
        }
		        
        currentTunedSpec.values.push_back(0.0); currentTunedSpec.values.push_back(0.0); currentTunedSpec.values.push_back(0.0); // upper edge
        vector<float> runningmean = SpecialConvolution(currentTunedSpec.values,hw);
        vector<float> runningstd;
        for (int i = 0; i < nNote; i++) { // first step: squared values into vector (variance)
            runningstd.push_back((currentTunedSpec.values[i] - runningmean[i]) * (currentTunedSpec.values[i] - runningmean[i]));
        }
        runningstd = SpecialConvolution(runningstd,hw); // second step convolve
        for (int i = 0; i < nNote; i++) { 
            runningstd[i] = sqrt(runningstd[i]); // square root to finally have running std
            if (runningstd[i] > 0) {
                // currentTunedSpec.values[i] = (currentTunedSpec.values[i] / runningmean[i]) > thresh ? 
                // 		                    (currentTunedSpec.values[i] - runningmean[i]) / pow(runningstd[i],m_whitening) : 0;
                currentTunedSpec.values[i] = (currentTunedSpec.values[i] - runningmean[i]) > 0 ?
                    (currentTunedSpec.values[i] - runningmean[i]) / pow(runningstd[i],m_whitening) : 0;
            }
            if (currentTunedSpec.values[i] < 0) {
                cerr << "ERROR: negative value in logfreq spectrum" << endl;
            }
        }
        tunedSpec.push_back(currentTunedSpec);
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
        cerr << "[Chordino Plugin] Mapping to semitone spectrum and chroma ... ";
    } else {
        cerr << "[Chordino Plugin] Performing NNLS and mapping to chroma ... ";
    }

	    
    vector<vector<double> > chordogram;
    vector<vector<int> > scoreChordogram;
    vector<float> chordchange = vector<float>(tunedSpec.size(),0);
    count = 0;

    FeatureList chromaList;
    
    

    for (FeatureList::iterator it = tunedSpec.begin(); it != tunedSpec.end(); ++it) {
        Feature currentTunedSpec = *it; // logfreq spectrum
        Feature currentChromas; // treble and bass chromagram

        currentChromas.hasTimestamp = true;
        currentChromas.timestamp = currentTunedSpec.timestamp;    

        float b[nNote];
	
        bool some_b_greater_zero = false;
        float sumb = 0;
        for (int i = 0; i < nNote; i++) {
            // b[i] = m_dict[(nNote * count + i) % (nNote * 84)];
            b[i] = currentTunedSpec.values[i];
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
                int dictsize = nNote*signifIndex.size();
                float *curr_dict = new float[dictsize];
                for (unsigned iNote = 0; iNote < signifIndex.size(); ++iNote) {
                    for (unsigned iBin = 0; iBin < nNote; iBin++) {
                        curr_dict[iNote * nNote + iBin] = 1.0 * m_dict[signifIndex[iNote] * nNote + iBin];
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
        currentChromas.values = chroma;
 
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
                    currentChromas.values[i] /= chromanorm[2];
                }
            }
        }

        chromaList.push_back(currentChromas);

        // local chord estimation
        vector<double> currentChordSalience;
        double tempchordvalue = 0;
        double sumchordvalue = 0;
	        
        for (int iChord = 0; iChord < nChord; iChord++) {
            tempchordvalue = 0;
            for (int iBin = 0; iBin < 12; iBin++) {
                tempchordvalue += m_chorddict[24 * iChord + iBin] * chroma[iBin];                
            }
            for (int iBin = 12; iBin < 24; iBin++) {
                tempchordvalue += m_chorddict[24 * iChord + iBin] * chroma[iBin];
            }
            if (iChord == nChord-1) tempchordvalue *= .7;
            if (tempchordvalue < 0) tempchordvalue = 0.0;
            tempchordvalue = pow(1.3,tempchordvalue);
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
		

    // bool m_useHMM = true; // this will go into the chordino header file.
	if (m_useHMM == 1.0) {
        cerr << "[Chordino Plugin] HMM Chord Estimation ... ";
        int oldchord = nChord-1;
        double selftransprob = 0.99;
	    
        // vector<double> init = vector<double>(nChord,1.0/nChord);
        vector<double> init = vector<double>(nChord,0); init[nChord-1] = 1;
        
        double *delta;
        delta = (double *)malloc(sizeof(double)*nFrame*nChord);                
        
        vector<vector<double> > trans;
        for (int iChord = 0; iChord < nChord; iChord++) {
            vector<double> temp = vector<double>(nChord,(1-selftransprob)/(nChord-1));            
            temp[iChord] = selftransprob;
            trans.push_back(temp);
        }
        vector<int> chordpath = ViterbiPath(init, trans, chordogram, delta);


        Feature chord_feature; // chord estimate
        chord_feature.hasTimestamp = true;
        chord_feature.timestamp = timestamps[0];
        chord_feature.label = m_chordnames[chordpath[0]];
        fsOut[m_outputChords].push_back(chord_feature);
        
        chordchange[0] = 0;
        for (int iFrame = 1; iFrame < chordpath.size(); ++iFrame) {
            // cerr << chordpath[iFrame] << endl;
            if (chordpath[iFrame] != oldchord ) {
                Feature chord_feature; // chord estimate
                chord_feature.hasTimestamp = true;
                chord_feature.timestamp = timestamps[iFrame];
                chord_feature.label = m_chordnames[chordpath[iFrame]];
                fsOut[m_outputChords].push_back(chord_feature);
                oldchord = chordpath[iFrame];         
            }
            /* calculating simple chord change prob */            
            for (int iChord = 0; iChord < nChord; iChord++) {
                chordchange[iFrame-1] += delta[(iFrame-1)*nChord + iChord] * log(delta[(iFrame-1)*nChord + iChord]/delta[iFrame*nChord + iChord]);
            }
        }
        
        // cerr << chordpath[0] << endl;
	} else {
        /* Simple chord estimation
           I just take the local chord estimates ("currentChordSalience") and average them over time, then
           take the maximum. Very simple, don't do this at home...
        */
        cerr << "[Chordino Plugin] Simple Chord Estimation ... ";
        count = 0; 
        int halfwindowlength = m_inputSampleRate / m_stepSize;
        vector<int> chordSequence;
        for (vector<Vamp::RealTime>::iterator it = timestamps.begin(); it != timestamps.end(); ++it) { // initialise the score chordogram
            vector<int> temp = vector<int>(nChord,0);
            scoreChordogram.push_back(temp);
        }
        for (vector<Vamp::RealTime>::iterator it = timestamps.begin(); it < timestamps.end()-2*halfwindowlength-1; ++it) {		
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
            if (bestchordL != bestchordR) {
                chordchange[maxindex+count] += (halfwindowlength - abs(maxindex-halfwindowlength)) * 2.0 / halfwindowlength;
            }
            count++;	
        }
        // cerr << "*******  agent finished   *******" << endl;
        count = 0;
        for (vector<Vamp::RealTime>::iterator it = timestamps.begin(); it != timestamps.end(); ++it) { 
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
            count++;
        }


        // mode filter on chordSequence
        count = 0;
        string oldChord = "";
        for (vector<Vamp::RealTime>::iterator it = timestamps.begin(); it != timestamps.end(); ++it) {
            Feature chord_feature; // chord estimate
            chord_feature.hasTimestamp = true;
            chord_feature.timestamp = *it;
            // Feature currentChord; // chord estimate
            // currentChord.hasTimestamp = true;
            // currentChord.timestamp = currentChromas.timestamp;

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
            // cerr << chordchange[count] << endl;            
            if (oldChord != maxChord) {
                oldChord = maxChord;
                chord_feature.label = m_chordnames[maxChordIndex];
                fsOut[m_outputChords].push_back(chord_feature);
            }
            count++;
        }
    }
    Feature chord_feature; // last chord estimate
    chord_feature.hasTimestamp = true;
    chord_feature.timestamp = timestamps[timestamps.size()-1];
    chord_feature.label = "N";
    fsOut[m_outputChords].push_back(chord_feature);
    cerr << "done." << endl;
    
    for (int iFrame = 0; iFrame < nFrame; iFrame++) {
        Feature chordchange_feature;
        chordchange_feature.hasTimestamp = true;
        chordchange_feature.timestamp = timestamps[iFrame];
        chordchange_feature.values.push_back(chordchange[iFrame]);
        // cerr << chordchange[iFrame] << endl;
        fsOut[m_outputHarmonicChange].push_back(chordchange_feature);
    }
    
    // for (int iFrame = 0; iFrame < nFrame; iFrame++) cerr << fsOut[m_outputHarmonicChange][iFrame].values[0] << endl;
    
    
    return fsOut;     
}
