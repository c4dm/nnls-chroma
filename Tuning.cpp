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

#include "Tuning.h"

#include "chromamethods.h"

#include <cstdlib>
#include <fstream>
#include <cmath>

#include <algorithm>

const bool debug_on = false;

const vector<float> hw(hammingwind, hammingwind+19);

Tuning::Tuning(float inputSampleRate) :
    NNLSBase(inputSampleRate)
{
    if (debug_on) cerr << "--> Tuning" << endl;
}

Tuning::~Tuning()
{
    if (debug_on) cerr << "--> ~Tuning" << endl;
}

size_t 
Tuning::getPreferredStepSize() const
{
    if (debug_on) cerr << "--> getPreferredStepSize" << endl;
    return 2048*4; 
}

string
Tuning::getIdentifier() const
{
    if (debug_on) cerr << "--> getIdentifier" << endl;
    return "tuning";
}

string
Tuning::getName() const
{
    if (debug_on) cerr << "--> getName" << endl;
    return "Tuning";
}

string
Tuning::getDescription() const
{
    // Return something helpful here!
    if (debug_on) cerr << "--> getDescription" << endl;
    return "This plugin provides a number of features derived from a log-frequency amplitude spectrum of the DFT: some variants of the log-frequency spectrum, including a semitone spectrum derived from approximate transcription using the NNLS algorithm; based on this semitone spectrum, chroma features and a simple chord estimate.";
}

Tuning::ParameterList
Tuning::getParameterDescriptors() const
{
    if (debug_on) cerr << "--> getParameterDescriptors" << endl;
    ParameterList list;

    ParameterDescriptor d0;
    d0.identifier = "rollon";
    d0.name = "spectral roll-on";
    d0.description = "The bins below the spectral roll-on quantile will be set to 0.";
    d0.unit = "";
    d0.minValue = 0;
    d0.maxValue = 0.05;
    d0.defaultValue = 0;
    d0.isQuantized = true;
	d0.quantizeStep = 0.005;
    list.push_back(d0);


    return list;
}

Tuning::OutputList
Tuning::getOutputDescriptors() const
{
    if (debug_on) cerr << "--> getOutputDescriptors" << endl;
    OutputList list;
    
    int index = 0;

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
    m_outputTuning = index++;
	
    OutputDescriptor d10;
    d10.identifier = "localtuning";
    d10.name = "Local Tuning";
    d10.description = "Tuning based on the history up to this timestamp.";
    d10.unit = "Hz";
    d10.hasFixedBinCount = true;
    d10.binCount = 1;
    d10.hasKnownExtents = true;
    d10.minValue = 427.47;
    d10.maxValue = 452.89;
    d10.isQuantized = false;
    d10.sampleType = OutputDescriptor::FixedSampleRate;
    d10.hasDuration = false;
    // d10.sampleRate = (m_stepSize == 0) ? m_inputSampleRate/2048 : m_inputSampleRate/m_stepSize;
    list.push_back(d10);
    m_outputLocalTuning = index++;
  
    return list;
}


bool
Tuning::initialise(size_t channels, size_t stepSize, size_t blockSize)
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
Tuning::reset()
{
    if (debug_on) cerr << "--> reset";
    NNLSBase::reset();
}

Tuning::FeatureSet
Tuning::process(const float *const *inputBuffers, Vamp::RealTime timestamp)
{   
    if (debug_on) cerr << "--> process" << endl;

    NNLSBase::baseProcess(inputBuffers, timestamp);
	
    Feature f10; // local tuning
    f10.hasTimestamp = true;
    f10.timestamp = timestamp;
    float normalisedtuning = m_localTuning[m_localTuning.size()-1];
    float tuning440 = 440 * pow(2,normalisedtuning/12);
    f10.values.push_back(tuning440);
	
    FeatureSet fs;
    fs[m_outputLocalTuning].push_back(f10);
    return fs;	
}

Tuning::FeatureSet
Tuning::getRemainingFeatures()
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
		    
    char buffer0 [50];
		
    sprintf(buffer0, "estimated tuning: %0.1f Hz", cumulativetuning);
		    
    // push tuning to FeatureSet fsOut
    Feature f0; // tuning
    f0.hasTimestamp = true;
    f0.timestamp = Vamp::RealTime::frame2RealTime(0, lrintf(m_inputSampleRate));;
    f0.label = buffer0;
    fsOut[m_outputTuning].push_back(f0);  
		    
    return fsOut;     

}

