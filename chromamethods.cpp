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

#include "chromamethods.h"

#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <boost/tokenizer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/lexical_cast.hpp>

#include "chorddict.cpp"

using namespace std;
using namespace boost;


/** Special Convolution
    special convolution is as long as the convolvee, i.e. the first argument. in the valid core part of the 
    convolution it contains the usual convolution values, but the pads at the beginning (ending) have the same values
    as the first (last) valid convolution bin.
**/

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
    int oversampling = 80;
	
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

    int nFFT = fft_f.size();
	
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
            if (cq_f[iCQ] * pow(2.0, 0.084) + fft_width > fft_f[iFFT] && cq_f[iCQ] * pow(2.0, -0.084 * 2) - fft_width < fft_f[iFFT]) { // within a generous neighbourhood
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

void dictionaryMatrix(float* dm) {
    int binspersemitone = 3; // this must be 3
    int minoctave = 0; // this must be 0
    int maxoctave = 7; // this must be 7
    float s_param = 0.7;
	
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

    float curr_f;
    float floatbin;
    float curr_amp;
    // now for every combination calculate the matrix element
    for (unsigned iOut = 0; iOut < 12 * (maxoctave - minoctave); ++iOut) {
        // cerr << iOut << endl;
        for (unsigned iHarm = 1; iHarm <= 20; ++iHarm) {
            curr_f = 440 * pow(2,(minMIDI-69+iOut)*1.0/12) * iHarm;
            // if (curr_f > cq_f[nNote-1])  break;
            floatbin = ((iOut + 1) * binspersemitone + 1) + binspersemitone * 12 * log2(iHarm);
            // cerr << floatbin << endl;
            curr_amp = pow(s_param,float(iHarm-1));
            // cerr << "curramp" << curr_amp << endl;
            for (unsigned iNote = 0; iNote < nNote; ++iNote) {
                if (abs(iNote+1.0-floatbin)<2) {
                    dm[iNote  + 256 * iOut] += cospuls(iNote+1.0, floatbin, binspersemitone + 0.0) * curr_amp;
                    // dm[iNote + nNote * iOut] += 1 * curr_amp;
                }
            }
        }
    }


}

static
std::vector<std::string>
getPluginPath()
{
    //!!! This is duplicated from PluginHostAdapter::getPluginPath,
    //!!! which is not available to us in the plugin (only to the
    //!!! host)

    std::vector<std::string> path;
    std::string envPath;

    char *cpath = getenv("VAMP_PATH");
    if (cpath) envPath = cpath;

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#define DEFAULT_VAMP_PATH "%ProgramFiles%\\Vamp Plugins"
#else
#define PATH_SEPARATOR ':'
#ifdef __APPLE__
#define DEFAULT_VAMP_PATH "$HOME/Library/Audio/Plug-Ins/Vamp:/Library/Audio/Plug-Ins/Vamp"
#else
#define DEFAULT_VAMP_PATH "$HOME/vamp:$HOME/.vamp:/usr/local/lib/vamp:/usr/lib/vamp"
#endif
#endif

    if (envPath == "") {
        envPath = DEFAULT_VAMP_PATH;
        char *chome = getenv("HOME");
        if (chome) {
            std::string home(chome);
            std::string::size_type f;
            while ((f = envPath.find("$HOME")) != std::string::npos &&
                    f < envPath.length()) {
                envPath.replace(f, 5, home);
            }
        }
#ifdef _WIN32
        char *cpfiles = getenv("ProgramFiles");
        if (!cpfiles) cpfiles = (char *)"C:\\Program Files";
        std::string pfiles(cpfiles);
        std::string::size_type f;
        while ((f = envPath.find("%ProgramFiles%")) != std::string::npos &&
               f < envPath.length()) {
            envPath.replace(f, 14, pfiles);
        }
#endif
    }

    std::string::size_type index = 0, newindex = 0;

    while ((newindex = envPath.find(PATH_SEPARATOR, index)) < envPath.size()) {
	path.push_back(envPath.substr(index, newindex - index));
	index = newindex + 1;
    }
    
    path.push_back(envPath.substr(index));

    return path;
}

vector<string> chordDictionary(vector<float> *mchorddict) {

    typedef tokenizer<char_separator<char> > Tok;
    char_separator<char> sep(",; ","=");

    string chordDictBase("chord.dict");
    string chordDictFilename;

    vector<string> ppath = getPluginPath();
    for (int i = 0; i < ppath.size(); ++i) {
	chordDictFilename = ppath[i] + "/" + chordDictBase;
	cerr << "Looking for chord.dict in " << chordDictFilename << "..." << endl;
	if (iostreams::stream<iostreams::file_source>(chordDictFilename.c_str())
	    .is_open()) {
	    cerr << "(Success)" << endl;
	    break;
	}
    }

    iostreams::stream<iostreams::file_source> chordDictFile(chordDictFilename);
    string line;
    int iElement = 0;
    int nChord = 0;
	
    vector<string> loadedChordNames;
    vector<float> loadedChordDict;
    if (chordDictFile.is_open()) {
        while (std::getline(chordDictFile, line)) { // loop over lines in chord.dict file		
            // first, get the chord definition
            string chordType;
            vector<float> tempPCVector;			
            // cerr << line << endl;
            if (!line.empty() && line.substr(0,1) != "#") {
                Tok tok(line, sep);			
                for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter) { // loop over line elements
                    string tempString = *tok_iter;
                    // cerr << tempString << endl;
                    if (tok_iter == tok.begin()) { // either the chord name or a colon
                        if (tempString == "=") {
                            chordType = "";
                        } else {
                            chordType = tempString;
                            tok_iter++; // is this cheating ? :)
                        }
                    } else {
                        tempPCVector.push_back(lexical_cast<float>(*tok_iter));
                    }
                }
					
                // now make all 12 chords of every type
                for (unsigned iSemitone = 0; iSemitone < 12; iSemitone++) {				
                    // add bass slash notation
                    string slashNotation = "";
                    for (unsigned kSemitone = 1; kSemitone < 12; kSemitone++) {
                        if (tempPCVector[(kSemitone) % 12] > 0.99) {
                            slashNotation = bassnames[iSemitone][kSemitone];
                        }
                    }
                    for (unsigned kSemitone = 0; kSemitone < 12; kSemitone++) { // bass pitch classes
                        // cerr << ((kSemitone - iSemitone + 12) % 12) << endl;
                        float bassValue = 0;
                        if (tempPCVector[(kSemitone - iSemitone + 12) % 12]==1) {
                            bassValue = 1;
                        } else {
                            if (tempPCVector[((kSemitone - iSemitone + 12) % 12) + 12] == 1) bassValue = 0.5;
                        }
                        loadedChordDict.push_back(bassValue);
                    }
                    for (unsigned kSemitone = 0; kSemitone < 12; kSemitone++) { // chord pitch classes
                        loadedChordDict.push_back(tempPCVector[((kSemitone - iSemitone + 12) % 12) + 12]);
                    }
                    ostringstream os;				
                    if (slashNotation.empty()) {
                        os << notenames[12+iSemitone] << chordType;
                    } else {
                        os << notenames[12+iSemitone] << chordType << "/" << slashNotation;
                    }
                    // cerr << os.str() << endl;
                    loadedChordNames.push_back(os.str());
                }
            }
        }
        // N type
        loadedChordNames.push_back("N");
        for (unsigned kSemitone = 0; kSemitone < 12; kSemitone++) loadedChordDict.push_back(0.5);
        for (unsigned kSemitone = 0; kSemitone < 12; kSemitone++) loadedChordDict.push_back(1.0);
	
        // normalise
        float sum = 0;
        for (int i = 0; i < loadedChordDict.size(); i++) {
            sum += pow(loadedChordDict[i],2);
            if (i % 24 == 23) {
                float invertedsum = 1.0/sqrt(sum);
                for (int k = 0; k < 24; k++) {
                    loadedChordDict[i-k] *= invertedsum; 
                }
                sum = 0;
            }
		
        }
	

        nChord = 0;
        for (int i = 0; i < loadedChordNames.size(); i++) {
            nChord++;
        }
        chordDictFile.close();


        // mchorddict = new float[nChord*24];
        for (int i = 0; i < nChord*24; i++) {
            mchorddict->push_back(loadedChordDict[i]);			
        }
			
    } else {// use default from chorddict.cpp
        // mchorddict = new float[nChorddict];
        for (int i = 0; i < nChorddict; i++) {
            mchorddict->push_back(chorddict[i]);
        }
		
        nChord = nChorddict/24;
        // mchordnames = new string[nChorddict/24];
        char buffer1 [50];
        for (int i = 0; i < nChorddict/24; i++) {
            if (i < nChorddict/24 - 1) {
                sprintf(buffer1, "%s%s", notenames[i % 12 + 12], chordtypes[i]);
            } else {
                sprintf(buffer1, "N");
            }
            ostringstream os;
            os << buffer1;
            loadedChordNames.push_back(os.str());

        }
		
    }
    // cerr << "before leaving" << chordnames[1] << endl;
    return loadedChordNames;
}
