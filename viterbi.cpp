
#include "viterbi.h"
#include <iostream>

std::vector<int> ViterbiPath(std::vector<double> init, std::vector<vector<double> > trans, std::vector<vector<double> > obs, double *delta, vector<double> *scale) {
    
    unsigned nState = init.size();
    unsigned nFrame = obs.size();
    
    // check for consistency
    if (trans[0].size() != nState || trans.size() != nState || obs[0].size() != nState) {
        cerr << "ERROR: matrix sizes inconsistent." << endl;
    }

    vector<vector<int> > psi; //  "matrix" of remembered indices of the best transitions
    vector<int> path = vector<int>(nFrame, nState-1); // the final output path (current assignment arbitrary, makes sense only for Chordino, where nChord-1 is the "no chord" label)
    
    double deltasum = 0;
    
    /* initialise first frame */
    for (unsigned iState = 0; iState < nState; ++iState) {
        delta[iState] = init[iState] * obs[0][iState];
        deltasum += delta[iState];
    }
    for (unsigned iState = 0; iState < nState; ++iState) delta[iState] /= deltasum; // normalise (scale)
    scale->push_back(1.0/deltasum);
    psi.push_back(vector<int>(nState,0));
    
    /* rest of the forward step */
    for (unsigned iFrame = 1; iFrame < nFrame; ++iFrame) {
        deltasum = 0;
        psi.push_back(vector<int>(nState,0));
        /* every state wants to know which previous state suits it best */
        for (unsigned jState = 0; jState < nState; ++jState) {            
            int bestState = nState - 1;
            double bestValue = 0;
            if (obs[iFrame][jState] > 0) {
                for (unsigned iState = 0; iState < nState; ++iState) {
                    double currentValue = delta[(iFrame-1) * nState + iState] * trans[iState][jState];
                    if (currentValue > bestValue) {
                        bestValue = currentValue;
                        bestState = iState;
                    }
                }
            }

            delta[iFrame * nState + jState] = bestValue * obs[iFrame][jState];
            deltasum += delta[iFrame * nState + jState];
            psi[iFrame][jState] = bestState;
        }
        if (deltasum > 0) {
            for (unsigned iState = 0; iState < nState; ++iState) {            
                delta[iFrame * nState + iState] /= deltasum; // normalise (scale)
            }
            scale->push_back(1.0/deltasum);
        } else {
            for (unsigned iState = 0; iState < nState; ++iState) {            
                delta[iFrame * nState + iState] = 1.0/nState;
            }
            scale->push_back(1.0);
        }
        
    }
    
    /* initialise backward step */
    double bestValue = 0;
    for (unsigned iState = 0; iState < nState; ++iState) {
        double currentValue = delta[(nFrame-1) * nState + iState];
        if (currentValue > bestValue) {
            bestValue = currentValue;            
            path[nFrame-1] = iState;
        }
    }

    /* rest of backward step */
    for (int iFrame = nFrame-2; iFrame > -1; --iFrame) {
        path[iFrame] = psi[iFrame+1][path[iFrame+1]];
    }    
    
    return path;
}
