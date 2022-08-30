/*
  ==============================================================================

  ==============================================================================
*/

#pragma once
#include<JuceHeader.h>
#include <cmath>

#include "math_const.h"

#ifdef __cplusplus 
extern "C" {
#endif
#include "transforms.h"

#ifdef __cplusplus 
}
#endif

/*
Simpson integrator: integrates time and frequency using "Simpson rule".
Just create one and use the two public get methods to get average time and frequency
*/

class SimpsonIntegrator {
public:

    SimpsonIntegrator(double* hilbertTransform, int samplingFrequency, int length, int numChannels, int notCeiledLength);
    SimpsonIntegrator(double* hilbertTransform, int samplingFrequency, int length, int numChannels, int notCeiledLength, float freqshift); //integrate shifted signal
    ~SimpsonIntegrator();
    float getAverageFrequency();
    float getAverageTime();

private:

    void computeAverageFrequency(double* hilbertTransform);
    void computeAverageTime(double* hilbertTransform);
    void computeAverageFrequency(double* hilbertTransform, float freqShift);
    void computeAverageTime(double* hilbertTransform, float freqShift);

    int samplingFrequency;
    int length;
    int notCeiledLength;
    int numChannels;

    float averageFrequency;
    float averageTime;

    double* hilbertSpectrum;
};