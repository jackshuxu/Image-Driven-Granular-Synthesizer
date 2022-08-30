/*
  ==============================================================================

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "GrainEnvelope.h"
#include "SimpsonIntegrator.h"
#include "math_const.h"

enum class EnvType {raisedCosineBell, gaussian, trapezoidal};

class Grain {

public:
    Grain(int grainDuration, int startPos, bool highreSolution, 
          float freqShift, int envelopeType, float envelopeWidth,
          int hostRate, int direction, juce::dsp::LinkwitzRileyFilter<float> * hiPass,
          double* hilbert, int ceiledLen);
    ~Grain();
    float getCurrentSample(int channel);    //Return the current sample playing on the given channel
    void updateIndex();                     //Increment the current sample playing index. Set finish to true if the grain is finished
    bool isFinished();
    juce::AudioBuffer<float>* getBuffer();
    void applyCrossFade(int crossfade, bool atStart);
    int remainingLife();

private:
    int startPosition;      //in the loaded audio file
    int length;
    int currentPosition;    //current playing sample in the grain
    bool finished;
    double* hilbertTransform; //hilbert transform for each channel

    juce::AudioBuffer<float>* buffer;  
    void channelFreqShift(float freqShift, int channel, int envType, float envWidth, int hostRate); //shifts a channel of freqshift [Hz]
    void channelFreqShift(float freqShift, int channel, int envType, float envWidth, int hostRate, juce::AudioBuffer<float>* analiticSignal);
    int bufferHilbertIndex(int channel, int index);

    int ceiledLength;
};


