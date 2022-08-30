/*
  ==============================================================================

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Grain.h"


class Granulator
{
public:
    Granulator();
    ~Granulator();
    void initialize(int portionLength);
    void process(juce::AudioBuffer<float>& outputBuffer, int numSamples);
    void setProcessorSampleRate(double processorSampleRate);

    //===================================
    float getEnvWidth();
    int getSectionSize();
    float getDensity();
    int getEnvIndex();
    int getGrainSize();
    float getSpeedModule();
    bool getIsPlaying();
    void setIsPlaying(bool val);
    void setHasLoadedFile(bool hasDone, int fileLength);
    bool getHasLoadedFile();
    void setSampleRate(double sampleRate);
    int getSpeedDirection();
    bool getInit();
    void setInit(bool val);

    void setReadPosition(int readPosition);
    juce::Array<std::pair<float, float>>* getxyPlane();

    int getReadPosition();
    int getRealPosition();
    void setRealPosition(int newPos);
    int getFilePos();
    int getxyArrayPosition();
    std::pair<float, float> getCurrentxyPosition();
    float getCurrentFrequencyShift();

    int getCurrentTime();
    float getCurrentGain();
    int getSpread();
    bool randomize();

    int getSampleRate();
    void setSampleRate(int sampleRate);

    float getHighCutFreq() { return 20000.f - highCutFreq; }
    float getLowCutFreq() { return 20000.f - lowCutFreq; }

    float getPitch() { return pitch / 10.f; }

    int get_num_samples() { return num_samples; }
    int get_wave_len() { return wave_len; }
    int get_num_channels() { return 1; }


    juce::Array<Grain*> activeGrains;     
    int nextOnset;      
    int position;       
    int portionLength;
    double processorSampleRate;
    juce::dsp::LinkwitzRileyFilter<float> hiPass;

    //====================================
    juce::Random r_spread;
    juce::Random r_random;
    float filePos; // [0,1]
    float sectionSize; // [0,1]
    float envWidth;
    float density;
    int envIndex;
    float grainSize;
    float speedModule;
    bool isPlaying;
    bool hasLoadedFile;
    int fileLength;
    int speedDirection;
    int readposition;
    int realPosition;
    juce::Array<std::pair<float, float>> xyPlane;
    bool init;
    float currentGain;
    int sampleRate;

    int spread;
    float random;

    float highCutFreq;
    float lowCutFreq;

    float pitch;

    int num_samples;
    int wave_len;

    double* hilbertTransform;
    int ceiledLength;
};