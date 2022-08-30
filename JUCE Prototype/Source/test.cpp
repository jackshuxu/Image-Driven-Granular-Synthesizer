/*
  ==============================================================================

    test.cpp
    Created: 3 Aug 2022 8:24:12pm
    Author:  tangjiwen

  ==============================================================================
*/

#include "test.h"

int getNum()
{
    juce::Array<int> test_array = { 1, 1 };
    test_array.add(1);
    return test_array[0] + test_array[1] + test_array[2];
}


Granulator* granulator = new Granulator();
void initial_gran(float* input,
                  double* sample_rate,
                  long long* wave_len, 
                  float* file_pos, 
                  float* section_size, 
                  float* env_index, 
                  float* env_width, 
                  float* density, 
                  float* gran_size, 
                  float* speed, 
                  float* gain, 
                  float* random, 
                  float* random_spread, 
                  float* pitch)
{
    int sr = *sample_rate;
    int size = *wave_len;
    int spread = *random_spread;
    granulator->sampleRate = sr;
    granulator->wave_len = size;
    granulator->filePos = *file_pos;
    granulator->sectionSize = *section_size;
    granulator->envIndex = *env_index;
    granulator->envWidth = *env_width;
    granulator->density = *density;
    granulator->grainSize = *gran_size;
    granulator->speedModule = abs(*speed);
    granulator->speedDirection = *speed == 0 ? 0 : (*speed / abs(*speed));
    granulator->currentGain = *gain;
    granulator->random = *random;
    granulator->spread = spread * size / 2 * 0.5;
    granulator->pitch = *pitch;

    juce::AudioBuffer<float> buffer(&input, 1, size);

    int ceiledLength = pow(2, ceil(log2(buffer.getNumSamples())));
    double* hilbertTransform = (double*)calloc((size_t)(buffer.getNumChannels() * (size_t)2 * ceiledLength), sizeof(double));

    for (int i = 0; i < buffer.getNumChannels(); i++)
    {
        for (int j = 0; j < buffer.getNumSamples(); j++)
        { //apply envelope
            if (hilbertTransform != NULL)
            {
                hilbertTransform[i * 2 * ceiledLength + j * 2] = buffer.getSample(i, j);
                hilbertTransform[i * 2 * ceiledLength + j * 2 + 1] = 0; //real signal ----> a value every two set to zero
                //printf("==========%f==========", buffer.getSample(i, j));
            }
        }

        if (hilbertTransform != NULL)
            hilbert(&hilbertTransform[i * ceiledLength], ceiledLength); //Transform file
    }

    granulator->hilbertTransform = hilbertTransform;
    granulator->ceiledLength = ceiledLength;
    granulator->setInit(true);
    granulator->setIsPlaying(true);
    granulator->setHasLoadedFile(true, size);

    granulator->initialize(size);
}

void release_gran()
{
    delete granulator;
}

float* juce_process(float* output, int* sample_num)
{
    //if (*sample_num) do
    //{
    //    *output++ = (float)(0.01f * ((float)rand() / (float)RAND_MAX * 2.f - 1.f) + *(input + *input_sample_cnt));
    //    *input_sample_cnt += 1;
    //    if (*input_sample_cnt >= *wave_len)
    //    {
    //        *input_sample_cnt = 0;
    //    }
    //} while((*sample_num)--);
    juce::AudioSampleBuffer buffer(&output, 1, *sample_num);
    buffer.clear();

    granulator->process(buffer, *sample_num);
    buffer.applyGain(juce::Decibels::decibelsToGain(granulator->getCurrentGain()));
    
    return output;
}
