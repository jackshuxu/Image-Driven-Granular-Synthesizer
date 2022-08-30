/*
  ==============================================================================

    test.h
    Created: 3 Aug 2022 8:24:04pm
    Author:  tangjiwen

  ==============================================================================
*/

#pragma once
#include "JuceHeader.h"
#include "Granulator.h"

extern "C" __declspec(dllexport) void release_gran();
extern "C" __declspec(dllexport) void initial_gran(float*,
                                                   double*,
                                                   long long*,
                                                   float*,
                                                   float*,
                                                   float*,
                                                   float*,
                                                   float*,
                                                   float*,
                                                   float*, 
                                                   float*, 
                                                   float*, 
                                                   float*,
                                                   float*);

extern "C" __declspec(dllexport) float* juce_process(float*, int*);

extern "C" __declspec(dllexport) int getNum();

