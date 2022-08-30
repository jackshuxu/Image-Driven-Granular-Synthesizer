/*
  ==============================================================================

  ==============================================================================
*/

#pragma once
#include"JuceHeader.h"
#include "math_const.h"


class GrainEnvelope
{

public:
	
    static float getEnvelopeValue(int index, int selectedEnv, int length, float mainlobewidth); //sample index, selected envelope type, grain length, main lobe width  
    
protected:    

private:
    static float getGaussian(int index, int length, float mainLobeWidth);
    static float getTrapezoidal(int index, int length, float mainLobeWidth);
    static float getRaisedCosine(int index, int length, float mainLobeWidth);
    
};

