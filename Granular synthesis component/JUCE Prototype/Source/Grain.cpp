/*
  ==============================================================================

  ==============================================================================
*/

#include "Grain.h"

Grain::Grain(int grainDuration, int startPos, bool highreSolution, 
             float freqShift, int envelopeType, float envelopeWidth, 
             int hostRate, int direction, juce::dsp::LinkwitzRileyFilter<float>* hiPass,
             double* hilbert, int ceiledLen) :
    length(grainDuration), startPosition(startPos)
{

    this->buffer = new juce::AudioBuffer<float>(1, this->length);
    this->hilbertTransform = hilbert;
    this->ceiledLength = ceiledLen;
    this->currentPosition = 0;
    this->finished = false;

    if (freqShift > 0) { 
    
        juce::ScopedPointer<juce::AudioBuffer<float>> analiticSignal = new juce::AudioBuffer<float>(1, this->length);
        for (auto channel = 0; channel < 1; channel++) {
            for (auto i = 0; i < length; i++) {
                analiticSignal->setSample(channel, i, hilbertTransform[bufferHilbertIndex(channel, i) + 1]);
                this->buffer->setSample(channel, i, hilbertTransform[bufferHilbertIndex(channel, i)]);
            }
        }
        hiPass->reset();
        hiPass->setType(juce::dsp::LinkwitzRileyFilterType::highpass);
        juce::dsp::ProcessSpec spec{ hostRate, static_cast<juce::uint32> (grainDuration), static_cast<juce::uint32> (1)};
        hiPass->prepare(spec);
        hiPass->setCutoffFrequency(freqShift + hostRate / (2 * length));
        juce::dsp::AudioBlock<float> block(*this->buffer);
        juce::dsp::ProcessContextReplacing<float> context(block);
        hiPass->process(context);

        juce::dsp::AudioBlock<float> analiticBlock(*analiticSignal);
        juce::dsp::ProcessContextReplacing<float> analiticContext(analiticBlock);
        hiPass->process(analiticContext);

        for (int i = 0; i < 1; i++)
            channelFreqShift(freqShift, i, envelopeType, envelopeWidth, hostRate, analiticSignal);
    }
    else {
        for (int i = 0; i < 1; i++)
            channelFreqShift(freqShift, i, envelopeType, envelopeWidth, hostRate);

    }

    //If the speed is negative reverse the buffer
    if (direction == -1)
        this->buffer->reverse(0, this->length);

}

Grain::~Grain()
{
    delete buffer;
}

void Grain::channelFreqShift(float freqShift, int channel, int envType, float envWidth, int hostRate)
{
    for (int i = 0; i < length; i++) {//freq shift  --->    ref links
        float phaseInc = freqShift * i / hostRate;
        phaseInc -= ceil(phaseInc);

        float theta = TWOPI * phaseInc; //angle
        float newValue = hilbertTransform[bufferHilbertIndex(channel, i)] * juce::dsp::FastMathApproximations::cos(theta) -
            hilbertTransform[bufferHilbertIndex(channel, i) + 1] * juce::dsp::FastMathApproximations::sin(theta); //rotation
        newValue *= GrainEnvelope::getEnvelopeValue(i, envType, this->length, envWidth);
        this->buffer->setSample(channel, i, newValue); //rewrite channel

    }
}

void Grain::channelFreqShift(float freqShift, int channel, int envType, float envWidth, int hostRate, juce::AudioBuffer<float>* analiticSignal)
{
    for (int i = 0; i < length; i++) {//freq shift  --->    ref links
        float phaseInc = freqShift * i / hostRate;
        phaseInc -= floor(phaseInc); //handle phase

        float theta = TWOPI * phaseInc; //angle
        float newValue = this->buffer->getSample(channel, i) * cos(theta) -
            analiticSignal->getSample(channel, i) * sin(theta); //rotation
        newValue *= GrainEnvelope::getEnvelopeValue(i, envType, this->length, envWidth);
        this->buffer->setSample(channel, i, newValue); //rewrite channel
    }
}


float Grain::getCurrentSample(int channel)
{
    return buffer->getSample(channel % 1, this->currentPosition);
}

void Grain::updateIndex()
{
    this->currentPosition++;
    if (this->currentPosition == this->length) {
        this->finished = true;
    }
}

int Grain::bufferHilbertIndex(int channel, int index)
{
   return 2 * (channel * ceiledLength + startPosition + index);
}

bool Grain::isFinished()
{
    return this->finished;
}

juce::AudioBuffer<float>* Grain::getBuffer()
{
    return this->buffer;
}

void Grain::applyCrossFade(int crossfade, bool atStart)
{
    if (atStart) {     
        crossfade = std::min(crossfade, this->length - 1);
        this->buffer->applyGainRamp(0, crossfade, 0, 1);
    }
    else {
        this->buffer->applyGainRamp(this->length - crossfade, crossfade, 1, 0);
    }
}

int Grain::remainingLife()
{
    return this->length - this->currentPosition;
}

