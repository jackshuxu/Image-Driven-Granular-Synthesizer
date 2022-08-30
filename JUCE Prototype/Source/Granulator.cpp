/*
  ==============================================================================

  ==============================================================================
*/
#include "Granulator.h"

Granulator::Granulator()
{
    this->activeGrains = juce::Array<Grain*>();
    this->nextOnset = 0;
    this->position = 0;
    this->processorSampleRate = 0;
    this->portionLength = 0;
    this->hiPass.setType(juce::dsp::LinkwitzRileyFilterType::highpass);



    this->isPlaying = false;
    this->hasLoadedFile = false;
    this->filePos = 0;
    this->envIndex = 1;
    this->envWidth = 0.125;
    this->sectionSize = 0.5;
    this->density = 25.;
    this->grainSize = 25.;
    this->speedModule = 1;
    this->sampleRate = 0;
    this->fileLength = 0;
    this->speedDirection = 1;
    this->readposition = 0;
    this->init = false;
    this->currentGain = 0.75;
    this->spread = 0;
    this->highCutFreq = 19999.f;
    this->lowCutFreq = 1.f;
    this->pitch = 0.f;
    this->num_samples = 0;
    this->wave_len = 0;
}

Granulator::~Granulator()
{
   /* for (auto grain : this->activeGrains) {       
        delete grain;
        grain = nullptr;
    }*/
    for (int i = 0; i < this->activeGrains.size(); ++i)
    {
        auto g = this->activeGrains.removeAndReturn(i);
        delete g;
    }
}

//Initialize everything, add the first grain and take the onSet time for the next grain. Called when PLAY button is clicked
void Granulator::initialize(int portionLength)
{
    for (int i = 0; i < this->activeGrains.size(); ++i)
    {
        auto g = this->activeGrains.removeAndReturn(i);
        delete g;
    }

    this->position =getFilePos();
    setReadPosition(0); 

    this->activeGrains.add(new Grain(getGrainSize(), 
                                     this->position, 
                                     false, 
                                     getPitch(),
                                     getEnvIndex(),
                                     getEnvWidth(),
                                     this->sampleRate,
                                     getSpeedDirection(),&this->hiPass,
                                     hilbertTransform,
                                     ceiledLength));
    
    this->nextOnset = round(this->sampleRate/getDensity());
    this->portionLength = portionLength;


}


// Process the sound
void Granulator::process(juce::AudioBuffer<float>& outputBuffer, int numSamples)
{
    for (int samplePos = 0; samplePos < numSamples; samplePos++) {
        float toExtract = 0;

        
        if (this->activeGrains.isEmpty()) {
            for (int i = 0; i < outputBuffer.getNumChannels(); i++) {
                outputBuffer.addSample(i, samplePos, 0);

            }
        }
        else 
        {
            //printf("==============%i============", this->activeGrains.size());
            for (auto grain : this->activeGrains) 
            {
                for (int i = 0; i < outputBuffer.getNumChannels(); i++) 
                {
                    outputBuffer.addSample(i, samplePos, grain->getCurrentSample(i));   
                    toExtract += grain->getCurrentSample(i);
                    
                }
                grain->updateIndex();
                if (grain->isFinished()) 
                {
                    this->activeGrains.remove(this->activeGrains.indexOf(grain));
                    delete grain;
                }
            }            
        }

        
        this->position++;


        this->nextOnset--;
        if (this->nextOnset == 0) {
            int sectionSize = juce::jmin(getSectionSize(), 
                get_wave_len() - (int)getGrainSize() -getFilePos());
            sectionSize = juce::jmax(sectionSize, (int)getGrainSize()); 
            int timePassed = this->position - getFilePos(); 
            float readPositionShift = getSpeedModule() * timePassed;  
            int circularShift = (int)readPositionShift % sectionSize; 
            int readPosition = getFilePos() + circularShift * getSpeedDirection(); //initial position + shift * direction

            if (getSpeedDirection() < 0) 
                readPosition += sectionSize;

            if (readPosition == getFilePos()) 
                this->position = getFilePos();

            setReadPosition(readPosition);
            readPosition = getCurrentTime();
            if(readPosition > get_wave_len() - getGrainSize())
                readPosition = readPosition + getFilePos() - get_wave_len();
           
            juce::dsp::ProcessSpec spec{ this->sampleRate, static_cast<juce::uint32> (getGrainSize()), get_num_channels() };
            this->hiPass.prepare(spec);

            if (randomize()) {
                readPosition += getSpread();
                if (readPosition < 0)
                    readPosition += get_wave_len();
                if (readPosition > get_wave_len() - getGrainSize()) {
                    readPosition += getGrainSize();
                    readPosition = readPosition % get_wave_len();
                }
            }
            
            setRealPosition(readPosition);
   
            Grain* toAdd = new Grain(getGrainSize(),
                                     readPosition,
                                     false,
                                     getPitch(),
                                     getEnvIndex(),
                                     getEnvWidth(),
                                     this->sampleRate,
                                     getSpeedDirection(),&this->hiPass,
                                     hilbertTransform,
                                     ceiledLength);

            if (!this->activeGrains.isEmpty()) 
            {
                int crossfade = this->activeGrains.getLast()->remainingLife();
                toAdd->applyCrossFade(crossfade, true);
                this->activeGrains.getLast()->applyCrossFade(crossfade, false);               
            }
            this->nextOnset = this->sampleRate / getDensity();
            this->activeGrains.add(toAdd);             
        }
    }
}

void Granulator::setProcessorSampleRate(double processorSampleRate)
{
    this->processorSampleRate = processorSampleRate;
}

int Granulator::getFilePos()
{
    return round(this->filePos * this->fileLength);
}

float Granulator::getEnvWidth()
{
    return this->envWidth;
}

int Granulator::getSectionSize()
{
    return round(this->sectionSize * this->fileLength);
}

float Granulator::getDensity()
{
    return this->density;
}

int Granulator::getEnvIndex()
{
    return this->envIndex;
}

int Granulator::getGrainSize()
{
    return pow(10, -3) * this->grainSize * this->sampleRate;
}

float Granulator::getSpeedModule()
{
    return this->speedModule;
}

bool Granulator::getIsPlaying()
{
    return this->isPlaying;
}

void Granulator::setIsPlaying(bool val)
{
    this->isPlaying = val;
}

void Granulator::setHasLoadedFile(bool hasDone, int fileLength)
{
    this->hasLoadedFile = hasDone;
    this->fileLength = fileLength;
}

bool Granulator::getHasLoadedFile()
{
    return this->hasLoadedFile;
}

void Granulator::setSampleRate(double sampleRate)
{
    this->sampleRate = sampleRate;
}

int Granulator::getSpeedDirection()
{
    return this->speedDirection;
}

bool Granulator::getInit()
{
    return this->init;
}

void Granulator::setInit(bool val)
{
    this->init = val;
}

int Granulator::getReadPosition()
{
    return this->readposition;
}

int Granulator::getRealPosition()
{
    return this->realPosition;
}

void Granulator::setRealPosition(int newPos)
{
    this->realPosition = newPos;
}

void Granulator::setReadPosition(int readPosition)
{
    this->readposition = readPosition;
}

juce::Array<std::pair<float, float>>* Granulator::getxyPlane()
{
    return &xyPlane;
}

int Granulator::getxyArrayPosition()
{
    if (!getHasLoadedFile() || !getIsPlaying())
        return 0;

    float position = (float)abs(readposition - filePos * get_wave_len()) * //value to map
        (float)xyPlane.size() /  //new range
        ((sectionSize) * (float)get_wave_len()); //old range
    return position;
}

std::pair<float, float> Granulator::getCurrentxyPosition()
{
    int pos = getxyArrayPosition();
    auto ret = xyPlane[pos];
    return ret;
}

float Granulator::getCurrentFrequencyShift()
{
    if (xyPlane.isEmpty())
        return 0;
    int pos = getxyArrayPosition();
    float freq = juce::jmax(0.0f, juce::jmin(xyPlane[pos].first, 1.0f)) * 400.0f - 200.0f;
    return freq;
}

int Granulator::getCurrentTime()
{
    if (xyPlane.isEmpty())
        return this->readposition;

    int pos = getxyArrayPosition();
    float currentTime = (this->filePos * get_wave_len())
        + juce::jmax(0.0f, juce::jmin(xyPlane[pos].first, 1.0f)) * this->sectionSize
        * get_wave_len();
    return (int)juce::jmin(currentTime, (float)fileLength - 1); // in samples
}

float Granulator::getCurrentGain()
{
    return this->currentGain;
}

int Granulator::getSpread()
{
    if (this->spread == 0) return 0;
    auto r = r_spread.nextInt(juce::Range<int>(-this->spread, this->spread));
    //printf("==============%i==============", this->spread);
    return r;
}

bool Granulator::randomize()
{
    return ((r_random.nextFloat() * 100) < this->random);
}


int Granulator::getSampleRate()
{
    return this->sampleRate;
}

void Granulator::setSampleRate(int sampleRate)
{
    this->sampleRate = sampleRate;
}