usage: 1. load "testtjw.lsp"
       2. load "image" or other image data
       3. use gran-sound(filename, duration, image-data) to generate sound
gran-sound is a encapsulation of gran-t and gran-t is a encapsulation of snd-gran
before you run the gran-sound function, you should use snd-juce-import to load dsp.dll first.
after generating sound, you should use snd-juce-free to unload dll.
you can find snd-gran's definition in the header sndfnintdefs.h and dsp source code in the juce_dll_text folder.