# AntelopeIT_Astrolib

A library of astronomy functions written for the Arduino platform adapted from the [IDLAstro](https://idlastro.gsfc.nasa.gov/contents.html) library and the book Astronomical Algorithms (second edition) by Jean Meeus.

## Overview

The library is a small subset of functions from the IDLAstro library adapted for use on an Arduino for incorporation into a specific, larger project. Necessarily the implementation has been tailoried to meet requirements of the project and modified to accomodate the limitations of the
Arduino platform. The IDLAstro library features more capable implementations than was required in this project e.g. a function in this library handles a single set of data whereas the IDLAstro version of the function will process an array of data sets at a time. 

The examples have been written to demonstrate how the library can be used and also as a check to ensure that the functions return the same results as predicted by Meeus in his book Astronomical Algorithms. The book seems to serve as the inspiration for the IDLAstro library as well being referenced in the code examples for that library.

The library has been incorporated into a larger astronomical hardware [project](https://github.com/Antelope-IT/StarPointer). For reference the project uses an [Arduino MKR1010](https://docs.arduino.cc/hardware/mkr-wifi-1010) and appears to be working well. 

Although the library has been written primarily for use on the Arduino Platform there is nothing in its implementation specific to the Arduino. As such it should be posible to extract the `.cpp` and `.h` files for use in your own project.

## Documentation

The library has been documented using Doxygen by following this [guide](https://learn.adafruit.com/the-well-automated-arduino-library/doxygen) from [Adafruit](https://learn.adafruit.com/).

The project documentation can be found [here](https://antelope-it.github.io/AntelopeIT_Astrolib/).

## How To Use

1. Download the AntelopeIT_Astrolib.zip file from this repository
2. Extract the contents of the zip file to the directory AntelopeIT_Astrolib in the `libraries` sub-directory of your sketchbook location `home\<user>\Arduino\libraries`  
3. Add the following to your sketch

```C
#include <AntelopeIT_Astrolib.h>
```

It should now be possible to load the example files in the Arduino IDE and call the Astrolib functions from within your own code 

## Health Warning

You are welcome to use what you see here and to learn from my mistakes but I make no guarantees regards the quality, functionality or suitability of this project for your application. You are free to use and adapt the information, ideas and designs contained in this project but you do so entirely at your own risk.
