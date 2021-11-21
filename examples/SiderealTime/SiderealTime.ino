/*****************************************************************************
 *  Antelope IT Astronomy Library - Sidereal Time Example.
 * 
 *  A simple example to demonstrate how the library can be used
 *  to convert a time and date in the gregorian calendar to its equivalent
 *  expressed in sidereal time. 
 * 
 *  The example given uses the values from example 22a on Pg.84 from 
 *  the book Astronomical Algorithms by J Meeus.
 * 
 *  Written by Mark Hargrove for Antelope IT Ltd.
 * 
 *  GPLv3 license, all text above must be included in any redistribution
 * 
 *****************************************************************************/

#include <Arduino.h>
#include <AntelopeIT_Astrolib.h>

void setup() {
    while (!Serial);        // wait until Serial is ready
    Serial.begin(115200);   // connect at 115200

    Serial.println("============     Sidereal Time Example     ============" ); 
    Serial.println("        Taken from Meeus Example 11.a Page 84.");
    Serial.println("");
    double lst = convertToSidereal(0.0,1987,4,10 ,0,0, 0.0);
    Serial.print("Calculated Mean Sidereal Time: = " );
    Serial.println(String(lst,5));  
    Serial.println("Expected Mean Sidereal Time  : = 13.17954633 Hrs => 13Hrs 10min 46.3668s."  ); 

    const JulianDate jd10Apr1987 = convertToJulianDate2000(1987,4,10 ,0, 0, 0.0); 

    Nutation nut10Apr1987 = nutate(jd10Apr1987);
    double eps10Apr1987 = calculateTrueObliquity(jd10Apr1987, nut10Apr1987);
    Serial.print("True Obliquity eps = "); 
    Serial.print(String(eps10Apr1987,5));
    Serial.println(" radians" );
    double last = convertToApparentSiderealTime(lst, nut10Apr1987, eps10Apr1987);
    Serial.print("Calculated Apparent Sidereal Time: = ");
    Serial.println(String( last, 8));
    Serial.println("Expected Apparent Sidereal Time  : = 13.17948197 Hrs => 13Hrs 10min 46.1351s."  );
    Serial.println("============  End Sidereal Time Example    ============" );
    Serial.println();
}

void loop() {

}