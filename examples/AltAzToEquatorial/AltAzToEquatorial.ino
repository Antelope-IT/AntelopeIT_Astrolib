/*****************************************************************************
 *  Antelope IT Astronomy Library - AltAzToEquatorial Example.
 * 
 *  A simple example to demonstrate how the library can be used
 *  to convert a set of Alt-Azimuth values to their Equatorial
 *  equivalent.
 * 
 *  Written by Mark Hargrove for Antelope IT Ltd.
 * 
 *  GPLv3 license, all text above must be included in any redistribution
 * 
 *****************************************************************************/

#include <Arduino.h>
#include <AntelopeIT_Astrolib.h>

void printArray (String name, double arr[],int length) {
  Serial.print(name + ":, ");
  for (auto i = 0; i< length; i++) {
    Serial.print(String(arr[i], 4) +  ", ") ; 
  }
  Serial.println();
}

void setup() {
    while (!Serial);        // wait until Serial is ready
    Serial.begin(115200);   // connect at 115200

    Serial.println("============  AltAzToEquatorial Example.  ============");
    JulianDate jd25Dec2041 = convertToJulianDate2000(2041,12,26,5,0,0.0);
    Equatorial starCoords = convertToEquatorial(37.91138889, 264.9183333,jd25Dec2041, 31.963, -111.6, 2096,273.0,781.0);
    DegMinSecs dmsRa = toDMS(starCoords.ra / hoursToDegrees);
    DegMinSecs dmsDec = toDMS(starCoords.dec);
    Serial.println("Calculated RA: " + String(starCoords.ra / hoursToDegrees,10) + " => " + String(dmsRa.degrees) +"Hrs " + String(dmsRa.minutes) + "mins " + String(dmsRa.seconds,2)+" sec." );
    Serial.println("Expected   RA: 0.2206111111 => 0Hrs 13mins 14.20 sec." );
    Serial.println();
    Serial.println("Calculated DEC: " + String(starCoords.dec,10) + " => " + String(dmsDec.degrees) + "deg " + String(dmsDec.minutes) + "mins " + String(dmsDec.seconds,2) + " sec.");
    Serial.println("Expected   DEC: 15.18361111   => 15deg 11mins  1.00 sec.");
    Serial.println("============ End AltAzToEquatorial Example.============");
    Serial.println();
}

void loop() {

}