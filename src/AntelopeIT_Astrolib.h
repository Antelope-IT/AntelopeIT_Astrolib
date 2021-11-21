/**
 *  @file AntelopeIT_Astrolib.h
 *
 *  This is a library of Astronomy related functions inspired by the IDLAstro Library 
 *  ------> https://idlastro.gsfc.nasa.gov/contents.html
 * 
 *  and the book Astronomical Algorithms (second edition) by Jean Meeus.
 *    
 *  This library is not a complete implementation of the entire IDLAstro library or
 *  all of the algorithms in the book, it is just a selection of methods and
 *  functions to meet the needs of a particular project.
 *
 *  Written by Mark Hargrove for Antelope IT Ltd.
 * 
 *  GPLv3 license, all text above must be included in any redistribution
 */

#ifndef __ANTELOPEIT_ASTROLIB_H__
#define __ANTELOPEIT_ASTROLIB_H__

#include <math.h>

static const double sidereal_second  = 1.00273790971;       //< solar second
static const double degreesToRadians = M_PI / 180.0;
static const double secondsToRadians = degreesToRadians/3600;
static const double hoursToDegrees = 15.0;
static const double JD2000 = 2451545.0;

/**
 * @brief Struct for holding a Julian Date
 * 
 * Having converted a date to a Julian date for use in the calculations below the Julian Date is passed between 
 * functions but often what is required is the number of centuries since J2000 or occaisionally the number of days
 * this leaves each function having to calculate these two numbers for each function that requires this result.
 * Therefore as an optimisation these values are calculated once and passed with the Julian Date for use as required
 */
typedef struct JulianDate {
  double julianDate;      /**< Date expressed as expressed as the number of days since the start of the Julian Period (January 1, 4713 BC). */
  double julianCenturies; /**< Number of centuries since J2000 expressed as a decimal. */
  double julianDays;      /**< Number of days since J2000 expressed as a decimal. */
} JulianDate;

/**
 * @brief Struct for holding the result from a Nutation calculation
 */
typedef struct Nutation {
  double longitude; /**< nutation in longitude. */
  double oblique;   /**< nutation in obliquity. */ 
} Nutation;

/**
 * @brief Struct for holding the result from a calculation returning coordinates as an hour angle and declination 
 * 
 * This is identical to Equatorial but distinguishes the fact the the RA is expressed as an hour angle 
 */
typedef struct HourAngleDeclination {
  double ha;   /**< Hour angle expressed in degrees (decimal). */
  double dec;  /**< declination expressed in degrees (decimal). */
} HourAngleDeclination;

/**
 * @brief Struct for holding a coordinate value expressed in degrees minutes and seconds
 * 
 */
typedef struct DegMinSecs {
  int degrees; /**< Number of degrees 0 - 360.  */
  int minutes; /**< Number of arc minutes 0 - 60.  */
  double seconds; /**< Number of arc seconds expressed as a decimal.  */
} DegMinSecs;

/**
 * @brief Struct for holding the result from a calculation returning coordinates in Right Ascension and Declination
 * 
 */
typedef struct Equatorial {
  double ra;  /**< Right Ascension expressed in degrees (decimal). */
  double dec; /**< declination expressed in degrees (decimal). */
  /**
   * @brief Adds one Equatorial value to this value
   * 
   * @param rhs The Equatorial value to add to this.
   * @return struct Equatorial& This the Equatorial value to accumulate the value in.
   */
  struct Equatorial& operator+=(const Equatorial& rhs) { ra += rhs.ra; dec += rhs.dec; return *this; }
    /**
   * @brief Subtracts one Equatorial value to this value
   * 
   * @param rhs The Equatorial value to Subract from this.
   * @return struct Equatorial& This the Equatorial value to accumulate the value in.
   */
  struct Equatorial& operator-=(const Equatorial& rhs) { ra -= rhs.ra; dec -= rhs.dec; return *this; } 
} Equatorial;

/**
 * @brief returns the sum of two Equatorial coordinates
 * 
 * @param x  First set of coordinates to add
 * @param y  Second set of coordinates to add
 * @return Equatorial The sum of the two Equatorial co-ordinates
 */
static Equatorial operator +(const Equatorial& x, const Equatorial& y) {
  return { x.ra + y.ra, x.dec + y.dec};
}

/**
 * @brief returns the difference between two Equatorial coordinates
 * 
 * @param x  First set of coordinates
 * @param y  Second set of coordinates
 * @return Equatorial The difference between two Equatorial co-ordinates
 */
static Equatorial operator -(const Equatorial& x, const Equatorial& y) {
  return { x.ra - y.ra, x.dec - y.dec};
}

/**
 * @brief Divides the Equatorial coordinates by a scalar
 * 
 * @param x  First set of coordinates which will receive the result
 * @param y  Second set of coordinates to add
 */
static Equatorial operator /(const Equatorial& x, const double& y) {
  double r = x.ra/y;
  double d = x.dec/y;
  return {r, d};
}

/**
 * @brief Struct for returning the position of the Sun 
 * 
 */
typedef struct SunPosition {
  Equatorial apparent;  /**< The apparent coordinates of the Sun Ra, Dec . */ 
  double eclipticLongitude; /**< Mean Ecliptic Longitude for date. */
  double oblique; /**< The obliquity of the ecliptic for date. */
} SunPosition;

SunPosition calculateSunPosition(const double julianDate, const bool inRadians = false) ;
double calculateTrueObliquity(JulianDate julian, Nutation nut);
double convertToApparentSiderealTime(double siderealTime, Nutation nut, double eps);
double convertToApparentSiderealTime(double siderealTime, JulianDate julian, Nutation nut) ;
Equatorial convertToEquatorial(double alt, double az, JulianDate julian, double latitude, double longitude, double heightAboveSealevel = 0.0, double temperature = 0.0, double pressure = -1.0, bool azWestFromSouth = false, bool applyPrecess = true, bool applyNutation = true, bool applyRefraction = true, bool applyAberration = true);
double convertToJulianDate(const int year,const int month,const int day,const int hour, const int min, const double sec );
JulianDate convertToJulianDate2000(const int year,const int month,const int day,const int hour, const int min, const double sec );
double convertToSidereal(double longitudeEast,const int year,const int month,const int day,const int hour, const int min, const double sec );
double convertToSidereal(double longitudeEast, JulianDate julian);
Nutation nutate(JulianDate julian);
DegMinSecs toDMS(const double& angle);

#endif
