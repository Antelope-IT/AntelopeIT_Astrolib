/**
 * @file AntelopeIT_Astrolib.cpp
 * @mainpage Antelope IT Astronomy Library
 *
 * @section intro_sec Introduction
 *
 * This is the documentation for Antelope IT's Astronomy Library for the
 * Arduino platform.  It is inspired by the IDLAstro Library 
 * ---> https://idlastro.gsfc.nasa.gov/contents.html
 * 
 *  and the book Astronomical Algorithms (second edition) by Jean Meeus.
 *    
 *  This library is not a complete implementation of the entire IDLAstro library or
 *  all of the algorithms in the book, it is just a selection of methods and
 *  functions to meet the needs of a particular project.
 *
 * @section dependencies Dependencies
 *
 * This library has no dependencies other than math.h.
 *
 * @section author Author
 *
 * Antelope IT (contact@antelope-it.co.uk)
 * 
 * Copyright Antelope IT Ltd (c) 2021
 *
 * @section license License
 *
 * GPLv3 license, all text here must be included in any redistribution.
 *
 */

#include "AntelopeIT_Astrolib.h"

/**
 * @brief Converts an angle expresed in decimal degrees to a structure with Degrees, Minutes and Seconds
 * 
 * @param angle angle expressed in decimal degrees nnn.nnnnn... to convert
 * @return DegMinSecs 
 */
DegMinSecs toDMS(const double& angle) {
  int d = static_cast<int>(angle);
  double mr = (angle-d)*60;
  int m = static_cast<int>(mr);
  double s = (mr - m) * 60;

  return {d,abs(m),abs(s)};
}


/**
 * @brief   Division, rounding down (rather than towards zero).
 * From C++11 onwards, integer division is defined to round towards zero, so we
 * can rely on that when implementing this.  This is only used with denominator b
 * > 0, so we only have to treat negative numerator, a, specially.
 * 
 * @param a 
 * @param b 
 * @return int 
 */
static inline int intdiv(int a, int b)
{
    return (a - (a < 0 ? b - 1 : 0)) / b;
}

/**
 * @brief Calculate a polynomial of the form
 *        y = A + Bx + Cx^2 + Dx^3... Nx^(n-1) using powers.
 * 
 * @param x The value for `x` to use in the polynomial
 * @param coefficients An array of coefficents [A,B,C,D...N]
 * @param len The number of coefficients (n)
 * @return double The result
 */
double calculate_polynomial(double x, const double coefficients[], int len) {
  double sum = 0.0;
  for (int i = 0; i < len; i++)
  {
    sum+=coefficients[i] * pow(x, i);
  }
  return sum;
}

/**
 * @brief Uses Horner's method to calculate a polynomial of the form
 *        y = A + Bx + Cx^2 + Dx^3... Nx^(n-1) 
 *        From: Meeus Astronomical Algorithms Chap. 1 pg. 10-11.
 * 
 * @param x The value for `x` to use in the polynomial
 * @param coefficients An array of coefficents [A,B,C,D...N]
 * @param len The number of coefficients (n)
 * @return double The result
 */
double calculate_hornerPolynomial(double x, const double coefficients[], int len) {
  double result = 0.0;
  for(int i = len; i > 0; i--) 
  {
    result = x * result + coefficients[i-1];
  }

  return result;
}

/**
 * @brief Returns an equivalent angle within the range of 0 - 360 degrees or 2 pi radians 
 * 
 * @param angle double angle to transpose to the range of a full circle
 * @param isRadians bool Is the angle in degrees or radian (default = degrees)
 * @return double Equivalent angle in range 0 - 360 degrees 2 pi radians  
 */
double circle_range(double angle, bool isRadians = false) {
  double limit = isRadians ? M_PI * 2 : 360.0;
  angle = fmod(angle,limit);
  if(angle < 0) {
    angle+=limit;
  }
  return angle;
};

/**
 * @brief Converts a Gregorian calender date to a Julian date
 * 
 * @param year Calendar Year 
 * @param month Calendar Month
 * @param day Calendar Day
 * @param hour Hour 
 * @param min Minutes
 * @param sec Seconds
 * @return double Julian Date 
 */
double convertToJulianDate(const int year,const int month,const int day,const int hour, const int min, const double sec ) {
  double julianTime = (hour / 24.0) + (min / 1440.0) + (sec / 86400.0) - 0.5;
  int leapYearAdjust = intdiv(14 - month, 12);
  unsigned int y = (unsigned int)year + 4800 - leapYearAdjust;
  int m = month + 12 * leapYearAdjust - 3;
  double julianDate = static_cast<double>(day + intdiv(153 * m + 2, 5) + 365 * y + intdiv(y, 4) - intdiv(y, 100) + intdiv(y, 400) - 32045);
  julianDate+=julianTime;

  return julianDate;
}

/**
 * @brief Converts a Gregorian calender date to a Julian date
 * 
 * @param year Calendar Year 
 * @param month Calendar Month
 * @param day Calendar Day
 * @param hour Hour 
 * @param min Minutes
 * @param sec Seconds
 * @return JulianDate Structure containing date, number of centuries and number of days since J2000
 */
JulianDate convertToJulianDate2000(const int year,const int month,const int day,const int hour, const int min, const double sec ) {
  double julianTime = (hour / 24.0) + (min / 1440.0) + (sec / 86400.0) - 0.5;
  
  int leapYearAdjust = 0;
  if(month ==1 || month ==2 ) {
    leapYearAdjust = -1;
  }
  double julianDate = static_cast<double>(day - 32075l + intdiv(1461l* (year+4800l+leapYearAdjust),4) + intdiv(367l*(month - 2-leapYearAdjust*12),12) - intdiv(3*intdiv(year+4900l+leapYearAdjust,100),4));
  julianDate+=julianTime;

  // Compute the number of centuries since J2000
  double julianDays = julianDate - JD2000;
  double julianCenturies = julianDays/36525.0;  
 
  return {julianDate, julianCenturies, julianDays};
}

/**
 * @brief Converts a Local Julian Date to Local mean Sidereal time 
 * 
 * @param longitudeEast The longitude of the location in degrees measured east of Greenwich
 * @param julian The JulianDate to conver to Sidereal time
 * @return double 
 */
double convertToSidereal(double longitudeEast, JulianDate julian) {
  static const double constants[] = {280.46061837e0, 360.98564736629e0, 0.000387933e0, 38710000.0};
  
  double t0 = julian.julianDays;
  double t = julian.julianCenturies;

  double theta = constants[0]+ (constants[1]*t0 ) + pow(t,2)*(constants[2] -t/constants[3]);
  double lst = (theta + longitudeEast) / hoursToDegrees;
  if(lst < 0) {
    lst = 24.0 + (fmod(lst,24.0));
  }
  lst = fmod(lst, 24.0);
  return lst;
}

/**
 * @brief Converts a Local civillian time UT to a Local mean Sidereal time 
 * 
 * @param longitudeEast  The longitude of the location in degrees measured east of Greenwich
 * @param year Calendar Year 
 * @param month Calendar Month
 * @param day Calendar Day
 * @param hour Hour 
 * @param min Minutes
 * @param sec Seconds
 * @return double The local mean sidereal time for the given date / time
 */
double convertToSidereal(double longitudeEast,const int year,const int month,const int day,const int hour, const int min, const double sec ) {
  JulianDate jd = convertToJulianDate2000(year, month,day,hour, min, sec );
  return convertToSidereal(longitudeEast, jd);
}

/**
 * @brief Calculates the true obliquity of the ecliptic given the Julian date and the nutation in longitude and obliquity on that date
 * 
 * @param julian JulianDate The Julian Date to calculate the true obliquity for
 * @param nut Nutation A nutation result calculated for the Julian Date in julian
 * @return double The True obliquity of the ecliptic in radians (eps - epsilon in Meeus)
 */
double calculateTrueObliquity(JulianDate julian, Nutation nut) {           
  static const double coeffs[] = {23.439291111*3600.0, -46.8150, -0.00059, 0.001813};
  double eps0 = calculate_hornerPolynomial(julian.julianCenturies, coeffs,4 );   // mean obliquity of the ecliptic in degrees
  return (eps0 + nut.oblique)/3600.0*degreesToRadians;                            // true obliquity of the ecliptic in radians
}

/**
 * @brief Convert a Local Mean Sidereal Time (LMST) to the equivalent Local Apparent Sidereal Time (LAST)
 * 
 * Where the Nutation and true obliquity of the ecpliptic are already known then this method can be used 
 * as an optimisation - saves recalculating the true obliquity value
 * 
 * @param siderealTime The Local Mean Sidereal Time (LMST) to convert
 * @param nut Nutation result for the given LMST
 * @param eps The true obliquity of ecliptic in radians calculated using the Nutation result calculated for the given LMST 
 * @return double The Local Apparent Aidereal Time for the given LMST
 */
double convertToApparentSiderealTime(double siderealTime, Nutation nut, double eps) {
  return (3600.0 * siderealTime * hoursToDegrees + nut.longitude * cos(eps)) / (hoursToDegrees * 3600.0);
}

/**
 * @brief Convert a Local Mean Sidereal Time (LMST) to the equivalent Local Apparent Sidereal Time (LAST)
 * 
 * @param siderealTime The Local Mean Sidereal Time (LMST) to convert
 * @param julian JulianDate Julian date used to calculate the local mean sidereal time
 * @param nut Nutation result for the given LMST
 * @return double The Local Apparent Aidereal Time for the given LMST
 */
double convertToApparentSiderealTime(double siderealTime, JulianDate julian, Nutation nut) {
  double eps = calculateTrueObliquity(julian, nut);
  return convertToApparentSiderealTime(siderealTime, nut, eps);
}

/**
 * @brief Calculates a correction factor for atmospheric refaction given either temperature, pressure or altitude
 * 
 * @param altitude Observed Elevation measured in degrees
 * @param temperature Temperature in degrees K (optional default = 283.0)
 * @param pressure Pressure milibars (optional default = 1010)
 * @return double Correction in degrees to apply to the Observed elevation
 */
double calculateAtmosphericRefractionCorrection(double altitude, double temperature = 283.0, double pressure = 1010.0) {
  double refraction = altitude < 15 ? 3.569*(0.1594 + 0.0196*altitude + 0.00002*pow(altitude,2))/(1.0 + 0.505*altitude+0.0845*pow(altitude,2)) : 0.0166667/tan((altitude + 7.31/(altitude+4.4))*degreesToRadians);
  double tempcorrection = pressure/1010.0 * 283.0/temperature;
  refraction = tempcorrection * refraction;
  return refraction;
}

/**
 * @brief Given the observed elevation calculates the apparent 
 * 
 * @param altitude Observed Elevation measured in degrees
 * @param heightAboveSeaLevel Altitude above sea level in meters only used if Temperature or Pressure are not specified
 * @param temperature Temperature in degrees K (optional default = 283.0)
 * @param pressure Pressure milibars (optional default = 1010)
 * @return double The apparent altitude/elevation measured in degrees.
 */
double correctForAtmosphericRefraction(double altitude, double heightAboveSeaLevel, double temperature = 0.0, double pressure = - 1.0) {
  static const double alpha = 0.0065;  // temp lapse rate [deg c per meter]
  pressure = pressure > 0 ? pressure : 1010 *pow((1-6.5/288000*heightAboveSeaLevel),5.255);
  temperature = temperature > 0 ? temperature : (heightAboveSeaLevel>11000) ? 211.5 : 283.0 - alpha*altitude;
  return altitude - calculateAtmosphericRefractionCorrection(altitude, temperature,  pressure);
}

/**
 * @brief Returns the nutation in longitude and obliquity for a given Julian date.
 * 
 * @param julian Julian Date to do the nutation calculation for.
 * @return Nutation The calculated nutation in longitude and obliquity
 */
Nutation nutate(JulianDate julian) {
  // Julian Centuries from 2000
  const double jC = julian.julianCenturies;

  // Calculate Mean Elongation of the Moon
  static const double coeff1[] = {297.85036,  445267.111480, -0.0019142, 1.0/189474 };
  double d = calculate_hornerPolynomial(jC, coeff1, 4);
  d = circle_range(d) * M_PI / 180;

  // Calculate the Sun's Mean Anomaly
  static const double coeff2[] = {357.52772, 35999.050340, -0.0001603, -1.0/3e5 };
  double m = calculate_hornerPolynomial(jC, coeff2, 4);
  m = circle_range(m) * M_PI / 180;

  // Calculate the Moons's Mean Anomaly
  static const double coeff3[] = {134.96298, 477198.867398, 0.0086972, 1.0/5.625e4 };
  double mprime = calculate_hornerPolynomial(jC, coeff3, 4);
  mprime = circle_range(mprime) *  M_PI / 180;

  // Calculate the Moon's argument of Latitude
  static const double coeff4[] = {93.27191, 483202.017538, -0.0036825, -1.0/3.27270e5 };
  double F = calculate_hornerPolynomial(jC, coeff4, 4);
  F = circle_range(F) *  M_PI / 180;

  // Calculate Longitude of the ascending node of the Moon's mean orbit on the ecliptic, measured from the mean equinox of the date
  static const double coeff5[] = {125.04452, -1934.136261, 0.0020708, 1./4.5e5 };
  double omega = calculate_hornerPolynomial(jC, coeff5, 4);
  omega = circle_range(omega) * M_PI / 180;

  static const int d_lng[] = {0,-2,0,0,0,0,-2,0,0,-2,-2,-2,0,2,0,2,0,0,-2,0,2,0,0,-2,0,-2,0,0,2,
    -2,0,-2,0,0,2,2,0,-2,0,2,2,-2,-2,2,2,0,-2,-2,0,-2,-2,0,-1,-2,1,0,0,-1,0,0, 
    2,0,2};

  static const int m_lng[] ={ 0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,1,0,-1,0,0,0,1,1,-1,0, 
    0,0,0,0,0,-1,-1,0,0,0,1,0,0,1,0,0,0,-1,1,-1,-1,0,-1};

  static const int mp_lng[] = {0,0,0,0,0,1,0,0,1,0,1,0,-1,0,1,-1,-1,1,2,-2,0,2,2,1,0,0,-1,0,-1, 
    0,0,1,0,2,-1,1,0,1,0,0,1,2,1,-2,0,1,0,0,2,2,0,1,1,0,0,1,-2,1,1,1,-1,3,0};
  
  static const int f_lng[] = {0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,2,2,2,2,0,0,2,0,0, 
    0,-2,2,2,2,0,2,2,0,2,2,0,0,0,2,0,2,0,2,-2,0,0,0,2,2,0,0,2,2,2,2};

  static const int om_lng[] = {1,2,2,2,0,0,2,1,2,2,0,1,2,0,1,2,1,1,0,1,2,2,0,2,0,0,1,0,1,2,1, 
    1,1,0,1,2,2,0,2,1,0,2,1,1,1,0,1,1,1,1,1,0,0,0,0,0,2,0,0,2,2,2,2};

  static const int sin_lng[] = {-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217, 
    -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 26, -22, 
    21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, -7, 
    6,6,6,-6,-6,5,-5,-5,-5,4,4,4,-4,-4,-4,3,-3,-3,-3,-3,-3,-3,-3 };

  static const double sdelt[] = {-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0.0, -0.5, 0.0, 0.1, 
    0.0, 0.0, 0.1, 0.0,-0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.1, 0.0, 0.1, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const int cos_lng[] = { 92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95,0,-70,-53,0, 
    -33, 26, 32, 27, 0, -24, 16,13,0,-12,0,0,-10,0,-8,7,9,7,6,0,5,3,-3,0,3,3,
    0,-3,-3,3,3,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

  static const double cdelt[] = {8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0  };

  Nutation nut = {};

  for(size_t i = 0; i < 63; i++) {
    double arg = d*d_lng[i] + m*m_lng[i] +mprime*mp_lng[i] + F*f_lng[i] +omega*om_lng[i];
    double sarg = sin(arg);
    double carg = cos(arg);
    nut.longitude+=(sdelt[i]*jC + sin_lng[i])*sarg;
    nut.oblique+=(cdelt[i]*jC + cos_lng[i])*carg;
  }

  nut.longitude*=0.0001;
  nut.oblique*=0.0001;

  return nut;
}

/**
 * @brief Calculates the correction for a set of observed coordinates due to nutation
 * 
 * Where the True obliquity of the ecliptic is already known can be used as a optimisation
 * saving re-calculating the the true obliquity value
 * 
 * @param nut Nutation in longitude and obliquity
 * @param eps True obliquity of teh ecliptic in degrees
 * @param coords The coordinates to calculate the correction for.
 * @return Equatorial 
 */
Equatorial correctForNutation(Nutation nut, double eps, Equatorial coords) {
  double ce = cos(eps);
  double se = sin(eps);
  double x = cos(coords.ra*degreesToRadians) * cos(coords.dec*degreesToRadians);
  double y = sin(coords.ra*degreesToRadians) * cos(coords.dec*degreesToRadians);
  double z = sin(coords.dec*degreesToRadians);
  double x2 = x - (y*ce + z*se)*nut.longitude * secondsToRadians;
  double y2 = y + (x*ce*nut.longitude - z*nut.oblique) * secondsToRadians;
  double z2 = z + (x*se*nut.longitude + y*nut.oblique) * secondsToRadians;
  double r = sqrt(pow(x2,2) + pow(y2,2) + pow(z2,2));
  double xyproj = sqrt(pow(x2,2) + pow(y2,2));
  double ra2 = x2 * 0.0;
  double dec2= x2 * 0.0;

  if ((xyproj==0) && (z!=0)) {
    // places where xyproj=0 (point at NCP or SCP)
    dec2 = asin(z2/r);
    ra2 = 0.0;
  }                            
  
  if(xyproj!=0) {
    // places other than NCP or SCP
    ra2= atan2(y2,x2);
    dec2 = asin(z2/r);
  }
  
  // convert to DEGREES
  ra2 = ra2 /degreesToRadians;
  dec2 = dec2 /degreesToRadians;

  if(ra2<0.0) {
    ra2 = ra2+ 360.0;
  }

  double d_ra = (ra2 - coords.ra) * 3600.0;
  double d_dec = (dec2 - coords.dec) * 3600.0;

  return Equatorial {d_ra, d_dec};
}

/**
 * @brief Calculates the correction for a set of observed coordinates due to nutation
 * 
 * @param jd JulianDate Julian date to calculate the correction for.
 * @param coords The coordinates to calculate the correction for.
 * @return Equatorial 
 */
Equatorial correctForNutation(JulianDate jd, Equatorial coords) { 
  Nutation nut = nutate(jd);
  // Calulate the true Obliquity
  double eps = calculateTrueObliquity(jd, nut);
  return correctForNutation(nut, eps, coords);
}

/**
 * @brief Calculate the RA and Dec of the Sun at a given Julian date.
 * 
 * @param julianDate Julian Date of the day and time
 * @param inRadians Indicates the coordinates are given in radians
 * @return SunPosition The apparent position of the Sun for the given date with obliquity and longitude of the ecliptic 
 */
SunPosition calculateSunPosition(const double julianDate, const bool inRadians) {
  // Julian Centuries from 1900
  const double jC = (julianDate - 2415020.0)/36525.0;

  // Suns mean longitude
  double L = (279.696678e0+(fmod((36000.768925e0*jC), 360.0e0)))*3600.0e0;
  // Earths mean anomaly
  double me = 358.475844e0 + (fmod((35999.049750e0*jC), 360.0e0));
  // Correction for Earths Eliptic orbit
  double meRadians = me*degreesToRadians;
  double ellcor = (6910.1e0 - 17.2e0*jC)*sin(meRadians) + 72.3e0*sin(2.0e0*meRadians);
  L += ellcor;

  // Venus' mean anomaly
  double mv = 212.603219e0 + (fmod((58517.803875e0*jC) , 360.0e0)); 
  // Correction for Venus' Perturbations
  double vencorr =
    4.8e0 * cos((299.1017e0 + mv - me)*degreesToRadians) + 
      5.5*cos((148.3133e0 + 2.0e0 * mv - 2.0e0 * me )*degreesToRadians) + 
      2.5*cos((315.9433e0 + 2.0e0 * mv - 3.0e0 * me )*degreesToRadians) + 
      1.6*cos((345.2533e0 + 3.0e0 * mv - 4.0e0 * me )*degreesToRadians) + 
      1.0*cos((318.15e0   + 3.0e0 * mv - 5.0e0 * me )*degreesToRadians);
  
  L += vencorr;

  // Mars' mean anomaly
  double mm = 319.529425e0  +  (fmod(( 19139.858500e0 * jC) , 360.0));
  // Correction for Mars' perturbations
  double marscorr =
    2*cos((343.8883e0 -2*mm  +  2*me)*degreesToRadians) + 
    1.8*cos((200.4017e0 -  2*mm  + me)*degreesToRadians);
  
  L += marscorr;

  // Jupiter's mean anomaly
  double mj = 225.328328e0  +  (fmod(( 3034.6920239e0 * jC)  , 360.0e0));
  // Correction for Jupiter's perturbations
  double jupcorr = 7.2e0 * cos(( 179.5317e0 - mj + me )* degreesToRadians) + 
    2.6e0 * cos((263.2167e0  -  mj ) *degreesToRadians) + 
      2.7e0 * cos(( 87.1450e0  -  2.0e0 * mj  +  2.0e0 * me ) * degreesToRadians) + 
        1.6e0 * cos((109.4933e0  -  2.0e0 * mj  +  me ) * degreesToRadians);
  L += jupcorr;

  // Correct for the Moons perturbations using the mean elongation of the Moon from the Sun (d)
  double d = 350.7376814e0  + (fmod(( 445267.11422e0 * jC) , 360.0e0 ));
  double mooncorr  = 6.5e0 * sin(d*degreesToRadians);
  L += mooncorr;

  // Correct for long period terms
  double longterm  = + 6.4e0 * sin((231.19e0 + 20.20e0 * jC )*degreesToRadians);
  L  +=  longterm;
  L = fmod(( L + 2592000.0e0) , 1296000);
  double longmed = L/3600.0e0;

  // Correct for Aberration
  L  -=  20.5e0;

  // Correct for Nutation using the longitude of the Moons mean node OMEGA
  double omega = 259.183275e0 - (fmod(( 1934.142008e0 * jC ) , 360.0e0));
  double omegaRadians = omega*degreesToRadians;
  L -=  17.2e0 * sin(omegaRadians);

  // Calculate the true obliquity
  double oblt  = 23.452294e0 - 0.0130125e0*jC + (9.2e0*cos(omegaRadians))/3600.0;
  double obltRadians = oblt*degreesToRadians;
  
  // Calculate the RA and DEC
  L = L/3600;
  double lradians = L*degreesToRadians;
  double ra  = atan2(sin(lradians) * cos(obltRadians) , cos(lradians) );

  if(ra<0.0e0) {
    ra += 2.0* M_PI;
  } 

  double dec = asin(sin(lradians) * sin(obltRadians));
  
  if(inRadians){
    oblt = obltRadians;
    longmed = longmed*degreesToRadians;
  } else {
    ra = ra/degreesToRadians;
    dec = dec/degreesToRadians;
  }

  return {{ra, dec},longmed,oblt};
}

/**
 * @brief Calculate the correction due to Aberration for a set of coordinates
 * 
 * @param julian The Julian date to calculate the correction for.
 * @param coords The coordinates RA, Dec in degrees to calculate the correction for 
 * @param oblique The true obliquity of the ecliptic
 * @return Equatorial The correction to apply to the observed coordinates. 
 */
Equatorial correctForAberration(JulianDate julian, Equatorial coords, double oblique) {
  const double k = 20.49552;                             // Constant of aberration, in arcseconds
  double jC = julian.julianCenturies;                    // Julian centuries from J2000 of jd.
  SunPosition sunPos = calculateSunPosition(julian.julianDate, true);

  double jCSqr = pow(jC,2);
  double e = 0.016708634 - 0.000042037*jC - 0.0000001267*jCSqr ;
  double pi = 102.93735 + 1.71946*jC + 0.00046*jCSqr;

  double decRadians = coords.dec*degreesToRadians;
  double raRadians = coords.ra*degreesToRadians;
  double piRadians = pi*degreesToRadians;
  double cd = cos(decRadians); 
  double sd = sin(decRadians);
  double ce = cos(oblique); 
  double te = tan(oblique);
  double cp = cos(piRadians);
  double sp = sin(piRadians);
  double cs = cos(sunPos.eclipticLongitude);
  double ss = sin(sunPos.eclipticLongitude);
  double ca = cos(raRadians); 
  double sa = sin(raRadians);
  double term1 = (ca*cs*ce+sa*ss)/cd;
  double term2 = (ca*cp*ce+sa*sp)/cd;
  double term3 = (cs*ce*(te*cd-sa*sd)+ca*sd*ss);
  double term4 = (cp*ce*(te*cd-sa*sd)+ca*sd*sp);
  double d_ra = -k * term1 + e*k * term2;
  double d_dec = -k * term3 + e*k * term4;

  return {d_ra,d_dec};
}

/**
 * @brief Given a 3x3 matrix it populates it with the precession coefficients which can be 
 * used to precess a set of coordinates from one equinox to another.
 * 
 * 
 * @param equinoxStart Julian Equinox to start the precession calculation from  
 * @param euqinoxEnd  Julian Equinox to end the precession calculation on
 * @param matrix The populated matrix. 
 */
 void populatePrecessionMatrix( double equinoxStart, double euqinoxEnd, double matrix[3][3]) {
  double t = 0.001e0*( euqinoxEnd - equinoxStart);
  double st = 0.001e0*( equinoxStart - 2000.0);
  double a = secondsToRadians * t * (23062.181e0 + st*(139.656e0 +0.0139e0*st) + t*(30.188e0 - 0.344e0*st+17.998e0*t));
  double b = secondsToRadians * t * t * (79.280e0 + 0.410e0*st + 0.205e0*t) + a;
  double c = secondsToRadians * t * (20043.109e0 - st*(85.33e0 + 0.217e0*st) + t*(-42.665e0 - 0.217e0*st -41.833e0*t));

  double sina = sin(a);
  double sinb = sin(b);
  double sinc = sin(c);
  double cosa = cos(a);
  double cosb = cos(b);
  double cosc = cos(c);

  double premat[3][3] = {
    {cosa*cosb*cosc-sina*sinb, sina*cosb+cosa*sinb*cosc, cosa*sinc},
    {-cosa*sinb-sina*cosb*cosc, cosa*cosb-sina*sinb*cosc, -sina*sinc},
    {-cosb*sinc, -sinb*sinc, cosc}
  };

  for (size_t r = 0; r < 3; r++)
  {
    for (size_t c = 0; c < 3; c++)
    {
      matrix[r][c] = premat[r][c];
    }
  }
};

/**
 * @brief Precess a set of coordinates from one equinox to another.
 * 
 * @param coords The coordinates RA, Dec in dregrees to precess
 * @param equinoxStart The Julian Equinox for the starting coordinates
 * @param equinoxEnd  The Julian Equinox to precess the coordinates to. 
 * @return Equatorial The coordinates precessed to the destination equinox.
 */
Equatorial precess(Equatorial coords,double equinoxStart, double equinoxEnd) {
  double ra_rad = coords.ra*degreesToRadians;   
  double dec_rad = coords.dec*degreesToRadians;
  double a = cos(dec_rad); 
  double x[3] = {a*cos(ra_rad), a*sin(ra_rad), sin(dec_rad)}; 

  double r [3][3] = {
    {0.0,0.0,0.0},
    {0.0,0.0,0.0},
    {0.0,0.0,0.0},
  };
  
  populatePrecessionMatrix(equinoxStart,equinoxEnd, r);

  double x2[3] = {0.0,0.0,0.0};

  for (size_t i = 0; i < 3; i++)
  {
    x2[i] = r[0][i] * x[0]; 
    x2[i]+= r[1][i] * x[1];
    x2[i]+= r[2][i] * x[2]; 
  }
    
  ra_rad = atan2(x2[1],x2[0]);
  dec_rad = asin(x2[2]);

  double ra = ra_rad / degreesToRadians;

  if(ra < 0.0) {
    ra+= 360.0;
  }

  double dec = dec_rad / degreesToRadians;

  return {ra , dec};
}

/**
 * @brief Given a set of altitude, azimuth coordinates calculate the equivalent equatorial coordinates expressed as Hour angle, Declination   
 * 
 * @param altitudeRadians Altitude / Elevation in radians
 * @param azimuthRadians  Azimuth in radians
 * @param latitudeRadians Latitude of the observer expressed in radians
 * @return HourAngleDeclination Coordinates expressed as Ha in degrees, Declination in degrees for the given coordinates.
 */
HourAngleDeclination convertToHourAngleDeclination(double altitudeRadians, double azimuthRadians, double latitudeRadians) {
  double sinAlt = sin(altitudeRadians);
  double cosAlt = cos(altitudeRadians);
  double sinAz = sin(azimuthRadians);
  double cosAz = cos(azimuthRadians);
  double sinLat = sin(latitudeRadians);
  double cosLat = cos(latitudeRadians);

  double hourAngle = atan2(-sinAz * cosAlt, - cosAz * sinLat * cosAlt + sinAlt * cosLat);
  double hourAngleDegrees = hourAngle / degreesToRadians;
  
  if(hourAngleDegrees < 0.0){
    hourAngleDegrees += 360.0;
    hourAngleDegrees = fmod(hourAngleDegrees , 360.0);
  }
    
  double sinDeclination = sinLat * sinAlt + cosLat * cosAlt * cosAz;
  double declination = asin(sinDeclination) / degreesToRadians;

  return {hourAngleDegrees, declination};
}

/**
 * @brief Convert a set of coordinates Hour Angle, Declination in degrees to RA, Dec
 * 
 * @param coords Coordinates to convert
 * @param localApparentSiderealTimeDegrees The Local Apparent Sidereal Time (LAST)  
 * @return Equatorial 
 */
Equatorial convertHaDecToRaDec(HourAngleDeclination coords, double localApparentSiderealTimeDegrees) {
  double hourAngleDegrees = coords.ha;
  double dec = coords.dec;

  double ra = fmod((localApparentSiderealTimeDegrees - hourAngleDegrees + 360.0) , 360.0);
  return {ra, dec};
}

/**
 * @brief Given an observed local set of coordinates Altitude and Azimuth calculate the equivalent equatorial
 * coordinates RA, Dec in degrees.
 * 
 * @param alt Altitude / Elevation in degrees
 * @param az  Azimuth angle in degrees measured East from North (default azWestFromSouth = false) 
 * @param julian Julian Date of the observation 
 * @param latitude North geodetic latitude of the observation location in degrees
 * @param longitude East longitude of location in degrees (Specify West longitude with a negative sign.)
 * @param heightAboveSealevel Height above sea level in meters of observation location
 * @param temperature Temperature degrees Kelvin at the Observation location 
 * @param pressure  Pressure milibars at the Observation location
 * @param azWestFromSouth Set to true if the Azimuth value is measured West from South  (default = false)
 * @param applyPrecess Correct the observed coordinates for Precession (default = true)
 * @param applyNutation Correct the observed coordinates for Nutation (default = true)
 * @param applyRefraction Correct the observed coordinates for Atmospheric refraction Re(default = true)
 * @param applyAberration Correct the observed coordinates for Aberration (default = true)
 * @return Equatorial The set of coordinates RA, Dec expressed in degrees.
 */
Equatorial convertToEquatorial(double alt, double az, JulianDate julian, double latitude, double longitude, double heightAboveSealevel, double temperature, double pressure, bool azWestFromSouth, bool applyPrecess, bool applyNutation, bool applyRefraction, bool applyAberration)  {
 
  if(applyRefraction) {
    alt = correctForAtmosphericRefraction(alt, heightAboveSealevel, temperature, pressure);
  }

  if(azWestFromSouth) {
    az-= 180.0; 
  }

  Nutation nut = nutate(julian);
  double localMeanSiderealTime = convertToSidereal(longitude,julian);
  double epsilon = calculateTrueObliquity(julian, nut);
  double localApparentSiderealTime = convertToApparentSiderealTime(localMeanSiderealTime, nut, epsilon);
  HourAngleDeclination coordsHaDec = convertToHourAngleDeclination(alt * degreesToRadians, az * degreesToRadians ,latitude*degreesToRadians);

  Equatorial coords = convertHaDecToRaDec(coordsHaDec, localApparentSiderealTime * hoursToDegrees);

  if(applyNutation || applyAberration) {
    Equatorial nutationAdjust = applyNutation ? correctForNutation(nut,epsilon,coords) : Equatorial {};
    Equatorial aberrationAdjust = applyAberration ? correctForAberration(julian,coords,epsilon): Equatorial {};

    Equatorial totalAdjust = aberrationAdjust + nutationAdjust;
    coords= coords - totalAdjust/3600;
  }

  if(applyPrecess) {
    double currentEquinox = (julian.julianDate - JD2000)/365.25 + 2000.0; // compute current equinox
    coords = precess(coords, currentEquinox, 2000.0);
  }

  return coords;
}
