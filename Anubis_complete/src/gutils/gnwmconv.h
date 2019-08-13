
#ifndef GNWMCONV_H
#define GNWMCONV_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: meteorological data conversion conversion utilities
  Version: $ Rev: $

  2011-04-26 /PV: created
  2012-04-06 /JD: extracted coordinate system conversion utilities only

-*/

#include <string>

using namespace std;

enum ESATUR { HAAN, KRAUS };
   
namespace gnut {
   
double esat(double t, ESATUR etyp=HAAN);  // temperature [K]
double  q2e(double q, double p);          // Pa or hPa (output = input units)
double  e2q(double e, double p);          // Pa or hPa (output = input units)
   
} // namespace
   
#endif // GMETCONV_H
