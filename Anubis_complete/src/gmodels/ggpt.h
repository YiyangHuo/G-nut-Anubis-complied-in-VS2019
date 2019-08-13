
#ifndef GGPT_H
#define GGPT_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: Global Pressure Temperature model 
  Version: $ Rev: $

  Reference:
  J. Boehm, R. Heinkelmann, H. Schuh, Short Note: A Global Model of Pressure 
  and Temperature for Geodetic Applications, Journal of Geodesy, 
  doi:10.1007/s00190-007-0135-3, 2007.
 
  input data
  ----------
  dmjd: modified julian date
  dlat: ellipsoidal latitude in radians
  dlon: longitude in radians
  dhgt: ellipsoidal height in m

  output data
  -----------
  pres: pressure in hPa
  temp: temperature in Celsius
  undu: Geoid undulation in m (from a 9x9 EGM based model)

  Johannes Boehm, 2006 June 12
  rev 2006 June 16: geoid undulation is accounted for
  rev 2006 July 17: barometric height equation used instead of Berg (1948) for
                    the vertical extrapolation of the pressure
  ref 2006 Aug. 14: recursions for Legendre polynomials (O. Montenbruck)
  ref 2011 Jul. 21: latitude -> ellipsoidal latitude (J. Boehm)

  Adaptation to Bernese GPS Software:  U. Hugentobler (14-08-2006)
  Recursions for Legendre polynomials: O. Montenbruck (14-08-2006)
  Adaptation to G-Nut Software:        J. Dousa       (28-11-2011)

  2011-11-28 /JD: created

-*/

#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;

namespace gnut {   

class t_gpt {

 public:
  t_gpt(){};
 ~t_gpt(){};
 
  // dlat, dlon --> RADIANS !
  int gpt_v1( double  dmjd, double  dlat, double  dlon, double  dhgt,
  	      double& pres, double& temp, double& undu );

 protected:
   static double a_geoid[55];
   static double b_geoid[55];
   static double ap_mean[55]; 
   static double bp_mean[55];
   static double ap_amp[55];
   static double bp_amp[55];
   static double at_mean[55];
   static double bt_mean[55];
   static double at_amp[55];
   static double bt_amp[55];

 private:
   
};

} // namespace

#endif
