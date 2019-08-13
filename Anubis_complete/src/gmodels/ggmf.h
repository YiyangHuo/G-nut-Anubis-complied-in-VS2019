
#ifndef GGMF_H
#define GGMF_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: Global Mapping Function (GMF)
  Version: $ Rev: $

  Reference:
  J. Boehm, A.E. Niell, P. Tregoning, H. Schuh (2006),
  Global Mapping Functions (GMF): A new empirical mapping function based on numerical 
  weather model data, Geoph. Res. Letters, Vol. 33, L07304, doi:10.1029/2005GL025545.

  A.E. Niell (1996), Global Mapping Functions for the Atmosphere Delay
  of Radio Wavelengths, Journal of Geophysical Research, Vol 101, pp. 3227-3246.
 
  input data
  ----------
  dmjd: modified julian date
  dlat: ellipsoidal latitude in radians
  dlon: longitude in radians
  dhgt: height in m
  zd:   zenith distance in radians

  output data
  -----------
  gmfh,dgmfh: hydrostatic mapping function and derivative wrt z
  gmfw,dgmfw: wet mapping function and derivative wrt z
  
  Johannes Boehm, 2005 August 30
  rev 2006 Aug. 14: derivatives (U. Hugentobler)
  rev 2006 Aug. 14: recursions for Legendre polynomials (O. Montenbruck)

  2012-04-26 /JD: created

-*/

#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>

#include "gutils/gconst.h"

using namespace std;

namespace gnut {   

class t_gmf {

 public:
  t_gmf(){};
 ~t_gmf(){};

  // dlat, dlon --> RADIANS !
  int gmf( double  dmjd, double  dlat, double   dlon, double   dhgt, double zd,
           double& gmfh, double& gmfw, double& dgmfh, double& dgmfw );

 protected:
   static double ah_mean[55];
   static double bh_mean[55];
   static double ah_amp[55]; 
   static double bh_amp[55];
   static double aw_mean[55];
   static double bw_mean[55];
   static double aw_amp[55];
   static double bw_amp[55];

 private:
   
};

} // namespace

#endif
