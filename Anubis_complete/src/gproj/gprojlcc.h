
#ifndef GPROJLCC_H
#define GPROJLCC_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: Spherical to Lambert Conformal Conic projection
  Version: $ Rev: $

  2015-02-19 /JD: created

-*/

#include <string>
#include <iomanip>
#include <iostream>

#include "gutils/gpair.h"
#include "gproj/gproj.h"

using namespace std;

namespace gnut {
 
class t_gprojlcc : public t_gproj {
	
 public:
  t_gprojlcc( float truelat1, float truelat2, float moad_cen_lat,
	      float cen_lat,  float cen_lon,  float std_lon );
  virtual ~t_gprojlcc(){};

  virtual t_gpair ll2proj(const t_gpair& latlon,int round = 0 ); // from spherical LL[deg] to LCI XY[idx]  (idx ge 0!)
  virtual t_gpair proj2ll(const t_gpair& xy_idx);                // from LCI XY[idx] to spherical LL[deg]  (idx ge 0!)

  bool operator==(const t_gprojlcc& prj) const;
   
  float truelat1()const{ return _truelat1; }
  float truelat2()const{ return _truelat2; }
  float mcenlat()const{  return _mcenlat;  }
  float cen_lat()const{  return _cen_lat;  }
  float cen_lon()const{  return _cen_lon;  }
  float std_lon()const{  return _std_lon;  }

 protected:

  float   _truelat1;
  float   _truelat2;
  float   _mcenlat;
  float   _cen_lat;
  float   _cen_lon;
  float   _std_lon;
};

} // namespace

#endif
