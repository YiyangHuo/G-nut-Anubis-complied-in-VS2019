
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.
 
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <math.h>
#include <iostream>

#include "gutils/gcommon.h"
#include "gproj/gprojlcc.h"

using namespace std;

namespace gnut {

// Constructor
// ---------   
t_gprojlcc::t_gprojlcc( float truelat1, float truelat2, float moad_cen_lat,
   	                float cen_lat,  float cen_lon,  float std_lon )
{
   _name = "Lambert Conic Conformal, local coordinates [m]";
   _id   = PROJ_LCC;

   _truelat1 = truelat1;
   _truelat2 = truelat2;
   _mcenlat  = moad_cen_lat;
   _cen_lat  = cen_lat;
   _cen_lon  = cen_lon;
   _std_lon  = std_lon;
};


// Conversion from LatLon [deg] sphere to Lambert Conform Conic projection XY [m] !
// ---------- 
t_gpair t_gprojlcc::ll2proj(const t_gpair& latlon, int dec )
{
  gtrace("t_gprojlcc::ll2proj(t_gpair)");

  double xn, psi0, psi1, xc, yc, flp, r, psx, xloc, yloc, ylon;
  float sign = 1.0;
  double conv = 57.29578;
  double pole = 90.0;
   
  double a = 6370.0;
  if( _mcenlat < 0.0 )  sign = -1.0;

  if( fabs(_truelat1) > 90.0 ) //  if( fabs(_truelat1 > 90.0 )) // fixed by JD/2015-08-12
  {
    _truelat1  = 60.0;
    _truelat2  = 30.0;
    _truelat1 *= sign;
    _truelat2 *= sign;
  }
   
  if( _truelat1 == _truelat2 ){
     xn = sin( _truelat2 / conv ); 
  }
  else{
     xn  =   log10( cos( _truelat1 / conv ) )
           - log10( cos( _truelat2 / conv ) );
     xn /=   log10( tan(( 45.0 - sign * _truelat1 / 2.0 ) / conv ))
           - log10( tan(( 45.0 - sign * _truelat2 / 2.0 ) / conv ));
  };
  psi1  = 90.0 - sign*_truelat1;
  psi1 /= conv;

  if( _mcenlat < 0.0 ){ psi1 = - psi1; pole = -pole; }
  
  psi0 = (pole - _mcenlat ) / conv;
  xc = 0.0;
  yc = -a / xn * sin( psi1 ) * pow(  tan(psi0/2.0) / tan(psi1/2.0), xn );
  
  // actual computation for the specified location
  ylon = latlon[1] - _std_lon;
  if( ylon >  180.0 ) ylon -= 360.0;
  if( ylon < -180.0 ) ylon += 360.0;
   
  flp = xn * ylon / conv;
  psx = ( pole - latlon[0] ) / conv;
  r =  -a / xn * sin( psi1 ) * pow( tan(psx/2.0) / tan(psi1/2.0), xn );
   
  if( _mcenlat < 0.0 ){
    xloc =   r* sin(flp);
    yloc =   r* cos(flp);
  }else{
    xloc = - r* sin(flp);
    yloc =   r* cos(flp);
  }
	
  xloc = (xloc - xc);   // km
  yloc = (yloc - yc);   // km
  
  t_gpair xy(xloc, yloc);

  // decimate request !
  if( dec ){
    xy[0] = round(xy[0]*dec)/dec;
    xy[1] = round(xy[1]*dec)/dec;
  }
  
  return xy;
}; 


// from LCI XY[idx] to spherical LL[deg]  (idx ge 0!)
t_gpair t_gprojlcc::proj2ll(const t_gpair& xy_idx)
{
  gtrace("t_gprojlcc::proj2ll(t_gpair)");

  cerr << "t_gprojlcc:proj2ll(t_gpair) - NOT IMPLEMENTED YET\n";
  t_gpair ll(0.0,0.0);

  return ll;
}
   


// equal operator
// ----------
bool t_gprojlcc::operator==(const t_gprojlcc& prj) const
{
  return ( _truelat1 == prj.truelat1() ) &&
         ( _truelat2 == prj.truelat2() ) &&
         ( _mcenlat  == prj.mcenlat()  ) &&
         ( _cen_lat  == prj.cen_lat()  ) &&
         ( _cen_lon  == prj.cen_lon()  ) &&
         ( _std_lon  == prj.std_lon()  );
}

   
} // namespace
