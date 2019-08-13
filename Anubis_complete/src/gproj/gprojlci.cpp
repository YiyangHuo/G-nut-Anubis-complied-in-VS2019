
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
#include "gproj/gprojlci.h"

using namespace std;

namespace gnut {

// Constructor
// ---------   
t_gprojlci::t_gprojlci( double lon0, double lat0, double mesh,
	                double xoff, double yoff, double radi )
{
  _name = "Lambert Conic Conformal, local indexes [0..long]";
  _id   = PROJ_LCI;
/*
  _lon0 = lon0;
  _lat0 = lat0;
  _xoff = xoff;
  _yoff = yoff;
  _mesh = mesh;
  _radi = radi;
*/
  _xoff =  287;                    // offset [idx]
  _yoff = 1444;                    // offset [idx]
  _lon0 = 17.0;                    // lon0 [deg]
  _lat0 = 46.24470064;             // lat0 [deg]
  _radi = 6371229.0;               // radius [m]
  _mesh = 4710.621094;             // mesh distance [m]

  _K    =             sin(D2R * _lat0);
  _R0   = _radi * pow(cos(D2R * _lat0), 1.0 - _K) * pow(1.0 + _K, _K) / _K;

};


// Conversion from LatLon [deg] sphere to Lambert Conform Conic projection XY [idx] !
// ---------- 
t_gpair t_gprojlci::ll2proj(const t_gpair& latlon, int dec)
{
  gtrace("t_gprojlci::ll2proj(t_gpair)");

  double lat = latlon[0];  
  double lon = latlon[1];
   
  double R   = _R0 * pow((1.0 - sin( D2R*lat))/(1.0 + sin( D2R*lat)), _K/2 );
  double zid = R/_mesh;
  double xid = zid * sin( _K * (lon - _lon0)*D2R );
  double yid = zid * cos( _K * (lon - _lon0)*D2R );

  // decimate request !
  if( dec ){
    xid = round(xid*dec)/dec;
    yid = round(yid*dec)/dec;
  }

  t_gpair xy(_xoff+xid-1, _yoff-yid-1);
   
  return xy;
}; 


// Conversion from Lambert Conform Conic projection XY [idx] LatLon [deg] sphere !
// ---------- 
t_gpair t_gprojlci::proj2ll(const t_gpair& xy_idx)   // xy = 0..X,0..X (starts with ZERO)
{
  gtrace("t_gprojlci::proj2ll(t_gpair)");

  double xidx = 1 + xy_idx[0];
  double yidx = 1 + xy_idx[1];

  double R = _mesh * sqrt( pow(_yoff-yidx, 2.0) + pow(_xoff-xidx, 2.0) );
  double A = pow(_R0/R, 1.0/ _K);

  double lon = _lon0 - R2D*(atan2((_xoff-xidx),(_yoff-yidx)) / _K);
  double lat = R2D*asin( (A*A - 1.0)/( A*A + 1.0) );

  t_gpair ll(lat, lon);
   
  return ll;
}; 


// equal operator
// ----------
bool t_gprojlci::operator==(const t_gprojlci& prj) const
{
  return ( _lon0 == prj.lon0() ) &&
         ( _lat0 == prj.lat0() ) &&
         ( _xoff == prj.xoff() ) &&
         ( _yoff == prj.yoff() ) &&
         ( _mesh == prj.mesh() ) &&
         ( _radi == prj.radi() );
}
   
} // namespace
