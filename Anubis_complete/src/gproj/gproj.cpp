
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

#include "gproj/gproj.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {


// Conversion from LatLon [deg] sphere to generic project projection, keeping the HEIGHT!
// ---------- 
t_gtriple t_gproj::ll2proj(const t_gtriple& latlon, int dec )
{
  gtrace("t_gproj::ll2proj(t_gtriple)");

  t_gpair ll( latlon.gpair() );
  t_gpair xx = ll2proj(ll, dec);
   
  t_gtriple xxh(xx[0], xx[1], latlon[2]);

#ifdef DEBUG
  cout << " XXH: " << xx[0] << " " << xx[1]
       << " LLH: " << ll[0] << " " << ll[1] << endl;
#endif

  return xxh;
}

// Conversion from generic project projection to LatLon [deg] sphere, keeping the HEIGHT!
// ---------- 
t_gtriple t_gproj::proj2ll(const t_gtriple& xy_idx)
{
  gtrace("t_gproj::proj2ll(t_gtriple)");

  t_gpair xy( xy_idx.gpair() );
  t_gpair xx = proj2ll(xy);

  t_gtriple llh(xx[0], xx[1], xy_idx[2]);

#ifdef DEBUG
  cout << " XX: " <<  xx[0] << " " <<  xx[1]
       << " LL: " << llh[0] << " " << llh[1] << " " << llh[2] << endl;
#endif

  return llh;
}

   
} // namespace
