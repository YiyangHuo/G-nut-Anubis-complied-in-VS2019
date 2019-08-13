
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

#include <cmath>
#include <iostream>
#include <iomanip>

#include "gutils/gconst.h"
#include "gutils/gnwmconv.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {


// saturated WV pressure
// ----------
double esat(double t, ESATUR etyp)  // temperature [K]
{
  if( t < 140 ) t = t + TC2TK;      // deg C -> K (check)

  switch( etyp ){
    case HAAN  : return (6.1070 * exp(17.38*(t-TPOINT)/(t-34.16)) );        // hPa (info: de Haan, 2013, KNMI)
    case KRAUS : return (6.1078 * exp(17.1 *(t-TC2TK)/(235 + t-TC2TK)) );   // hPa (book: J Boehm, 2013, TUW)
  }

  return 0.0;
}

// SHUM to EPRE
// ----------
double q2e(double q, double p)   // Pa or hPa (output = input units)
{
   return ( q * p/(Eps + q*(1 - Eps)));
}
   

// EPRE to SHUM
// ----------
double e2q(double e, double p)   // Pa or hPa (output = input units)
{
  return ( Eps * e/(p - (1-Eps)*e));
}


} // namespace
