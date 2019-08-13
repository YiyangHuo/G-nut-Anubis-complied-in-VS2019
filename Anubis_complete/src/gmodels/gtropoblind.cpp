
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
#include <iomanip>

#include "gmodels/gtropoblind.h"

using namespace std;

namespace gnut {   

// ----------------------------- BLIND MODELS ------------------------------------

// ---------
// MOPS ZHD interface
// ---------
double t_blindmops::getZHD(const t_gtriple& ell, const t_gtime& epo) // RAD
{
  gtrace("t_blindmops::getZHD");
      
  double zhd, zwd;
  double ele = 90.0;
  double doy = double( epo.doy() );
  double lat = ell[0]*R2D;
  double hel = ell[2];

  _model->mops(lat, hel, ele, doy, zhd, zwd);  // DEG

#ifdef DEBUG
  cout << " MOPS: " << fixed << setprecision(3)
       << "  DOY: " << setw(12) << doy
       << "  ELE: " << setw(12) << ele
       << "  ZHD: " << setw(12) << zhd
       << endl;           
#endif

  return zhd;
}

// ---------
// MOPS ZWD interface
// ---------
double t_blindmops::getZWD(const t_gtriple& ell, const t_gtime& epo) // RAD
{
  gtrace("t_blindmops::getZWD");
      
  double zhd, zwd;
  double ele = 90.0;
  double doy = double( epo.doy() );
  double lat = ell[0]*R2D;
  double hel = ell[2];

  _model->mops(lat, hel, ele, doy, zhd, zwd); // DEG

#ifdef DEBUG
  cout << " MOPS: " << fixed << setprecision(3)
       << "  DOY: " << setw(12) << doy
       << "  ELE: " << setw(12) << ele
       << "  ZWD: " << setw(12) << zwd
       << endl;           
#endif

  return zwd;
}

} // namespace
