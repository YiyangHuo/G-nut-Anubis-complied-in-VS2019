
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

#include <iostream>
#include <iomanip>

#include "gutils/gconst.h"

using namespace std;

namespace gnut {

// static map initializer
// ---------
t_map_refr REFR_COEF()
{
  t_map_refr m;
   
  m["ESSEN"]    = { K1_ESS, K2_ESS, K3_ESS };
  m["SMITH"]    = { K1_SMW, K2_SMW, K3_SMW };
  m["BEVIS"]    = { K1_BEV, K2_BEV, K3_BEV };
  m["THAYER"]   = { K1_THA, K2_THA, K3_THA };
  m["RUEGER"]   = { K1_RUE, K2_RUE, K3_RUE };
  m["FOELSCHE"] = { K1_FOE, K2_FOE, K3_FOE };

  return m;
}

} // namespace
