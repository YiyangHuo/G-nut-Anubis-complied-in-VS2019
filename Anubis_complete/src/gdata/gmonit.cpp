
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
#include <stdlib.h>

#include "gdata/gmonit.h" 
 
using namespace std;

namespace gnut {

/* ----------
 * constructor
 */
t_gmonit::t_gmonit( string id )
  : _moni_id(id) 
{}


/* ----------
 * destructor
 */
t_gmonit::~t_gmonit()
{}


/* ----------
 * basic class for monitoring purpose
 */
void t_gmonit::show( ostringstream& os, int verb )
{   
   os << _moni_id << " - method not implemented\n";
   return;
}

} // namespace
