
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

#include "gutils/gmatrixconv.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"
#include "gmodels/geop.h"

using namespace std;

namespace gnut {   

//Constructor
t_geop::t_geop()
{
}

//Destructor
t_geop::~t_geop()
{
}

// Nutation Matrix 
// --------------------------------------
Matrix t_geop::nutMatrix(double mjd)
{
//  boost::mutex::scoped_lock lock(_mutex);      

  Matrix I(3,3);
  I << 1 << 0 << 0
    << 0 << 1 << 0
    << 0 << 0 << 1;
     
  return I;
 
}

// Precession Matrix
// ------------------------------------------
Matrix t_geop::precMatrix (double mjd_1)
{
//  boost::mutex::scoped_lock lock(_mutex);
   
  Matrix I(3,3);
  I << 1 << 0 << 0
    << 0 << 1 << 0
    << 0 << 0 << 1;
     
  return I;   
}

// Normalize angle into interval 0 - 2pi
// -----------------------------------
double t_geop::_normangle(double x)
{
   double norm;
   norm = fmod(x, 2*G_PI);
   if (norm < 0 ) norm += 2*G_PI;
   
   return norm;
}

} // namespace
