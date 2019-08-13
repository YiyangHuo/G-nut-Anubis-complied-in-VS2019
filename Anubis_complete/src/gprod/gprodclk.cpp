
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

#include "gprod/gprodclk.h"

using namespace std;

namespace gnut {
   
// constructor
// ----------
t_gprodclk::t_gprodclk(const t_gtime& t, shared_ptr<t_gobj> pt)
 : t_gprod(t, pt),
  _clk(0.0),
  _clk_rms(0.0)
{ 
   gtrace("t_gprodclk::t_gprodclk");
   
   id_type(CLK);
   id_group(GRP_PRODUCT);
}

// destructor
// --------------
t_gprodclk::~t_gprodclk()
{
}


// add clk
// --------------------------
void t_gprodclk::clk(const double& val, const double& rms)
{
   _gmutex.lock();
   
   _clk     = val;
   _clk_rms = rms;

   _gmutex.unlock();
   return;
}

// get clk
// --------------------------
double t_gprodclk::clk()
{
  return _clk;
}

// get clk rms
// -------------------------
double t_gprodclk::clk_rms()
{
  return _clk_rms;
}

   
// add ICB
// --------------------------
void t_gprodclk::icb(const double& val, const double& rms)
{
   _gmutex.lock();
   
   _icb     = val;
   _icb_rms = rms;

   _gmutex.unlock();
   return;
}

// get ICB
// --------------------------
double t_gprodclk::icb()
{
  return _icb;
}

// get ICB rms
// -------------------------
double t_gprodclk::icb_rms()
{
  return _icb_rms;
}


} // namespace