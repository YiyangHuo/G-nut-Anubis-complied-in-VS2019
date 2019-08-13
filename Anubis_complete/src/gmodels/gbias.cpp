
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

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "gmodels/gbias.h"

using namespace std;

namespace gnut {   
#define DIFF_DCB (86400*31)   // 31 days

// constructor
// ----------
t_gbias::t_gbias()
{
  _beg = FIRST_TIME;
  _end = LAST_TIME;
  
  id_type(t_gdata::BIAS);
  id_group(t_gdata::GRP_MODEL);
}


// destructor
// ----------
t_gbias::~t_gbias(){}


// get single code bias
// ----------
double t_gbias::bias(bool meter)
{
  if(meter) return _val * 1e-9 * CLIGHT;  
  else return _val;
}

// add bias
// ----------   
void t_gbias::set(t_gtime beg, t_gtime end, double val, GOBS obs1, GOBS obs2)
{
  #ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
  _gmutex.lock();
  
  _beg = beg;
  _end = end;
  
  _gobs = obs1;
  _ref  = obs2;
  
  _val   = val;
  
  _gmutex.unlock(); return;
}

// add bias
// ----------
void t_gbias::set( double val, GOBS obs1, GOBS obs2)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
  _gmutex.lock();
  
  _gobs = obs1;
  _ref  = obs2;
  
  _val   = val;
  
  _gmutex.unlock(); return;

}

// set reference signal
void t_gbias::ref(GOBS ref){
  _gmutex.lock();

  _ref = ref;

  _gmutex.unlock(); return;
}
  
// test validity   
bool t_gbias::valid(const t_gtime& epo)
{
  _gmutex.lock();
  
  bool ret = true;;

  if(epo < _beg || epo > _end) ret = false;
  else ret = true;
  
  _gmutex.unlock();
  return ret;
}

   
} // namespace
