
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
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

#include "gdata/gnav.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gnav::t_gnav() : t_geph()
{
  id_type(t_gdata::EPH);
  id_group(t_gdata::GRP_EPHEM);
}


// destructor
// ----------
t_gnav::~t_gnav(){}


// get GNSS NAV validity time [s]
int t_gnav::nav_validity( GSYS gs )
{
  switch( gs ){
    case GPS : return MAX_GPS_TIMEDIFF;
    case GLO : return MAX_GLO_TIMEDIFF;
    case GAL : return MAX_GAL_TIMEDIFF;
    case BDS : return MAX_BDS_TIMEDIFF;
    case QZS : return MAX_QZS_TIMEDIFF;
    case SBS : return MAX_SBS_TIMEDIFF;
    case IRN : return MAX_IRN_TIMEDIFF;
    case GNS : return MAX_NAV_TIMEDIFF;
  }
  return MAX_NAV_TIMEDIFF;
}


// healthy check
// ----------
bool t_gnav::healthy() const
{   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  bool tmp = this->_healthy();
   
  _gmutex.unlock(); return tmp;
}


} // namespace
