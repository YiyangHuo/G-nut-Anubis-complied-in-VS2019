
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

#include "gdata/geph.h"
#include "gutils/gtypeconv.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_geph::t_geph()
 : t_gdata(),
   _sat(""),
   _epoch(t_gtime::GPS),
   _validity(true),
   _gio_ptr(0)
{
  _type  = EPHGPS;
}


// destructor
// ----------
t_geph::~t_geph(){
}


// gsys
// ----------
GSYS t_geph::gsys() const
{   
   gtrace("t_geph::gsys");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  GSYS tmp = GNS;
  _gmutex.lock();
  if( _valid() ) tmp = t_gsys::char2gsys( _sat[0] );
  _gmutex.unlock();
  return tmp;
}


// gsat
// ----------
string t_geph::gsat() const
{   
   gtrace("t_geph::gsat");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  string tmp = "";
  _gmutex.lock();
  if( _valid() ) tmp = _sat;
  _gmutex.unlock();
  return tmp;
}


// get parameter value
// ----------
t_timdbl t_geph::param( const NAVDATA& n )
{
  t_timdbl tmp;
  return tmp;
}


// get parameter value
// ----------
int t_geph::param( const NAVDATA& n, double val )
{
  return 0;
}


// get parameter value
// ----------
bool t_geph::param_cyclic( const NAVDATA& n )
{
  if( n == 2 || n == 6 ) return true;
  return false;
}

   
// valid
// ----------
bool t_geph::valid()
{
   gtrace("t_geph::valid");
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  bool tmp = this->_valid();
  _gmutex.unlock();  
  return tmp;
}


// clean data
// ----------
void t_geph::clear()
{ 
   gtrace("t_geph::clear");   
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  this->_clear();
  _gmutex.unlock();
   return;
}


// clean internal function
// ----------
void t_geph::_clear()
{
   gtrace("t_geph::_clear");   
   
  _sat.clear();
  _epoch = FIRST_TIME;
}


// clean internal function
// ----------
bool t_geph::_valid() const
{   
   gtrace("t_geph::_valid");   
   
  if( !_validity   ||
      _sat.empty() || 
      _sat   == "" ||
      _epoch == FIRST_TIME ) return false;

  return true; 
}


// type set
// ----------
//void t_geph::set_types( unsigned int types ){
//  _types = types;
//}


// return type
// ----------
//unsigned int t_geph::types() const
//{
// if( (_types & EPH_NAV ) == EPH_NAV  ) cout << "_types.nav   \n";
// if( (_types & EPH_RTCM) == EPH_RTCM ) cout << "_types.rtcm  \n";
// if( (_types & EPH_SP3 ) == EPH_SP3  ) cout << "_types.sp3   \n";
// if( (_types & EPH_CLK ) == EPH_CLK  ) cout << "_types.clk   \n";
// return _types;
//}

} // namespace
