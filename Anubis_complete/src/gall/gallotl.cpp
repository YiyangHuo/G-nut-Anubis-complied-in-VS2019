
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

#include "gall/gallotl.h"

namespace gnut {  

//Constructor
t_gallotl::t_gallotl()
{
//   cout << "vytvarim t_gallotl" << endl;
   id_type(  t_gdata::ALLOTL );
}

//Destructor
t_gallotl::~t_gallotl()
{
   _mapotl.clear();
}

int t_gallotl::data(Matrix& otldata, const string& site)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( _mapotl.find(site) == _mapotl.end() || _mapotl.size() == 0 ){   
//    cout << "t_gallotl: " << site << " not found " << _mapotl.size() << endl;
     _gmutex.unlock(); return -1; 
  }else otldata = _mapotl[site].data();

  _gmutex.unlock(); return 1;
}

double t_gallotl::lat(const string& site)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  double tmp = 0.0;
  if( _mapotl.find(site) != _mapotl.end() ) tmp = _mapotl[site].lat();
   
  _gmutex.unlock(); return tmp;
}

double t_gallotl::lon(const string& site)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  double tmp = 0.0;
  if( _mapotl.find(site) != _mapotl.end() ) tmp = _mapotl[site].lon();

  _gmutex.unlock(); return tmp;
}


// Add element to map container
// -----------------------------
void t_gallotl::add(t_gotl& otl)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  string site = otl.site();
  _mapotl[site] = otl;
   
  _gmutex.unlock(); return;
}

// Print all data
// --------------------
void t_gallotl::print()
{
   map<string, t_gotl>::iterator it;
   for (it = _mapotl.begin(); it != _mapotl.end(); it++)
   {
      t_gotl otl = it->second;
      cout << "Site: " << otl.site() << " Lon: " << otl.lon() << " Lat: " << otl.lat() << endl;
      cout << "Data: \n" << otl.data() << endl;
   }
   
}

} // namespace
