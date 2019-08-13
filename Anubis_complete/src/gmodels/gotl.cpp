
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

#include "gmodels/gotl.h"

namespace gnut {   

//Constructor
t_gotl::t_gotl()
{
  id_type(t_gdata::OTL);   
}

//Destructor
t_gotl::~t_gotl()
{
}

string t_gotl::site()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   string tmp = _site;
   _gmutex.unlock(); return tmp;
}

Matrix t_gotl::data()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   Matrix tmp = _data;
   _gmutex.unlock(); return tmp;
}

double t_gotl::lat()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   double tmp = _lat;
   _gmutex.unlock(); return tmp;
}

double t_gotl::lon()
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   double tmp = _lon;
   _gmutex.unlock(); return tmp;
}

// Set data
// -------------
void t_gotl::setdata(const string& site, const double& lon, const double& lat, const Matrix& data)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
   _site = site;
   _lon  = lon;
   _lat  = lat;
   _data = data;
   _gmutex.unlock();
   return;
}

} // namespace
