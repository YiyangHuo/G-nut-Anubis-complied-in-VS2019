
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
#include <algorithm>

#include "gio/glog.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// constructor
// ---------
t_glog::t_glog( string mask )
  : t_giof( mask ),
    _verb(0),
    _size(CACHE_LINES)
{}


// destructor
// ---------
t_glog::~t_glog(){}


// comment
// ---------
void t_glog::comment(int l, const string& str)
{ 
   if( fabs(l) > _verb ) return;
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_log_gmutex);
#endif
  _log_gmutex.lock();

   string text, text1;
   // negative will not print date !
   if( l < 0 ){
     text  = str + "\n";
     text1 = str + "\n"; 
   }else{
     text1 = text = t_gtime::current_time(t_gtime::LOC).str_ymdhms() + " " + str + "\n";

//   text1.replace(9,1," ");     // cache tens of seconds
     text1.replace(8,2,"  ");    // cache full minutes
//   text1.replace(6,4,"    ");  // cache tens of minutes
//   text1.replace(5,5,"     "); // cache full hours
   }
   
   // print repeating message maximally every minut
   if( find(_cache.begin(), _cache.end(), text1) == _cache.end() ){
       
     this->write( text.c_str(), text.size());
     this->flush(); // IMPORTANT ???

     // maintain cache
     _cache.push_back(text1);
     if( (int)_cache.size() > _size ) _cache.erase(_cache.begin());
   }
 _log_gmutex.unlock();
   return;
}


// comment & identificator
// ---------
void t_glog::comment(int l, const string& ide, const string& str)
{
  comment(l,"[" + ide + ":" + int2str(l) + "] " + str);
}

} // namespace