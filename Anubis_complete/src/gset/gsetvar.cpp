
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
#include <sstream>
#include <algorithm>

#include "gset/gsetvar.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetvar::t_gsetvar()
 : t_gsetbase()
{
   _set.insert(XMLKEY_VAR);
}


// Destructor
// -----------
t_gsetvar::~t_gsetvar()
{}


// Return value
// ----------
string t_gsetvar::get_attr(string key)
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_VAR).attribute(key.c_str()).as_string();

  _gmutex.unlock(); return tmp;
}


// Set value
// ----------
void t_gsetvar::set_attr(string key, string val)
{
  _gmutex.lock();

  // check/set attribute
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node   = _default_node(parent, XMLKEY_VAR);

  _default_attr(node, key.c_str(), val);

  _gmutex.unlock();
}


} // namespace
