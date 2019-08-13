
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

//#include <cstring>
//#include <iostream>

#include "gio/gnote.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// constructor
// ---------
t_gnote::t_gnote(t_note n, string s)
{
  _stat = n;
  _text = s;
}


// destructor
// ---------
t_gnote::~t_gnote()
{}


// overloading << operator
// -----------------------------
ostream& operator<<(ostream& os, const t_gnote& n)
{
  os << n.str();
  return os;
}


// overloading == operator
// -----------------------------
bool t_gnote::operator==(const t_gnote& n) const
{
//  boost::mutex::scoped_lock lock(_mutex_triple);
  return ( n.status() == _stat &&
	   n.text()   == _text   
	 );
}


// overloading < operator
// -----------------------------
bool t_gnote::operator<(const t_gnote& n) const
{
//  boost::mutex::scoped_lock lock(_mutex_triple);
  return ( n.status() < _stat &&
	   n.text()   < _text
	 );
}



// get string
// ---------
string t_gnote::_str()const
{ 
  string head;
//  switch( _stat ){
//    case MESSAGE : head = "Message:"; break;
//    case WARNING : head = "Warning:"; break;
//  case ERROR   : head = "Error:";   break;
//  }
   
  return head + _text;
}


} // namespace