
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

#include "gutils/gobs.h"

using namespace std;

namespace gnut {

// -------------------------------------------------------------------------------------------
// class T_GATTR
// -------------------------------------------------------------------------------------------
bool t_gattr::valid() const
{
  return ( _gattr != ATTR );
}

// set attr
// ----------
void t_gattr::attr( const GOBSATTR& a )
{
  _gattr = a;
}

// get attr
// -----------
GOBSATTR t_gattr::attr()const
{ 
  return _gattr;
}

// operator
// ----------
bool t_gattr::operator==(const t_gattr& g) const
{
  return ( _gattr == g.attr() );
}


// -------------------------------------------------------------------------------------------
// class T_GBAND
// -------------------------------------------------------------------------------------------
bool t_gband::valid() const
{
  return ( t_gattr::valid() && _gband != BAND );
}

// set band
// ----------
void t_gband::band( const GOBSBAND& b )
{
  _gband = b;
}

// get band
// -----------
GOBSBAND t_gband::band()const
{ 
  return _gband;
}

// set attr
// ----------
void t_gband::gattr( const t_gattr& g )
{
  _gattr = g.attr();
  _gband = BAND;
}

// get attr
// -----------
t_gattr t_gband::gattr()const
{
  t_gattr g(_gattr);
  return g;
}

// operators
// ----------
bool t_gband::operator==(const t_gband& g) const
{
  return ( _gband == g.band() &&
	   _gattr == g.attr()
         );
}

// -------------------------------------------------------------------------------------------
// class T_GOBS
// -------------------------------------------------------------------------------------------

// valid ?
// -----------
bool t_gobs::valid() const
{
  return ( t_gband::valid() && _gtype != TYPE );
}

// set type
// ----------
void t_gobs::type( const GOBSTYPE& t )
{
  _gtype = t;
}

// get type
// -----------
GOBSTYPE t_gobs::type()const
{
  return _gtype;
}

// set gband
// ----------
void t_gobs::gband( const t_gband& g )
{
  _gattr = g.attr();
  _gband = g.band();
  _gtype = TYPE;
}

// get gband
// -----------
t_gband t_gobs::gband()const
{
  t_gband g(_gband,_gattr);
  return g;
}

// operator
// ----------
bool t_gobs::operator==(const t_gobs& g) const
{
  return ( _gtype == g.type() &&
	   _gband == g.band() &&
           _gattr == g.attr()
	 );
}

// set from GOBS
// -----------
int t_gobs::gobs(const GOBS& g)
{
  string s = gobs2str( g );
  _gtype = str2gobstype( s );
  _gband = str2gobsband( s );
  _gattr = str2gobsattr( s );

  return 1;
}

// set from string
// -----------
int t_gobs::gobs(const string& s)
{
  _gtype = str2gobstype( s );
  _gband = str2gobsband( s );
  _gattr = str2gobsattr( s );

  return 1;
}

// get GOBS enum
// -----------
GOBS t_gobs::gobs()const
{
  string s = gobstype2str(_gtype) + 
             gobsband2str(_gband) + 
             gobsattr2str(_gattr);
  return str2gobs( s );
}

// get true if code observation
// -----------
bool t_gobs::is_code()const
{
  return( _gtype == TYPE_C || _gtype == TYPE_P );
}


// get true if phase observation
// -----------
bool t_gobs::is_phase()const
{
  return( _gtype == TYPE_L );
}

} // namespace
