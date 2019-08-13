
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
#include <cmath>

#include "gmodels/gfunclin.h"

using namespace std;

namespace gnut {   

// function (linH0)
// ---------
int t_gfunc_linH0::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linH0)\n"; return -1;
  }
   
  for( int i=1; i<=dat.size(); ++i ){
   val(i) = fit(1) - fit(2)*dat(i);
  }
  return 0;
}
// ---------
int t_gfunc_linH0::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linH0)\n"; return -1;
  }
  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) =   1;
    val(i,2) = - dat(i);
  }
  return 0;
}


// function (linHS)
// ---------
int t_gfunc_linHS::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( _coeff.find("H0") == _coeff.end() || fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linHS)\n"; return -1;
  }

  double H0 = _coeff["H0"];

  for( int i=1; i<=dat.size(); ++i ){
   val(i) = fit(2) - fit(1)*(dat(i) - H0);
  }
  return 0;
}
// ---------
int t_gfunc_linHS::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( _coeff.find("H0") == _coeff.end() || fit.size() != 2 ){
    cout << "# warning: not properly initialized function (linHS)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = - dat(i);
    val(i,2) =   1;
  }
  return 0;
}


// function (linH1)
// ---------
int t_gfunc_linH1::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( _coeff.find("a0") == _coeff.end() ||
      _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
  {
    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
  }

  double a0 = _coeff["a0"];
  double H0 = _coeff["H0"];
   
  for( int i=1; i<=dat.size(); ++i ){
    val(i) = a0 - fit(1)*(dat(i)-H0);
  }
  return 0;
}
// ---------
int t_gfunc_linH1::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
  {
    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
  }

  double H0 = _coeff["H0"];
      
  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = - (dat(i)-H0);
  }
  return 0;
}

// function (linT)
// ---------
int t_gfunc_linT::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
//  if( _coeff.find("a0") == _coeff.end() ||
//      _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
//  {
//    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
//  }

//  double a0 = _coeff["a0"];
//  double H0 = _coeff["H0"];
   for( int i=1; i<=dat.size(); ++i )
     {
	val(i) = fit(1)+fit(2)*dat(i);
     }
   return 0;
}
// ---------
int t_gfunc_linT::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
//  if( _coeff.find("H0") == _coeff.end() || fit.size() != 1 )
//{
//    cout << "# warning: not properly initialized function (linH1)\n"; return -1;
//  }

//  double H0 = _coeff["H0"];
   for( int i=1; i<=dat.size(); ++i )
     {
	val(i,1) = 1;
	val(i,2) = dat(i);
     }    
  return 0;
}  
   
// function (linCS)
// ---------
int t_gfunc_linCS::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( fit.size() != 2 )
  {
    cout << "# warning: not properly initialized function (linCS)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i) = fit(1)*cos(dat(i)) + fit(2)*sin(dat(i));
  }

  return 0;
}
// ---------
int t_gfunc_linCS::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( fit.size() != 2 )
  {
    cout << "# warning: not properly initialized function (linCS)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = cos(dat(i));
    val(i,2) = sin(dat(i));
  }

  return 0;
}


// function (linCSH)
// ---------
int t_gfunc_linCSH::value( const ColumnVector& dat, const ColumnVector& fit, ColumnVector& val )
{
  if( fit.size() != 3 )
  {
    cout << "# warning: not properly initialized function (linCSH)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i) = ( fit(1)*cos(dat(i)) + fit(2)*sin(dat(i)) ) / fit(3);
  }

  return 0;
}
// ---------
int t_gfunc_linCSH::deriv( const ColumnVector& dat, const ColumnVector& fit, Matrix& val )
{
  if( fit.size() != 3 )
  {
    cout << "# warning: not properly initialized function (linCSH)\n"; return -1;
  }

  for( int i=1; i<=dat.size(); ++i ){
    val(i,1) = cos( dat(i) )/fit(3); // fit(3) = H_tropo
    val(i,2) = sin( dat(i) )/fit(3); // fit(3) = H_tropo
    val(i,3) = -(fit(1)*cos(dat(i)) + fit(2)*sin(dat(i))) /fit(3) /fit(3); // fit(3) = H_tropo
  }

  return 0;
}

} // namespace
