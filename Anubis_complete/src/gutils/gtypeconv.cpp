
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

#include <algorithm>
#include <sstream>
#include <cmath>
#include <limits>
#include <stdlib.h>

#include "gutils/gtypeconv.h"
#include "gutils/gconst.h"

using namespace std;

namespace gnut {

// double equivalence according to machine capability
// ----------
bool double_eq(const double& a, const double& b)
{
  if ( abs(a-b)< 1e-30) return true;   // must be editeed with machine epsilon (FOLLOWING FLOAT BELLOW!)
  else return false;
}

// float equivalence according to machine capability
// ----------   
bool float_eq(const float& a, const float& b)
{
  if ( abs(a-b)< std::numeric_limits<float>::epsilon() ) return true;
  else return false;
}

// round double
// ----------
double dround(double d)
{
  return floor(d+0.5);
}
	

// double to string conversion (width/digits by iomanip can be added !!!)
// ----------
string dbl2str(const double& d, int prec )
{
  ostringstream out;
  out << fixed << setprecision(prec) << " " << setw(0)  << d;
  return out.str();
}

#ifndef STR2DBL
// string to double conversion (avoiding blanks)
// ----------
double str2dbl(const string& s)
{
  return strtod( s.c_str(), NULL ); 
   
  // 2x slower alternative using stringstream !

//  double i;
//  istringstream istr(s);
//  istr >> i;
//  return i;
}
#else
// http://tinodidriksen.com/2011/05/28/cpp-convert-string-to-double-speed/
// http://pastebin.com/dHP1pgQ4
// bool naive(T & r, const char *p)
double str2dbl(const string& s)
{ 
  return str2dbl( s.c_str() );
}
// string to double conversion (avoiding blanks)      // http://tinodidriksen.com/2011/05/28/cpp-convert-string-to-double-speed/
// ----------
double str2dbl(const char* p)
{ 

#define white_space(c) ((c) == ' ' || (c) == '\t')
#define valid_digit(c) ((c) >= '0' && (c) <= '9')
#define null_double    0.0;

  //Skip leading white space, if any.
  while( white_space(*p) ) p += 1;
 	  
  double r = 0.0;
  int    c = 0; // counter to check how many numbers we got!
	  
  // Get the sign!
  bool neg = false;
  if(      *p == '-' ){ neg = true;  ++p; }
  else if( *p == '+' ){ neg = false; ++p; }
	  
  // Get the digits before decimal point
  while( valid_digit(*p) ){ r = (r*10.0) + (*p - '0'); ++p; ++c; }
	  
  // Get the digits after decimal point
  if( *p == '.' ){ 
    double f = 0.0;
    double scale = 1.0;
    ++p;
    while( *p >= '0' && *p <= '9' ){
        f = (f*10.0) + (*p - '0');
        ++p;
        scale*=10.0;
        ++c;
    }
    r += f / scale;
  }
	  
  // FIRST CHECK:
  if( c==0 ){ return null_double; } // we got no decimal places! this cannot be any number!
	  
  // Get the digits after the "e"/"E" (exponent)
  if( *p == 'e' || *p == 'E' || *p == 'd' || *p == 'D' ){
    int e = 0;

    bool negE = false;
    ++p;
    if(      *p == '-'){ negE = true;  ++p; }
    else if( *p == '+'){ negE = false; ++p; }

    // Get exponent
    c = 0;
    while( valid_digit(*p) ){
      e = (e*10) + (*p - '0');
      ++p; ++c;
    }
    
    if( !neg && e>std::numeric_limits<double>::max_exponent10 ){ e = std::numeric_limits<double>::max_exponent10; }
    else if(    e<std::numeric_limits<double>::min_exponent10 ){ e = std::numeric_limits<double>::max_exponent10; }

    // SECOND CHECK:
    if( c==0 ) return null_double; // we got no  exponent! this was not intended!!
	  
    double scaleE = 1.0;
   
    // Calculate scaling factor.	  
    while( e >= 50 ){ scaleE *= 1E50; e -= 50; }
  //while (e >=  8) { scaleE *= 1E8;  e -=  8; }
    while (e >   0) { scaleE *= 10.0; e -=  1; }
 
    if( negE ){ r /= scaleE; }
    else{       r *= scaleE; }
  }
	  
  // POST CHECK:
  // skip post whitespaces
  while( white_space(*p) ) ++p;

  if( *p != '\0' ){ return null_double; } // if next character is not the terminating character
	  
  // Apply sign to number
  if( neg ) r = -r;
	  
  return r;
}
#endif


// string (scientific) to double conversion (avoiding blanks)
// ----------
double strSci2dbl(const string& s)
{   
  double i = 0.0;
  string str(s);
  substitute(str,"d","E");
  substitute(str,"D","E");
  istringstream istr(str);
  istr >> i;
  return i;
}


// integer to string conversion (widtht can be added !!!)
// ----------
string int2str(const int& i)
{
  ostringstream out;
  out << i;
  return out.str();
}


// string to integer conversion (avoiding blanks)
// ----------
int str2int(const string& s)
{
  int i = 0;
  istringstream istr(s);
  istr >> i;
  return i;
}


// return trimmed trailing spaces and changes the content
// ----------
string rtrim(const string& s)
{   
  string str;
  size_t endpos = s.find_last_not_of(" \t");
  if( string::npos != endpos ) str = s.substr( 0, endpos+1 );
  return str;
}

   
// return trimmed string leading spaces
// ----------
string ltrim(const string& s)
{   
  string str;
  size_t startpos = s.find_first_not_of(" \t");
  if( string::npos != startpos ) str = s.substr( startpos );
  return str;
}


// return trimmed leading and tailing spaces
// ----------
string trim(const string& s)
{
  return ltrim( rtrim(s) );
}


// string substitute
// ----------
size_t substitute(string& line, const string& a, const string& b, bool caseSensitive )
{   
   size_t n = 0;
     
   if( caseSensitive ){
     string tmp(line);
     while( (n = line.find(a)) != string::npos ){
       tmp = line.substr(0, n) + b + line.substr( n + a.length() );
       line = tmp;
     }

   }else{
     string lineLC(line);
     string findLC(a);
     transform(lineLC.begin(), lineLC.end(), lineLC.begin(), ::tolower);     
     transform(findLC.begin(), findLC.end(), findLC.begin(), ::tolower);

     while( (n = lineLC.find(findLC)) != string::npos ){
       lineLC = line.substr(0, n) + b + line.substr( n + a.length() );
       line   = lineLC;
     }
   }

   return n+b.length();
}

// Frac part of double
// ------------------------
double frac (double x)
{
    return x-floor(x);
}

// Get Float without Leading Zero
// ------------------------
string floatWoutZero(double x, int l)
{
   stringstream ss;
   ss << std::setw(1+l) << std::setprecision(l);
   ss << std::fixed << x;

   string str = ss.str();
   if(x > 0.f && x < 1.f)
      return str.substr(1, str.size()-1);
    else if(x < 0.f && x > -1.f)
      return "-" + str.substr(2, str.size()-1);
   
  return str;
}

} // namespace
