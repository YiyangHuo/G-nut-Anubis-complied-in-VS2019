
#ifndef GTYPECONV_H
#define GTYPECONV_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: type conversion utilities
  Version: $ Rev: $

  2011-04-26 /PV: created
  2012-04-06 /JD: extracted type conversion utilities from previous utils.h

-*/

#include <cmath>
#include <string>
#include <iomanip>
#include <iostream>

using namespace std;

namespace gnut {
   
bool double_eq(const double&, 
               const double&);          // double equivalence according to machine capability
bool float_eq( const float&, 
               const float&);           // float equivalence according to machine capability
double dround(double d);                // round double 

string dbl2str(const double&, int p=3); // double to string conversion (width/digits by iomanip can be added !!!)
double str2dbl(const string&);          // string to double conversion (avoiding blanks)
double strSci2dbl(const string&);       // string (Scientific) to double conversion (including convert d,D to E!)
#ifdef STR2DBL
double str2dbl(const char*);            // faster char* to double conversion (avoiding blanks)
#endif

string int2str(const int&);             // integer to string conversion (widtht can be added !!!)
int    str2int(const string&);          // string to integer conversion (avoiding blanks)

size_t substitute(string& line, const string& a, const string& b, bool caseSensitive = true);

string ltrim(const string&);            // return trimmed string leading spaces
string rtrim(const string&);            // return trimmed string trailing spaces
string  trim(const string&);            // returm trimmed string leading and trailing spaces

double frac (double x);
string floatWoutZero(double x, int i = 3);
   
} // namespace
   
#endif