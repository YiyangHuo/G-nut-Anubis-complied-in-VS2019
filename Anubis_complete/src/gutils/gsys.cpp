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
#include <iostream>
#include <iomanip>
#include <sstream>

#include "gutils/gsys.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

// ---------------------------------------------------------------------------------
// class GSYS
// ---------------------------------------------------------------------------------

// constructor
// ---------
t_gsys::t_gsys(GSYS sys) : _gsys(sys)
{}

// constructor
// ---------
t_gsys::t_gsys(string sys)
 : _gsys(GNS)
{
 if( sys.size() > 0 ) _gsys = str2gsys(sys);
}

// constructor
// ---------
t_gsys::t_gsys(char c)
{
  _gsys = char2gsys(c);
}


// --------------------------------------------------------
// STATIC FUNCTIONS
// --------------------------------------------------------
// get system band
// ----------
GOBSBAND t_gsys::gfrq2band(GSYS gs, GFRQ gfreq )
{
  t_map_freq m = GNSS_FREQ_PRIORITY();
  if( gs == GNS || gfreq == LAST_GFRQ ) return BAND;

  for( size_t i=0; i < m[gs].size(); ++i ){
    if( m[gs].at(i) == gfreq ) return band_priority(gs, FREQ_SEQ(i));
  }
  return BAND;
}

// get system freq
// ----------
GFRQ t_gsys::band2gfrq(GSYS gs, GOBSBAND gband )
{
  t_map_band m = GNSS_BAND_PRIORITY();
  if( gs == GNS || gband == BAND ) return LAST_GFRQ;

  for( size_t i=0; i < m[gs].size(); ++i ){
    if( m[gs].at(i) == gband ) return freq_priority(gs, FREQ_SEQ(i));
  }
  return LAST_GFRQ;
}

// get system sequence
// ----------
FREQ_SEQ t_gsys::gfrq2freq(GSYS gs, GFRQ gfreq )
{
  t_map_freq m = GNSS_FREQ_PRIORITY();
  if( gs == GNS || gfreq == LAST_GFRQ ) return FREQ_X;

  for( size_t i=0; i < m[gs].size(); ++i ){
    if( m[gs].at(i) == gfreq ) return FREQ_SEQ(i);
  }
  return FREQ_X;
}

// get system freq sequence
// ----------
FREQ_SEQ t_gsys::band2freq(GSYS gs, GOBSBAND gband )
{
  t_map_band m = GNSS_BAND_PRIORITY();
  if( gs == GNS || gband == BAND ) return FREQ_X;

  for( size_t i=0; i < m[gs].size(); ++i ){
    if( m[gs].at(i) == gband ) return FREQ_SEQ(i);
  }
  return FREQ_X;
}

   
// get band selection
// ----------
GOBSBAND t_gsys::band_priority(GSYS gs, FREQ_SEQ iseq)  // FREQ_SEQ priority band
{
  t_map_band m = GNSS_BAND_PRIORITY();

  if( gs == GNS || iseq >= m[gs].size() ) return BAND;

  return m[gs][iseq];
}

// get freq selection
// ----------
GFRQ t_gsys::freq_priority(GSYS gs, FREQ_SEQ iseq)  // FREQ_SEQ priority frequency
{
  t_map_freq m = GNSS_FREQ_PRIORITY();

  if( gs == GNS || iseq >= m[gs].size() ) return LAST_GFRQ;

  return m[gs][iseq];
}

// get attr selection
// ----------
GOBSATTR t_gsys::attr_priority(GSYS gs, GOBSBAND gb, GOBSTYPE gt, unsigned int iseq)  // iseq priority sequence
{
  t_map_gnss m = GNSS_DATA_PRIORITY();

  if( gs == GNS  ||
      gt == TYPE ||
      gb == BAND || iseq > m[gs][gb][gt].size() ) return ATTR;

  return m[gs][gb][gt][iseq];
}


// convert GSYS enum to GSYS string
// ---------- 
string t_gsys::gsys2str( GSYS sys )
{
  gtrace("t_gsys::gsys2str");

// PROBLEM WITH STATIC FUNCTION
//  boost::mutex::scoped_lock lock(_mutex);

  switch( sys ){
   case GPS  : return "GPS";
   case GLO  : return "GLO";
   case GAL  : return "GAL";
   case BDS  : return "BDS";
   case SBS  : return "SBS";
   case QZS  : return "QZS";
   case IRN  : return "IRN";
   case GNS  : return "GNS";
     
   default : { cout << "*** warning: unknown GNSS system!\n"; cout.flush(); }
   
        // cannot be while static function !!!
        // if( _log ) _log->comment(0,"gobsgnss","warning: unknown frequency code ");
        //     else 
  }
   
  return "GNS";
}


// convert GSYS enum to GSYS char
// ---------- 
char t_gsys::gsys2char( GSYS sys )
{
  gtrace("t_gsys::gsys2char");

// PROBLEM WITH STATIC FUNCTION
//  boost::mutex::scoped_lock lock(_mutex);

  switch( sys ){
   case GPS  : return 'G';
   case GLO  : return 'R';
   case GAL  : return 'E';
   case BDS  : return 'C';
   case SBS  : return 'S';
   case QZS  : return 'J';
   case IRN  : return 'I';
   case GNS  : return 'X';
     
   default : { cout << "*** warning: unknown GNSS system \n"; cout.flush(); }
        // cannot be while static function !!!
        // if( _log ) _log->comment(0,"gobsgnss","warning: unknown frequency code ");
        //     else 
  }
   
  return 'X';
}


// convert GSYS string to GSYS enum
// ---------- 
GSYS t_gsys::str2gsys( string s )
{
  gtrace("t_gsys::str2gsys [" + s + "]" );
   
// PROBLEM WITH STATIC FUNCTION
//  boost::mutex::scoped_lock lock(_mutex);

  if( s.size() == 0 ){
    cout << "*** warning: not defined GNSS system code [NULL]\n"; cout.flush();
    return GNS;
  }

  transform(s.begin(), s.end(), s.begin(), ::toupper);
  
  if(      s == "G"    || s == "GPS"  || s == "NAVSTAR"   ) return GPS;
  else if( s == "R"    || s == "GLO"  || s == "GLONASS"   ) return GLO;
  else if( s == "E"    || s == "GAL"  || s == "GALILEO"   ) return GAL;
  else if( s == "C"    || s == "COMP" || s == "COMPASS"   ) return BDS;
  else if( s == "C"    || s == "BDS"  || s == "BEIDOU"    ) return BDS;
  else if( s == "S"    || s == "SBS"  || s == "EGNOS"     ) return SBS;
  else if( s == "S"    || s == "SBAS"                     ) return SBS;
  else if( s == "J"    || s == "QZS"  || s == "JAXA"      ) return QZS;
  else if( s == "J"    || s == "QZSS"                     ) return QZS;
  else if( s == "I"    || s == "IRN"  || s == "IRNSS"     ) return IRN;
  else if( s == "X"    || s == "GNS"  || s == "GNSS"      ) return GNS;
  else if( s == "M"                                       ) return GNS;
  else{ cout << "*** warning: not defined GNSS system code [" << s[0] << "]\n"; cout.flush(); }

  return GNS;
}


// convert GSYS string to GSYS enum
// ---------- 
char t_gsys::str2char(string s)
{
  if( s.size() > 0 ) return gsys2char( str2gsys(s) );
  return 'X';
}	


// convert GSYS char to GSYS enum
// ---------- 
GSYS t_gsys::char2gsys(char c)
{
  gtrace("t_gsys::char2gsys");

  if(      c == 'G' ) return GPS;
  else if( c == 'R' ) return GLO;
  else if( c == 'E' ) return GAL;
  else if( c == 'C' ) return BDS;
  else if( c == 'S' ) return SBS;
  else if( c == 'J' ) return QZS;
  else if( c == 'I' ) return IRN;
  else if( c == 'M' ) return GNS;
  else if( c == 'X' ) return GNS;
  else{ cout << "*** warning: not defined GNSS system char [" << c << "]\n"; cout.flush(); }

  return GNS;
}


// convert GSYS char to GSYS enum
// ---------- 
string t_gsys::char2str(char c)
{
  return gsys2str( char2gsys(c) );
}	


// convert satellite name
// ---------- 
string t_gsys::eval_sat( string sat, GSYS sys )
{
  istringstream is( sat ); // .substr(0,3) NE!
  GSYS gnss = sys;
  int  svn  =  0;
  char chr  = 'G';        // if empty, use GPS
  size_t l  = is.str().length();
   
  if( l == 0 ) return "X00";

  if( l < 3 || sat[0] == ' ' )  is        >> svn;
  else                          is >> chr >> svn;

  if( is.fail() ){ return "X00"; }
  if( chr != 'G' ) gnss = char2gsys( chr );

  switch( gnss ){
    case GPS       : chr = 'G'; break;
    case GLO       : chr = 'R'; break;
    case GAL       : chr = 'E'; break;
    case BDS       : chr = 'C'; break;
    case SBS       : chr = 'S'; break;
    case QZS       : chr = 'J'; break;
    case IRN       : chr = 'I'; break;
    case GNS       : chr = 'X'; break;
    default        : chr = 'X'; break;
  }

  if( svn > QZS_OFFSET ){ svn -= QZS_OFFSET; }
  if( svn > SBS_OFFSET ){ svn -= SBS_OFFSET; }

  ostringstream os;
  os << setw(1) << chr << setfill('0') << setw(2) << svn;

#ifdef DEBUG3
  cout << "  ARGM: [" << sat << "] " << gsys2str( sys )
       << "  EVAL: " << gsys2str( gnss ) << " [" << chr << "] " << svn
       << "  RESU: " << os.str() << " " << l << endl;  cout.flush();
#endif
  return os.str();
}


// convert satellite name
// ---------- 
string t_gsys::eval_sat( int svn, GSYS sys )
{
  char chr;
  switch( sys ){
    case GPS       : chr = 'G'; break;
    case GLO       : chr = 'R'; break;
    case GAL       : chr = 'E'; break;
    case BDS       : chr = 'C'; break;
    case SBS       : chr = 'S'; break;
    case QZS       : chr = 'J'; break;
    case IRN       : chr = 'I'; break;
    case GNS       : chr = 'X'; break;
    default        : chr = 'X'; break;
  }
  if( svn > QZS_OFFSET ){ svn -= QZS_OFFSET; }
  if( svn > SBS_OFFSET ){ svn -= SBS_OFFSET; }
    
  ostringstream os;
  os << setw(1) << chr << setfill('0') << setw(2) << svn;

  return os.str();
}

// set GSYS from string
// ---------- 
void t_gsys::from_string( string sys )
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gsys = str2gsys( sys );
}


// set GSYS from enum
// ---------- 
void t_gsys::from_gsys( GSYS sys )
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gsys = sys;
}


// overloaded equivalence operator
// ----------
bool t_gsys::operator==(const string& sys) const
{
//  transform(sys.begin(), sys.end(), sys.begin(), ::toupper);
  if( _gsys == str2gsys(sys)  ) return true;
  return false;
}


// overloaded equivalence operator
// ----------
bool t_gsys::operator==(const GSYS& sys) const
{
//  if( _gsys == sys.gsys() ) return true;
  return false;
}

/*
// decode band number to frequency number
// -----------------------------------
FREQ_SEQ t_gsys::band2freqID(GSYS gs, GOBSBAND band)
{
  gtrace("t_gsys::band2freqID");

   int freq = 0;
   
   switch (gs) {
    case GPS :
      if(band == BAND_1) freq = 1;
      if(band == BAND_2) freq = 2;
      if(band == BAND_5) freq = 3;
      break;		          
    case GLO :		          
      if(band == BAND_1) freq = 1;
      if(band == BAND_2) freq = 2;
      if(band == BAND_3) freq = 3;
      if(band == BAND_7) freq = 3;     
      if(band == BAND_5) freq = 4;
      break;      	          
    case GAL:		          
      if(band == BAND_1) freq = 1;     
      if(band == BAND_6) freq = 2;
      if(band == BAND_7) freq = 3;
      if(band == BAND_8) freq = 4;
      if(band == BAND_5) freq = 5;
      break;		          
    case BDS:		          
      if(band == BAND_2) freq = 1;
      if(band == BAND_6) freq = 2;
      if(band == BAND_7) freq = 3;
      break;		          
    case QZS:		          
      if(band == BAND_1) freq = 1;
      if(band == BAND_6) freq = 2;
      if(band == BAND_2) freq = 3;
      if(band == BAND_5) freq = 4;
      break;		          
    case SBS:		          
      if(band == BAND_1) freq = 1;
      if(band == BAND_2) freq = 2;
      if(band == BAND_5) freq = 3;
      break;
    default:
      break;
   }   
   return freq;
}

// decode frequency number to band number
// -----------------------------------
GOBSBAND t_gsys::freq2band(GSYS gs, int freq)
{
  gtrace("t_gsys::freq2band");

//   GSYS gs = this->gsys();
   
   GOBSBAND band = BAND;
   
   switch (gs) {
    case GPS :
      if(freq == 1) band = BAND_1;
      if(freq == 2) band = BAND_2;
      if(freq == 3) band = BAND_5;
      break;
    case GLO :
      if(freq == 1) band = BAND_1;
      if(freq == 2) band = BAND_2;
      if(freq == 3) band = BAND_3;
      if(freq == 4) band = BAND_5;
      break;      
    case GAL:
      if(freq == 1) band = BAND_1;
      if(freq == 2) band = BAND_6;
      if(freq == 3) band = BAND_7;
      if(freq == 4) band = BAND_8;
      if(freq == 5) band = BAND_5;
      break;
    case BDS:
      if(freq == 1) band = BAND_2;
      if(freq == 2) band = BAND_6;
      if(freq == 3) band = BAND_7;
      break;
    case QZS:
      if(freq == 1) band = BAND_1;
      if(freq == 2) band = BAND_6;
      if(freq == 3) band = BAND_2;
      if(freq == 4) band = BAND_5;
      break;
    case SBS:
      if(freq == 1) band = BAND_1;
      if(freq == 2) band = BAND_2;
      if(freq == 3) band = BAND_5;      
      break;
    default:
      break;
   }   
   return band;
}
*/


// ---------------------------------------------------------------------------------
// class FREQ
// ---------------------------------------------------------------------------------

// get GFRQ enum from string
// ----------
GFRQ t_gfreq::str2gfreq( string freq )
{
  gtrace("t_gsys::str2gfreq");
   
// PROBLEM WITH STATIC FUNCTION
//  boost::mutex::scoped_lock lock(_mutex);

  transform(freq.begin(), freq.end(), freq.begin(), ::toupper);
   
  if(      freq == "G01" ) return G01;
  else if( freq == "G02" ) return G02;
  else if( freq == "G05" ) return G05;
  else if( freq == "R01" ) return R01;
  else if( freq == "R02" ) return R02;
//else if( freq == "R01" ) return R01_CDMA;
//else if( freq == "R02" ) return R02_CDMA;
  else if( freq == "R03" ) return R03_CDMA;
  else if( freq == "R05" ) return R05_CDMA;
  else if( freq == "E01" ) return E01;
  else if( freq == "E05" ) return E05;
  else if( freq == "E07" ) return E07;
  else if( freq == "E08" ) return E08;
  else if( freq == "E06" ) return E06;
  else if( freq == "C02" ) return C02;
  else if( freq == "C06" ) return C07;
  else if( freq == "C07" ) return C06;
  else if( freq == "J01" ) return J01;
  else if( freq == "J02" ) return J02;
  else if( freq == "J05" ) return J05;
  else if( freq == "J06" ) return J06;
  else if( freq == "S01" ) return S01;
  else if( freq == "S05" ) return S05;
  else if( freq == "I05" ) return I05;
//else if( freq == "I09" ) return I09;
//else if( freq.empty()        ) return LAST_GFRQ;
  else if( freq == "LAST_GFRQ" ) return LAST_GFRQ;
  
//  cout << "*** warning: not defined frequency code [" + freq + "]\n"; cout.flush();

  // cannot be while static function !!!
  //  if( _log ) _log->comment(0,"gobsgnss","warning: not defined frequency code " + freq);
  //  else                      cout << "*** warning: not defined frequency code " + freq << endl;

  return LAST_GFRQ;
}


// get string from GFRQ enum
// ---------- 
string t_gfreq::gfreq2str( GFRQ freq )
{
  gtrace("t_gsys::gfreq2str");

// PROBLEM WITH STATIC FUNCTION
//  boost::mutex::scoped_lock lock(_mutex); 
   
  switch( freq ){
   case G01 : return "G01";
   case G02 : return "G02";
   case G05 : return "G05";
   case R01 : return "R01";
   case R02 : return "R02";
   case R01_CDMA : return "R01";
   case R02_CDMA : return "R02";
   case R03_CDMA : return "R03";
   case R05_CDMA : return "R05";
   case E01 : return "E01";
   case E05 : return "E05";
   case E07 : return "E07";
   case E08 : return "E08";
   case E06 : return "E06";
   case C02 : return "C02";
   case C06 : return "C07";
   case C07 : return "C06";
   case J01 : return "J01";
   case J02 : return "J02";
   case J05 : return "J05";
   case J06 : return "J06";
   case S01 : return "S01";
   case S05 : return "S05";
   case I05 : return "I05";
// case I09 : return "I09";          
   case LAST_GFRQ : return "LAST_GFRQ";
     
   default : {  cout << "*** warning: unknown frequency code \n"; cout.flush(); }
   
        // cannot be while static function !!!
        // if( _log ) _log->comment(0,"gobsgnss","warning: unknown frequency code ");
        //     else 
  }
   
  return "LAST_GFRQ";
}
   
// get true if BDS GEO satellite
// ------------------------------
bool t_gsys::bds_geo(const string& sat)
{
   set<string> geo;                                  // prns of geost. satellites
   geo.insert("C01");
   geo.insert("C02");
   geo.insert("C03");
   geo.insert("C04");
   geo.insert("C05");      
   if( geo.find(sat) != geo.end() ) return true;
     
   return false;
}
   
} // namespace
