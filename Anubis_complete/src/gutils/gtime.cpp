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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <thread>
#include <chrono>
#include <algorithm>

#if  defined _WIN32 || defined _WIN64
#include <io.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {


// constructor 
// ----------
t_gtime::t_gtime(const t_tsys& ts)
  : _mjd(0), _sod(0), _dsec(0.0), _tsys(UTC)
{
  from_time(time(NULL), 0.0, true);  // input UTC time-system
  _tsys = ts;                        // now switch to required TS
}


// constructor 
// ----------
t_gtime::t_gtime(const time_t& tt, const double& ds, const t_tsys& ts )
  : _mjd(0), _sod(0), _dsec(0.0), _tsys(UTC)
{
  from_time(tt, ds, true); // input UTC time-system
  _tsys = ts;              // now switch to required TS
}


// constructor 
// ----------
t_gtime::t_gtime(const int& yr, const int& mn, 
                 const int& dd, const int& hr, 
                 const int& mi, const int& sc, 
                 const double& ds, const t_tsys& ts )
  : _mjd(0), _sod(0), _dsec(0.0), _tsys(ts)
{
  from_ymdhms(yr,mn,dd,hr,mi,sc + ds,true);
  // old style
  // from_ymd(yr,mn,dd,static_cast<int>(hr*3600 + mi*60 + sc),ds,true);
}


// constructor
// ----------
t_gtime::t_gtime(const int& gw, const int& dw, const int& sd, 
                 const double& ds, const t_tsys& ts)
  : _mjd(0), _sod(0), _dsec(0.0), _tsys(ts)
{
  from_gwd(gw, dw, sd, ds, true);
}


// constructor
// ----------
t_gtime::t_gtime(const int& gw, const double& sow, const t_tsys& ts)
  : _mjd(0), _sod(0), _dsec(0.0), _tsys(ts)
{
  from_gws(gw, sow, true);
}


// constructor
// ----------
t_gtime::t_gtime(const int& mjd, const int& sd,
                 const double& ds, const t_tsys& ts)
  : _mjd(0), _sod(0), _dsec(0.0), _tsys(ts)
{ 
 from_mjd(mjd, sd, ds, true);
}


// destructor
// ----------
t_gtime::~t_gtime()
{}


// convert string to tsys
// ----------
string t_gtime::tsys2str(t_tsys ts)
{
  switch( ts ){
   case USER : return "USER";
   case TAI  : return "TAI";
   case UTC  : return "UTC";
   case LOC  : return "LOC";
   case GPS  : return "GPS";
   case GLO  : return "GLO";
   case GAL  : return "GAL";
   case BDS  : return "BDS";
     
   default : { cout << "*** warning: unknown time system!\n"; cout.flush(); }
  }
  return "TAI";
}
  
   
// convert tsys to string
// ----------
t_gtime::t_tsys t_gtime::str2tsys(string s)
{
  if( s.size() == 0 ){
    cout << "*** warning: not defined time system\n"; cout.flush();
    return TAI;
  }

  transform(s.begin(), s.end(), s.begin(), ::toupper);
  
  if(      s == "USER" ) return USER;
  else if( s == "TAI"  ) return TAI;
  else if( s == "UTC"  ) return UTC;
  else if( s == "LOC"  ) return LOC;
  else if( s == "GPS"  ) return GPS;
  else if( s == "GLO"  ) return GLO;
  else if( s == "GAL"  ) return GAL;
  else if( s == "BDS"  ) return BDS;
  else{ cout << "*** warning: not defined correct time system [" << s[0] << "]\n"; cout.flush(); }

  return TAI;
}


// set current time // acccuracy of 1 sec!
// ----------
t_gtime t_gtime::current_time(const t_tsys& ts)
{
  t_gtime tmp(UTC);
  tmp.from_time(time(NULL), 0.0, true);
  tmp.tsys( ts );
  return tmp;
}


// set from time
// -----------
int t_gtime::from_time(const time_t& tt, const double& ds, const bool& conv )
{
  struct tm *tm = gmtime( &tt );
  from_ymd( tm->tm_year+1900,  tm->tm_mon+1,   tm->tm_mday,
	          tm->tm_hour*3600 + tm->tm_min*60 + tm->tm_sec,
	          ds, conv );

//  if(conv) _to_tai(); --> ALREADY CONVERTED
  _norm_dsec();
  _norm_sod();
  return 0;
}


// set from GPS week and second of week
// -----------
int t_gtime::from_gws(const int& gw, const double& sow, const bool& conv)
{ 
  int dw = (int)floor(sow/86400);
  _mjd = 44244 + 7*gw + dw;
  _sod = (int)floor(sow - dw*86400);
  _dsec = sow - dw*86400 - _sod;

  if(conv) _to_tai(); // CONVERT
  _norm_dsec();
  _norm_sod();
  return 0;
}


// set from GPS week and dow of week
// -----------
int t_gtime::from_gwd(const int& gw, const int& dw, const int& sd, 
                      const double& ds, const bool& conv )
{
  _mjd = 44244 + 7*gw + dw;
  _sod = sd;
  _dsec = ds;

  if(conv) _to_tai(); // CONVERT
  _norm_dsec();
  _norm_sod();
  return 0;
}


// set from year, month, day, second of day, fractional seconds
// -----------
int t_gtime::from_ymd(const int& yr, const int& mn, const int& dd, 
                      const int& sd, const double& ds, const bool& conv )
{
  int year(yr);
  _norm_year(year);
  _mjd  = _ymd_mjd(year,mn,dd);
  _sod  = sd;
  _dsec = ds;

  if(conv) _to_tai(); // CONVERT
  _norm_dsec();
  _norm_sod();
  return 0;
}


// set from year, month, day, hour, minute, float seconds
// -----------
int t_gtime::from_ymdhms(const int& yr, const int& mn, 
                         const int& dd, const int& h, 
                         const int& m, const double& s, 
                         const bool& conv )
{
  int year(yr);
  _norm_year(year);
  _mjd  = _ymd_mjd(year,mn,dd);
   
  _sod  = int((h*3600.0 + m*60.0) + floor(s));
  _dsec = s - 1.0*(floor(s));

  if(conv) _to_tai(); // CONVERT
  _norm_dsec();
  _norm_sod();
  return 0;
}


// set from MJD
// -----------
int t_gtime::from_mjd(const int& mjd, const int& sd, 
                      const double& ds, const bool& conv)
{
  _mjd  = mjd;
  _sod  = sd;
  _dsec = ds;

  if(conv) _to_tai(); // CONVERT
  _norm_dsec();
  _norm_sod();
  return 0;
}


// from string
// ----------
int t_gtime::from_str(const string& ifmt, const string& idat, const bool& conv )
{
  int cYMD = _ymd;
  int cHMS = _hms;
  int y,b,d,h,m,s;   y = b = d = h = m = s = 0;
  size_t idx;
  size_t fmtadd = 0;
  size_t datadd = 0;
  size_t fmtpos = 0;
  size_t datpos = 0;

  // search and process keys from left to right
  while( (fmtpos < ifmt.length()) && ( idx = ifmt.find('%',fmtpos) ) != string::npos ){
    
    fmtpos++; // a priori encrease
     
    for( int i=0; i < MAX_DT; ++i ){
      string tmp = TD[i];
 
      // too short ifmt format string, skip!
      if( ifmt.length() - idx - 1 < tmp.length() ) continue;

      // dat string not identified, skipped !
      if( ifmt.substr(idx+1,tmp.length()).compare( tmp ) != 0 ) continue;

      fmtpos  = idx + tmp.length() + 1; // end of format reading
      datpos  = idx - fmtadd + datadd;  // sum of all idat characters
      fmtadd +=       tmp.length() + 1; // sum of all ifmt characters

      if(      !tmp.compare("Y")){ y =   str2int(idat.substr(datpos,4));  datadd+=4; cYMD+=t_tymd(_year);
      }else if(!tmp.compare("y")){ y =   str2int(idat.substr(datpos,2));  datadd+=2; cYMD+=t_tymd(_year);
      }else if(!tmp.compare("b")){ b =       mon(idat.substr(datpos,3));  datadd+=3; cYMD+=t_tymd(_mon);
      }else if(!tmp.compare("m")){ b =   str2int(idat.substr(datpos,2));  datadd+=2; cYMD+=t_tymd(_mon);
      }else if(!tmp.compare("d")){ d =   str2int(idat.substr(datpos,2));  datadd+=2; cYMD+=t_tymd(_day);
      }else if(!tmp.compare("j")){ b = 1;                                            cYMD+=t_tymd(_mon);
	                                 d =   str2int(idat.substr(datpos,3));  datadd+=3; cYMD+=t_tymd(_day);
      }else if(!tmp.compare("H")){ h =   str2int(idat.substr(datpos,2));  datadd+=2; cHMS+=t_thms(_hour);
      }else if(!tmp.compare("M")){ m =   str2int(idat.substr(datpos,2));  datadd+=2; cHMS+=t_thms(_min);
      }else if(!tmp.compare("S")){ s =   str2int(idat.substr(datpos,2));  datadd+=2; cHMS+=t_thms(_sec);
      }else if(!tmp.compare("W")){ y = 1980;                                         cYMD+=t_tymd(_year);
                                   b = 1;                                            cYMD+=t_tymd(_mon);
  	                               d+= 6+str2int(idat.substr(datpos,4))*7;datadd+=4;
      }else if(!tmp.compare("w")){ d+=   str2int(idat.substr(datpos,1));  datadd+=1; cYMD+=t_tymd(_day);
      }else if(!tmp.compare("v")){ s+=   str2int(idat.substr(datpos,6));  datadd+=6; cYMD+=t_tymd(_day);
	                                                                                   cHMS+=t_thms(_hour)
	                                                                                        +t_thms(_min)
                                                                                          +t_thms(_sec);
      }else if(!tmp.compare("s")){ s =   str2int(idat.substr(datpos,5));  datadd+=5; cHMS+=t_thms(_hour)
	                                                                                        +t_thms(_min)
	                                                                                        +t_thms(_sec);
      }else if(!tmp.compare("I")){ y = 2000;                                         cYMD+=t_tymd(_year);
	                             b = 1;                                                cYMD+=t_tymd(_mon);
	                             d = 1-51544+(str2int(idat.substr(datpos,11)));
	                                                                                   cYMD+=t_tymd(_day);
	                             s = 86400*(  str2int(idat.substr(datpos,11))
					                                -(str2int(idat.substr(datpos,11))) );
	                                                                        datadd+=11;cHMS+=t_thms(_hour)
	                                                                                        +t_thms(_min) 
	                                                                                        +t_thms(_sec);
      }else if(!tmp.compare("J")){ y = 2000;                                         cYMD+=t_tymd(_year);
 	                                 b = 1;                                            cYMD+=t_tymd(_mon);
                                   d = 1-51544 + str2int(idat.substr(datpos,5));
	                                                                      datadd+=5; cYMD+=t_tymd(_day);
      }else if(!tmp.compare("T")){ cerr << "input time-system conversion not yet supported here!\n";

      }else{ cerr << " warning : t_gtime - unknown date/time identifier [" << tmp << "]\n";
      }
    }
  }

  if( cYMD == ( _year | _mon | _day ) && cHMS == ( _hour | _min | _sec ) ){
     
    from_ymd( y, b, d, h*3600+m*60+s, 0.0, conv );
     
#ifdef DEBUG
  cout << " from_str = " << y
                  << " " << b
                  << " " << d
                  << " " << h
                  << " " << m 
                  << " " << s
                  << " " << _mjd
                  << " " << _sod
                  << "\n";
#endif
     
    return 0;
  }

  if( cYMD == ( _year | _mon | _day ) ){

#ifdef DEBUG
  cout << " from_str = " << y
                  << " " << b
                  << " " << d
                  << " " << h
                  << " " << m 
                  << " " << s
                  << "\n";
#endif

    from_ymd( y, b, d, 0, 0.0, conv );
    return 0;
  }

  return 1;
}

// reset_dsec - reset dsec only!
// ---------- 
int t_gtime::reset_dsec(){
  _norm_dsec();
  _dsec = 0.0; 
  return 0;
}


// reset_sod - reset sod + dsec only!
// ----------
int t_gtime::reset_sod(){
  _norm_sod();
  _tai_to();
  _sod = 0;
  _to_tai();
  reset_dsec();
  return 0;
}


// time zone difference (LOC - UTC) [sec]
// ----------
int t_gtime::tzdiff() const
{
  switch( _tsys ){
    case USER  : return 0;
    case TAI   : return 0;
    case UTC   : return 0;
    case GPS   : return 0;
    case GLO   : return 3*3600;
    case GAL   : return 0;
    case BDS   : return 0;
    case LOC   : { time_t local, utc;
                   local =   time(NULL);
                   utc   = mktime(gmtime(&local));
                   return int(local - utc);
                 }
     default    : cerr << " warning : t_tsys not recognized !\n";
  }
  return 0;
}


// day saving time (summer - winter time, e.g. CEST - CET) [sec]
// ----------
int t_gtime::dstime() const
{  
  switch( _tsys ){
    case USER  : return 0;
    case TAI   : return 0;
    case UTC   : return 0;
    case GPS   : return 0;
    case GLO   : return 0;
    case BDS   : return 0;
    case GAL   : return 0;
    case LOC   : { 
                   struct tm tm;
                   ymd(tm.tm_year, tm.tm_mon, tm.tm_mday, false); // TAI !
                   hms(tm.tm_hour, tm.tm_min, tm.tm_sec,  false); // TAI !

                   tm.tm_year  -=  1900;
                   tm.tm_mon   -=  1;
                   tm.tm_sec   -=  leapsec() + tzdiff();  // convert tai to LOC !
                   tm.tm_isdst  = -1;

                   mktime(&tm);
                   return int(3600.0 * tm.tm_isdst); // [sec]
                 }
    default    : cerr << " warning : t_tsys not recognized !\n";
  }
  return 0;
}


// get TAI-UTC leap seconds [sec]
// ----------
int t_gtime::leapsec() const
{  
  int leap = 0;
  for( int i = sizeof(leapseconds)/sizeof(*leapseconds) -1; i >= 0; --i ){
    double mjd_time = _mjd + _sod/86400.0 - leapseconds[i]._leap/86400.0;
    double mjd_leap = static_cast<double>( _ymd_mjd( leapseconds[i]._year,
			                             leapseconds[i]._mon,
		 	                             leapseconds[i]._day  ));
    if( mjd_leap <= mjd_time ){
      leap = leapseconds[i]._leap;
      break;
    }
  } 
  return leap;
}


// getleap year
// ----------
int t_gtime::leapyear(const int& y) const
{
  if( y%4   == 0.0 ) return 1; //     leap year !
  if( y%400 == 0.0 ) return 1; //     leap year !
  if( y%100 == 0.0 ) return 0; // not leap year
   
  return 0;
}

   
   
// get modified julian date
// b = true (convert to TSYS)
// ----------
int t_gtime::mjd(const bool& conv) const
{
  // do not convert
  if( !conv || _tsys == USER || _tsys == TAI) return _mjd;

  // do convert
  int mjd = static_cast<int>(_mjd + (_sod+tai_tsys(_tsys))/86400.0); // 86400.0  have to be float!
  return mjd;
}


// get double modified jujian date
// -----------------------------------
double t_gtime::dmjd(const bool& conv) const
{
  double ret;
  ret = (double)mjd(conv) +
        (double)sod(conv)/24.0/60.0/60.0 + dsec(conv)/24.0/60.0/60.0;
   
  return ret;
}


// get seconds of day
// b = true (convert to TSYS)
// ----------
int t_gtime::sod(const bool& conv) const
{
  // do not convert
  if( !conv || _tsys == USER || _tsys == TAI ) return _sod;

  // do convert
  int sod = _sod + tai_tsys(_tsys);
  
  while( sod >= 86400 ){ sod -= 86400; }
  while( sod <  0     ){ sod += 86400; }

  return sod;
}


// get dsec
// b = true (convert to TSYS)
// ----------
double t_gtime::dsec(const bool& conv) const
{
  return _dsec;
}


// get time_t
// ----------
time_t t_gtime::tim() const
{
  int y,b,d,h,m,s;
  ymd(y,b,d,true);  // convert to UTC ?
  hms(h,m,s,true);  // convert to TAI UTC ?

  struct tm t;
  t.tm_year  = y - 1900;
  t.tm_mon   = b - 1;
  t.tm_mday  = d;
  t.tm_hour  = h;
  t.tm_min   = m;
  t.tm_sec   = s;
  t.tm_isdst = -1;

  return mktime(&t);
}


// get day of year
// ----------
int t_gtime::doy() const
{
  int y, m, d;
  ymd(y,m,d,true);

  int doy = d;
  for( int i=0; i<m && i<13; ++i)  doy += monthdays[i];
  if( m > 2 )                      doy += leapyear(y);   // count for leap year
  return doy;
}


// get BDS week
// ----------
int t_gtime::bwk() const
{
  return static_cast<int>((mjd(true)-44244.0)/7.0) - CONV_BWK2GWK;  // return true TS
}

// get GPS week
// ----------
int t_gtime::gwk() const
{
  return static_cast<int>((mjd(true)-44244.0)/7.0);  // return true TS
}


// get day of GPS(BDS) week
// ----------
int t_gtime::dow() const
{
  return static_cast<int>( mjd(true) - 44244.0 - gwk()*7 );  // no problem gwk for BDS here
}


// get seconds of GPS(BDS) week
// ----------
int t_gtime::sow() const
{
  return static_cast<int>( dow()*86400.0 + sod() ); // no problem gwk for BDS here
}


// get year (4-char)
// ----------
int t_gtime::year() const
{  
  int y, m, d;
  ymd( y, m, d, true);

  return static_cast<int>( y );
}

// get year (2-char) (1980-2079 only)
// ----------
int t_gtime::yr() const
{  
  int y, m, d;
  ymd( y, m, d, true);
  int yr = this->yr(y);
  return yr;
}


// get month
// ----------
int t_gtime::mon() const
{  
  int y, m, d;
  ymd( y, m, d, true);

  return static_cast<int>( m );
}


// get day
// ----------
int t_gtime::day() const
{  
  int y, m, d;
  ymd( y, m, d, true);

  return static_cast<int>( d );
}


// get hour
// ----------
int t_gtime::hour() const
{ 
  return static_cast<int>( (sod(true)%86400)/3600.0 );
}


// get minutes
// ----------
int t_gtime::mins() const
{ 
  return static_cast<int>( (sod(true)%3600)/60.0 );
}


// get seconds
// ----------
int t_gtime::secs() const
{ 
  return static_cast<int>( sod(true)%60 );
}


// get hour
// ----------
void t_gtime::hms( int& h, int& m, int& s, bool conv ) const
{
  h = static_cast<int>( ( sod(conv)%86400 )/3600.0 );
  m = static_cast<int>( ( sod(conv)%3600  )/60.0 );
  s = static_cast<int>( ( sod(conv)%60    ) );
}

   
// MJD -> YMD
// ----------
void t_gtime::ymd( int& y, int& m, int& d, bool conv ) const
{
  int    jj,mm,dd;
  long   ih, ih1, ih2 ;
  double t1, t2,  t3, t4;
  double mjd = 0.0;

  mjd = this->mjd( conv ) + (this->sod( conv ))/86400.0;

  //  DO NOT USE THIS DUE TO UNDERFLOW
  //mjd = this->mjd( conv ) + (this->sod( conv ) + this->dsec( conv ))/86400; 
   
  t1  =  1.0 + mjd - fmod( mjd, 1.0 ) + 2400000.0;
  t4  =  fmod( mjd, 1.0 );
  ih  =  int( (t1 - 1867216.25)/36524.25 );
  t2  =  t1 + 1 + ih - ih/4;
  t3  =  t2 - 1720995.0;
  ih1 =  int( (t3 - 122.1)/365.25 );
  t1  =  ih1*365.25 - fmod( ih1*365.25, 1.0 );
  ih2 =  int( (t3 - t1)/30.6001 );
  dd  =  int(t3 - t1 - int( ih2*30.6001 ) + t4);
  mm  =  ih2 - 1;
  if( ih2 > 13 ) mm = ih2 - 13;
  jj  = ih1;
  if( mm <= 2 ) jj = jj + 1;
  y=jj;
  m=mm;
  d=dd;
}


// to string
// ----------
string t_gtime::str(const string& ofmt, const bool& conv ) const
{
  t_gtime gt(*this);

  char cstr[12];
  string str = ofmt; // copy of requested format
  size_t idx =  0;
  int y,b,d,h,m,s;
  gt.ymd( y, b, d, conv );
  gt.hms( h, m, s, conv );
  int y2 = gt.yr(y);

  // search and process all % identificators
  while( ( idx = str.find('%') ) != string::npos && idx+1 <= str.length() ){
   
    bool replace = false;
    for( int i=0; i < MAX_DT; ++i ){
      string tmp = TD[i];

      // dat string not identified, skipped !
      if( str.substr(idx+1,tmp.length()).compare( tmp ) != 0 ) continue;

      if(      !tmp.compare("Y")){ sprintf(cstr,"%04i",y);
      }else if(!tmp.compare("y")){ sprintf(cstr,"%02i",y2);
      }else if(!tmp.compare("b")){ sprintf(cstr,"%3.3s", gt.mon(b).c_str());
      }else if(!tmp.compare("m")){ sprintf(cstr,"%02i",b);
      }else if(!tmp.compare("d")){ sprintf(cstr,"%02i",d);
      }else if(!tmp.compare("j")){ sprintf(cstr,"%03i",gt.doy());
      }else if(!tmp.compare("H")){ sprintf(cstr,"%02i",h);
      }else if(!tmp.compare("M")){ sprintf(cstr,"%02i",m);
      }else if(!tmp.compare("S")){ sprintf(cstr,"%02i",s);
      }else if(!tmp.compare("W")){ sprintf(cstr,"%04i",gt.gwk());
      }else if(!tmp.compare("w")){ sprintf(cstr,"%01i",gt.dow());
      }else if(!tmp.compare("v")){ sprintf(cstr,"%6i", gt.sow());
      }else if(!tmp.compare("s")){ sprintf(cstr,"%5i", gt.sod( conv ));
      }else if(!tmp.compare("J")){ sprintf(cstr,"%5i", gt.mjd( conv ));
      }else if(!tmp.compare("I")){ sprintf(cstr,"%11.5f", gt.mjd( conv )
					                + gt.sod( conv )/86400.0
					                + gt.dsec( conv )/86400.0);
      }else if(!tmp.compare("T")){ sprintf(cstr,"%3s", sys().c_str() );
      }else{ cerr << " warning : t_gtime - unknown date/time identifier [" << tmp << "]\n";
      }

      str.replace(idx,2,cstr);
      replace = true;
    }
    // replace unknown % occurance
    if( !replace ) str.replace(idx,1,"");
  }
  return str;
}


// to string
// ----------
string t_gtime::str_ymd(const string& str, const bool& conv) const
{
  char cstr[12];
  int y=0, b=0, d=0;
  this->ymd( y, b, d, conv );
  sprintf(cstr,"%04i-%02i-%02i",y,b,d);
  return str + string(cstr);
}


// to string
// ----------
string t_gtime::str_hms(const string& str, const bool& conv) const
{
  char cstr[12];
  int h=0, m=0, s=0;
  this->hms( h, m, s, conv );
  sprintf(cstr,"%02i:%02i:%02i",h,m,s);
  return str + string(cstr);
}


// to string
// ----------
string t_gtime::str_ymdhms(const string& str, 
			   const bool& ts,
			   const bool& conv) const
{
  char cstr[25];
  int y=0, b=0, d=0, h=0, m=0, s=0;
  this->ymd( y, b, d, conv );
  this->hms( h, m, s, conv );
  if( ts ) sprintf(cstr,"%04i-%02i-%02i %02i:%02i:%02i[%3s]",y,b,d,h,m,s, sys().c_str());
  else     sprintf(cstr,"%04i-%02i-%02i %02i:%02i:%02i",y,b,d,h,m,s);
  return str + string(cstr);
}


// get yr (1980-2079 only)
// ----------
int t_gtime::yr(const int& y) const
{  
  if( y <= 2079 && y >= 2000 ) return y-2000;
  if( y <  2000 && y >= 1980 ) return y-1900;
   
  return -1;
}


// get month 3-char string
// ----------
string t_gtime::mon(const int& m) const
{ 
  if( m < 1 || m > 12 ) return "XXX";
  return MON[m];
}


// get month number
// ----------
int t_gtime::mon(const string& str) const
{ 
  for(int i = 0; i<12; ++i){
     if( str.compare(MON[i]) == 0 )
       return i+1; 
  }
  return 0;
}


// get system identifier
// ----------
string t_gtime::sys() const
{ 
  string ts( t_tstr[_tsys] );
  return ts;
}


// operator reduce [sec] [TAI]
// ----------
double t_gtime::operator-(const t_gtime& t) const
{
  return ( this->diff(t) );
}


// operator minus [TAI]
// ----------
t_gtime t_gtime::operator-(const double& sec) const
{
  t_gtime tmp(*this);
  tmp.add_secs( - static_cast<int>(sec));
  tmp.add_dsec( - static_cast<double>( sec - static_cast<int>(sec) ));  // already normed
  return tmp;
}


// operator add [TAI]
// ----------
t_gtime t_gtime::operator+(const double& sec) const
{
  t_gtime tmp(*this);
  tmp.add_secs( + static_cast<int>(sec));
  tmp.add_dsec( + static_cast<double>( sec - static_cast<int>(sec) ));  // already normed
  return tmp;
}


// operator comp [TAI]
// ----------
bool t_gtime::operator<(const t_gtime& t) const
{
  return (    ( mjd(true)  < t.mjd(true) )
           || ( mjd(true) == t.mjd(true) && sod(true)  <  t.sod(true) )
	   || ( mjd(true) == t.mjd(true) && sod(true) ==  t.sod(true) && dsec(true) < t.dsec(true) )
	 );
}


// operator comp [TAI]
// ----------
bool t_gtime::operator>(const t_gtime& t) const
{
  return (    ( mjd(true)  > t.mjd(true) )
           || ( mjd(true) == t.mjd(true) && sod(true)  >  t.sod(true) )
	   || ( mjd(true) == t.mjd(true) && sod(true) ==  t.sod(true) && dsec(true) > t.dsec(true) )
	 );
}


// operator equiv [TAI]
// ----------
bool t_gtime::operator==(const t_gtime& t) const	        
{
    return ( mjd(true) == t.mjd(true) && sod(true) == t.sod(true) && dsec(true) == t.dsec(true) );
}

// operator not equiv [TAI]
// ----------
bool t_gtime::operator!=(const t_gtime& t) const	        
{
    return ( mjd(true) != t.mjd(true) || sod(true) != t.sod(true) || dsec(true) != t.dsec(true) );
}

// operator =
// ----------
t_gtime t_gtime::operator=(const t_gtime& t)
{
  _mjd  = t.mjd(false);
  _sod  = t.sod(false);
  _dsec = t.dsec(false);
  tsys( t.tsys() );

  return *this;
}

// add seconds
// ----------
void t_gtime::add_secs(const int& sec)
{
  _sod += sec;
  _norm_sod();
}


// add dseconds
// ----------
void t_gtime::add_dsec(const double& dsec)
{
  _dsec += dsec;
  _norm_dsec();
  _norm_sod();
}


// OBSOLETE, because via operator - exists
// time difference (this - t) [s]
// ----------
double t_gtime::diff(const t_gtime& t) const
{   
  bool b = false; // compare in TAI
   
  return (    mjd(b)*86400.0 +   sod(b)  +  dsec(b)
	        - t.mjd(b)*86400.0 - t.sod(b) - t.dsec(b) );
}


// -----------------------------------------------------------------------------------
// private functions
// -----------------------------------------------------------------------------------

// convert to TAI from requested time
// ----------
void t_gtime::_to_tai()
{
  _norm_dsec();
  _norm_sod();
  switch(_tsys){
   case USER  : break;;
   case TAI   : break;;
   case UTC   : _sod += leapsec(); break;;
   case GPS   : _sod += TAI_GPS; break;;
   case GLO   : _sod -= tzdiff();
                _sod += leapsec(); break;;
   case GAL   : _sod += TAI_GAL; break;;
   case BDS   : _sod += TAI_BDS; break;;
   case LOC   : _sod -= tzdiff();
                _sod -= dstime(); 
                _sod += leapsec(); break;;
   default    : break;;
  }
  _norm_dsec();
  _norm_sod();
}


// convert to tsys and return
// -----------
t_gtime t_gtime::_tai_to()
{
  _sod += tai_tsys(_tsys);
  _norm_dsec();
  _norm_sod();

  return *this;
}

// get TAI-tsys difference [sec]
// ----------
int t_gtime::tai_tsys(const t_tsys& ts) const
{
  double sec = 0.0;
  switch(ts){
   case USER  : break;;
   case TAI   : break;;
   case UTC   : sec -= leapsec(); break;;
   case GPS   : sec -= TAI_GPS; break;;
   case GLO   : sec -= leapsec();
                sec += tzdiff(); break;;
   case GAL   : sec -= TAI_GAL; break;;
   case BDS   : sec -= TAI_BDS; break;;
   case LOC   : sec -= leapsec();
                sec += dstime();
                sec += tzdiff(); break;;
   default    : break;;
  }

  return sec;
}

// normalize dsec
// ----------
void t_gtime::_norm_dsec()
{
  while( _dsec >= 1.0 ){ _dsec -= 1.0; _sod += 1; }
  while( _dsec <  0.0 ){ _dsec += 1.0; _sod -= 1; }
}


// normalize sod
// ----------
void t_gtime::_norm_sod()
{
  while( _sod >= 86400 ){ _sod -= 86400; _mjd += 1; }
  while( _sod <  0     ){ _sod += 86400; _mjd -= 1; }
}


// normalize year
// ----------
void t_gtime::_norm_year(int& year) const
{
  if( year < 100.0 ) year += year < 80 ? 2000 : 1900;
}


// Year,Mon,Day -> MJD
// ----------
int t_gtime::_ymd_mjd(const int& yr, const int& mn, const int& dd) const
{
  int year(yr);
  int mon(mn);
  _norm_year(year);

  double mjd = 0.0;
  if( mon <= 2 ){ mon = mon + 12; year = year - 1; }
   
  mjd  = 365.25*year - fmod(365.25*year, 1.0) - 679006.0;
  mjd += floor(30.6001*(mon + 1)) + 2.0 - floor(year / 100) + floor(year / 400) + dd;
  return int(mjd);
}

// multiplatform ms sleep function
void t_gtime::gmsleep(unsigned int ms)         // milisecond
{

#if    defined __linux__ 
   usleep(ms*1000);
#elif  defined __APPLE__
   usleep(ms*1000);
#elif  defined _WIN32 || defined _WIN64
   Sleep(ms);
#else
   clock_t goal = ms + clock();
   while (goal > clock());
#endif 
}


// multiplatform us sleep function
void t_gtime::gusleep(unsigned int us)
{
// this_thread::sleep_for(chrono::microseconds(us));
#if    defined __linux__ 
   usleep(us);
#elif  defined __APPLE__
   usleep(us);
#elif  defined _WIN32 || defined _WIN64
   this_thread::sleep_for(chrono::microseconds(us)); // c++11 feature
//#else
//   clock_t goal = us + clock();
//   while (goal > clock());
#endif

}

} // namespace
