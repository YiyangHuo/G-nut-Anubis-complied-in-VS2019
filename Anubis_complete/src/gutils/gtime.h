
#ifndef GTIME_H
#define GTIME_H 

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements gtime class (day and precise time)
  Version: $ Rev: $

  2011-01-10 /JD: created
  2018-07-13 /JD: updated (fixed current_time, revised FIRST/LAST time, speeded etc.)

 * Internal representation is in TAI time-system (TS)
 * If input is in other TS, it needs to be set via constructor
 * and all inputs are then converted from TS to TAI.
 * All outputs are implicitely in TAI thus use of
 * gps(), utc(), local() duplicating function is required.
 * Input TS can be changed only via convert(TS) function.
-*/

#include <time.h>
#include <string>
#include <string.h>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>  // Win/Lin mutex
#endif

#include "gutils/gmutex.h"

using namespace std;

namespace gnut {
   
#if   defined GTIME_TAI
#define                 DEFAULT_TIME TAI   // default gtime TAI
#elif defined GTIME_UTC
#define                 DEFAULT_TIME UTC   // default gtime UTC
#elif defined GTIME_USER
#define                 DEFAULT_TIME USER  // default gtime USER
#elif defined GTIME_LOC
#define                 DEFAULT_TIME LOC   // default gtime LOCAL
#elif defined GTIME_GPS
#define                 DEFAULT_TIME GPS   // default gtime GPS
#else
#define                 DEFAULT_TIME GPS   // default gtime GPS
#endif

#define MAX_DT  16  // number of time string identifiers
#define MAXLEAP 50  // maximum number of leap second in table
#define TAI_GPS 19  // seconds of which GPS is ahead of TAI at 6.1.1980
#define TAI_BDS 33  // seconds of which BDS is ahead of TAI at 1.1.2006
#define TAI_GAL 19  // coincide with GPS time

#define CONV_BWK2GWK  1356 // convert BDS to GPS week

const static string MON[12] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

static const int monthdays[13] = { 0,31,28,31,30,31,30,31,31,30,31,30,31 };

struct leapsec_table
{
  int    _year;
  int     _mon;
  int     _day;
  int    _leap;
};

const static struct leapsec_table leapseconds[] = {
  { 1971, 12, 31, 11 },
  { 1972, 12, 31, 12 },
  { 1973, 12, 31, 13 },
  { 1974, 12, 31, 14 },
  { 1975, 12, 31, 15 },
  { 1976, 12, 31, 16 },
  { 1977, 12, 31, 17 },
  { 1978, 12, 31, 18 },
  { 1979, 12, 31, 19 },
  { 1981, 06, 30, 20 },
  { 1982, 06, 30, 21 },
  { 1983, 06, 30, 22 },
  { 1985, 06, 30, 23 },
  { 1987, 12, 31, 24 },
  { 1989, 12, 31, 25 },
  { 1990, 12, 31, 26 },
  { 1992, 06, 30, 27 },
  { 1993, 06, 30, 28 },
  { 1994, 06, 30, 29 },
  { 1995, 12, 31, 30 },
  { 1997, 06, 30, 31 },
  { 1998, 12, 31, 32 },
  { 2005, 12, 31, 33 },
  { 2008, 12, 31, 34 },
  { 2011, 06, 30, 35 },
  { 2015, 06, 30, 36 },
  { 2017, 01, 01, 37 }
};

// supported io keys of fixed size !
const static string TD[MAX_DT] = { 
   "Y",    // [4] year            (1900..2099)                             
   "y",	   // [2] 2-dig year      (80..99,00..79)                          
   "m",	   // [2] month           (1..12+)                                
   "b",	   // [3] month           (Jan..Dec)                                
   "d",	   // [2] day of month    (1..31+)                                 
   "j",	   // [3] day of year     (0..365+)                                
   "H",	   // [2] hours           (0..23+)                                 
   "M",	   // [2] minutes         (0..59+)                                 
   "S",	   // [2] seconds         (0..59+)                                 
   "W",	   // [4] GPS week        (0000..XXXX)                             
   "w",	   // [1] day of week     (0..6)                                   
   "v",	   // [6] seconds of week (0..604800+)                             
   "s",	   // [5] seconds of day  (0..86400+)                              
   "J",    // [5]  modified julian date = integer only ! (e.g. 55725)      
   "I",	   // [11] modified julian date                  (e.g. 55725.81712)
   "T" 	   // [3] time-system string
};

// conversion to GAL + GLO are not implemented
const static char t_tstr[8][4] = { "USR", "TAI", "UTC", "LOC", "GPS", "GAL", "GLO", "BDS" };

 //static t_gmutex  _gmutex;
#ifdef BMUTEX
 static boost::mutex   _static_mutex;
#endif
   
class t_gtime {

 public:

   // conversion to GAL + GLO are not implemented
   enum t_tsys { USER, TAI, UTC, LOC, GPS, GAL, GLO, BDS };
   static string tsys2str(t_tsys ts);
   static t_tsys str2tsys(string s);

   t_gtime(const t_tsys& ts=DEFAULT_TIME); // initiate with current time !
   
   t_gtime(const time_t& tt, 
           const double& dsec=0,
           const t_tsys& ts=DEFAULT_TIME);
   
   t_gtime(const int& yr, const int& mn,
           const int& dd, const int& hr,
           const int& mi, const int& sc,
           const double& ds=0.0, const t_tsys& ts=DEFAULT_TIME);
   
   t_gtime(const int& gw, const int& dw,
           const int& sd, const double& ds=0.0, 
           const t_tsys& ts=DEFAULT_TIME);
   
   t_gtime(const int& gw, const double& sow,
           const t_tsys& ts=DEFAULT_TIME);
   
   t_gtime(const int& mjd, const int& sd,
           const double& ds=0.0, const t_tsys& ts=DEFAULT_TIME);
   
  ~t_gtime( );

   // returning only static functions
   static t_gtime current_time(const t_tsys& ts=DEFAULT_TIME);   // get current time

   // set functions
   int from_time(const time_t& tt,
                 const double& ds=0.0, const bool& conv=true);   // set from time_t,dsec
   
   int from_gws(const int& gw,
                const double& sow, const bool& conv=true);       // set from GPSwk,sow,dsec
   
   int from_gwd(const int& gw, const int& dw, const int& sd,
                const double& ds=0.0, const bool& conv=true);    // set from GPSwk+dw,sod,dsec
   
   int from_ymd(const int& yr, const int& mn, const int& dd,
                const int& sd, const double& ds=0.0,
                const bool& conv=true);                          // set from yr+mn+dd,sod,dsec

   int from_ymdhms(const int& yr, const int& mn, const int& dd,
                   const int& h, const int& m, const double& s,
                   const bool& conv=true);                       // set from yr+mn+dd + H+M+S
   
   int from_mjd(const int& mjd, const int& sd,
                const double& ds=0.0, const bool& conv=true);    // set from mjd,sod,dsec
   
   int from_str(const string& ifmt, const string& dat, 
                const bool& conv=true);                          // set from defined string
   
   int reset_dsec();                                             // reset dsec only (in TAI)!
   int reset_sod();                                              // reset sod + dsec only (in TAI)!

   int tsys(const t_tsys& ts){ _tsys = ts; return 0; }           // set new io time-system

   // get functions (default includes conversion to _tsys, if conv=false return TAI)
   int     mjd(const bool& conv=true) const;       // get MJD
   double dmjd(const bool& conv=true) const;       // get MJD
   int     sod(const bool& conv=true) const;       // get seconds of day
   double dsec(const bool& conv=true) const;       // get fractional sec
   				             	                                        
   // get functions (always includes conversion to _tsys, for TAI use only mjd, sod and dsec)
   time_t  tim() const;                            // get time_t
   int     gwk() const;                            // get GPS week
   int     bwk() const;                            // get BDS week
   int     dow() const;                            // get integer day of GPS week
   int     sow() const;                            // get integer Sec of Week
   int    year() const;                            // get integer 4-char year
   int      yr() const;                            // get integer 2-char year (1980-2079 only)
   int     doy() const;                            // get integer day of year
   int     day() const;                            // get integer day of month
   int     mon() const;                            // get integer day
   int    hour() const;                            // get integer hour
   int    mins() const;                            // get integer minute
   int    secs() const;                            // get integer seconds
   string  sys() const;                            // get time system identifier
   			                      			                      
   void   hms(int& h,int& m,int& s,bool conv=true) const; // get hour, minute, seconds
   void   ymd(int& y,int& m,int& d,bool conv=true) const; // get year, month, day
   string str(const string& ofmt="%Y-%m-%d %H:%M:%S",
              const bool& conv=true) const;               // get any string
   string str_ymd(const string& str="",
                  const bool& conv=true) const;           // get fixed date string
   string str_hms(const string& str="",
                  const bool& conv=true) const;           // get fixed time string
   string str_ymdhms(const string& str="",
		     const bool& ts=false,
                     const bool& conv=true) const;        // get fixed date and time string

   // general conversion functions
   int    yr(const int& y) const;                  // get integer 2-char year (1980-2079 only)
   string mon(const int& m) const;                 // get month string
   int    mon(const string& str) const;            // get month number
   
   int tzdiff() const;                             // get UTC-LOC difference [sec]
   int dstime() const;                             // get day-saving time [sec]
   int leapsec() const;                            // get TAI-UTC leap seconds [sec]
   int leapyear(const int& y) const;               // get YEAR with leapsec [0,1]
   int tai_tsys(const t_tsys& ts) const;           // get TAI-tsys difference [sec]
   t_tsys tsys() const {return _tsys; }            // get required inp/out time-system

   void   add_secs(const int& sec);                // date/time with added X-seconds
   void   add_dsec(const double& dsec);            // date/time with added X-dseconds
   double diff(const t_gtime& t) const;            // time difference (this - t) [s]

   bool    operator<(const t_gtime& t) const;
   bool    operator>(const t_gtime& t) const;
   bool    operator==(const t_gtime& t) const;
   bool    operator!=(const t_gtime& t) const;   
   double  operator-(const t_gtime& t) const;      // [sec]
   t_gtime operator-(const double& sec) const;     // gtime
   t_gtime operator+(const double& sec) const;     // gtime
   t_gtime operator=(const t_gtime& t);

   static void gmsleep(unsigned int ms);           // multiplatform msleep function [ms]
   static void gusleep(unsigned int us);           // multiplatform usleep function [us]
   
 protected:   

 private:
   int     _mjd;                                   // integer MJD [TAI]
   int     _sod;                                   // seconds of day [TAI]
   double  _dsec;                                  // fractional seconds
   t_tsys  _tsys;                                  // time system
   t_gtime _tai_to();                              // conversion from TAI to TS with changing gtime !
   void    _to_tai();                              // conversion from TS to TAI with changing gtime !

   void    _norm_dsec();                           // normalize dsec (range in 0-1 sec)
   void    _norm_sod();                            // normalize sod  (range in 0-86399 sec)
   void    _norm_year(int& year) const;            // normalize year (range 1980 - 2079)
   int     _ymd_mjd(const int& yr,
                    const int& mn, 
                    const int& dd ) const;
   
   enum t_tymd { _ymd = 0, _year = 1, _mon = 2, _day = 4 };   // internal counter
   enum t_thms { _hms = 0, _hour = 1, _min = 2, _sec = 4 };   // internal counter

};

const static t_gtime FIRST_TIME(44239,     0, 0.0, t_gtime::TAI);   // first valid time  "1980-01-01 00:00:00"
const static t_gtime  LAST_TIME(80763, 86399, 0.0, t_gtime::TAI);   // last  valid time  "2079-12-31 23:59:59"


} // namespace

#endif
