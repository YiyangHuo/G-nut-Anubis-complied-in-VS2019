
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
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
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>

#include "gdata/gnavsbs.h"
#include "gmodels/gephplan.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

/* --- */
t_gnavsbs::t_gnavsbs()
  : t_gnav(),
    _toc(t_gtime::GPS)
{
  gtrace("t_gnavsbs::t_gnavsbs");
   
  id_type(t_gdata::EPHSBS);
  id_group(t_gdata::GRP_EPHEM);
   
}

/* --- */
t_gnavsbs::~t_gnavsbs(){
  gtrace("t_gnavsbs::~t_gnavsbs");
}



// check message // NEED TO IMPROVE !
// ----------
int t_gnavsbs::chk(set<string>& msg)
{
  gtrace("t_gnavsbs::chk");
  _gmutex.lock();
	 
  if( ! _healthy() ){
    msg.insert("Issue: "+_sat+" nav unhealthy satellite " + _toc.str_ymdhms());
    _validity = false;
  }

  if( _tki >= 9.9999e8 ){
    msg.insert("Issue: "+_sat+" nav too large value [tki] " + dbl2str(_tki));
    _validity = false;
  }

  _gmutex.unlock();  return 0;
}


/* --- */
/*int t_gnavsbs::ura( double acc ) const
{
   int i; 
   for( i = 0; i < (int)sizeof(ura_eph); i++)
     if(acc>ura_eph[i]) break;
   
   return i;
}
*/

// local function
// ----------
int t_gnavsbs::pos( const t_gtime& t, double  xyz[], double var[], double vel[], bool chk_health )
{
  gtrace("t_gnavsbs::_pos");  

  if( sat().empty() ) return -1; // not valid !!!
  if( chk_health && _healthy() == false ){
     if( _log ) _log->comment(3,"gephsbs","not healthy sat " + sat() + " excluded from pos calculation " + t.str_ymdhms());
     return -1;  // HEALTH NOT OK   
  }
   
  xyz[0] = xyz[1] = xyz[2] = 0.0;
  if( var ) var[0] = var[1] = var[2] = 0.0;
  if( vel ) vel[0] = vel[1] = vel[2] = 0.0;   

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif

  _gmutex.lock();   

//  cout << "t = " << t.str_hms() << " toc = " << _toc.str_hms("", true) << endl;
  double Tk = t.diff(_toc);                    // T - toe difference
  
//  cout << "t: " << t.str_ymdhms() << " " << t.sys() << " toe: " << _toc.str_ymdhms() << " " << _toc.sys() << endl;
//  cout << "Tk: " << Tk << endl;
   
  if ( fabs(Tk) > MAX_SBS_TIMEDIFF  )
  {
     if (_log) _log->comment(2,"gnavsbs","CRD " + _sat + ": The user time and Glonass ephemerides epoch differ too much. TOE: " 
			                  + _toc.str_ymdhms() + " T: " + t.str_ymdhms());
     else cerr << "gnavsbs: CRD: The user time and Glonass ephemerides epoch differ too much" << endl;
     _gmutex.unlock(); return -1;
  }
   
   if (double_eq(_x,0.0) || double_eq(_y,0.0) || double_eq(_z,0.0)) {_gmutex.unlock(); return -1;}

   xyz[0] = _x + _x_d*Tk + _x_dd*Tk*Tk/2.0;
   xyz[1] = _y + _y_d*Tk + _y_dd*Tk*Tk/2.0;
   xyz[2] = _z + _z_d*Tk + _z_dd*Tk*Tk/2.0;   

//  cout << sat() << " " << t.str_hms() << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;     
   
   // velocity at positon t
   if (vel){
     vel[0] = _x_dd;
     vel[1] = _y_dd;
     vel[2] = _z_dd;
   }
   
  _gmutex.unlock();   
  return 0;
}


/* ---
   t is glo time of transmission, 
   i.e. glo time corrected for transit time (range/speed of light)

*/
int t_gnavsbs::clk( const t_gtime& t, double*    clk, double*   var, double*  dclk, bool chk_health )
{
  gtrace("t_gnavsbs::clk");   

  if( sat().empty() ) return -1; // not valid !!!

  if( chk_health && _healthy() == false ){
     if( _log ) _log->comment(3,"gephsbs","not healthy sat " + sat() + " excluded from clk calculation " + t.str_ymdhms());   
     return -1;  // HEALTH NOT OK
  }
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

//  cout << "gnavsbs: CLK " << _sat << ": The user time and SBAS ephemerides epoch differs too much. TOE: " << _toc.str_ymdhms() << " T: " << t.str_ymdhms() << endl;

  double Tk = t.diff(_toc);                    // T - toe difference
   
//cout << "Tk = " << Tk << " maxEphAge = " << MAX_SBS_TIMEDIFF << endl;
  if ( fabs(Tk) > MAX_SBS_TIMEDIFF )
  {
     if (_log) _log->comment(2,"gnavsbs","CLK " + _sat + ": The user time and SBAS ephemerides epoch differs too much. TOE: "
			                  + _toc.str_ymdhms() + " T: " + t.str_ymdhms());     
     else cerr << "gnavsbs: CLK: The user time and SBAS ephemerides epoch differs too much" << endl;
     _gmutex.unlock();
     return -1;          
  }

  for (int i=0;i<2;i++) {   	
      Tk -= _f0 + _f1*Tk;
  }

  *clk = _f0 + _f1*Tk;
   
  _gmutex.unlock(); return 0;
}

int t_gnavsbs::data2nav( string sat, const t_gtime& ep, const t_gnavdata& data ){

  gtrace("t_gnavsbs::data2nav");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  if( sat.substr(0,1) == "S" )  _sat =     sat;
  else                          _sat = "S"+sat;
   
   _epoch = ep;
   _toc   = ep;
   _f0    = data[0];  
   _f1    = data[1];
   _tki   = data[2];
    if (_tki < 0) _tki += 86400;

   _x     = data[3] * 1.e3;
   _x_d   = data[4] * 1.e3;
   _x_dd  = data[5] * 1.e3;
   
   _health = data[6]; 
     
   _y     = data[7] * 1.e3;
   _y_d   = data[8] * 1.e3;
   _y_dd  = data[9] * 1.e3;
   
   _C_rms = data[10];
   
   _z     = data[11] * 1.e3;
   _z_d   = data[12] * 1.e3;
   _z_dd  = data[13] * 1.e3;
   
   _iod   = data[14];      
   
   _gmutex.unlock(); return 0;
}


// convert gnav_sbs element to general gnavdata
// ----------
int t_gnavsbs::nav2data( t_gnavdata& data ){
  
  gtrace("t_gnavsbs::nav2data");
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

   data[0]  = _f0;    
   data[1]  = _f1;
   data[2]  = _tki;   //  if (_tki < 0) _tki += 86400;

   data[3]  = _x     / 1.e3;
   data[4]  = _x_d   / 1.e3;
   data[5]  = _x_dd  / 1.e3;

   data[6]  = _health;
     
   data[7]  = _y     / 1.e3;
   data[8]  = _y_d   / 1.e3;
   data[9]  = _y_dd  / 1.e3;
   
   data[10] = _C_rms;
   
   data[11] = _z     / 1.e3;
   data[12] = _z_d   / 1.e3;
   data[13] = _z_dd  / 1.e3;
   
   data[14] = _iod;

  _gmutex.unlock(); return 0;
}

   
// get parameter value
// ----------
t_timdbl t_gnavsbs::param( const NAVDATA& n )
{
  _gmutex.lock();

  t_timdbl tmp;
  switch (n){
   case NAV_X     : tmp = make_pair(_toc,_x      *1e0);  break; // meters
   case NAV_XD    : tmp = make_pair(_toc,_x_d    *1e0);  break; // meters
   case NAV_XDD   : tmp = make_pair(_toc,_x_dd   *1e0);  break; // meters
   case NAV_Y     : tmp = make_pair(_toc,_y      *1e0);  break; // meters
   case NAV_YD    : tmp = make_pair(_toc,_y_d    *1e0);  break; // meters
   case NAV_YDD   : tmp = make_pair(_toc,_y_dd   *1e0);  break; // meters
   case NAV_Z     : tmp = make_pair(_toc,_z      *1e0);  break; // meters
   case NAV_ZD    : tmp = make_pair(_toc,_z_d    *1e0);  break; // meters
   case NAV_ZDD   : tmp = make_pair(_toc,_z_dd   *1e0);  break; // meters
    
   case NAV_IOD   : tmp = make_pair(_toc, _iod   *1e0);  break; //
   case NAV_HEALTH: tmp = make_pair(_toc,_health *1e0);  break; //

   default : break;
  }

  _gmutex.unlock(); return tmp;
}


// set parameter value
// ----------
int t_gnavsbs::param( const NAVDATA& n, double val )
{
  _gmutex.lock();

  switch (n){       // SELECTED only, ! use the same MULTIPLICATOR as in param()

   case NAV_IOD   : _iod    = val/1.e0;  break;
   case NAV_HEALTH: _health = val/1.e0;  break;

   default : break;
  }

  _gmutex.unlock(); return 0;
}


// print function
// ----------
string t_gnavsbs::line() const
{
  gtrace("t_gnavsbs::line");   
   
  int w = 20;
  ostringstream tmp;

  tmp  << " " << setw(3) << sat()
       << " " << _toc.str("%Y-%m-%d %H:%M:%S")
       << scientific << setprecision(12)
       << setw(w) << _f0
       << setw(w) << _f1
       << setw(w) << _tki
       << setw(w) << _x
       << setw(w) << _x_d
       << setw(w) << _x_dd
       << setw(w) << _health
       << setw(w) << _y
       << setw(w) << _y_d
       << setw(w) << _y_dd
       << setw(w) << _C_rms
       << setw(w) << _z
       << setw(w) << _z_d
       << setw(w) << _z_dd
       << setw(w) << _iod;     

  return tmp.str();      
}


// print function
// ----------
string t_gnavsbs::linefmt() const
{
  gtrace("t_gnavsbs::linefmt");
   
  ostringstream tmp;

  tmp  << " " << setw(3) << sat() << fixed
       << " " << _toc.str("%Y-%m-%d %H:%M:%S")
//     << " " << _toe.str("%Y-%m-%d %H:%M:%S")
//     << " " << _tot.str("%Y-%m-%d %H:%M:%S")
       << fixed    << setprecision(0)
  
       << setw( 8) << 0
       << setw( 8) << _C_rms
       << setw( 4) << _iod
       << setw( 4) << _health
       << " |"
       << setw(12) << setprecision(3) << _f0     *1e9  //   [sec]
       << setw( 8) << setprecision(3) << _f1     *1e9  //   [sec]
       << setw( 8) << setprecision(0) << _tki    *1e0  //
       << " |"
       << setw(14) << setprecision(3) << _x      *1e0  // 1 [km]
       << setw(11) << setprecision(3) << _x_d    *1e0  // 2 [km/s]
       << setw( 9) << setprecision(3) << _x_dd   *1e6  // 3 [km/s^2]
       << " |"
       << setw(14) << setprecision(3) << _y      *1e0  // 1 [km]
       << setw(11) << setprecision(3) << _y_d    *1e0  // 2 [km/s]
       << setw( 9) << setprecision(3) << _y_dd   *1e6  // 3 [km/s^2]
       << " |"
       << setw(14) << setprecision(3) << _z      *1e0  // 1 [km]
       << setw(11) << setprecision(3) << _z_d    *1e0  // 2 [km/s]
       << setw( 9) << setprecision(3) << _z_dd   *1e6  // 3 [km/s^2]
    ;

  return tmp.str();
}


// healthy check
bool t_gnavsbs::_healthy() const
{
   if( _health == 0 ) return true;
   return false;
}


} // namespace
