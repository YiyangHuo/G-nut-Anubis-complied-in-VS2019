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

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "gutils/gcommon.h"
#include "gutils/gconst.h"
#include "gmodels/gnwmsurf.h"
#include "gmodels/ggpt.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_gnwmsurf::t_gnwmsurf()
: t_gnwmbase(),
   _log(0),
   _decay_min(DECAY_MIN),
   _decay_max(DECAY_MAX),
   _vert_temp(dT_LINEAR),
   _vert_zwd(dW_DE),
   _surf_zwd(W_INTEGR)
{
  gtrace("t_gnwmsurf::constructor");
  id_type(t_gdata::NWMSURF);
  id_group(t_gdata::GRP_PRODUCT);
  _epoch.tsys(t_gtime::GPS);
}

// constructor
// ----------
t_gnwmsurf::t_gnwmsurf(set<MET_DATA> data, set<MET_SURF> surf)
: t_gnwmbase(),
   _log(0),
   _decay_min(DECAY_MIN),
   _decay_max(DECAY_MAX),
   _vert_temp(dT_LINEAR),
   _vert_zwd(dW_DE),
   _surf_zwd(W_INTEGR)
{
  gtrace("t_gnwmsurf::constructor(data,surf)");
  id_type(t_gdata::NWMSURF);
  id_group(t_gdata::GRP_PRODUCT);

  _set_dt = data;
  _set_sf = surf;
  _epoch.tsys(t_gtime::GPS);
}


// destructor
// ----------
t_gnwmsurf::~t_gnwmsurf(){}


// add value
// ----------
int t_gnwmsurf::add_surf(MET_SURF typ, double val, bool reset )
{
  double max = 999999;
  double min = -max;

  switch( typ )
  {	
   case LAT    : min = -90.0; max =  90.0;              // deg
                 while( val < min ) val+=(max-min);
                 while( val > max ) val-=(max-min);
                 break;
   case LON    : min =   0.0; max = 360.0;              // deg
                 while( val < min ) val+=(max-min);
                 while( val > max ) val-=(max-min);
                 break;
   case XPRJ   : 
   case YPRJ   : min = -999999999.0; max = -min; break; // m NO LIMITS
   case HEL    : min =   -250.0; max = 250000.0; break; // m
   case OROGR  : min =   -250.0; max =  10000.0; break; // m
   case GEOID  : min =   -150.0; max =    150.0; break; // m
   case T2M    : min =    100.0; max =    400.0; break; // K
   case TSK    : min =    100.0; max =    400.0; break; // K
   case MSURF  : min =    100.0; max =    400.0; break; // K
   case TSURF  : min =    100.0; max =    400.0; break; // K
   case WSURF  : min =      0.0; max =   1000.0; break; // mm
   case ESURF  : min =      0.0; max =   1000.0; break; // hPa
   case QSURF  : min =      0.0; max =   1000.0; break; // kg/kg
   case MSL    : min =    900.0; max =   1100.0; break; // hPa
   case LSM    : min =      0.0; max =      1.0; break; // 0/1
   case TCWV   : min =      0.0; max =   1000.0; break; // mm
   case ME_Ph  : min =     -0.5; max =      0.5; break; //  --> SHOULD BE 0 (CONTROL) !
     
   case eqH : case eqW : case eqE : min = 1.0; max = 100.0; break; // km
   case scH : case scW : case scE : min = 1.0; max =  25.0; break; // km 
     
   case dT    : // GENERAL
   case dM    : // GENERAL
   case MS_Tp : // two-parameter fit
   case MP_Tp : case  SP_Tp : case  AP_Tp : case ML_Th : case SL_Th : case  AL_Th :  // NEW
   case MP_Mp : case  SP_Mp : case  AP_Mp : case ML_Mh : case SL_Mh : case  AL_Mh :  // NEW
       
                 // min = -25.0; max = 500.0; // K/km
                 // if( val < min ) val = 0.0;
                 min = _decay_min;
                 max = _decay_max;
                 break;

   case dE    : // GENERAL
   case MS_Ep : // two-parameter fit
   case MP_Ep : case SP_Ep : case  AP_Ep : case ML_Ep : case SL_Ep : case AL_Ep : // NEW
   case MP_Eh : case SP_Eh : case  AP_Eh : case ML_Eh : case SL_Eh : case AL_Eh : // NEW
   case MP_Fp : case SP_Fp : case  AP_Fp : case ML_Fp : case SL_Fp : case AL_Fp : // NEW
   case MP_Fh : case SP_Fh : case  AP_Fh : case ML_Fh : case SL_Fh : case AL_Fh : // NEW
   case ME_Eh :

   case LR_Ep : case ER_Ep : case LL_Ep : case LM_Ep : // OLD
   case LR_Eh : case ER_Eh : case LL_Eh : case LM_Eh : // OLD
   case LR_Fp : case ER_Fp : case LL_Fp : case LM_Fp : // OLD
   case LR_Fh : case ER_Fh : case LL_Fh : case LM_Fh : // OLD
     
                 // min = -25.0; max = 2500.0;  // -
                 min = _decay_min;
                 max = _decay_max;
                 break;
     
   case dW    : // GENERAL
   case dI    : // GENERAL
   case MP_Wp : case SP_Wp : case AP_Wp : case ML_Wp : case SL_Wp : case AL_Wp : // NEW
   case MP_Wh : case SP_Wh : case AP_Wh : case ML_Wh : case SL_Wh : case AL_Wh : // NEW
   case ME_Wh :

   case MP_Ip : case SP_Ip : case AP_Ip : case ML_Ip : case SL_Ip : case AL_Ip : // NEW
   case MP_Ih : case SP_Ih : case AP_Ih : case ML_Ih : case SL_Ih : case AL_Ih : // NEW

   case MP_WE : case ML_WE : // NEW
  
   case MP_Wi : case SP_Wi : case ML_Wi : case SL_Wi : // NEW
   case MP_Fi : case SP_Fi : case ML_Fi : case SL_Fi : // NEW
       
   case LR_Wp : case ER_Wp : case LL_Wp : case LM_Wp : // OLD
   case LR_Wh : case ER_Wh : case LL_Wh : case LM_Wh : // OLD
   case ER_Wi : // OLD
   case ER_Fi : // OLD
   case LR_WE : // OLD

                 // min = -2.0; max = 250.0;  // -
                 min = _decay_min;
                 max = _decay_max;
                 break;
     
   case dR    : // GENERAL
   case MP_Rt : case ML_Rt :  // NEW
   case SP_Rt : case SL_Rt :  // NEW
   case LR_Rt : case ER_Rt :  // OLD
     
                 // min = -2.0; max = 100.0;  // -
                 min = _decay_min;
                 max = _decay_max;
                 break;
     
  };
   
  if( val < min || val > max ){
    ostringstream os;
    os << "# Warning: " << fixed << setw(6) << met_surf_str(typ) << setprecision(2)
       << "  surf: "    << val
       << "  range: "   << min << ":" << max
       << "  epo: "     << _epoch.str_ymdhms()
       << "  blh: "     << " " << setw(7) << get_surf(LAT)
                        << " " << setw(7) << get_surf(LON) 
	                << " " << setw(7) << get_surf(HEL);
    if( _log ) _log->comment(0,"gnwmsurf",os.str());
    else       cerr << os.str() << endl;
     
    if( !reset ) return 1;

    // reset
    if( val < min ){ _map_surf[typ] = min; }  // reset to minimum
    if( val > max ){ _map_surf[typ] = max; }  // reset to maximum
  }else{             _map_surf[typ] = val; }  // true value
  
  _set_sf.insert(typ); // if( _set_sf.find(typ) != _set_sf.end() ) _set_sf.insert(typ);

  return 0;
}


// add value
// ----------
int t_gnwmsurf::add_surf(MET_DATA typ, double val, bool reset )
{
  double max = 999999;
  double min = -max;

  // check values
  switch( typ )
  {	
   case GEOP : min = -250.0; max = 250000.0; break; // m
   case GEOM : min = -250.0; max = 250000.0; break; // m
   case GRAV : min =    9.0; max =     10.0; break; // m/s2
   case PRES : min =    0.0; max = 110000.0; break; // Pa
   case EPRE : min =    0.0; max = 100000.0; break; // Pa     // changed for geoid values
   case TEMP : min =  100.0; max =    500.0; break; // K      // changed for geoid values
   case TM   : min =  100.0; max =    500.0; break; // K      // changed for geoid values
   case SHUM : min =    0.0; max =    100.0; break; // g/kg
   case ZTD  : min =    0.0; max =    999.0; break; // m      // changed for geoid values
   case ZHD  : min =    0.0; max =    999.0; break; // m      // changed for geoid values
   case ZWD  : min =    0.0; max =    999.0; break; // m      // changed for geoid values
   case IWV  : min =    0.0; max =    999.0; break; // kg/m2
   case NGRD : min =   -9.9; max =      9.0; break; // mm
   case EGRD : min =   -9.9; max =      9.0; break; // mm
   case NGH  : min =   -9.9; max =      9.0; break; // mm
   case NGW  : min =   -9.9; max =      9.0; break; // mm
   case EGH  : min =   -9.9; max =      9.0; break; // mm
   case EGW  : min =   -9.9; max =      9.0; break; // mm
   case RATL : min =    0.0; max =    100.0; break; // -
   case RATE : min =    0.0; max =    100.0; break; // -
   case RATW : min =    0.0; max =    100.0; break; // -
   case RATX : min =    0.0; max =    100.0; break; // -
  }
   
  if( val < min || val > max ){
    ostringstream os;
    os << "# Warning: " << fixed << setw(6) << met_data_str(typ) << setprecision(2)
       << "  data: "    << val 
       << "  range: "   << min << ":" << max
       << "  epo: "     << _epoch.str_ymdhms()
       << "  blh: "     << " " << setw(7) << get_surf(LAT)
                        << " " << setw(7) << get_surf(LON)
	                << " " << setw(7) << get_surf(HEL);
    if( _log ) _log->comment(0,"gnwmsurf",os.str());
    else       cerr << os.str() << endl;
     
    if( !reset ) return 1;

    // reset
    if( val < min ){ _map_data[typ] = min; }  // reset to minimum
    if( val > max ){ _map_data[typ] = max; }  // reset to maximum
  }else{             _map_data[typ] = val; }  // true value

  _set_dt.insert(typ); //  if( _set_dt.find(typ) != _set_dt.end() ) _set_dt.insert(typ);

  return 0;
}


// convert typ
// ----------
int t_gnwmsurf::add_surfdata(MET_DATA typ, double val)
{
  switch ( typ ){
    case TM   : return add_surf(MSURF,val);
    case TEMP : return add_surf(TSURF,val);
    case SHUM : return add_surf(QSURF,val);
    case EPRE : return add_surf(ESURF,val);
    case ZWD  : return add_surf(WSURF,val);
    default   : return add_surf(typ,val);
  }
  return add_surf(typ, val);
}


// add value rms
// ----------
int t_gnwmsurf::add_rms(MET_SURF typ, double val)
{
  _map_rms[typ] = val;   
  return 0;
}


// add value rms
// ----------
int t_gnwmsurf::add_rms(MET_DATA typ, double val)
{
  _map_drms[typ] = val;
  return 0;
}

// add iterations
// ----------
int t_gnwmsurf::add_itr(MET_SURF typ, int val)
{
  _map_iter[typ] = val;
  return 0;
}


// get value
// ----------
double t_gnwmsurf::get_surf(MET_SURF typ)const
{
  map<MET_SURF, double>::const_iterator it = _map_surf.find( typ );   
  if( it == _map_surf.end() ) return NWM_UNKNOWN; // undefined value!
  return it->second;
}


// get value
// ----------
double t_gnwmsurf::get_surf(MET_DATA typ)const
{
  map<MET_DATA, double>::const_iterator it = _map_data.find( typ );
  if( it == _map_data.end() ) return NWM_UNKNOWN; // undefined value!
  return it->second;
}
   
   
// convert typ
// ----------
double t_gnwmsurf::get_surfdata(MET_DATA typ)const
{
  switch ( typ ){
    case TM   : return get_surf(MSURF); // use model value!
    case TEMP : return get_surf(TSURF); // use model value!
    case SHUM : return get_surf(QSURF); // use model value!
    case EPRE : return get_surf(ESURF); // use model value!
    case ZWD  : return get_surf(WSURF); // use model value!
    default   : return get_surf(typ);
  }
  return get_surf(typ);
}


// get value 
// ----------
double t_gnwmsurf::get_in_hgp(MET_DATA typ, double hgp)
{
  double val = NWM_UNKNOWN;
  double inp = get_surfdata(typ); if( inp == NWM_UNKNOWN ) return val;
   
  if(      get_surf(LAT)   == NWM_UNKNOWN ) cerr << "#warning: unknown LAT !\n";
  else if( get_surf(HEL)   == NWM_UNKNOWN ) cerr << "#warning: unknown HEL !\n";
  else if( get_surf(GEOID) == NWM_UNKNOWN ) cerr << "#warning: unknown GEOID !\n";
  else if( get_surf(TSURF) == NWM_UNKNOWN ) cerr << "#warning: unknown TSURF !\n";
  else{

    val = vert_scaling( typ, hgp,
			     geom2geop( get_surf(LAT), get_surf(HEL) - get_surf(GEOID)),
			     inp,
  	  	             get_surfdata(TEMP), // --> TSURF
			     get_surfdata(PRES) );
  }
   
  return val;
}

// get value 
// ----------
double t_gnwmsurf::get_in_hms(MET_DATA typ, double hms)
{
  double val = NWM_UNKNOWN;
  double inp = get_surfdata(typ); if( inp == NWM_UNKNOWN ) return val;

  if(      get_surf(LAT)   == NWM_UNKNOWN ) cerr << "#warning: unknown LAT !\n";
  else if( get_surf(HEL)   == NWM_UNKNOWN ) cerr << "#warning: unknown HEL !\n";
  else if( get_surf(GEOID) == NWM_UNKNOWN ) cerr << "#warning: unknown GEOID !\n";
  else if( get_surf(TSURF) == NWM_UNKNOWN ) cerr << "#warning: unknown TSURF !\n";
  else{

    val = vert_scaling( typ, geom2geop( get_surf(LAT), hms ),
 	  	  	     geom2geop( get_surf(LAT), get_surf(HEL) - get_surf(GEOID)),
   	                     inp,
			     get_surfdata(TEMP), // --> TSURF
			     get_surfdata(PRES) );
  }

  return val;
}

// get value 
// ----------
double t_gnwmsurf::get_in_hel(MET_DATA typ, double hel)
{
  double val = NWM_UNKNOWN;
  double inp = get_surfdata(typ); if( inp == NWM_UNKNOWN ) return val;
   
  if(      get_surf(LAT)   == NWM_UNKNOWN ) cerr << "#warning: unknown LAT !\n";
  else if( get_surf(HEL)   == NWM_UNKNOWN ) cerr << "#warning: unknown HEL !\n";
  else if( get_surf(GEOID) == NWM_UNKNOWN ) cerr << "#warning: unknown GEOID !\n";
  else if( get_surf(TSURF) == NWM_UNKNOWN ) cerr << "#warning: unknown TSURF !\n";
  else{
     
    val = vert_scaling( typ, geom2geop( get_surf(LAT), hel           - get_surf(GEOID) ),
			     geom2geop( get_surf(LAT), get_surf(HEL) - get_surf(GEOID) ),
   	                     inp,
			     get_surfdata(TEMP), // --> TSURF
			     get_surfdata(PRES) );
  }

#ifdef DEBUG
  cout << met_data_id(typ)
       << " HEL: "  << get_surf(HEL)
       << " GEOM: " << get_surf(GEOM) 
       << " N: "    << get_surf(GEOID) 
       << " T: "    << get_surfdata(TEMP)
       << " req: "  << geom2geop( get_surf(LAT), hel           - get_surf(GEOID) )
       << " ref: "  << geom2geop( get_surf(LAT), get_surf(HEL) - get_surf(GEOID) )
       << " srf: "  << get_surfdata(typ)
       << " val: "  << val
       << endl;
#endif

  return val;
}


// get value 
// ----------
double t_gnwmsurf::get_interp(MET_DATA typ, double hel)
{
  if(_log) _log->comment(1,"gnwmsurf", "Warning - direct vertical interpolation not supported.");

  return get_in_hel(typ,hel);
}


// get value rms
// ----------
double t_gnwmsurf::get_rms( MET_SURF typ )const
{
  map<MET_SURF, double>::const_iterator it = _map_rms.find( typ );
  if( it == _map_rms.end() ) return NWM_UNKNOWN; // undefined value!
  return it->second;
}


// get value rms
// ----------
double t_gnwmsurf::get_rms( MET_DATA typ )const
{
  map<MET_DATA, double>::const_iterator it = _map_drms.find( typ );
  if( it == _map_drms.end() ) return NWM_UNKNOWN; // undefined value!
  return it->second;
}


// get iteration
// ----------
double t_gnwmsurf::get_itr( MET_SURF typ )const
{
  map<MET_SURF, double>::const_iterator it = _map_iter.find( typ );
  if( it == _map_iter.end() ) return NWM_UNKNOWN; // undefined value!
  return it->second;
}


// vertical scaling
// ----------
double t_gnwmsurf::vert_scaling(MET_DATA type,
                                  double geopX,
                                  double geop,
                                  double data,
                                  double temp,
                                  double pres )
{ 
  gtrace("t_gnwmsurf::vert_scaling");
   
  if( data == NWM_UNKNOWN || geopX == NWM_UNKNOWN || 
      geop == NWM_UNKNOWN   )  return NWM_UNKNOWN;

  double res  = NWM_UNKNOWN;
  double lpsT = get_surf(dT);
  double lpsM = get_surf(dM);
  double lpsE = get_surf(dE);
  double lpsW = get_surf(dW);
  double lpsI = get_surf(dI);
  double lpsR = get_surf(dR);
  double sclH = get_surf(scH); // km

  // to be SURE
  if( pres == NWM_UNKNOWN ) pres = get_surf(PRES);  // assume to be at orography
  if( temp == NWM_UNKNOWN ) temp = get_surf(TEMP);  // assume to be at orography
  if( temp == NWM_UNKNOWN ) temp = get_surf(TSURF); // assume to be at orography
  if( lpsT == NWM_UNKNOWN ) lpsT = 6.5;
  if( lpsM == NWM_UNKNOWN ) lpsM = 6.5;
  if( sclH == NWM_UNKNOWN ) sclH = Rd*temp/G_WMO;   // troposphere scale height [km] approx 8km

  // temperature dependence
  if( _vert_temp == dT_CONST ){ lpsT = lpsM = 6.5; } //  temp = 288.15 - lpsT*geop2geom(lat,geop); } // USSA 1970 lapse rate
  if( _vert_temp == dT_ZERO  ){ lpsT = lpsM = 0.0; }

  switch( type )
  {
   case TEMP: //  temperature [K]
     res = data - lpsT/1000*(geopX-geop);
     break;
     
   case TM:   //  mean temperature [K]
     res = data - lpsM/1000*(geopX-geop);
     break;

   case PRES: //  pressure [hPa]
     if( _vert_temp == dT_ZERO ){ res = data * exp(  - (geopX-geop)/1000/sclH ); }
     else{                        res = data * pow(1 - lpsT/1000*(geopX-geop)/temp, G_WMO/(Rd*lpsT/1000)); }
     break;
     
   case EPRE: //  water vapour pressure [hPa]
     if( lpsE == NWM_UNKNOWN ) break;
     if( _vert_temp == dT_ZERO ){ res = data * exp(  - (geopX-geop)*(lpsE+1)/1000/sclH ); }
     else{                        res = data * pow(1 - lpsT/1000*(geopX-geop)/temp, (lpsE+1) * G_WMO/(Rd*lpsT/1000)); }
     break;

   case SHUM: //  specific humidity [g/kg]
     if( lpsE == NWM_UNKNOWN ) break;
     // use mixing ratio instead of spec.humidity (~to specific humidity scaling)
     // w = q/(1-q) ;  q = w/(1+w)
     data /= 1000; // [g/kg] --> [kg/kg]
     if( _vert_temp == dT_ZERO ){ res = data/(1-data) * exp(  - (geopX-geop)*(lpsE+1)/1000/sclH  ); }
     else{                        res = data/(1-data) * pow(1 - lpsT/1000*(geopX-geop)/temp, (lpsE) * G_WMO/(Rd*lpsT/1000)); }    
     res  = res/(1+res); // back conversion to spec.humidity
     res *= 1000; // [g/kg] --> [kg/kg]
     break;

   case ZHD: //  zenith hydrostatic delay (integrated) [m]
     // impact of g due to dH is negligible
     // -> dH = 1 km --> dZHD = ~ 1mm !
     // -> dH = 10km --> dZHD = ~ 1mm !
     if( _vert_temp == dT_ZERO ){ res = data * exp(  - (geopX-geop)/100/sclH ); }
     else{                        res = data * pow(1 - lpsT/1000*(geopX-geop)/temp, G_WMO/(Rd*lpsT/1000)); }

     // NOT NECESSARY, BUT CORRECT
     // ==========================
//     double pX   = pres;
//     double lat  = get_surf(LAT);
//     if( _vert_temp == dT_ZERO ){ pX *= exp(  - (geopX-geop)/100/sclH ); }
//     else{                        pX *= pow(1 - lpsT/1000*(geopX-geop)/temp, G_WMO/(Rd*lpsT/1000)); }
//     // Saastamoinen at req. geopX
//     res = 0.002277 * pX / (1.0 - 0.00266 * cos(2.0*lat*D2R) - 0.00000028 * (geopX));
          

#ifdef DEBUG
     cout << "\nZHD(vcor): " << fixed << setprecision(3)
          << " hel: "  << get_surf(HEL)
          << " gp0: "  << geop 
          << " gpX: "  << geopX
          << " ps: "   << get_surf(PRES)
          << " p0: "   << pres
          << " pX: "   << pX
          << " T:  "   << temp
          << " orig: " << data 
          << " newX: " << res
          << " laps: " << lpsT
          << " corr: " << pow(1 - lpsT/1000*(geopX-geop)/temp, G_WMO/(Rd*lpsT/1000))
          << " corr: " << pow(1 - lpsT/1000*(geopX-geop)/temp, G_WMO/(Rd*lpsT/1000))
          << endl;
#endif
     break;
     
   case ZWD: //  zenith wet delay [m]
     if( lpsW == NWM_UNKNOWN ) break;
     if( lpsE == NWM_UNKNOWN ) break;
     if( _vert_temp == dT_ZERO        ){ res = data * exp(  - (geopX-geop)*(lpsW+1)/1000/sclH ); }
     else if( _vert_zwd == dW_DE      ){ res = data * pow(1 - lpsT/1000*(geopX-geop)/temp, (lpsW+1) * G_WMO/(Rd*lpsT/1000)    ); }
     // MOPS-related using dE (Lambda)
     else if( _vert_zwd == dW_DE_MOPS ){ res = data * pow(1 - lpsT/1000*(geopX-geop)/temp, (lpsE+1) * G_WMO/(Rd*lpsT/1000)    ); }
     else if( _vert_zwd == dW_MOPS    ){ res = data * pow(1 - lpsT/1000*(geopX-geop)/temp, (lpsE+1) * G_WMO/(Rd*lpsT/1000) - 1); }
     break;

   case IWV: //  integrated water vapour [kg/m3]
     if( lpsI == NWM_UNKNOWN ) break;
     if( _vert_temp == dT_ZERO        ){ res = data * exp(  - (geopX-geop)*(lpsW+1)/1000/sclH ); }
     else{                               res = data * pow(1 - lpsT/1000*(geopX-geop)/temp, (lpsI+1) * G_WMO/(Rd*lpsT/1000)); }
     break;

   case RATX: case RATL : case RATE : case RATW: // gamma/lambda ratio [%]
     
     // DEFAULT (FIX) (already in gnwmprof.cpp)
//     if( lpsR == NWM_UNKNOWN || lpsR <= 0.0001 ) lpsR = RATIO_DH; // USE FIX: UNKNOWN OR ZERO
//     if( data == NWM_UNKNOWN                   ) data = RATIO_EW; // USE FIX: UNKNOWN 

     res = data - lpsR*(geopX-geop);
     if( res <   0.0 ) res =   0.0;
     if( res > 100.0 ) res = 100.0;
//   cout << "RATIO: " << data << " lpsR:" << lpsR << "  res: " << res << "  geop->geopX:  " << geop << " " << geopX << endl;
     break;

   default: cerr << "# --> warning - gnwsurf: unknown type for vertical scaling\n";
  }

  return res;
}


// calculate surface values
// ----------
int t_gnwmsurf::calc_surface()
{
  gtrace("t_gnwmsurf::_calc_surface");

  double lat  = get_surf(LAT);    // latitude [deg]
  double lon  = get_surf(LON);    // latitude [deg]
  double hel  = get_surf(HEL);    // elipsoidal height [m]
  double geom = get_surf(GEOM);   // ortometric height [m]
  double geop = get_surf(GEOP);   // geopotential height [gpm]

  if( lat  == NWM_UNKNOWN ||
      lon  == NWM_UNKNOWN ) return 1; // no Latitude !

  // add GEOID (undulation) // calculated at AT ANY EPOCH
  t_gpt gpt; 
  double pres, temp, undu;

  gpt.gpt_v1( 51544.0, lat*D2R, lon*D2R, geop, pres, temp, undu );
  add_surf( GEOID, undu );
   
  // calculate/add gravity at GEOM
  double r_lat =  radius_wgs84(lat);
  double g_lat = gravity_wgs84(lat);
  add_surf( GRAV, g_lat * pow( r_lat / (r_lat + geom ), 2.0) );

  if( geop == NWM_UNKNOWN &&
      geom == NWM_UNKNOWN &&
       hel == NWM_UNKNOWN ) return 1; // no Height - ANY !

  // calculate/add geopotential from geometric
  if( geop == NWM_UNKNOWN ){
    if( geom != NWM_UNKNOWN ) add_surf( GEOP, geom2geop(lat, geom) );
    if(  hel != NWM_UNKNOWN ) add_surf( GEOP, geom2geop(lat,  hel - undu) );
  }

  // calculate/add geometric from geopotential height
  if( geop != NWM_UNKNOWN ){
    if( geom == NWM_UNKNOWN ) add_surf( GEOM, geop2geom(lat, geop) );
    if(  hel == NWM_UNKNOWN ) add_surf(  HEL, geop2geom(lat, geop) + undu );
  }

  // copy geometric heights
  if( geom != NWM_UNKNOWN && hel  == NWM_UNKNOWN )  add_surf(  HEL, geom + undu );
  if(  hel != NWM_UNKNOWN && geom == NWM_UNKNOWN )  add_surf( GEOM,  hel - undu );
   
#ifdef DEBUG
  cout << " calc_surf: " << fixed << setprecision(3)
       << "  HEL: "  << setw(8) << get_surf( HEL )  << " [" << setw(8) << hel  << "]"
       << "  GEOM: " << setw(8) << get_surf( GEOM ) << " [" << setw(8) << geom << "]"
       << "  GEOP: " << setw(8) << get_surf( GEOP ) << " [" << setw(8) << geop << "]"
       << endl; cout.flush();
#endif

  return 0;
}


// WGS84 Earth radius(lat), Mahoney (2001)  // [m]
// ----------
double t_gnwmsurf::radius_wgs84(const double& lat) const
{
  double r_lat = A_WGS / (1 + F_WGS + G_RATIO - 2 * F_WGS * pow(sin(lat*D2R), 2.0));
  return r_lat;
}


// WGS84 gravity(lat), Mahoney (2001)        // [m/s2]
// ----------
double t_gnwmsurf::gravity_wgs84(const double& lat) const
{
  double g_lat = G_EQUA * (1 +          KS_SOM * pow(sin(lat*D2R), 2.0))
                   / (sqrt(1 - pow(E_WGS, 2.0) * pow(sin(lat*D2R), 2.0)));
  return g_lat;
}


// geop [gpm] -> geom [m], Mahoney (2001)
// ----------
double t_gnwmsurf::geop2geom(const double& lat, const double& geop) const
{
  double r_lat = radius_wgs84(lat);
  double g_lat = gravity_wgs84(lat);

//double geom = (A_WGS * geop) / (A_WGS - geop);               // simple form (acc roughly 20cm)
  double geom = (r_lat * geop) / ((g_lat/G_WMO)*r_lat - geop); // Mahoney form (acc ~2mm up to 20km!)
  return geom;
}


// geom [m] -> geop [gpm], Mahoney (2001)
// ----------
double t_gnwmsurf::geom2geop(const double& lat, const double& geom) const
{
  double r_lat = radius_wgs84(lat);
  double g_lat = gravity_wgs84(lat);

//double geop = (A_WGS * geom) / (A_WGS + geom);               // simple form (acc roughly 20cm)
  double geop = g_lat/G_WMO * (r_lat * geom)/(r_lat + geom);   // Mahoney form (acc ~2mm up to 20km!)
  return geop;
}


// printf surface values
void t_gnwmsurf::print_surface() const
{
  gtrace("t_gnwmsurf::print_surface");
   
  t_map_metsurf::iterator it;

//  double zhd01 = (0.002277 * get_surf(PRES0))
//               / (1.0 - 0.00266 * cos(2.0*get_surf(LAT)*D2R) - 0.00000028 * get_surf(GEOM0));
//  double zhd02 = (0.0022768 * get_surf(PRES0))
//               / (1.0 - 0.0026  * cos(2.0*get_surf(LAT)*D2R) - 0.00000028 * get_surf(GEOM0));

  cout << "#beg_surface\n" << fixed     << setprecision(3)
       << "#site: " << "NONAME"                   << endl;
  cout << "#date: " << _epoch.str_ymdhms()        << endl;
  cout << "#lat:  " << setw(9) << get_surf(LAT)   << endl;
  cout << "#lon:  " << setw(9) << get_surf(LON)   << endl;
  cout << "#hel:  " << setw(9) << get_surf(HEL)   << endl;

  if( _set_dt.find(GEOM) != _set_dt.end() )  cout << "#geom: " << setw(9) << get_surf(GEOM)  << endl;
  if( _set_dt.find(GEOP) != _set_dt.end() )  cout << "#geop: " << setw(9) << get_surf(GEOP)  << endl;
  if( _set_sf.find(GEOID)!= _set_sf.end() )  cout << "#geoid:" << setw(9) << get_surf(GEOID) << endl;
  if( _set_dt.find(PRES) != _set_dt.end() )  cout << "#pres: " << setw(9) << get_surf(PRES)  << endl;
  if( _set_dt.find(EPRE) != _set_dt.end() )  cout << "#epre: " << setw(9) << get_surf(EPRE)  << endl;
  if( _set_dt.find(SHUM) != _set_dt.end() )  cout << "#shum: " << setw(9) << get_surf(SHUM)  << endl;
  if( _set_dt.find(TEMP) != _set_dt.end() )  cout << "#temp: " << setw(9) << get_surf(TEMP)  << endl;
  if( _set_dt.find(TM)   != _set_dt.end() )  cout << "#tm:   " << setw(9) << get_surf(TM)    << endl;
  if( _set_sf.find(T2M)  != _set_sf.end() )  cout << "#t2m:  " << setw(9) << get_surf(T2M)   << endl;
  if( _set_sf.find(TSK)  != _set_sf.end() )  cout << "#tsk:  " << setw(9) << get_surf(TSK)   << endl;
  if( _set_dt.find(ZTD)  != _set_dt.end() )  cout << "#ztd:  " << setw(9) << get_surf(ZHD)+get_surf(ZWD)   << endl;
  if( _set_dt.find(ZHD)  != _set_dt.end() )  cout << "#zhd:  " << setw(9) << get_surf(ZHD)   << endl;
  if( _set_dt.find(ZWD)  != _set_dt.end() )  cout << "#zwd:  " << setw(9) << get_surf(ZWD)   << endl;
  if( _set_dt.find(IWV)  != _set_dt.end() )  cout << "#iwv:  " << setw(9) << get_surf(IWV)   << endl;
     
  if( _set_sf.find(TCWV) != _set_sf.end() )  cout << "#tcwv: " << setw(9) << get_surf(TCWV)  << endl;
  if( _set_sf.find(MSL)  != _set_sf.end() )  cout << "#msl:  " << setw(9) << get_surf(MSL)   << endl;
  if( _set_sf.find(LSM)  != _set_sf.end() )  cout << "#lsm:  " << setw(9) << get_surf(LSM)   << endl;

  if( _set_dt.find(RATX) != _set_dt.end() )  cout << "#ratX: " << setw(9) << get_surf(RATX)  << endl;
  if( _set_dt.find(RATL) != _set_dt.end() )  cout << "#ratL: " << setw(9) << get_surf(RATL)  << endl;
  if( _set_dt.find(RATE) != _set_dt.end() )  cout << "#ratE: " << setw(9) << get_surf(RATE)  << endl;
  if( _set_dt.find(RATW) != _set_dt.end() )  cout << "#ratW: " << setw(9) << get_surf(RATW)  << endl;

  if( _set_sf.find(TSURF)!= _set_sf.end() )  cout << "#Ts:   " << setw(9) << get_surf(TSURF) << endl;  // or get_surfdata(TEMP)
  if( _set_sf.find(ESURF)!= _set_sf.end() )  cout << "#Es:   " << setw(9) << get_surf(ESURF) << endl;  // or get_surfdata(EPRE)
  if( _set_sf.find(QSURF)!= _set_sf.end() )  cout << "#Qs:   " << setw(9) << get_surf(QSURF) << endl;  // or get_surfdata(IWV)
  if( _set_sf.find(WSURF)!= _set_sf.end() )  cout << "#Ws:   " << setw(9) << get_surf(WSURF) << endl;  // or get_surfdata(ZWD)
  if( _set_sf.find(MSURF)!= _set_sf.end() )  cout << "#Ms:   " << setw(9) << get_surf(MSURF) << endl;  // or get_surfdata(TM)

  cout << setprecision(2);

  // general lapse/decay rates
  if( _set_sf.find(dT)    != _set_sf.end() ) cout << "#dT:   " << setw(9) << get_surf(dT)
                                                               << setw(9) << get_rms( dT)
	                                                       << setw(9) << get_itr( dT) << endl;
  if( _set_sf.find(dM)    != _set_sf.end() ) cout << "#dM:   " << setw(9) << get_surf(dM)
                                                               << setw(9) << get_rms( dM)
	                                                       << setw(9) << get_itr( dM) << endl;
  if( _set_sf.find(dE)    != _set_sf.end() ) cout << "#dE:   " << setw(9) << get_surf(dE)
                                                               << setw(9) << get_rms( dE)
	                                                       << setw(9) << get_itr( dE) << endl;
  if( _set_sf.find(dW)    != _set_sf.end() ) cout << "#dW:   " << setw(9) << get_surf(dW)
                                                               << setw(9) << get_rms( dW)
	                                                       << setw(9) << get_itr( dW) << endl;
  if( _set_sf.find(dI)    != _set_sf.end() ) cout << "#dI:   " << setw(9) << get_surf(dI)
                                                               << setw(9) << get_rms( dI)
	                                                       << setw(9) << get_itr( dI) << endl;
  if( _set_sf.find(dR)    != _set_sf.end() ) cout << "#dR:   " << setw(9) << get_surf(dR)
                                                               << setw(9) << get_rms( dR)
	                                                       << setw(9) << get_itr( dR) << endl;

  if( _set_sf.find(eqH)   != _set_sf.end() ) cout << "#eqH:  " << setw(9) << get_surf(eqH)
                                                               << setw(9) << get_rms( eqH)
	                                                       << setw(9) << get_itr( eqH) << endl;
  if( _set_sf.find(eqW)   != _set_sf.end() ) cout << "#eqW:  " << setw(9) << get_surf(eqW)
                                                               << setw(9) << get_rms( eqW)
	                                                       << setw(9) << get_itr( eqW) << endl;
  if( _set_sf.find(eqE)   != _set_sf.end() ) cout << "#eqE:  " << setw(9) << get_surf(eqE)
                                                               << setw(9) << get_rms( eqE)
	                                                       << setw(9) << get_itr( eqE) << endl;

  if( _set_sf.find(scH)   != _set_sf.end() ) cout << "#eqH:  " << setw(9) << get_surf(scH)
                                                               << setw(9) << get_rms( scH)
	                                                       << setw(9) << get_itr( scH) << endl;
  if( _set_sf.find(scW)   != _set_sf.end() ) cout << "#scW:  " << setw(9) << get_surf(scW)
                                                               << setw(9) << get_rms( scW)
	                                                       << setw(9) << get_itr( scW) << endl;
  if( _set_sf.find(scE)   != _set_sf.end() ) cout << "#scE:  " << setw(9) << get_surf(scE)
                                                               << setw(9) << get_rms( scE)
	                                                       << setw(9) << get_itr( scE) << endl;

  // control (
  if( _set_sf.find(ME_Ph) != _set_sf.end() ) cout << "#mePh: " << setw(9) << get_surf(ME_Ph)
                                                               << setw(9) << get_rms( ME_Ph)
	                                                       << setw(9) << get_itr( ME_Ph) << endl;

  // T lapse rate
  if( _set_sf.find(MS_Tp) != _set_sf.end() ) cout << "#msTp: " << setw(9) << get_surf(MS_Tp)
                                                               << setw(9) << get_rms( MS_Tp)
	                                                       << setw(9) << get_itr( MS_Tp) 
	                      << setw(20) << "adj_surf_TEMP: " << setw(9) << get_surf(TSURF) << endl;
   
  if( _set_sf.find(MP_Tp) != _set_sf.end() ) cout << "#mpTp: " << setw(9) << get_surf(MP_Tp)
                                                               << setw(9) << get_rms( MP_Tp)
	                                                       << setw(9) << get_itr( MP_Tp) << endl;
  if( _set_sf.find(SP_Tp) != _set_sf.end() ) cout << "#spTp: " << setw(9) << get_surf(SP_Tp)
                                                               << setw(9) << get_rms( SP_Tp)
	                                                       << setw(9) << get_itr( SP_Tp) << endl;
  if( _set_sf.find(AP_Tp) != _set_sf.end() ) cout << "#apTp: " << setw(9) << get_surf(AP_Tp)
                                                               << setw(9) << get_rms( AP_Tp)
	                                                       << setw(9) << get_itr( AP_Tp) << endl;
  if( _set_sf.find(ML_Th) != _set_sf.end() ) cout << "#mlTh: " << setw(9) << get_surf(ML_Th)
                                                               << setw(9) << get_rms( ML_Th)
	                                                       << setw(9) << get_itr( ML_Th) << endl;
  if( _set_sf.find(SL_Th) != _set_sf.end() ) cout << "#slTh: " << setw(9) << get_surf(SL_Th)
                                                               << setw(9) << get_rms( SL_Th)
	                                                       << setw(9) << get_itr( SL_Th) << endl;
  if( _set_sf.find(AL_Th) != _set_sf.end() ) cout << "#alTh: " << setw(9) << get_surf(AL_Th)
                                                               << setw(9) << get_rms( AL_Th)
	                                                       << setw(9) << get_itr( AL_Th) << endl;
   
  // Tm lapse rate
  if( _set_sf.find(MP_Mp) != _set_sf.end() ) cout << "#mpMp: " << setw(9) << get_surf(MP_Mp)
                                                               << setw(9) << get_rms( MP_Mp)
	                                                       << setw(9) << get_itr( MP_Mp) << endl;
  if( _set_sf.find(SP_Mp) != _set_sf.end() ) cout << "#spMp: " << setw(9) << get_surf(SP_Mp)
                                                               << setw(9) << get_rms( SP_Mp)
	                                                       << setw(9) << get_itr( SP_Mp) << endl;
  if( _set_sf.find(AP_Mp) != _set_sf.end() ) cout << "#alMp: " << setw(9) << get_surf(AP_Mp)
                                                               << setw(9) << get_rms( AP_Mp)
	                                                       << setw(9) << get_itr( AP_Mp) << endl;
  if( _set_sf.find(ML_Mh) != _set_sf.end() ) cout << "#mlMh: " << setw(9) << get_surf(ML_Mh)
                                                               << setw(9) << get_rms( ML_Mh)
	                                                       << setw(9) << get_itr( ML_Mh) << endl;
  if( _set_sf.find(SL_Mh) != _set_sf.end() ) cout << "#slMh: " << setw(9) << get_surf(SL_Mh)
                                                               << setw(9) << get_rms( SL_Mh)
	                                                       << setw(9) << get_itr( SL_Mh) << endl;
  if( _set_sf.find(AL_Mh) != _set_sf.end() ) cout << "#alMh: " << setw(9) << get_surf(AL_Mh)
                                                               << setw(9) << get_rms( AL_Mh)
	                                                       << setw(9) << get_itr( AL_Mh) << endl;

//  cout << "#lrTp: " << setw(9) << get_surf(LR_Tp) << setw(9) << get_rms( LR_Tp) << setw(9) << get_itr(LR_Tp) << endl;
//  cout << "#lmTp: " << setw(9) << get_surf(LM_Tp) << setw(9) << get_rms( LM_Tp) << setw(9) << get_itr(LM_Tp) << endl;
//  cout << "#erTp: " << setw(8) << get_surf(ER_Tp) << setw(9) << get_rms( ER_Tp) << setw(9) << get_itr(ER_Tp) << endl;
//  cout << "#llTp: " << setw(8) << get_surf(LL_Tp) << setw(9) << get_rms( LL_Tp) << endl;
//  cout << "#lrMp: " << setw(8) << get_surf(LR_Mp) << setw(9) << get_rms( LR_Mp) << setw(9) << get_itr(LR_Mp) << endl;
//  cout << "#lmMp: " << setw(8) << get_surf(LM_Mp) << setw(9) << get_rms( LM_Mp) << setw(9) << get_itr(LM_Mp) << endl;
//  cout << "#erMp: " << setw(8) << get_surf(ER_Mp) << setw(9) << get_rms( ER_Mp) << setw(9) << get_itr(ER_Mp) << endl;
//  cout << "#llMp: " << setw(8) << get_surf(LL_Mp) << setw(9) << get_rms( LL_Mp) << endl;

   
  // WV decay 
  if( _set_sf.find(MS_Ep) != _set_sf.end() ) cout << "#msEp: " << setw(9) << get_surf(MS_Ep)
                                                               << setw(9) << get_rms( MS_Ep)
	                                                       << setw(9) << get_itr( MS_Ep) 
	                      << setw(20) << "adj_surf_EPRE: " << setw(9) << get_surf(ESURF) << endl;
   
  if( _set_sf.find(MP_Ep) != _set_sf.end() ) cout << "#mpEp: " << setw(9) << get_surf(MP_Ep)
                                                               << setw(9) << get_rms( MP_Ep)
	                                                       << setw(9) << get_itr( MP_Ep) << endl;
  if( _set_sf.find(SP_Ep) != _set_sf.end() ) cout << "#spEp: " << setw(9) << get_surf(SP_Ep)
                                                               << setw(9) << get_rms( SP_Ep)
	                                                       << setw(9) << get_itr( SP_Ep) << endl;
  if( _set_sf.find(AP_Ep) != _set_sf.end() ) cout << "#apEp: " << setw(9) << get_surf(AP_Ep)
                                                               << setw(9) << get_rms( AP_Ep)
	                                                       << setw(9) << get_itr( AP_Ep) << endl;
  if( _set_sf.find(ML_Ep) != _set_sf.end() ) cout << "#mlEp: " << setw(9) << get_surf(ML_Ep)
                                                               << setw(9) << get_rms( ML_Ep)
	                                                       << setw(9) << get_itr( ML_Ep) << endl;
  if( _set_sf.find(SL_Ep) != _set_sf.end() ) cout << "#slEp: " << setw(9) << get_surf(SL_Ep)
                                                               << setw(9) << get_rms( SL_Ep)
	                                                       << setw(9) << get_itr( SL_Ep) << endl;
  if( _set_sf.find(AL_Ep) != _set_sf.end() ) cout << "#alEp: " << setw(9) << get_surf(AL_Ep)
                                                               << setw(9) << get_rms( AL_Ep)
	                                                       << setw(9) << get_itr( AL_Ep) << endl;
  if( _set_sf.find(MP_Eh) != _set_sf.end() ) cout << "#mpEh: " << setw(9) << get_surf(MP_Eh)
                                                               << setw(9) << get_rms( MP_Eh)
	                                                       << setw(9) << get_itr( MP_Eh) << endl; 
  if( _set_sf.find(SP_Eh) != _set_sf.end() ) cout << "#spEh: " << setw(9) << get_surf(SP_Eh)
                                                               << setw(9) << get_rms( SP_Eh)
	                                                       << setw(9) << get_itr( SP_Eh) << endl;
  if( _set_sf.find(AP_Eh) != _set_sf.end() ) cout << "#apEh: " << setw(9) << get_surf(AP_Eh)
                                                               << setw(9) << get_rms( AP_Eh)
	                                                       << setw(9) << get_itr( AP_Eh) << endl;
  if( _set_sf.find(ML_Eh) != _set_sf.end() ) cout << "#mlEh: " << setw(9) << get_surf(ML_Eh)
                                                               << setw(9) << get_rms( ML_Eh)
	                                                       << setw(9) << get_itr( ML_Eh) << endl;
  if( _set_sf.find(SL_Eh) != _set_sf.end() ) cout << "#slEh: " << setw(9) << get_surf(SL_Eh)
                                                               << setw(9) << get_rms( SL_Eh)
	                                                       << setw(9) << get_itr( SL_Eh) << endl;
  if( _set_sf.find(AL_Eh) != _set_sf.end() ) cout << "#alEh: " << setw(9) << get_surf(AL_Eh)
                                                               << setw(9) << get_rms( AL_Eh)
	                                                       << setw(9) << get_itr( AL_Eh) << endl;
   
  // OLD
  if( _set_sf.find(LR_Eh) != _set_sf.end() ) cout << "#lrEh: " << setw(9) << get_surf(LR_Eh)
                                                               << setw(9) << get_rms( LR_Eh)
	                                                       << setw(9) << get_itr( LR_Eh) << endl;
  if( _set_sf.find(ER_Eh) != _set_sf.end() ) cout << "#erEh: " << setw(9) << get_surf(ER_Eh)
                                                               << setw(9) << get_rms( ER_Eh)
	                                                       << setw(9) << get_itr( ER_Eh) << endl;
  if( _set_sf.find(LL_Eh) != _set_sf.end() ) cout << "#llEh: " << setw(9) << get_surf(LL_Eh)
                                                               << setw(9) << get_rms( LL_Eh)
	                                                       << setw(9) << get_itr( LL_Eh) << endl;
  if( _set_sf.find(LM_Eh) != _set_sf.end() ) cout << "#lmEh: " << setw(9) << get_surf(LM_Eh)
                                                               << setw(9) << get_rms( LM_Eh)
	                                                       << setw(9) << get_itr( LM_Eh) << endl;
  if( _set_sf.find(LR_Ep) != _set_sf.end() ) cout << "#lrEp: " << setw(9) << get_surf(LR_Ep)
                                                               << setw(9) << get_rms( LR_Ep)
	                                                       << setw(9) << get_itr( LR_Ep) << endl;
  if( _set_sf.find(ER_Ep) != _set_sf.end() ) cout << "#erEp: " << setw(9) << get_surf(ER_Ep)
                                                               << setw(9) << get_rms( ER_Ep)
	                                                       << setw(9) << get_itr( ER_Ep) << endl;
  if( _set_sf.find(LL_Ep) != _set_sf.end() ) cout << "#llEp: " << setw(9) << get_surf(LL_Ep)
                                                               << setw(9) << get_rms( LL_Ep)
	                                                       << setw(9) << get_itr( LL_Ep) << endl;
  if( _set_sf.find(LM_Ep) != _set_sf.end() ) cout << "#lmEp: " << setw(9) << get_surf(LM_Ep)
                                                               << setw(9) << get_rms( LM_Ep)
	                                                       << setw(9) << get_itr( LM_Ep) << endl;
   
  if( _set_sf.find(LR_Fh) != _set_sf.end() ) cout << "#lrFh: " << setw(9) << get_surf(LR_Fh)
                                                               << setw(9) << get_rms( LR_Fh)
	                                                       << setw(9) << get_itr( LR_Fh) << endl;
  if( _set_sf.find(ER_Fh) != _set_sf.end() ) cout << "#erFh: " << setw(9) << get_surf(ER_Fh)
                                                               << setw(9) << get_rms( ER_Fh)
	                                                       << setw(9) << get_itr( ER_Fh) << endl;
  if( _set_sf.find(LL_Fh) != _set_sf.end() ) cout << "#llFh: " << setw(9) << get_surf(LL_Fh)
                                                               << setw(9) << get_rms( LL_Fh)
	                                                       << setw(9) << get_itr( LL_Fh) << endl;
  if( _set_sf.find(LM_Fh) != _set_sf.end() ) cout << "#lmFh: " << setw(9) << get_surf(LM_Fh)
                                                               << setw(9) << get_rms( LM_Fh)
	                                                       << setw(9) << get_itr( LM_Fh) << endl;
  if( _set_sf.find(LR_Fp) != _set_sf.end() ) cout << "#lrFp: " << setw(9) << get_surf(LR_Fp)
                                                               << setw(9) << get_rms( LR_Fp)
	                                                       << setw(9) << get_itr( LR_Fp) << endl;
  if( _set_sf.find(ER_Fp) != _set_sf.end() ) cout << "#erFp: " << setw(9) << get_surf(ER_Fp)
                                                               << setw(9) << get_rms( ER_Fp)
	                                                       << setw(9) << get_itr( ER_Fp) << endl;
  if( _set_sf.find(LL_Fp) != _set_sf.end() ) cout << "#llFp: " << setw(9) << get_surf(LL_Fp)
                                                               << setw(9) << get_rms( LL_Fp)
	                                                       << setw(9) << get_itr( LL_Fp) << endl;
  if( _set_sf.find(LM_Fp) != _set_sf.end() ) cout << "#lmFp: " << setw(9) << get_surf(LM_Fp)
                                                               << setw(9) << get_rms( LM_Fp)
	                                                       << setw(9) << get_itr( LM_Fp) << endl;
   
  // ZWD decay
  if( _set_sf.find(ME_Wh) != _set_sf.end() ) cout << "#meWh: " << setw(9) << get_surf(ME_Wh)
                                                               << setw(9) << get_rms( ME_Wh)
	                                                       << setw(9) << get_itr( ME_Wh) << endl;
  if( _set_sf.find(MP_Wp) != _set_sf.end() ) cout << "#mpWp: " << setw(9) << get_surf(MP_Wp)
                                                               << setw(9) << get_rms( MP_Wp)
	                                                       << setw(9) << get_itr( MP_Wp) << endl;
  if( _set_sf.find(SP_Wp) != _set_sf.end() ) cout << "#spWp: " << setw(9) << get_surf(SP_Wp)
                                                               << setw(9) << get_rms( SP_Wp)
	                                                       << setw(9) << get_itr( SP_Wp) << endl;
  if( _set_sf.find(AP_Wp) != _set_sf.end() ) cout << "#apWp: " << setw(9) << get_surf(AP_Wp)
                                                               << setw(9) << get_rms( AP_Wp)
	                                                       << setw(9) << get_itr( AP_Wp) << endl;
  if( _set_sf.find(ML_Wp) != _set_sf.end() ) cout << "#mlWp: " << setw(9) << get_surf(ML_Wp)
                                                               << setw(9) << get_rms( ML_Wp)
	                                                       << setw(9) << get_itr( ML_Wp) << endl;
  if( _set_sf.find(SL_Wp) != _set_sf.end() ) cout << "#slWp: " << setw(9) << get_surf(SL_Wp)
                                                               << setw(9) << get_rms( SL_Wp)
	                                                       << setw(9) << get_itr( SL_Wp) << endl;
  if( _set_sf.find(AL_Wp) != _set_sf.end() ) cout << "#alWp: " << setw(9) << get_surf(AL_Wp)
                                                               << setw(9) << get_rms( AL_Wp)
	                                                       << setw(9) << get_itr( AL_Wp) << endl;
  if( _set_sf.find(MP_Wh) != _set_sf.end() ) cout << "#mpWh: " << setw(9) << get_surf(MP_Wh)
                                                               << setw(9) << get_rms( MP_Wh)
	                                                       << setw(9) << get_itr( MP_Wh) << endl;
  if( _set_sf.find(SP_Wh) != _set_sf.end() ) cout << "#spWh: " << setw(9) << get_surf(SP_Wh)
                                                               << setw(9) << get_rms( SP_Wh)
	                                                       << setw(9) << get_itr( SP_Wh) << endl;
  if( _set_sf.find(AP_Wh) != _set_sf.end() ) cout << "#apWh: " << setw(9) << get_surf(AP_Wh)
                                                               << setw(9) << get_rms( AP_Wh)
	                                                       << setw(9) << get_itr( AP_Wh) << endl;
  if( _set_sf.find(ML_Wh) != _set_sf.end() ) cout << "#mlWh: " << setw(9) << get_surf(ML_Wh)
                                                               << setw(9) << get_rms( ML_Wh)
	                                                       << setw(9) << get_itr( ML_Wh) << endl;
  if( _set_sf.find(SL_Wh) != _set_sf.end() ) cout << "#slWh: " << setw(9) << get_surf(SL_Wh)
                                                               << setw(9) << get_rms( SL_Wh)
	                                                       << setw(9) << get_itr( SL_Wh) << endl;
  if( _set_sf.find(AL_Wh) != _set_sf.end() ) cout << "#alWh: " << setw(9) << get_surf(AL_Wh)
                                                               << setw(9) << get_rms( AL_Wh)
	                                                       << setw(9) << get_itr( AL_Wh) << endl;

  if( _set_sf.find(LR_Wh) != _set_sf.end() ) cout << "#lrWh: " << setw(9) << get_surf(LR_Wh)
                                                               << setw(9) << get_rms( LR_Wh)
	                                                       << setw(9) << get_itr( LR_Wh) << endl;
  if( _set_sf.find(ER_Wh) != _set_sf.end() ) cout << "#erWh: " << setw(9) << get_surf(ER_Wh)
                                                               << setw(9) << get_rms( ER_Wh)
	                                                       << setw(9) << get_itr( ER_Wh) << endl;
  if( _set_sf.find(LL_Wh) != _set_sf.end() ) cout << "#llWh: " << setw(9) << get_surf(LL_Wh)
                                                               << setw(9) << get_rms( LL_Wh)
	                                                       << setw(9) << get_itr( LL_Wh) << endl;
  if( _set_sf.find(LM_Wh) != _set_sf.end() ) cout << "#lmWh: " << setw(9) << get_surf(LM_Wh)
                                                               << setw(9) << get_rms( LM_Wh)
	                                                       << setw(9) << get_itr( LM_Wh) << endl;
  if( _set_sf.find(LR_Wp) != _set_sf.end() ) cout << "#lrWp: " << setw(9) << get_surf(LR_Wp)
                                                               << setw(9) << get_rms( LR_Wp)
	                                                       << setw(9) << get_itr( LR_Wp) << endl;
  if( _set_sf.find(ER_Wp) != _set_sf.end() ) cout << "#erWp: " << setw(9) << get_surf(ER_Wp)
                                                               << setw(9) << get_rms( ER_Wp)
	                                                       << setw(9) << get_itr( ER_Wp) << endl;
  if( _set_sf.find(LL_Wp) != _set_sf.end() ) cout << "#llWp: " << setw(9) << get_surf(LL_Wp)
                                                               << setw(9) << get_rms( LL_Wp)
	                                                       << setw(9) << get_itr( LL_Wp) << endl;
  if( _set_sf.find(LM_Wp) != _set_sf.end() ) cout << "#lmWp: " << setw(9) << get_surf(LM_Wp)
                                                               << setw(9) << get_rms( LM_Wp)
	                                                       << setw(9) << get_itr( LM_Wp) << endl;

  // IWV decay
  if( _set_sf.find(MP_Ip) != _set_sf.end() ) cout << "#mpIp: " << setw(9) << get_surf(MP_Ip)
                                                               << setw(9) << get_rms( MP_Ip)
	                                                       << setw(9) << get_itr( MP_Ip) << endl;
  if( _set_sf.find(SP_Ip) != _set_sf.end() ) cout << "#spIp: " << setw(9) << get_surf(SP_Ip)
                                                               << setw(9) << get_rms( SP_Ip)
	                                                       << setw(9) << get_itr( SP_Ip) << endl;
  if( _set_sf.find(AP_Ip) != _set_sf.end() ) cout << "#apIp: " << setw(9) << get_surf(AP_Ip)
                                                               << setw(9) << get_rms( AP_Ip)
	                                                       << setw(9) << get_itr( AP_Ip) << endl;
  if( _set_sf.find(ML_Ip) != _set_sf.end() ) cout << "#mlIp: " << setw(9) << get_surf(ML_Ip)
                                                               << setw(9) << get_rms( ML_Ip)
	                                                       << setw(9) << get_itr( ML_Ip) << endl;
  if( _set_sf.find(SL_Ip) != _set_sf.end() ) cout << "#slIp: " << setw(9) << get_surf(SL_Ip)
                                                               << setw(9) << get_rms( SL_Ip)
	                                                       << setw(9) << get_itr( SL_Ip) << endl;
  if( _set_sf.find(AL_Ip) != _set_sf.end() ) cout << "#alIp: " << setw(9) << get_surf(AL_Ip)
                                                               << setw(9) << get_rms( AL_Ip)
	                                                       << setw(9) << get_itr( AL_Ip) << endl;
  if( _set_sf.find(MP_Ih) != _set_sf.end() ) cout << "#mpIh: " << setw(9) << get_surf(MP_Ih)
                                                               << setw(9) << get_rms( MP_Ih)
	                                                       << setw(9) << get_itr( MP_Ih) << endl;
  if( _set_sf.find(SP_Ih) != _set_sf.end() ) cout << "#spIh: " << setw(9) << get_surf(SP_Ih)
                                                               << setw(9) << get_rms( SP_Ih)
	                                                       << setw(9) << get_itr( SP_Ih) << endl;
  if( _set_sf.find(AP_Ih) != _set_sf.end() ) cout << "#apIh: " << setw(9) << get_surf(AP_Ih)
                                                               << setw(9) << get_rms( AP_Ih)
	                                                       << setw(9) << get_itr( AP_Ih) << endl;
  if( _set_sf.find(ML_Ih) != _set_sf.end() ) cout << "#mlIh: " << setw(9) << get_surf(ML_Ih)
                                                               << setw(9) << get_rms( ML_Ih)
	                                                       << setw(9) << get_itr( ML_Ih) << endl;
  if( _set_sf.find(SL_Ih) != _set_sf.end() ) cout << "#slIh: " << setw(9) << get_surf(SL_Ih)
                                                               << setw(9) << get_rms( SL_Ih)
	                                                       << setw(9) << get_itr( SL_Ih) << endl;
  if( _set_sf.find(AL_Ih) != _set_sf.end() ) cout << "#alIh: " << setw(9) << get_surf(AL_Ih)
                                                               << setw(9) << get_rms( AL_Ih)
	                                                       << setw(9) << get_itr( AL_Ih) << endl;

  // TEST
  if( _set_sf.find(ER_Wi) != _set_sf.end() ) cout << "#erWi: " << setw(9) << get_surf(ER_Wi)
                                                               << setw(9) << get_rms( ER_Wi)
	                                                       << setw(9) << get_itr( ER_Wi) << endl;
  if( _set_sf.find(ER_Fi) != _set_sf.end() ) cout << "#erFi: " << setw(9) << get_surf(ER_Fi)
                                                               << setw(9) << get_rms( ER_Fi)
	                                                       << setw(9) << get_itr( ER_Fi) << endl;

  if( _set_sf.find(LR_WE) != _set_sf.end() ) cout << "#lrWE: " << setw(9) << get_surf(LR_WE)
                                                               << setw(9) << get_rms( LR_WE)
	                                                       << setw(9) << get_itr( LR_WE) << endl;

  // Ratio lapse rate
  if( _set_sf.find(ML_Rt) != _set_sf.end() ) cout << "#mlRt: " << setw(9) << get_surf(ML_Rt)
                                                               << setw(9) << get_rms( ML_Rt)
                                                               << setw(9) << get_itr( ML_Rt) << endl;
  if( _set_sf.find(MP_Rt) != _set_sf.end() ) cout << "#mpRt: " << setw(9) << get_surf(MP_Rt)
                                                               << setw(9) << get_rms( MP_Rt)
                                                               << setw(9) << get_itr( MP_Rt) << endl;
  if( _set_sf.find(LR_Rt) != _set_sf.end() ) cout << "#lrRt: " << setw(9) << get_surf(LR_Rt)
                                                               << setw(9) << get_rms( LR_Rt)
                                                               << setw(9) << get_itr( LR_Rt) << endl;
  if( _set_sf.find(ER_Rt) != _set_sf.end() ) cout << "#erRt: " << setw(9) << get_surf(ER_Rt)
                                                               << setw(9) << get_rms( ER_Rt)
                                                               << setw(9) << get_itr( ER_Rt) << endl;   

  cout << "#end_surface\n";

  return;
}

} // namespace
