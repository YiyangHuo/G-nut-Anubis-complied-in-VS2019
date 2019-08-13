
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

#include "gprod/gprodmet.h"
#include "gutils/gnwmconv.h"

using namespace std;

namespace gnut {
   
// constructor
// ----------
t_gprodmet::t_gprodmet(const t_gtime& t, shared_ptr<t_gobj> pt)
 : t_gprod(t, pt)
{ 
   gtrace("t_gprodmet::t_gprodmet");
   
   id_type(MET);
   id_group(GRP_PRODUCT);
}

// destructor
// --------------
t_gprodmet::~t_gprodmet()
{
}

// get parameter list
// ------------------
set<METEO_ID> t_gprodmet::get_params()
{
  set<METEO_ID> params;
  map<METEO_ID,double>::const_iterator it;

  for(it = _map_val.begin(); it != _map_val.end(); ++it){
    switch( it->first ){
       case METEO_RHUM :
       case METEO_EPRE : params.insert(METEO_RHUM); params.insert(METEO_EPRE); break;
       default : params.insert(it->first); break;
    };
  }
  return params;
}
   
   
   
// add parameter
// -------------
void t_gprodmet::met( const METEO_ID key, const double& val, const double& rms)
{
  switch( key ){
   case METEO_PRES    :   pres(val,rms); return;
   case METEO_EPRE    :   epre(val,rms); return;
   case METEO_RHUM    :   rhum(val,rms); return;
   case METEO_SHUM    :   shum(val,rms); return;
   case METEO_TEMP    :   temp(val,rms); return;
   case METEO_MTEMP   :  mtemp(val,rms); return;
   case METEO_TDEW    :   tdew(val,rms); return;
   case METEO_ZTD     :    ztd(val,rms); return;
   case METEO_ZHD     :    zhd(val,rms); return;
   case METEO_ZWD     :    zwd(val,rms); return;
   case METEO_IWV     :    iwv(val,rms); return;
   case METEO_DT      :  dtemp(val,rms); return;
   case METEO_DM      : dmtemp(val,rms); return;
   case METEO_DE      :  depre(val,rms); return;
   case METEO_DW      :   dzwd(val,rms); return;
   case METEO_DI      :   diwv(val,rms); return;
   case METEO_UNKNOWN : return;
  }	
  return;
}

// get parameter
// -------------
double t_gprodmet::met( const METEO_ID key)
{
  switch( key ){
   case METEO_PRES    : return      pres();
   case METEO_EPRE    : return      epre(); 
   case METEO_RHUM    : return      rhum(); 
   case METEO_SHUM    : return      shum();
   case METEO_TEMP    : return      temp(); 
   case METEO_MTEMP   : return     mtemp();
   case METEO_TDEW    : return      tdew(); 
   case METEO_ZTD     : return       ztd();
   case METEO_ZHD     : return       zhd(); 
   case METEO_ZWD     : return       zwd(); 
   case METEO_IWV     : return       iwv();
   case METEO_DT      : return     dtemp();
   case METEO_DM      : return    dmtemp();
   case METEO_DE      : return     depre();
   case METEO_DW      : return      dzwd();
   case METEO_DI      : return      diwv();
   case METEO_UNKNOWN : return NWM_UNKNOWN;
  }	
  return NWM_UNKNOWN;
}

// get rms of parameter
// --------------------
double t_gprodmet::rms( const METEO_ID key)
{
  switch( key ){
   case METEO_PRES    : return   pres_rms();
   case METEO_EPRE    : return   epre_rms(); 
   case METEO_RHUM    : return   rhum_rms(); 
   case METEO_SHUM    : return   shum_rms();
   case METEO_TEMP    : return   temp_rms(); 
   case METEO_MTEMP   : return  mtemp_rms();
   case METEO_TDEW    : return   tdew_rms(); 
   case METEO_ZTD     : return    ztd_rms(); 
   case METEO_ZHD     : return    zhd_rms(); 
   case METEO_ZWD     : return    zwd_rms(); 
   case METEO_IWV     : return    iwv_rms(); 
   case METEO_DT      : return  dtemp_rms();
   case METEO_DM      : return dmtemp_rms();
   case METEO_DE      : return  depre_rms();
   case METEO_DW      : return   dzwd_rms();
   case METEO_DI      : return   diwv_rms();
   case METEO_UNKNOWN : return  NWM_UNKNOWN;
  }	
  return NWM_UNKNOWN;
}

      
   
// add atmospheric pressure [hPa]
// ------------------------
void t_gprodmet::pres(const double& val, const double& rms)
{
  _gmutex.lock();

  _map_val[METEO_PRES] = val;
  _map_rms[METEO_PRES] = rms;
  _gmutex.unlock();
  return;
}

// get atmospheric pressure [hPa]
// ------------------------
double t_gprodmet::pres()
{
  if( _map_val.find(METEO_PRES) != _map_val.end() ) return _map_val[METEO_PRES];
  return NWM_UNKNOWN;
}

// get rms of atmospheric pressure [hPa]
// -------------------------------
double t_gprodmet::pres_rms()
{
  if( _map_rms.find(METEO_PRES) != _map_rms.end() ) return _map_rms[METEO_PRES];
  return NWM_UNKNOWN;
}


   
   
// add partial water vapour presssure [hPa]
// ----------------------------------
void t_gprodmet::epre(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_EPRE] = val;
  _map_rms[METEO_EPRE] = rms;

  _gmutex.unlock();
  return;
}

// get partial water vapour pressure [hPa]
// ---------------------------------
double t_gprodmet::epre()
{
  if( _map_val.find(METEO_EPRE) != _map_val.end() ) return _map_val[METEO_EPRE];
  if( _map_val.find(METEO_TEMP) != _map_val.end() &&
      _map_val.find(METEO_RHUM) != _map_val.end() ) return _map_val[METEO_RHUM]*esat(_map_val[METEO_TEMP])/100;
  return NWM_UNKNOWN;
}

// get rms of partial water vapour presssure [hPa]
// -----------------------------------------
double t_gprodmet::epre_rms()
{
  if( _map_rms.find(METEO_EPRE) != _map_rms.end() ) return _map_rms[METEO_EPRE];
  return NWM_UNKNOWN;
}



// add relative humidity [%]
// ---------------------
void t_gprodmet::rhum(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_RHUM] = val;
  _map_rms[METEO_RHUM] = rms;

  _gmutex.unlock();
  return;
}

// get relative humidity [%]
// ---------------------
double t_gprodmet::rhum()
{
  if( _map_val.find(METEO_RHUM) != _map_val.end() ) return _map_val[METEO_RHUM];
  if( _map_val.find(METEO_TEMP) != _map_val.end() &&
      _map_val.find(METEO_EPRE) != _map_val.end() ) return 100*_map_val[METEO_EPRE]/esat(_map_val[METEO_TEMP]);
  return NWM_UNKNOWN;
}

// get rms of relative humidity [%]
// ----------------------------
double t_gprodmet::rhum_rms()
{
  if( _map_rms.find(METEO_RHUM) != _map_rms.end() ) return _map_rms[METEO_RHUM];
  return NWM_UNKNOWN;
}


   

// add specific humidity [g/kg]
// ---------------------
void t_gprodmet::shum(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_SHUM] = val;
  _map_rms[METEO_SHUM] = rms;

  _gmutex.unlock();
  return;
}

// get specific humidity [g/kg]
// ---------------------
double t_gprodmet::shum()
{
  if( _map_val.find(METEO_SHUM) != _map_val.end() ) return _map_val[METEO_SHUM];
//if( _map_val.find(METEO_PRES) != _map_val.end() &&
//    _map_val.find(METEO_EPRE) != _map_val.end() ) return 1000.0*e2p(_map_val[METEO_EPRE]*100.0,_map_val[METEO_PRES]*100.0);
  return NWM_UNKNOWN;
}

// get rms of specific humidity [g/kg]
// ----------------------------
double t_gprodmet::shum_rms()
{
  if( _map_rms.find(METEO_SHUM) != _map_rms.end() ) return _map_rms[METEO_SHUM];
  return NWM_UNKNOWN;
}


      
// add dry temperature [K]
// -------------------
void t_gprodmet::temp(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_TEMP] = val;
  _map_rms[METEO_TEMP] = rms;

  if( val < 150 ) _map_val[METEO_TEMP] += TC2TK;      // deg C -> K (check)

  _gmutex.unlock();
  return;
}

// get dry temperature [K]
// -------------------
double t_gprodmet::temp()
{
  if( _map_val.find(METEO_TEMP) != _map_val.end() ) return _map_val[METEO_TEMP];
  return NWM_UNKNOWN;
}

// get rms of dry temperature [K]
// --------------------------
double t_gprodmet::temp_rms()
{
  if( _map_rms.find(METEO_TEMP) != _map_rms.end() ) return _map_rms[METEO_TEMP];
  return NWM_UNKNOWN;
}

   
   

// add mean temperature [K]
// -------------------
void t_gprodmet::mtemp(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_MTEMP] = val;
  _map_rms[METEO_MTEMP] = rms;

  if( val < 150 ) _map_val[METEO_MTEMP] += TC2TK;      // deg C -> K (check)

  _gmutex.unlock();
  return;
}

// get mean temperature [K]
// -------------------
double t_gprodmet::mtemp()
{
  if( _map_val.find(METEO_MTEMP) != _map_val.end() ) return _map_val[METEO_MTEMP];
  return NWM_UNKNOWN;
}

// get rms of mean temperature [K]
// --------------------------
double t_gprodmet::mtemp_rms()
{
  if( _map_rms.find(METEO_MTEMP) != _map_rms.end() ) return _map_rms[METEO_MTEMP];
  return NWM_UNKNOWN;
}



      
// add dew temperature [K]
// -------------------
void t_gprodmet::tdew(const double& val, const double& rms)
{
  _gmutex.lock();

  _map_val[METEO_TDEW] = val;
  _map_rms[METEO_TDEW] = rms;

  if( val < 150 ) _map_val[METEO_TDEW] += TC2TK;      // deg C -> K (check)

  _gmutex.unlock();
  return;
}

// get dew temperature [K]
// -------------------
double t_gprodmet::tdew()
{
  if( _map_val.find(METEO_TDEW) != _map_val.end() ) return _map_val[METEO_TDEW];
  return NWM_UNKNOWN;
}

// get rms of dew temperature [K]
// --------------------------
double t_gprodmet::tdew_rms()
{
  if( _map_rms.find(METEO_TDEW) != _map_rms.end() ) return _map_rms[METEO_TDEW];
  return NWM_UNKNOWN;
}
   

   

// add zenith total delay [mm]
// ----------------------
void t_gprodmet::ztd(const double& val, const double& rms)
{
  _gmutex.lock();
      
  _map_val[METEO_ZTD] = val;
  _map_rms[METEO_ZTD] = rms;
   
  if( val < 3 ) _map_val[METEO_ZTD] *= 1000;      // m -> mm (check)
  if( rms < 3 ) _map_rms[METEO_ZTD] *= 1000;      // m -> mm (check)

  _gmutex.unlock();
  return;
}

// get zenith total delay [mm]
// ----------------------
double t_gprodmet::ztd()
{
  if( _map_val.find(METEO_ZTD) != _map_val.end() ) return _map_val[METEO_ZTD];
  return NWM_UNKNOWN;
}

// get rms of zenith total delay [mm]
// -----------------------------
double t_gprodmet::ztd_rms()
{
  if( _map_rms.find(METEO_ZTD) != _map_rms.end() ) return _map_rms[METEO_ZTD];
  return NWM_UNKNOWN;
}
   

   

// add zenith hydrostatic delay [mm]
// ----------------------------
void t_gprodmet::zhd(const double& val, const double& rms)
{
  _gmutex.lock();

  _map_val[METEO_ZHD] = val;
  _map_rms[METEO_ZHD] = rms;

  if( val < 3 ) _map_val[METEO_ZHD] *= 1000;      // m -> mm (check)
  if( rms < 3 ) _map_rms[METEO_ZHD] *= 1000;      // m -> mm (check)

  _gmutex.unlock();
  return;
}

// get zenith hydrostatic delay [mm]
// ----------------------------
double t_gprodmet::zhd()
{
  if( _map_val.find(METEO_ZHD) != _map_val.end() ) return _map_val[METEO_ZHD];
  return NWM_UNKNOWN;
}

// get rms of zenith hydrostatic delay [mm]
// -----------------------------------
double t_gprodmet::zhd_rms()
{
  if( _map_rms.find(METEO_ZHD) != _map_rms.end() ) return _map_rms[METEO_ZHD];
  return NWM_UNKNOWN;
}



   
// add zenith wet delay [mm]
// --------------------
void t_gprodmet::zwd(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_ZWD] = val;
  _map_rms[METEO_ZWD] = rms;

  if( val < 1 ) _map_val[METEO_ZWD] *= 1000;      // m -> mm (check)
  if( rms < 1 ) _map_rms[METEO_ZWD] *= 1000;      // m -> mm (check)

  _gmutex.unlock();
  return;
}

// get zenith wet delay [mm]
// --------------------
double t_gprodmet::zwd()
{
  if( _map_val.find(METEO_ZWD) != _map_val.end() ) return _map_val[METEO_ZWD];
  return NWM_UNKNOWN;
}

// get rms of zenith wet delay [mm]
// ---------------------------
double t_gprodmet::zwd_rms()
{
  if( _map_rms.find(METEO_ZWD) != _map_rms.end() ) return _map_rms[METEO_ZWD];
  return NWM_UNKNOWN;
}

   
   
   
// add integrated water vapour [kg/m2]
// --------------------
void t_gprodmet::iwv(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_IWV] = val;
  _map_rms[METEO_IWV] = rms;

  _gmutex.unlock();
  return;
}

// get integrated water vapour [kg/m2]
// --------------------
double t_gprodmet::iwv()
{
  if( _map_val.find(METEO_IWV) != _map_val.end() ) return _map_val[METEO_IWV];
  return NWM_UNKNOWN;
}

// get rms of integrated water vapour [kg/m2]
// ---------------------------
double t_gprodmet::iwv_rms()
{
  if( _map_rms.find(METEO_IWV) != _map_rms.end() ) return _map_rms[METEO_IWV];
  return NWM_UNKNOWN;
}

   

   
// add temperature lapse-rate [K/km]
// --------------------
void t_gprodmet::dtemp(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_DT] = val;
  _map_rms[METEO_DT] = rms;

  if( val < 1 ) _map_val[METEO_DT] *= 1000;      // K/m -> K/km (check)
  if( rms < 1 ) _map_rms[METEO_DT] *= 1000;      // K/m -> K/km (check)

  _gmutex.unlock();
  return;
}

// get temperature lapse-rate [K/km]
// --------------------
double t_gprodmet::dtemp()
{
  if( _map_val.find(METEO_DT) != _map_val.end() ) return _map_val[METEO_DT];
  return NWM_UNKNOWN;
}

// get rms of temperature lapse-rate [K/km]
// ---------------------------
double t_gprodmet::dtemp_rms()
{
  if( _map_rms.find(METEO_DT) != _map_rms.end() ) return _map_rms[METEO_DT];
  return NWM_UNKNOWN;
}


   
// add mean temperature lapse-rate [K/km]
// --------------------
void t_gprodmet::dmtemp(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_DM] = val;
  _map_rms[METEO_DM] = rms;

  if( val < 1 ) _map_val[METEO_DM] *= 1000;      // K/m -> K/km (check)
  if( rms < 1 ) _map_rms[METEO_DM] *= 1000;      // K/m -> K/km (check)

  _gmutex.unlock();
  return;
}

// get mean temperature lapse-rate [K/km]
// --------------------
double t_gprodmet::dmtemp()
{
  if( _map_val.find(METEO_DM) != _map_val.end() ) return _map_val[METEO_DM];
  return NWM_UNKNOWN;
}

// get rms of mean temperature lapse-rate [K/km]
// ---------------------------
double t_gprodmet::dmtemp_rms()
{
  if( _map_rms.find(METEO_DM) != _map_rms.end() ) return _map_rms[METEO_DM];
  return NWM_UNKNOWN;
}
   

// add zwd decay-rate [-]
// --------------------
void t_gprodmet::dzwd(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_DW] = val;
  _map_rms[METEO_DW] = rms;

  _gmutex.unlock();
  return;
}

// get zwd decay-rate [-]
// --------------------
double t_gprodmet::dzwd()
{
  if( _map_val.find(METEO_DW) != _map_val.end() ) return _map_val[METEO_DW];
  return NWM_UNKNOWN;
}

// get rms of zwd decay-rate [-]
// ---------------------------
double t_gprodmet::dzwd_rms()
{
  if( _map_rms.find(METEO_DW) != _map_rms.end() ) return _map_rms[METEO_DW];
  return NWM_UNKNOWN;
}


// add iwv decay-rate [-]
// --------------------
void t_gprodmet::diwv(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_DI] = val;
  _map_rms[METEO_DI] = rms;

  _gmutex.unlock();
  return;
}

// get iwv decay-rate [-]
// --------------------
double t_gprodmet::diwv()
{
  if( _map_val.find(METEO_DI) != _map_val.end() ) return _map_val[METEO_DI];
  return NWM_UNKNOWN;
}

// get rms of iwv decay-rate [-]
// ---------------------------
double t_gprodmet::diwv_rms()
{
  if( _map_rms.find(METEO_DI) != _map_rms.end() ) return _map_rms[METEO_DI];
  return NWM_UNKNOWN;
}


// add E decay-rate [-]
// --------------------
void t_gprodmet::depre(const double& val, const double& rms)
{
  _gmutex.lock();
   
  _map_val[METEO_DE] = val;
  _map_rms[METEO_DE] = rms;

  _gmutex.unlock();
  return;
}

// get E decay-rate [-]
// --------------------
double t_gprodmet::depre()
{
  if( _map_val.find(METEO_DE) != _map_val.end() ) return _map_val[METEO_DE];
  return NWM_UNKNOWN;
}

// get rms of E decay-rate [-]
// ---------------------------
double t_gprodmet::depre_rms()
{
  if( _map_rms.find(METEO_DE) != _map_rms.end() ) return _map_rms[METEO_DE];
  return NWM_UNKNOWN;
}


} // namespace
