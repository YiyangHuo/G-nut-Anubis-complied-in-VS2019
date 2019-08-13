
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

#include <cmath>

#include "gprod/gprodtrp.h"

using namespace std;

namespace gnut {
   
// constructor
// ----------
t_gprodtrp::t_gprodtrp(const t_gtime& t, shared_ptr<t_gobj> pt)
 : t_gprodmet(t, pt),
  _grd_N(0.0),
  _grd_E(0.0),
  _grd_N_rms(0.0),
  _grd_E_rms(0.0),
  _lock_ztd(false)
{ 
   gtrace("t_gprodtrp::t_gprodtrp");
   
   _ztd_mf = GMF;
   _grd_mf = CHEN_HERRING;
   
   id_type(TRP);
   id_group(GRP_PRODUCT);
   
   _map_val[METEO_ZHD] = 0.0;  _map_rms[METEO_ZHD] = 0.0;
   _map_val[METEO_ZWD] = 0.0;  _map_rms[METEO_ZWD] = 0.0;
}

// destructor
// --------------
t_gprodtrp::~t_gprodtrp()
{
}


// add zhd
// --------------------------
void t_gprodtrp::zhd(const double& val, const double& rms)
{
   _gmutex.lock();
   if(!_lock_ztd){
     _map_val[METEO_ZHD] = val;
     _map_rms[METEO_ZHD] = rms;
   }
   
   _gmutex.unlock();
   return;
}

// get zhd
// --------------------------
double t_gprodtrp::zhd()
{
  return _map_val[METEO_ZHD];
}

// get zhd rms
// -------------------------
double t_gprodtrp::zhd_rms()
{
  return _map_rms[METEO_ZHD];
}

// add zwd
// --------------------------
void t_gprodtrp::zwd(const double& val, const double& rms)
{
   _gmutex.lock();
   
   if(!_lock_ztd){
     _map_val[METEO_ZWD] = val;
     _map_rms[METEO_ZWD] = rms;
   }
   
   _gmutex.unlock();
   return;
}

// get zwd
// --------------------------
double t_gprodtrp::zwd()
{
  return _map_val[METEO_ZWD];
}

// get zwd rms
// -------------------------
double t_gprodtrp::zwd_rms()
{
  return _map_rms[METEO_ZWD];
}
 
// add ztd
// --------------------------
void t_gprodtrp::ztd(const double& val, const double& rms)
{
   _gmutex.lock();

   if( _map_val.find(METEO_ZHD) == _map_val.end() ||
       _map_val[METEO_ZHD]      == 0.0
   ){
      _lock_ztd = true;  // no more zhd/zwd can be overwritten
      _map_val[METEO_ZHD] = 2.3;
      _map_rms[METEO_ZHD] = 1.0;
   }

   if( _map_val.find(METEO_ZWD) == _map_val.end() ||
       _map_val[METEO_ZWD]      == 0.0
   ){
      _lock_ztd = true;  // no more zhd/zwd can be overwritten
      _map_val[METEO_ZWD] = val - _map_val[METEO_ZHD];
      _map_rms[METEO_ZWD] = rms;
   }else if(_map_rms.find(METEO_ZWD) == _map_rms.end() ||
	    _map_rms[METEO_ZWD]      == 0.0){
      _map_rms[METEO_ZWD] = rms;
   }
	   

   _gmutex.unlock();
   return;
}   
   
// get ztd
// --------------------------
double t_gprodtrp::ztd()
{
  return _map_val[METEO_ZHD] + _map_val[METEO_ZWD];
}

// get zwd rms
// -------------------------
double t_gprodtrp::ztd_rms()
{
  return sqrt(_map_rms[METEO_ZHD]*_map_rms[METEO_ZHD]
	    + _map_rms[METEO_ZWD]*_map_rms[METEO_ZWD]);
}

// add grd
// --------------------------
void t_gprodtrp::grd_N(const double& grd, const double& rms)
{
   _gmutex.lock();
   
  _grd_N     = grd;
  _grd_N_rms = rms;
   
   _gmutex.unlock();
   return;
}

// get grd
// --------------------------
double t_gprodtrp::grd_N() const
{
  return _grd_N;
}

// get grd rms
// -------------------------
double t_gprodtrp::grd_N_rms() const
{
  return _grd_N_rms;
}

// add grd
// --------------------------
void t_gprodtrp::grd(const double& grdN, const double& grdE,
		     const double& rmsN, const double& rmsE )
{
   _gmutex.lock();
   
  _grd_N     = grdN;
  _grd_N_rms = rmsN;
  _grd_E     = grdE;
  _grd_E_rms = rmsE;
   
   _gmutex.unlock();
   return;
}


// add grd
// --------------------------
void t_gprodtrp::grd_E(const double& grd, const double& rms)
{
   _gmutex.lock();
   
  _grd_E     = grd;
  _grd_E_rms = rms;
   
   _gmutex.unlock();
   return;
}

// get grd
// --------------------------
double t_gprodtrp::grd_E() const
{
   return _grd_E;
}

// get grd rms
// -------------------------
double t_gprodtrp::grd_E_rms() const
{
   return _grd_E_rms;
}

// add ZTD mapping function
// ------------------------
void t_gprodtrp::ztd_mf(const ZTDMPFUNC& mf)
{
   _gmutex.lock();
   
   _ztd_mf = mf;
   
   _gmutex.unlock();
   return;
}   
   
// add GRD mapping function
// ------------------------
void t_gprodtrp::grd_mf(const GRDMPFUNC& mf)
{
   _gmutex.lock();
   
   _grd_mf = mf;
   
   _gmutex.unlock();
   return; 
}   
   
} // namespace