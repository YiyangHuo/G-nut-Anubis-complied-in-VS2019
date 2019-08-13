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

*/

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include "gall/gallsurf.h"
#include "gset/gsetgen.h"
#include "gset/gsetnwm.h"
#include "gset/gsetout.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gtropo.h"
#include "gmodels/gfunclin.h"
#include "gmodels/gfitlmm.h"
#include "gproj/gprojlcc.h"
#include "gproj/gprojlci.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallsurf::t_gallsurf() 
 : t_ginterp(),
   _ver("1.0"),
   _overwrite(false),
   _reg_auto(true),
   _min_lat(MIN_LAT),
   _max_lat(MAX_LAT),
   _min_lon(MIN_LON),
   _max_lon(MAX_LON),

   _grid_ll(false),
   _grid_lat(0.0),
   _grid_lon(0.0),
   _dist_lat(0.0),
   _dist_lon(0.0),
   _vert_dst(0.0),

   _vert_temp(dT_LINEAR),
   _vert_zwd(dW_DE),
   _surf_zwd(W_INTEGR),
   _vert_fit(POW),
   _vert_adj(LMM),

   _set(0),
   _out(0),
   _sum(0),
   _grd(0),
   _srf(0),
   _fit(0)
{
  id_type(t_gdata::ALLSURF);
  id_group(t_gdata::GRP_PRODUCT);

  string str_DATA[] = { "TEMP", "PRES", "ZWD", "ZHD" };  // default meteorological parameters
  set<string> param( str_DATA, str_DATA + 4);

  _interp_1d = SPLINE;
  _interp_3d = VER2HOR;
  _interp_ht = SCALE;

  _set_param( param );
}



// constructor
// ----------
t_gallsurf::t_gallsurf(t_gsetbase* gset) // use only if NWM settings available !
 : t_ginterp(),
   _ver("1.0"),
   _overwrite(false),
   _reg_auto(true),

   _grid_lat(0.0),
   _grid_lon(0.0),
   _dist_lat(0.0),
   _dist_lon(0.0),
   _vert_dst(0.0),

   _set(gset),
   _out(0),
   _sum(0),
   _grd(0),
   _srf(0),
   _fit(0)
{
  gtrace("t_gallsurf::construct - gset");

  id_type(t_gdata::ALLSURF);
  id_group(t_gdata::GRP_PRODUCT);

  double minLat = dynamic_cast<t_gsetnwm*>(gset)->min_lat();
  double maxLat = dynamic_cast<t_gsetnwm*>(gset)->max_lat();
  double minLon = dynamic_cast<t_gsetnwm*>(gset)->min_lon();
  double maxLon = dynamic_cast<t_gsetnwm*>(gset)->max_lon();

  set<string> param = dynamic_cast<t_gsetnwm*>(_set)->param();
   
  _set_param( param );
  _set_output();
  _set_region( minLat, maxLat, minLon, maxLon );

  _vert_temp    = dynamic_cast<t_gsetnwm*>(gset)->vert_temp();
  _vert_zwd     = dynamic_cast<t_gsetnwm*>(gset)->vert_zwd(); 
  _surf_zwd     = dynamic_cast<t_gsetnwm*>(gset)->surf_zwd();
  _vert_fit     = dynamic_cast<t_gsetnwm*>(gset)->vert_fit();   
  _vert_adj     = dynamic_cast<t_gsetnwm*>(gset)->vert_adj();

  _interp_1d = str_to_interp_1d( dynamic_cast<t_gsetnwm*>(_set)->interp_time()  );
  _interp_3d = str_to_interp_3d( dynamic_cast<t_gsetnwm*>(_set)->interp_space() );
  _interp_ht = str_to_interp_ht( dynamic_cast<t_gsetnwm*>(_set)->interp_vert() );
}


// constructor
// ----------
t_gallsurf::t_gallsurf(double minLat, double maxLat, double minLon, double maxLon)
 : t_ginterp(),
   _ver("1.0"),
   _overwrite(false),
   _reg_auto(true),

   _grid_lat(0.0),
   _grid_lon(0.0),
   _dist_lat(0.0),
   _dist_lon(0.0),
   _vert_dst(0.0),

   _set(0),
   _out(0),
   _sum(0),
   _grd(0),
   _srf(0),
   _fit(0)
{
  _set_region( minLat, maxLat, minLon, maxLon );
//  if( _min_lat == 0.0 && _max_lat == 0.0 && _min_lon == 0.0 && _max_lon == 0.0 ) _reg_auto = false;

  string str_DATA[] = { "PRES", "TEMP", "ZWD", "ZHD" }; // default meteorological parameters
  set<string> param( str_DATA, str_DATA + 4);
  _set_param( param );
}


// destructor
// ----------
t_gallsurf::~t_gallsurf()
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_map_met::iterator itEPO;
  for( itEPO = _map_met.begin(); itEPO != _map_met.end(); ++itEPO ){
     
    t_map_dat::iterator itGRID;
    for( itGRID = itEPO->second.begin(); itGRID != itEPO->second.end(); ++itGRID ){

#ifdef DEBUG
      cerr << " deleting profile " << fixed << setprecision(3)
           << setw(10) << itGRID->first.crd(0)
           << setw(10) << itGRID->first.crd(1)
	   << itEPO->first.str_ymdhms(" t_gallsurf: ")
	   <<  " " << (itGRID->second)->get_surf(LAT)
	   <<  " " << (itGRID->second)->get_surf(LON)
           << endl;
#endif
    }
    itEPO->second.clear();
  }
  _map_met.clear(); // NOT NECESSARY FOR VALGRIND
   
  if( _grd ){ if( _grd->is_open() ){ _grd->close(); }; delete _grd; }
  if( _srf ){ if( _srf->is_open() ){ _srf->close(); }; delete _srf; }
  if( _fit ){ if( _fit->is_open() ){ _fit->close(); }; delete _fit; }
  if( _out ){ if( _out->is_open() ){ _out->close(); }; delete _out; }
  if( _sum ){ if( _sum->is_open() ){ _sum->close(); }; delete _sum; }
   
  _gmutex.unlock(); return;
}

   
// set projection (only if not set before)
// ----------
int t_gallsurf::gproj(const shared_ptr<t_gproj>& prj)
{  

  if( _proj ){ //  != make_shared<t_gproj>() ){
    if( *_proj == *prj ){ 
      if( _log ) _log->comment(0,"gallsurf","Specific projection already identified in NWM container" );
      return 1;
    }else{
    if( _log ) _log->comment(0,"gallsurf","Different projections identified in NWM container" );
    else               cerr << "gallsurf - Different projections identified in NWM container\n";
    return -1;
    }
  }

  if( _log ) _log->comment(3,"gallsurf","Setting/using specific projection in NWM container" );

  _proj = prj;

  return 1;
}


// add gnwm profile
// ----------
int t_gallsurf::add(shared_ptr<t_gnwmsurf> met, bool proj)
{
  gtrace("t_gallsurf::add");
 
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
   
  _gmutex.lock();

  t_gpair LatLon( met->get_surf(LAT), met->get_surf(LON) );
  t_gpair xy_prj( met->get_surf(YPRJ), met->get_surf(XPRJ) );

  if( LatLon[0] < _min_lat || LatLon[0] > _max_lat ||
      LatLon[1] < _min_lon || LatLon[1] > _max_lon    )
  {
    if( _log && _log->verb() > 1 )
       _log->comment(0,"gallsurf","Lat/Lon: "
					       +dbl2str(LatLon[0])+" /"+dbl2str(LatLon[1])
					       +dbl2str(xy_prj[0])+" /"+dbl2str(xy_prj[1])
					       +"   out of limits: "
					       +dbl2str(_min_lat) +" /"+dbl2str(_max_lat)+" "
					       +dbl2str(_min_lon) +" /"+dbl2str(_max_lon));
    _gmutex.unlock(); return 1;
  }

  // only if projection applied
  if( proj && _proj ){
     
    // SWITCH X/Y to Y/X ---> HERE as LAT/LON (Y/X)
    t_gpair LatLonYX( met->get_surf(YPRJ), met->get_surf(XPRJ) );
#ifdef DEBUG
    cout << "add profile with PROJECTION:  Lat/Lon: " << LatLon[0]   << " "  << LatLon[1]
                                     << "  Y/X_PRJ: " << LatLonYX[0] <<  " " << LatLonYX[1] << endl;
#endif
    LatLon = LatLonYX;
  }

  if( LatLon[0] == NWM_UNKNOWN || LatLon[1] == NWM_UNKNOWN ){
    cerr << "*** warn: met - not defined lat/lon of the profile!" << fixed
         << setw(10) << LatLon[0]
         << setw(10) << LatLon[1]
         << endl;
    _gmutex.unlock(); return -1;
  }

  // profile does not exist or _overwrite (--> delete & new)
  t_gtime epo( met->epoch() );
  if( _map_met[epo].find( LatLon ) == _map_met[epo].end() ){
     
    _map_met[epo][LatLon] = met;
    this->_new_element(met); // allocate in heap (legacy, settings in gallprof used only)
    if( _log && _log->verb() >= 1 ) _log->comment(1,"gallsurf",epo.str_ymdhms("new profile: ")+dbl2str(LatLon[0])+" "+dbl2str(LatLon[1]) );

  }else{

    if( _map_met.lower_bound(epo+60) == _map_met.end() ){ // LAST EPOCH ?
          epo.add_secs(+1); } // plus  1 sec (asume  forwards list of NWM files)
    else{ epo.add_secs(-1); } // minus 1 sec (asume backwards list of NWM files)

    if( _log ) _log->comment(1,"gallsurf","warning - cannot overwrite, epoch shifted: "+epo.str_ymdhms() );
    else               cerr << "gallsurf:  warning - cannot overwrite, epoch shifted: "+epo.str_ymdhms() << endl;     

    _map_met[epo][LatLon] = met;
    this->_new_element(met); // allocate in heap (legacy, settings in gallprof used only)
//  _gmutex.unlock(); return -1;
  }

#ifdef DEBUG
  cout << " adding profile " << fixed << setprecision(3)
       << epo.str_ymdhms(" ")
       << setw(10) << met->get_surf(LAT)
       << setw(10) << met->get_surf(LON)
       << setw(10) << _map_met[epo][LatLon]
       << endl;
#endif
  _gmutex.unlock(); return 0;
}


// add param
// ----------
int t_gallsurf::add_param(MET_DATA par)
{
  gtrace("t_gallsurf::add_param(MET_DATA)");
  _gmutex.lock();

  bool found = false;
  vector<MET_DATA>::const_iterator it;
  for( it = _vec_par.begin(); it < _vec_par.end(); ++it )
    if( *it == par ) found = true;

  if( !found ) _vec_par.push_back(par);

  _gmutex.unlock(); return 0;
}

   
// add param
// ----------
int t_gallsurf::add_param(MET_SURF par)
{
  gtrace("t_gallsurf::add_param(MET_SURF)");
  _gmutex.lock();

  bool found = false;
  vector<MET_SURF>::const_iterator it;
  for( it = _vec_srf.begin(); it < _vec_srf.end(); ++it )
    if( *it == par ) found = true;

  if( !found ) _vec_srf.push_back(par);
  
  _gmutex.unlock(); return 0;
}


   
// get all epochs
// ----------
set<t_gtime> t_gallsurf::epochs()
{
  gtrace("t_gallsurf::epochs");
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  set<t_gtime> s_epochs;
   
  t_map_met::const_iterator itEPO = _map_met.begin();
  while( itEPO != _map_met.end() ){
    s_epochs.insert( itEPO->first );
    ++itEPO;
  }

  _gmutex.unlock();
  return s_epochs;
}


// get all pairs
// ----------
set<t_gpair> t_gallsurf::pairs( const t_gtime& epo )
{
  gtrace("t_gallsurf::pairs");
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  set<t_gpair> s_pairs;

  // not correct time
  if( _map_met.find(epo) == _map_met.end() ){
    _gmutex.unlock(); return s_pairs;
  }
     
  t_map_dat::const_iterator itGRID = _map_met[epo].begin();
  while( itGRID != _map_met[epo].end() ){
    s_pairs.insert( itGRID->first );
    ++itGRID;
  }

  _gmutex.unlock();
  return s_pairs;
}


// -----------------------------
// dLat, dLon grid step              ---->>>> POZOR, NEFUNGUJE SPRAVNE PRO VYBER STANIC!!!!!
// -----------------------------
int t_gallsurf::_get_grid_step()
{
   gtrace("t_gallsurf::get_grid_step");

   if( ! double_eq(_grid_lat, 0.0) ) return 0;
   if( ! double_eq(_grid_lon, 0.0) ) return 0;
     
   double Lon = NWM_UNKNOWN;
   double Lat = NWM_UNKNOWN;
   double StepLat = 9999;
   double StepLon = 9999;

   for(t_map_met::iterator iEP = _map_met.begin(); iEP != _map_met.end(); ++iEP)
   {
     t_map_dat dat = iEP->second;
     for(t_map_dat::iterator itD = dat.begin(); itD != dat.end(); ++itD)
     {
       if( Lat == NWM_UNKNOWN ) Lat = itD->first.crd(0);	
       if( Lon == NWM_UNKNOWN ) Lon = itD->first.crd(1);

       double dLat = fabs(itD->first.crd(0) - Lat);
       double dLon = fabs(itD->first.crd(1) - Lon);

       if( dLat > 0.0 && dLat < StepLat ) StepLat = fabs(itD->first.crd(0)-Lat);
       if( dLon > 0.0 && dLon < StepLon ) StepLon = fabs(itD->first.crd(1)-Lon);
     }
   }
   
#ifdef DEBUG
   cout << " dLat dLon  ... " <<   StepLat   << " : " <<   StepLon   << endl;
#endif
   _grid_lat = StepLat;
   _grid_lon = StepLon;

   return 0;
}


// -------------------------------------
// Get epochs for interpolation (single epoch)
// -------------------------------------
int t_gallsurf::_get_node_epochs(const t_gtime& epo, set<t_gtime>& set_dt)
{
   gtrace("t_gallsurf::get_node_epochs (single epoch)");
   
   t_map_met::iterator itEP, itEPforw, itEPback;
   itEP = _map_met.lower_bound(epo);
     
   // exact time (1 value)
   if( itEP->first == epo ){
      set_dt.insert(epo);
//      cerr << "EPOCHS -> " << epo.str_ymdhms() <<  " " << itEP->first.str_ymdhms() << endl;
//      if( _log ) _log->comment(0,"gallsurf","Warning - single epoch found: " + epo.str_ymdhms() + " " + int2str(set_dt.size()) );
//      else               cerr << "gallsurf:  Warning - single epoch found: " <<epo.str_ymdhms() <<" " <<set_dt.size() << endl;      
      return 0; 
   }
 
   itEPforw = itEPback = itEP;
   if( itEP == _map_met.begin() || itEP == _map_met.end() ){
     if( _log ) _log->comment(0,"gnwmsurf","no data for temporal interpolation");
     return -1;
   }else{
    // linear (2 values)
    --itEPback;
    set_dt.insert(itEPforw->first);
    set_dt.insert(itEPback->first);
    ++itEPforw;
   }
   
   if( itEPback == _map_met.begin() || itEPforw == _map_met.end() || _interp_1d != SPLINE ){
    --itEPforw;
   }else{
     // spline (4 values) 
     --itEPback;
     set_dt.insert(itEPforw->first);
     set_dt.insert(itEPback->first);
   }
#ifdef DEBUG
   cout << " node epochs found = "
        << " " << interp_1d_to_str( _interp_1d )
        << " " << interp_2d_to_str( _interp_2d )
        << " " << interp_3d_to_str( _interp_3d )
        << " " << interp_ht_to_str( _interp_ht )
        << " " << epo.str_ymdhms()
        << itEPback->first.str_ymdhms(" beg: ")
        << itEPforw->first.str_ymdhms(" end: ")
        << " " << set_dt.size() << endl;
#endif
   
#ifdef DEBUG
   cout << " epochs : ";
   for(itEP = _map_met.begin(); itEP != _map_met.end(); ++itEP ) cout << itEP->first.str_ymdhms(" ");
   cout << endl;
#endif

  return 0;
}


// -------------------------------------
// Get epochs for interpolation (interval)
// ------------------------------------- 
int t_gallsurf::_get_node_epochs(const t_gtime& beg, const t_gtime& end, set<t_gtime>& set_dt)
{
  gtrace("t_gallsurf::get_node_epochs (interval)");

#ifdef DEBUG
  cerr << " map_met epoch:"
       << beg.str_ymdhms(" ")
       << end.str_ymdhms(" ") << endl;
#endif

  // interface to a single epoch
  if( beg == end ) return _get_node_epochs(beg, set_dt);
     
  t_map_met::iterator itEPO;
  for( itEPO = _map_met.begin() ; itEPO != _map_met.end(); ++itEPO )
  {
//  cerr << " map_met epoch : " << itEPO->first.str_ymdhms() << endl;
    if( beg        == itEPO->first  ||  itEPO->first == end      ||
      ( beg-6*3600  < itEPO->first  &&  itEPO->first  < end + 6*3600) ) // 6h reserve
    {
      set_dt.insert( itEPO->first );
       
#ifdef DEBUG
      cout << " node epochs found = " 
           << " " << interp_1d_to_str( _interp_1d )
           << " " << interp_2d_to_str( _interp_2d )
           << " " << interp_3d_to_str( _interp_3d )
           << " " << interp_ht_to_str( _interp_ht )
	   << itEPO->first.str_ymdhms(" epo: ")
           << beg.str_ymdhms(" beg: ")
           << end.str_ymdhms(" end: ")
           << " " << set_dt.size() << endl;
#endif
    }
  }   

  if( set_dt.size() == 0 ) return 1;

  if( _map_met.size() < 2 && beg != end ){
    if( _log ) _log->comment(0,"gallsurf","Warning - single epoch found. No temporal interpolation possible!" );
    else               cerr << "gallsurf:  Warning - single epoch found. No temporal interpolation possible!\n";
    return -1;
  }

  return 0;
}

// ---------------------------------------------
// get data points
// ---------------------------------------------
int t_gallsurf::_get_grid_points( const t_gpair& pt,
				  const set<t_gtime>& set_dt,
			                set<t_gpair>& set_pt )
{
  gtrace("t_gallsurf::get_grid_points");

  if( set_dt.size() <= 0   ||
    ( double_eq(_grid_lat, 0.0) ||
      double_eq(_grid_lon, 0.0) )){ set_pt.clear(); return -1; }

  // without cache
  unsigned int nPT = 0;
  t_gtime tt( *(set_dt.begin()) );

  t_map_met::iterator itMET = _map_met.find(tt);

  t_map_dat data = itMET->second;
  t_map_dat::reverse_iterator itDAT = data.rbegin(); // combine forw/reverse iterator
   
  for( ; itDAT != data.rend(); ++itDAT)
  {     
    t_gpair node = itDAT->first;

#ifdef DEBUG
   if( pt[0] - node[0] < _grid_lat )
    cout << " test = " <<   pt.crd(0)  << " " <<   pt.crd(1)
         <<        " " << node.crd(0)  << " " << node.crd(1)
	 << " grd_lat: " << _grid_lat 
	 << " grd_lon: " << _grid_lon
	 << " grd_ll: "  << _grid_ll
         << " proj: "    << _proj
         << endl;
#endif
   
    // 1 point
    if( pt == node )
    {     
//    _grid_points.insert(node);     // cache
      set_pt.insert(node);
      nPT = 1;  break;       
    }

    bool okLat(false), okLon(false);     
    if( fabs( pt[0] - node[0] ) < _grid_lat ) okLat = true;
    if( fabs( pt[1] - node[1] ) < _grid_lon ) okLon = true;
     
    // cyclic boundary (LatLon only)
    if( _grid_ll &&  ( fabs( pt[1] - node[1]-360.0 ) < _grid_lon ||
                       fabs( pt[1] - node[1]+360.0 ) < _grid_lon ) ) okLon = true;

    // 2 points
    if( ( pt[0] == node[0] && okLon ) ||
	( pt[1] == node[1] && okLat ) )
    {
//    _grid_points.insert(node);     // cache
      set_pt.insert(node);
      nPT = 2;
      if( set_pt.size() == 2 )  break;       
      continue; 
    }

    if( okLat && okLon )
    {	   
//    _grid_points.insert(node);     // cache
      set_pt.insert(node);
      nPT = 4;
      if( set_pt.size() == 4 )  break;
      continue; 
    }
  }
    
  if( nPT == 0 || nPT != set_pt.size() ){
    set_pt.clear();
     
    if( _log ){
      ostringstream os;
      os << fixed << setprecision(3)
	 << "No grid points " << set_pt.size()
	 << " found for interpolation: " << setw(8) << pt.crd(0) << " " << setw(8) << pt.crd(1)
	 << " .. skip user point";
      _log->comment(0,"gallsurf",os.str());
    }
    return -1;
  }

#ifdef DEBUG     
  if( _log && _log->verb() >= 1 ){
    for( set<t_gpair>::iterator it = set_pt.begin(); it != set_pt.end(); ++it )
    {	  
      ostringstream os;
      os << fixed << setprecision(3)
         << " user: " << setw(8)  <<  pt.crd(0) << " " << setw(8) <<  pt.crd(1)
         << " grid: " << setw(8)  << it->crd(0) << " " << setw(8) << it->crd(1)
         << " " << set_pt.size();
      _log->comment(1,"gallsurf",os.str());
    }
  }
#endif
   
  if( _log && _log->verb() >= 1 ){
    ostringstream os;
    os <<  "inc_lat: " << _grid_lat
       << " inc_lon: " << _grid_lon
       << " lat/lon: " << pt[0] << " " << pt[1]
       << " pt_size: " << set_pt.size()
       << " ep_size: " << set_dt.size()
       << " tot_epo: " << _map_met.size();
   _log->comment(1,"gallsurf", os.str());
  }

  return 0;
}


// ---------------------------------------------
// Create maps of GRID/REQ.PT for all NWM epochs
// ---------------------------------------------
int t_gallsurf::_get_grid_interp(const t_gtriple& pt,
			         const set<t_gtime>& set_epo,
			         const set<t_gpair>& set_pt,
				 const vector<MET_SURF> vec_surf,
			         t_map_met& map_grid,          // map of NWM epochs/pairs/profiles
			         t_map_met& map_prof)          // map of NWM epochs/reqpt/profiles req.point)
{
  gtrace("t_gallsurf::_get_grid_interp");

  // loop over NWM epochs (not requested epochs)
  set<t_gtime>::iterator itEPO;
  for( itEPO = set_epo.begin(); itEPO != set_epo.end(); ++itEPO )
  {  
     
//  cout << "interpolate: NWM Epoch = " << itEPO->str_ymdhms()  << endl;
     
    t_gtime tt( *itEPO );
    t_map_dat::iterator itDAT;
    t_map_met::iterator itMET = _map_met.find( tt );
    if( itMET == _map_met.end() ){ cerr << "# --> warning: should not occure [interpolate]!\n"; continue; }
    int irc = 0;

    // GRID POINTS (all epochs): create map of gnwmsurf pointers (epochs/grid_pt)
    set<t_gpair>::iterator itPT;
    for( itPT = set_pt.begin(); itPT != set_pt.end(); ++itPT )
    {

      if( (itDAT = itMET->second.find( *itPT )) != itMET->second.end() ){
             map_grid[tt][itDAT->first] = itDAT->second;
      }else{ cerr << "error: missing point for bilinear interpolation !\n";
   	     map_grid.erase(*itEPO);
	     irc++; continue;
      }
    }

    if( irc > 0 ) continue; // skip problematic

    // REQ.POINT (all epochs): create map of NEW gnwmsurf pointers (epochs/req_pt)
    map_prof[tt][pt.gpair()] = make_shared<t_gnwmsurf>();
    map_prof[tt][pt.gpair()]->epoch(tt);
    map_prof[tt][pt.gpair()]->add_surf(LAT, pt[0]);
    map_prof[tt][pt.gpair()]->add_surf(LON, pt[1]);
    map_prof[tt][pt.gpair()]->add_surf(HEL, pt[2]);
    map_prof[tt][pt.gpair()]->calc_surface();
//  cerr << tt.str_ymdhms("epoch: ") << " " << pt[0] << " " << pt[1] << " " << pt[2] << endl;

    // interpolate basic parameters (surface)
    vector<MET_SURF>::const_iterator itSURF;
    for( itSURF = vec_surf.begin(); itSURF != vec_surf.end(); ++itSURF )
    {
      // PREPARE DATA for BILINEAR: interpolate basic surface parameters
      map<t_gpair, double> map_space;
      map<t_gtime, double> map_gtime;
      for( itPT = set_pt.begin(); itPT != set_pt.end(); ++itPT ){
        map_space[ *itPT ] = itMET->second[ *itPT ]->get_surf( *itSURF );
	 
	// COPY PROFILE SETTINGS (probably not necessary, but better do)
        map_prof[tt][pt.gpair()]->decay_min( itMET->second[ *itPT ]->decay_min() );
        map_prof[tt][pt.gpair()]->decay_max( itMET->second[ *itPT ]->decay_max() );
        map_prof[tt][pt.gpair()]->vert_temp( itMET->second[ *itPT ]->vert_temp() );
        map_prof[tt][pt.gpair()]->vert_zwd(  itMET->second[ *itPT ]->vert_zwd() );
        map_prof[tt][pt.gpair()]->surf_zwd(  itMET->second[ *itPT ]->surf_zwd() );

#ifdef DEBUG
        cout << "interp_aux1: "  << tt.str_ymdhms() << fixed << setprecision(5)
	     << setw(8) << t_gnwmbase::met_surf_str( *itSURF ) << " " << (*itPT)[0] << " " << (*itPT)[1]
	     << "  input: " << map_space[ *itPT ] << endl;
#endif
      }
      double res = NWM_UNKNOWN;
      if( bilinear( map_space, pt.gpair(), res ) == 0 ){ 
	map_prof[tt][pt.gpair()]->add_surf(*itSURF, res);
	 
	// save also basic parameters
//        if( m_surf ) (*m_surf)[tt][*itSURF] = res; // OBSOLETE !!!

#ifdef DEBUG	 
        cout << "interp_aux2: " << tt.str_ymdhms() << fixed << setprecision(5)
	     << setw(8) << map_prof[tt][pt.gpair()]->met_surf_str( *itSURF ) << " " << pt[0] << " " << pt[1]
	     << " output: " << res 
	     << endl;
#endif
      }
      map_space.clear();
    }
  }
  return 0;
}

// temporal interpolation (linear/spline)
// -----------------
int t_gallsurf::_interp_temporal(const t_gtime& beg, 
			 	 const t_gtime& end,
	 	                 const double& sampling,
	  		 	 const map<t_gtime,double>& inp,
	 	                       map<t_gtime,double>& out)
{ 
  gtrace("t_gallsurf::_interp_temporal");

  map<t_gtime, double> map_gtime;   // local map
  map<t_gtime, double>::const_iterator itEPO,itE;

  t_gtime epo(beg);
  t_gtime epo_BEF0(FIRST_TIME);
  t_gtime epo_BEF( FIRST_TIME);
  t_gtime epo_AFT( FIRST_TIME);
  t_gtime epo_AFT1(FIRST_TIME);

  // loop over NWM epochs (shifting epochs map by 1)
  for( itEPO = inp.begin(); itEPO != inp.end(); ++itEPO )
  {
    t_gtime tt = itEPO->first; //  cerr << "tt = " << tt.str_ymdhms() << endl;
     
    // add current epoch (if not available)
    if( map_gtime.find(tt) == map_gtime.end() ) map_gtime[tt] = itEPO->second;
     
    // SINGLE EPOCH INTERPOLATION
    if( tt == beg && tt == end ){ out[tt] = inp.at(tt); return 0;}

    // skip the first epoch
    if( epo_BEF == FIRST_TIME ){ epo_BEF = tt; continue; }

    epo_AFT = tt;

    // support SPLINE if data available at epoBEF0 & epo_AFT1
    if( _interp_1d == SPLINE && map_gtime.size() > 2 ){
      itE = itEPO;

      // more data after available
      if( ++itE != inp.end() ){
         epo_AFT1 = itE->first;
	 double tmp = inp.at(epo_AFT1);
         map_gtime[epo_AFT1] = tmp;
//       cout << epo_AFT1.str_ymdhms("adding epoch AFT1: ") << endl;
      }
      // continue with interpolation
      else if( map_gtime.find(epo_BEF0) != map_gtime.end() ){ 
        map_gtime.erase(epo_BEF0); 
//      cout << epo_BEF0.str_ymdhms("removing epoch BEF0: ") << endl;
      }
       
      // check model forecast change (i.e. dT<10s) to avoid using spline in model-boundary and use linear instead
      t_gtime tst = FIRST_TIME ;
      for( itE = map_gtime.begin(); itE != map_gtime.end(); ++itE ){
	 
        if( fabs(itE->first - tst) > 10 ){ tst = itE->first; continue; } // OK for SPLINE (assume the same NWM model)

	if( _log ) _log->comment(1,"gallsurf","NWM change identified (force LINEAR interpolation) "
				                   +tst.str_ymdhms()+" "+itE->first.str_ymdhms());
	map_gtime.erase(map_gtime.begin()->first);
	map_gtime.erase(epo_AFT1);
	break;                                      // OK for LINEAR only (don't allow SPLINE)
      }
    }

    // SHOULD NOT HAPPEN!
    if( map_gtime.size()%2 != 0 )
      cerr << tt.str_ymdhms("Warn - not correct map size: ") << " " << map_gtime.size() << " !\n";
     
#ifdef DEBUG
    cout << fixed << setprecision(4) << tt.str_ymdhms("#EPO: ") << " " << map_gtime.size() << ": ";
    for( map<t_gtime, double>::iterator it = map_gtime.begin(); it != map_gtime.end(); ++it)
      cout << " ... " << it->first.str_ymdhms();
    cout << endl;
#endif

    // sync begin if necessray
    if( epo < epo_BEF ) epo = epo_BEF;

    // interpolate btw 2 NWP epochs
    while( epo == epo_BEF || epo == epo_AFT || ( epo > epo_BEF && epo < epo_AFT ) )
    {
      double res = 0.0;
      if( map_gtime.size() == 2 && linear( map_gtime, epo, res ) == 0 ){ out[epo] = res; }
      if( map_gtime.size() == 4 && spline( map_gtime, epo, res ) == 0 ){ out[epo] = res; }
//      cout << epo.str_ymdhms("#EPO: ") << setw(12) << res << endl;
      epo = epo + sampling;
    }

    // shift by single NWM (source) epoch
    if( map_gtime.size() == 4 && _interp_1d == SPLINE ){ map_gtime.erase(epo_BEF0); }
    if( map_gtime.size() == 2 && _interp_1d == LINEAR ){ map_gtime.erase(epo_BEF ); }

    epo_BEF0 = epo_BEF;
    epo_BEF  = tt;

  } // loop over NWM epochs
  return 0;
}


// point interpolation - spatial + temporal
// -----------------
int t_gallsurf::interpolate_point( const t_gtriple& ell,
 	                           const t_gtime& beg,
			           const t_gtime& end,
			           map<t_gtime, t_map_metdata>          *m_data,
				   map<t_gtime, t_map_metsurf>          *m_surf,
				   map<t_gtime, shared_ptr<t_gnwmsurf>> *m_prof,
			           double sampling )
{ 
  gtrace("t_gallsurf::interpolate_point");

  t_gtriple yxh(ell);

  // USE THE PROJECTION (and switch the system to GRID (X=Lat, Y=Lon, but projetion: X=Lon, Y=Lat)
  if( _proj ){ t_gtriple xxx = _proj->ll2proj( ell ); yxh[0]=xxx[1]; yxh[1]=xxx[0]; yxh[2]=xxx[2]; }
  t_gpair req_pt = yxh.gpair();   

  if( _log && _log->verb() >= 1 ){
     _log->comment(1,"gallsurf","Req crd[ELL]: "+dbl2str(ell[0])+" "+dbl2str(ell[1])+" "+dbl2str(ell[2]));
     _log->comment(1,"gallsurf","Req crd[YXH]: "+dbl2str(yxh[0])+" "+dbl2str(yxh[1])+" "+dbl2str(yxh[2]));
  }

  _vert_dst = 0.0;

  if( m_data == 0 && m_surf == 0 && m_prof == 0 )  return 1;
  if( sampling <= 0 ) sampling = DEFAULT_SAMPLING;

  set<MET_SURF>    set_surf;
  vector<MET_SURF> vec_surf;  vector<MET_SURF>::iterator itSURF;
  vector<MET_DATA> vec_par;   vector<MET_DATA>::iterator itPAR;
   
  // user-requested parameters (from _vec_surf)
  for( itSURF = _vec_srf.begin(); itSURF != _vec_srf.end(); ++itSURF ){
    vec_surf.push_back(*itSURF);
    set_surf.insert(*itSURF);    // just to  list
  }

  // mandatory background surface parameters (for vert_scaling)
//if( set_surf.find(TSURF) == set_surf.end() ) vec_surf.push_back(TSURF);
//if( set_surf.find(ESURF) == set_surf.end() ) vec_surf.push_back(ESURF);
  if( set_surf.find(dT)    == set_surf.end() ) vec_surf.push_back(dT);
  if( set_surf.find(dM)    == set_surf.end() ) vec_surf.push_back(dM);
  if( set_surf.find(dE)    == set_surf.end() ) vec_surf.push_back(dE);
  if( set_surf.find(dW)    == set_surf.end() ) vec_surf.push_back(dW);
  if( set_surf.find(dI)    == set_surf.end() ) vec_surf.push_back(dI);

//if( set_surf.find(dR)    == set_surf.end() ) vec_surf.push_back(dR);  // DO NOT SUPPORT NOW
  if( _interp_3d           == HOR2VER        ) vec_surf.push_back(HEL); // interp. HEL of req.pt: 1.space 2.vert 3.time

  set<t_gpair> set_pt;         set<t_gpair>::iterator itPT;             // source grid pointers
  set<t_gtime> set_epo;        set<t_gtime>::iterator itEPO;            // source epochs

  t_map_met map_prof;          t_map_met::iterator itPROF;              // map (profiles) for interpolated point (fixed)
  t_map_met map_grid;          t_map_met::iterator itMET;               // map (profiles) for NWM grid points (BILINEAR:4)
                               t_map_dat::iterator itDAT;

  map<t_gpair, double> map_space;  
  map<t_gtime, double> map_gtime;
  map<t_gtime, double> map_results;
  map<t_gtime, double>::const_iterator itE;

  int irc = 0;
  irc += _get_node_epochs(beg, end, set_epo);                     // get epoch nodes (1,2,4)
  irc += _get_grid_step();                                        // get grid step (if not set)  --> NEFUNGUJE SPRAVNE PRO VYBER STANIC!!
  irc += _get_grid_points(req_pt, set_epo, set_pt);               // get nearby grid points (1,2,4)

#ifdef DEBUG    // INTERPOLATING EPOCHS
  for( itEPO = set_epo.begin(); itEPO != set_epo.end(); ++itEPO )
    cerr << (*itEPO).str_ymdhms("interp_epochs: ") 
         << beg.str_ymdhms(" beg:") << end.str_ymdhms(" end:") << endl;
#endif

  if( irc != 0 ){
     if( _log ) _log->comment(0,"gallsurf","Warning: initialization problem [interpolate]");
     else               cerr << "gallsurf - Warning: initialization problem [interpolate]\n";
     return -1;
  }

  // fill maps of epochs/data for GRID(4 x nEpo) & REQ_PT(1 x nEpo)
  irc += _get_grid_interp(yxh, set_epo, set_pt, vec_surf, map_grid, map_prof);

  // temporal for basic parameters (only if requested --> if( m_surf || m_prof )
  if( m_surf || m_prof ){
    map_gtime.clear();
    for( itSURF = vec_surf.begin(); itSURF != vec_surf.end(); ++itSURF )
    {
#ifdef DEBUG
    cerr << "interpolate surf: " << t_gnwmbase::met_surf_str(*itSURF,false) << endl;
#endif
      for( itEPO = set_epo.begin(); itEPO != set_epo.end(); ++itEPO ){
        map_gtime[*itEPO] = map_prof[*itEPO][yxh.gpair()]->get_surf(*itSURF);
//        cerr << " map_gtime:  " << map_gtime.size()  << " " << itEPO->str_ymdhms() << endl; cerr.flush();
      }

      _interp_temporal( beg, end, sampling, map_gtime, map_results );

      for( itE = map_results.begin(); itE != map_results.end(); ++itE )
      { t_gtime tt = itE->first;
        if( ! (tt<beg || tt>end) ) (*m_surf)[tt][*itSURF] = map_results[tt]; // cut requested period
      }

      map_gtime.clear();
    }
  }

  // loop over PARAMETERS (LEVEL data)
  // =================================
  for( itPAR = _vec_par.begin(); itPAR != _vec_par.end(); ++itPAR )
  {
#ifdef DEBUG
    cout << "interpolate param: " << t_gnwmbase::met_data_str(*itPAR,false) << endl;
#endif
    map_gtime.clear();

    if( *itPAR == ZTD ){ continue; }

    // loop over NWM original epochs
    for( itMET = map_grid.begin(); itMET != map_grid.end(); ++itMET )
    {
      t_gtime tt(itMET->first);
	     
      // loop over grid points
      for( itDAT  = itMET->second.begin();
	   itDAT != itMET->second.end() && *itPAR != ZTD; 
	 ++itDAT )
      {
	shared_ptr<t_gnwmsurf> gsurf = itDAT->second; // profile data
	double val = gsurf->get_surfdata(*itPAR);     // individual data
	double hel = yxh[2];                          // reference height
	
	if( fabs(hel - gsurf->get_surf(HEL)) > _vert_dst ) _vert_dst = fabs(hel - gsurf->get_surf(HEL));
		   
	if( _interp_3d == VER2HOR ){                                                         // 1.vert, 2.horiz, 3.time
	  if(      _interp_ht == SCALE ){       val = gsurf->get_in_hel( *itPAR, yxh[2] ); } //  .. approximate vert.profile
	  else if( _interp_ht == INTERPOLATE ){ val = gsurf->get_interp( *itPAR, yxh[2] ); } //  .. interpolate vert.profile
        }else if( _interp_3d == HOR2VER ){      hel = map_prof[tt][req_pt]->get_surf(HEL);   // 1.horiz, 2.vert, 3.time
	                                        val = gsurf->get_in_hel( *itPAR, hel );    } //  .. mean grid altitude
		  
	// add corrected/uncorrected VALUE for spatial interpolation
	map_space[itDAT->first] = val;
		  
	if( _log && _log->verb() > 2 ){
          ostringstream os;
          os << fixed    << setprecision(2)
	     << "#GRD: " << tt.str_ymdhms()
             <<  setw(7) << itDAT->first.crd(0)
             <<  setw(7) << itDAT->first.crd(1)
		         << setprecision(3)
             <<  setw(9) << itDAT->second->met_data_id(  *itPAR )
             <<  setw(9) << itDAT->second->get_surfdata( *itPAR );    // orig value
          if(!double_eq( val, itDAT->second->get_surfdata( *itPAR )) ){ 
	    os << " -> " << setw(9) << val                          // new  value
	                 << setprecision(3)
	       << "    " << setw(9) << gsurf->get_surf( HEL )       // from height
	       << " -> " << setw(9) << hel;                         // to   height
	  }
          _log->comment(3,"gallsurf",os.str());
         }

      }      
      // spatial interpolation for (un)corrected values
      double res = 0.0;
      if( bilinear( map_space, req_pt, res ) == 0 ){
	map_prof[tt][req_pt]->add_surf(    *itPAR, res); // original data (TEMP, EPRE, ..)
	map_prof[tt][req_pt]->add_surfdata(*itPAR, res); // surface  data (TSURF<=>TEMP, ESURF<=>EPRE, .. )

        if( _interp_3d == HOR2VER ){                                // 1.horiz, 2.vert, 3.time
          res = map_prof[tt][req_pt]->get_in_hel( *itPAR, yxh[2] ); // cout << "#[HOR2VER] .. VERT SCALE \n";
        }	 
        map_gtime[tt] = res;

        if( _log && _log->verb() > 2 ){
	  ostringstream os;
          os << fixed    << setprecision(4)
             << "#INT: " << tt.str_ymdhms()
	     << " "      << yxh[0]
	     << " "      << yxh[1]
	     << " "      << yxh[2]
             << " "      << map_prof[tt][req_pt]->met_data_id(  *itPAR )
	     << " "      << map_prof[tt][req_pt]->get_surfdata( *itPAR );
          if(!double_eq( res, map_prof[tt][req_pt]->get_surfdata( *itPAR )) ){
  	  os << " -> " << res
	     << "    " << map_prof[tt][req_pt]->get_surf( HEL )
	     << " -> " << yxh[2];
	  }
          os << " "      << map_gtime.size();
          _log->comment(3,"gallsurf",os.str());
	}

      }else{
        if(_log) _log->comment(0,"gallsurf", "Warning - bilinear interpolation failed");
      }
      map_space.clear();
    }

    // temporal interpolation
    _interp_temporal( beg, end, sampling, map_gtime, map_results );
     
    for( itE = map_results.begin(); itE != map_results.end(); ++itE )
    { 
      t_gtime tt = itE->first;
      if( ! (tt<beg || tt>end) )
	(*m_data)[tt][*itPAR] = map_results[tt]; // cut requested period

    }     

    map_gtime.clear();

  } // loop over PARAMETERS

  // HANDLE ZTD (cummulative parameter)
  for( itPAR = _vec_par.begin(); itPAR != _vec_par.end(); ++itPAR ){
    if( *itPAR == ZTD ){
      for( itE = map_results.begin(); itE != map_results.end(); ++itE ){
        t_gtime tt = itE->first;
        if( ! (tt<beg || tt>end) )
	   (*m_data)[tt][*itPAR] = (*m_data)[tt][ZHD] + (*m_data)[tt][ZWD]; // cummulate + cut ZTD
      }
    }
  }

  // clean local data
  for( itPROF = map_prof.begin(); itPROF != map_prof.end(); ++itPROF ){
    t_gtime tt = itPROF->first;
    itDAT = itPROF->second.find(req_pt);
    if( m_prof && itDAT != itPROF->second.end() && ! (tt<beg || tt>end) ){ 
      (*m_prof)[tt] = itPROF->second[req_pt];
    }
  }
  map_prof.clear();

  return 0;
}


// ---------------------------------------------
// estimate gradients
// ---------------------------------------------
int t_gallsurf::analyse_gradient( const t_gtriple& ell,                         // ELL:BLH of point of interest
				  const double& cut_off,                        // [deg] elevation cut-off
                                  const double sampling,                        // [sec] sampling interval
				  map<t_gtime, map<MET_DATA, double> >& rf_data)// ref.point interpolated data
{   
  gtrace("t_gallsurf::analyse_gradient");

  map<t_gtime, map<MET_DATA, double> >::const_iterator itEPO, itXXX;

  map<t_gtime,double> inp_egrd, inp_ngrd;
  map<t_gtime,double> out_egrd, out_ngrd;
   
#ifdef DEBUG  // WRITE INPUT DATA (ZHD, ZWD, ZTD)
  if( _proj ) cout << "USING THE LCC PROJECTION ! <---- \n";
  cout << "REQ POINT: " << fixed << setprecision(3) << ell[0] << " " << ell[1] << " " << ell[2] << endl;
  for( itEPO = rf_data.begin(); itEPO != rf_data.end(); ++itEPO)
  {
    cout << itEPO->first.str_ymdhms() << " ";
     
    map<MET_DATA, double>::const_iterator itDAT;
    for(itDAT = itEPO->second.begin(); itDAT != itEPO->second.end(); ++itDAT)
       cout << "  " << t_gnwmbase::met_data_id(itDAT->first) << ": " << itDAT->second;
     
    cout << endl;
  }
#endif

  map<t_gpair, double> plane_datH;  map<t_gpair, double> plane_datW;  map<t_gpair, double> plane_datT;
  map<t_gpair, double> plane_eleH;  map<t_gpair, double> plane_eleW;  map<t_gpair, double> plane_eleT;
  map<t_gpair, double> plane_azim;  map<t_gpair, double> plane_dist;  map<t_gpair, double> plane_wght;
   
  map<t_gpair, double>::const_iterator  itA, itH, itW, itT, itX, itELH, itELW, itELT, itDST;
  
  bool temp_interp = false;
  // loop over epochs (user point)
  for( itEPO = rf_data.begin(); itEPO != rf_data.end(); ++itEPO)
  {
    t_gtime tt(itEPO->first);
    t_map_met::iterator itM = _map_met.find(tt);  // support only gradients in the main container!
     
    if( itM == _map_met.end() ){ 
      // cout << "*** Warning: analys_gradients - no data in map, skip epoch: " << tt.str_ymdhms() << endl;
      temp_interp = true;
      continue;
    }
    double rf_ztd = rf_data[tt][ZTD]; // [m]
    double rf_zhd = rf_data[tt][ZHD]; // [m]
    double rf_zwd = rf_data[tt][ZWD]; // [m]
     
    int count = 0;
    // loop over grid points
    for(t_map_dat::iterator itD = itM->second.begin(); itD != itM->second.end(); ++itD)
    {
      shared_ptr<t_gnwmsurf> gsurf = static_pointer_cast<t_gnwmsurf>(itD->second);
       
      t_gtriple grd( gsurf->get_surf(LAT),gsurf->get_surf(LON),gsurf->get_surf(HEL) ); // ALWAYS USE ELL

#ifdef DEBUG
      if( _proj ){ t_gtriple xyz = _proj->ll2proj( ell ); // USE THE PROJECTION!
                   cout << "MSG: USING THE LCC PROJECTION !" << fixed << setprecision(3)
	                << setw(14) << ell[0] << setw(14) << ell[1] << setw(14) << ell[2] << " --> "
	                << setw(14) << xyz[0] << setw(14) << xyz[1] << setw(14) << xyz[2] << endl;
      }
#endif

      // relative coordinates w.r.t. point of interes, given in [m]
      t_gpair ellgrd = grd.gpair();                                   // [deg] always LatLon !
      t_gpair ellrel(( ellgrd[0] - ell[0] )*D2R*A_WGS,                // [m]
                     ( ellgrd[1] - ell[1] )*D2R*A_WGS );              // [m]
       
      double t0     = gsurf->get_surfdata( TEMP );
      double dw     = gsurf->get_surf( dW );
      double hScZHD = Rd * t0 / G_WMO            / 1000.0; // [km] dry scale height
      double hScZWD = Rd * t0 / G_WMO / (dw + 1) / 1000.0; // [km] wet scale height
       
      double dist = sqrt(pow(ellrel[0],2) + pow(ellrel[1],2))/1000.0;  // [km]
      double azim = 90.0 - R2D*atan2(ellrel[0],ellrel[1]);             // [deg]
      if( azim < 0 ) azim += 360.0;
              
      double zwd_norm = gsurf->get_in_hel( ZWD, ell[2] );
      double zhd_norm = gsurf->get_in_hel( ZHD, ell[2] );
       
      double zhd_diff = zhd_norm            - rf_zhd; double eleH = R2D*(atan2(hScZHD,dist)-2*atan2(dist,2*R_SPHERE));
      double zwd_diff =            zwd_norm - rf_zwd; double eleW = R2D*(atan2(hScZWD,dist)-2*atan2(dist,2*R_SPHERE));
      double ztd_diff = zhd_norm + zwd_norm - rf_ztd; double eleT = R2D*(atan2(hScZHD,dist)-2*atan2(dist,2*R_SPHERE));
      
      if( eleH >= cut_off ){                              // [deg]                // selected points only
        plane_dist[ellrel] = dist;                        // [km]
        plane_eleH[ellrel] = eleH;                        // [deg]
        plane_eleW[ellrel] = eleW;                        // [deg]
        plane_eleT[ellrel] = eleT;                        // [deg]
        plane_azim[ellrel] = azim*D2R;                    // [deg]
        plane_datH[ellrel] = zhd_diff*1000*tan(eleH*D2R); // [mm]  ZHD gradient
        plane_datW[ellrel] = zwd_diff*1000*tan(eleW*D2R); // [mm]  ZWD gradient
        plane_datT[ellrel] = ztd_diff*1000*tan(eleT*D2R); // [mm]  ZTD gradient
        plane_wght[ellrel] = pow(sin(eleT*D2R),2.0);      // observation weight
              
#ifdef DEBUG
       double hel = gsurf->get_surf(HEL);
       double zwd = gsurf->get_surfdata(ZWD);
       double zhd = gsurf->get_surfdata(ZHD);
       cout << fixed << setprecision(0)
            << setw(3)   << count
	    << "  CRD:"  << setw(8) << ellrel[0] << " "     << setw(8) << ellrel[1] << setprecision(1)
	    << " "       << setw(7) << hel       << " ->"   << setw(7) << ell[2]    << setprecision(1)
	    << "  dst:"  << setw(6) << dist  
	    << "  sky:"  << setw(5) << eleH      << " / "   << setw(5) << azim      << setprecision(3)
	    << "  ZHD:"  << setw(6) << zhd       << " ->"   << setw(6) << zhd_norm
            << "  ZWD:"  << setw(6) << zwd       << " ->"   << setw(6) << zwd_norm
            << "  REF:"  << setw(6) << rf_ztd    << " "     << setw(8) << (zhd_norm + zwd_norm - rf_ztd)*1000
            << "  TIM:"  << tt.str_ymdhms()
	    << endl;
#endif
	count++;
      }
    }

#ifdef DEBUG
    count = 0;
    cout << "\n FIT for REQ PT: " << fixed << ell[0] << " " << ell[1] << " " << ell[2] << endl;
    for( itA   = plane_azim.begin(), itT   = plane_datT.begin(), itX  = plane_wght.begin(),
	 itH   = plane_datH.begin(), itW   = plane_datW.begin(), 
	 itELH = plane_eleH.begin(), itELW = plane_eleW.begin(),
	 itDST = plane_dist.begin(), itELT = plane_eleT.begin();
	 itA  != plane_azim.end();   ++itA,  ++itH, ++itW, ++itT, ++itX, ++itELH, ++itELW, ++itELT, ++itDST )
    {   count++;
        cout << fixed << setw(3)  << count << setprecision(0)
	              << "  crd:" << setw(7) << itA->first[0] 
	                          << setw(7) << itA->first[1]
	              << "  dst:" << setw(7) << itDST->second   // km
	              << "  elH:" << setw(6) << itELH->second   // deg
	              << "  elW:" << setw(6) << itELW->second   // deg
	              << "  elT:" << setw(6) << itELT->second   // deg
	              << "  azi:" << setw(6) << itA->second*R2D // deg
	                          << setprecision(3)
                      << "  gH:"  << setw(7) << itH->second     // mm
                      << "  gW:"  << setw(7) << itW->second     // mm
                      << "  gT:"  << setw(7) << itT->second     // mm
                      << "  wgt"  << setw(7) << itX->second     // weight
	              << endl;
    }
#endif

    if( plane_azim.size() < 4 ){ 
      if( _log ) _log->comment(0,"gallsurf","Error: gradient analysis has a few data. Skipped");
      else               cerr << "gallsurf - Error: gradient analysis has a few data. Skipped\n";
      return -1;
    }

    ColumnVector fitH(2); fitH << 0.0 << 0.0;  ColumnVector rmsH(2); rmsH << 0.0 << 0.0;
    ColumnVector fitW(2); fitW << 0.0 << 0.0;  ColumnVector rmsW(2); rmsW << 0.0 << 0.0;
    ColumnVector fitT(2); fitT << 0.0 << 0.0;  ColumnVector rmsT(2); rmsT << 0.0 << 0.0;
     
    ColumnVector   A(plane_azim.size());  itA = plane_azim.begin();
    ColumnVector   H(plane_datH.size());  itH = plane_datH.begin();
    ColumnVector   W(plane_datW.size());  itW = plane_datW.begin();
    ColumnVector   T(plane_datT.size());  itT = plane_datT.begin();
    DiagonalMatrix P(plane_wght.size());  itX = plane_wght.begin(); // Weighting

    map<string, double> cf;
    string msgH, msgW, msgT;

    int itrH, itrW, itrT; 
    t_gfunc_linCS funcCS(cf); 
    t_gfitlmm gfit(&funcCS);

    for(unsigned int i=0; i < plane_azim.size(); ++i, ++itA, ++itH, ++itW, ++itT, ++itX  ){
      A[i] = itA->second;  T[i] = itT->second;
      H[i] = itH->second;  W[i] = itW->second;
      P[i] = itX->second;
    }

    if( gfit.fit( A, H, P, fitH, rmsH, itrH, msgH ) < 0 ){ cerr << "Error: fitting ZHD gradients\n"; return -1; }
    if( gfit.fit( A, W, P, fitW, rmsW, itrW, msgW ) < 0 ){ cerr << "Error: fitting ZWD gradients\n"; return -1; }
    if( gfit.fit( A, T, P, fitT, rmsT, itrT, msgT ) < 0 ){ cerr << "Error: fitting ZTD gradients\n"; return -1; }

     
#ifdef DEBUG
	  cout << fixed    << setw(0) << itEPO->first.str_ymdhms()  <<  setprecision(2) 
//                         << setw(4) << itEPO->first.doy()
//	                   << setw(6) << itEPO->first.sod()
	       << "  ZHD:" << setw(5) << rf_zhd  << setw(2) << itrH
               << " "      << setw(6) << fitH(1) << setw(5) << rmsH(1)
	       << " "      << setw(6) << fitH(2) << setw(5) << rmsH(2)
	       << "  ZWD:" << setw(5) << rf_zwd  << setw(2) << itrW
	       << " "      << setw(6) << fitW(1) << setw(5) << rmsW(1)
	       << " "      << setw(6) << fitW(2) << setw(5) << rmsW(2)
	       << "  ZTD:" << setw(5) << rf_ztd  << setw(2) << itrT
	       << " "      << setw(6) << fitT(1) << setw(5) << rmsT(1)
	       << " "      << setw(6) << fitT(2) << setw(5) << rmsT(2)
  	       << " "      << setw(6) << fitH(1) + fitW(1)
	       << " "      << setw(6) << fitH(2) + fitW(2)
               << endl;
#endif

    inp_ngrd[tt] = rf_data[tt][NGRD] = (fitH(1) + fitW(1))/1000.0; // [mm] -> [m]
    inp_egrd[tt] = rf_data[tt][EGRD] = (fitH(2) + fitW(2))/1000.0; // [mm] -> [m]
       
    plane_dist.clear();  plane_azim.clear();  plane_wght.clear();
    plane_eleH.clear();  plane_eleW.clear();  plane_eleT.clear();
    plane_datH.clear();  plane_datW.clear();  plane_datT.clear();
  }

  // optional loop over epochs (user point) - interpolate if data were not available in _met_dat
  if( temp_interp && rf_data.size() > 0 ){
    t_gtime beg = inp_ngrd.begin()->first; 
    t_gtime end = inp_ngrd.rbegin()->first;
    if( beg != end ){
      _interp_temporal( beg, end, sampling, inp_ngrd, out_ngrd );
      _interp_temporal( beg, end, sampling, inp_egrd, out_egrd );
     
      for( itEPO = rf_data.begin(); itEPO != rf_data.end(); ++itEPO){
        t_gtime tt(itEPO->first);
	rf_data[tt][NGRD] = out_ngrd[tt];
	rf_data[tt][EGRD] = out_egrd[tt];
      }
    }
  }
  return 0;
}


// new element allocated in heap
// ----------
int t_gallsurf::_new_element(shared_ptr<t_gnwmsurf> met)
{
  return 0;
}


// set region of interest
// ----------
int t_gallsurf::_set_region(const double minLat,
	 	 	    const double maxLat,
	 		    const double minLon,
	 		    const double maxLon
			   )
{
  gtrace("t_gallsurf::_set_region");

  _min_lat = minLat;
  _max_lat = maxLat;
  _min_lon = minLon;
  _max_lon = maxLon;
   
  if( _min_lat != MIN_LAT && _max_lat != MAX_LAT && 
      _min_lon != MIN_LON && _max_lon != MAX_LON ) _reg_auto = false;

  return 0;
}


// Set param
// -----------------------------
void t_gallsurf::_set_param( const set<string>& param )
{
  gtrace("t_gallsurf::_set_param");

  // the parameters requested in the interpolation
  if( param.find("TEMP")  != param.end() ) _vec_par.push_back(TEMP); // MUST BE FIRST
  if( param.find("T")     != param.end() ) _vec_par.push_back(TEMP);
  if( param.find("PRES")  != param.end() ) _vec_par.push_back(PRES); // MUST BE SECOND
  if( param.find("P")     != param.end() ) _vec_par.push_back(PRES);
  if( param.find("EPRE")  != param.end() ) _vec_par.push_back(EPRE);
  if( param.find("E")     != param.end() ) _vec_par.push_back(EPRE);
   
  if( param.find("IWV")   != param.end() ) _vec_par.push_back(IWV);
  if( param.find("ZHD")   != param.end() ) _vec_par.push_back(ZHD);
  if( param.find("ZWD")   != param.end() ) _vec_par.push_back(ZWD);
  if( param.find("ZTD")   != param.end() )
  {	
   if( param.find("ZHD")  == param.end() ){_vec_par.push_back(ZHD);}
   if( param.find("ZWD")  == param.end() ){_vec_par.push_back(ZWD);}
                                           _vec_par.push_back(ZTD);  // MUST BE AFTER ZHD+ZWD
  }
   
  if( param.find("SHUM")  != param.end() ) _vec_par.push_back(SHUM);
  if( param.find("Q")     != param.end() ) _vec_par.push_back(SHUM);
  if( param.find("TM")    != param.end() ) _vec_par.push_back(TM);
  if( param.find("GP")    != param.end() ) _vec_par.push_back(GEOP);
  if( param.find("G")     != param.end() ) _vec_par.push_back(GRAV);

  // not yet supported (due to resulting map from the interpolation)
  if( param.find("DM")    != param.end() ) _vec_srf.push_back(dM);
  if( param.find("DT")    != param.end() ) _vec_srf.push_back(dT);
  if( param.find("DE")    != param.end() ) _vec_srf.push_back(dE);
  if( param.find("DW")    != param.end() ) _vec_srf.push_back(dW);
  if( param.find("DI")    != param.end() ) _vec_srf.push_back(dI);
//if( param.find("DR")    != param.end() ) _vec_srf.push_back(dR);

  if( param.find("SCH")   != param.end() ) _vec_srf.push_back(scH);
  if( param.find("SCW")   != param.end() ) _vec_srf.push_back(scW);
  if( param.find("SCE")   != param.end() ) _vec_srf.push_back(scE);
   
  if( param.find("EQH")   != param.end() ) _vec_srf.push_back(eqH);
  if( param.find("EQW")   != param.end() ) _vec_srf.push_back(eqW);
  if( param.find("EQE")   != param.end() ) _vec_srf.push_back(eqE);

//if( param.find("GEOM")  != param.end() ) _vec_par.push_back(GEOM); // not SURF, but DATA
  if( param.find("HGT")   != param.end() ||
      param.find("HEL")   != param.end() ||
      param.find("H")     != param.end() ) _vec_srf.push_back(HEL);

  if( param.find("RATIO") != param.end() ||
      param.find("RATX")  != param.end() ) _vec_par.push_back(RATX);

  return;
}

// Set output
// -----------------------------
void t_gallsurf::_set_output()
{
  gtrace("t_gallsurf::_set_output");

  if( ! _set ) return;

  string tmp;
  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("out");
  if( ! tmp.empty() ){
//    substitute(tmp,"$(rec)",_site, false);
    _out = new t_giof;
    _out->tsys();
    _out->mask(tmp);
    _out->append( dynamic_cast<t_gsetout*>(_set)->append() );
  }
   
  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("sum");
  if( ! tmp.empty() ){
//    substitute(tmp,"$(rec)",_site, false);
    _sum = new t_giof;
    _sum->tsys();
    _sum->mask(tmp);
    _sum->append( dynamic_cast<t_gsetout*>(_set)->append() );
  }

  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("grd");
  if( ! tmp.empty() ){
    _grd = new t_giof;
    _grd->tsys();
    _grd->mask(tmp);
    _grd->append( dynamic_cast<t_gsetout*>(_set)->append() );
  }


  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("srf");
  if( ! tmp.empty() ){
    _srf = new t_giof;
    _srf->tsys();
    _srf->mask(tmp);
    _srf->append( dynamic_cast<t_gsetout*>(_set)->append() );
  }

  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("fit");
  if( ! tmp.empty() ){
    _fit = new t_giof;
    _fit->tsys();
    _fit->mask(tmp);
    _fit->append( dynamic_cast<t_gsetout*>(_set)->append() );
  }

  return;
}

} // namespace
