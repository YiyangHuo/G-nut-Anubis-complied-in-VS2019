
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
#include <iostream>
#include <iomanip>
#include <sstream>

#include "gall/gallobs.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallobs::t_gallobs() 
  : t_gdata(),
    _set(0),
    _nepoch(0),
    _overwrite(false)
{
  gtrace("t_gallobs::constructor");
  id_type(  t_gdata::ALLOBS );
  id_group( t_gdata::GRP_OBSERV );
}


// destructor
// ----------
t_gallobs::~t_gallobs()
{
  gtrace("t_gallobs::destructor");
  _mapobj.clear();
  _filter.clear();
}


// settings
// ----------
void t_gallobs::gset(t_gsetbase* gset)
{  
  _set = gset;
  _sys = dynamic_cast<t_gsetgen*>(_set)->sys();
  _smp = dynamic_cast<t_gsetgen*>(_set)->sampling();
  _scl = dynamic_cast<t_gsetgen*>(_set)->sampling_scalefc(); // scaling 10^decimal-digits

//  _diff_sec = NOMINAL_DIFF_SEC;
//  if( _sampling < 1.0 ) _diff_sec *= _sampling;

  return;
}


// return set of available GNSS systems
// --------------------------
set<GSYS> t_gallobs::sys(string site)
{
  gtrace("t_gallobs::sys");
   t_map_oref::const_iterator itEpo;
   t_map_osat::iterator itSat;
   set<GSYS> temp;
   t_gtime t;
   
   for( itEpo = _mapobj[site].begin(); itEpo != _mapobj[site].end(); ++itEpo) // loop over epochs
   {
     t = itEpo->first;
     for( itSat = _mapobj[site][t].begin(); itSat != _mapobj[site][t].end(); ++itSat) // loop over satellites
     {
	temp.insert(itSat->second->gsys());
     }
   }   
   return temp;
}


// return used systems for epoch
// -------------------------------
set<GSYS> t_gallobs::sys(const string& site, const t_gtime& t)
{
  gtrace("t_gallobs::sys");
   
  _gmutex.lock();
  set<GSYS> gnss = _gsys(site, t);
  _gmutex.unlock(); return gnss;
}
   
   
// protected
// ---------
set<GSYS> t_gallobs::_gsys( const string& site, const t_gtime& t)
{
  set<GSYS> gnss;

  vector<t_gsatdata> satdata = _gobs(site, t);
 
  vector<t_gsatdata>::iterator it;
  for (it = satdata.begin(); it != satdata.end(); it++){
//    cout << it->epoch().str_hms() << " " << it->sat() << endl;
    GSYS gsys = it->gsys();
    gnss.insert(gsys);
  }

  return gnss;
}
   

// return list of available stations
// ----------
set<string> t_gallobs::stations()
{
  gtrace("t_gallobs::stations");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  set<string> all_sites;
  t_map_oobj::const_iterator itSITE = _mapobj.begin();
   
  while( itSITE != _mapobj.end() ){
    all_sites.insert(itSITE->first);
    itSITE++;
  }
  _gmutex.unlock(); return all_sites;
}


// return list of available satellites
// ----------
set<string> t_gallobs::sats(const string& site, 
                            const t_gtime& t, GSYS gnss)
{
  gtrace("t_gallobs::sats");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  set<string> all_sats = _sats(site,t,gnss);
  _gmutex.unlock(); return all_sats;
}


// return list of available satellites
// ----------
set<string> t_gallobs::_sats(const string& site,
                             const t_gtime& t, GSYS gnss)
{
  set<string> all_sats;

  if( _mapobj.find(site)    == _mapobj.end() ||
      _mapobj[site].find(t) == _mapobj[site].end() ){ return all_sats; }

  t_map_osat::const_iterator itSAT = _mapobj[site][t].begin();  
  while( itSAT != _mapobj[site][t].end() ){
    GSYS sys = itSAT->second->gsys();
    if(gnss == GNS) all_sats.insert(itSAT->first);
    else if(gnss == GPS && sys == GPS) all_sats.insert(itSAT->first);
    else if(gnss == GLO && sys == GLO) all_sats.insert(itSAT->first);
    else if(gnss == GAL && sys == GAL) all_sats.insert(itSAT->first);
    else if(gnss == BDS && sys == BDS) all_sats.insert(itSAT->first);
    else if(gnss == QZS && sys == QZS) all_sats.insert(itSAT->first);
    else if(gnss == SBS && sys == SBS) all_sats.insert(itSAT->first);     
    itSAT++;
  }
  return all_sats;
}


// return list of available satellites
// ----------
vector<t_gsatdata> t_gallobs::obs(const string& site,
                                  const t_gtime& t)
{
  gtrace("t_gallobs::obs");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  vector<t_gsatdata> all_obs = _gobs(site,t);

  _gmutex.unlock(); return all_obs;
}


// protected
// ---------
vector<t_gsatdata> t_gallobs::_gobs(const string& site,
                                    const t_gtime& t)
{
  vector<t_gsatdata> all_obs;
  t_gtime tt(t_gtime::GPS);

  if (_find_epo(site, t, tt) < 0) {
     return all_obs;
  }

  t_map_osat::iterator itSAT = _mapobj[site][tt].begin();
  while( itSAT != _mapobj[site][tt].end() ){    
//    all_obs.push_back( itSAT->second ); 
    all_obs.push_back( *itSAT->second );
    itSAT++;
  }

  return all_obs;
}


// TEMPORARY !!!
// return list of available satellites (POINTERS!)
// ----------
vector<t_spt_gobs> t_gallobs::obs_pt(const string& site, const t_gtime& t)
{
  gtrace("t_gallobs::obs_pt");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<t_spt_gobs> all_obs;
  t_gtime tt(t_gtime::GPS);
   
  if( _find_epo(site, t, tt) < 0 ){
    _gmutex.unlock(); return all_obs;
  }
               
  t_map_osat::iterator itSAT = _mapobj[site][tt].begin();
  while( itSAT != _mapobj[site][tt].end() ){
     
    string sat = itSAT->first;
     
    // TESTING NEW METHOD
    all_obs.push_back( dynamic_pointer_cast<t_gobsgnss>(itSAT->second) );
 
    itSAT++;
  }

  _gmutex.unlock(); return all_obs;
}


// return single difference between satellites
// ---------
/*t_gsdbs t_gallobs::sdbs(const string& site,
                        const t_gtime& t)
{
   t_gsdbs diff;
   
   vector<t_gsatdata> gobs = this->obs(site, t);   
   diff.pivot( *gobs.begin() );
      
   for (unsigned int i = 1; i < gobs.size(); i++){
      gobs[i] = gobs[i] - gobs[0];      
   }        
   
   gobs.erase(gobs.begin()); // erase reference observation

   diff.data(gobs);
   
   return diff;
}
*/


// return list of available satellites
// ----------
vector<t_gtime> t_gallobs::epochs(const string& site)
{
  gtrace("t_gallobs::epochs");
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<t_gtime> all_epochs;

  if( _mapobj.find(site) == _mapobj.end() ){
    // cout << "site not found: " << site << t.str_ymdhms() << endl; cout.flush();
     _gmutex.unlock(); return all_epochs;
  }
     
  t_map_oref::iterator it = _mapobj[site].begin();

  while( it != _mapobj[site].end() ){
    all_epochs.push_back( it->first );
    it++;
  }
  
  _gmutex.unlock(); return all_epochs;
}


// return first position for satellite
// ----------
t_gtime t_gallobs::beg_obs(const string& site, double smpl)
{
  gtrace("t_gallobs::beg_obs");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = LAST_TIME;

  if( _mapobj.find(site)    != _mapobj.end() &&
      _mapobj[site].begin() != _mapobj[site].end() ) tmp = _mapobj[site].begin()->first;

  // get first synchronized obs
  if( smpl > 0.0 ){
    auto itEpoB = _mapobj[site].begin();
    auto itEpoE = _mapobj[site].end();
    int sod = static_cast<int>(dround(itEpoB->first.sod()+itEpoB->first.dsec()));

    while( sod % static_cast<int>(smpl) != 0 && ++itEpoB != itEpoE ){
      sod = static_cast<int>(dround(itEpoB->first.sod()+itEpoB->first.dsec()));
      tmp = itEpoB->first;
      tmp.reset_sod();
      tmp.add_secs(sod); 
//      cerr << " sync sod = " << fixed << setprecision(4) << sod << "  " << itEpoB->first.str_ymdhms() << "  " << tmp.str_ymdhms() << endl;
    }
  }
   
  _gmutex.unlock(); return tmp;
}


// return last position for satellite
// ----------
t_gtime t_gallobs::end_obs(const string& site)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = FIRST_TIME;
  if( _mapobj.find(site)    != _mapobj.end() &&
      _mapobj[site].begin() != _mapobj[site].end() )
    tmp = _mapobj[site].rbegin()->first;

  _gmutex.unlock(); return tmp;
}


// add observations ( both P and L in meters !!!)
// ----------
int t_gallobs::addobs(t_spt_gobs obs)
{
  gtrace("t_gallobs::addobs");
   
  t_gtime t(obs->epoch()),   tt = t;
  string site = obs->site();
  string sat  = obs->sat();
//  cout << "ADDING site: " << site << endl;

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  int  epo_found = _find_epo(site, t, tt);
  auto itSAT     = _mapobj[site][tt].find(sat);

//  auto itEPO = _mapobj[site].find(tt);
//  auto itSAT = _mapobj[site][tt].find(sat);

  // add new observations (or overwrite)
  // ===================================
  if( _overwrite || epo_found < 0                     // epoch exists (smart search)
                 || itSAT == _mapobj[site][tt].end()  // satellite exists
//                 || itEPO == _mapobj[site].end()
//                 || itSAT == _mapobj[site][tt].end()
  ){
#ifdef DEBUG
    cout << " ADDING " << site << " " <<  _mapobj[site].size() << " " << _nepoch // << " " <<  nepochs(site) // MUTEX !!!
         << " sat "    <<  sat << " " << tt.str_ymdhms() << " " << "\n";cout.flush();
#endif

    // delete old if exists
    // ====================
    if( _overwrite && epo_found > 0                     // epoch exists (smart search)
                   && itSAT != _mapobj[site][tt].end()  // satellite exists    
//    if( _overwrite && itEPO != _mapobj[site].end()
//                   && itSAT != _mapobj[site][tt].end()
    ){
        _mapobj[site][tt].erase(sat);
        if( _log && _log->verb() >= 2 ) _log->comment(2, "gallobs", site + tt.str_ymdhms(" obs replaced ") + sat);
    }

    // too many epochs (remove old)
    // ============================
    if( _nepoch > 0 && _mapobj[site].size() > _nepoch + 10 ){  // +10 .. reduce number of removals to every tenths

      auto itEPO = _mapobj[site].begin();
      while( itEPO != _mapobj[site].end() ){

        if( _mapobj[site].size() <= _nepoch ){ break; }
	
        t_gtime t = itEPO->first;
        if( _log && _log->verb() >= 2 )
        {
           string satells(" ");
           for(auto itSAT = _mapobj[site][t].begin(); itSAT != _mapobj[site][t].end(); ++itSAT ){
             satells = satells + " " + itSAT->first;
           }
       
           _log->comment(2, "gallobs", site + t.str_ymdhms(" obs removed ") + satells);
        }

        _mapobj[site][t].erase(_mapobj[site][t].begin(),_mapobj[site][t].end());
        _mapobj[site].erase(itEPO++);
      }
    }

    // distinguish GNSS specific message
    // =================================
    if( obs->id_type() == t_gdata::OBSGNSS ){

      _mapobj[site][tt][sat] = obs;

    }else{
       if( _log ){  _log->comment(0,"gallobs","warning: t_gobsgnss record not identified!"); }
       _gmutex.unlock(); return 1;
    }
  }else{
    if( _log && _log->verb() >= 2 ) _log->comment(2,"gallobs", site + tt.str_ymdhms(" skipped ") + " " + sat);
    _gmutex.unlock(); return 0;
  }

  // comments
  if( _log && _log->verb() >= 3 ){

    vector<GOBS> v_obs = obs->obs();
    vector<GOBS>::iterator itOBS = v_obs.begin();
    ostringstream lg; lg << fixed << setprecision(3);
    for( ; itOBS != v_obs.end(); ++itOBS ) 
      lg << " " << gobs2str( *itOBS ) << ":" << obs->getobs( *itOBS );

    _log->comment(3,"gallobs",site + tt.str_ymdhms(" add obs: ") + " " + sat + lg.str() );
  }
  _gmutex.unlock(); return 0;
}


// clean function
// ----------
void t_gallobs::clean_outer( const  string& obj,
                             const t_gtime& beg, 
                             const t_gtime& end )
{
  gtrace("t_gallobs::clean_outer");

  if( end < beg ) return;

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  // loop over all object [default]
  t_map_oobj::const_iterator itOBJ  = _mapobj.begin();
  t_map_oobj::const_iterator itOBJ2 = _mapobj.end();
    
  if( _log ){
    _log->comment(2,"gallobs","obs clean request: " + obj + beg.str_ymdhms(": ") + end.str_ymdhms(" - "));
  }

  // loop over a single object (if obj defined !)
  if( ! obj.empty() ) itOBJ = itOBJ2 = _mapobj.find(obj);
  if( itOBJ2 != _mapobj.end() ) itOBJ2++;

  for( ; itOBJ != itOBJ2; ++itOBJ ){

    string site = itOBJ->first;

    // find and CLEAN all data (epochs) out of the specified period !
    t_map_oref::iterator it;
    t_map_oref::iterator itFirst = _mapobj[site].begin();
    t_map_oref::iterator itLast  = _mapobj[site].end();
    t_map_oref::iterator itBeg   = _mapobj[site].lower_bound(beg);  // greater|equal
    t_map_oref::iterator itEnd   = _mapobj[site].upper_bound(end);  // greater only!

    t_map_osat::iterator itOBS;

    // delete before BEGIN request
    if (itBeg != itLast&& itBeg != itFirst){
      if( _log && _log->verb() >= 3  ){
        string lg("obs removed before: " + itBeg->first.str_ymdhms() + " " + site);
        _log->comment(3,"gallobs",lg);
      }

      for( it = itFirst; it != itBeg; ++it ){	 
        t_gtime tt(it->first);
	 
        _mapobj[site][tt].erase(_mapobj[site][tt].begin(),_mapobj[site][tt].end());
      }
       
      _mapobj[site].erase(itFirst,itBeg); // erase all but itBeg
    }

    // delete after END request
    if( itEnd != itLast ){ // && ++itEnd != itLast ){
      if( _log && _log->verb() >= 3 ){
        string lg("obs removed after : " + itEnd->first.str_ymdhms() + " " + site);
        _log->comment(3,"gallobs",lg);
      }

      for( it = itEnd; it != itLast; ++it ){
	    
        t_gtime tt(it->first);
	 
        _mapobj[site][tt].erase(_mapobj[site][tt].begin(),_mapobj[site][tt].end());	 
      }                   
      _mapobj[site].erase(itEnd,itLast);
    }
//    itOBJ++; // WHILE LOOP ONLY
  }
  _gmutex.unlock(); return;
}


// return list of available satellites
// ----------
t_gallobs::t_map_osat t_gallobs::find( string site, const t_gtime& t)
{
  gtrace("t_gallobs::find");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tt(t_gtime::GPS);
  t_map_osat  tmp;
  if(_find_epo(site, t, tt) < 0){
     _gmutex.unlock(); return tmp;
  }

//  shared_ptr<t_map_osat> tmp = &_mapobj[site][tt];
  _gmutex.unlock(); return _mapobj[site][tt]; // tmp;
}


// return list of available satellites
// ----------   
double t_gallobs::find( string site, const t_gtime& t, const string& prn, const GOBS& gobs)
{
   gtrace("t_gallobs::find");

#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif
   _gmutex.lock();
   
   double obs = 0.0;   
   
   t_gtime tt(t_gtime::GPS);
   t_map_osat  tmp;
   if(_find_epo(site, t, tt) < 0){
      _gmutex.unlock(); return obs;
   }   
   
   t_map_osat mdata = _mapobj[site][tt];
   t_map_osat::iterator it = mdata.find(prn);

   if( it != mdata.end() ){
      t_spt_gobs satdata = it->second;
      vector<GOBS> vobs = satdata->obs(); 
      for(vector<GOBS>::iterator it2 = vobs.begin(); it2 != vobs.end(); it2++){
        if(*it2 == gobs) obs = satdata->getobs(*it2);
        else continue;
      }
   }
   
   _gmutex.unlock(); return obs;
}


// get number of occurance of individual signals
// ----------
t_gallobs::t_map_frq t_gallobs::frqobs(string site)
{
   gtrace("t_galloqc::frqobs");
   
   _gmutex.lock();
   
   t_map_frq mfrq;
   
   if( _mapobj.find(site) == _mapobj.end() ){ _gmutex.unlock(); return mfrq; }   
   
   t_map_oref::const_iterator itEpo = _mapobj[site].begin();
   t_map_osat::iterator itSat;

   t_gtime t;
   
   while (itEpo != _mapobj[site].end()) {
      t = itEpo->first;
      for(itSat = _mapobj[site][t].begin(); itSat != _mapobj[site][t].end(); itSat++) { // loop over satellites
	 string prn     = itSat->first;
	 t_spt_gobs obs = itSat->second;
	 
	 vector<GOBS> vgobs = obs->obs();
	 for(vector<GOBS>::iterator it = vgobs.begin(); it != vgobs.end(); it++){
	    GOBS gobs = *it;
	    int bnd = gobs2band(gobs);
	    GOBSBAND bend = int2gobsband(bnd);
	    mfrq[prn][bend][gobs]++;
	 }
	 
      }
      itEpo++;
   }
   
   _gmutex.unlock();
   return mfrq;
}
   

// add site-specific filtered data/epochs
// ----------
void t_gallobs::xdata(string site, string file, t_xfilter xflt)
{
  gtrace("t_gallobs::xdata");

  _gmutex.lock();
  
  _filter[site][file].xdat = xflt.xdat;
  _filter[site][file].beg  = xflt.beg;
  _filter[site][file].end  = xflt.end;

  _gmutex.unlock();
}


// get site-specific filtered data/epochs
// ---------
t_gallobs::t_xfilter t_gallobs::xdata(string site, string file)
{
  gtrace("t_gallobs::xdata");

  _gmutex.lock();

  if( _filter.find(site) != _filter.end() ){
    if( _filter[site].find(file) != _filter[site].end() ){
      _gmutex.unlock(); return _filter[site][file];
    }
  }

  t_xfilter tmp;
  tmp.xdat[XDATA_BEG] = 0;
  tmp.xdat[XDATA_END] = 0;
  tmp.xdat[XDATA_SMP] = 0;
  tmp.xdat[XDATA_SYS] = 0;

  tmp.beg = LAST_TIME;
  tmp.end = FIRST_TIME;

  _gmutex.unlock(); return tmp;
}


// number of epochs for station
// ----------
unsigned int t_gallobs::nepochs(const string& site)
{
  gtrace("t_gallobs::nepochs");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( _mapobj.find(site) == _mapobj.end() ){ _gmutex.unlock(); return 0; }

  unsigned int tmp = _mapobj[site].size();
/*
  if( site == "STAN" ){
    t_map_oref::iterator itEPO = _mapobj[site].begin();
    while( itEPO != _mapobj[site].end() ){
       
      t_gtime t = itEPO->first;
      cout << site << " --- time: " << t.str_ymdhms() << " " << t.dsec() << "  size: " << tmp << endl;
      itEPO++;
    }
  }
*/

  _gmutex.unlock(); return tmp;
}


// number of epochs for station, interval and sampling, list of OBS types
// ----------
unsigned int t_gallobs::nepochs(const string& site, const t_gtime& beg, const t_gtime& end, double sampl, map<GSYS, pair<int,int> >& n)
{
  gtrace("t_gallobs::nepochs");

  if( _mapobj.find(site) == _mapobj.end() ) return 0;   

  set<GSYS>::iterator itS;
  map<GSYS, pair<int, int> >::iterator it;

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  n.clear();

  t_gtime tt = beg;
  while( tt < end || tt == end ){
    set<GSYS> gnss = _gsys(site, tt);
    for(itS = gnss.begin(); itS != gnss.end(); itS++) n[*itS].second += 1;
//  tt = tt + sampl;
    tt.add_dsec(sampl);
  }

  for( it = n.begin(); it != n.end(); it++ )
    it->second.first = (int)floor( (end - beg) / sampl + DIFF_SEC(sampl)) + 1;
   
  _gmutex.unlock(); return 1;
}


// return true if map contains the site and false if not
//-----------------------------------------------------
bool t_gallobs::isSite(string site)
{
  gtrace("t_gallobs::isSite");

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
 
  bool tmp = false;
  if ( _mapobj.find(site) != _mapobj.end() ) tmp = true;

  _gmutex.unlock(); return tmp;
}


// get all t_gsatdata for prn from epoch interval
// ----------------------------------------------
vector<t_spt_gobs> t_gallobs::obs_prn_pt(const string& site, const string& prn, 
			  	 	 const t_gtime& beg, const t_gtime& end)
{
  gtrace("t_gallobs::obs_prn_pt");

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
 
  vector<t_spt_gobs> all_obs;
  t_gtime tt(t_gtime::GPS);

  if( _mapobj.find(site) == _mapobj.end() ){
#ifdef DEBUG     
      cout << "site not found: " << site << t.str_ymdhms() << endl; cout.flush();
#endif
    _gmutex.unlock(); return all_obs;
  }       

  t_map_oref::iterator it1 = _mapobj[site].lower_bound(beg);  // greater || equal
  t_map_oref::iterator it2 = _mapobj[site].lower_bound(end);  // greater || equal   

  for ( ;it1 != it2; ++it1)
  {
     tt = it1->first;
     if( _mapobj[site].find(tt) != _mapobj[site].end() )
     {
       auto itSAT = _mapobj[site][tt].find(prn);
       if( itSAT != _mapobj[site][tt].end() ) all_obs.push_back( itSAT->second );
     }
  }
  _gmutex.unlock(); return all_obs;
}


// find epoch from the map w.r.t. DIFF_SEC
// ---------------------------------------
int t_gallobs::_find_epo(const string& site, const t_gtime& epo, t_gtime& tt)
{
  gtrace("t_gallobs::_find_epo");

  if( _mapobj.find(site) == _mapobj.end() ){
     if( _log && _log->verb() >= 2 ){
        _log->comment(3,"gallobs","site not found: " + site + " " + epo.str_ymdhms());
      }
      return -1;
  }
   
  t_map_oref::iterator it1 = _mapobj[site].lower_bound(epo);  // greater || equal
  t_map_oref::iterator it0 = it1;                             // previous value

  if( it0 != _mapobj[site].begin() ) it0--;                               // set previous value  
  if( it1 == _mapobj[site].end() ) it1 = it0;

  if( it1 == _mapobj[site].end() && it0 == _mapobj[site].end() ){
    cout << site << " " << epo.str_ymdhms() << " " << epo.dsec() << " no observations !\n";
    return -1;
  }
     
#ifdef DEBUG
    cout << site << " " << epo.str_ymdhms() << " " << epo.dsec() << fixed << setprecision(6)
         << " dif_sec: " << DIFF_SEC(_smp)
         << " t1_diff: " << setw(14) << it1->first.diff(epo)
         << " t0_diff: " << setw(14) << it0->first.diff(epo) << endl;
#endif

  if(      fabs(it1->first - epo) <= DIFF_SEC(_smp) ) tt = it1->first;  // it = it1;                 // set closest value
  else if( fabs(it0->first - epo) <= DIFF_SEC(_smp) ) tt = it0->first;  // it = it0;                 // set closest value
  else{
#ifdef DEBUG
     if( _log ){
       ostringstream lg;
       lg << site << tt.str_ymd()  << tt.dsec()
          << " it0 = " << it0->first.str_ymdhms() << " diff: " << it1->first.diff(epo)
          << " it1 = " << it1->first.str_ymdhms() << " diff: " << it0->first.diff(epo)
          << " not close observations found!";
       _log->comment(0,"gallobs",lg.str());
     }     
#endif
     return -1;                                     // not found !
  }
   
   return 1;
}

   
} // namespace
