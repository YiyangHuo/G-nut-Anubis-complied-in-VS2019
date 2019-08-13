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
#include <sstream>
#include <iomanip>
#include <cmath>

#include "gall/gallnav.h" 
#include "gutils/gsys.h"
#include "gutils/gtimesync.h"
#include "gutils/gtypeconv.h"
#include "gutils/gstat.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallnav::t_gallnav() 
  : t_gdata(),
    _com(false),
    _offset(0),
    _nepoch(t_gtime::GPS),
    _multimap(false),
    _overwrite(false),
    _chk_health(true),
    _chk_navig(true),
    _chk_tot(false)
{
  gtrace("t_gallnav::constructor");
  id_type(  t_gdata::ALLNAV );
  id_group( t_gdata::GRP_EPHEM);
}


// destructor
// ----------
t_gallnav::~t_gallnav()
{
  gtrace("t_gallnav::destructor");

  _mapsat.clear();
}

   
// return gnav element
// ----------
shared_ptr<t_geph> t_gallnav::find( string sat, const t_gtime& t, bool chk_mask )
{	
  gtrace("t_gallnav::find");

  _gmutex.lock();

  shared_ptr<t_geph> tmp = t_gallnav::_find( sat, t, _chk_health && chk_mask );

  _gmutex.unlock(); return tmp;
};

// return gnav elements
// ----------
vector<shared_ptr<t_geph>> t_gallnav::find_mult( string sat, const t_gtime& t )
{
  gtrace("t_gallnav::find_mult");

  _gmutex.lock();

  vector<shared_ptr<t_geph>> vec = this->_find_mult(sat, t);

  _gmutex.unlock(); 
  return vec; 
}


// return position
// ----------
int t_gallnav::pos( string sat, const t_gtime& t, double  xyz[],
		                               double  var[], double  vel[], bool chk_mask ) // [m]
{
  gtrace("t_gallnav::pos");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  shared_ptr<t_geph> tmp = t_gallnav::_find( sat, t, _chk_health && chk_mask );

  if( tmp == _null ){
    for(int i = 0; i<3; i++){
	                     xyz[i] = 0.0;
	           if( var ) var[i] = 0.0;
	           if( vel ) vel[i] = 0.0;
     }
     _gmutex.unlock(); return -1;
  }         

  int irc = tmp->pos( t, xyz, var, vel, _chk_health && chk_mask );
      
  _gmutex.unlock(); return irc;
   
//  return find( sat, t )->pos( t, xyz, var, vel, _chk_health && chk_mask );
}


// return satellite health
// ----------
bool t_gallnav::health( string sat, const t_gtime& t )
{
  gtrace("t_gallnav::health");
  
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  
  shared_ptr<t_geph> tmp = t_gallnav::_find( sat, t, false );

  if( tmp == _null ){ _gmutex.unlock(); return false; }
  
  bool status = tmp->healthy();

 _gmutex.unlock(); return status;
}


// return position aka navigation
// ----------
int t_gallnav::nav( string sat, const t_gtime& t, double  xyz[],
		                               double  var[], double  vel[], bool chk_mask ) // [m]
{
  gtrace("t_gallnav::nav");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  shared_ptr<t_geph> tmp = t_gallnav::_find( sat, t, _chk_health && chk_mask );

  if( tmp == _null ){	
    for(int i = 0; i<3; i++){
	                     xyz[i] = 0.0;
	           if( var ) var[i] = 0.0;
	           if( vel ) vel[i] = 0.0;
     }
     _gmutex.unlock(); return -1;
  }     
  int irc = tmp->nav( t, xyz, var, vel, _chk_health && chk_mask );

  _gmutex.unlock(); return irc;
   
//  return find( sat, t )->pos( t, xyz, var, vel, _chk_health && chk_mask );
}


// return clock corrections
// ----------
int t_gallnav::clk( string sat, const t_gtime& t, double*  clk,
		                                double*  var, double* dclk, bool chk_mask ) // [s]
{
  gtrace("t_gallnav::clk");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  shared_ptr<t_geph> tmp = t_gallnav::_find( sat, t, _chk_health && chk_mask );

  if( tmp == _null ){
              *clk  = 0.0;
    if( var ) *var  = 0.0;
    if(dclk ) *dclk = 0.0;
    _gmutex.unlock(); return -1;
  }

  int irc = tmp->clk( t, clk, var, dclk, _chk_health && chk_mask );
  _gmutex.unlock(); return irc;
   
//  return this->find( sat, t )->clk( t, clk, var, dclk, _chk_health && chk_mask );
}


// print function
// -------------------
void t_gallnav::print(string sat, const t_gtime& t)
{
  gtrace("t_gallnav::print");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  shared_ptr<t_geph> tmp = t_gallnav::_find( sat, t );
  tmp->print(); 
   
  _gmutex.unlock(); return;
}


// return list of available constellation
// ----------
set<GSYS> t_gallnav::systems()
{
  gtrace("t_gallnav::systems");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  set<GSYS> all_sys;
  auto itPRN = _mapsat.begin();
   
  while( itPRN != _mapsat.end() ){
    GSYS gsys = t_gsys::str2gsys(itPRN->first.substr(0, 1));
    if( all_sys.find(gsys) == all_sys.end() ) all_sys.insert(gsys);
    itPRN++;
  }
  _gmutex.unlock(); return all_sys;
}

// return list of available satellites
// ----------
set<string> t_gallnav::satellites()
{
  gtrace("t_gallnav::satellites");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  set<string> all_sat;
  auto itPRN = _mapsat.begin();
   
  while( itPRN != _mapsat.end() ){
    if( all_sat.find(itPRN->first) == all_sat.end() ) all_sat.insert(itPRN->first);
    itPRN++;
  }
  _gmutex.unlock(); return all_sat;
}


// return first position for satellite
// ----------
t_gtime t_gallnav::beg_gnav( string prn )
{
 gtrace("t_gallnav::beg_gnav");
   
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = LAST_TIME;

  if( ! prn.empty() ){

    if( _mapsat.find(prn) != _mapsat.end() && _mapsat[prn].size() > 0 ){
      tmp = _mapsat[prn].begin()->first;
    }
    
  }else{
    for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
        if( _mapsat[itSAT->first].begin()->first < tmp ){
          tmp = _mapsat[itSAT->first].begin()->first;
        }
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// return last position for satellite
// ----------
t_gtime t_gallnav::end_gnav( string prn )
{
 gtrace("t_gallnav::end_gnav");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime tmp = FIRST_TIME;

  if( ! prn.empty() ){

    if( _mapsat.find(prn) != _mapsat.end() &&  _mapsat[prn].size() > 0  ){
      tmp = _mapsat[prn].rbegin()->first;
    }

  }else{
    for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT ){
      for( auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it ){
        if( _mapsat[itSAT->first].rbegin()->first > tmp ){ 
          tmp = _mapsat[itSAT->first].rbegin()->first; 
        }
      }
    }
  }

  _gmutex.unlock(); return tmp;
}


// add navigation message
// ----------
int t_gallnav::add( shared_ptr<t_gnav> nav )
{   
  gtrace("t_gallnav::add");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  t_gtime ep(nav->epoch());
  string sat = nav->sat();
  
  if( _chk_navig ){
    set<string> msg; nav->chk(msg);
    if( _log ){
      for(auto it = msg.begin(); it != msg.end(); it++)
        _log->comment(2, "gallnav", *it);
    }
  }

  // test navigation type
  bool add = false;
  if( _multimap                                        ){ add = true; }    // multimap
  else if( _mapsat[sat].find(ep) == _mapsat[sat].end() ){ add = true; }    // non-existent
  else if( nav->id_type() == t_gdata::EPHGAL           ){                  // check nav-type
    auto itB = _mapsat[sat].lower_bound(ep);
    auto itE = _mapsat[sat].upper_bound(ep);
    while( itB != itE ){
      if( dynamic_pointer_cast<t_gnav>(itB->second)->gnavtype() == nav->gnavtype() ){ 
        add = false; break; } // exclude the message and skip!
      else{ add = true;         } // ok
      ++itB;
    }
  }
  if(!nav->valid()) add = false;    // validity test  
  
  if( add ){
    if( _log && _log->verb() >= 3 ){
      ostringstream lg;
      lg << "add sat [" << nav->str_type() << "]: " << sat << " " << ep.str_ymdhms()
         << " iod: " << fixed << setw(3) <<  nav->iod()
         << " flg: " << fixed << setw(3) <<  nav->healthy();
      _log->comment(3,"gallnav",lg.str());
    }
    _mapsat[sat].insert(make_pair(ep,nav));
    
  }else if( _log && _log->verb() >= 3 ){
    ostringstream lg;
    lg << "skip sat [" << nav->str_type() << "]: " << sat << " " << ep.str_ymdhms()
       << " iod: " << fixed << setw(3) <<  nav->iod()
       << " flg: " << fixed << setw(3) <<  nav->healthy();
    _log->comment(3,"gallnav",lg.str());
  }	

  _gmutex.unlock(); return 0;
}


// return number of epochs
// ----------
unsigned int t_gallnav::nepochs( const string& prn )
{
 gtrace("t_gallnav::nepochs");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  unsigned int tmp = 0;
  if( _mapsat.find(prn) != _mapsat.end() ) tmp = _mapsat[prn].size();
   
  _gmutex.unlock(); return tmp;
}


// list of epochs
// ----------
vector<t_gtime> t_gallnav::vec_epo( string prn )
{
 gtrace("t_gallnav::list_epo");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  vector<t_gtime> v_epo;

  auto itFIRST = _mapsat.begin();
  auto itLAST  = _mapsat.end();

  if( !prn.empty() ) itLAST++ = itFIRST = _mapsat.find( prn );

#ifdef DEBUG   
  if( itFIRST == itLAST          ) cout << "FIRST = LAST  \n";
  if( itFIRST == _mapsat.begin() ) cout << "FIRST = begin \n";
  if( itLAST  == _mapsat.end()   ) cout << "FIRST = end   \n";
#endif

  for( auto itSAT = itFIRST; itSAT != itLAST; ++itSAT )
    for( auto itEPO = itSAT->second.begin(); itEPO != itSAT->second.end(); ++itEPO )
      v_epo.push_back( itEPO->first );
   
  _gmutex.unlock(); return v_epo;
}

// list of epochs
// ----------
set<t_gtime> t_gallnav::vec_epo( GSYS gsys )
{
  _gmutex.lock();
   
  set<t_gtime> v_epo;

  auto itFIRST = _mapsat.begin();
  auto itLAST  = _mapsat.end();  

#ifdef DEBUG   
  if( itFIRST == itLAST          ) cout << "FIRST = LAST  \n";
  if( itFIRST == _mapsat.begin() ) cout << "FIRST = begin \n";
  if( itLAST  == _mapsat.end()   ) cout << "FIRST = end   \n";
#endif

  for( auto itSAT = itFIRST; itSAT != itLAST; ++itSAT ){
    string prn = itSAT->first;
    GSYS gs = t_gsys::str2gsys(prn.substr(0, 1));
    if (gs != gsys)
      continue;
    for ( auto itEPO = itSAT->second.begin();
         itEPO != itSAT->second.end(); ++itEPO) {
      v_epo.insert(itEPO->first);
    }
  }

  _gmutex.unlock(); return v_epo;
}

// list of nav messages
// ----------
vector<shared_ptr<t_geph>> t_gallnav::vec_nav( string prn )
{
 gtrace("t_gallnav::list_nav");

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<shared_ptr<t_geph>> v_nav;

  auto itFIRST = _mapsat.begin();
  auto itLAST  = _mapsat.end();

  if( !prn.empty() ) itLAST++ = itFIRST = _mapsat.find( prn );

  for( auto itSAT = itFIRST; itSAT != itLAST; ++itSAT ){
    for( auto itEPO = itSAT->second.begin(); itEPO != itSAT->second.end(); ++itEPO ){
      v_nav.push_back( itEPO->second );
    }
  }
   
  _gmutex.unlock(); return v_nav;
}


// list of calculated crd
// ----------
map<string, t_gtriple> t_gallnav::map_xyz(set<string> prns, const t_gtime& epo)
{
  map<string, t_gtriple> m_xyz;

  t_gtriple xyz;

  for ( auto itPRN = prns.begin(); itPRN != prns.end(); itPRN++) {
    string prn = *itPRN;

    double sat_xyz[3] = {0.0, 0.0, 0.0};
    double sat_var[3] = {0.0, 0.0, 0.0};
    double sat_vel[3] = {0.0, 0.0, 0.0};

    if (pos(prn, epo, sat_xyz, sat_var, sat_vel) >= 0) {
      xyz.set(0, sat_xyz[0]);
      xyz.set(1, sat_xyz[1]);
      xyz.set(2, sat_xyz[2]);

      m_xyz[prn] = xyz;

    } else {
      continue;
    }
  }

  return m_xyz;
}   

// list of calculated crd and clk   
map<string, double> t_gallnav::map_clk(set<string> prns, const t_gtime& epo)
{
   map<string, double> m_clk;      
      
   double sat_clk  = 0;
   double var  = 0;
   double dclk = 0;
      
   for( auto itPRN =  prns.begin(); itPRN != prns.end(); itPRN++){
      if( clk(*itPRN, epo, &sat_clk, &var, &dclk) >= 0 ) {
	 m_clk[*itPRN] = sat_clk;
      }
   }	      

   return m_clk;
}    

// list of calculated pos from all redundant navig. messages   
map<string, map<shared_ptr<t_geph>, t_gtriple> > t_gallnav::multi_xyz(set<string> prns, const t_gtime& epo)
{
   map<string, map<shared_ptr<t_geph>, t_gtriple>> m_xyz;
      
   double sat_xyz[3] = {0.0, 0.0, 0.0};
   double sat_var[3] = {0.0, 0.0, 0.0};
   double sat_vel[3] = {0.0, 0.0, 0.0};   
   t_gtriple xyz;
   
   for( auto itPRN =  prns.begin(); itPRN != prns.end(); itPRN++) {
      vector<shared_ptr<t_geph>> vec_geph = this->_find_mult(*itPRN, epo);

      for( auto itEPH = vec_geph.begin(); itEPH != vec_geph.end(); itEPH++) {
	 if( *itEPH && (*itEPH)->pos(epo, sat_xyz, sat_var, sat_vel) >= 0 ) {	      
	    xyz.set(0, sat_xyz[0]);
	    xyz.set(1, sat_xyz[1]);
	    xyz.set(2, sat_xyz[2]);
	    m_xyz[*itPRN][*itEPH] = xyz;
	 }else continue;
      
      }
                  
   }	      
   
   return m_xyz;   
}
   
// list of calculated clk from all redundant navig. messages   
map<string, map<shared_ptr<t_geph>, double> > t_gallnav::multi_clk(set<string> prns, const t_gtime& epo)
{
   map<string, map<shared_ptr<t_geph>, double>> m_clk;      
      
   double sat_clk  = 0;
   double var  = 0;
   double dclk = 0;
   
   for( auto itPRN =  prns.begin(); itPRN != prns.end(); itPRN++) {
      vector<shared_ptr<t_geph>> vec_geph = this->_find_mult(*itPRN, epo);

      for( auto itEPH = vec_geph.begin(); itEPH != vec_geph.end(); itEPH++) {
        if( *itEPH && (*itEPH)->clk(epo, &sat_clk, &var, &dclk) >= 0 ) m_clk[*itPRN][*itEPH] = sat_clk;
      }
                  
   }	      
   
   return m_clk;
}
   

// clean invalid messages
// ----------
void t_gallnav::clean_invalid()
{
  gtrace("t_gallnav::clean_invalid");

  _gmutex.lock();

  for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT ) { 
    string prn = itSAT->first;
     
    auto itEPO = itSAT->second.begin();     
    while( itEPO != itSAT->second.end() ){
      if( ! itEPO->second->valid() ){
        if( _log ) _log->comment(2,"gallnav","Del NAV invalid: "+prn
				 +" "+itEPO->second->epoch().str_ymdhms()
				 +" "+itEPO->second->gio()->path());
        _mapsat[prn].erase(itEPO++);
      }else{ itEPO++; }
    }
  }
  
  _gmutex.unlock();
}

// clean redundant messages
// ----------
void t_gallnav::clean_duplicit()
{
  gtrace("t_gallnav::clean_duplicit");

  _gmutex.lock();

  for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT )
  { 
    string prn = itSAT->first;
    map<t_gtime,set<int>> list_epo;

    for( auto itEPO  = itSAT->second.begin();
              itEPO != itSAT->second.end(); )
    {
      t_gtime epo = itEPO->first;
      int src = itEPO->second->src();
      
      if( list_epo.find(epo)      == list_epo.end()  ||
          list_epo[epo].find(src) == list_epo[epo].end()
      ){
        list_epo[epo].insert( src );
        if( _log ) _log->comment(2,"gallnav","OK  NAV unique  : "+prn
		  	       +" "+itEPO->second->epoch().str_ymdhms()
		  	       +" "+int2str(itEPO->second->src())
			       +" "+itEPO->second->gio()->path());
        ++itEPO;
      }else{ // remove redundant
        if( _log ) _log->comment(2,"gallnav","Del NAV multiple: "+prn
		  	       +" "+itEPO->second->epoch().str_ymdhms()
		  	       +" "+int2str(itEPO->second->src())
			       +" "+itEPO->second->gio()->path());

        _mapsat[prn].erase(itEPO++);
      }
    }
  }

  _gmutex.unlock();
}


// clean function
// ----------
void t_gallnav::clean_outer( const t_gtime& beg, const t_gtime& end )
{
 gtrace("t_gallnav::clean_outer");

  if( end < beg ) return;
  if( beg == FIRST_TIME ) return;
  if( end ==  LAST_TIME ) return;

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

#ifdef DEBUG   
  if( _log ) _log->comment(1,"gallnav","NAV cleaned interval: "+beg.str_ymdhms()+" - "+end.str_ymdhms());
#endif

  // loop over all satellites
  auto itPRN = _mapsat.begin();
  while( itPRN != _mapsat.end() ){
    string prn = itPRN->first;

    // find and CLEAN all data (epochs) out of the specified period !
    auto itFirst = _mapsat[prn].begin();
    auto itLast  = _mapsat[prn].end();
    auto itBeg   = _mapsat[prn].lower_bound(beg);  // greater only   old: // greater|equal
    auto itEnd   = _mapsat[prn].upper_bound(end);  // greater only!
     
//  cout << "distance = " << distance(itBeg,itEnd) << " prn " << prn << endl;
   
    // remove before BEGIN request
    if( itBeg != itFirst ){
      if( _log && _log->verb() >= 2 ){
        _log->comment(2,"gallnav","Del NAV before: "+itFirst->first.str_ymdhms()+" "+prn);
      }       
      auto it = itFirst;

      while( it != itBeg && it != itLast){
	if( (it->second).use_count() == 1 ){ _mapsat[prn].erase(it++); }else{ it++; }
      }       
    }
     
    // remove after END request
    if( itEnd != itLast ){
      if( _log && _log->verb() >= 2 ){
        _log->comment(2,"gallnav","Del NAV after : "+itEnd->first.str_ymdhms()+" "+prn );
      }
      auto it = itEnd;

      while( it != itLast ){
        if( (it->second).use_count() == 1 ){ _mapsat[prn].erase(it++); }else{ it++; }
      }
    }
    itPRN++;
  }
#ifdef BMUTEX   
  lock.unlock();
#endif
  _gmutex.unlock();
}


// get number of satellites
// ----------
int t_gallnav::nsat(GSYS gs)
{
  gtrace("t_gallnav::nsat(GSYS)");

  int nsatell = 0;
  for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT )
  {
    string sat = itSAT->first;
    if( t_gsys::char2gsys(sat[0]) == gs || gs == GNS ) nsatell++;
  }   
  return nsatell;
}


// get interval between messages   
// ----------
int t_gallnav::intv(GSYS gs)
{
  gtrace("t_gallnav::intv(GSYS)");
   
  int interval = 0;  if( gs == GNS ) return interval;

  // collect sampling intervals
  map<int,int> m_intv;
  for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT )
  {
    string sat = itSAT->first; if( t_gsys::char2gsys(sat[0]) != gs ) continue;

    t_gtime lstEPO;
    for( auto itEPO = _mapsat[sat].begin(); itEPO != _mapsat[sat].end(); ++itEPO )
    {
      if( itEPO != _mapsat[sat].begin() ){
        int diff = int(floor(itEPO->first - lstEPO));
        m_intv[diff]++;
      }
      lstEPO = itEPO->first; 
    }
  }

  // auto-detect sampling interval
  int count = 0;
  for( auto it = m_intv.begin(); it != m_intv.end(); ++it ){
	
//    cout << " " << it->first << ":" << it->second;
    if( it->second > count ){
      interval = it->first;
      count    = it->second;
    }
//    cout << endl;
  }
//  cout << t_gsys::gsys2str(gs) << " : " << interval << endl;

  return interval;
}
   

// get existing number of messages
// ----------
int t_gallnav::have(GSYS gs, const t_gtime& beg, const t_gtime& end)
{
  gtrace("t_gallnav::have(GSYS,beg,end)");

  int existing = 0;
  for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT )
  {
    string sat = itSAT->first; if( t_gsys::char2gsys(sat[0]) != gs ) continue;

    for( auto itEPO = _mapsat[sat].begin(); itEPO != _mapsat[sat].end(); ++itEPO ){
      if( itEPO->first < beg-900 || itEPO->first > end+900 ) continue;
      else existing++;
    }
  }   
  return existing;
}


// get expected number of messages
// ----------
int t_gallnav::expt(GSYS gs, const t_gtime& beg, const t_gtime& end)
{
  gtrace("t_gallnav::expt(GSYS,beg,end)");

  if(end < beg) return 0;
  int xint = intv(gs);
  int xsat = nsat(gs);
  int diff = int(floor(end-beg));
  if( xint == 0 || xsat == 0 || diff == 0.0 ) return 0;
  return int(floor( diff/xint * xsat ));
}


// get excluded number of messages
// ----------
int t_gallnav::excl(GSYS gs, const t_gtime& beg, const t_gtime& end)
{
  gtrace("t_gallnav::excl(GSYS,beg,end)");
  int exclude = 0;
  return exclude;
}

   
// return list of available satellites
// ----------
int t_gallnav::consolidate( double cfdi )
{
  gtrace("t_gallnav::consolidate");

  if( cfdi <= 0.0 ) cfdi = 10.0; // DEFAULT

  map<shared_ptr<t_geph>,double> m_penalty;

  t_gtime t_first = FIRST_TIME;
  t_gtime t_last  = LAST_TIME;

  for( auto itSAT = _mapsat.begin(); itSAT != _mapsat.end(); ++itSAT )
  {
    string prn = itSAT->first;

    NAVDATA str_KPL[] = { NAV_A, NAV_E, NAV_M, NAV_I, NAV_IDOT, NAV_OMEGA, NAV_OMG, NAV_OMGDOT, NAV_DN, NAV_F0 };
    NAVDATA str_XYZ[] = { NAV_X, NAV_XD, NAV_XDD, NAV_Y, NAV_YD, NAV_YDD, NAV_Z, NAV_ZD, NAV_ZDD };

    vector<NAVDATA> param;
    if( prn[0] == 'R' || prn[0] == 'S' ){
          vector<NAVDATA> xxx( str_XYZ, str_XYZ +  9); param = xxx; }
    else{ vector<NAVDATA> xxx( str_KPL, str_KPL + 10); param = xxx; }

    for( auto itPAR = param.begin(); itPAR != param.end(); ++itPAR )
    {       
      vector<shared_ptr<t_geph>> v_ptr;
      vector<shared_ptr<t_geph>>::iterator itPTR;
       
      vector<t_gtime> v_tim;
      vector<t_gtime>::const_iterator itTIM;

      vector<double> v_val, v_timdif, v_valdif;
      vector<double>::const_iterator itVAL,itDIF;

      // prepare values & differences vectors for statistics
      for( auto itEPH = _mapsat[prn].begin(); itEPH != _mapsat[prn].end(); ++itEPH )
      {       
        t_timdbl tmp = itEPH->second->param(*itPAR);
        double mindif = 60.0; // minimum timedif applied for differences
	 
        if( !v_tim.empty() ){ // itEPH != _mapsat[prn].begin() ){
          // check selected parameters (differences) for periodic switch
          if( itEPH->second->param_cyclic( *itPAR ) ){
            if( tmp.second - v_val.back() >  G_PI ) tmp.second -= 2*G_PI;
            if( tmp.second - v_val.back() < -G_PI ) tmp.second += 2*G_PI;
            mindif = 1.0; // periodic
          }
          v_timdif.push_back( tmp.first  - v_tim.back() );
          v_valdif.push_back( fabs(v_timdif.back()) > mindif // relax a short time difference to a minute
                              ? (tmp.second - v_val.back()) / (tmp.first - v_tim.back())
                              : (tmp.second - v_val.back()) / mindif
                            );
        }
        
        v_ptr.push_back( itEPH->second ); // save pointer
        v_tim.push_back( tmp.first     );
        v_val.push_back( tmp.second    );

#ifdef DEBUG
        cout << fixed << setprecision(7)
             << "prn: " << prn
             << " " << v_tim.back().str_ymdhms()
             << " " << setw(12) << v_val.back();
        if( itEPH != _mapsat[prn].begin() ){ 
          cout << " " << setw(12) << v_valdif.back()
               << " " << setw(12) << v_timdif.back();
        }
        cout << endl;
#endif
      }

      if( v_val.size() <= 1 ){ continue; } // no values for v_val or v_valdif

      // calculate values & differences statistics
      t_gstat stt_val(v_val);   stt_val.calc_stat(); stt_val.calc_median(); stt_val.calc_minmax();
      t_gstat stt_dif(v_valdif);stt_dif.calc_stat(); stt_dif.calc_median(); stt_dif.calc_minmax();

      ostringstream oSTT;
      oSTT << fixed << setprecision(5)
           << prn << " stat: " << setw(12) << "parameter"          << setw(12) << "difference"         << endl
           << prn << " medi: " << setw(12) << stt_val.get_median() << setw(12) << stt_dif.get_median() << endl
           << prn << " mean: " << setw(12) << stt_val.get_mean()   << setw(12) << stt_dif.get_mean()   << endl
           << prn << " sdev: " << setw(12) << stt_val.get_std()    << setw(12) << stt_dif.get_std()    << endl
           << prn << " rms : " << setw(12) << stt_val.get_rms()    << setw(12) << stt_dif.get_rms()    << endl
           << prn << " min : " << setw(12) << stt_val.get_min()    << setw(12) << stt_dif.get_min()    << endl
           << prn << " max : " << setw(12) << stt_val.get_max()    << setw(12) << stt_dif.get_max()    << endl;

      if(_log) _log->comment(-3,oSTT.str());

      for( itPTR  = v_ptr.begin(),itVAL  = v_val.begin(),itTIM  = v_tim.begin(),itDIF  = v_valdif.begin();
           itPTR != v_ptr.end(),  itVAL != v_val.end();
           ++itPTR,             ++itVAL,               ++itTIM,				  ++itDIF
      ){
        bool badV = ( fabs( *itVAL - stt_val.get_median() ) > fabs(cfdi * stt_val.get_std()) && stt_val.get_std() > 0 );
        bool badD = false;
        if( itDIF != v_valdif.end()){		 
			badD = ( fabs( *itDIF - stt_dif.get_median() ) > fabs(cfdi * stt_dif.get_std()) && stt_dif.get_std() > 0 );
		}
		else {
			itDIF--;
		}
        if( itPTR == v_ptr.begin() ){

          ostringstream os;
          os << fixed << setprecision(9) << prn << itTIM->str_ymdhms(" ")
             << "   rmsV: " << setw(20) << cfdi * stt_val.get_std()
             << "   rmsD: " << setw(20) << cfdi * stt_dif.get_std();
          if( _log ) _log->comment(-3,os.str());
        }

        ostringstream os;
        os << fixed << setprecision(6) 
           << prn << itTIM->str_ymdhms(" ")
           <<      " " << setw( 2) << *itPAR
           << " val: " << setw(18) << *itVAL
           <<      " " << setw(18) << *itVAL - stt_val.get_median()
	         <<" badO: " << setw( 1) << badV;

        if( itDIF != v_valdif.end() ){
          os << setprecision(12)
             << " dif: " << setw(20) << *itDIF
             <<      " " << setw(20) << *itDIF - stt_dif.get_median()
  	         <<" badD: " << setw( 1) << badD;
        }else{ os << endl; }
        if( _log ) _log->comment(-3,os.str());
        
        if( badV ){ m_penalty[*itPTR]     += 1.00; }
        if( badD ){ m_penalty[*itPTR]     += 0.25;
        if((itPTR + 1) != v_ptr.end())  m_penalty[*(itPTR+1)] += 0.25; }
      }
    }
  }

  for(auto it = m_penalty.begin(); it != m_penalty.end(); ++it ){
    if( it->second > 3 ){ // criteria
      shared_ptr<t_geph> itP = it->first; itP->valid(false); // SET INVALID

      if( _log ) _log->comment(1,"gallnav","Set NAV invalid: "+itP->sat()+" "+itP->epoch().str_ymdhms()
			       +" [penalty"+dbl2str(it->second,1)+"] "+itP->gio()->path());
    }
  }

  return 0;
}


// return list of available satellites
// ----------
shared_ptr<t_geph> t_gallnav::_find( string sat, const t_gtime& t, bool chk_mask )
{
  gtrace("t_gallnav::_find sat/time");

  if( _mapsat.find(sat) == _mapsat.end() ) return _null; // make_shared<t_geph>();

  GSYS gnss = t_gsys::str2gsys(sat.substr(0,1));

  if( _mapsat[sat].size() == 0 ){
    if( _log ) _log->comment(2,"gallnav", sat+" no gnav elements found");
    return _null;
  }
  
  auto it = _mapsat[sat].lower_bound(t);  // greater|equal  (can be still end())
  if( it == _mapsat[sat].end() ) it--;                   // size() > 0 already checked above

  double maxdiff = t_gnav::nav_validity(gnss)*1.1; 

  if (gnss != GLO){
    
    // 2 conditions for all cases when NAV is in fugure and the iterator should be moved back (it--)
    if( fabs(t - it->second->epoch()) > maxdiff ||  // too far navigation message in future!
        (_chk_tot && !it->second->chktot(t))        // check ToT for using past messages only
    ){
      if(_mapsat[sat].size() > 0 && it != _mapsat[sat].begin()) it--; // one more step back
    }

  }else{    
    t_gtime tt = t - maxdiff;  // time span of Glonass ephemerides is 15 min
    
    auto itt = _mapsat[sat].lower_bound(tt);  // greater|equal
    if( itt == _mapsat[sat].end() ){                         // size() > 0 already checked above
      if( _log ) _log->comment(2,"gallnav", sat+" gnav element not found: " + t.str_ymdhms() );
      return _null;
    }

    double dt1 = itt->first - t;
    double dt2 =  it->first - t;
    if( fabs(dt1) < fabs(dt2) ) it = itt;
  }

  if(fabs(t - it->second->epoch()) > maxdiff ||
     (_chk_tot && !it->second->chktot(t)) ) {   // simulation real-time epoch search 
    if(_log){ 
      string lg(sat+" gnav element not found: " + t.str_ymdhms());
      _log->comment(2,"gallnav",lg);
    }
    return _null;
  }

  if(_chk_health && chk_mask){
    if(it->second->healthy()) return it->second;
    else                      return _null;
  }
  return it->second;
  
}   
   
// return list of available satellites
// ----------
vector<shared_ptr<t_geph>> t_gallnav::_find_mult( string sat, const t_gtime& t )
{
  gtrace("t_gallnav::_find vec sat/time");

  vector<shared_ptr<t_geph>> vec_geph;
   
  if( _mapsat.find(sat) == _mapsat.end() ) {
    vec_geph.push_back(_null);
    return vec_geph;
  }
   
  GSYS gnss = t_gsys::str2gsys(sat.substr(0,1));

  auto it = _mapsat[sat].lower_bound(t);  // greater|equal

  double maxdiff = t_gnav::nav_validity(gnss)*1.1;  

  if (gnss != GLO){
    if( it == _mapsat[sat].end() ){
      t_gtime tt = t - maxdiff;
      auto itt = _mapsat[sat].lower_bound(tt);  // greater|equal
      if( itt == _mapsat[sat].end() ){       
        if(_log){ string lg(sat+" gnav element not found: "+ t.str_ymdhms());
          _log->comment(3,"gallnav",lg);
        }       
        vec_geph.push_back(_null); return vec_geph;
      } else it = itt;
    }
  }else{    
    t_gtime tt = t - maxdiff;  // time span of Glonass ephemerides is 15 min
    auto itt = _mapsat[sat].lower_bound(tt);  // greater|equal     
    if( it == _mapsat[sat].end() ){
      if( itt == _mapsat[sat].end() ) { 
        vec_geph.push_back(_null); return vec_geph;
      }else it = itt;
    }else{
      if( itt == _mapsat[sat].end() ) { 
        vec_geph.push_back(_null); return vec_geph;
      }
      double dt1 = itt->first - t;
      double dt2 =  it->first - t;
      if ( fabs(dt1) < fabs(dt2) ) it = itt;
    }
  }
  
  t_gtime epo = it->first;
  unsigned int cnt = _mapsat[sat].count(epo);

  if(cnt >= 1) {    
    for( auto itMULT = _mapsat[sat].equal_range(epo).first; itMULT != _mapsat[sat].equal_range(epo).second; ++itMULT) {      
      vec_geph.push_back(itMULT->second);
    }
  }
  
  return vec_geph;
}   
   
   
// return list of available satellites
// ----------
shared_ptr<t_geph> t_gallnav::_find( string sat, int iod, const t_gtime& t )
{
  gtrace("t_gallnav::_find sat/iod");

  if( _mapsat.find(sat) == _mapsat.end() ) return _null; // make_shared<t_geph>();

  auto it = _mapsat[sat].begin();
  while( it != _mapsat[sat].end() ){
    shared_ptr<t_gnav> pt_nav;
    if( ( pt_nav = dynamic_pointer_cast<t_gnav>(it->second) ) != NULL ){
//       cout << "nav->epoch() " << pt_nav->epoch().str_ymdhms() << " nav->iod() = " << pt_nav->iod()  << " rtcm iod: " << iod << endl;
      if( pt_nav->iod() == iod && abs(pt_nav->epoch().diff(t)) <= MAX_GPS_TIMEDIFF ) {
        break; // found
      }
    }
    it++;
  }

  // not found !
  if( it == _mapsat[sat].end() ){
    if( _log && _log->verb() >= 3 ){
      ostringstream lg;
      lg << sat + " gnav element not found "
         << " [iod: " << iod << "]";
      _log->comment(3,"gallnav",lg.str());
    }
    return _null; // make_shared<t_geph>();
  }

  return it->second;
}

} // namespace
