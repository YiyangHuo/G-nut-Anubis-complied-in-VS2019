
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

#include "gall/gallprod.h" 
#include "gprod/gprodcrd.h"
#include "gprod/gprodtrp.h"
#include "gprod/gprodion.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallprod::t_gallprod( )
  : t_gdata()
{
  gtrace("t_gallprod::constructor");
  id_type(  t_gdata::ALLPROD);
  id_group( t_gdata::GRP_PRODUCT );
}


// destructor
// ----------
t_gallprod::~t_gallprod()
{
  gtrace("t_gallprod::destructor");
  this->clear();
}


// clean
// --------------------
void t_gallprod::clear()
{
  gtrace("t_gallprod::clear");

  _gmutex.lock();
  _map_prod.clear();
  _gmutex.unlock();
}


// add product
// ---------------------
int t_gallprod::add(shared_ptr<t_gprod> prod, string site) // str only if prod->obj_pt = 0
{
  gtrace("gallprod::add");
   
  ID_TYPE type = prod->id_type();
   
  string id = prod->obj_id(); // if( id.empty() ){ return -1; } // no, use id instead!
   
  if( id.empty() ) id = site;
  if( id.empty() && ! site.empty() && id != site )
    cerr << "warning [gallprod] - adding product with both [and different] non-zero string and object id!\n";
   
  _gmutex.lock();
  _map_prod[id][type][prod->epoch()] = prod;  // allocate in heap

#ifdef DEBUG   
  cout << "adding pt: " <<  prod->epoch().str_ymdhms(" ") << ":" << prod->epoch().dsec() << " type: " << type << endl;
#endif

  _gmutex.unlock(); return 0;
}


// get product
// --------------------
shared_ptr<t_gprod> t_gallprod::get(const string& site, ID_TYPE type, const t_gtime& t)
{
  gtrace("gallprod::get(str)");
   
  _gmutex.lock(); 
   
  shared_ptr<t_gprod> obj_pt = _find(site, type, t);

#ifdef DEBUG
  cout << "geting pt: " << obj_pt << " " << obj_pt.use_count() << " " << _map_prod.size() << " ... "
       << t.str_ymdhms(" ") << ":" << t.dsec() << endl;
#endif

  _gmutex.unlock(); return obj_pt;
}

// list of sites
// --------------------
set<string> t_gallprod::prod_sites()
{
  gtrace("t_gallprod::prod_sites");

  set<string> site_list;
  _gmutex.lock();
  t_map_prd::const_iterator it;  
  for( it = _map_prod.begin(); it != _map_prod.end(); ++it ){
     site_list.insert( it->first );
  }
  _gmutex.unlock(); return site_list;
}
      
// list of product types
// --------------------
set<t_gdata::ID_TYPE> t_gallprod::prod_types(const string& site)
{
  gtrace("t_gallprod::prod_types " + site );

  set<t_gdata::ID_TYPE> type_list;

  _gmutex.lock();

  t_map_id::const_iterator itTYPE;
  t_map_prd::const_iterator it = _map_prod.find( site );

  if( it == _map_prod.end() ) { _gmutex.unlock(); return type_list;}
   
  for( itTYPE = _map_prod[site].begin(); itTYPE != _map_prod[site].end(); ++itTYPE ){
     type_list.insert( itTYPE->first );
  }

  _gmutex.unlock(); return type_list;
}

// list of product types
// --------------------
set<t_gtime> t_gallprod::prod_epochs(const string& site, ID_TYPE type)
{
  gtrace("t_gallprod::prod_epochs " + site);

  set<t_gtime> epo_list;

  _gmutex.lock();

  t_map_epo::const_iterator itEPO;
  t_map_id::const_iterator itTYPE;
  t_map_prd::const_iterator it = _map_prod.find(site);

  if (it == _map_prod.end()) {
    _gmutex.unlock();
    return epo_list;
  } else {
    itTYPE = it->second.find(type);
    if (itTYPE == it->second.end()) {
      _gmutex.unlock();
      return epo_list;
    } else {
      for (itEPO = itTYPE->second.begin(); itEPO != itTYPE->second.end(); ++itEPO) {
        epo_list.insert(itEPO->first);
      }
    }
  }

  _gmutex.unlock();
  return epo_list;
}   

// clean data out of the interval
// ----------------------------   
void t_gallprod::clean_outer( const t_gtime& beg, const t_gtime& end )
{
  gtrace("t_gallprod::clean_outer");

  if( end < beg ) return;

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

#ifdef DEBUG   
  if( _log ) _log->comment(1,"gallprod","products clean request: "
			              + beg.str("%Y-%m-%d %H:%M:%S[%T] - ")
                                      + end.str("%Y-%m-%d %H:%M:%S[%T]"));
#endif      

   // loop over all sites
   t_map_prd::iterator itKEY = _map_prod.begin();
   while( itKEY != _map_prod.end() ){
      string key = itKEY->first;

      t_map_id::iterator itID = itKEY->second.begin();
      while( itID != itKEY->second.end() ){
	 string id = t_gdata::type2str(itID->first);
	 
	 t_map_epo::iterator it;
	 t_map_epo::iterator itFirst = itID->second.begin();
	 t_map_epo::iterator itLast  = itID->second.end();	 
	 t_map_epo::iterator itBeg   = itID->second.lower_bound(beg);
	 t_map_epo::iterator itEnd   = itID->second.upper_bound(end);	 
	 
	 
	 // remove before BEGIN request
	 if( itBeg != itFirst ){
       
            it = itFirst;

	    // begin is last
	    if( itBeg == itLast ){
	       itBeg--;
	       if( _log ) _log->comment(2,"gallprod","products removed before: " + itFirst->first.str("%Y-%m-%d %H:%M:%S[%T] ") + key + " " + id);
	       
  	       while( it != itLast ) // itID->second.erase(it++);
	         if( (it->second).use_count() == 1 ){ itID->second.erase(it++); }
	         else{ it++; }
	       
	    // begin is not last
	    }else{
	       if( _log ) _log->comment(2,"gallprod","products removed before: " + itBeg->first.str("%Y-%m-%d %H:%M:%S[%T] ") + key + " " + id);
	       
    	       while( it != itBeg ) // itID->second.erase(it++);
	         if( (it->second).use_count() == 1 ){ itID->second.erase(it++); }
 	         else{ it++; }		       
	    }
	 }
	 
         // remove after END request
	 if( itEnd != itLast ){
	    if( _log ) _log->comment(2,"gallprod","products removed after : " + itEnd->first.str("%Y-%m-%d %H:%M:%S[%T] ") + key + " " + id);
	    it = itEnd;
	    
	    while( it != itLast ) // itID->second.erase(it++);
	      if( (it->second).use_count() == 1 ){ itID->second.erase(it++); }
	      else{ it++; }
	 }
	 itID++;
      }
      itKEY++;
   }   
   
#ifdef BMUTEX   
  lock.unlock();
#endif
  _gmutex.unlock();
  return; 
}

// remove appropriate element t_gprod*
// ---------------------
void t_gallprod::rem(const string& site, ID_TYPE type, const t_gtime& t)
{
   gtrace("t_gallprod::rem(str)");
  
   _gmutex.lock();

   t_map_epo::iterator itEP;
   t_map_id::iterator itID;
   t_map_prd::iterator it = _map_prod.find( site );
   if( it != _map_prod.end() ){
      itID = it->second.find(type);
      if( itID != it->second.end() ){
	 itEP = itID->second.find(t);
	 if( itEP != itID->second.end() ) _map_prod[site][type].erase(itEP);
      }     
   }

   _gmutex.unlock();
   return;
}

// get ionospheric delay for particular position and epoch
double t_gallprod::iono(const double& lat, const double& lon, const t_gtime& epo, bool interp)
{
  _gmutex.lock();
  double val = 0.0;

  double tec1, tec2;
  tec1 = tec2 = 0.0;

  t_gtime epo1, epo2;
  epo1 = FIRST_TIME;
  epo2 = LAST_TIME;  
  
  shared_ptr<t_gprodion> p_ion1, p_ion2;
  p_ion1 = p_ion2 = nullptr;
  
//  double lon_rot = 0.0; 
  
  t_map_id::iterator itID;
  t_map_prd::iterator it = _map_prod.find("");
  if (it != _map_prod.end()) {
    itID = it->second.find(t_gdata::ION);
    if (itID != it->second.end()) {            

      t_map_epo::iterator itEP2 = itID->second.upper_bound(epo);
      if (itEP2 != itID->second.end()){
        epo2 = itEP2->second->epoch();
        p_ion2 = static_pointer_cast<t_gprodion>( itEP2->second );

      }
      t_map_epo::iterator itEP1 = itEP2;
      if(itEP1 != itID->second.begin()) {                
        itEP1--;
        epo1 = itEP1->second->epoch();
        p_ion1 = static_pointer_cast<t_gprodion>( itEP1->second );

      }
    }
  }

//  lon_rot = lon + OMEGA * (epo - epo1);

  if(p_ion1 != nullptr){
    tec1 = p_ion1->tec(lat, lon);
  }
  if(p_ion2 != nullptr){    
    tec2 = p_ion2->tec(lat, lon);   
  }

  if (interp) {
    if (tec1 > 0 && tec2 > 0 && epo1 != FIRST_TIME && epo2 != LAST_TIME) {
      val = ((epo2 - epo) / (epo2 - epo1)) * tec1 + ((epo - epo1) / (epo2 - epo1)) * tec2;
    } else {
      if (_log) _log->comment(1, "t_gallprod", "Warning: TEC temporal extrapolation not allowed!");
    }
  }else{
    double dt1 = epo - epo1;
    double dt2 = epo2 - epo;
    if(dt1 <= dt2 && tec1 > 0) val = tec1;
    if(dt1 >  dt2 && tec2 > 0) val = tec2;
  }

  #ifdef DEBUG
    cout << "IONEX interpolation: " << epo1.str_ymdhms() << " " << epo2.str_ymdhms() << " " << epo.str_ymdhms() << " " << lat << " " << lon << " " << val << " " << tec1 << " " << tec2 << endl;  
    //int ooo; cin >> ooo;
  #endif

  _gmutex.unlock();
  return val;
}

// find appropriate element t_gprod*
// ---------------------
shared_ptr<t_gprod> t_gallprod::_find(const string& site, ID_TYPE type, const t_gtime& t)
{
  gtrace("t_gallprod::_find(str)");

  t_map_epo::iterator itEP;
  t_map_id::iterator itID;
  t_map_prd::iterator it = _map_prod.find(site);
  if (it != _map_prod.end()) {
    itID = it->second.find(type);
    if (itID != it->second.end()) {
      itEP = itID->second.find(t);
      if (itEP != itID->second.end())
        return itEP->second;
    }
  }

  shared_ptr<t_gprod> tmp;
  return tmp;
}

} // namespace
