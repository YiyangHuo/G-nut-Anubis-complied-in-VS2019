
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

#include "gall/gallobj.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallobj::t_gallobj()
  : _gpcv(0),
    _gotl(0)
{
  gtrace("t_gallobj::constructor");
  id_type(  t_gdata::ALLOBJ );
  id_group( t_gdata::GRP_OBJECT );
  _aloctrn();  
}


// constructor
// ----------
t_gallobj::t_gallobj(t_gallpcv* pcv, t_gallotl* otl)
  : _gpcv(pcv),
    _gotl(otl)
{
  gtrace("t_gallobj::constructor");
  id_type(  t_gdata::ALLOBJ );
  id_group( t_gdata::GRP_OBJECT );
  _aloctrn();
}


// destructor
// ----------
t_gallobj::~t_gallobj()
{
  gtrace("t_gallobj::destructor");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _mapobj.clear();
  _gmutex.unlock(); return;
}

// add object
// ----------
int t_gallobj::add(shared_ptr<t_gobj> obj)
{
  gtrace("t_gallobj::add");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  string s = obj->id();

  // object empty
  if( obj->id().empty() ){ _gmutex.unlock();  return -1; } 

  // object does not exist
  // or _overwrite (--> remove & addnew)
  if( _mapobj.find(s) == _mapobj.end() ){
     
    _mapobj[s] = obj;
     
    if( _log ){
      _mapobj[s]->glog(_log);
      if( _log->verb() >= 1 ) _log->comment(1,"gallobj","add new obj " + s);
    }

  }else{
    // FINALLY COULD ADD RECORDS FROM OBJ to existing OBJ (merge/overwrite mode)
    if( _log && _log->verb() >= 0 ) _log->comment(0,"gallobj","warning - cannot overwrite object: " + s);
    else                                    cerr << "gallobj:  warning - cannot overwrite object: " + s << endl;
    _gmutex.unlock(); return -1;
  }

  // individual get/sync PCV for new object
  // --------------------------------------------------------------------------
  // HOWEVER, IT MUST BE GUARANTED THAT WHENEVER GOBJ IS MODIFIED PCVs are SYNC AGAIN!
  // --------------------------------------------------------------------------

   _mapobj[s]->sync_pcv(_gpcv);
       
   _gmutex.unlock(); return 0;
}

   
// get single object
// -------------------------------
shared_ptr<t_gobj> t_gallobj::obj(string s)
{
  gtrace("t_gallobj::obj");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  shared_ptr<t_gobj> p_obj;
   
  t_map_obj::iterator  it = _mapobj.find(s);
  if( it == _mapobj.end() ){
//    cout << "object " << _mapobj.size() <<  " " << s << " not found !\n";
     if( _log && _log->verb() >= 2 ) _log->comment(2,"gallobj","object: " + s + " not found.");
    _gmutex.unlock(); return p_obj;
  }else{
    p_obj = it->second;
  }

  _gmutex.unlock(); return p_obj;
}

// get object number
// -----------
int t_gallobj::obj_num()
{
   return _mapobj.size();
}

/*
// get all object names
// ----------
vector<string> t_gallobj::objects_id(t_gdata::ID_TYPE id)
{
  gtrace("t_gallobj::objects_id");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<string> all_obj;
  t_map_obj::const_iterator itOBJ = _mapobj.begin();
        
  while( itOBJ != _mapobj.end() ){
    if( id != NONE && id != itOBJ->second->id_type() ){ ++itOBJ; continue; } // filter objects
    all_obj.push_back( itOBJ->first );
    itOBJ++;
  }
   
  _gmutex.unlock(); return all_obj;
}
*/

// get all object elements
// ----------
map<string, shared_ptr<t_gobj>> t_gallobj::objects(t_gdata::ID_TYPE id)
{
  gtrace("t_gallobj::objects");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  map<string,shared_ptr<t_gobj>> all_obj;
  t_map_obj::const_iterator itOBJ = _mapobj.begin();

  while( itOBJ != _mapobj.end() ){
    if( id != NONE && id != itOBJ->second->id_type() ){ ++itOBJ; continue; } // filter objects
    string site = itOBJ->second->id();
    all_obj[ site ] =  itOBJ->second;
    itOBJ++;
  } 

  _gmutex.unlock(); return all_obj;
}


// set PCV
// ------------------------
void t_gallobj::setPCV(t_gallpcv* pcv)
{
  gtrace("t_gallobj::setPCV");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _gpcv = pcv;
  _gmutex.unlock(); return;
}


// set OTL
// ------------------------
void t_gallobj::setOTL(t_gallotl* otl)
{
  gtrace("t_gallobj::setOTL");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _gotl = otl;
  _gmutex.unlock(); return;
}


// get all object elements
// ----------
void t_gallobj::sync_pcvs()
{
  gtrace("t_gallobj::sync_pcvs");
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<shared_ptr<t_gobj>> all_obj;
  t_map_obj::const_iterator itOBJ = _mapobj.begin();
    
  while( itOBJ != _mapobj.end() ){
#ifdef DEBUG     
    cout << "GALLOBJ syncing: " << itOBJ->second->id() << endl;
#endif     
    itOBJ->second->sync_pcv( _gpcv );
    itOBJ++;
  }

#ifdef DEBUG
  cout << " PCVs synchronized for all objects !\n";
#endif
  _gmutex.unlock(); return;
}


// alocate all t_gtrn objects for all GNSS
// ------------------
void  t_gallobj::_aloctrn()
{
   gtrace("t_gallobj::_aloctrn");
   int numsatG = 32;
   int numsatR = 24;
   int numsatE = 30;
   int numsatC = 35;
   int numsatQ =  5;
// int numsatS =  5;   
   
   for (int i = 0; i != GNS; i++){
      GSYS sys = static_cast<GSYS>(i);
      int num;

      switch (sys){
       case GPS : 
	 num = numsatG;
	 break;
       case GAL : 
	 num = numsatE;
	 break;
       case GLO :
	 num = numsatR;
	 break;
       case BDS :
	 num = numsatC;
	 break;	 
       case QZS :
	 num = numsatQ;
	 break;	 
//       case SBS :           // NOT in sequence! FOR THIS CLASS MERGE WITH gnss.h FUNCTIONS!
//	 num = numsatS;
//	 break;	 
       default : 
	 num = 0;
	 break;
      }
      for (int j = 1; j <= num; j++){
 	 string satID;
	 if (j <= 9) satID = string(1,t_gsys::gsys2char(sys)) + "0" + int2str(j);
	 else satID = string(1,t_gsys::gsys2char(sys)) + int2str(j);

	 shared_ptr<t_gtrn> trn_new = make_shared<t_gtrn>();
         trn_new->id( satID );
//         trn_new->ant( satID, FIRST_TIME );    // add antena name
          _mapobj[satID] = trn_new;
      }
   }
     
}

// Read Satelite info file
// ------------------
// !!!!!!! MUST be enhanced for reading whole satelite info file !!!!!!!!!!!!
void t_gallobj::read_satinfo(t_gtime& epo)
{
  gtrace("t_gallobj::read_satinfo");
  t_map_obj::iterator it;
   for (it = _mapobj.begin(); it != _mapobj.end(); it++){
      if (it->second->id_type() == TRN) {
        string ID = it->first;
        it->second->ant(ID, epo);
      }
      
   }
   
}

} // namespace
