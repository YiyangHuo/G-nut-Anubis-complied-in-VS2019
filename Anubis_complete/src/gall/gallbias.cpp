
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

#include "gall/gallbias.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallbias::t_gallbias()
   : _overwrite(false)
{
  id_type(  t_gdata::ALLBIAS );
  id_group( t_gdata::GRP_MODEL);  
}


// destructor
// ----------
t_gallbias::~t_gallbias(){
  _mapbias.clear();
}


// add satellite bias
// ----------
void t_gallbias::add(const string& ac, const t_gtime& epo, const string& obj, t_spt_bias pt_cb)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if(pt_cb == nullptr) {
    _gmutex.unlock();
    return;
  }

  if(pt_cb->ref() == X){     // When ref GOBS is X => bias is expressed in absolute sense
    _mapbias[ac][epo][obj][pt_cb->gobs()] = pt_cb;
  }else{                     // Differential biases need to be splitted
    if(_mapbias[ac][epo][obj].size() == 0){    
      shared_ptr<t_gbias> pt_ref = make_shared<t_gbias>();                            // create first reference bias 
      pt_ref->set(pt_cb->beg(), pt_cb->end(), 0.0, pt_cb->ref(), pt_cb->ref());       // reference bias is set up to zero
      _mapbias[ac][epo][obj][pt_ref->gobs()] = pt_ref;                                // store new bias (reference)
      _mapbias[ac][epo][obj][pt_cb->gobs()] =  pt_cb;                                 // store new bias
    }else{
      t_spt_bias pt_obs1 = _find(ac, epo, obj, pt_cb->gobs());
      t_spt_bias pt_obs2 = _find(ac, epo, obj, pt_cb->ref());
      if(pt_obs1 != nullptr && pt_obs2 == nullptr){         // connection with first signal
        _connect_first(pt_obs1, pt_cb);
        _mapbias[ac][epo][obj][pt_cb->gobs()] = pt_cb;      // store modified bias
      }else if(pt_obs1 == nullptr && pt_obs2 != nullptr){   // connection with second signal
        _connect_second(pt_obs2, pt_cb);
        _mapbias[ac][epo][obj][pt_cb->gobs()] = pt_cb;      // store modified bias
      }else if(pt_obs1 != nullptr && pt_obs2 != nullptr){   // connectin two groups with different reference signal
                                                            // WARNING!!! - this case has not been tested (not happen in tested files)
        // connection with first signal
        _connect_first(pt_obs1, pt_cb);
        // all biases connected with second signal need to be consolidated
        _consolidate(ac, obj, pt_cb, pt_obs2);
      }
    }
  }
  
  _gmutex.unlock(); return;
}

// get single fcb bias for phase
// -------------------------------
double t_gallbias::fcbbias(const string prd, const string& prn, const t_gtime& epo, const GOBS gobs1, const GOBS gobs2){

#ifdef BMUTEX   
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();

	double val = 999.0;
	double valid = 900;	//seconds
	if (prd=="WSB") valid = 24 * 3600;

	auto itAC = _mapbias.find(prd);
	if (itAC != _mapbias.end()) {
		auto itEPO = itAC->second.upper_bound(epo);
		if (itEPO == itAC->second.end() && itAC->second.size() != 0) itEPO--;         // no epochs
		
		auto itEOP0 = itEPO;
		if (itEOP0 != itAC->second.end() && fabs(itEOP0->first - epo) < fabs(itEPO->first - epo)) itEPO = itEOP0;

		if (itEPO != itAC->second.end()) {
			auto itSAT = itEPO->second.find(prn);
			if (itSAT != itEPO->second.end() && itSAT->second.find(gobs1) != itSAT->second.end()
				&& itSAT->second.find(gobs2) != itSAT->second.end()) {
				t_spt_bias pobs1 = itSAT->second.find(gobs1)->second;
				t_spt_bias pobs2 = itSAT->second.find(gobs2)->second;
				val = pobs1->bias(false) - pobs2->bias(false);
			}
		}
	}

	_gmutex.unlock();
	return val;
}

// get single code bias
// -------------------------------
double t_gallbias::get(const t_gtime& epo, const string& obj, const GOBS& gobs1, const GOBS& gobs2){

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  double dcb = 0.0;

  auto itAC = _mapbias.begin();
  string ac = "";
  if(itAC != _mapbias.end()){
    ac = itAC->first;
  }else {
    _gmutex.unlock();
    return dcb;
  }
    
  t_spt_bias pobs1 = _find(ac, epo, obj, gobs1);
  t_spt_bias pobs2 = _find(ac, epo, obj, gobs2);  

  bool found = false;

  if(pobs1 != nullptr && pobs2 != nullptr && pobs1->ref() == pobs2->ref()){
    found = true;
    dcb = pobs1->bias() - pobs2->bias();
  }

  if(!found && _log){
    string dcb_str = "(" + gobs2str(gobs1) + "-" + gobs2str(gobs2) + ")";
    _log->comment(1, "gallbias", "WARNING: code bias " + dcb_str + " not found!");
  }
  
  _gmutex.unlock();
  return dcb;
}
  
// get single code bias
// -------------------------------
t_spt_bias t_gallbias::_find(const string& ac, const t_gtime& epo, const string& obj, const GOBS& gobs)
{
  t_spt_bias pt_bias = nullptr;

  auto itAC = _mapbias.find(ac);
  if(itAC != _mapbias.end()){
    auto itEPO = itAC->second.upper_bound(epo);
    if(itEPO != itAC->second.begin() && itEPO != itAC->second.begin()) itEPO--;  // between epochs
    if(itEPO == itAC->second.end() && itAC->second.size() != 0) itEPO--;         // no epochs
    
    if(itEPO != itAC->second.end()){
      auto itOBJ = itEPO->second.find(obj);
      if(itOBJ != itEPO->second.end()){
        auto itGOBS = itOBJ->second.find(gobs);
        if(itGOBS != itOBJ->second.end()){
          if(itGOBS->second->valid(epo)){
              pt_bias = itGOBS->second;
          }
        } 
      } 
    }
  }
  
  return pt_bias;
}

// get all biases with particular reference singal
// -------------------------------
vector<t_spt_bias> t_gallbias::_find_ref(const string& ac, const t_gtime& epo, const string& obj, const GOBS& ref)
{
  vector<t_spt_bias> vec_bias;

  auto itAC = _mapbias.find(ac);
  if(itAC != _mapbias.end()){
    auto itEPO = itAC->second.find(epo);
    if(itEPO != itAC->second.end()){
      auto itOBJ = itEPO->second.find(obj);
      if(itOBJ != itEPO->second.end()){
        for(auto itGOBS = itOBJ->second.begin(); itGOBS != itOBJ->second.end(); itGOBS++){
          if(itGOBS->second->ref() == ref) vec_bias.push_back(itGOBS->second);
        }
      } 
    } 
  }
  
  return vec_bias;
}
  
// Connect DCB pt_cb2 with first GOBS
void t_gallbias::_connect_first(const t_spt_bias pt_cb1, t_spt_bias pt_cb2)
{
  double newval = pt_cb1->bias(false) - pt_cb2->bias(false);
  pt_cb2->set(newval, pt_cb2->ref(), pt_cb1->ref());
}  

// Connect DCB pt_cb2 with second GOBS
void t_gallbias::_connect_second(const t_spt_bias pt_cb1, t_spt_bias pt_cb2)
{
  double newval = pt_cb1->bias(false) + pt_cb2->bias(false);
  pt_cb2->set(newval, pt_cb2->gobs(), pt_cb1->ref());
}

// Consolidate all biases with reference signal of pt_cb2
void t_gallbias::_consolidate(const string& ac, const string& obj, const t_spt_bias pt_cb1, t_spt_bias pt_cb2)
{
  double diff = pt_cb2->ref() - pt_cb1->ref();
  t_gtime epo = pt_cb1->beg();
  vector<t_spt_bias> vec = _find_ref(ac, epo, obj, pt_cb2->ref());

  for(auto itSPT = vec.begin(); itSPT != vec.end(); itSPT++){
    double newval = (*itSPT)->bias() + diff;
    GOBS   gobs   = (*itSPT)->gobs();
    GOBS   newref = pt_cb1->ref();
    (*itSPT)->set(newval, gobs, newref);
  }
  
}
  
} // namespace
