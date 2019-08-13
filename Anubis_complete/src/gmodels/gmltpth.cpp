
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

#include <iostream>
#include <iomanip>

#include "gmodels/gmltpth.h"
#include "gutils/gnss.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_gmltpth::t_gmltpth(t_gallobs* gobs, t_gsetbase* gset, string site)
{
   _gobs = gobs;
   _set  = gset;
   _site = site;
   
   _sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();
   _mp_limit = 5;

}


// destructor
// ----------
t_gmltpth::~t_gmltpth(){}


// calculate multipath
// ----------------
void t_gmltpth::multipath(const string& prn, GSYS gsys, const t_gtime& epo,
		       const unsigned int& nepo, t_map_mpobs& m_obs)
{
  gtrace("t_gmltpth::multipath");

  string gs = t_gsys::gsys2str(gsys);
  string ep = epo.str_ymdhms();   
     
  // get all observations for the satellite
  vector<shared_ptr<t_gobsgnss>> vobs_prn = _gobs->obs_prn_pt(_site, prn, epo, epo + 2*nepo*_sampling );
  if( vobs_prn.size() == 0 ) return;                    // nothing to do

  // L-bands selection for GNSS (priority)
  t_gobs code;
//  t_gobs gbL1(TYPE_L, BAND, ATTR);
//  t_gobs gbL2(TYPE_L, BAND, ATTR);

//  gbL1.band( t_gsys::band_priority(gsys,0) );           // get first band
//  gbL2.band( t_gsys::band_priority(gsys,1) );           // get second band

  // NEED ONLY ONCE! ===> SIGNAL/ATTR SELECTION BY QUANTITATIVE PRIORITY 
  // --> can be done at _list_mult + saved in member data!
//  gbL1.attr( _gobs->select_attr( _stat_obs, gsys, TYPE_L, gbL1.band() ) );
//  gbL2.attr( _gobs->select_attr( _stat_obs, gsys, TYPE_L, gbL2.band() ) );
   
  vector<double> vMP;
  vector<double>::const_iterator itV;
  double obs, obs_0, mean, rms, difMP;
  double rmsMP = 1.0; // 1m a priori rms
  double sdif2 = rmsMP*rmsMP;
  unsigned int count;
  bool firstCOD = true;

  //select phase GOBS
  t_map_mpobs::iterator itGOBS;  
  set<GOBS> phases;   
  for( itGOBS = m_obs.begin(); itGOBS != m_obs.end(); ++itGOBS )
  {
     GOBS go = itGOBS->first;
     if( gobs_phase(go) ) phases.insert(go);
  }  
  if(phases.size() < 2) {     
     return;
  }  
  t_gobs gbL1( *phases.begin()  );
  t_gobs gbL2( *phases.rbegin() );   
   
   
  // loop GOBS
//  t_gallobs::t_map_stt_obs::const_iterator itGOBS;
  for( itGOBS = m_obs.begin(); itGOBS != m_obs.end(); ++itGOBS )
  {
    GOBS gobs  = itGOBS->first;
     
    if( ! gobs_code(gobs) ||           // code only
        vobs_prn.size() == 0 ){ continue; }       // nothing to do (data removed due to cycle-slips)

    string go  = gobs2str( gobs );
    code.gobs(gobs);

#ifdef DEBUG
//    if( nepo == DAY_EPOCHS )	  
    cout << ep
         << "  "    << prn
         << "  go:" << go 
         << "  b1:" << gobsband2str(gbL1.band()) 
	 << "  b2:" << gobsband2str(gbL2.band())
         << "  a1:" << gobsattr2str(gbL1.attr()) 
	 << "  a2:" << gobsattr2str(gbL2.attr())
         << "  o1:" << gobs2str( gbL1.gobs() )
         << "  o2:" << gobs2str( gbL2.gobs() )	    
         << "  o1:" << (*vobs_prn.begin())->obs_L(gbL1)
         << "  o2:" << (*vobs_prn.begin())->obs_L(gbL2)
         << endl; cout.flush();
#endif

    // rmsMP, sdif2 (the same as previous)
    vMP.clear();
    obs_0 = 0.0;
    mean  = 0.0;
    count = 0;

    // estimate constant term (ambig.) for required time interval
    vector<shared_ptr<t_gobsgnss>>::iterator it = vobs_prn.begin();

    for( it = vobs_prn.begin(); it != vobs_prn.end(); ++it){

      obs = (*it)->MP(code,gbL1,gbL2);

#ifdef DEBUG       
      cout << prn << " " << go << " DATA: "
           << "  mp:" << obs
           << "  ac:" << gobsattr2str(code.attr())
           << "  co:" << (*it)->obs_C(code)
           << "  l1:" << (*it)->obs_L(gbL1)
           << "  l2:" << (*it)->obs_L(gbL2) 
           << endl;
#endif
      if( double_eq(obs,0.0) ){ continue; } // no data

      // cycle-slip detection, vMP preparation
      difMP = obs_0 - obs;
      obs_0 = obs;

      // NOT NECESSARY to repeat cycle-slip detection for 2nd, 3rd obs-type (vobs_prn prepared)
      // if 'cycle-slip' to be done individually, reset rmsMP & sdif2 !!!!
      if( ! firstCOD ){
        mean += obs;
	vMP.push_back(obs);
        continue;
      }

      double chk_limit = _mp_limit * rmsMP;
      if( vMP.size() > 0 && fabs(difMP) < chk_limit ){   // klouzava RMS
        sdif2 += difMP*difMP;
	rmsMP  = sqrt(sdif2 / ++count);
      }

#ifdef DEBUG
//      if( nepo == DAY_EPOCHS ){
        cout << (*it)->epoch().str_ymdhms() 
	     << " " << fixed << W(4) << prn << W(4) << go
	     << " " << setprecision(0) << W(3)  << vMP.size()
	     << " " << setprecision(0) << W(3)  << count
	     << " " << setprecision(3) << W(16) << obs
	     << " " << setprecision(3) << W(16) << difMP
	     << " " << setprecision(1) << W(6)  << rmsMP;
        if( vMP.size() > 0 && fabs(difMP) < chk_limit )  // klouzava RMS
	                            cout << "\n";
	else if( vMP.size() == 0 )  cout << " -> OBS initiate ! \n";
        else                        cout << " -> CS ! \n";
//      }
#endif	 

      if( vMP.size() == 0 || fabs(difMP) < _mp_limit * rmsMP )
      { 
        mean += obs;
	vMP.push_back(obs);
	if( vMP.size() >= nepo ) break;                      // data already collected

      }else if( vMP.size() < nepo ){                         // check minimum epochs
#ifdef DEBUG
	cout << ep << " " << prn << fixed << setprecision(3) 
	     << " " << go << " vMP: " << W(2) << vMP.size()
	     << " MP slip:" << W(8) << difMP << " MP rms:" << W(8) << rmsMP << endl;
#endif
        mean = 0.0;
        vMP.clear();
        vobs_prn.erase(vobs_prn.begin(), it);                    // reduce initial vector of observation
	it = vobs_prn.begin();
        continue;

      }else break;                                           // data collected
       
    } // loop over observations

     
    if( vMP.size() < nepo ) continue;                        // check minimum # epochs
    if(   it != vobs_prn.end() &&
	++it != vobs_prn.end()) vobs_prn.erase(it,vobs_prn.end()); // erase rest of data

#ifdef DEBUG
  cout << ep << " " << prn
       << " go: " << go 
       << " bL1:" << gobsband2str(gbL1.band())
       << " bL2:" << gobsband2str(gbL2.band())
       << " o1: " << gobs2str( (*vobs_prn.begin())->id_phase(gbL1.band()) )
       << " o2: " << gobs2str( (*vobs_prn.begin())->id_phase(gbL2.band()) )
       << " .....: " << vMP.size()
       << endl; cout.flush();
#endif

    mean /= vMP.size();
    rms = 0.0;
    count = 0;
    for( itV = vMP.begin(); itV != vMP.end(); ++itV )
    {     
      // more robust statistics (mainly for 2nd, 3rd, .. code)
      if( fabs(*itV-mean) < _mp_limit * rmsMP ){
	 count++;
	 rms += (*itV-mean)*(*itV-mean);

#ifdef DEBUG	 
      }else{
      cout << fixed << ep
 	   << " " << prn
	   << " " << go 
	   << " " << setprecision(0) << W(3) << vMP.size()
	   << " " << setprecision(3) << W(8) << *itV
	   << " " << setprecision(3) << W(8) << *itV-mean
	   << " " << setprecision(2) << W(8) << rmsMP
	   << " MP outlier !\n";
#endif
      }
    }

    firstCOD = false;
    t_gtime tlast(LAST_TIME);

    rms = sqrt(rms/(vMP.size()-1)) * MP_UNITS;
     
    if( rms > 999.0 ){      
      // max printable value   
      m_obs[gobs][epo][prn].first  = 999.0;       // DOUBLE (RMS)
      m_obs[gobs][epo][prn].second = vMP.size();  // COUNT
       
    }else{
      m_obs[gobs][epo][prn].first  = rms;         // DOUBLE (RMS)
      m_obs[gobs][epo][prn].second = vMP.size();  // COUNT
      m_obs[gobs][epo][""].first  += rms;         // DOUBLE (Sat mean RMS)
      m_obs[gobs][epo][""].second++;              // COUNT

      m_obs[gobs][tlast][prn].first += rms;       // DOUBLE (Epo mean RMS)
      m_obs[gobs][tlast][prn].second++;           // COUNT
      m_obs[gobs][tlast][""].first  += rms;       // DOUBLE (Epo/Sat mean RMS)
      m_obs[gobs][tlast][""].second++;            // COUNT
    }     
  }
}
   
   

} // namespace
