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

#include "gall/galloqc.h"
#include "gall/gallprec.h"
#include "gset/gsetinp.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_galloqc::t_galloqc()
   : t_gallobs(),
    _beg(FIRST_TIME),
    _end(LAST_TIME),
    _cutoff(0.0),
    _useHealth(DEF_HEALTH)

{
  gtrace("t_galloqc::constructor");
  id_type(  t_gdata::ALLOBS );
  id_group( t_gdata::GRP_OBSERV );  
}


// destructor
// ----------
t_galloqc::~t_galloqc()
{
  gtrace("t_galloqc::destructor");
}


// settings
// ----------
void t_galloqc::gset(t_gsetbase* gset)
{  
  t_gallobs::gset(gset);
  
  _cutoff    = dynamic_cast<t_gsetproc*>(_set)->minimum_elev();
  _useHealth = dynamic_cast<t_gsetqc*>(_set)->useHealth();

  return;
}
   
// Attribute/signal selection (quantitative priority)
// --------------------------------------------------
GOBSATTR t_galloqc::select_attr( const t_map_stt_sys& stt, GSYS gs, GOBSTYPE gt, GOBSBAND gb )
{   
  gtrace("t_galloqc::select_attr");

  GOBSATTR ga = ATTR;

  if( stt.find(gs) == stt.end() ) return ga;

  map<GOBSATTR,int> signal_sel;

  t_map_stt_obs::const_iterator itGOBS;
  t_map_stt_sat::const_iterator itSAT;

  for( itGOBS = stt.at(gs).begin(); itGOBS != stt.at(gs).end(); ++itGOBS )
  {
    GOBS gobs = itGOBS->first;
    string go = gobs2str(gobs);
     
    if( gb != str2gobsband(go) ||
        gt != str2gobstype(go)   ) continue;

    for( itSAT = stt.at(gs).at(gobs).begin(); itSAT != stt.at(gs).at(gobs).end(); ++itSAT )
    {
      signal_sel[str2gobsattr(go)] += itSAT->second.first;
    }     
  }

  // maximize
  int max_count = 0;
  for( map<GOBSATTR,int>::const_iterator it = signal_sel.begin(); it != signal_sel.end(); ++it )
  { 
    if( it->second > max_count ){ 
      ga = it->first;
      max_count = it->second;
//    cout << " minimize 1 : " << max_count << " " << gobsattr2str( ga ) << " " << gobsband2str( gb ) << endl;

    }
  }
   
  return ga;
}


// return first position for satellite
// ----------
t_gtime t_galloqc::beg_obs(const string& site, double smpl)
{
  gtrace("t_gallobs::beg_obs");

  t_gtime tt = t_gallobs::beg_obs(site, smpl);
  
  if( tt < _beg ) return _beg;

  return tt;
}


// return last position for satellite
// ----------
t_gtime t_galloqc::end_obs(const string& site)
{
  t_gtime tt = t_gallobs::end_obs(site);

  if( tt > _end ) return _end;

  return tt;
}


// Observation statistics
// ----------------------
t_galloqc::t_map_stt_sys t_galloqc::stat_obs( string site,
                                              double* smp_est,
                                              double* min_ele,
                                              double  cut_off,
//                                            const t_gtime& beg,
                                              const t_gtime& end,
                                              t_map_sky_sat* sky_sat,
                                              t_map_stt_smp* smp_stt,
                                              t_map_sky_sys* sky_stt
                                            )
{
  gtrace("t_galloqc::stat_obs");

  _gmutex.lock();

  string oo;
  
  t_map_stt_smp smp_stt_loc;

  double sample = 0.0;
  if( smp_est ) sample = (*smp_est);
  
  t_map_stt_sys stat_obs;
  
  if( _mapobj.find(site) == _mapobj.end() ){ _gmutex.unlock(); return stat_obs; }

  t_map_stt_sys::const_iterator itGSYS;
  t_map_stt_sat::const_iterator itPRN;
  t_map_stt_obs::const_iterator itGOBS;
  
  t_map_oref::const_iterator itEpo = _mapobj[site].begin(); // start
  t_map_osat::const_iterator itSat;
  
  vector<GOBS>::iterator itVOBS;
  t_gtime tt     = itEpo->first;
  t_gtime tt_sav = itEpo->first;
   
  for( ; itEpo != _mapobj[site].end(); ++itEpo ) // loop over epochs in interval
  {
    tt = itEpo->first;
    
    if( tt < _beg ) continue;
    if( tt > _end || tt > end ) break;  // here it combines two end requests ???

    double dint = fabs(tt-tt_sav);
    // collect sampling histogram
    if( dint > 0.0 ){
      // carefull counting for > 1Hz
      double tmp = round(dint*1000)/1000;           // round to scaled digits (seconds <1Hz, round to msec)
      if( _scl > 1.0 ) tmp = round(dint*_scl)/_scl; // round to scaled digits (high-rate)

      // if(dint > 0.5)
      //      cerr << tt.str_ymdhms("temp: ")+dbl2str(tt.dsec())  << " " << tmp << " " << dint << " " << _scl << endl;
      if(tt_sav != tt){ smp_stt_loc[tmp] += 1; }
      tt_sav = tt;
    }

    for( itSat = _mapobj[site][tt].begin(); itSat != _mapobj[site][tt].end(); ++itSat ) // loop over satellites
    {
      GSYS   sys = itSat->second->gsys();
      string sat = itSat->second->sat();

      if( !itSat->second->health() && _useHealth >= STT_HEALTH ){
//        cout << sat << " is unhealthy (skip obs_stat) " << _useHealth << " \n";
        continue;
      } // skip unhealthy

#ifdef DEBUG
      cout << " ok sat [" << sat << "] " << tt.str_ymdhms() << " dint:" << dint << " sample:" << sample << " DIFF:" << DIFF_SEC(dint) << endl;
#endif
      
      // handle requested sampling (consider missing epochs with no elevation support)
      if( sample != 0.0 && fabs(dint - sample) > DIFF_SEC(dint) ) // POTENTIAL PROBLEM (only 1Hz now) NEED A CHANGE IF >1Hz
      {
        if( dint < sample ){                                      // skip the epoch
#ifdef DEBUG
          if( sat[0] == 'G' ){
            cout << " skip this epoch due to sample: "
                 << fixed << setprecision(15)
                 << setw(20) << dint << setw(20) << sample << setw(20) << fabs(dint - sample) << endl;
          }
#endif         
          continue;
        }
      }

      // existing GOBS in the epoch
      vector<GOBS> obs = itSat->second->obs();
      for( unsigned int i=0; i < obs.size(); ++i )
      {
        if( itSat->second->getobs(obs[i]) != 0.0 )
        {
          double ele = itSat->second->getele();
          
#ifdef DEBUG         
          if( ele < *min_ele ){ 
            *min_ele = ele; cout << sat << tt.str_ymdhms(" ") << " " << ele << " --> "<< *min_ele << endl;
          }else{            cout << sat << tt.str_ymdhms(" ") << " " << ele << endl; }
#endif

          if( ele >= 0.0 && ele < *min_ele )    // only if status HEALTH (checked already above)
          {
            oo = site+tt.str_ymdhms(" Min.Elev ")+" ["+sat+"]:"+dbl2str(ele)+" deg [<"+trim(dbl2str(*min_ele))+" deg]";
            if( _log && _log->verb() >= 3 ){ if( _log ) _log->comment(3,"galloqc",oo); }

            *min_ele = ele;
          }

          stat_obs[sys][obs[i]][sat].first += 1;
          if( ele >= cut_off ){
            stat_obs[sys][obs[i]][sat].second += 1; // count above cut-off mask ! (instead of expected observations)
          }

          //      if(ele < 0) cerr << sat << " " << ele << " " << cut_off << " " << *ELEVATIONS.begin() << endl;
          // 
          // add histogram of elevations if requested
          if( sky_stt && ele > *ELEVATIONS.begin() ){
            int e = static_cast<int>(*(prev(ELEVATIONS.lower_bound(ele))));
            (*sky_stt)[sys][obs[i]][e] += 1;
          }
        }
      }
    }
  }
  
  // identify sampling interval (auto)
  if( smp_est ){
    double count = 0;
    for( t_map_stt_smp::iterator it = smp_stt_loc.begin(); it != smp_stt_loc.end(); ++it ){

      if( smp_stt ){ (*smp_stt)[it->first] = it->second; } // copy to output parameter
      
      if( it->second > count ){

        *smp_est = it->first;
        count   = it->second;        
      }
    }
  }

  if( _log ) _log->comment(2,"galloqc",oo);

  _gmutex.unlock(); return stat_obs;
}
  

// gaps number
// ---------------------
int t_galloqc::gaps(string site, int gap, map<t_gtime, int>& len, double& min, double& max)
{
  gtrace("t_galloqc::gaps");

  if( _mapobj.find(site) == _mapobj.end() ) return 0;

  _gmutex.lock();

   min = 86400.0;
   max =     0.0;
   
   t_map_oref::const_iterator itEpo;
   t_gtime t;
   t_gtime t_sav = _mapobj[site].begin()->first;

   for( itEpo = _mapobj[site].begin(); itEpo != _mapobj[site].end(); ++itEpo ) // loop over epochs
   {
     if( itEpo == _mapobj[site].begin() ) continue;
     
     t = itEpo->first;     
     if( t < _beg ) continue;
     if( t > _end ) break;
     
     if( abs(t-t_sav) > gap ) len[t_sav] = (int)abs(t-t_sav);
     if( abs(t-t_sav) < min ) min = abs(t-t_sav);
     if( abs(t-t_sav) > max ) max = abs(t-t_sav);
    
     t_sav = t;
   }
   
   _gmutex.unlock(); return 1;
}


// short block of data
// ----------------------------
int t_galloqc::block(string site, int b, int g, map<t_gtime, int>& len)
{
  gtrace("t_galloqc::block");
   
  if( _mapobj.find(site) == _mapobj.end() ) return 0;

  _gmutex.lock();
  
  t_map_oref::const_iterator itEpo;
  t_gtime t;
  t_gtime t_sav = _mapobj[site].begin()->first;
  t_gtime sb = _mapobj[site].begin()->first;
  t_gtime eb;
  
  for( itEpo = _mapobj[site].begin(); itEpo != _mapobj[site].end(); ++itEpo ) // loop over epochs
  {
    t = itEpo->first;
    if( t < _beg ) continue;
    if( t > _end ) break;
    
    if (abs(t-t_sav) > g){
      eb = t_sav;      
      if (abs(eb-sb) < b){
        len[sb] = (int)abs(eb-sb);
        sb = t;
      }
      else sb = t;
    }
    
    t_sav = t;
  }
   
  _gmutex.unlock(); return 1;
}


// number of code/phase bands
// ----------------------------------------
int t_galloqc::nbands(string site, t_map_bsys& nb, double sampl)
{
  gtrace("t_galloqc::nbands");

  if (_mapobj.find(site) == _mapobj.end()) return -1;

  _gmutex.lock();
 
//  t_map_oref::const_iterator itEpoB = _mapobj[site].begin();
//  t_map_oref::const_iterator itEpoE;

  t_gtime t,tt;

  // all epochs are involved
  if( double_eq(sampl, 0.0) )
  {
    for( auto itEpoB = _mapobj[site].begin(); itEpoB != _mapobj[site].end(); ++itEpoB ) // loop over epochs
//    while( itEpoB != _mapobj[site].end() )
    {
      t = itEpoB->first;
      if( t < _beg ) continue;
      if( t > _end ) break;

      for( auto itSat = _mapobj[site][t].begin(); itSat != _mapobj[site][t].end(); ++itSat )// loop over satellites
      { 
        GSYS sys = itSat->second->gsys();
        string prn = itSat->first;
        if( itSat->second->getele() > _cutoff ){
          if( _useHealth >= ALL_HEALTH ){
            if( itSat->second->health() ){ itSat->second->nbands(nb[sys][t][prn]); }
          }else{                           itSat->second->nbands(nb[sys][t][prn]); }
        }else{                             itSat->second->nbands(nb[sys][t][prn]); }
#ifdef __DEBUG__
        cout << " bands ok = "  << prn
             << " " << t.str_ymdhms(" epo: ") 
             << " " << nb.size()
             << " " << nb[sys].size()
             << " " << nb[sys][t].size()
             << " " << nb[sys][t][prn].first
             << " " << nb[sys][t][prn].second
             << " " << fixed << setprecision(6) << t.dsec()
             << endl;
#endif
      }
//      itEpoB++;
    }
    _gmutex.unlock(); return 1; 
  }
  
  // just sampled epochs are involved
//  auto itEpoB = _mapobj[site].begin();
  auto itEpoE = _mapobj[site].end(); itEpoE--;
//  t = itEpoB->first;
  
// // sync
// while( t.sod()%int(sampl) != 0 && itEpoB != itEpoE ){ t = (++itEpoB)->first; }

  vector<t_gsatdata> v_obs;
  for( auto itEpoB = _mapobj[site].begin(); itEpoB != _mapobj[site].end(); ++itEpoB ) // loop over epochs
  { 
    t = itEpoB->first;

    if( t < _beg ) continue;
    if( t > _end ) break;
    
    if( t.sod()%int(sampl) != 0 ) continue; // sync first sample

//  vector<t_gsatdata> v_obs;
//  while( t < itEpoE->first || t == itEpoE->first )
//  {
//    if( t < _beg ) continue;
//    if( t > _end ) break;

    // handle clock drift!!
    _gmutex.unlock();
    v_obs = obs(site, t);
    _gmutex.lock();
      
    if( v_obs.size() > 0 ) 
    {
      tt = v_obs[0].epoch();
      for( auto itSat = _mapobj[site][tt].begin(); itSat != _mapobj[site][tt].end(); ++itSat )  // loop over satellites
      {
        GSYS sys = itSat->second->gsys();
        string prn = itSat->first;
        
        if( itSat->second->getele() > _cutoff ){
          if( _useHealth >= ALL_HEALTH ){
            if( itSat->second->health() ){ itSat->second->nbands(nb[sys][t][prn]); }
          }else{                           itSat->second->nbands(nb[sys][t][prn]); }
        }else{                             itSat->second->nbands(nb[sys][t][prn]); }
      }
    }
    t.add_dsec(sampl);   //  relative drift (tt+sampl) - not good solution, but should work also for absolute drift (t+sampl)
  }  
  
  _gmutex.unlock(); return 1;
}

   
// est elevations
// ===============
//    if( t < _beg ) continue;     // not necessary, better applied from arguments!
//    if( t > _end ) break;        // not necessary, better applied from arguments!
//---------------
int t_galloqc::est_elevations(string site,
                              t_gtriple xyz_rec,
                              t_gallnav* gnav,
                              const t_gtime& beg,                        // should be setup reasonably !
                              const t_gtime& end,                        // should be setup reasonably !
                              const double& mask,
                              t_map_sky_sat& sat_zero,
                              t_map_sky_sat& sat_mask )
{
  gtrace("t_galloqc::est_elevations");
   
  if( !gnav || xyz_rec.zero() || _mapobj.find(site) == _mapobj.end() ) return -1;

  _gmutex.lock();

  t_map_osat::iterator itSAT;
  t_map_oref::iterator itEPO;
  t_map_oobj::iterator itOBJ;

  map<string,t_gtime> epo_data, epo_mask, epo_zero;  // active satellite (first observation)
     
  double r_radi, s_radi, s_topo;
  double zero   = 0.0;
  int    sample = 300; // GLONASS rather sensitive to this settings (larger better)

  t_gtriple top_sat;
  t_gtriple xyz_sat;

  t_gtime   epo  = beg;
  t_gtime   epoB = _mapobj[site].begin()->first;
  t_gtime   epoE = _mapobj[site].rbegin()->first;
//  cout << epoB.str_ymdhms("===epB: ") <<  epoE.str_ymdhms(" epE: ") << " " << _mapobj[site].size() << endl;

  // all satellites
  set<GSYS> sys_obs = sys(site);
  set<string> sats;
  set<string>::const_iterator itPRN;
  t_map_sats gnss_sats = GNSS_SATS();
  for( t_map_sats::const_iterator itGNS = gnss_sats.begin(); itGNS != gnss_sats.end(); ++itGNS ){
    if( sys_obs.find(itGNS->first) != sys_obs.end() ){             // ACTUAL DATA
      for( itPRN = itGNS->second.begin(); itPRN != itGNS->second.end(); ++itPRN ){
        sats.insert(*itPRN);
      }
    }
  }

  // loop over all interval via SAMPLE step
  for(; ! (epo > end+sample); epo.add_secs(sample) )
  {
//    cout << epo.str_ymdhms("-->epo: ") <<   end.str_ymdhms(" end: ") << endl;
    // loop over all satellites
    for(itPRN = sats.begin(); itPRN != sats.end(); ++itPRN )
    {
      string prn = *itPRN; // cout << "PRN: " << prn << epo.str_ymdhms(" ") << endl;
      t_gtime nav_beg = gnav->beg_gnav( prn );
      t_gtime nav_end = gnav->end_gnav( prn );
      int nav_validity = t_gnav::nav_validity(t_gsys::char2gsys(prn[0]));
//      cout << prn << nav_beg.str_ymdhms(" beg: ") << nav_beg.str_ymdhms(" end: ") << " " << nav_validity << endl;

      // EXCEPTION IF NO DATA
      if( _mapobj[site].find(epo) == _mapobj[site].end() ){
//        cout << epo.str_ymdhms("-->epo: ") <<   end.str_ymdhms(" end: ") << " initiate elev "+prn << endl;
        epo_data[prn] = epo;
      }

      // observed satellite at EPO only !
      else if( (itSAT = _mapobj[site][epo].find(prn)) != _mapobj[site][epo].end() )
      {
        // initialize epo_zero; search the 1st observation (backward)
        if( epo_data.find(prn) == epo_data.end() ){
          epo_data[prn] = epo;
          for( itEPO = _mapobj[site].lower_bound(epo - sample);
               (beg < itEPO->first && epo != itEPO->first && itEPO != _mapobj[site].end()); 
               ++itEPO )
          {
//            cout << itEPO->first.str_ymdhms("     search data .. "+prn+": ") << endl;
            if( itEPO->second.find(prn) != itEPO->second.end() ){
              epo_data[prn] = itEPO->first; 
              break; // 1st observation found
            }
          }
//          cout << epo.str_ymdhms("   data initiate PRN: "+prn+": ") << endl;
        }
//        cout << epo.str_ymdhms("-->epo: ") <<   end.str_ymdhms(" end: ") << " data initiated: "+prn << endl;
      }
      else{ 
//        cout << epo.str_ymdhms("-->epo: ") <<   end.str_ymdhms(" end: ") << " NO initiate "+prn << endl;
        epo_data[prn] = epo;
      }
      
      if( epo_data.find(prn) == epo_data.end() ) continue;

//      cout << "EPO: " << epo_data[prn].str_ymdhms() << " size:" << epo_data.size() << endl;

      // ======================================
      //  ONLY ACTIVE SATELLITES HERE AND LATER
      // ======================================
       
      bool search_back = true;
      vector<double> elev;
      vector<t_gtime> tepo, xepo;
      t_gtime tt = epo, tc = epo;
      double dele, ele, tdif;
      bool health = true;
      
      while( search_back )
      {
        if( epo == beg ) break;
        if(  tt  < beg ){ search_back = false; tt = beg; } // if( idx == 0 ) tt = epo-sample; else tt = epo;
        
        if(  tt  < nav_beg - nav_validity ){ // message before NAV
//             cout << prn+tt.str_ymdhms(" skip [NAV-2*3600]: ")+nav_beg.str_ymdhms(" ") << endl; 
             if( epo_zero.find(prn) != epo_zero.end() ){ sat_zero[prn][epo_zero[prn]] = nav_beg - nav_validity;
                                                         epo_zero.erase(prn); }
             if( epo_mask.find(prn) != epo_mask.end() ){ sat_mask[prn][epo_mask[prn]] = nav_beg - nav_validity;;
                                                         epo_mask.erase(prn); }
             break; 
        }
        
        if(  tt  > nav_end + nav_validity ){ // message before NAV
//             cout << prn+tt.str_ymdhms(" skip [NAV+2*3600]: ")+nav_end.str_ymdhms(" ") << endl;
             if( epo_zero.find(prn) != epo_zero.end() ){ sat_zero[prn][epo_zero[prn]] = nav_end + nav_validity;;
                                                         epo_zero.erase(prn); }
             if( epo_mask.find(prn) != epo_mask.end() ){ sat_mask[prn][epo_mask[prn]] = nav_end + nav_validity;;
                                                         epo_mask.erase(prn); }
             break; 
        }
        if(  tt  < epo - nav_validity  ){ // last EPO too far
//             cout << prn+tt.str_ymdhms(" skip [EPO-2*3600]: ")+    epo.str_ymdhms(" ") << endl;
             if( epo_zero.find(prn) != epo_zero.end() ){ sat_zero[prn][epo_zero[prn]] = epo - nav_validity;;
                                                         epo_zero.erase(prn); }
             if( epo_mask.find(prn) != epo_mask.end() ){ sat_mask[prn][epo_mask[prn]] = epo - nav_validity;;
                                                         epo_mask.erase(prn); }
             break;
        }

        if( gnav->nav( prn, tt, xyz_sat.crd_array(), NULL, NULL, false ) < 0 || xyz_sat.zero() ){
//          cout << "INIT cannot find navigation message: "+prn+tt.str_ymdhms(" ") << endl;

          health = gnav->health(prn, tt);
          xepo.push_back(tt);
          auto last_epo = xepo.rbegin();

          //          && _useHealth >= STT_HEALTH
          if( !health && ++last_epo != xepo.rend() ) // fill only for unhealthy (healthy are implicite!)
          {
//            cout << prn << " is unhealthy (set health) " << _useHealth << " \n";
            for( itEPO  = _mapobj[site].lower_bound(tt);
               ( itEPO != _mapobj[site].end() && itEPO->first < *last_epo ); ++itEPO )
            {
              tc  = itEPO->first;
              if( _mapobj[site][tc].find(prn) != _mapobj[site][tc].end() ){
                _mapobj[site][tc][prn]->health(health);
//                cout << tc.str_ymdhms("         hlt -> fill: "+prn+": ") << fixed << setprecision(3) << setw(8) << -9.0 << " " << health << endl; 
              }
            }
            search_back = false;
          }

        }else{
          health  = gnav->health(prn, tt);
          top_sat = xyz_sat - xyz_rec;
          r_radi  = xyz_rec.crd_cvect().norm_Frobenius()/1000.0; // km
          s_radi  = xyz_sat.crd_cvect().norm_Frobenius()/1000.0; // km
          s_topo  = top_sat.crd_cvect().norm_Frobenius()/1000.0; // km
     
          double e = -90 + acos( ( pow(s_topo, 2.0) + pow(r_radi, 2.0) - pow(s_radi, 2.0) )
                                 / ( 2 * r_radi * s_topo ) ) / G_PI*180;
          elev.push_back(e);
          tepo.push_back(tt);

//          cout << tt.str_ymdhms("        estimate ELE: "+prn+": ") << fixed << setprecision(3)
//               << setw(8) << e << setw(15) << top_sat[0] << setw(15) << top_sat[1] << setw(15) << top_sat[2] << endl;
          
          if( elev.size() < 2 ){ 
            tt.add_secs(-sample); 
            continue; 
          }

          size_t prev = elev.size() - 1;
          size_t next = elev.size() - 2;
          tdif = tepo[next] - tepo[prev];
          dele = elev[next] - elev[prev];

          for( itEPO  = _mapobj[site].lower_bound(tt);
             ( itEPO != _mapobj[site].end() && itEPO->first < tepo[next]+1 ); ++itEPO )
          {
            tc  = itEPO->first;
            ele = elev[prev] + dele * (tc-tepo[prev])/(tdif);
            if( _mapobj[site][tc].find(prn) != _mapobj[site][tc].end() ){
              _mapobj[site][tc][prn]->addele(ele);
              _mapobj[site][tc][prn]->health(health);
//              cout << tc.str_ymdhms("         ele -> fill: "+prn+": ") << fixed << setprecision(3) << setw(8) << ele << " " << health << endl;
            }
          }
          search_back = false; // cout << prn << "  search_back:false\n";
          double DELE_LIMIT = 0.1;
          // mask zero cut-off
          if( elev[prev] < zero && elev[next] < zero ){
//            cout << tt.str_ymdhms("         ele -> skip: "+prn+": ")+dbl2str(elev[prev])+dbl2str(elev[next]) << endl; 
            break;
          }else{ 
//            cout << tt.str_ymdhms("         ele -> zero: "+prn+": ")+dbl2str(elev[prev])+dbl2str(elev[next]) << endl;
          }
          
          if( !health && _useHealth >= STT_HEALTH ){
//            cout << prn << " is unhealthy (skip elev/azi) " << health << " " << _useHealth << " \n";
            continue;
         } // skip all unhealthy

          if( epo_zero.find(prn) == epo_zero.end() ) // search horizon (START/STOP)
          { 
            tc = tepo[prev]; // + tdif/2.0;
            if( fabs(dele) > DELE_LIMIT && tepo[prev] > beg ) tc = tepo[prev] + (zero - elev[prev]) * tdif/dele;
            if( tc < tepo[next]  ){
              epo_zero[prn] = (tc > beg) ? tc : beg; // descend
              if( elev[prev] > mask )
              epo_mask[prn] = (tc > beg) ? tc : beg; // descend
//              cout << tc.str_ymdhms(" .. INIT ele -> zero: "+prn+": ")+dbl2str(dele)+dbl2str(elev[prev])+dbl2str(elev[next]) << endl;
            }else{ break; } // stop
          }else if( elev[prev] < zero || elev[next] < zero ){
            tc = tepo[next]; // tepo[prev] + tdif/2.0;
            if( fabs(dele) > DELE_LIMIT ){ tc = tepo[prev] + (zero - elev[prev]) * tdif/dele; }
            if( tc > tepo[prev] ){
              sat_zero[prn][epo_zero[prn]] = (tc < beg) ? beg : tc; // ascend
              epo_zero.erase(prn);
//              cout << tc.str_ymdhms(" .. STOP ele -> zero: "+prn+": ")+dbl2str(dele)+dbl2str(elev[prev])+dbl2str(elev[next]) << endl;
            } // next iteration
          }
          
          // mask user cut-off
          if( elev[prev] < mask && elev[next] < mask ){
                 // cout << tt.str_ymdhms("         ele -> skip: "+prn+": ")+dbl2str(elev[prev])+dbl2str(elev[next]) << endl;
            break;
          }else{ // cout << tt.str_ymdhms("         ele -> mask: "+prn+": ")+dbl2str(elev[prev])+dbl2str(elev[next]) << endl;
          }
          
          if( epo_mask.find(prn) == epo_mask.end() ){  // search mask (START/STOP)
            if( elev[prev] < mask && elev[next] < mask ) break;
            tc = tepo[prev]; // + tdif/2.0;
            if( fabs(dele) > DELE_LIMIT && tepo[prev] > beg ) tc = tepo[prev] + (mask - elev[prev]) * tdif/dele;
            if( tc < tepo[next] ){
              epo_mask[prn] = (tc > beg) ? tc : beg; // initiate
//               cout << tc.str_ymdhms(" .. INIT ele -> mask: "+prn+": ")+dbl2str(dele)+dbl2str(elev[prev])+dbl2str(elev[next]) << endl;
            }else break; // stop
          }else if( elev[prev] < mask || elev[next] < mask ){
            tc = tepo[next]; // tepo[prev] + tdif/2.0;
            if( fabs(dele) > DELE_LIMIT ) tc = tepo[prev] + (mask - elev[prev]) * tdif/dele;
            if( tc > tepo[prev] ){
              sat_mask[prn][epo_mask[prn]] = (tc < beg) ? beg : tc;
              epo_mask.erase(prn);
//               cout << tc.str_ymdhms(" .. STOP ele -> mask: "+prn+": ")+dbl2str(dele)+dbl2str(elev[prev])+dbl2str(elev[next]) << endl;
            } // next iteration
          }
        }       
        tt.add_secs(-sample);
      }
    }
  }

  map<string,t_gtime>::const_iterator it;
  for( it = epo_zero.begin(); it != epo_zero.end(); ++it ) sat_zero[it->first][it->second] = end;
  for( it = epo_mask.begin(); it != epo_mask.end(); ++it ) sat_mask[it->first][it->second] = end;

  _gmutex.unlock(); return 1;
}


// add elevations
//---------------
int t_galloqc::add_elevations(string site, t_gtriple xyz_rec, t_gallnav* gnav )
{
  gtrace("t_galloqc::add_elevations");
   
  if( !gnav || xyz_rec.zero() || _mapobj.find(site) == _mapobj.end() ) return -1;

  _gmutex.lock();

  t_map_osat::iterator itSAT;
  t_map_oref::iterator itEPO;
  t_map_oobj::iterator itOBJ;

  t_gtriple top_sat;
  t_gtriple xyz_sat;
  t_gtime   epo;
  double r_radi, s_radi, s_topo, ele;
  bool health = true;
  
  // loop over all observations for given object
  for(itEPO = _mapobj[site].begin(); itEPO != _mapobj[site].end(); ++itEPO )
  {
    epo = itEPO->first;
    if( epo < _beg ) continue;
    if( epo > _end ) break;

    for(itSAT = _mapobj[site][epo].begin(); itSAT != _mapobj[site][epo].end(); ++itSAT )
    {
      string prn = itSAT->first;

      if( gnav->nav( prn, epo, xyz_sat.crd_array(), NULL, NULL, false ) < 0 || xyz_sat.zero() ){}
      else
      { 
        health = gnav->health(prn, epo);
        top_sat = xyz_sat - xyz_rec;

        r_radi = xyz_rec.crd_cvect().norm_Frobenius()/1000.0; // km
        s_radi = xyz_sat.crd_cvect().norm_Frobenius()/1000.0; // km
        s_topo = top_sat.crd_cvect().norm_Frobenius()/1000.0; // km
     
        ele = -90 + acos( ( pow(s_topo, 2.0) + pow(r_radi, 2.0) - pow(s_radi, 2.0) )
                               / ( 2 * r_radi * s_topo ) ) / G_PI * 180;

        _mapobj[site][epo][prn]->addele(ele);
        _mapobj[site][epo][prn]->health(health);
      }

#ifdef DEBUG
cout << fixed << setprecision(3) << " ELE " << prn << " calc:"
         << " " << epo.str_ymdhms()
         << " " << ele
         << " " << top_sat[0] << " " << top_sat[1] << " " << top_sat[2] 
//       << " " << xyz_sat[0] << " " << xyz_sat[1] << " " << xyz_sat[2] 
//       << " " << xyz_rec[0] << " " << xyz_rec[1] << " " << xyz_rec[2]
              << endl;
#endif
    }
  }

  _gmutex.unlock(); return 1;
}


// add elevations
//---------------
int t_galloqc::apr_elevations(string site, t_gtriple xyz_rec, t_gallnav* gnav )
{
  gtrace("t_galloqc::apr_elevations");
   
  if( !gnav || xyz_rec.zero() || _mapobj.find(site) == _mapobj.end() ) return -1;

  _gmutex.lock();

  t_map_osat::iterator itSAT;
  t_map_oref::iterator itEPO;
  t_map_oobj::iterator itOBJ;

  t_gtriple top_sat,xyz_sat;
  double r_radi, s_radi, s_topo, ele;
  bool health = true; 
  
  t_gtime beg = _mapobj[site].begin()->first;
  t_gtime end = _mapobj[site].rbegin()->first;
  t_gtime epo = beg;
  double  len = fabs(end - beg);
  double  sec = 86400.0;
  double  smp = 64.0;
  t_gtime ref = beg + len/2.0;
  t_gpoly fit_poly;
  vector<double> X, Y;

#ifdef DEBUG
  cout << " len = " << beg.str_ymdhms("") <<  " - " << end.str_ymdhms("")
       <<     " : " << ref.str_ymdhms("") << " -> " << fixed << len << endl;
#endif

  // loop over available satellites
  set<string> sats = gnav->satellites();
  for( set<string>::iterator itPRN = sats.begin(); itPRN != sats.end(); ++itPRN )
  {
    epo = beg;
    string prn = *itPRN;
    fit_poly.reset(); X.clear(); Y.clear();

    // loop over sampled observations for given object
    while( epo < end+1 )
    { 
      if( epo < _beg ){ epo.add_dsec(len/smp); continue; }
      if( epo > _end ){ break; }

      if( gnav->nav( prn, epo, xyz_sat.crd_array(), NULL, NULL, false ) < 0 || xyz_sat.zero() ){}
      else
      { 
        health = gnav->health(prn, epo);
        
        top_sat = xyz_sat - xyz_rec;

        r_radi = xyz_rec.crd_cvect().norm_Frobenius()/1000.0; // km
        s_radi = xyz_sat.crd_cvect().norm_Frobenius()/1000.0; // km
        s_topo = top_sat.crd_cvect().norm_Frobenius()/1000.0; // km
     
        ele = -90 + acos( ( pow(s_topo, 2.0) + pow(r_radi, 2.0) - pow(s_radi, 2.0) )
                                           / ( 2 * r_radi * s_topo ) ) / G_PI * 180;

        X.push_back(epo.diff(ref)/sec);
        Y.push_back(ele);
    
//      cout << fixed << setprecision(2) << epo.str_ymdhms() << " fit_prn:" << prn << setw(8) << ele << endl;
      }
      epo.add_dsec(len/smp);
    }

    if( X.size()/smp *100.0 < 65 ){ if( _log ) _log->comment(1,"galloqc","few elevation samples for estimating polynomials, PRN: "+prn); continue; }

    // calculate polynomials
    int deg = X.size() - 8;
    if( deg > 13 ) deg = 13;
    fit_poly.fitpolynom(  X, Y, deg, sec, ref);

    // fill all elevations - loop over all observations for given object
    for( t_map_oref::iterator itEPO = _mapobj[site].begin(); itEPO != _mapobj[site].end(); ++itEPO )
    {
      epo = itEPO->first;
      double ele = 0.0; fit_poly.evaluate( epo.diff(ref)/sec, 0, ele );
//    double ele = 0.0, de = 0.0; int_poly.interpolate( X, Y, epo.diff(ref)/sec, ele, de);
//      cout << fixed << setprecision(2) << epo.str_ymdhms() << " ele_prn:" << prn << setw(8) << ele << endl;

      // fill specific satellite
      t_map_osat::iterator itOBJ = itEPO->second.find(prn);
      if( itOBJ != itEPO->second.end() ) {
        itEPO->second[prn]->health(health);
        itEPO->second[prn]->addele(ele);
      }
    }
  }

  _gmutex.unlock(); return 1;
}

   
} // namespace
