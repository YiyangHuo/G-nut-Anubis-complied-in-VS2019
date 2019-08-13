
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
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

#include <iomanip>
#include <memory>

#include "gproc/gpreproc.h"
#include "gutils/gcommon.h"
#include "gutils/gtimesync.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gbancroft.h"

using namespace gnut;

namespace gnut {

// Constructor
// ----------
t_gpreproc::t_gpreproc(t_gallobs* obs, t_gsetbase* settings)
: t_gmonit("gpreproc"),
  _log(0)
{
  _obs = obs;
  _set = settings;   
  _msoffset.clear();
  _beg_end = true;

  _sys = dynamic_cast<t_gsetgen*>(_set)->sys();
  _scl = dynamic_cast<t_gsetgen*>(_set)->sampling_scalefc();
  _sat = dynamic_cast<t_gsetgnss*>(_set)->sat();
  _sigCode      = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(GPS);
  _sigPhase     = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(GPS);
  _sigCode_GLO  = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(GLO);
  _sigPhase_GLO = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(GLO);
  _sigCode_GAL  = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(GAL);
  _sigPhase_GAL = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(GAL);
  _sigCode_BDS  = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(BDS);
  _sigPhase_BDS = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(BDS);
  _sigCode_QZS  = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(QZS);
  _sigPhase_QZS = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(QZS);
  _sigCode_IRN  = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(IRN);
  _sigPhase_IRN = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(IRN);   
   
  for (set<string>::iterator it = _sys.begin(); it != _sys.end(); it++){
    _dI[*it] = 0.0;
  }
   
//   if( _set ){ // AVOID gsetproc !
//     _sigCode  = dynamic_cast<t_gsetproc*>(_set)->sigma_C();
//     _sigPhase = dynamic_cast<t_gsetproc*>(_set)->sigma_L();
//   }
}

// Destructor
// ----------
t_gpreproc::~t_gpreproc()
{ 
}


// run preprocessing from time beg until time end (main method)
// --------------------------------------------------------------
int t_gpreproc::preprocess(string site, const t_gtime& beg_r, const t_gtime& end_r, double sampl, bool sync, bool save)
{
  gtrace("t_gpreproc::preprocess");

  int sign = 1;
  t_gtime beg;
  t_gtime end;   
  
  if(!_beg_end){
    beg = end_r;
    end = beg_r;
    sign = -1;      
    if(_log) _log->comment(2, "t_gpreproc", "Preprocessing in end -> begin direction!");
  }else{
    beg = beg_r;
    end = end_r;
    if(_log) _log->comment(2, "t_gpreproc", "Preprocessing in begin -> end direction!");
  }   
   
  this->_site = site;
  //   cerr << "start prep for " << _site << endl;
  if ( _obs->nepochs(_site) <= 1) return -1;   

  bool firstEpo = true;
  vector<shared_ptr<t_gobsgnss>> epoData;
  vector<t_gobsgnss> epoDataPre;

  double CJ = 0.0;   // integer ms clk jump   

  double subint = 0.1;
  if( _scl > 0) subint = 1.0/_scl;
  if(sampl > 1) subint = pow(10,floor(log10(sampl)));
// cerr << "subint 2:" << subint << " _scl: " << _scl << endl;

  bool time_loop = true;   
  t_gtime epoch(beg); 
  while( time_loop )
  {            
    if(_beg_end && (epoch < end || epoch == end)) time_loop = true;
    else if(_beg_end && epoch > end) time_loop = false;
        
    if(!_beg_end && (epoch > end || epoch == end)) time_loop = true;
    else if(!_beg_end && epoch < end) time_loop = false;
      
    //     cout << epoch.str_ymdhms("epoch: ")+dbl2str(epoch.dsec()) << " sampl:" << sampl << endl;

     // synchronization
    if(sync){
      if( ! time_sync(epoch, sampl, _scl, _log) ){
        //        cerr << epoch.str_ymdhms("synchronize: ")+dbl2str(epoch.dsec()) << " add: " << subint/100 << " smpl:" << sampl << endl;
        epoch.add_dsec( subint/100 ); // add_dsec used for synchronization!
        continue;
      }
    }
      
    if( sampl >= 1 ) epoch.reset_dsec(); //  LOW-RATE (>=1Hz), i.e. if not HIGH-RATE !!
    //    cout << "Preprocessing: " << epoch.str_ymdhms() << " " << epoch.dsec() << " sync: " << sync << endl;
    epoData = _obs->obs_pt(_site, epoch);
  
    if( epoData.size() == 0 ){
      if( sampl >= 1 )epoch.add_secs( sign*(int)sampl ); // =<1Hz data
      else            epoch.add_dsec( sign*sampl ); //  >1Hz data
      continue;
    }
                  
    if( firstEpo ){
      for(unsigned int i = 0; i < epoData.size(); i++) epoDataPre.push_back(*epoData[i]);
      firstEpo = false;
    }       

    // test data gap
    shared_ptr<t_gobsgnss> g = *(epoData.begin());
    double diffEpo = g->epoch() - epoDataPre.begin()->epoch();
    bool gap = !firstEpo && diffEpo > sampl + DIFF_SEC(sampl);
    if(gap && _log){
      _log->comment(2, "t_gpreproc", "Warning: Data gap" + dbl2str(diffEpo,0) + "s (sampling =" + dbl2str(sampl,0) + "s, beg = " +
                                     epoDataPre.begin()->epoch().str_ymdhms() + ", end = " + g->epoch().str_ymdhms() + ")");
    }
    if(save){
      if(gap){
        _gapReport(epoData);
      }
      _compare(epoDataPre, epoData);     
    }      
      
    double Ns = epoData.size();   // total number of viewed satellites
    double Na = Ns;                // number of satellite involved into clk jump detection
    double n  = 0;                 // number of detected jumps in particular epoch
    _sumS = 0;
      
    shared_ptr<t_gobsgnss> obs2;
      
    for (vector<shared_ptr<t_gobsgnss>>::iterator it = epoData.begin(); it != epoData.end(); ++it) {
      if (*it == 0) {
//     cout << "*it = 0" << endl;
        Na--; continue;
      } 
    
      obs2 = *it;
      t_gobsgnss* obs1 = 0; 
  
      vector<t_gobsgnss>::iterator itPre;
      for (itPre = epoDataPre.begin(); itPre != epoDataPre.end(); ++itPre){
        if ( (*it)->sat() == (itPre)->sat() ) {
          obs1 = &(*itPre);
        }     
      }

      if (obs1 == 0) {
        //     cout << "obs1 = 0" << endl;
        Na--; continue;
      } // sat not found in previous epoch
      
      t_gsys gsys(obs1->gsys());  
      if (gsys.gsys() != GPS) Na--;   // clk jump is determined only from GPS       

      string satname  =  obs1->sat();
      if ( _sat.find(satname) == _sat.end() )  continue;  // sat excluded in config

      // Cycle slip detection and repair
      if( _slip(obs1, obs2) > -1 ) _transform(obs2, save);

      // Clock jump detection and repair
      int irc_jmp = _jumps(obs1, obs2);
      if ( irc_jmp ==  1 ) n++;
      if ( irc_jmp == -1 ) Na--;  
    }      // end sats loop
    
    epoDataPre.clear();
      
    for(unsigned int i = 0; i < epoData.size(); i++) epoDataPre.push_back(*epoData[i]);

    if (n == Na) {                                         // clk jump effect at all sats in the same way
      double M = (1e3 * _sumS) / (CLIGHT * n);
      double Js = abs( M-floor(M+0.5) );

      if (Js <= 1e-5) {                                   // it is clock jump instead of common slip
        CJ += floor(M+0.5);                               // ms jump
        if (save) _mbreaks[epoch] = (int)CJ;                             // store for logging
        _remove_slip(epoData);     
      }
    }
      
    if( CJ != 0.0 ){
      _repair(epoData, 1e-3*CJ*CLIGHT);                      // corrected carrier phase
      //        cout << "Opravil jsem epoData: " << epoch.str_hms() << " o hodnotu " << CJ << " ms" << endl;
    }

    if( sampl >= 1 )epoch.add_secs( sign*(int)sampl ); // =<1Hz data
    else            epoch.add_dsec( sign*sampl );      //  >1Hz data
  } 
   
//   cerr << "Preprocess for " << site << " finished" << endl;
  return 1;   
}
   
// Compute threshold for cycle slip
// -------------------------------------
double t_gpreproc::_slipThreshold(string LC, t_spt_gobs obs, t_gband& band1, t_gband& band2, double sampl)
{  
  gtrace("t_gpreproc::_slipThreshold");

  double iono = 0.4/60;    // maximum iono jump/sec  ??? set via XML ???

  double f1 = obs->frequency( band1.band() ); // first  band frequency
  double f2 = obs->frequency( band2.band() ); // second band frequency
  double k1 = f1/(f1-f2);
  double k2 = f2/(f1-f2);
  double k3 = f1/(f1+f2);
  double k4 = f2/(f1+f2);

//  double alfa = (f1*f1)/(f1*f1 - f2*f2);
  
  double sigCode = 0;
  double sigPhase = 0;
   
  GSYS gs = obs->gsys();
  switch (gs){
  case GPS: sigCode  = _sigCode;     sigPhase = _sigPhase;
    break;
  case GLO: sigCode = _sigCode_GLO; sigPhase = _sigPhase_GLO;
    break;
  case GAL: sigCode = _sigCode_GAL; sigPhase = _sigPhase_GAL;
    break;
  case BDS: sigCode = _sigCode_BDS; sigPhase = _sigPhase_BDS;
    break;    
  case QZS: sigCode = _sigCode_QZS; sigPhase = _sigPhase_QZS;
    break;
  case IRN: sigCode = _sigCode_IRN; sigPhase = _sigPhase_IRN;
    break;      
  case SBS: sigCode = _sigCode;     sigPhase = _sigPhase;
    break;
  case GNS: sigCode = _sigCode;     sigPhase = _sigPhase;
    break;
  }
  
   
#ifdef DEBUG
   cout << "sigma phase = " << sigPhase << endl;
   cout << "iono*sample = " << iono*sampl << endl;
   cout << "sample = " << sampl << endl;
#endif   

   if (LC.compare("L4") == 0) {
     return 2.5*sqrt(4*sigPhase*sigPhase) + iono*sampl;
   } else if (LC.compare("MW") == 0) {
     double bias = 0.05;
     double mltp = 0.3;
     return 2.5*sqrt( 2*k1*k1*sigPhase*sigPhase +
                      2*k2*k2*sigPhase*sigPhase +
                      2*k3*k3*sigCode *sigCode  +
                      2*k4*k4*sigCode *sigCode  + bias + mltp);
   } else if (LC.compare("WL") == 0) {
     double bias = 0.05;
     return 2.5*sqrt( 2*k1*k1*sigPhase*sigPhase +
                      2*k2*k2*sigPhase*sigPhase + bias);
   } else return 0.0;      
}
  
    
// check coherency between range and phase caused by clock jump
// --------------------------------------------------------------
int t_gpreproc::_jumps(t_gobsgnss* gobs1, t_spt_gobs gobs2)
{              
  gtrace("t_gpreproc::_jumps");
   
  string prn = gobs2->sat();
  GSYS GS = gobs1->gsys();
   
  if( GS != GPS ) return -1;  
   
  set<GOBSBAND> freq_1 = gobs1->band_avail();
  set<GOBSBAND> freq_2 = gobs2->band_avail();   
   
  if(freq_1.size() < 2 || freq_2.size() < 2) {return -1;}
   
  GOBSBAND band;
  set<GOBSBAND>::reverse_iterator it = freq_1.rbegin();
  // int b1_1 = t_gsys::freq2band(gobs1->gsys(), *it);
  band = *it;  // t_gsys::freq2band(gobs1->gsys(), *it);
  t_gband b1_1( band, ATTR);  it++;
  // int b2_1 = t_gsys::freq2band(gobs1->gsys(), *it);
  band = *it;  // t_gsys::freq2band(gobs1->gsys(), *it);
  t_gband b2_1( band, ATTR);
   
  it = freq_2.rbegin();
  // int b1_2 = t_gsys::freq2band(gobs2->gsys(), *it);
  band = *it; // t_gsys::freq2band(gobs2->gsys(), *it);
  t_gband b1_2( band, ATTR);  it++; 
  // int b2_2 = t_gsys::freq2band(gobs2->gsys(), *it);
  band = *it; // t_gsys::freq2band(gobs2->gsys(), *it);
  t_gband b2_2( band, ATTR);
  
  // j-epoch
  double L1j = gobs2->obs_L(b1_2); // [m]
  double L2j = gobs2->obs_L(b2_2); // [m]
  double P1j = gobs2->obs_C(b1_2); // [m]
  double P2j = gobs2->obs_C(b2_2); // [m] 

  // i-epoch
  double L1i = gobs1->obs_L(b1_1); // [m]
  double L2i = gobs1->obs_L(b2_1); // [m]
  double P1i = gobs1->obs_C(b1_1); // [m]
  double P2i = gobs1->obs_C(b2_1); // [m]

  int    nl = 0;
  double dl = 0.0;
  if (!double_eq(L1j,0.0) && !double_eq(L1i,0.0)){ dl += L1j - L1i; nl++; }
  if (!double_eq(L2j,0.0) && !double_eq(L2i,0.0)){ dl += L2j - L2i; nl++; }
  if( nl == 0 ) return -1;
  dl /= nl;

  int    np = 0;
  double dp = 0.0;
  if (!double_eq(P1j,0.0) && !double_eq(P1i,0.0)){ dp += P1j - P1i; np++; }
  if (!double_eq(P2j,0.0) && !double_eq(P2i,0.0)){ dp += P2j - P2i; np++; }
  if( np == 0 ) return -1;
  dp /= np;

  double S = dp - dl;              // jump detection observable
  double sig = 5;                  // sigma of S
  double k = 1e-3 * CLIGHT - 3*sig; // jump threshold

#ifdef DEBUG
   cout << gobs2->epoch().str_hms() << " " << gobs2->sat() << " S: " << S << " dp: " << dp << " dl: " << dl << " k: " << k << endl;
#endif
   
   if (abs(S) >= k) {   // candidate of clock jump
     _sumS += S; 
     return 1;   
   }
   return 0;
}
 
// Set navigation file (rinexn or sp3)
// --------------------------
void t_gpreproc::setNav(t_gallnav* nav)
{
  this->_nav = nav;
}

// Repair phase due to clk jump
// -----------------------------------
void t_gpreproc::_repair(vector<t_spt_gobs> epoData, double dL)
{
  gtrace("t_gpreproc:_repair");

  for (vector<t_spt_gobs>::iterator it = epoData.begin(); it != epoData.end(); it++)
  {
    (*it)->mod_L(dL, X);  // modify all phases by dL [m]
  }
   
}

// Set site name
// -------------------
void t_gpreproc::setSite(string site)
{
  this->_site = site;
}

// Get map with cysle slips
// ---------------------------------
map<t_gtime, map<string, map<GOBS, double> > > t_gpreproc::getSlips()
{
  return _mslips;
}
   
// Get map with cysle slips : gap version
// ---------------------------------
map<t_gtime, map<string, map<GOBS, int> > > t_gpreproc::getSlipsGap()
{
  return _mslipsGap;
}

// Get map with clk jumps
// ---------------------------------
map<t_gtime, int> t_gpreproc::getClkjump()
{
  return _mbreaks;
}

// Transform dN to slips on individual band
// ----------------------------------
int t_gpreproc::_transform(t_spt_gobs gobs,  bool save)
{   
  gtrace("t_gpreproc:_transform");

  int d = _v_lcslp.size();
  if (d <= 1) return -1;
   
  unsigned int maxit = 1;
  for(t_vec_slp::iterator it = _v_lcslp.begin(); it != _v_lcslp.end(); it++){
    if(it->second.size() > maxit) maxit = it->second.size();
  }

  for(unsigned int tmp = 1; tmp <= maxit; tmp++){
    Matrix M(d,d); M = 0;
    ColumnVector S(d); S = 0;
    map<GOBS, double> m_orig;      
   
    bool trans = false;

    for(int i = 1; i < d; i++){
      // find narr-lane slp
      if(_v_lcslp.find(i+10) != _v_lcslp.end()) {
        S(d) = - _v_lcslp[i+10].begin()->val;
        if( !double_eq(S(d), 0.0) ) trans = true;  
        m_orig[_v_lcslp[i+10].begin()->obs1.gobs()] = 0;
        m_orig[_v_lcslp[i+10].begin()->obs2.gobs()] = 0;
      
        if(_v_lcslp[i+10].size() >= 2) {
          _v_lcslp[i+10].erase(_v_lcslp[i+10].begin());
        }
        
        M(d,1) = 2;
        M(d,2) = -1;
      }
     
      // find wide-lane slps
      if(_v_lcslp.find(i) != _v_lcslp.end()){
            
        M(i,d) = -1;
        M(i,d-i) = 1;

        S(i)    = - _v_lcslp[i].begin()->val;
      
        m_orig[_v_lcslp[i].begin()->obs1.gobs()] = 0;
        m_orig[_v_lcslp[i].begin()->obs2.gobs()] = 0;
        if(_v_lcslp[i].size() >= 2) {
          _v_lcslp[i].erase(_v_lcslp[i].begin());
        }
        
        if( !double_eq(S(d-i), 0.0) ) trans = true;
        
      }
    }      

//cout << gobs->sat() << " " << gobs->epoch().str_hms() << endl;      
//cout << "M = " << M << endl;
//cout << "S = " << S << endl;

    ColumnVector O(S); O = 0;
      
    if(trans){
      O = M.i()*S;
      set<GOBSBAND> bands;
      for(map<GOBS, double>::iterator it = m_orig.begin(); it != m_orig.end(); it++){
        t_gobs s;
        s.gobs(it->first);
        //      int f = t_gsys::band2freq(gobs->gsys(), s.band());
        bands.insert(s.band());
      }

      vector<GOBSBAND> vec_bnd = sort_band(gobs->gsys(), bands);
         
      int i = 1;
      for(auto it = vec_bnd.begin(); it != vec_bnd.end(); it++, i++){
        //      GOBSBAND band = t_gsys::freq2band(gobs->gsys(), *it);
        for(map<GOBS, double>::iterator it2 = m_orig.begin(); it2 != m_orig.end(); it2++){
          t_gobs s;
          s.gobs(it2->first);
          
          if(s.band() == *it) {
            it2->second = O(i);
            // save lli to t_gobsgnss
            gobs->addlli( it2->first, 1 );
            gobs->addslip(it2->first, int(it2->second));
            
            if (_log ) {
              ostringstream msg;
              msg << _site << ": cycle slip estimated and stored: " << gobs->epoch().str_hms() << " " << gobs->epoch().dsec() << " " <<
                gobs->sat() << " " << gobs2str( it2->first ) << " " << dbl2str( it2->second, 0 );
              _log->comment(2, "gpreproc", msg.str() );
            }
          }
        }     
      }   
      
      if (save) _save(gobs, m_orig);  
    }
    
  }
  
  return 1;   
}


// Save estimated cycle slips
// ----------------------------
void t_gpreproc::_save(t_spt_gobs gobs, const map<GOBS,double>& slips)
{
  t_gtime epo = gobs->epoch();
  string  prn = gobs->sat();   
  
  for(map<GOBS, double>::const_iterator it = slips.begin(); it != slips.end(); it++){
    GOBS g  = it->first;
    double slp = it->second;
    //      cout << "preproc test: " << gobs->epoch().str_hms() << " " << gobs->sat() << " " << gobs2str(g) << "  " << slp << endl;
    if( !double_eq(slp, 0.0) ) _mslips[epo][prn][g] = slp;      
  }      
}

// remove slips from gobsgnss* and _mslips
// ---------------------------------------
void t_gpreproc::_remove_slip(vector<t_spt_gobs> gobs)
{
  gtrace("t_gpreproc:_remove_slip");

  t_gtime epo = gobs[0]->epoch();

  map<t_gtime, map<string, map<GOBS, double> > >::iterator itEPO;
  itEPO = _mslips.find(epo);  
  if (itEPO != _mslips.end()) _mslips.erase(itEPO);
   
  for(vector<t_spt_gobs>::iterator itG = gobs.begin(); itG != gobs.end(); itG++){
    set<GOBSBAND> freq = (*itG)->band_avail();
    for( auto itF = freq.begin(); itF != freq.end(); itF++){
      //   int band = t_gsys::freq2band( (*itG)->gsys(), *itF);
      //   GOBSBAND band = t_gsys::freq2band( (*itG)->gsys(), *itF);
      t_gband( *itF, ATTR );
      //   GOBS g = (*itG)->pha_id(*itF);
      GOBS g = (*itG)->id_phase(*itF);
      (*itG)->addlli( g, 0 );   
    }      
  }
}

// only common items remain
// -----------------------------------
void t_gpreproc::_common(set<GOBSBAND>& set1, set<GOBSBAND>& set2)
{
  gtrace("t_gpreproc:_common");

  for( auto it1 = set1.begin(); it1 != set1.end(); ){
    auto it2 = set2.find(*it1);
    if(it2 == set2.end()) {
      auto tmp = it1;
      ++it1;
      set1.erase(tmp);
    }else ++it1;
  }
   
  for(auto it1 = set2.begin(); it1 != set2.end(); ){
    auto it2 = set1.find(*it1);
    if(it2 == set1.end()) {
      auto tmp = it1;
      ++it1;
      set2.erase(tmp);
    }else ++it1;
  }   
}

// only common items remain
// -----------------------------------
void t_gpreproc::_common(set<GOBS>& set1, set<GOBS>& set2)
{   
  gtrace("t_gpreproc:_common");
   
  for(set<GOBS>::iterator it1 = set1.begin(); it1 != set1.end(); ){
    set<GOBS>::iterator it2 = set2.find(*it1);
    if(it2 == set2.end()) {
      set<GOBS>::iterator tmp = it1;
      ++it1;
      set1.erase(tmp);
    }else ++it1;
  }
  
  for(set<GOBS>::iterator it1 = set2.begin(); it1 != set2.end(); ){
    set<GOBS>::iterator it2 = set1.find(*it1);
    if(it2 == set1.end()) {
      set<GOBS>::iterator tmp = it1;
      ++it1;
      set2.erase(tmp);
    }else ++it1;
  }  
  
}

// check phase cycl slips
// --------------------------------
int t_gpreproc::_slip(t_gobsgnss* gobs1, t_spt_gobs gobs2)
{  
  gtrace("t_gpreproc:_slip");   
   
  bool slip = false;

  set<GOBSBAND> bands_t1 = gobs1->band_avail();
  set<GOBSBAND> bands_t2 = gobs2->band_avail();  

  this->_common(bands_t1, bands_t2);      

  if(bands_t1.size() != bands_t2.size()) {
    cerr << "ERROR: problem in t_gpreprocc::_common" << endl;
    return -1; 
  }      

  if(bands_t2.size() <= 1){
    if( _log ) _log->comment(2,"t_gpreproc", gobs2->epoch().str_ymdhms( "Not enough bands available: " ) + " " + dbl2str(gobs2->epoch().dsec()) + " " + gobs2->sat() );
    return -1;
  }  

  int nfreq = bands_t1.size();

  t_gobs s1;
  t_gobs s2;
  t_gobs sF;
  t_gobs s_narr;

  // sort according to wavelenght
  vector<GOBSBAND> sorted_t1 = sort_band(gobs1->gsys(), bands_t1);
  vector<GOBSBAND> sorted_t2 = sort_band(gobs2->gsys(), bands_t2);

  if(sorted_t1.size() < 2 || sorted_t2.size() < 2) return -1;
  
  vector<GOBSBAND>::reverse_iterator itFRQ = sorted_t1.rbegin();
  //   int b = t_gsys::freq2band(gobs1->gsys(), *itFRQ);

  set<GOBS> gf2_t1 = gobs1->obs_phase(*itFRQ);
  set<GOBS> gf2_t2 = gobs2->obs_phase(*itFRQ);

  this->_common(gf2_t1, gf2_t2);     // signals for wide-lane

  //second freq: for narrow-lane
  vector<GOBSBAND>::iterator it_narr = sorted_t1.begin();
  if(nfreq >= 2) it_narr++;
  //   int b_narr = t_gsys::freq2band(gobs1->gsys(), *it_narr);
  set<GOBS> gNL_t1 = gobs1->obs_phase(*it_narr);
  set<GOBS> gNL_t2 = gobs2->obs_phase(*it_narr);   
  this->_common(gNL_t1, gNL_t2);     // signals for wide-lane         

  _m_lcslp.clear();
  _v_lcslp.clear();   

  for(int i = 1; i < nfreq; i++){
    ++itFRQ;

//    b = t_gsys::freq2band(gobs1->gsys(), *itFRQ);
    set<GOBS> gf1_t1 = gobs1->obs_phase(*itFRQ);
    set<GOBS> gf1_t2 = gobs2->obs_phase(*itFRQ);

    this->_common(gf1_t1, gf1_t2);     // reference signal - last one

    //shift reference band if no common signals
    if(gf2_t1.size() == 0){
      //       int b = t_gsys::freq2band(gobs1->gsys(), *itFRQ);
      gf2_t1 = gobs1->obs_phase(*itFRQ);
      gf2_t2 = gobs2->obs_phase(*itFRQ);  
      this->_common(gf2_t1, gf2_t2);     // signals for wide-lane
      i--; nfreq--; continue;
    }

    // find freq for extra-wide-lane
    vector<GOBSBAND>::reverse_iterator itFF;
    set<GOBS> gFF_t1;
    set<GOBS> gFF_t2;
    if(i>1){
      vector<GOBSBAND>::reverse_iterator itFF = itFRQ; --itFF;
      //     b = t_gsys::freq2band(gobs1->gsys(), *itFF);
      gFF_t1 = gobs1->obs_phase(*itFF);
      gFF_t2 = gobs2->obs_phase(*itFF);   
      this->_common(gFF_t1, gFF_t2);     // freq for extra-widelane
    }      

    set<GOBS>::iterator itGOBSFF = gFF_t1.begin();
    set<GOBS>::iterator itGOBSf2 = gf2_t1.begin();
    set<GOBS>::iterator itGOBSNL = gNL_t1.begin();            
                  
    bool endFF = true;
    bool endNL = true;

    for(set<GOBS>::iterator itGOBSf1 = gf1_t1.begin(); ; ){
      if(gf1_t1.size() == 0 || gf2_t1.size() == 0) {
        _v_lcslp.clear();
        break;
      }

      s1.gobs(*itGOBSf1);
      s2.gobs(*itGOBSf2);
   
      t_gobs_pair gobs_pair(s1, s2);

      double diff = 0;
      double lam = CLIGHT / gobs1->frequency_lc( s1.band(), 1, s2.band(), -1 );

      if(i == 1){
        //        cout << "MW test " << gobs2->epoch().str_hms() << " " << gobs2->sat() << endl;      
        double lct1 = gobs1->MW(s1, s2);
        double lct2 = gobs2->MW(s1, s2);
        diff = (lct2 - lct1) / lam ;
#ifdef DEBUG        
         cout << gobs1->epoch().str_hms() << " " << gobs2->epoch().str_hms() << " " << gobs1->sat() << " "
              << gobs2str(s1.gobs()) << "  " << gobs2str(s2.gobs()) 
              << fixed << setprecision(3) << " MW_t1: " << lct1 << " MW_t2: " << lct2 << " diff: " << diff << endl;
#endif         
        t_gband band1 = s1.gband();
        t_gband band2 = s2.gband();  
//      double thr = _slipThreshold("MW", gobs2, band1, band2) / lam;
        double wlSlp = 0;
        if(fabs(diff) > 2){//thr) {
          slip = true;
          wlSlp = round(diff);
        }
        gobs_pair.val = wlSlp;
        _v_lcslp[i].push_back(gobs_pair);
        _m_lcslp[i][gobs_pair] = wlSlp;
      }else{
        //       cout << "WL test " << endl;     
      // Compute wide-lane time differenced observations
        double lct1 = gobs1->LWL(s1, s2);
        double lct2 = gobs2->LWL(s1, s2);  
        double dWL = lct2 - lct1;
        
//cout << gobs1->sat() << " " << gobs2str(s1.gobs()) << "  " << gobs2str(s2.gobs()) << " ";
//cout << fixed << setprecision(3) << " WL_t1: " << lct1 << " WL_t2: " << lct2 << " dWL: " << dWL << endl;
    
        t_gband band1 = s1.gband();
        t_gband band2 = s2.gband();     
//      double thrWL = _slipThreshold("WL", gobs2, band1, band2) / lam;
      
      // Compute extra-wide-lane time differenced observatins 
        sF.gobs(*itGOBSFF);
        lct1 = gobs1->LWL(sF, s2);
        lct2 = gobs2->LWL(sF, s2);      
        double dEWL = lct2 - lct1;
        gobs_pair.obs1 = sF;
    
      // find extra-wide-lane cycle slip from previous cascade
        double ewlSlp = 0;
        int x = i-1;
        ewlSlp = _findSlp(x, gobs_pair);
        double elam   = 0;
        elam = CLIGHT / gobs1->frequency_lc(gobs_pair.obs1.band(), 1, gobs_pair.obs2.band(), -1);     
        t_gband bandF = sF.gband();
//      double thrEWL = _slipThreshold("WL", gobs2, bandF, band2) / elam;
        diff = (dEWL - dWL + ewlSlp*elam) / lam;
        //        cout << gobs1->sat() << " " << gobs2str(sF.gobs()) << "  " << gobs2str(s2.gobs()) << fixed << setprecision(3) << " EWL_t1: " << lct1 << " EWL_t2: " << lct2 << " diff: " << diff << endl;
//      double thr = thrEWL + thrWL;
        gobs_pair.obs1 = s1;
        double wlSlp = 0;
        if(fabs(diff) > 2){//thr) {
          slip = true;
          wlSlp = round(diff);
        }
        gobs_pair.val = wlSlp;
        _v_lcslp[i].push_back(gobs_pair);     
        _m_lcslp[i][gobs_pair] = wlSlp;
      }   
      
     // narrow-lane slip
      if(i == nfreq - 1){
        //        cout << "NL test " << endl;
        if(gf1_t1.size() == 0 || gNL_t1.size() == 0) break;
        double lct1 = gobs1->LWL(s1, s2);
        double lct2 = gobs2->LWL(s1, s2);     
        double dWL = lct2 - lct1;
      
      // find wide-lane slip from last cascade
        double wlSlp = 0;
        wlSlp = _findSlp(i, gobs_pair);
      
        s_narr.gobs(*itGOBSNL);
        lct1 = gobs1->LNL(s1, s_narr);
        lct2 = gobs2->LNL(s1, s_narr);
        double dNL = lct2 - lct1;
        gobs_pair.obs1 = s1;
        gobs_pair.obs2 = s_narr;      
        double nlam = CLIGHT / gobs1->frequency_lc(gobs_pair.obs1.band(), 2, gobs_pair.obs2.band(), -1);
        double disf = _disf(gobs1, gobs2, s1, s_narr);
        if (!slip) _iono(gobs1, gobs2, s1, s_narr);
        string prn = gobs1->sat();
      
        diff = (dWL - dNL - disf*_dI[prn] + wlSlp*lam) / nlam;
        //     cout << gobs1->epoch().str_hms() << " " << gobs2->epoch().str_hms() << " " << gobs1->sat() << " " << gobs2str(s1.gobs()) << "  " << gobs2str(s_narr.gobs()) 
        //         << " diff: " << diff << " " << nlam << "  " << dWL << "  " << dNL << " " << wlSlp << "  " << disf*_dI[prn] << " " << disf << " " << _dI[prn] << endl;
        double nlSlp = 0;
        if(fabs(diff) > 2) {   
          nlSlp = round(diff);
        }
        //cout << gobs1->sat() << " " << gobs2str(s1.gobs()) << "  " << gobs2str(s2.gobs()) << fixed << setprecision(3) << " NL_t1: " << lct1 << " NL_t2: " << lct2 << " diff: " << diff << endl;
        gobs_pair.val = nlSlp;
        _v_lcslp[i+10].push_back(gobs_pair);      
        _m_lcslp[i+10][gobs_pair] = nlSlp;
      }
      
      if(i>1) {
        ++itGOBSFF;
        if(itGOBSFF == gFF_t1.end()) endFF = true;
        else endFF = false;
      }
      if(i==nfreq-1) {
        ++itGOBSNL;
        if(itGOBSNL == gNL_t1.end()) endNL = true;
        else endNL = false;
      }   

      ++itGOBSf1; ++itGOBSf2;
      if(itGOBSf1 == gf1_t1.end() && itGOBSf2 == gf2_t1.end() && endFF && endNL) break;
      if(itGOBSf1 == gf1_t1.end()) --itGOBSf1;
      if(itGOBSf2 == gf2_t1.end()) --itGOBSf2;
      if(i>1 && endFF) --itGOBSFF;
     if(i==nfreq-1 && endNL) --itGOBSNL;
    }
   }
  
  return 1;
}

// Costructor for t_gobs_pair
// ----------------------------
t_gobs_pair::t_gobs_pair(t_gobs& gobs1, t_gobs& gobs2)
: obs1(gobs1),
  obs2(gobs2)
{   
}


// operator< for t_gobs_pair
// --------------------------
bool t_gobs_pair::operator<(const t_gobs_pair& t) const
{
  return ( this->obs1.attr() < t.obs1.attr() ||
      
           this->obs2.attr() < t.obs2.attr() );   
}

// ionosphere scale factor: wl - nl
// --------------------------------------
double t_gpreproc::_disf(t_gobsgnss* gobs1, t_spt_gobs gobs2, t_gobs& s1, t_gobs& s2)
{
  double isf_wl = gobs1->isf_lc(s1.band(), 1, s2.band(), -1);
  double isf_nl = gobs1->isf_lc(s1.band(), 2, s2.band(), -1);
  double disf = isf_wl - isf_nl;
  
  return disf;
}

// ionodphere
// -----------------------------------------
void t_gpreproc::_iono(t_gobsgnss* gobs1, t_spt_gobs gobs2, t_gobs& s1, t_gobs& s2)
{
  double k  = pow(gobs1->frequency(s1.band()),2) / pow(gobs1->frequency(s2.band()),2);

  double Lb1_1 = gobs1->obs_L(s1);
  double Lb2_1 = gobs1->obs_L(s2);
  double Lb1_2 = gobs2->obs_L(s1);
  double Lb2_2 = gobs2->obs_L(s2);
  double dL = Lb1_2 - Lb2_2 - Lb1_1 + Lb2_1;   
  
  double dI = -dL/(k-1);
  
  string prn = gobs1->sat();
   
#ifdef DEBUG   
   cout << "iono: " << dL << " " << Lb1_2 << " " << Lb2_2 << " " << Lb1_1 << " " << Lb2_1 << endl;
#endif   
   
   _dI[prn] = dI;
   
//   ddI = _dI[prn] - dI;      
}

// find wl slip
// ----------------------------------
double t_gpreproc::_findSlp(int& i, t_gobs_pair& gpair)
{
  double wlSlp = 0;
  if( _m_lcslp.find(i) != _m_lcslp.end() ){
    if( _m_lcslp[i].find(gpair) != _m_lcslp[i].end() ) {
      wlSlp = _m_lcslp[i][gpair];
    }        
  }      
  
  return wlSlp;
}

// report epo data as slips due to epoch data gap
// ------------------------------------   
void t_gpreproc::_gapReport(vector<shared_ptr<t_gobsgnss>> epoData)
{   
  gtrace("t_gpreproc:_gapReport");

  for(vector<shared_ptr<t_gobsgnss>>::iterator it = epoData.begin(); it != epoData.end(); it++ ){
    t_gtime epo = (*it)->epoch();
    string prn  = (*it)->sat(); 
    vector<GOBS> allobs = (*it)->obs();      
    for(vector<GOBS>::iterator itGOBS = allobs.begin(); itGOBS != allobs.end(); itGOBS++){
      if(gobs_phase(*itGOBS)) _mslipsGap[epo][prn][*itGOBS] = 1;
    }
  }
}
   
// report epo data as slips due to satellite data gap
// ------------------------------------   
void t_gpreproc::_gapReport(shared_ptr<t_gobsgnss> data)
{   
  gtrace("t_gpreproc:_gapReport");

  t_gtime epo = data->epoch();
  string prn  = data->sat();
  vector<GOBS> allobs = data->obs();      
  for(vector<GOBS>::iterator itGOBS = allobs.begin(); itGOBS != allobs.end(); itGOBS++){
    if(gobs_phase(*itGOBS)) _mslipsGap[epo][prn][*itGOBS] = 2;
  }
}   

// compare two epoch data
// --------------------------------
void t_gpreproc::_compare(vector<t_gobsgnss> data1, vector<shared_ptr<t_gobsgnss>> data2)
{
  gtrace("t_gpreproc:_compare");

  for(vector<shared_ptr<t_gobsgnss>>::iterator it2 = data2.begin(); it2 != data2.end(); it2++){
    t_gtime epo = (*it2)->epoch();
    string prn  = (*it2)->sat();  
    bool foundPRN = false;
    for(vector<t_gobsgnss>::iterator it1 = data1.begin(); it1 != data1.end(); it1++){
      if( (*it2)->sat() == it1->sat() ) {
        foundPRN = true;        
        vector<GOBS> vgobs1 = it1->obs();
        vector<GOBS> vgobs2 = (*it2)->obs();
        for(vector<GOBS>::iterator itGOBS2 = vgobs2.begin(); itGOBS2 != vgobs2.end(); itGOBS2++){
          bool foundGOBS = false;
          for(vector<GOBS>::iterator itGOBS1 = vgobs1.begin(); itGOBS1 != vgobs1.end(); itGOBS1++){
            if(*itGOBS2 == *itGOBS1) {foundGOBS = true; break;}
          }
          if(!foundGOBS && gobs_phase(*itGOBS2)) _mslipsGap[epo][prn][*itGOBS2] = 3;
        }
        break;
      }     
    }
    if(!foundPRN) _gapReport(*it2);
  }
}

   
} // namespace
