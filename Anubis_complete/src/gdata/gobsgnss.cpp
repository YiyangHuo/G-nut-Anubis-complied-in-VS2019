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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "gdata/gobsgnss.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"
#include "gdata/gpppdata.h"

using namespace std;

namespace gnut {

// constructor
// -----------
t_gobsgnss::t_gobsgnss()
  : t_gdata(),
    _apr_ele(-1),
    _channel(DEF_CHANNEL),
    _rtcm_end(2),
    _health(true)
{ 
  id_type(t_gdata::OBSGNSS);
  id_group(t_gdata::GRP_OBSERV);
}

t_gobsgnss::t_gobsgnss(const string& site, const string& sat, const t_gtime& t)
  : t_gdata(),
     _staid(site),
     _satid(sat),
     _gsys(t_gsys::char2gsys( _satid[0] )),
     _epoch(t),
     _apr_ele(-1),
     _channel(DEF_CHANNEL),
     _rtcm_end(2),
     _health(true)
 { 
  id_type(t_gdata::OBSGNSS);
  id_group(t_gdata::GRP_OBSERV);
}


// destructor
// -----------
t_gobsgnss::~t_gobsgnss()
{
#ifdef DEBUG
  cout << "OBS-destruct: " << site() << " " << sat() << " " 
       << epoch().str("  %Y-%m-%d %H:%M:%S[%T] ") << fixed << setprecision(3);
   
  vector<GOBS> v_obs = this->getobs();
  vector<GOBS>::iterator itOBS = v_obs.begin();
  for( ; itOBS != v_obs.end(); ++itOBS ) 
      cout << " " << gobs2str( *itOBS ) << ":" << this->getobs( *itOBS );
  cout  << endl; cout.flush();   
#endif
}


// add observation - special map of doubles
// -----------
void t_gobsgnss::addobs(const GOBS& obs, const double& d)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
//   cout << "addobs: " << sat() << " " << gobs2str(obs) << " " << d << endl;
  _gobs[obs] = d;
  _gmutex.unlock();
  return;
}


// add lost-of-lock indicators (LLI) - special map of integers
// -----------
void t_gobsgnss::addlli(const GOBS& obs, const int& i)
{  
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  if( i != 0 ) _glli[obs] = i;
  _gmutex.unlock();
   return;
}
   
// add an estimated cycle slip
// -----------
void t_gobsgnss::addslip(const GOBS& obs, const int& i)
{  
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _gslip[obs] = i;
  _gmutex.unlock();
   return;
}   


// get vector of observations types (GOBS)
// -----------
vector<GOBS> t_gobsgnss::obs()
{  
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  vector<GOBS> tmp;
  map<GOBS,double>::iterator it = _gobs.begin();  
  while( it != _gobs.end() ){
    tmp.push_back( it->first );
    it++;
  }
   
  _gmutex.unlock(); return tmp;
}


// get vector of observations types (GOBS) for a badn
// -----------
set<GOBS> t_gobsgnss::obs_phase(const int& band)
{  
  _gmutex.lock();
   
  set<GOBS> tmp;
  map<GOBS,double>::iterator it = _gobs.begin();  
  while( it != _gobs.end() ){
     GOBS gobs = it->first;
     int b = gobs2band( gobs );
     if(band == b && gobs_phase(gobs)) tmp.insert( gobs );
     it++;
  }
   
  _gmutex.unlock(); return tmp;
}

// add OBS - code [m], phase [whole cycles], dopler[cycles/sec], snr [DBHZ], ...  observations (double)
// -----------
double t_gobsgnss::getobs(const GOBS& obs)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  double tmp = NULL_GOBS;
  if( _gobs.find(obs) != _gobs.end() ) tmp = _gobs[obs];
   
  _gmutex.unlock(); return tmp;
}


// add OBS - code [m], phase [whole cycles], dopler[cycles/sec], snr [DBHZ], ...  observations (double)
// -----------
double t_gobsgnss::getobs(const string& obs)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  double tmp = NULL_GOBS;
  if( _gobs.find(str2gobs(obs)) != _gobs.end() ){ tmp = _gobs[str2gobs(obs)]; }

  _gmutex.unlock(); return tmp;
}


// get LLI - lost of lock indicator (int)
// -----------
int t_gobsgnss::getlli(const GOBS& obs)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  int tmp = 0;
  if( _glli.find(obs) != _glli.end() ) tmp = _glli[obs];
   
  _gmutex.unlock(); return tmp;
}
   
// get an estimated cycle slip
// -----------
int t_gobsgnss::getslip(const GOBS& obs)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  int tmp = 0;
  if( _gslip.find(obs) != _gslip.end() ) tmp = _gslip[obs];
   
  _gmutex.unlock(); return tmp;
}   


// get LLI - lost of lock indicator (int)
// -----------
int t_gobsgnss::getlli(const string& obs)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  int tmp = 0;
  if( _glli.find(str2gobs(obs)) != _glli.end() ) tmp = _glli[str2gobs(obs)];
   
  _gmutex.unlock(); return tmp;
}


// add observation - special map of doubles
// -----------
void t_gobsgnss::addele(double d)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _apr_ele = d;
  _gmutex.unlock();
  return;
}


// add observation - special map of doubles
// -----------
double t_gobsgnss::getele()
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  double d = _apr_ele;
  _gmutex.unlock(); return d;
}


// set chanel number for Glonass satellites
// -----------------
void t_gobsgnss::channel(int ch)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _channel = ch;
  _gmutex.unlock();
   return;
}


// get chanel number for Glonass satellites
// -----------------
int t_gobsgnss::channel() const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  int tmp = _channel;
   
  _gmutex.unlock(); return tmp;
}


// return # of observations
size_t t_gobsgnss::size() const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  size_t tmp = _gobs.size();
   
  _gmutex.unlock(); return tmp;
}


// // get system-specific frequency for band i
// ---------- 
// TEMPORARY !!!!
double t_gobsgnss::frequency( int b )
{
  GOBSBAND gb = GOBSBAND(b);
  return this->frequency( gb );
}

// // get system-specific frequency for band i
// ---------- 
double t_gobsgnss::frequency( GOBSBAND b )
{
//  return t_gsys::frequency( _gsys, b ); // CHYBI CHANNEL FOR GLONASS !!

  switch( _gsys ){
     
   case GPS  : switch( b ){
                 case BAND_1 : return G01_F;
                 case BAND_2 : return G02_F;
                 case BAND_5 : return G05_F;
                default : return 0.0; }
     
   case GLO  : switch( b ){
                 case BAND_1 : return R01_F(_channel); 
                 case BAND_2 : return R02_F(_channel); 
                 case BAND_3 : return R03_F_CDMA;
                 case BAND_5 : return R05_F_CDMA;
                default : return 0.0; }
     
   case GAL  : switch( b ){
                 case BAND_1 : return E01_F;
                 case BAND_5 : return E05_F;
                 case BAND_6 : return E06_F;
                 case BAND_7 : return E07_F;
                 case BAND_8 : return E08_F;
                default : return 0.0; }
     
   case BDS  : switch( b ){
                 case BAND_2 : return C02_F;
                 case BAND_6 : return C06_F;
                 case BAND_7 : return C07_F;
                default : return 0.0; }
     
   case QZS  : switch( b ){
                 case BAND_1 : return J01_F;
                 case BAND_2 : return J02_F;
                 case BAND_5 : return J05_F;
                 case BAND_6 : return J06_F;
                default : return 0.0; }
   
   case IRN  : switch( b ){
                 case BAND_5 : return I05_F;
                default : return 0.0;}
	  
     
   case SBS  : switch( b ){
                 case BAND_1 : return S01_F;
                 case BAND_5 : return S05_F;
                default : return 0.0; }
     
   case GNS  : return 0.0;
  }

  return 0.0;
}


// get system-specific wavelength for band i
// ---------- 
double t_gobsgnss::wavelength(GOBSBAND b)
{
  double frq = this->frequency( b );
  if( frq != 0.0 ) return CLIGHT/frq;
     
  return 0.0;
}


// get wavelength for iono-fee LC band i,j
// ---------- 
double t_gobsgnss::wavelength_L3(GOBSBAND b1, GOBSBAND b2)
{
  double c1, c2; 
  _coef_ionofree(b1, c1, b2, c2);

  double f1 = this->frequency( b1 );
  double f2 = this->frequency( b2 );
   
  if ( !double_eq(f1, 0.0) && !double_eq(f2, 0.0) ){
     double lamb_L3 = c1 * CLIGHT/f1 + c2 * CLIGHT/f2;     
     return lamb_L3;
  }
  return 0.0;   
}


// get GNSS system from satellite IDOC
// ----------
GSYS t_gobsgnss::gsys()
{
  return _gsys;
}

// code observations [m]
// ----------
double t_gobsgnss::obs_C(const t_gobs& go)
{ 
  double tmp = this->_obs_range( go );      
/*  
  //apply DCB correction for code  
  GOBS g;  
  if(go.attr() == ATTR) g = _id_range(go.band());
  else g = go.gobs();  

  string sat = this->sat();
  if (_gsys == GPS) {
    int prn = str2int(sat.substr(1, 2));
    if (g == C1 || g == C1C) {tmp += t_gbiasDCB::dcbcorr(prn, 1);}
    if (g == C2 || g == C2C) {tmp += t_gbiasDCB::dcbcorr(prn, 2);}    
  }
*/    
  return tmp;
}
   
// code observations [m]
// ----------
double t_gobsgnss::obs_C(const t_gband& gb)
{
  double tmp = this->_obs_range( gb );
/*  
  t_gobs go(TYPE, gb.band(), gb.attr() );
  
  //apply DCB correction for code  
  GOBS g;  
  if(go.attr() == ATTR) g = _id_range(go.band());
  else g = go.gobs();  

  string sat = this->sat();
  if (_gsys == GPS) {
    int prn = str2int(sat.substr(1, 2));
    if (g == C1 || g == C1C) {tmp += t_gbiasDCB::dcbcorr(prn, 1);}
    if (g == C2 || g == C2C) {tmp += t_gbiasDCB::dcbcorr(prn, 2);}    
  }
*/  
  return tmp;
}


// code observations [m]
// ----------
double t_gobsgnss::_obs_range(const t_gobs& go)
{
  // Glonass has 255 channel number 
  if ( ! _valid_obs() ) return NULL_GOBS;

  // AUTO SELECTION
  if( go.attr() == ATTR )
  {
    GOBS gobs = _id_range(go.band());
    if( gobs == X ) return  NULL_GOBS;
    else            return _gobs[gobs];
    return NULL_GOBS;
  }

  // FIXED SELECTION
  // here avoid a PROBLEM WITH the legacy TYPE_P (only when NULL_ATTR)
  // -----------------------------
  map<GOBS,double>::const_iterator it = _gobs.end();

  if( go.attr() == ATTR_NULL &&
      ( go.type() == TYPE ||
        go.type() == TYPE_P )
    )
    {   
      t_gobs gobs(TYPE_P, go.band(), go.attr() );
      it = _gobs.find( gobs.gobs() );
    }  
  
  if( it == _gobs.end() && 
      go.type() != TYPE_P
    ){
    t_gobs gobs(TYPE_C, go.band(), go.attr() );
    it = _gobs.find( gobs.gobs() );
  }

  if( it == _gobs.end() ) return NULL_GOBS;

  return it->second;
}

// code observations [m]
// ----------
double t_gobsgnss::_obs_range(const t_gband& gb)
{
  t_gobs gobs(TYPE, gb.band(), gb.attr() );
   
#ifdef DEBUG
   cout << _epoch.str_hms() << " sat: " << sat() << " band: " << gobsband2str(gb.band()) << " attr: " << gobsattr2str(gb.attr()) << endl;
#endif

  return _obs_range( gobs );
}


// phase observations [m]
// ----------
double t_gobsgnss::obs_L(const t_gobs& go)
{
  double tmp = this->_obs_phase( go.gband() );
  return tmp;
}


// phase observations [m]
// ----------
double t_gobsgnss::obs_L(const t_gband& gb)
{
  double tmp = this->_obs_phase( gb );
  return tmp;
}


// phase observations [m]
// ----------
double t_gobsgnss::_obs_phase(const t_gband& gb)
{
  // Glonass has 255 channel number 
  if ( ! _valid_obs() ) return NULL_GOBS; 

  // AUTO SELECTION
  if( gb.attr() == ATTR )
  {
    GOBS gobs = _id_phase(gb.band());
    if( gobs == X ) return  NULL_GOBS;
    else            return _gobs[gobs]*this->wavelength( gb.band() ); // transfer from whole cycles to meters!;
  }

  // FIXED SELECTION
  t_gobs go(TYPE_L, gb.band(), gb.attr() );

  map<GOBS,double>::const_iterator it;
  it = _gobs.find( go.gobs() );
   
  if( it == _gobs.end() ) return NULL_GOBS;
  return it->second*this->wavelength( gb.band() ); // transfer from whole cycles to meters!
}


// modify phase observations [m]
// ----------
int t_gobsgnss::mod_L(const double& dL, const GOBS& gobs, const int i)
{   
  map<GOBS,double>::iterator it;

  if( gobs == X )
  {	
    for(it = _gobs.begin(); it != _gobs.end(); ++it )
    {
      // phase only !!!!  
      string gobs_str = gobs2str(it->first);
      if (gobs_str.compare(0,1,"L") != 0) continue;
        
      // TEMPORARY CONVERSION !
      GOBSBAND gb1 = GOBSBAND( gobs2band(it->first) );

      if (i == 1) it->second += dL/this->wavelength(gb1);  // transfer back to cycles
      else if (i == 0) it->second += dL; 
      else cout << "Observations not modified!" << endl;
//      cout << "! POZOR ! modifikuji: " << gobs2str(it->first) << " o " << dL << " band " << band << " " << epoch().str_hms() << endl;
    }

  // modify requested only
  }else{
    it = _gobs.find(gobs);
    if( it == _gobs.end() ) return -1;

    // TEMPORARY CONVERSION !
    GOBSBAND gb1 = GOBSBAND( gobs2band(gobs) );
     
    if (i == 1) it->second += dL/this->wavelength(gb1);  // transfer back to cycles
    else if (i == 0) it->second += dL; 
    else cout << "Observations not modified!" << endl;     
//    cout << "! POZOR ! modifikuji: " << gobs2str(gobs) << " o " << dL << " band " << band << endl;

  }
  return 0; // not found
}


// return value of carrier-phase linear commbination frequency
// for 2 or 3 bands with given coefficients
// ----------
double t_gobsgnss::frequency_lc(const int& band1, const double& coef1,
				const int& band2, const double& coef2,
				const int& band3, const double& coef3 )
{  
  double f1 = 0.0;
  double f2 = 0.0;
  double f3 = 0.0;

  // TEMPORARY CONVERSION !
  GOBSBAND gb1 = GOBSBAND( band1 );
  GOBSBAND gb2 = GOBSBAND( band2 );
  GOBSBAND gb3 = GOBSBAND( band3 );

  if( coef1 != 0.0 ){ f1 = this->frequency( gb1 ); if( double_eq(f1,NULL_GOBS) ) return NULL_GOBS; }
  if( coef2 != 0.0 ){ f2 = this->frequency( gb2 ); if( double_eq(f2,NULL_GOBS) ) return NULL_GOBS; }
  if( coef3 != 0.0 ){ f3 = this->frequency( gb3 ); if( double_eq(f3,NULL_GOBS) ) return NULL_GOBS; }

  double freq = ( coef1*f1 + coef2*f2 + coef3*f3 );

#ifdef DEBUG
    cout << _satid << " " << _staid
         << " frequency_lc BAND(" << band1 << "," << band2 << "," << band3 << ") "
         <<         "  COEF(" << coef1 << "," << coef2 << "," << coef3 << ") "
         << fixed << setprecision(3) << setw(16) << f1 << setw(16) << f2 << setw(16) << f3 << "  freq = " << setw(16) << freq << endl;
#endif
  return freq;
}


// return value of carrier-phase linear combination ionosphere scale factor
// for 2 or 3 bands with given coefficients
// ----------
double t_gobsgnss::isf_lc(const int& band1, const double& coef1,
		  	  const int& band2, const double& coef2 )
{  

   double f1 = 1.0;
   double f2 = 1.0;
   
   // TEMPORARY CONVERSION !
   GOBSBAND gb1 = GOBSBAND( band1 );
   GOBSBAND gb2 = GOBSBAND( band2 );

   f1 = this->frequency( gb1 ); if( double_eq(f1,NULL_GOBS) ) return NULL_GOBS;
   f2 = this->frequency( gb2 ); if( double_eq(f2,NULL_GOBS) ) return NULL_GOBS;
   
   double denom = coef1*f1 + coef2*f2;
   double num   = f1*f1*(coef1/f1 + coef2/f2);
   double isf   = num/denom;
   
#ifdef DEBUG
    cout << _satid << " " << _staid
         << " isf_lc BAND(" << band1 << "," << band2 << ") "
         <<         "  COEF(" << coef1 << "," << coef2 << ") "
         << fixed << setprecision(3) << setw(16) << f1 << setw(16) << f2 << setw(16) << "  isf = " << setw(16) << isf << endl;
#endif
   return isf;
}

// return value of carrier-phase linear combination noise factor
// for 2 or 3 bands with given coefficients
// ----------
double t_gobsgnss::pnf_lc(const int& band1, const double& coef1,
		  	  const int& band2, const double& coef2,
			  const int& band3, const double& coef3 )
{  
  double f1 = 1.0;
  double f2 = 1.0;
  double f3 = 1.0;

  // TEMPORARY CONVERSION !
  GOBSBAND gb1 = GOBSBAND( band1 );
  GOBSBAND gb2 = GOBSBAND( band2 );
  GOBSBAND gb3 = GOBSBAND( band3 );

  if( coef1 != 0.0 ){ f1 = this->frequency( gb1 ); if( double_eq(f1,NULL_GOBS) ) return NULL_GOBS; }
  if( coef2 != 0.0 ){ f2 = this->frequency( gb2 ); if( double_eq(f2,NULL_GOBS) ) return NULL_GOBS; }
  if( coef3 != 0.0 ){ f3 = this->frequency( gb3 ); if( double_eq(f3,NULL_GOBS) ) return NULL_GOBS; }

  double denom = pow(coef1*f1 + coef2*f2 + coef3*f3,2);
  double num   = f1*f1*coef1*coef1 + f2*f2*coef2*coef2 + f3*f3*coef3*coef3;
  double pnf   = num/denom;

#ifdef DEBUG
    cout << _satid << " " << _staid
         << " pnf_lc BAND(" << band1 << "," << band2 << "," << band3 << ") "
         <<         "  COEF(" << coef1 << "," << coef2 << "," << coef3 << ") "
         << fixed << setprecision(3) << setw(16) << f1 << setw(16) << f2 << setw(16) << f3 << "  pnf = " << setw(16) << sqrt(pnf) << endl;
#endif
  return sqrt(pnf);
}


// return value of general code linear combination
// for 2 or 3 bands with given coefficients
// ----------
double t_gobsgnss::_lcf_range(const t_gobs* gobs1, const int& coef1,
	 	              const t_gobs* gobs2, const int& coef2,
		              const t_gobs* gobs3, const int& coef3 )
{
  double C1 = NULL_GOBS, f1 = 0.0;
  double C2 = NULL_GOBS, f2 = 0.0;
  double C3 = NULL_GOBS, f3 = 0.0;
   
  if( gobs1 && coef1 != 0.0 ){ f1 = frequency(gobs1->band()); C1 = obs_C(*gobs1); if( double_eq(C1,NULL_GOBS) ) return NULL_GOBS; }
  if( gobs2 && coef2 != 0.0 ){ f2 = frequency(gobs2->band()); C2 = obs_C(*gobs2); if( double_eq(C2,NULL_GOBS) ) return NULL_GOBS; }
  if( gobs3 && coef3 != 0.0 ){ f3 = frequency(gobs3->band()); C3 = obs_C(*gobs3); if( double_eq(C3,NULL_GOBS) ) return NULL_GOBS; }

  double lc = NULL_GOBS;

  if (fabs(C1 - C2) > 100)
  {
	  if (_log)
		  _log->comment(2, "t_gobsgnss", "Inconsistencies of the code observations " + _epoch.str_ymdhms());
	  return lc;
  }
  
  string sat = this->sat();
  
  if( gobs3 ) lc = ( coef1*f1*C1 + coef2*f2*C2 + coef3*f3*C3 )
                 / ( coef1*f1    + coef2*f2    + coef3*f3    );
  else        lc = ( coef1*f1*C1 + coef2*f2*C2 )
                 / ( coef1*f1    + coef2*f2    );

#ifdef DEBUG
    cout << fixed << setprecision(3)
         << " _lcf_range C1(" << gobs1->band() << ":" << coef1 << "):" << setw(16) << C1
         <<           "  C2(" << gobs2->band() << ":" << coef2 << "):" << setw(16) << C2;
    if(gobs3) cout << "  C3(" << gobs3->band() << ":" << coef3 << "):" << setw(16) << C3;
    cout << " lc = " << lc << endl;
#endif
   
  return lc;
}


// return value of general phase linear combination
// for 2 or 3 bands with given coefficients
// ----------
double t_gobsgnss::_lcf_phase(const t_gobs* gobs1, const int& coef1,
	 	              const t_gobs* gobs2, const int& coef2,
		              const t_gobs* gobs3, const int& coef3 )
{
  double L1 = NULL_GOBS, f1 = 0; 
  double L2 = NULL_GOBS, f2 = 0;
  double L3 = NULL_GOBS, f3 = 0;
   
  if( gobs1 && coef1 != 0.0 ){ f1 = frequency(gobs1->band()); L1 = obs_L(*gobs1); if( double_eq(L1,NULL_GOBS) ) return NULL_GOBS; }
  if( gobs2 && coef2 != 0.0 ){ f2 = frequency(gobs2->band()); L2 = obs_L(*gobs2); if( double_eq(L2,NULL_GOBS) ) return NULL_GOBS; }
  if( gobs3 && coef3 != 0.0 ){ f3 = frequency(gobs3->band()); L3 = obs_L(*gobs3); if( double_eq(L3,NULL_GOBS) ) return NULL_GOBS; }
   
  double lc = NULL_GOBS;

  if( gobs3 ) lc = ( coef1*f1*L1 + coef2*f2*L2 + coef3*f3*L3 )
                 / ( coef1*f1    + coef2*f2    + coef3*f3    );
  else        lc = ( coef1*f1*L1 + coef2*f2*L2 )
                 / ( coef1*f1    + coef2*f2    );

#ifdef DEBUG
    cout << fixed << setprecision(3)
         << " _lcf_phase L1(" << gobs1->band() << ":" << coef1 << "):" << setw(16) << L1
         <<           "  L2(" << gobs2->band() << ":" << coef2 << "):" << setw(16) << L2;
    if(gobs3) cout << "  L3(" << gobs3->band() << ":" << coef3 << "):" << setw(16) << L3;
    cout << " lc = " << lc << endl;
#endif
   
  return lc;
}


// return coefficients (c1,c2) of the code multipath linear combination
// ----------
int t_gobsgnss::_coef_multpath(GOBSBAND bC, double& cC,
		 	       GOBSBAND b1, double& c1,
                               GOBSBAND b2, double& c2)
{
  if( b1 == b2 ){ c1 = c2 = 0.0; return -1; } // phase should be different!

  double pfacC = pow(this->frequency( bC ), 2);
  double pfac1 = pow(this->frequency( b1 ), 2);
  double pfac2 = pow(this->frequency( b2 ), 2);

  double denominator = (pfac1-pfac2)*pfacC;

  cC = 1.0;
  c1 = - (pfac1*(pfacC+pfac2))/denominator;
  c2 = + (pfac2*(pfacC+pfac1))/denominator;

  return 0;
}


// Ionosphere-linear code combination [m]
// ----------
double t_gobsgnss::P3(GOBSBAND b1, GOBSBAND b2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  double coef1, coef2;
   
  t_gobs g1( TYPE, b1, ATTR );
  t_gobs g2( TYPE, b2, ATTR );

  _coef_ionofree( g1.band(), coef1, g2.band(), coef2 );
  
  double lc = _lc_range( &g1, coef1, &g2, coef2 ); // still dual-frequency only

//  cout << "t_gobsgnss::P3() " << sat() << " " << fixed << setprecision(3) << lc << endl;

  _gmutex.unlock(); return lc;
}


// Ionosphere-linear code combination [m]
// ----------
double t_gobsgnss::P3(const t_gobs& g1, const t_gobs& g2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  double coef1, coef2;
   
  _coef_ionofree( g1.band(), coef1, g2.band(), coef2 );

  double lc = _lc_range( &g1, coef1, &g2, coef2 ); // still dual-frequency only

//  cout << "t_gobsgnss::P3() " << sat() << " " << fixed << setprecision(3) << lc << endl;

  _gmutex.unlock(); return lc;
}


// get Geometry-free phase combination [m]
// ----------
double t_gobsgnss::P4(GOBSBAND b1, GOBSBAND b2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  double coef1, coef2;
  
  t_gobs g1( TYPE, b1, ATTR );
  t_gobs g2( TYPE, b2, ATTR );

  _coef_geomfree( g1.band(), coef1, g2.band(), coef2 );

  double lc = _lc_range( &g1, coef1, &g2, coef2 ); // still dual-frequency only
  _gmutex.unlock(); return lc;
}


// get Geometry-free phase combination [m]
// ----------
double t_gobsgnss::P4(const t_gobs& g1, const t_gobs& g2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  double coef1, coef2;
     
  _coef_geomfree( g1.band(), coef1, g2.band(), coef2 );

  double lc = _lc_range( &g1, coef1, &g2, coef2 ); // still dual-frequency only
  _gmutex.unlock(); return lc;
}


// Ionosphere-free linear phase combination [m] !
// ----------
double t_gobsgnss::L3(GOBSBAND b1, GOBSBAND b2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

   if( _gsys == GLO && double_eq(_channel, DEF_CHANNEL)) {
//      cout << epoch().str_hms() << " " << sat() << " " << _channel << endl;      
      _gmutex.unlock(); 
      return 0.0;
   }   
   
  double coef1, coef2;
   
  t_gobs g1( TYPE_L, b1, ATTR );
  t_gobs g2( TYPE_L, b2, ATTR );
   
  _coef_ionofree( b1, coef1, b2, coef2 );

  double lc = _lc_phase( &g1, coef1, &g2, coef2 ); // still dual-frequency only

//  cout << "t_gobsgnss::L3() " << sat() << " " << fixed << setprecision(3) << lc << endl;

  _gmutex.unlock(); return lc;
}


// Ionosphere-free linear phase combination [m] !
// ----------
double t_gobsgnss::L3(const t_gobs& g1, const t_gobs& g2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

   if( _gsys == GLO && double_eq(_channel, DEF_CHANNEL)) {
//      cout << epoch().str_hms() << " " << sat() << " " << _channel << endl;      
      _gmutex.unlock(); 
      return 0.0;
   }   
   
  double coef1, coef2;

  _coef_ionofree( g1.band(), coef1, g2.band(), coef2 );

  double lc = _lc_phase( &g1, coef1, &g2, coef2 ); // still dual-frequency only

//  cout << "t_gobsgnss::L3() " << sat() << " " << fixed << setprecision(3) << lc << endl;

  _gmutex.unlock(); return lc;
}


// get Geometry-free combination [m] !
// ----------
double t_gobsgnss::L4(GOBSBAND b1, GOBSBAND b2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
   if( _gsys == GLO && double_eq(_channel, DEF_CHANNEL)) {
      _gmutex.unlock(); 
      return 0.0;
   }

  t_gobs g1( TYPE_L, b1, ATTR );
  t_gobs g2( TYPE_L, b2, ATTR );

  double coef1, coef2;

  _coef_geomfree( b1, coef1, b2, coef2 );

  double lc = _lc_phase( &g1, coef1, &g2, coef2 ); // still dual-frequency only

#ifdef DEBUG   
  cout << " L4 = " << fixed << setprecision(3)
       << setw(16)      <<  lc
       << setw(16)      << _lc_phase(&g1, 1, &g2, 0)  // first  band
       << setw(16)      << _lc_phase(&g1, 0, &g2, 1)  // second band
       << endl;
#endif   

  _gmutex.unlock(); return lc;
}


// get Geometry-free combination [m] !
// ----------
double t_gobsgnss::L4(const t_gobs& g1, const t_gobs& g2) // const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
   if( _gsys == GLO && double_eq(_channel, DEF_CHANNEL)) {
      _gmutex.unlock(); 
      return 0.0;
   }
   
  double coef1, coef2;

  _coef_geomfree( g1.band(), coef1, g2.band(), coef2 );

//double lc =  phase_lc( &g1, coef1, &g2, coef2 ); // still dual-frequency only
  double lc = _lc_phase( &g1, coef1, &g2, coef2 ); // still dual-frequency only

#ifdef DEBUG   
  cout << " L4 = " << fixed << setprecision(3)
       << setw(16)      <<  lc
       << setw(16)      << _lc_phase(&g1, 1, &g2, 0)  // first  band
       << setw(16)      << _lc_phase(&g1, 0, &g2, 1)  // second band
       << endl;
#endif   

  _gmutex.unlock(); return lc;
}


// Melbourne-Wuebenna combination phase & code [m] !
// ----------
double t_gobsgnss::MW(const t_gobs& g1, const t_gobs& g2, bool phacod_consistent )
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( gsys() == GLO && double_eq(_channel, DEF_CHANNEL) ){
    _gmutex.unlock(); 
    return 0.0;
  }

  double wideL = NULL_GOBS;
  double narrP = NULL_GOBS;

  t_gobs L1( TYPE_L, g1.band(), g1.attr() ); // should always be as input 1st argument
  t_gobs L2( TYPE_L, g2.band(), g2.attr() ); // should always be as input 2nd argument

  t_gobs C1( TYPE, g1.band(), ATTR ); // TYPE universal to handle P+C code
  t_gobs C2( TYPE, g2.band(), ATTR ); // TYPE universal to handle P+C code

  if( phacod_consistent ){  // overwrite attribute
     C1.gattr( g1.gattr() );
     C2.gattr( g2.gattr() );
  }

  wideL = _lcf_phase( &L1, 1, &L2, -1 );
  narrP = _lcf_range( &C1, 1, &C2,  1 );

  if( double_eq(narrP, NULL_GOBS) ){ _gmutex.unlock(); return NULL_GOBS; }
  if( double_eq(wideL, NULL_GOBS) ){ _gmutex.unlock(); return NULL_GOBS; }

  double lc = narrP - wideL;  

#ifdef DEBUG
  cout << _epoch.str_ymdhms() << " " << _satid << " " << gobs2str(L1.gobs()) << " " << gobs2str(C1.gobs())<< " " << gobs2str(L2.gobs())<< " " << gobs2str(C2.gobs())  << endl;
  cout << "wideL = " <<  wideL  << endl;
  cout << "narrP = " <<  narrP  << endl;
  cout << "MW = " <<  lc        << endl;
#endif

  _gmutex.unlock(); return lc;
}


// get wide-lane combination for phase [m]!
// ---------------------------------------
double t_gobsgnss::LNL(const t_gobs& g1, const t_gobs& g2)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
   if( _gsys == GLO && double_eq(_channel, DEF_CHANNEL)) {
      _gmutex.unlock(); 
      return 0.0;
   }

  double lc = _lcf_phase( &g1, 2, &g2, -1 ); // still dual-frequency only

#ifdef DEBUG   
  cout << " LNL = " << fixed << setprecision(3)
       << setw(16)  << lc
       << endl;
#endif   

  _gmutex.unlock(); return lc;
}


// get wide-lane combination for phase [m]!
// ---------------------------------------
double t_gobsgnss::LWL(const t_gobs& g1, const t_gobs& g2)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  if( _gsys == GLO && double_eq(_channel, DEF_CHANNEL)) {
     _gmutex.unlock(); 
     return 0.0;
  }
   
  double lc = _lcf_phase( &g1, 1, &g2, -1 ); // still dual-frequency only

#ifdef DEBUG   
  cout << " LWL = " << fixed << setprecision(3)
       << setw(16)  << lc
       << endl;
#endif   

  _gmutex.unlock(); return lc;
}


//get mulitpath LC
//---------------------
double t_gobsgnss::MP(const t_gobs& code, const t_gobs& L1, const t_gobs& L2)
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( !code.is_code() ) return NULL_GOBS;
  if( !L1.is_phase()  ) return NULL_GOBS;
  if( !L2.is_phase()  ) return NULL_GOBS;

  double coef1, coef2, coefC;
  double C = _obs_range( code );

  if( double_eq(C, NULL_GOBS) ){ _gmutex.unlock(); return NULL_GOBS; }

  _coef_multpath( code.band(), coefC, L1.band(), coef1, L2.band(), coef2 );

  double L = _lc_phase( &L1, coef1, &L2, coef2 );
  if( double_eq(L, NULL_GOBS) ){ _gmutex.unlock(); return NULL_GOBS; }

  string sat = this->sat();

  double lc = coefC*C + L;

#ifdef DEBUG
  cout << "MP code " << gobs2str(code.gobs())
       << fixed << setprecision(3)
       <<    ": coef1 = " << setw(7) << coef1
       << "     coef2 = " << setw(7) << coef2
       << "         C = " << setw(7) << C
       << "         L = " << setw(7) << L
       << "        lc = " << setw(7) << lc
       << endl;
#endif

  _gmutex.unlock(); return lc;
}


// valid
// ----------
bool t_gobsgnss::valid()
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  bool tmp = this->_valid();
   
  _gmutex.unlock(); return tmp;
   
}


// clean data
// ----------
void t_gobsgnss::clear()
{ 
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  this->_clear();
   
  _gmutex.unlock();
   return;
}


// clean internal function
// ----------
void t_gobsgnss::_clear()
{
  _staid.clear();
  _epoch = FIRST_TIME;
}


// clean internal function
// ----------
bool t_gobsgnss::_valid() const
{   
  if( _staid.empty() || 
      _staid == ""   || 
      _epoch == FIRST_TIME ) return false;

  return true; 
}

// Band of code/phase
// ------------------------------
void t_gobsgnss::nbands(pair<int,int>& nb)
{
  for (int b=1; b<=8; b++){
    if( _cod_id(b) != X ) nb.first++;
    if( _pha_id(b) != X ) nb.second++;      
  }
}

// Public version of _freq_avail()
// --------------------------------
set<GOBSBAND> t_gobsgnss::band_avail(bool phase)
{
  _gmutex.lock();

  set<GOBSBAND> s_f;
     
  if(phase) s_f = _band_avail();
  else      s_f = _band_avail_code();
  
  _gmutex.unlock();
  return s_f;
}

// Public version of _freq_avail()
// --------------------------------
set<GFRQ> t_gobsgnss::freq_avail()
{
  _gmutex.lock();
  set<GFRQ> s_f = _freq_avail();   
  
  _gmutex.unlock();
  return s_f;
}

// find out wether GFRQ is included   
bool t_gobsgnss::contain_freq(FREQ_SEQ freq)
{
  _gmutex.lock();
  
  bool res = false;
  
  set<GFRQ> frqs = _freq_avail();
  for(auto it = frqs.begin(); it != frqs.end(); it++){
    FREQ_SEQ fs = t_gsys::gfrq2freq(_gsys, *it);
    if(fs == freq) res = true;
  }
  
  _gmutex.unlock();
  return res;
}
   
// Get available frequences
// --------------------------------
set<GFRQ> t_gobsgnss::_freq_avail()
{
  set<GFRQ> s_f;
  for(map<GOBS, double>::iterator it = _gobs.begin(); it != _gobs.end(); it++){
    t_gobs gobs; gobs.gobs( it->first );
    if( double_eq(it->second, 0.0) || !gobs.is_phase() ) continue;
    else {s_f.insert( t_gsys::band2gfrq( _gsys, gobs.band()) );
    }
  }
  return s_f;
}

// Get available bands for phase
// --------------------------------
set<GOBSBAND> t_gobsgnss::_band_avail()
{
   set<GOBSBAND> s_f;
   for(map<GOBS, double>::iterator it = _gobs.begin(); it != _gobs.end(); it++){
      t_gobs gobs; gobs.gobs( it->first );

      if( double_eq(it->second, 0.0) || !gobs.is_phase() ) continue;
      else s_f.insert( gobs.band() );
   }
   return s_f;
}
   
// Get available bands for phase
// --------------------------------
set<GOBSBAND> t_gobsgnss::_band_avail_code()
{
   set<GOBSBAND> s_f;
   for(map<GOBS, double>::iterator it = _gobs.begin(); it != _gobs.end(); it++){
      t_gobs gobs; gobs.gobs( it->first );
      
      if( double_eq(it->second, 0.0) || gobs.is_phase() ) continue;
      else s_f.insert( gobs.band() );
   }
   return s_f;
}

// Validity test for GLONASS data
// ---------------
bool t_gobsgnss::_valid_obs() const
{
   if( _gsys == GLO && _channel == 255)
        return false;
   else return true;
}


// return value of general carrier-phase linear combination
// for 2 or 3 bands with given coefficients
// ----------
double t_gobsgnss::_lc_range(const t_gobs* g1, const double& coef1,
	                     const t_gobs* g2, const double& coef2,
	                     const t_gobs* g3, const double& coef3 )
{
  double C1 = NULL_GOBS;
  double C2 = NULL_GOBS;
  double C3 = NULL_GOBS;

  if( g1 && coef1 != 0.0 ){ C1 = obs_C( *g1 ); if( double_eq(C1,NULL_GOBS) ) return NULL_GOBS; }
  if( g2 && coef2 != 0.0 ){ C2 = obs_C( *g2 ); if( double_eq(C2,NULL_GOBS) ) return NULL_GOBS; }
  if( g3 && coef3 != 0.0 ){ C3 = obs_C( *g3 ); if( double_eq(C3,NULL_GOBS) ) return NULL_GOBS; }

  if (fabs(C1 - C2) > 100)
  {
	  if (_log)
		  _log->comment(2, "t_gobsgnss", "Inconsistencies of the code observations " + _epoch.str_ymdhms() );
	  return 0.0;
  }

  string sat = this->sat();

  double lc = ( coef1*C1 + coef2*C2 + coef3*C3 );

#ifdef DEBUG
   int b1, b2, b3;
   b1 = b2 = b3 = 0;
   if(g1) b1 = g1->band();    if(g2) b2 = g2->band();    if(g2) b2 = g2->band();   
    cout << _satid << " " << _staid
         << " _lc_range BAND(" << b1   << "," <<  b2 << " ," << b3 << ") "     
         <<          "  COEF(" << coef1 << "," << coef2 << "," << coef3 << ") "
         << fixed << setprecision(3)
         << setw(16) << C1 << setw(16) << C2 << setw(16) << C3 
         << "  lc = " << setw(16) << lc << endl;
#endif
  return lc;
}

// return value of general carrier-phase linear combination
// for 2 or 3 bands with given coefficients
// ----------
double t_gobsgnss::_lc_phase(const t_gobs* g1, const double& coef1,
	                     const t_gobs* g2, const double& coef2,
	                     const t_gobs* g3, const double& coef3 )
{
  double L1 = NULL_GOBS;
  double L2 = NULL_GOBS;
  double L3 = NULL_GOBS;

  if( g1 && coef1 != 0.0 ){ L1 = obs_L( *g1 ); if( double_eq(L1,NULL_GOBS) ) return NULL_GOBS; }
  if( g2 && coef2 != 0.0 ){ L2 = obs_L( *g2 ); if( double_eq(L2,NULL_GOBS) ) return NULL_GOBS; }
  if( g3 && coef3 != 0.0 ){ L3 = obs_L( *g3 ); if( double_eq(L3,NULL_GOBS) ) return NULL_GOBS; }

  double lc = ( coef1*L1 + coef2*L2 + coef3*L3 );

#ifdef DEBUG
   int b1, b2, b3;
   b1 = b2 = b3 = 0;
   if(g1) b1 = g1->band();    if(g2) b2 = g2->band();    if(g2) b2 = g2->band();
    cout << _satid << " " << _staid
         << " _lc_phase BAND(" << b1   << "," <<  b2 << " ," << b3 << ") "
         <<          "  COEF(" << coef1 << "," << coef2 << "," << coef3 << ") "
         << fixed << setprecision(3)
         << setw(16) << L1/wavelength(g1->band()) << setw(16) << L2/wavelength(g2->band()) << setw(16) << L3 
         << "  lc = " << setw(16) << lc << endl;
#endif
  return lc;
}


// get selected code (GOBS)
// ----------
GOBS t_gobsgnss::_id_range(GOBSBAND b)
{ 
  map<GOBS,double>::const_iterator it = _gobs.begin();

  while( it != _gobs.end() ){
    GOBS gobs = it->first;
    if( gobs_code( gobs ) ){
      GOBSBAND tmp = str2gobsband( gobs2str(gobs) );
      if( tmp == b ) return gobs;
    }
    it++;
  }
  return X;
}
   

// get selected phase (GOBS)
// ----------
GOBS t_gobsgnss::_id_phase(GOBSBAND b)
{ 
  map<GOBS,double>::const_iterator it = _gobs.begin();

  while( it != _gobs.end() ){
    GOBS gobs = it->first;
    if( gobs_phase( gobs ) ){
      GOBSBAND tmp = str2gobsband( gobs2str(gobs) );
      if( tmp == b ) return gobs;
    }
    it++;
  }
  return X;
}
   

// get selected code (GOBS)
// ----------
GOBS t_gobsgnss::_cod_id(const int& band) // , const GOBS& obs)
{    
  size_t sz = sizeof(code_choise[band]) / sizeof(GOBS);
  map<GOBS,double>::iterator it; // = _gobs.find( obs );

  GOBS tmp = X; // obs;

  // choices 
  for( size_t i=0; i<sz; ++i ){
    tmp = code_choise[band][i];
    
    it = _gobs.find( tmp );
    if( tmp == X || it == _gobs.end() ) continue;
       
    if( ! double_eq(it->second, NULL_GOBS ) ) return tmp;
//    { cout << " found  code ID: [" << band << "] " << gobs2str( tmp ) << endl; return tmp; }
  }

  return X;
}

// get selected phase (GOBS)
// ----------
GOBS t_gobsgnss::_pha_id(const int& band) // , const GOBS& obs)
{    
  size_t sz = sizeof(phase_choise[band]) / sizeof(GOBS);
  map<GOBS,double>::iterator it; // = _gobs.find( obs );
 
  GOBS tmp = X; // obs;

  // choices
  for( size_t i=0; i<sz; ++i ){
    tmp = phase_choise[band][i];
       
    it = _gobs.find( tmp );
    if( tmp == X || it == _gobs.end() ) continue;

    if( ! double_eq( it->second, NULL_GOBS) ) return tmp;
//    { cout << " found phase ID: [" << band << "] " << gobs2str( tmp ) << endl; return tmp; }
  }

  return X;
}

   
// return coefficients (c1,c2) of the ionosphere-free linear combination
// ----------
int t_gobsgnss::_coef_ionofree(GOBSBAND b1, double& c1, 
                               GOBSBAND b2, double& c2 )
{
  if( b1 == b2 ){ c1 = c2 = 0.0; return -1; }
   
  double fac1 = this->frequency( b1 );
  double fac2 = this->frequency( b2 );
  
  c1 =   fac1*fac1 / (fac1*fac1 - fac2*fac2);
  c2 = - fac2*fac2 / (fac1*fac1 - fac2*fac2);
   
#ifdef DEBUG
  cout << " coefs (D) = " << fixed << setprecision(3)
       << setw(15) << fac1
       << setw(15) << fac2
       << setw(15) << c1
       << setw(15) << c2
       << endl;
#endif

  return 0;
}

// return coefficients (c1,c2) of the ionosphere-free linear combination
// ----------
int t_gobsgnss::coef_ionofree(GOBSBAND b1, double& c1, 
                              GOBSBAND b2, double& c2 )
{
  return _coef_ionofree(b1,c1,b2,c2);
}


// return coefficients (c1,c2) of the geometry-free linear combination
// ----------
int t_gobsgnss::_coef_geomfree(GOBSBAND b1, double& c1, 
                               GOBSBAND b2, double& c2)
{   
  if( b1 == b2 ){ c1 = c2 = 0.0; return -1; }
  c1 =  1.0;
  c2 = -1.0;
  
  return 0;
}

// return coefficients (c1,c2) of the geometry-free linear combination
// ----------
int t_gobsgnss::coef_geomfree(GOBSBAND b1, double& c1, 
                              GOBSBAND b2, double& c2)
{
  return _coef_geomfree(b1,c1,b2,c2);
}


// return coefficients (c1,c2) of the narrow-lane linear combination
// ----------
int t_gobsgnss::_coef_narrlane(GOBSBAND b1, double& c1,
                               GOBSBAND b2, double& c2)
{   
  if( b1 == b2 ){ c1 = c2 = 0.0; return -1; }

  double fac1 = this->frequency( b1 );
  double fac2 = this->frequency( b2 );
  
  c1 =   fac1 / (fac1 + fac2);
  c2 =   fac2 / (fac1 + fac2);
  
  return 0;
}

// return coefficients (c1,c2) of the narrow-lane linear combination
// ----------
int t_gobsgnss::coef_narrlane(GOBSBAND b1, double& c1,
                              GOBSBAND b2, double& c2)
{
  return _coef_narrlane(b1,c1,b2,c2);
}

// return coefficients (c1,c2) of the wide-lane linear combination
// ----------
int t_gobsgnss::_coef_widelane(GOBSBAND b1, double& c1, 
                               GOBSBAND b2, double& c2)
{   
  if( b1 == b2 ){ c1 = c2 = 0.0; return -1; }

  double fac1 = this->frequency( b1 );
  double fac2 = this->frequency( b2 );
  
  c1 =   fac1 / (fac1 - fac2);
  c2 = - fac2 / (fac1 - fac2);
  
  return 0;
}

// return coefficients (c1,c2) of the wide-lane linear combination
// ----------
int t_gobsgnss::coef_widelane(GOBSBAND b1, double& c1,
                              GOBSBAND b2, double& c2)
{
  return _coef_widelane(b1,c1,b2,c2);
}


// -------------------------------------------------------------------------------------------
// t_obscmb class 
// -------------------------------------------------------------------------------------------
   
// operator< for t_obscmb
// --------------------------
bool t_obscmb::operator<(const t_obscmb& t) const
{
  return ( this->first.type() < t.first.type() &&
           this->first.band() < t.first.band() &&
           this->first.attr() < t.first.attr() );
}

} // namespace
