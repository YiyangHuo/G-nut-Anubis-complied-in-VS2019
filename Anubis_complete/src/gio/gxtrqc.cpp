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

#include <set>
#include <numeric>
#include <iostream>

#include "gio/gxtrqc.h"
#include "gio/gfile.h"
#include "gutils/gnss.h"
#include "gutils/gobs.h"
#include "gutils/gsys.h"
#include "gutils/gstat3d.h"
#include "gutils/gfileconv.h"
#include "gproc/gsppflt.h"
#include "gset/gsetqc.h"
#include "gset/gsetinp.h"
#include "gset/gsetproc.h"
#include "gprod/gprodclk.h"

using namespace std;
using namespace pugi;

namespace gnut {

double rad2deg = 180/G_PI;
const char NOT_YET[] = "NOT YET IMPLEMENTED";

// Constructors
// ---------
 t_gxtrqc::t_gxtrqc()
{
}
   
t_gxtrqc::t_gxtrqc(t_gsetbase* set,
             const string& pgm, const t_gtime& dt)
 : t_gxtr(set, pgm, dt),
   _gobs(0),
   _xtr(0),
   _xyz_est(0.0,0.0,0.0), // +=  used !
   _xyz_rep(0.0,0.0,0.0), // +=  used !
   _smp_est(0.0),
   _min_ele(90.0),
   _cut_ele(0.0),
   _hours(0.0),
   _minInt(0.0),
   _maxInt(0.0),
   _pos_cut(0.0),
   _pos_int(0),
   _xml(false),
   _kinematic(true),
   _ele_new(true),
   _ele_app(false),
   _ngap(0),
   _npcs(0),
   _summ(1),
   _head(1),
   _stat(1),
   _band(1),
   _gaps(1),
   _prep(1),
   _elev(1),
   _mult(1),
   _calc(1),
   _stnr(1),
   _sinf(1),
   _gkpi(0),
   _step(STP),
   _tgap(GAP),
   _tpcs(PCS),
   _mp_nepochs(MP_NEPOCHS),
//   _mp_all(MP_ALL),          // via settings only
   _mp_limit(MP_LIMIT),
   _auto_band(false),
   _ClkSync(0)
{
  gtrace("t_gxtrqc::t_gxtrqc");

  _get_settings();
  _qcdata = make_shared<t_gqcdata>();

}

// Destructor
// ---------
t_gxtrqc::~t_gxtrqc()
{
  gtrace("t_gxtrqc::desctructor");
  if( _xtr ){ if( _xtr->is_open() ) _xtr->close(); delete _xtr; }
}


// set data & gnav
// ---------
void t_gxtrqc::setDAT(t_galloqc* gobs, t_gallnav* gnav)
{   
  gtrace("t_gxtrqc::setDAT");
  _gobs = gobs;
  _gnav = gnav;
}       
   
// Extract data from galloqc
// ----------------------------
void t_gxtrqc::summary(string site)
{
  gtrace("t_gxtrqc::summary");
  if( site.empty() && _site.empty() ) return;

#ifdef DEBUG   
  t_gtime epoch(t_gtime::GPS); epoch.from_str("%Y-%m-%d %H:%M:%S", "2013-04-15 17:30:00");
  vector<t_gsatdata> v = _gobs->obs("GOP7", epoch);
  cout << "START\n"; cout.flush();
  for(auto itV=v.begin(); itV!=v.end(); ++itV){ cout << itV->sat() << itV->epoch().str_ymdhms(" ") << endl;}
  cout << "END\n"; cout.flush(); exit(1);
#endif   

  _clear();
   
  if( ! site.empty() ) _site = site;

  ostringstream osHDR, osCAL, osSUM, osOBS, osGAP, osBND, osPRE, osELE, osSNR, osMPH, osINF, osKPI;

  // get first sync obs
  _obs_sync = _gobs->beg_obs(_site, _step);
  _obs_beg  = _gobs->beg_obs(_site);
  _obs_end  = _gobs->end_obs(_site);
  _beg      = ( _set_beg != FIRST_TIME ) ? _set_beg : _obs_beg;
  _end      = ( _set_end != LAST_TIME  ) ? _set_end : _obs_end;
  _ref      = _beg; //  + (_beg - _end)/2;   // ref time

  this->_setOut();
  string mask = ""; if( _xtr ){ mask = _xtr->mask(); }
  if( _log ){ _log->comment(0,"gxtrqc", "SITE: " + _site + " " + mask); }
  if( _log && _obs_beg != _obs_sync ){
    _log->comment(1,"gxtrqc", "Sync begin "+_obs_beg.str_ymdhms()+_obs_sync.str_ymdhms(" -> ") );
  }

  if( _xml && _summ == 0 ) _summ = -1;
  if( _xml && _stat == 0 ) _stat = -1;
  if( _xml && _calc == 0 ) _calc = -1;
  if( _xml && _elev == 0 ) _elev = -1;
  if( _xml && _prep == 0 ) _prep = -1;
  if( _xml && _band == 0 ) _band = -1;
  if( _xml && _gaps == 0 ) _gaps = -1;
  if( _xml && _mult == 0 ) _mult = -1;
  if( _xml && _stnr == 0 ) _stnr = -1;

  _get_rnxhdr();

  // order is important!
  t_gtime t1,t2;
  t1 = t2 = t_gtime::current_time();

  _header(osHDR);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" header [" + int2str(_head) + "]: " + dbl2str(t2-t1) );

  if( _summ != 0 || _gkpi != 0 || _band != 0 )
  {
    if( _ele_new ){
      _gobs->est_elevations( _site, _pos, _gnav, _beg, _end, _cut_ele, _sat_view, _sat_mask );
    }else{
      sat_view( _sat_view, _pos, _gnav, _beg, _end,      0.0 );
      sat_view( _sat_mask, _pos, _gnav, _beg, _end, _cut_ele );
      if( _ele_app ) _gobs->apr_elevations( _site, _pos, _gnav );
      else           _gobs->add_elevations( _site, _pos, _gnav );
    }
  }
  
  if( _summ > 0 ) // needed for summary 1 & 2
  {
    t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" satview[" + int2str(_summ) + "]: " + dbl2str(t2-t1) );
    _stat_obs = _gobs->stat_obs( _site, &_smp_est, &_min_ele, _cut_ele, _end, &_sat_view, &_stat_smp, &_stat_ele );
    t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" obsview[" + int2str(_summ) + "]: " + dbl2str(t2-t1) );
  }else{
    _stat_obs = _gobs->stat_obs( _site, &_smp_est, &_min_ele, _cut_ele, _end, NULL, &_stat_smp, NULL );
    t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" obsstt [" + int2str(_summ) + "]: " + dbl2str(t2-t1) );
  }

  if( double_eq(_smp_est, 0.0) ){
    if( _log ) _log->comment(0,"gxtrqc","Warning: sampling could not be estimated!");
    else               cerr << "gxtrqc - Warning: sampling could not be estimated!\n";
    _smp_est = dynamic_cast<t_gsetgen*>(_set)->sampling_default();
  }

  _calcul(osCAL);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" estima [" + int2str(_calc) + "]: " + dbl2str(t2-t1) );
  _observ(osOBS);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" observ [" + int2str(_stat) + "]: " + dbl2str(t2-t1) );
  _nbands(osBND);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" nbands [" + int2str(_band) + "]: " + dbl2str(t2-t1) );
  _pieces(osGAP);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" pieces [" + int2str(_gaps) + "]: " + dbl2str(t2-t1) );
  _prepro(osPRE);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" prepro [" + int2str(_prep) + "]: " + dbl2str(t2-t1) );
  _skyplt(osELE);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" skyplt [" + int2str(_elev) + "]: " + dbl2str(t2-t1) );
  _mlpath(osMPH);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" mlpath [" + int2str(_mult) + "]: " + dbl2str(t2-t1) );
  _snoise(osSNR);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" snoise [" + int2str(_stnr) + "]: " + dbl2str(t2-t1) );
  _svinfo(osINF);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" svinfo [" + int2str(_sinf) + "]: " + dbl2str(t2-t1) );  
  _summar(osSUM);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtrqc",_site+" summar [" + int2str(_summ) + "]: " + dbl2str(t2-t1) );
  _grckpi(osKPI);

  if( _xtr ){
    _xtr->write("# " + _pgm + "\n"); // pgm version
    _xtr->write(osSUM.str()); // summary
    _xtr->write(osHDR.str()); // meta-data
    _xtr->write(osCAL.str()); // calc position 
    _xtr->write(osOBS.str()); // statistics
    _xtr->write(osBND.str()); // system bands
    _xtr->write(osGAP.str()); // gap & pieces
    _xtr->write(osPRE.str()); // preprocess
    _xtr->write(osELE.str()); // skyplot  
    _xtr->write(osMPH.str()); // multipath
    _xtr->write(osSNR.str()); // snr
    _xtr->write(osINF.str()); // sat info
    _xtr->write(osKPI.str()); // grc kpi
  }
   
  if( _xml ){
    this->_xml_meta();
    this->_xml_head();
    this->_xml_data();
    this->_xml_navi();
    this->_xml_full();
    make_path( _name );
    this->write( _name );
  }
}   


// ===================
// PROTECTED FUNCTIONS
// ===================
   
// Get settings
// -----------------------------
void t_gxtrqc::_get_settings()
{
  gtrace("t_gxtrqc::_get_settings");

  t_gxtr::_get_settings();
   
  _summ             = dynamic_cast<t_gsetqc*>(_set)->summ();
  _head             = dynamic_cast<t_gsetqc*>(_set)->head();
  _stat             = dynamic_cast<t_gsetqc*>(_set)->stat();
  _calc             = dynamic_cast<t_gsetqc*>(_set)->calc();
  _gaps             = dynamic_cast<t_gsetqc*>(_set)->gaps();
  _band             = dynamic_cast<t_gsetqc*>(_set)->band();
  _prep             = dynamic_cast<t_gsetqc*>(_set)->prep();
  _elev             = dynamic_cast<t_gsetqc*>(_set)->elev();
  _mult             = dynamic_cast<t_gsetqc*>(_set)->mult();
  _stnr             = dynamic_cast<t_gsetqc*>(_set)->stnr();   
  _sinf             = dynamic_cast<t_gsetqc*>(_set)->sinf();     
  _step             = dynamic_cast<t_gsetqc*>(_set)->step();
  _tgap             = dynamic_cast<t_gsetqc*>(_set)->tgap();
  _tpcs             = dynamic_cast<t_gsetqc*>(_set)->tpcs();
  _nsat             = dynamic_cast<t_gsetqc*>(_set)->nsat();
  _mp_limit         = dynamic_cast<t_gsetqc*>(_set)->mp_limit();
  _mp_all           = dynamic_cast<t_gsetqc*>(_set)->mp_all();
  _mp_nepochs       = dynamic_cast<t_gsetqc*>(_set)->mp_nepochs();
  _kinematic        = dynamic_cast<t_gsetqc*>(_set)->pos_kin();
  _ele_new          = dynamic_cast<t_gsetqc*>(_set)->ele_new();
  _ele_app          = dynamic_cast<t_gsetqc*>(_set)->ele_app();
  _cut_ele          = dynamic_cast<t_gsetqc*>(_set)->ele_cut();
  _pos_cut          = dynamic_cast<t_gsetqc*>(_set)->pos_cut();
  _pos_int          = dynamic_cast<t_gsetqc*>(_set)->pos_int();
  _sat_rec          = dynamic_cast<t_gsetqc*>(_set)->sat_rec();
  _useHealth        = dynamic_cast<t_gsetqc*>(_set)->useHealth();

   bool chkHealth   = dynamic_cast<t_gsetinp*>(_set)->chkHealth();

  _auto_band        = dynamic_cast<t_gsetproc*>(_set)->auto_band();
 
//  LOG still does not exist here
  if( chkHealth == false && _useHealth > POS_HEALTH ){
    cerr << "Warning: due to INP:chk_health=false the QC:use_health has been disabled\n";
  }
//if( chkHealth == false ){ _useHealth = POS_HEALTH; }
  
// DONT USE THIS, KEEP EVEN ZERO REQUEST FOR AUTO PROCESSING
//if( _smp_req == 0 ) _smp_req = dynamic_cast<t_gsetgen*>(_set)->sampling_default();  

  // pre-conditions
  if( !_stat           ) _stat = -1; // do always obs-statistics
  if( !_calc           ) _calc = -1; // do always calc for station position
  if( !_head           ) _head = -1; // do always header reading   
//if( !_stat && _summ > 0 ) _stat = -1; // always do obs-statistics for summary
//if( !_stat && _band > 0 ) _band = -1; // always do band summary
//if( !_stat && _gaps > 0 ) _gaps = -1; // always do gap summary
//if( !_stat && _mult > 0 ) _stat = -1; // always do obs-statistic for multipath
//if( !_prep && _mult > 0 ) _prep = -1; // always do preprocessing for multipath ! NOT NECESSARY ANYMORE!

  if( _pos.zero() ) {
    t_gtriple xyz;
    shared_ptr<t_grec> confObj = _rec->grec(_site);

    if( confObj != 0 ) xyz = confObj->crd(_obs_beg);
    if ( !xyz.zero() ) _pos = xyz;
  }   

  if( _mp_limit < 1 ){ _mp_limit = MP_LIMIT;
    cerr << "Warning: cycle-slip & outlier detection for MP too small [reset to " << MP_LIMIT << "]\n";
  }

  if( _mp_nepochs < 5 ){ _mp_nepochs = MP_NEPOCHS;
    cerr << "Warning: # of epochs for MP set too small [reset to " << MP_NEPOCHS << "]\n"; 
  }

#ifdef DEBUG
  cout << "settings: " << _summ << _head << _stat << _prep << _elev << _mult << endl;
#endif
}


// number of gaps
// ----------
void t_gxtrqc::_list_gap(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_gap");
  map<t_gtime, int> len;
  if( _gobs->gaps(_site, _tgap, len, _minInt, _maxInt) > 0 )
  {
    // FILL DATA MEMBER!
    _ngap = len.size();
     
    // verbosity
    if( _gaps < 1 ) return;
     
    _subsec(os,"LIST_GAP");
    _setkey(os,"GAPLST", '#', _ref.str_ymd() + " begTime");
    os << W(10) << "endTime" << W(8) << ">"+int2str(_tgap)+"s" << endl;

    // verbosity
    if( _gaps < 2 ) return;
     
    for(auto it = len.begin(); it != len.end(); ++it )
    {
      t_gtime tt( it->first );
      _setkey( os, "GAPLST", ' ',tt.str_ymdhms());
      tt.add_secs(it->second);
      os << W(10) << tt.str_hms() << W(8) << it->second << endl;
    }
  }
}


// number of blocks
// ----------
void t_gxtrqc::_list_pcs(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_pcs");
  map<t_gtime, int> len;
  if( _gobs->block(_site, _tpcs, _tgap, len) > 0 )
  {
    // FILL DATA MEMBER!
    _npcs = len.size();
     
    // verbosity
    if( _gaps < 1 ) return;

    _subsec(os,"LIST_PCS");
    _setkey(os,"PCSLST", '#', _ref.str_ymd() + " begTime");
    os << W(10) << "endTime" << W(8) << "<"+int2str(_tpcs)+"s" << endl;

    // verbosity
    if( _gaps < 2 ) return;

    for(auto it = len.begin(); it != len.end(); ++it )
    {
      t_gtime tt( it->first );
      _setkey( os, "PCSLST", ' ',tt.str_ymdhms());
      tt.add_secs(it->second);
      os << W(10) << tt.str_hms() << W(8) << it->second << endl;
    }
  }
}

// sampling histogram
// ----------
void t_gxtrqc::_list_smp(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_smp");

  _subsec(os,"LIST_SMP");
  _setkey(os,"SMPLST", '#');
  os << W(10) << "sampling" << W(8) << "count" << endl;

  for(auto it = _stat_smp.begin(); it != _stat_smp.end(); ++it )
  {
    _setkey(os,"SMPLST", ' ');
    os << W(10) << it->first << W(8) << it->second << endl;
  }
}


// list GNSS sys
// ----------
void t_gxtrqc::_list_gsys(ostringstream& os)
{   
  gtrace("t_gxtrqc::_list_gsys");
  if( _stat < 1 ) return;

  ostringstream osGNS, osSAT;

  _subsec(os,"LIST_GSYS");
  _setkey(osGNS,"GNSSYS", '=');
  osGNS << W(8) << int2str(_stat_obs.size());

  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS)
  {
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);

    osGNS << W(4) << gs;

    _setkey(osSAT, gs + "SAT", '=');

    // assume all signals has the full set of satellites (although sometimes zero observation)
    auto itGOBS = itGSYS->second.begin();

    osSAT << W(8) << int2str( itGOBS->second.size() );

    int offset = 0;
    if(     gsys==SBS) offset = SBS_OFFSET;
    else if(gsys==QZS) offset = QZS_OFFSET;
    for(int i = 1+offset; i<=_nsat+offset ; ++i )
    {
      string prn = t_gsys::eval_sat( i, gsys );
      if( itGOBS->second.find( prn) == itGOBS->second.end() )
      osSAT << W(4) << "-";
      else osSAT << W(4) << prn;
    }
    osSAT << endl;
  }
  os << osGNS.str() << endl << endl;
  os << osSAT.str();
}

// list GNSS types from header
// ----------
void t_gxtrqc::_list_typh(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_typh");
   
  os << endl;
  for(auto itHDR = _m_rnxHdr.begin(); itHDR != _m_rnxHdr.end(); ++itHDR )
  {
    t_rnxhdr rnxhdr = itHDR->second;
   
    t_rnxhdr::t_obstypes typehdr = rnxhdr.mapobs();

    for(auto itG = typehdr.begin(); itG != typehdr.end(); ++itG ){
      char c = itG->first[0];
      string gs_str = t_gsys::char2str(c);
      _setkey( os, gs_str + "HDR", '=');
      os << W(8) << itG->second.size();
      for(auto itB = itG->second.begin(); itB != itG->second.end(); ++itB ){
        os << W(4) << gobs2str(itB->first);
      }
      os << endl;
    }
  }
}


// list GNSS types
// ----------
void t_gxtrqc::_list_type(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_type");
  if( _stat > 0 ) _subsec(os,"LIST_TYPE");

  // special maps for data
  ostringstream osV, osVV;
   
  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS )
  {
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);
     
    _setkey( os, gs + "OBS", '=');
    os << W(8) << _stat_obs[gsys].size();

    map<string,vector<string>> sat_types;                // local or global ?

    // obstype-specific satellite list
    for(auto itGOBS = itGSYS->second.begin(); itGOBS != itGSYS->second.end(); ++itGOBS )
    {
      GOBS gobs = itGOBS->first;
      string go = gobs2str(gobs);
       
      os << W(4) << go;   // gobs2str( itGOBS->first );

      _setkey( osV, gs + go, ' ');
      osV << W(8) << ""; // reserve
       
      int offset = 0;
      if(     gsys==SBS) offset = SBS_OFFSET;
      else if(gsys==QZS) offset = QZS_OFFSET;
      for(int i = 1+offset; i<=_nsat+offset ; ++i )
      {
        string prn = t_gsys::eval_sat( i, gsys );

        auto itSAT = itGOBS->second.find( prn );
        if( itSAT != itGOBS->second.end() && itSAT->second.first > 0 )
        {
          sat_types[prn].push_back(go);
          osV << W(4) << prn;  
        }else{  osV << W(4) << "-"; }
      }
      osV << endl;
    }

    // sat-specific observation list
    for(auto itTYP = sat_types.begin(); itTYP != sat_types.end(); ++itTYP )
    {
      _setkey( osVV, itTYP->first + "OBS", ' ');
      osVV << W(8) << itTYP->second.size();
      for(auto it = itTYP->second.begin(); it != itTYP->second.end(); ++it ){ 
        osVV << W(4) << *it;
      }
      osVV << endl;
    }
    os << endl;
  }
   
  // extra verbosity
  if( _stat > 1 ) os << endl << osV.str();
  if( _stat > 2 ) os << endl << osVV.str();
}


// list bands
// ----------
void t_gxtrqc::_list_band(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_band");
  if( _band > 0 ){
    _subsec(os,"LIST_BAND");   
    _legend(os,"GNSxEP","FewBand" );
  }

  ostringstream osC,osL,osV,osX;
  _m_UnusEpo.clear();
  _m_UnusSat.clear();
  _m_FineEpo.clear();

  t_galloqc::t_map_bsat satUnus; // sat-specific map of pair<C,L> for unusable epochs
  t_galloqc::t_map_bsat satFine; // sat-specific map of pair<C,L> for total epochs

  // FILL DATA MEMBER!
  _gobs->nbands( _site, _qcdata->qc_bnd.data);
  // loop GSYS
  for(auto itBG  = _qcdata->qc_bnd.data.begin();
           itBG != _qcdata->qc_bnd.data.end(); ++itBG)
  {
    GSYS gsys = itBG->first;
    string gs = t_gsys::gsys2str(gsys);

    _m_UnusSat[gsys] = pair<int,int>(0,0); // (C/L-obs) unusable satellite epochs (GNSS-based)
    _m_UnusEpo[gsys] = pair<int,int>(0,0); // (C/L-obs) unusable epochs (GNSS-based)
    _m_FineEpo[gsys] = 0;                  // Total usable epochs (GNSS-based)
      
    int nEpoUnusC = 0;  // total unusable epochs (due to # C-bands)
    int nEpoUnusL = 0;  // total unusable epochs (due to # L-bands)

    for(auto itBT  = itBG->second.begin();
             itBT != itBG->second.end(); ++itBT) // loop EPOCH
    {
      int nSatTotal = itBT->second.size(); 
      int nSatUnusC = 0; // number of satellites with C-bands<2 in a single epoch
      int nSatUnusL = 0; // number of satellites with L-bands<2 in a single epoch
     
//     ---> PROBLEM WITH CLOCK DRIFT (not consistent for _nbands(site,t) && (_nbands,t,SAMPLE)
//     _m_nSatEpo[itBT->first][gsys] = nSatTotal; 
//     
      for(auto itBS  = itBT->second.begin();
               itBS != itBT->second.end(); ++itBS) // loop SAT
      {
        string prn( itBS->first );
        if( itBS->second.first  < 2 ){ satUnus[prn].first++;  nSatUnusC++; nEpoUnusC++; }else{ satFine[prn].first++;  } // code
        if( itBS->second.second < 2 ){ satUnus[prn].second++; nSatUnusL++; nEpoUnusL++; }else{ satFine[prn].second++; } // phase
      }

      // SPECIAL VERBOSITY for all incomplete C/L-band epochs!
      if( (nSatUnusC > 0 || nSatUnusL > 0) && itBG->first != SBS && itBG->first != QZS)
      {      
        _setkey( osV, gs+"XBN", ' ', itBT->first.str_ymdhms() );
        osV << W(4) << nSatUnusC << W(4) << nSatUnusL << W(4) << nSatTotal;

        for(auto itBS  = itBT->second.begin();
                 itBS != itBT->second.end(); ++itBS ) // loop SAT
        {
          if( itBS->second.first < 2 || itBS->second.second < 2 ){
            osV << W(4) << itBS->first;
          }
        }
        osV << endl;
      }

      if (itBG->first == SBS || itBG->first == QZS) _m_FineEpo[gsys]++;
      else if( nSatTotal - nSatUnusC < 4 || nSatTotal - nSatUnusL < 4){
           if( nSatTotal - nSatUnusC < 4 ) _m_UnusEpo[gsys].first++;
           if( nSatTotal - nSatUnusL < 4 ) _m_UnusEpo[gsys].second++; 
      }else _m_FineEpo[gsys]++;
       
      //Fill qcdata
      _qcdata->qc_kpi.DF_XX[gsys] = _m_UnusEpo[gsys].first;
      _qcdata->qc_kpi.DF_OK[gsys] = _m_FineEpo[gsys];

      // FILL DATA MEMBER!
      _m_UnusSat[gsys].first  += nSatUnusC;  // (C-obs) unusable satellite epochs (GNSS-based)
      _m_UnusSat[gsys].second += nSatUnusL;  // (L-obs) unusable satellite epochs (GNSS-based)

    } // EPO
     
    if( _band < 0 ) continue;

    _setkey(osC, gs+"CEP", '='); osC << W(8) << nEpoUnusC;
    _setkey(osL, gs+"LEP", '='); osL << W(8) << nEpoUnusL;
      
    int offset = 0;
    if(     gsys==SBS) offset = SBS_OFFSET;
    else if(gsys==QZS) offset = QZS_OFFSET;
    for(int i = 1+offset; i<=_nsat+offset ; ++i )
    {
      string prn = t_gsys::eval_sat( i, itBG->first);

      if( satFine.find(prn) != satFine.end() ){
        osC << W(4) << floor(   double(satFine[prn].first)
               / double(satUnus[prn].first  + satFine[prn].first  )*100.0);
        osL << W(4) << floor(   double(satFine[prn].second) 
               / double(satUnus[prn].second + satFine[prn].second )*100.0);
      }else if( satUnus.find(prn) != satUnus.end() ){
        osC << W(4) << 0;
        osL << W(4) << 0;
      }else{
        osC << W(4) << MISS;
        osL << W(4) << MISS;
      }
    }
    if( _band < 1 ) continue;
    os << osC.str() << endl; osC.str(""); osC.clear();
    os << osL.str() << endl; osL.str(""); osL.clear();

  } // GSYS
   
  // extra verbosity
  if( _band < 2 )  return;

  osC.str("");osC.clear();
  osL.str("");osL.clear();

  _qcdata->qc_bnd.data.clear();
  _gobs->nbands(_site, _qcdata->qc_bnd.data, _step); // different sampling from above!  Zero cutoff!

  // loop GSYS
  for(auto itBG  = _qcdata->qc_bnd.data.begin(); 
           itBG != _qcdata->qc_bnd.data.end(); ++itBG)
  {     
    GSYS gsys = itBG->first;
    string gs = t_gsys::gsys2str(gsys);

    // loop TIME
    for(auto itBT  = itBG->second.begin();
             itBT != itBG->second.end(); ++itBT)
    {
       // duplicite get 
       _m_nSatEpo[itBT->first][gsys] = itBT->second.size();
 
       _setkey( osC, gs+"CBN", ' ', itBT->first.str_ymdhms() ); osC << W(8) << itBT->second.size(); // << nSat;
       _setkey( osL, gs+"LBN", ' ', itBT->first.str_ymdhms() ); osL << W(8) << itBT->second.size(); // << nSat;

      // loop SAT (fixed)
      int offset = 0;
      if(     gsys==SBS) offset = SBS_OFFSET;
      else if(gsys==QZS) offset = QZS_OFFSET;
      for(int i = 1+offset; i<=_nsat+offset ; ++i )
      {   
        string prn = t_gsys::eval_sat( i, gsys);
//      cout << " " << prn;
        auto itTMP = itBT->second.find( prn );
        if( itTMP != itBT->second.end() ){
          osC << W(4) << itTMP->second.first;
          osL << W(4) << itTMP->second.second;

          // COPY TO NEW STRUCTURE !
          _qcdata->qc_tim.dat[gsys][itBT->first][prn].cbn = itTMP->second.first;
          _qcdata->qc_tim.dat[gsys][itBT->first][prn].lbn = itTMP->second.second;

        }else{
          osC << W(4) << MISS;
          osL << W(4) << MISS;
        }
      }
//      cout << endl;
      osC << endl;
      osL << endl;
    }
  }

  _legend( osX, "NxBAND","nSatell" );
  os << endl << osX.str() << osC.str() << osL.str();
  osX.str("");osX.clear();

  // single verbosity
  if( _band < 3 ) return;

  _setkey(osX, "GNSBND",'#');
  os << endl << osX.str()
     << W(4) << "C<2" 
     << W(4) << "L<2"
     << W(4) << "nSv"
     << W(7) << "SvList"
     << endl << osV.str();
  return;
}


// list cycle slips (verb 2)
// ---------------------
void t_gxtrqc::_list_cslp(ostringstream& os)
{     
  gtrace("t_gxtrqc::_list_cslp");

  if( _prep == 0 ) return;

  set<GOBS> gobs_pha;
  set<GSYS> gsys_avail;

  // loop GSYS
  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS )
  {
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);
    gsys_avail.insert(gsys);
     
    // loop GOBS
    for(auto itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS )
    {              
      GOBS gobs = itGOBS->first;
      string go = gobs2str( gobs );
      if( gobs_phase(gobs) ) gobs_pha.insert(gobs);
    }
  }

  // ALL slips listing
  for(auto itEP = _m_SlipsGap.begin(); itEP != _m_SlipsGap.end(); ++itEP )
  {     
    t_gtime t = itEP->first;
    map<string, map<GOBS, int> > mepo = itEP->second;

    for(auto itSAT = mepo.begin(); itSAT != mepo.end(); ++itSAT )
    {
      string gs = t_gsys::char2str( itSAT->first[0] );
      GSYS gsys = t_gsys::str2gsys( gs );

      for(auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it )
      {
        GOBS gobs = it->first;
        _m_totAll[gsys][gobs]++;

        if(      it->second == 1 ) _m_totEpo[gsys][gobs]++;
        else if( it->second == 2 ) _m_totSat[gsys][gobs]++;
        else if( it->second == 3 ) _m_totSig[gsys][gobs]++;
      }
    }
  }

  if( _prep > 0 ){   
    _subsec(os,"LIST_SLIPS");
    _setkey(os, "GNSSLP",'#');
    os << W(5) << "PRN";
    for(auto it = gobs_pha.begin(); it != gobs_pha.end(); ++it ){
      string go = gobs2str( *it );
      os << W(13) << trim(go);
    }
    os << endl;
  }

  for(auto itGSYS = gsys_avail.begin(); itGSYS != gsys_avail.end(); ++itGSYS )
  {
    for(auto itEPO = _m_Slips.begin(); itEPO != _m_Slips.end(); ++itEPO )
    {
      t_gtime t = itEPO->first;          
      map<string, map<GOBS, double> > mepo = itEPO->second;
      
      for(auto itSAT = mepo.begin(); itSAT != mepo.end(); ++itSAT )
      {
        string gs = t_gsys::char2str( itSAT->first[0] );
        GSYS gsys = t_gsys::str2gsys( gs );
   
        if( gsys != *itGSYS ) continue;
        
        for(auto it = itSAT->second.begin(); it != itSAT->second.end(); ++it )
        { 
          GOBS gobs = it->first;
          _m_totAll[gsys][gobs]++;
          _m_totSlp[gsys][gobs]++;
        }
      }
       
      if( _prep < 3 ) continue;

      for(auto itSAT = mepo.begin(); itSAT != mepo.end(); ++itSAT )
      {
        string prn = itSAT->first;
        string gs = t_gsys::char2str( prn[0] );
        GSYS gsys = t_gsys::str2gsys( gs );
        
        if( gsys != *itGSYS ) continue;
        
        _setkey(os, gs +"SLP", ' ', t.str_ymdhms());
        map<GOBS, double> mobs = itSAT->second;
        os << fixed << setprecision(0) << W(5) << prn;

        for(auto it = gobs_pha.begin(); it != gobs_pha.end(); ++it )
        {
          auto itGOBS = mobs.find(*it);
          if( itGOBS ==  mobs.end() ) os << W(13) << "-";
          else                        os << W(13) << itGOBS->second;
        }
        os << endl;
      }
    }
  }
}

// list receiver clk jumps
// --------------------------
void t_gxtrqc::_list_jump(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_jump");
  if( _prep < 2 ) return;

  _subsec(os,"LIST_JUMPS");
  _setkey( os, "CLKJMP", '#');
  os << W(11) << "Phase[ms]" << endl;
   
  if (_m_Breaks.size() == 0){
    _setkey(os, "CLKJMP", ' ');
    os << W(4) << "-" << endl;     
  }else{
    for(auto itBr = _m_Breaks.begin(); itBr != _m_Breaks.end(); itBr++)
    {
      t_gtime t = itBr->first;
      _setkey(os, "CLKJMP", ' ', t.str_ymdhms());
      os << W(4) << itBr->second << endl;
    }
  }
}


// list summ for preprocessing
// -----------------------------------
void t_gxtrqc::_list_prep(ostringstream& os)
{
  gtrace("t_gxtrqc::_list_prep");

  if( _prep == 0 ) return;

  GOBS gobs; GSYS gsys;

  for(auto itSYS = _m_totAll.begin(); itSYS != _m_totAll.end(); ++itSYS ){ gsys = itSYS->first;
    if( gsys == GNS ) continue;

    for(auto it = _m_totAll[gsys].begin(); it != _m_totAll[gsys].end(); ++it ){ 
      gobs = it->first; if( gobs == X ) continue;      
      _m_totAll[gsys][X] += it->second;
    }

    for(auto it = _m_totSlp[gsys].begin(); it != _m_totSlp[gsys].end(); ++it ){ 
      gobs = it->first; if( gobs == X ) continue;
      _m_totSlp[gsys][X] += it->second;
    }

    for(auto it = _m_totEpo[gsys].begin(); it != _m_totEpo[gsys].end(); ++it ){ 
      gobs = it->first; if( gobs == X ) continue;
      _m_totEpo[gsys][X] += it->second;
    }

    for(auto it = _m_totSat[gsys].begin(); it != _m_totSat[gsys].end(); ++it ){ 
      gobs = it->first; if( gobs == X ) continue;
      _m_totSat[gsys][X] += it->second;
    }

    for(auto it = _m_totSig[gsys].begin(); it != _m_totSig[gsys].end(); ++it ){ 
      gobs = it->first; if( gobs == X ) continue;
      _m_totSig[gsys][X] += it->second;
    }
  }
   
  if( _prep < 1 ) return;

  _subsec(os, "LIST_PREP_SUMM");
  _setkey(os, "GNSPRP", '#' );

  os << W(14) << "CS_Total"   << W(14) << "CS_Slip"
     << W(14) << "CS_Epoch"   << W(14) << "CS_Satell"
     << W(14) << "CS_Signal"  << endl;

  // now including GNS summary
  for(auto itSYS = _m_totAll.begin(); itSYS != _m_totAll.end(); ++itSYS )
  {
    GSYS gsys = itSYS->first;
    string gs = t_gsys::gsys2str(gsys);

    _setkey( os, gs + "PRP",'=');
    os << W(14) << _m_totAll[gsys][X]
       << W(14) << _m_totSlp[gsys][X]
       << W(14) << _m_totEpo[gsys][X]
       << W(14) << _m_totSat[gsys][X]
       << W(14) << _m_totSig[gsys][X]
       << endl;
  }

  os << endl;

  // SIGNAL-SPECIFIC OUTPUT
  _setkey( os, "GNSxxx", '#' );
  os << W(14) << "CS_Total"   << W(14) << "CS_Slip"
     << W(14) << "CS_Epoch"   << W(14) << "CS_Satell"
     << W(14) << "CS_Signal"  << endl;

  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS )
  {
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);

    for(auto itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS )
    {
      GOBS gobs = itGOBS->first;  
      string go = gobs2str(gobs);
       
      if( gobs_phase(gobs) ){     
        _setkey( os, gs+go, '=' );
        os << setprecision(0)
           << W(14) << _m_totAll[gsys][gobs]
           << W(14) << _m_totSlp[gsys][gobs]
           << W(14) << _m_totEpo[gsys][gobs]
           << W(14) << _m_totSat[gsys][gobs]
           << W(14) << _m_totSig[gsys][gobs]
           << endl;
      }
    }
  }
}


// list summ
// ----------
void t_gxtrqc::_list_summ(ostringstream& os)
{   
  gtrace("t_gxtrqc::_list_summ");
  if( _summ > 0 ) _subsec(os,"LIST_SUMM");

  // EPOCHS
  map<GSYS, pair<int,int> > nepo;
  _gobs->nepochs(_site, _obs_beg, _obs_end, _smp_est, nepo); // TIME CONSUMING

  // CLOCK JUMPS (clock synchronization)
  for(auto it = _stat_smp.begin(); it != _stat_smp.end(); ++it ){
    double tst = it->first - _smp_est;
    if( tst != 0.0 && fabs(tst) < _smp_est/100.0 ){ _ClkSync += it->second; } // unexpected interval
  }

  char tmp[3]; sprintf(tmp, "%02.0f", _cut_ele); string cut(tmp);
   
  double ElevDepCorrection = 1.0 - _cut_ele/90.0;

  // if High-rate (> 1Hz) values are normalized in Summary section
  if( _smp_est < 1.0 ){
     _setkey( os, "NOTICE", '#' ); os << " -----------------------------------------------------------\n";
     _setkey( os, "NOTICE", '#' ); os << " High-rate data - epoch-wise sumaries normalized to 1Hz data\n";
     _setkey( os, "NOTICE", '#' ); os << " -----------------------------------------------------------\n";
  }

  // Summary counts - expected observations
  for(auto itPRN = _sat_view.begin(); itPRN != _sat_view.end(); ++itPRN ){
    string prn = itPRN->first;
    GSYS gsys  = t_gsys::str2gsys( t_gsys::char2str( prn[0] ) );
    if( _stat_obs.find(gsys) == _stat_obs.end() ) continue;
    
    for(auto itBEG = _sat_view[prn].begin(); itBEG != _sat_view[prn].end(); ++itBEG ){
      int tmp = (int)floor(abs(itBEG->second - itBEG->first) / _smp_est)+1;

      // since 2.0.0: expected counts / signal (counting only those sats providing signal      
      // since 2.2.0: enabled also counting all the signals (thought not tracked by GNSS)
      bool sat_ok = false; // check below only for satellites providing signal(s)
      for(auto itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS )
      {
        GOBS gobs = itGOBS->first;
        auto itSAT = itGOBS->second.find( prn );
        if( _sat_rec == false || (itSAT != itGOBS->second.end() && itSAT->second.first > 0) )
        {
          _m_obsExp[gsys][gobs] += tmp;
          _m_obsExp[GNS][gobs]  += tmp;
          sat_ok = true;
        }
      }
      if( sat_ok )
      {
        _m_obsExp[gsys][X] += tmp;
        _m_obsExp[GNS][X]  += tmp;
      }
    }
  }

  // Summary counts - expected observations (elevation cut-off)
  for( auto itPRN = _sat_mask.begin(); itPRN != _sat_mask.end(); ++itPRN )
  {
    string prn = itPRN->first;
    GSYS gsys  = t_gsys::str2gsys( t_gsys::char2str( prn[0] ) );
    if( _stat_obs.find(gsys) == _stat_obs.end() ) continue;

    for(auto itBEG = _sat_mask[prn].begin(); itBEG != _sat_mask[prn].end(); ++itBEG )
    {
      int tmp = (int)floor(abs(itBEG->second - itBEG->first) / _smp_est) + 1;

      // since 2.0.0: expected counts / signal (counting only those sats providing signal
      // since 2.2.0: enabled also counting all the signals (thought not tracked by GNSS)
      bool sat_ok = false; // check below only for satellites providing signal(s)
      for( auto itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS )
      {
        GOBS gobs = itGOBS->first;
        auto itSAT = itGOBS->second.find( prn );
        if( _sat_rec == false || (itSAT != itGOBS->second.end() && itSAT->second.first > 0) )
        {
          _m_cutExp[gsys][gobs] += tmp;
          _m_cutExp[GNS][gobs]  += tmp;
          sat_ok = true;
        }
      }
      if( sat_ok )
      {
        _m_cutExp[gsys][X] += tmp;
        _m_cutExp[GNS][X]  += tmp;
      }
    }
  }
  
  // Summary counts - existent observations (elevation cut-off)
  _hours = 0.0;

  // System loop
  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS ){
    GSYS gsys = itGSYS->first; int slps = 0, zero=0, mask=0, xEle=0, oEle=0;

    if( _hours < nepo[gsys].second ) _hours = nepo[gsys].second;
     
    // Observation loop
    GOBS max = X;
    for(auto itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS ){
      GOBS gobs = itGOBS->first;

      // Elevation loop (existing elevations only !)
      for(auto itELE  = _stat_ele[gsys][gobs].begin();
      itELE != _stat_ele[gsys][gobs].end(); ++itELE ) _m_okElev[gsys][gobs] += itELE->second;

      // Satellite loop
      for(auto itSAT = _stat_obs[gsys][gobs].begin(); itSAT != _stat_obs[gsys][gobs].end(); ++itSAT ){
        if( itSAT->second.first > 0 ) _m_satHav[gsys][gobs] += 1;             // sum of all have satellites

        _m_obsHav[GNS][gobs]  += itSAT->second.first;
        _m_cutHav[GNS][gobs]  += itSAT->second.second;
        _m_obsHav[gsys][gobs] += itSAT->second.first;                         // sum of all have obs for SATs
        _m_cutHav[gsys][gobs] += itSAT->second.second;                        // sum of all have obs for SATs (cut-off)
      }
       
      _m_woElev[gsys][gobs]  = _m_obsHav[gsys][gobs] - _m_okElev[gsys][gobs]; // sum of all have obs wout Elevations (total)
       
      // select the highest # of phase observations (horizon)
      if( gobs_phase(gobs) )
        if( max == X || _m_obsHav[gsys][max] < _m_obsHav[gsys][gobs] ) max = gobs;
     
      // Add all wout elevation to user defined cut-off
      _m_cutHav[gsys][gobs] += int( _m_woElev[gsys][gobs] * ElevDepCorrection);
       
      // Add all signal-specific observations without elevations
      _m_obsExp[gsys][gobs] += int( _m_woElev[gsys][gobs]                    ); // zero mask
      _m_cutExp[gsys][gobs] += int( _m_woElev[gsys][gobs] * ElevDepCorrection); // user mask
    }

    slps = _m_totSlp[gsys][max];
    zero = _m_obsHav[gsys][max];
    mask = _m_cutHav[gsys][max];
    xEle = _m_woElev[gsys][max]; 
    oEle = _m_okElev[gsys][max];
     
    // Sum of highest counts among system-specific signals
    _m_obsHav[GNS][X]  += zero;                                   // horizon (all included)
    _m_cutHav[GNS][X]  += mask; // ALREADY CORRECTED!             // cut-off (xEle already corrected above)
    _m_obsExp[GNS][X]  += xEle;                                   // horizon,                   woutEle added
    _m_cutExp[GNS][X]  +=        int(xEle * ElevDepCorrection);   // cut-off,(only withEle), correction added
    _m_woElev[GNS][X]  += xEle;                                   // sum of all have obs wout Elevations (total over maxobs)
    _m_okElev[GNS][X]  += oEle;                                   // sum of all have obs with Elevations (total over maxobs)
    _m_totSlp[GNS][X]  += slps;                                   // sum of all cycle-slips (total over maxobs)
        
    // Add all wout elevation to sys-specific (Expected from highest count!)
    _m_obsHav[gsys][X] += zero;                                   // horizon (all included)
    _m_cutHav[gsys][X] += mask; // ALREADY CORRECTED!             // cut-off (xEle already corrected above)
    _m_obsExp[gsys][X] += xEle;                                   // horizon,                   woutEle added
    _m_cutExp[gsys][X] +=        int(xEle * ElevDepCorrection);   // cut-off (only withEle), correction added
    _m_woElev[gsys][X] += xEle;                                   // sum of all have obs wout Elevations
    _m_okElev[gsys][X] += oEle;                                   // sum of all have obs with Elevations

#ifdef DEBUG
    cerr << gsys << fixed << setprecision(0)
         << setw(6) << gobs2str( max ) << " : "
         << setw(6) << slps << "   "
         << setw(6) << zero << " -> "
         << setw(6) << mask << " + "
         << setw(6) << xEle << " -> " 
         << setw(6) << int(xEle*ElevDepCorrection)
         << " oE:" << setw(6) << _m_obsExp[gsys][X]
         << " oH:" << setw(6) << _m_obsHav[gsys][X]
         << " cE:" << setw(6) << _m_cutExp[gsys][X]
         << " cH:" << setw(6) << _m_cutHav[gsys][X]
         << "    "
         << " oE:" << setw(6) << _m_obsExp[GNS][X]
         << " oH:" << setw(6) << _m_obsHav[GNS][X]
         << " cE:" << setw(6) << _m_cutExp[GNS][X]
         << " cH:" << setw(6) << _m_cutHav[GNS][X]
         << endl;
#endif
  }

  // TOTAL SUMMARY (aka TEQC)
  _setkey( os, "TOTSUM", '#', "First_Epoch________" );
  os << W(20) << "Last_Epoch_________"
     << W(7)  << "Hours_"   << W(7) << "Sample"   << W(7) << "MinEle"
     << W(7)  << "#_Expt"   << W(7) << "#_Have"   << W(7) << "%Ratio"
     << W(7)  << "o/slps"   << W(7) << "woElev"
     << W(7)  << "Exp>"+cut << W(7) << "Hav>"+cut << W(7) << "%Rt>"+cut
     << endl;

  t_gtime beg = _obs_beg; // from DATA (already filtered)
  t_gtime end = _obs_end; // from DATA (already filtered)

  // collect filtered data ===> LIKE for XML-QC !!!!!
  for(auto itHDR = _m_rnxHdr.begin(); itHDR != _m_rnxHdr.end(); ++itHDR )
  {   
    t_rnxhdr rnxhdr = itHDR->second;
    t_gallobs::t_xfilter xflt = _gobs->xdata(_site, rnxhdr.path());
    if( xflt.beg < beg ) beg = xflt.beg; // use unfiltered BEG epoch
    if( xflt.end > end ) end = xflt.end; // use unfiltered END epoch    
    
  }

  _setkey( os, "TOTSUM", '=', beg.str_ymdhms() ); // _obs_beg.str_ymdhms() );
  os << fixed << setprecision(2)
     << W(20) << end.str_ymdhms()                 // _obs_end.str_ymdhms()
     << W(7)  << (_hours * _smp_est)/3600
     << W(7)  << _smp_est
     << W(7)  << _min_ele;

  double ratH = 100.0 * _m_obsHav[GNS][X] / _m_obsExp[GNS][X];
  double ratU = 100.0 * _m_cutHav[GNS][X] / _m_cutExp[GNS][X];
  if( ratH > 100.0 ){ 
//    cout << "RATIO: " << _m_obsHav[GNS][X] << " " << _m_obsExp[GNS][X] << endl;
    _log->comment(0,"gxtr","expected observations and have/expt ratio fixed (horizon ratio:"+dbl2str(ratH,1)+"%)");
    _m_obsExp[GNS][X] = _m_obsHav[GNS][X];
    ratH = 100.0;
  }
  if( ratU > 100.0 ){
//    cout << "RATIO: " << _m_cutHav[GNS][X] << " " << _m_cutExp[GNS][X] << endl;
    _log->comment(0,"gxtr","expected observations and have/expt ratio fixed (cut-off ratio:"+dbl2str(ratU,1)+"%)");
    _m_cutExp[GNS][X] = _m_cutHav[GNS][X];
    ratU = 100.0;
  }

  if( _m_obsExp[GNS][X] > 0 && _m_okElev[GNS][X] > 0 )
       os << W(7) << setprecision(0) << (( _smp_est<1 ) ? _m_obsExp[GNS][X]*_smp_est : _m_obsExp[GNS][X])
          << W(7) << setprecision(0) << (( _smp_est<1 ) ? _m_obsHav[GNS][X]*_smp_est : _m_obsHav[GNS][X])
          << W(7) << setprecision(2) << ratH;
  else os << W(7) << MISS
          << W(7) << setprecision(0) << (( _smp_est<1 ) ? _m_obsHav[GNS][X]*_smp_est : _m_obsHav[GNS][X])
          << W(7) << MISS;
   
  if( _prep )
       os << W(7) << setprecision(2) << (( _m_totSlp[GNS][X] > 0 ) ? _m_obsHav[GNS][X]/_m_totSlp[GNS][X] : _m_obsHav[GNS][X]);
  else os << W(7) << MISS;
  
  os << W(7) << _m_woElev[GNS][X];

  if( _m_cutExp[GNS][X] > 0 && _m_okElev[GNS][X] > 0 )
       os << W(7) << setprecision(0) << (( _smp_est<1 ) ? _m_cutExp[GNS][X]*_smp_est : _m_cutExp[GNS][X])
          << W(7) << setprecision(0) << (( _smp_est<1 ) ? _m_cutHav[GNS][X]*_smp_est : _m_cutHav[GNS][X])
          << W(7) << setprecision(2) << ratU;
  else os << W(7) << MISS
          << W(7) << MISS
          << W(7) << MISS;

  os << endl << endl;

  // GNSS-SPECIFIC SUMMARY
  _setkey( os, "GNSSUM", '#' );
  os << W(18) << " Epoch_Statistics_"
     << W(24) << " Excl_Epochs&Satellites_"
     << W(52) << " CycleSlips/Interruptions_And_Other_Discontinuities"
     << W(48) << " Code_Multipath_Mean_Statistics_Over_All_Signals"
     << endl;

  _setkey( os, "GNSSUM", '#' );   
  os << W(6) << "ExpEp"     << W(6) << "HavEp"     << W(6) << "UseEp"
     << W(6) << "xCoEp"     << W(6) << "xPhEp"     << W(6) << "xCoSv"     << W(6) << "xPhSv"
     << W(7) << "csAll"     << W(7) << "csEpo"     << W(7) << "csSat"     << W(7) << "csSig"
     << W(6) << "nSlp"      << W(6) << "nJmp"      << W(6) << "nGap"      << W(6) << "nPcs"
     << W(6) << "mp1"       << W(6) << "mp2"       << W(6) << "mp3"       << W(6) << "mpx" 
     << W(6) << "mp5"       << W(6) << "mp6"       << W(6) << "mp7"       << W(6) << "mp8"
     << endl;

  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS )
  {
    GSYS gsys = itGSYS->first;

    _setkey( os, t_gsys::gsys2str(gsys) + "SUM", '=' );

    // total epochs (expected, existing, usable)
    _nEpoExpect[gsys] = nepo[gsys].first;
    _nEpoExists[gsys] = nepo[gsys].second;
    os << fixed << setprecision(0)
       << W(6)  << ((_smp_est<1)?( _nEpoExpect[gsys] * _smp_est ):( _nEpoExpect[gsys] ))
       << W(6)  << ((_smp_est<1)?( _nEpoExists[gsys] * _smp_est ):( _nEpoExists[gsys] ));

    // total (C/L) fine epochs
    if( _m_FineEpo.find(gsys) != _m_FineEpo.end() ){
          os << W(6) << ((_smp_est<1)?( _m_FineEpo[gsys] * _smp_est ):( _m_FineEpo[gsys] )); }
    else{ os << W(6) <<  MISS; }

    // total (C/L) unused epochs
    if( _m_UnusEpo.find(gsys) != _m_UnusEpo.end() ){
          os << W(6) << ((_smp_est<1)?( _m_UnusEpo[gsys].first  * _smp_est ):( _m_UnusEpo[gsys].first  ))
             << W(6) << ((_smp_est<1)?( _m_UnusEpo[gsys].second * _smp_est ):( _m_UnusEpo[gsys].second )); }
    else{ os << W(6) <<  MISS << W(6) <<  MISS; }

    // total (C/L) unused satellite-specific epochs
    if( _m_UnusSat.find(gsys) != _m_UnusSat.end() ){
          os << W(6) << ((_smp_est<1)?( _m_UnusSat[gsys].first  * _smp_est ):( _m_UnusSat[gsys].first  ))
             << W(6) << ((_smp_est<1)?( _m_UnusSat[gsys].second * _smp_est ):( _m_UnusSat[gsys].second )); }
    else{ os << W(6) <<  MISS << W(6) <<  MISS; }

    // slips, jumps and other interruptions
    if( _prep ){
          os << W(7) << _m_totAll[gsys][X] // total  interruptions or cycle slips
             << W(7) << _m_totEpo[gsys][X] // epoch  interruptions
             << W(7) << _m_totSat[gsys][X] // satell interruptions
             << W(7) << _m_totSig[gsys][X] // signal interruptions
             << W(6) << _m_totSlp[gsys][X] // found  cycle slips
             << W(6) << _m_Breaks.size() + _ClkSync; } // both types!
    else{ os << W(7) <<  MISS
             << W(7) <<  MISS
             << W(7) <<  MISS
             << W(7) <<  MISS
             << W(6) <<  MISS
             << W(6) <<  MISS; }

    // gaps, pices
    if( _gaps ){
          os << W(6) << _ngap
             << W(6) << _npcs; }
    else{ os << W(6) <<  MISS 
             << W(6) <<  MISS; }

   // mean Multipath for all bands
   const unsigned int nband = 9;
   unsigned int band;

   if( _mult ){
     double cb[nband];
     double av[nband];
     int    cn[nband];
     t_gtime tt(LAST_TIME);

     for(band=1; band<nband; ++band) av[band] = cb[band] = cn[band] = 0;
 
     for(auto itOBS = _m_mp[gsys].begin(); itOBS != _m_mp[gsys].end(); ++itOBS )
     {
       GOBS gobs = itOBS->first;
       band = gobs2band(gobs);
   
       // search only mean characteristics (represented by LAST_TIME and empty PRN)
       if( _m_mp[gsys][gobs].find(tt)     == _m_mp[gsys][gobs].end()     ) continue;
       if( _m_mp[gsys][gobs][tt].find("") == _m_mp[gsys][gobs][tt].end() ) continue;

       cb[band] += _m_mp[gsys][gobs][tt][""].first / _m_mp[gsys][gobs][tt][""].second;
       cn[band]++;     
     }
     
     for(band = 1; band<nband; ++band ){
       switch( band ){
       case 0: continue;
       case 9: continue;
       default:
         if( cn[band] != 0 ){
           av[band] = cb[band]/cn[band]; // cm
           os << fixed << setprecision(1) << W(6) << av[band];
         }else
           os << W(6) << MISS;
       }
     }
   }else{ os << W(6) << MISS << W(6) << MISS << W(6) << MISS << W(6) << MISS
             << W(6) << MISS << W(6) << MISS << W(6) << MISS << W(6) << MISS;
   }     
   os << endl;
  }   

  if( _summ < 2 ) return;

  os << endl;

  _setkey( os, "GNSxxx", '#' );
  os << W(5) << "nSat"
     << W(7) << "ExpObs"   << W(7) << "HavObs"   << W(7) << "%Ratio"
     << W(7) << "Exp>"+cut << W(7) << "Hav>"+cut << W(7) << "%Rt>"+cut
     << W(7) << "wo/Ele";
   
  if( _summ > 2 ){
    auto it = ELEVATIONS.begin();
//  while( ++it != ELEVATIONS.end() ) os << W(7) << int2str(*(prev(it)))+"-"+int2str(*it);  // histogram
    while( ++it != ELEVATIONS.end() ) os << W(7) << "Ele>"+int2str(int(*(prev(it))));
  }
   
  os << endl;

  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS )
  {
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);

    for(auto itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS )
    {
      GOBS gobs = itGOBS->first;
      string go = gobs2str(gobs);
      _setkey( os, gs+go, '=' );

      os << fixed << setprecision(0) << W(5) << _m_satHav[gsys][gobs];

      // count observations with/wout elevations!
      int ok_Elev = _m_okElev[gsys][gobs];           // Existant_WithElevations!
      int wo_Elev = _m_woElev[gsys][gobs];           // Existant_WoutElevations!
//    int ExpZero = _m_obsExp[gsys][X];              // Expected_Horizon + wo_Elev_All (already included!)
      int ExpZero = _m_obsExp[gsys][gobs];           // Expected_Horizon + wo_Elev_All (already included!)
      int HavZero = _m_obsHav[gsys][gobs];           // Existant_Horizon only
//    int ExpMask = _m_cutExp[gsys][X];              // Expected_CutMask + wo_Elev_Weighted (already included!)
      int ExpMask = _m_cutExp[gsys][gobs];           // Expected_CutMask + wo_Elev_Weighted (already included!)
      int HavMask = _m_cutHav[gsys][gobs];           // Existant_CutMask only

      // >1Hz sampling
      if( _smp_est < 1 ){ 
        ExpZero = (int)(ExpZero*_smp_est);  HavZero = (int)(HavZero*_smp_est);  wo_Elev = (int)(wo_Elev*_smp_est);
        ExpMask = (int)(ExpMask*_smp_est);  HavMask = (int)(HavMask*_smp_est);  ok_Elev = (int)(ok_Elev*_smp_est);
      }

      double ratH = 100.0 * HavZero/ExpZero; if( ratH > 100.0 ){ ratH = 100.0; ExpZero=HavZero; _m_obsExp[gsys][gobs] = _m_obsHav[gsys][gobs]; }
      double ratU = 100.0 * HavMask/ExpMask; if( ratU > 100.0 ){ ratU = 100.0; ExpMask=HavMask; _m_cutExp[gsys][gobs] = _m_cutHav[gsys][gobs]; }

      // support for a few elevations only
      if( ok_Elev > 0 && wo_Elev < _m_obsHav[gsys][gobs] )
            os << W(7) << setprecision(0) << ExpZero
               << W(7) << setprecision(0) << HavZero
               << W(7) << setprecision(2) << ratH     // 100.0 * HavZero/ExpZero
               << W(7) << setprecision(0) << ExpMask
               << W(7) << setprecision(0) << HavMask
               << W(7) << setprecision(2) << ratU     // 100.0 * HavMask/ExpMask
               << W(7) << setprecision(0) << wo_Elev;
      else  os << W(7) << MISS
               << W(7) << setprecision(0) << HavZero
               << W(7) << MISS
               << W(7) << MISS
               << W(7) << MISS
               << W(7) << MISS
               << W(7) << setprecision(0) << wo_Elev;

      // histogram of sat elevation observations
      if( _summ > 2 ){
    
        map<double,int> m_ele;
    
        // use predefined ELE-BINs (_stat_ele might not exists)
        for(auto it = ELEVATIONS.begin(); it != prev(ELEVATIONS.end()); ++it ){
          auto itELE = _stat_ele[gsys][gobs].lower_bound( int(*it) );
          while( itELE != _stat_ele[gsys][gobs].end() ) m_ele[*it] += (itELE++)->second;
          
          if( _m_okElev[gsys][gobs] > 0 )
            os << W(7) << (( _smp_est<1 )?( m_ele[*it] * _smp_est ):( m_ele[*it] ));
          else os << W(7) << MISS;
        }

/*      // histogram 
   itELE = _stat_ele[gsys][gobs].begin();
   for( unsigned int i = 1; i < ELEVATIONS.size(); ++i ){ // use predefined ELE-BINs (_stat_ele might not exists)
     if( itELE != _stat_ele[gsys][gobs].end() ){
           os << W(7) << (( _smp_est<1 )?( itELE->second * _smp_est ):( itELE->second )); itELE++;
     }else os << W(7) << MISS;
   }
*/
      }
      os << endl;
    }
  }

  if( _summ < 4 ) return;

  os << endl;  _setkey(os, "SKYxxx", '#', "Ascending_Horizon__");
  os << W(20) << "Descending_Horizon_"
     << W(8)  << "Time[h]" 
     << W(8)  << "ExptObs"
     << endl;

  for(auto itPRN = _sat_view.begin(); itPRN != _sat_view.end(); ++itPRN ){
    string prn = itPRN->first;

    GSYS  gsys = t_gsys::char2gsys(prn[0]);
    if( _stat_obs.find(gsys) == _stat_obs.end() ) continue;
       
    for(auto itBEG = _sat_view[prn].begin(); itBEG != _sat_view[prn].end(); ++itBEG ){
      _setkey(os, "SKY"+prn, '=', itBEG->first.str_ymdhms() );
      os << fixed << setprecision(3)
         << W(20) <<  itBEG->second.str_ymdhms()
         << W(8)  << (itBEG->second - itBEG->first)/3600      // duration [h]
         << setprecision(0)
         << W(8)  << (itBEG->second - itBEG->first)/_smp_est  // expected observations
         << endl;
    }
    if( _summ > 4 ){ // special verbosity for MASK
      for(auto itBEG = _sat_mask[prn].begin(); itBEG != _sat_mask[prn].end(); ++itBEG ){
        _setkey(os, "MSK"+prn, '=', itBEG->first.str_ymdhms() );
        os << fixed << setprecision(3)
           << W(20) <<  itBEG->second.str_ymdhms()
           << W(8)  << (itBEG->second - itBEG->first)/3600      // duration [h]
           << setprecision(0)
           << W(8)  << (itBEG->second - itBEG->first)/_smp_est  // expected observations
           << endl;
      }
    }
  }
}

// Section: summary
// ----------------------------
void t_gxtrqc::_summar(ostringstream& os)
{
  gtrace("t_gxtrqc::_summar");
  if( _summ > 0 ) _section(os, "Summary statistics", _summ);
  else if( !_summ ) return;

  _list_summ(os);
}


// Section: Header observations
// ---------
void t_gxtrqc::_header(ostringstream& os)
{
  gtrace("t_gxtrqc::_header");
  if( _head > 0 ) _section(os, "Header information", _head);
  else if( !_head ) return;
 
  if( _head < -1 ) return;

  if( _m_rnxHdr.size() == 0 ){
    _pos = _xyz_approx();
    os << "=No header information available!" << endl;
    return;
  }
  
  for(auto itHDR = _m_rnxHdr.begin(); itHDR != _m_rnxHdr.end(); ++itHDR )
  {
    t_rnxhdr rnxhdr = itHDR->second;
     
    // check input from multiple files
    if( _m_rnxHdr.size() > 1 ){ 
      if( _log ) _log->comment(1,"gxtrqc","Warning: input from multiple files/headers identified!");
      os << "\n# Multiple file input --> multi-header report"
         << "\n# ============================================\n\n";
    }

    _pos = rnxhdr.aprxyz();              // APPROX POS check
    t_gtriple code_est = _xyz_approx();

    if( !code_est.zero() ) {
       _pos = code_est;
       t_gtriple diff = _pos - code_est;
       if(diff.norm() > 100) {
         if( _log ) _log->comment(0,"gxtrqc","Warning: RINEX header position differ too much from code etimation.");   
       }
    }

    if( _head == -1 ) return;

    // SETTINGS INFORMATION
    shared_ptr<t_grec> setrec = _rec->grec(_site);
    if( setrec == 0 ){
      if( _log ) _log->comment(1,"gxtrqc","Warning: no user receiver settings "+_site);
//     else               cerr << "gxtrqc - Warning: no user receiver settings "+_site << endl;
      setrec = make_shared<t_grec>();
    }

    string tt(_ref.str_ymdhms() + " ");
//  _obs_beg = _gobs->beg_obs(_site);
//  _obs_end = _gobs->end_obs(_site);

    string first = rnxhdr.first().str_ymdhms();
    string  last = rnxhdr.last().str_ymdhms();
    string rundt = rnxhdr.gtime().str_ymdhms();
    
    if( rnxhdr.first() == FIRST_TIME ) first = "xxxx-xx-xx xx:xx:xx";
    if( rnxhdr.last()  == LAST_TIME  )  last = "xxxx-xx-xx xx:xx:xx";
    if( rnxhdr.gtime() == FIRST_TIME ) rundt = "xxxx-xx-xx xx:xx:xx";

    _setkey(os, "RNXHDR", '#', tt);
    os << W(25) << left << "_RINEX_HEADER___________"
       << W(25) << left << "_RINEX_HEADER___________"
       << W(25) << left << "_RINEX_HEADER___________" << endl;
 
    _setkey(os, "RNXVER", '=', tt);
    !rnxhdr.rnxver().empty()   ? os << W(25) << left << rnxhdr.rnxver()
                               : os << W(25) << left << "-";
     rnxhdr.rnxsys() != ' '    ? os << W(25) << left << rnxhdr.rnxsys()
                               : os << W(25) << left << "-";
                                 os << W(25) << left << rundt << endl;

    _setkey(os, "RNXPGM", '=', tt);
    !rnxhdr.program().empty()  ? os << W(25) << left << rnxhdr.program()
                               : os << W(25) << left << "-";
    !rnxhdr.runby().empty()    ? os << W(25) << left << rnxhdr.runby()  << endl
                               : os << W(25) << left << "-"             << endl;

    _setkey(os, "RNXAGE", '=', tt);
    !rnxhdr.observer().empty() ? os << W(25) << left << rnxhdr.observer()
                               : os << W(25) << left << "-";
    !rnxhdr.agency().empty()   ? os << W(25) << left << rnxhdr.agency() << endl
                               : os << W(25) << left << "-"             << endl;

    os << endl;

    _setkey(os, "RNXHDR", '#', tt);
    os << W(48) << left << "_RINEX_HEADER___________________________________" << "  "
       << W(48) << left << "_USER_REQUEST___________________________________" << endl;

    _setkey(os, "BEGEND", '=', tt);
    os << W(24) << left << first
       << W(24) << left << last
       << W(2)  << ""
       << W(24) << _set_beg.str_ymdhms()
       << W(24) << _set_end.str_ymdhms() << endl;

    _setkey(os, "INTHDR", '=', tt);
    os << fixed << setprecision(3)
       << W(24) << left << rnxhdr.interval()
       << W(24) << left << " "
       << "  "
       << fixed << setprecision(3)
       << W(24) << left << _smp_req
       << W(24) << left << ""
       << endl;

    _setkey(os, "MARKER", '=', tt);
    os << W(48) << left << rnxhdr.markname()+" "+rnxhdr.marknumb()
       << W(2)  << ""
       << W(48) << left << setrec->name()+" "+setrec->domes() << endl;

    _setkey(os, "RECEIV", '=', tt);
    os << W(20) << left << rnxhdr.rectype()
       << W(20) << left << rnxhdr.recvers()
       << W(8)  << left << rnxhdr.recnumb().substr(0,8)
       << W(2)  << " "
       << W(48) << left << setrec->rec(_obs_beg) << endl;
     
    _setkey(os, "ANTENN", '=', tt);
    os << W(20) << left << rnxhdr.anttype()
       << W(20) << left << rnxhdr.antnumb()
       << W(8)  << ""
       << W(2)  << ""
       << W(48) << left << setrec->ant(_obs_beg) << endl;

    os << endl;

    _setkey(os, "RNXHDR", '#', tt);
    os << W(48) << left << "_RINEX_HEADER___________________________________" << "  "
       << W(48) << left << "_USER_REQUEST___________________________________" << endl;

    _setkey(os, "XYZAPR", '=', tt);
    os << fixed << setprecision(4)
       << W(16) << right << rnxhdr.aprxyz().crd(0)
       << W(16) << right << rnxhdr.aprxyz().crd(1)
       << W(16) << right << rnxhdr.aprxyz().crd(2)
       << W(2)  << ""
       << W(16) << right << setrec->crd(_obs_beg).crd(0)
       << W(16) << right << setrec->crd(_obs_beg).crd(1)
       << W(16) << right << setrec->crd(_obs_beg).crd(2) << endl;

    _setkey(os, "XYZECC", '=', tt);
    os << fixed << setprecision(4)
       << W(16) << right << rnxhdr.antxyz().crd(0)
       << W(16) << right << rnxhdr.antxyz().crd(1)
       << W(16) << right << rnxhdr.antxyz().crd(2)
       << W(2)  << ""
       << W(16) << right << setrec->eccxyz(_obs_beg).crd(0)
       << W(16) << right << setrec->eccxyz(_obs_beg).crd(1)
       << W(16) << right << setrec->eccxyz(_obs_beg).crd(2) << endl;
     
    _setkey(os, "ENUECC", '=', tt);
    os << fixed << setprecision(4)
       << W(16) << right << rnxhdr.antneu().crd(1) // E ! reverse with N
       << W(16) << right << rnxhdr.antneu().crd(0) // N
       << W(16) << right << rnxhdr.antneu().crd(2) // U
       << W(2)  << ""
       << W(16) << right << setrec->eccneu(_obs_beg).crd(1)          // E ! reverse with N
       << W(16) << right << setrec->eccneu(_obs_beg).crd(0)          // N 
       << W(16) << right << setrec->eccneu(_obs_beg).crd(2) << endl; // U
  }
}

// Section: Estimated
// ----------------------------
void t_gxtrqc::_calcul(ostringstream& os)
{
  if( _gnav == 0 || _gnav->satellites().size() == 0 ){ _calc = 0; return; }

  gtrace("t_gxtrqc::_calcul");

  if( _calc > 0 ) _section(os, "Estimated values", _calc);
  else if (!_calc) return;   

  double sample = _pos_int; // 15 min - don't synchronize with other Sampling (AZI/ELE/MPT/SNR..)
  string tt(_ref.str_ymdhms() + " ");
  t_gtriple xyz, blh;
  double clk;
  
  if( _calc > 0 ){   
    _setkey(os, "PERIOD", '=', tt);
    string first = _obs_beg.str_ymdhms();
    string  last = _obs_end.str_ymdhms();

    if( _obs_beg == FIRST_TIME ) first = "xxxx-xx-xx xx:xx:xx";
    if( _obs_end == LAST_TIME  )  last = "xxxx-xx-xx xx:xx:xx";
   
    os << W(24) << left << first
       << W(24) << left << last << endl;

    _setkey(os, "SAMPLE", '=', tt);
    os << fixed << setprecision(3)
       << W(24) << left << _smp_est << endl;
  }

  map<GSYS, vector<t_gtriple> > m_xyz;
  map<GSYS, vector<t_gtriple> > m_blh;
  map<GSYS, vector<t_gtime> >   m_time;
  map<GSYS, vector<double> >    m_clk;
  map<GSYS, vector<t_gpair>>    m_nsat;
  struct s_dop {double gdop, pdop, hdop, vdop;};
  map<GSYS, vector<s_dop>>      m_gdop;
  
// SPP version
  t_gsppflt gspp(_site, _set);  gspp.minsat(5); 
  gspp.setDAT(_gobs, _gnav);
  t_gallprod prod;
  gspp.setOUT(&prod); //link the _allprod with the prod, modyfy one equals to modify another
  gspp.setOBJ(_gobj);
  gspp.glog(_log);

  t_gtriple temp(0.0, 0.0, 0.0);
  for(t_gtime t = _obs_sync; t < _obs_end || t == _obs_end; t.add_dsec(sample) )
  {
    set<GSYS> sys = _gobs->sys(_site, t);
    
    for(auto itSYS = sys.begin(); itSYS != sys.end(); ++itSYS )
    {
      gspp.setgnss(*itSYS);
      prod.clear();
      gspp.processBatch(t,t);

    // use synchronized time for clock drift in RINEX since gsppflt reports rounded values
      t_gtime tt(t); //tt.reset_dsec();    
      shared_ptr<t_gprodcrd> prd_crd = static_pointer_cast<t_gprodcrd>( prod.get(_site, t_gdata::POS, tt) ); //crd refers to cordinate
      shared_ptr<t_gprodclk> prd_clk = static_pointer_cast<t_gprodclk>( prod.get(_site, t_gdata::CLK, tt) );

      if( prd_crd != nullptr && prd_clk != nullptr )//prd_crd refers to predict coordinate
      {
        xyz = prd_crd->xyz();
        clk = prd_clk->clk();
        int nsat = prd_crd->nSat();
        int nsat_excl = prd_crd->nSat_excl();
        t_gpair p_nsat(nsat, nsat_excl);
        s_dop dop;

        prd_crd->get_val("GDOP", dop.gdop);
        prd_crd->get_val("PDOP", dop.pdop);
        prd_crd->get_val("HDOP", dop.hdop);
        prd_crd->get_val("VDOP", dop.vdop);

        if( xyz.zero() ){ 
          _qcdata->qc_est.epo_excl[*itSYS]++;
          if( _log ) _log->comment(0,"gxtrqc", _site + " XYZ is zero at " + tt.str_ymdhms());
        }else if( xyz2ell(xyz, blh, false) != 0 ){
          _qcdata->qc_est.epo_excl[*itSYS]++;
          if( _log ) _log->comment(0,"gxtrqc", _site + " XYZ can not be converted to BLH");
        }else{
          m_xyz[*itSYS].push_back(xyz);
          m_blh[*itSYS].push_back(blh);
          m_clk[*itSYS].push_back(clk); 
          m_time[*itSYS].push_back(tt);
          m_nsat[*itSYS].push_back(p_nsat);
          m_gdop[*itSYS].push_back(dop);

          _qcdata->qc_est.epo_used[*itSYS]++;
          _qcdata->qc_est.pos_line[*itSYS][tt].xyz = xyz;
          _qcdata->qc_est.pos_line[*itSYS][tt].blh = blh;
          _qcdata->qc_est.pos_line[*itSYS][tt].clk = clk;
          _qcdata->qc_est.pos_line[*itSYS][tt].nsat = nsat;
          _qcdata->qc_est.pos_line[*itSYS][tt].xsat = nsat_excl;

          _qcdata->qc_est.pos_line[*itSYS][tt].gdop = dop.gdop;
          _qcdata->qc_est.pos_line[*itSYS][tt].pdop = dop.pdop; 
          _qcdata->qc_est.pos_line[*itSYS][tt].hdop = dop.hdop; 
          _qcdata->qc_est.pos_line[*itSYS][tt].vdop = dop.vdop; 
          
          t_gtriple diff;
          diff = xyz - temp;

          if( fabs(diff[0]) < 20 &&
              fabs(diff[1]) < 20 &&
              fabs(diff[2]) < 20 && _pos.zero() ){
            _pos = xyz;
            if( _calc < 1 ) return;
          }else{ temp = xyz; }
        }
      }else{
        _qcdata->qc_est.epo_excl[*itSYS]++;
      }
    }   
  }
  
// mean xyz
  int nsys = 0;
  for(auto itSYS = m_xyz.begin(); itSYS != m_xyz.end(); ++itSYS ){
     GSYS gsys = itSYS->first;
     t_gstat3d stat(itSYS->second);
     stat.calc_stat();
     _m_xyz_est[gsys] = stat.get_mean3d();
     _m_xyz_rep[gsys] = stat.get_std3d();
     int size = stat.get_size();
     int outl = stat.get_outl();

     if( _m_xyz_rep[gsys][0] > 9999 || _m_xyz_rep[gsys][1] > 9999 || _m_xyz_rep[gsys][2] > 9999 ){
       _m_xyz_est[gsys][0] = _m_xyz_est[gsys][1] = _m_xyz_est[gsys][2] =     0.0;
       _m_xyz_rep[gsys][0] = _m_xyz_rep[gsys][1] = _m_xyz_rep[gsys][2] = -9999.0;
     }else{
       _xyz_est[0] += _m_xyz_est[gsys][0];  _xyz_rep[0] += _m_xyz_rep[gsys][0];
       _xyz_est[1] += _m_xyz_est[gsys][1];  _xyz_rep[1] += _m_xyz_rep[gsys][1];
       _xyz_est[2] += _m_xyz_est[gsys][2];  _xyz_rep[2] += _m_xyz_rep[gsys][2];
       ++nsys;
     }

     if( _calc < 1 ) continue;

     string str = "XYZ" + t_gsys::gsys2str(itSYS->first);
     _setkey(os, str, '=',tt);
     os << fixed << setprecision(4)
        << W(16) << right << _m_xyz_est[gsys][0]
        << W(16) << right << _m_xyz_est[gsys][1]
        << W(16) << right << _m_xyz_est[gsys][2]
                 << setprecision(1)
        << W(8)  << right << _m_xyz_rep[gsys][0]
        << W(8)  << right << _m_xyz_rep[gsys][1]
        << W(8)  << right << _m_xyz_rep[gsys][2]
        << W(6)  << right << size
        << W(6)  << right << outl
        << endl;       
  }
  _xyz_est /= double(nsys);
  _xyz_rep /= double(nsys);

  if (_calc < 0) return; 

  // mean blh
  for(auto itSYS = m_blh.begin(); itSYS != m_blh.end(); itSYS++){
     t_gstat3d stat(itSYS->second);
     stat.calc_stat();
     blh = stat.get_mean3d();
     t_gtriple std = stat.get_std3d();
     int size = stat.get_size();
     int outl = stat.get_outl();
     std[0] *= A_WGS;     // convert to meters
     std[1] *= A_WGS;     // convert to meters
     if (std[0] > 9999 || std[1] > 9999 || std[2] > 9999 ) std[0] = std[1] = std[2] = -9999;
       
     string str = "BLH" + t_gsys::gsys2str(itSYS->first);
     _setkey(os, str, '=',tt);
     os << fixed << setprecision(9)
        << W(16) << right << blh[0]*R2D
        << W(16) << right << blh[1]*R2D 
                 << setprecision(4)
        << W(16) << right << blh[2] 
                 << setprecision(1)
        << W(8)  << right << std[0]
        << W(8)  << right << std[1]
                 << setprecision(1)       
        << W(8)  << right << std[2]      
        << W(6)  << right << size
        << W(6)  << right << outl       
        << endl;
  }

  if( _calc < 2 ) return;

  os << endl;
  _setkey(os, "POSGNS", '#', tt);
  os << fixed << setprecision(4)
     << W(16) << right <<   "X [m] "
     << W(16) << right <<   "Y [m] "
     << W(16) << right <<   "Z [m] "
              << setprecision(1)
     << W(16) << right << "B [deg] "
     << W(16) << right << "L [deg] "
     << W(10) << right <<   "H [m] "
     << W(6)  << right << "GDOP"
     << W(6)  << right << "PDOP"
     << W(6)  << right << "HDOP"
     << W(6)  << right << "VDOP"
     << W(16) << right << "REC_CLK[m]"
     << W(5)  << right << "#Sat"
     << W(6) << right << "#Excl"     
     << endl;
   
  double gdop, pdop, hdop, vdop;  
  gdop = pdop = hdop = vdop = -1;

  for(auto itSYS = m_xyz.begin(); itSYS != m_xyz.end(); itSYS++ ){
    vector<t_gtime> vtime = m_time[itSYS->first];
    vector<double> vclk   = m_clk[itSYS->first];
    vector<t_gpair> vnsat = m_nsat[itSYS->first];
    vector<s_dop>   vgdop = m_gdop[itSYS->first];
    vector<t_gtime>::const_iterator itT = vtime.begin();
    vector<double>::const_iterator itCLK = vclk.begin();
    vector<t_gpair>::const_iterator itNSAT = vnsat.begin();
    vector<s_dop>::const_iterator itDOP    = vgdop.begin();

    for(auto itPos = itSYS->second.begin(); itPos != itSYS->second.end(); itPos++ ){
      double clk = *itCLK;
      t_gpair nsat = *itNSAT;
      t_gtime epo = *itT;      
      t_gtriple xyz = *itPos;
      t_gtriple blh;
      xyz2ell(xyz, blh, false);
      string t(epo.str_ymdhms() + " ");
      string str_sys = t_gsys::gsys2str(itSYS->first);      
      gdop = (*itDOP).gdop;
      pdop = (*itDOP).pdop;
      hdop = (*itDOP).hdop;
      vdop = (*itDOP).vdop;            
      
      _setkey(os, "POS" + str_sys, ' ', t);
      os << fixed << setprecision(4)
      << W(16) << right << xyz[0]
      << W(16) << right << xyz[1]
      << W(16) << right << xyz[2]
      << fixed << setprecision(9)
      << W(16) << right << blh[0]*R2D
      << W(16) << right << blh[1]*R2D
      << fixed << setprecision(4)
      << W(10) << right << blh[2]
      << fixed << setprecision(1) 
      << W(6)  << right << gdop
      << W(6)  << right << pdop
      << W(6)  << right << hdop
      << W(6)  << right << vdop
      << W(16) << right << clk      
      << fixed << setprecision(0)
      << W(5)  << right << nsat[0] 
      << W(6)  << right << nsat[1] 
      << endl;
      itT++;
      itCLK++;
      itNSAT++;
      itDOP++;
    }
    os << endl;
  }
}


// Section: obs statistics
// ----------------------------
void t_gxtrqc::_observ(ostringstream& os)
{
  gtrace("t_gxtrqc::_observ");
  if( _stat > 0 ) _section(os, "Observation types", _stat);
  else if( !_stat ) return;
  
  ostringstream  osGNS, osSAT, osTYP, osOBS, osTYH;

  // ORDER IS IMPORTANT ! 
  _list_gsys(osGNS);
  _list_type(osTYP);
//_list_gobs(osOBS);
  _list_typh(osTYH);      

  if( _stat > 0 ) os << osGNS.str() << osSAT.str()
                     << osTYH.str() << osTYP.str();
//                   << osOBS.str();
                    
}

// Section: observations
// ----------------------------
void t_gxtrqc::_nbands(ostringstream& os)
{
  gtrace("t_gxtrqc::_nbands");
  if( _band > 0 ) _section(os, "Band available", _band);
  else if( !_band ) return;

  ostringstream osBND;   
  _list_band(osBND);
  if( _band > 0 ) os << osBND.str();
}

// Section: gaps & pieces
// ----------------------------
void t_gxtrqc::_pieces(ostringstream& os)
{
  gtrace("t_gxtrqc::_pieces");
  if( _gaps > 0 ) _section(os, "Gaps & Pieces", _gaps);
  else if( !_gaps ) return;

  ostringstream osGAP, osPCS, osSMP;
  _list_gap(osGAP);
  _list_pcs(osPCS);
  _list_smp(osSMP);
  if( _gaps > 0 ) os << osGAP.str() << endl
                     << osPCS.str() << endl
                     << osSMP.str();
}

// Section: preprocessing
// ---------------------------
void t_gxtrqc::_prepro(ostringstream& os)
{
  gtrace("t_gxtrqc::_prepro");
  if( _prep > 0 ) _section(os, "Preprocessing results", _prep);
  else if( !_prep ) return;

  t_gpreproc gpreproc(_gobs, _set);  gpreproc.glog(_log);
  gpreproc.setNav(_gnav);
  gpreproc.preprocess(_site, _obs_beg, _obs_end, _smp_est, true, true);

  _m_Slips     = gpreproc.getSlips();
  _m_SlipsGap  = gpreproc.getSlipsGap();
  _m_Breaks    = gpreproc.getClkjump();

  ostringstream osSLP, osJMP, osSUM;

  _list_cslp(osSLP);
  _list_jump(osJMP);
  _list_prep(osSUM);

  if( _prep > 0 ) os << osSUM.str() << endl << osJMP.str() << endl << osSLP.str();
}


// Section: skyplot (elevation and azimuth)
// ---------------------------
void t_gxtrqc::_skyplt(ostringstream& os)
{
  if( _gnav == 0 || _gnav->satellites().size() == 0 ){ _elev = 0; return; }
   
  gtrace("t_gxtrqc::_skyplt");
  if( _elev > 0 ) _section(os, "Elevation & Azimuth", _elev);
  else if( !_elev ) return;
   
  if(!_gnav ) return;
  if( _elev > 0 ) _subsec(os,"LIST_SKY");
  t_gtriple pos;
  ostringstream osAz,osEl;
     
  _legend( osEl, "GNSELE","Mean" );
  _legend( osAz, "GNSAZI","Mean" );
   
  t_gsppflt gspp(_site, _set);

  gspp.setDAT(_gobs, _gnav);
  t_gallprod prod;
  gspp.setOUT(&prod);
  gspp.setOBJ(_gobj);
  gspp.glog(_log);
   
  set<GSYS> sys = _gobs->sys(_site);
   
  for(auto itsys = sys.begin(); itsys != sys.end(); ++itsys)
  {
    GSYS gsys = *itsys;
    string gs = t_gsys::gsys2str( gsys );
     
    // currently only exact epoch implemented!
    t_gtime epo = _obs_sync;

    while(epo < _obs_end || epo == _obs_end )
    {
      t_gtriple crd(0,0,0);
      if (_kinematic || _pos.zero()){
        prod.clear();
        gspp.processBatch(epo,epo);
     
    // use rounded time for clock drift in RINEX since gsppflt reports rounded values
        t_gtime tt(epo);// tt.reset_dsec();
        shared_ptr<t_gprodcrd> p = static_pointer_cast<t_gprodcrd>( prod.get(_site, t_gdata::POS, tt) );
     
        if( p != 0 && ! p->xyz().zero() ) crd = p->xyz();
        else{        epo = epo + _step; continue; }
        
      }else crd = _pos;

      vector<t_gsatdata> satdata = _gobs->obs(_site, epo);

      for(auto itSAT = satdata.begin(); itSAT != satdata.end(); ){
        if( itSAT->addprd_nav(_gnav, false) < 0 ){
//   cout << "erasing due to addprd" << itSAT->sat() << " " << epo.str_hms() << endl;      
          itSAT = satdata.erase(itSAT);
          continue;
        }  
        if( itSAT->cmpVal(crd) < 0 ){
//   cout << "erasing due to cmpVal" << itSAT->sat() << " " << epo.str_hms() << endl;      
          itSAT = satdata.erase(itSAT);
          continue;
        }
        ++itSAT;
      }       
       
      double minEle = 0.0, maxEle = 0.0, aveEle = 0.0;
      double minAzi = 0.0, maxAzi = 0.0, aveAzi = 0.0;
      bool init = false;
      t_gtime t;

      // Min/Max Values
      int count = 0;
      for(auto itSAT = satdata.begin(); itSAT != satdata.end(); ++itSAT)
      {
        if( itSAT->gsys() != gsys ) continue;
        
        string sat = itSAT->sat();
        double ele = itSAT->ele()*rad2deg;
        double azi = itSAT->azi()*rad2deg;
        
        if( !itSAT->health() && _useHealth >= ALL_HEALTH ){
//           cout << sat << " is unhealthy (skip azi/ele) " << _useHealth << " " << itSAT->health() << "\n";
          continue;
        }

        count++; aveEle += ele; aveAzi += azi;
   
        if( ! init ){ minEle = maxEle = ele; minAzi = maxAzi = azi; init = true; }

        if( minEle > ele ) minEle = ele;
        if( maxEle < ele ) maxEle = ele;
        if( minAzi > azi ) minAzi = azi;
        if( maxAzi < azi ) maxAzi = azi;
        t = itSAT->epoch();

        _qcdata->qc_tim.dat[gsys][t][sat].ele = ele;                 // NEW STRUCTURE
        _qcdata->qc_tim.dat[gsys][t][sat].azi = azi;                 // NEW STRUCTURE
        _qcdata->qc_tim.dat[gsys][t][sat].health = itSAT->health();  // NEW STRUCTURE  
      }
      
      _setkey(osEl, gs+"ELE", ' ', epo.str_ymdhms());
      _setkey(osAz, gs+"AZI", ' ', epo.str_ymdhms());

      if( count > 0 ){
        aveEle = aveEle/count;
        aveAzi = aveAzi/count;
        osEl << W(8) << floor(aveEle);
        osAz << W(8) << floor(aveAzi);
      }else{ 
        osEl << W(8) << MISS;
        osAz << W(8) << MISS; 
      }

      // loop SAT (fixed)
      int offset = 0;
      if(     gsys==SBS) offset = SBS_OFFSET;
      else if(gsys==QZS) offset = QZS_OFFSET;
      for(int i = 1+offset; i<=_nsat+offset ; ++i )
      {
        string prn = t_gsys::eval_sat(i, gsys);

        if(    _qcdata->qc_tim.dat.find(gsys)         !=  _qcdata->qc_tim.dat.end()
            && _qcdata->qc_tim.dat[gsys].find(t)      !=  _qcdata->qc_tim.dat[gsys].end()
            && _qcdata->qc_tim.dat[gsys][t].find(prn) !=  _qcdata->qc_tim.dat[gsys][t].end()
//          && ! double_eq( _qcdata->qc_tim.dat[gsys][t][prn].ele, 0.0 )  // show only if real values exists
//          && ! double_eq( _qcdata->qc_tim.dat[gsys][t][prn].azi, 0.0 )  // show only if real values exists
        ){
          if( aveEle != 0 ){ // NAV available (if above is applied, this is irrelevant)
            osEl << W(4) << floor( _qcdata->qc_tim.dat[gsys][t][prn].ele );
            osAz << W(4) << floor( _qcdata->qc_tim.dat[gsys][t][prn].azi );
          }else{             // NAV unavailable (if above is applied, this is irrelevant)
            osEl << W(4) << 0;
            osAz << W(4) << 0;
          }
        }else{
          osAz << W(4) << MISS;
          osEl << W(4) << MISS;
        }
      }

      osEl << endl;
      osAz << endl;
      
      epo = epo + _step;
    }
  }
  
    if( _elev > 0 ) os << osEl.str() << endl << osAz.str();
}


// calculate multipath
// ----------------
void t_gxtrqc::multipath(const string& prn, GSYS gsys, const t_gtime& epo,
             const unsigned int& nepo, t_map_mpobs& m_obs)
{
  gtrace("t_gxtrqc::multipath");
   
  if( gsys > QZS ){ return; } // skip SBAS & regional   

  string gs = t_gsys::gsys2str(gsys);
  string ep = epo.str_ymdhms();
   
  // get all observations for the satellite
  vector<shared_ptr<t_gobsgnss>> vobs_prn = _gobs->obs_prn_pt(_site, prn, epo, epo + 2*nepo*_smp_est );
  if( vobs_prn.size() == 0 ) return;                    // nothing to do

  // L-bands selection for GNSS (priority)
  t_gobs code;
  t_gobs gbL1(TYPE_L, BAND, ATTR);
  t_gobs gbL2(TYPE_L, BAND, ATTR);

  if(!_auto_band){
    gbL1.band( t_gsys::band_priority(gsys,FREQ_1) );           // get first band
    gbL2.band( t_gsys::band_priority(gsys,FREQ_2) );           // get second band    
    
    // NEED ONLY ONCE! ===> SIGNAL/ATTR SELECTION BY QUANTITATIVE PRIORITY 
    // --> can be done at _list_mult + saved in member data!
    gbL1.attr( _gobs->select_attr( _stat_obs, gsys, TYPE_L, gbL1.band() ) );
    gbL2.attr( _gobs->select_attr( _stat_obs, gsys, TYPE_L, gbL2.band() ) );
  }
  
  vector<double> vMP;
  double obs, obs_0, mean, rms, difMP;
  double rmsMP = 1.0; // 1m a priori rms
  double sdif2 = rmsMP*rmsMP;
  unsigned int count;
  bool firstCOD = true;

  // loop GOBS
  t_galloqc::t_map_stt_obs::const_iterator itGOBS;
  for( itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS )
  {
    if( ! gobs_code(itGOBS->first) || vobs_prn.size() == 0 ){ return; } // code only || nothing to do (data removed due to cycle-slips)

    GOBS gobs  = itGOBS->first;
    string go  = gobs2str( gobs );
    code.gobs(gobs);

#ifdef DEBUG
    if( nepo == DAY_EPOCHS )    
    cout << ep
         << "  "    << prn
         << "  go:" << go 
         << "  b1:" << gobsband2str(gbL1.band()) 
         << "  b2:" << gobsband2str(gbL2.band())
         << "  a1:" << gobsattr2str(gbL1.attr()) 
         << "  a2:" << gobsattr2str(gbL2.attr())
         << "  o1:" << gobs2str(    gbL1.gobs())
         << "  o2:" << gobs2str(    gbL2.gobs())
         << "  o1:" << (*vobs_prn.begin())->obs_L(gbL1)
         << "  o2:" << (*vobs_prn.begin())->obs_L(gbL2)
         << "  oC:" << (*vobs_prn.begin())->obs_C(code)
         << "  sz:" << vobs_prn.size()
         << endl; cout.flush();
#endif

    // rmsMP, sdif2 (the same as previous)
    vMP.clear();
    obs_0 = 0.0;
    mean  = 0.0;
    count = 0;

    // estimate constant term (ambig.) for required time interval
    vector<shared_ptr<t_gobsgnss>>::iterator it = vobs_prn.begin();

    for( it = vobs_prn.begin(); it != vobs_prn.end(); ++it ){

      if(_auto_band){    // automatic dual band selection
        set<GOBSBAND> bands = (*it)->band_avail();
        auto itBAND = bands.begin();
        if(bands.size() < 2) continue;
        gbL1.band(*itBAND);
        itBAND++;
        gbL2.band(*itBAND);
    
        gbL1.attr( _gobs->select_attr( _stat_obs, gsys, TYPE_L, gbL1.band() ) );
        gbL2.attr( _gobs->select_attr( _stat_obs, gsys, TYPE_L, gbL2.band() ) );    
      }
          
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
    if(it != vobs_prn.end() && ++it != vobs_prn.end()) vobs_prn.erase(it,vobs_prn.end()); // erase rest of data

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
    for(auto itV = vMP.begin(); itV != vMP.end(); ++itV )
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
    if( rms > 999.0 ){ // cm
      if( _log ) _log->comment(2,"gxtrqc","Unreliable multipath RMS (>10m) "+prn+" "+gobs2str(gobs)+" "+epo.str_ymdhms());
      // max printable value   
//      m_obs[gobs][epo][prn].first  = 999.0;       // DOUBLE (RMS)
//      m_obs[gobs][epo][prn].second = vMP.size();  // COUNT
       
    }else{
      _qcdata->qc_tim.dat[gsys][epo][prn].mpt[gobs].first  = rms;          // NEW STRUCTURE
      _qcdata->qc_tim.dat[gsys][epo][prn].mpt[gobs].second = vMP.size();   // NEW STRUCTURE

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


// Section: multipath
// ----------------
void t_gxtrqc::_mlpath(ostringstream& os)
{
  gtrace("t_gxtrqc::_mlpath");
  if( _mult > 0 ){ _section(os, "Code multipath", _mult);
                    _subsec(os, "LIST_MP");
                    _legend(os, "GNSMxx","mean" ); 
  }else if( !_mult ) return;

  t_galloqc::t_map_stt_sys m_multSys;

  ostringstream osMP, osSUM, *osPT;

  // loop GSYS
  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS )
  {
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);

#ifdef DEBUG                                            // A FULL DAY MP TEST FOR ALL SATELLITES
    int offset = 0;
    if(     gsys==SBS) offset = SBS_OFFSET;
    else if(gsys==QZS) offset = QZS_OFFSET;
    for(int i = 1+offset; i<=_nsat+offset ; ++i ){ 
      string prn = t_gsys::eval_sat(i, gsys);
      multipath(prn, gsys, _obs_beg, DAY_EPOCHS, _m_mp[gsys]); // calc MP for all signals (1 SVN)
    }
#endif

    // loop time
    t_gtime epo = _obs_sync;
    while( epo < _obs_end || epo == _obs_end )
    {
//      if( _step < _mp_nepochs*_smp_est || _step < 30 ) cout << "multipath step prekrocen\n";
//      cout << epo.str_ymdhms("multipath: ") << " " << _step << " " << _mp_nepochs 
//    << " " <<_mp_nepochs*_smp_est << endl;
      int offset = 0;
      if(     gsys==SBS) offset = SBS_OFFSET;
      else if(gsys==QZS) offset = QZS_OFFSET;
      for(int i = 1+offset; i<=_nsat+offset ; ++i )
      {
        string prn = t_gsys::eval_sat(i, gsys);
        multipath(prn, gsys, epo, _mp_nepochs, _m_mp[gsys]); // calc MP for all signals (1 SVN)
      }      
      epo = epo + _step;
    } // while over epochs
  } // loop over systems
   
  for(auto itSys = _m_mp.begin(); itSys != _m_mp.end(); ++itSys )
  {
    GSYS gsys  = itSys->first;
    string gs  = t_gsys::gsys2str( gsys );
     
    for(auto itObs = _m_mp[gsys].begin(); itObs != _m_mp[gsys].end(); ++itObs )
    {
      GOBS gobs  = itObs->first;
      string go  = gobs2str( gobs );
      
      for(auto itEpo = _m_mp[gsys][gobs].begin(); itEpo != _m_mp[gsys][gobs].end(); ++itEpo )
      {
        t_gtime epo = itEpo->first;
        
        // label for EPOCH and all satellites
        string key = gs+"M"+go.substr(1,2);
        if( go[2] == ' ' ) key = gs+"M"+go.substr(0,2); // 2-char obs-type !
   
        osPT = &osMP;
        if( epo == LAST_TIME ){
          osPT = &osSUM;
          _setkey(*osPT, key, '='); // mean
        }else _setkey(*osPT, key, ' ', epo.str_ymdhms());

        // average
        if( // _m_mp[gsys][gobs][epo].size()      > 0 &&
           _m_mp[gsys][gobs][epo].find("")   != _m_mp[gsys][gobs][epo].end() &&
           _m_mp[gsys][gobs][epo][""].second  > 0
           ){
          double mean = _m_mp[gsys][gobs][epo][""].first /
                        _m_mp[gsys][gobs][epo][""].second;
          *osPT << fixed << setprecision(2) << W(8) << mean;
        }else *osPT << fixed << setprecision(2) << W(8) << MISS;
       
        // loop over sats/codes
        int offset = 0;
        if(     gsys==SBS) offset = SBS_OFFSET;
        else if(gsys==QZS) offset = QZS_OFFSET;
        for(int i = 1+offset; i<=_nsat+offset ; ++i )
        {
          string prn = t_gsys::eval_sat(i, gsys);

          if( _m_mp[gsys][gobs][epo].find(prn) != _m_mp[gsys][gobs][epo].end() &&
              _m_mp[gsys][gobs][epo][prn].second > 0 )
          {
            double val = _m_mp[gsys][gobs][epo][prn].first;
            if( epo == LAST_TIME )
              val /= _m_mp[gsys][gobs][epo][prn].second; // mean

            *osPT << setprecision(0) << W(4) << val;
          }else  *osPT << setprecision(0) << W(4) << MISS;
        }
        *osPT << endl;
      }
      osMP << endl;
    }
  }
  
  if( _mult < 1 ) return;

  os << osSUM.str() << endl;

  // extra verbosity 
  if( _mult < 2 ) return;

  os << osMP.str();
  return;
}


// Section: SNR
// ----------------
void t_gxtrqc::_snoise(ostringstream& os)
{
  gtrace("t_gxtrqc::_snoise");

  if( _stnr > 0 ){ _section(os, "Signal to noise ratio", _stnr); 
                    _legend(os, "GNSSxx","mean" );
  }else if( !_stnr ) return;

  ostringstream osv1, osv2;
     
  // loop GSYS
  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS ){
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);

    int offset = 0;
    if(     gsys==SBS) offset = SBS_OFFSET;
    else if(gsys==QZS) offset = QZS_OFFSET;     
     
    for(auto itGOBS = itGSYS->second.begin(); itGOBS != itGSYS->second.end(); ++itGOBS){
      GOBS gobs = itGOBS->first;
      string go = gobs2str(gobs);       
      if(go.compare(0,1,"S") != 0) continue;
      string key = gs+"S"+go.substr(1,2);
      if( go[2] == ' ' ) key = gs+"S"+go.substr(0,2); // 2-char obs-type !
   
      map<string, pair<double, int> > m_obs_stt;
      t_gtime epo = _obs_sync;
      while( epo < _obs_end || epo == _obs_end ){
        _setkey(osv2, key, ' ', epo.str_ymdhms());

        ostringstream osPRN;
        // loop over sats/codes
        double obs = 0.0, sum = 0.0, count = 0.0;
      
        for(int i = 1+offset; i<=_nsat+offset ; ++i ){
          string prn = t_gsys::eval_sat(i, gsys);
          obs = _gobs->find(_site, epo, prn, gobs);
          if(obs > 0){
            _qcdata->qc_tim.dat[gsys][epo][prn].snr[gobs] = obs; // NEW STRUCTURE

            osPRN << W(4) << fixed << setprecision(0) << obs;
            sum += obs; m_obs_stt[prn].first += obs;
            count++;    m_obs_stt[prn].second++;
          }else osPRN << W(4) << MISS;     
        }   
        
        if( count > 0 ) osv2 << W(8) << fixed << setprecision(2) << sum/count;
        else            osv2 << W(8) << fixed << setprecision(2) << MISS;
      
        osv2 << osPRN.str() << endl;
        osPRN.str("");
        osPRN.clear();
        epo = epo + _step;
      }
   
      _setkey(osv1, key, '=');
      ostringstream osPRN;
      double sum_prn  = 0.0;
      for( int i = 1+offset; i<=_nsat+offset ; ++i ){
        string prn = t_gsys::eval_sat(i, gsys);
        if( m_obs_stt.find(prn) != m_obs_stt.end() ){
          if( m_obs_stt[prn].second > 0 ){
            double mean_sys = m_obs_stt[prn].first / m_obs_stt[prn].second;
            sum_prn += mean_sys;
            osPRN << fixed << setprecision(0) << W(4) << mean_sys;
          }else osPRN << fixed << setprecision(0) << W(4) << MISS; 
          
        }else osPRN << W(4) << MISS;
      }
   
      if( m_obs_stt.size() > 0 ) osv1 << fixed << setprecision(2) << W(8) << sum_prn / m_obs_stt.size();
      else osv1 << fixed << setprecision(2) << W(8) << MISS;
   
      osv1 << osPRN.str() << endl;
      osv2 << endl;
    }               

  } // loop over systems   
     
  if(_stnr <  0) return;
  if(_stnr <  2) os << osv1.str();
  if(_stnr >= 2) os << osv1.str() << endl << osv2.str();
  
  return;
}     

// Section: Navigation indicators
// ----------------
void t_gxtrqc::_svinfo(ostringstream& os)
{
  gtrace("t_gxtrqc::_svinfo");

  if( _sinf > 0 ){ _section(os, "Satellite information", _sinf); 
                    _legend(os, "GNSIxx");
  }else if( !_sinf ) return;

  set<GSYS> set_sys = _gnav->systems();

  ostringstream os_v2;

  // extend window for NAV/SP3
  t_gtime beg = _beg - 7200;
  t_gtime end = _end + 7200;

  // loop GSYS
  for( auto itGSYS = set_sys.begin(); itGSYS != set_sys.end(); ++itGSYS )
  {
    string gs = t_gsys::gsys2str(*itGSYS);

    set<t_gtime> vec_epo = _gnav->vec_epo(*itGSYS); 
    string key = gs + "HLT";
        
    for(auto itEPO = vec_epo.begin(); itEPO != vec_epo.end(); itEPO++)
    {
      if( *itEPO < _beg || *itEPO > _end ) continue;
      
      _setkey(os, key, ' ', itEPO->str_ymdhms());  os << W(8) << "";
      
      if(_sinf >= 2){
        string key2 = gs + "RED";
        _setkey(os_v2, key2, ' ', itEPO->str_ymdhms());
        os_v2 << W(8) << "";
      }
      
      int offset = 0;
      if(     *itGSYS == SBS) offset = SBS_OFFSET;
      else if(*itGSYS == QZS) offset = QZS_OFFSET;

      for (int i = 1 + offset; i <= _nsat + offset; ++i) {
        string prn = t_gsys::eval_sat(i, *itGSYS);
        vector<shared_ptr<t_geph>> v_eph = _gnav->find_mult(prn, *itEPO);   // used find_mult because health status is not used
        if(*(v_eph.begin()) != nullptr) os << W(4) << (*v_eph.begin())->healthy();
        else                            os << W(4) << MISS;

        if(_sinf >= 2){
          if( v_eph.size() == 1 && *(v_eph.begin()) == nullptr) os_v2 << W(4) << MISS;
          else os_v2 << W(4) << v_eph.size();
        }
      }
      
      os << endl; 
      if( _sinf >=2 ) os_v2 << endl;
    }
    os << endl;
    if( _sinf >=2 ) os_v2 << endl;    
  }

  if(_sinf >= 2) os << os_v2.str();
}


// XML FILE META DATA store
// -----------------------------
void t_gxtrqc::_xml_meta()
{
  gtrace("t_gxtrqc::_xml_meta");
   
  t_gtime tt(t_gtime::UTC);        // UTC!
  string vers = "1.00";
  string syst = "";

  double     sampling = dynamic_cast<t_gsetgen*>(_set)->sampling();
//set<string> systems = dynamic_cast<t_gsetgen*>(_set)->sys();
  set<GSYS>   systems = _gobs->sys(_site);

  string appl = _pgm.substr(0,1+_pgm.find(']'));
   
  _default_node(_XFMT, "created", tt.str("%Y-%m-%dT%H:%M:%S").c_str(), true );
  _default_node(_XFMT, "program", trim(appl).c_str(), true );
   
  xml_node XSET = _XFMT.append_child("settings");
  
  _default_node(XSET, "time_beg", trim(_set_beg.str("%Y-%m-%dT%H:%M:%S")).c_str(), true); // QC data interval
  _default_node(XSET, "time_end", trim(_set_end.str("%Y-%m-%dT%H:%M:%S")).c_str(), true); // QC data interval
  _default_node(XSET, "data_int", trim(dbl2str(sampling,3)).c_str(), true);
  _default_node(XSET, "elev_min", trim(dbl2str(_cut_ele,2)).c_str(), true);
   
  for(auto itSYS = systems.begin(); itSYS != systems.end(); ++itSYS ){
    xml_node XSYS = XSET.append_child("system");
    _default_attr(XSYS, "type", t_gsys::gsys2str(*itSYS), true);
//  _default_attr(XSYS, "type", *itSYS, true);
  }   
}


// XML HEADER extract
// -----------------------------
void t_gxtrqc::_xml_head()
{
  gtrace("t_gxtrqc::_xml_head");
      
  int cnt = 1;
  for(auto itHDR = _m_rnxHdr.begin(); itHDR != _m_rnxHdr.end(); ++itHDR, ++cnt )
  {
    t_rnxhdr rnxhdr = itHDR->second;

    // to enable multi-header output
    xml_node _XHDR = _ROOT.append_child("head");
    _default_attr(_XHDR, "file_numb", int2str(cnt), true );

    string name = base_name( rnxhdr.path() );    
    t_gfile gf; gf.path(rnxhdr.path());
    string md5s   = gf.md5sum();
    string format = "RINEX" + trim(rnxhdr.rnxver()); substitute(format," ","");
    string radome;
    // safety reasons -> may happen anttype is shorter!
    if( rnxhdr.anttype().length() > 16 ) radome = rnxhdr.anttype().substr(16,4); 

    _default_node(_XHDR, "file_name",        trim(name).c_str(), true );
    _default_node(_XHDR, "file_md5sum",      trim(md5s).c_str(), true );
    _default_node(_XHDR, "file_format",      trim(format).c_str(), true);
    _default_node(_XHDR, "site_id",          trim(_site).c_str(), true);
    _default_node(_XHDR, "marker_name",      trim(rnxhdr.markname()).c_str(), true);
    _default_node(_XHDR, "marker_numb",      trim(rnxhdr.marknumb()).c_str(), true);
    _default_node(_XHDR, "receiver_type",    trim(rnxhdr.rectype()).c_str(), true);
    _default_node(_XHDR, "receiver_numb",    trim(rnxhdr.recnumb()).c_str(), true);
    _default_node(_XHDR, "antenna_type",     trim(rnxhdr.anttype().substr(0,15)).c_str(), true);
    _default_node(_XHDR, "antenna_dome",     trim(radome).c_str(), true); // DOME - may happen is not available!
    _default_node(_XHDR, "antenna_numb",     trim(rnxhdr.antnumb()).c_str(), true);
    _default_node(_XHDR, "software",         trim(rnxhdr.program()).c_str(), true);
    _default_node(_XHDR, "data_int",         trim(dbl2str(rnxhdr.interval())).c_str(), true);
   
    t_gtriple xyz = rnxhdr.aprxyz();
    t_gtriple neu = rnxhdr.antneu();

    xml_node POS = _XHDR.append_child("position");
    xml_node CRD =   POS.append_child("coordinate");
    xml_node ECC =   POS.append_child("eccentricity");
    _default_attr(POS, "type", string("GNS"), true);
    _default_attr(POS, "source", string("rinex_header"), true);
    string ecclabel[] = { "N", "E", "U" };

    xml_node XYZ = CRD.append_child("gml:Point");
    _default_attr(XYZ, "gml:id",  string("rinex_header"), true);
    _default_node(XYZ, "gml:pos", trim(dbl2str(xyz[0])+dbl2str(xyz[1])+dbl2str(xyz[2])).c_str(), true);

    for( int i = 0; i < 3; ++i ){
      xml_node XXX = ECC.append_child("axis");
      XXX.append_child(pugi::node_pcdata).set_value(trim(dbl2str(neu[i],4)).c_str());
      _default_attr(XXX, "name", ecclabel[i], true);
    }

    t_rnxhdr::t_obstypes typehdr = rnxhdr.mapobs();
    for(auto itG = typehdr.begin(); itG != typehdr.end(); ++itG ){
      xml_node SYS = _XHDR.append_child("system");
      _default_attr(SYS, "type", trim(t_gsys::char2str(itG->first[0])), true);

      for(auto itB = itG->second.begin(); itB != itG->second.end(); ++itB ){
        xml_node OBS = SYS.append_child("obs");
        _default_attr(OBS, "type", trim(gobs2str(itB->first)), true);
      }
    }
  }   
}


// XML DATA extract
// -----------------------------
void t_gxtrqc::_xml_data()
{
  gtrace("t_gxtrqc::_xml_data");

  // actual observation time !
  int xbeg = 0, xend = 0, xsmp = 0, xsys = 0;
  t_gtime beg = _obs_beg; // from DATA (already filtered)
  t_gtime end = _obs_end; // from DATA (already filtered)

  // collect filtered data
  for(auto itHDR = _m_rnxHdr.begin(); itHDR != _m_rnxHdr.end(); ++itHDR )
  {   
    t_rnxhdr rnxhdr = itHDR->second;
    t_gallobs::t_xfilter xflt = _gobs->xdata(_site, rnxhdr.path());

    xbeg += xflt.xdat[t_gallobs::XDATA_BEG];
    xend += xflt.xdat[t_gallobs::XDATA_END];
    xsmp += xflt.xdat[t_gallobs::XDATA_SMP];
    xsys += xflt.xdat[t_gallobs::XDATA_SYS];

    if( xflt.beg < beg ) beg = xflt.beg; // use unfiltered BEG epoch
    if( xflt.end > end ) end = xflt.end; // use unfiltered END epoch
  }

  _default_node(_XDAT, "time_beg", trim(beg.str("%Y-%m-%dT%H:%M:%S")).c_str(), true); // actual data interval UPDATE from DECODER?
  _default_node(_XDAT, "time_end", trim(end.str("%Y-%m-%dT%H:%M:%S")).c_str(), true); // actual data interval UPDATE from DECODER?
  _default_node(_XDAT, "data_int", trim(dbl2str(_smp_est,2)).c_str(),          true); // estimated sampling
  _default_node(_XDAT, "numb_epo", trim(dbl2str(_hours,0)).c_str(),            true); // # epochs
  _default_node(_XDAT, "numb_gap", trim(int2str(_ngap)).c_str(),               true); // # gaps
   
  xml_node TOT = _XDAT.append_child("total");
  _default_node(TOT, "expt",     trim(int2str(_m_obsExp[GNS][X])).c_str(),        true); // # expected obs
  _default_node(TOT, "have",     trim(int2str(_m_obsHav[GNS][X])).c_str(),        true); // # existing obs
  _default_node(TOT, "expt_usr", trim(int2str(_m_cutExp[GNS][X])).c_str(),        true); // # expected obs (user)
  _default_node(TOT, "have_usr", trim(int2str(_m_cutHav[GNS][X])).c_str(),        true); // # existing obs (user)
  _default_node(TOT, "elev_min", trim(dbl2str(_min_ele,2)).c_str(),               true); // minimum elevation
  _default_node(TOT, "cyc_slps", trim(int2str(_m_totSlp[GNS][X])).c_str(),        true); // # cycle slips
  _default_node(TOT, "clk_jmps", trim(int2str(_m_Breaks.size()+_ClkSync)).c_str(),true); // # clock jumps (both types)

  xml_node EXL = _XDAT.append_child("excluded");
  _default_node(EXL, "xbeg", trim(int2str(xbeg)).c_str(), true);
  _default_node(EXL, "xend", trim(int2str(xend)).c_str(), true);
  _default_node(EXL, "xint", trim(int2str(xsmp)).c_str(), true);
  _default_node(EXL, "xsys", trim(int2str(xsys)).c_str(), true);

  // GNSS-specific output
  t_gtime epo_mxsum(LAST_TIME);
  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS ){
    GSYS gsys = itGSYS->first;

    int nsat = _get_nsat(gsys);

    string xele = int2str(_m_woElev[gsys][X]);

    xml_node GNODE = _XDAT.append_child("system");
    _default_attr(GNODE, "type", trim(t_gsys::gsys2str(gsys)),             true); // DO NOT TRIM!     
    _default_attr(GNODE, "nsat", trim(int2str(nsat)),                      true);
    _default_attr(GNODE, "xele", trim(xele),                               true);

    xml_node EPO = GNODE.append_child("epo");
    _default_node(EPO, "expt",  trim(int2str(_nEpoExpect[gsys])).c_str(),  true); // expected number of epochs
    _default_node(EPO, "have",  trim(int2str(_nEpoExists[gsys])).c_str(),  true); // existing number of epochs
    _default_node(EPO, "dual",  trim(int2str(_m_FineEpo[gsys])).c_str(),   true); // dual-frq number of epochs
     
    xml_node AMB = GNODE.append_child("amb");
    _default_node(AMB, "nepo",  trim(int2str(_m_totEpo[gsys][X])).c_str(), true); // #interruptions due to epochs
    _default_node(AMB, "nsat",  trim(int2str(_m_totSat[gsys][X])).c_str(), true); // #interruptions due to satellites
    _default_node(AMB, "nsig",  trim(int2str(_m_totSig[gsys][X])).c_str(), true); // #interruptions due to signals
    _default_node(AMB, "nslp",  trim(int2str(_m_totSlp[gsys][X])).c_str(), true); // #interruptions due to cycle slips
          
    xml_node COD = GNODE.append_child("bnd");
    _default_attr(COD, "type", string("code"), true);
    _default_node(COD, "xepo", trim(int2str(_m_UnusEpo[gsys].first)).c_str(),  true);
    _default_node(COD, "xsat", trim(int2str(_m_UnusSat[gsys].first)).c_str(),  true);

    xml_node PHA = GNODE.append_child("bnd");
    _default_attr(PHA, "type", string("phase"), true);
    _default_node(PHA, "xepo", trim(int2str(_m_UnusEpo[gsys].second)).c_str(), true);
    _default_node(PHA, "xsat", trim(int2str(_m_UnusSat[gsys].second)).c_str(), true);
     
    ostringstream sat_list, obs_list, obs_sats, obs_have, obs_expt, cod_mpth, pha_cslp;

    // assume all signals has the full set of satellites (although sometimes zero observation)
    auto itGOBS = itGSYS->second.begin();

    int offset = 0;
    if(     gsys==SBS) offset = SBS_OFFSET;
    else if(gsys==QZS) offset = QZS_OFFSET;
    for(int i = 1+offset; i<=_nsat+offset ; ++i ){
      string prn = t_gsys::eval_sat( i, gsys );
      if( itGOBS->second.find(prn) != itGOBS->second.end() ) sat_list << W(4) << prn;
    }

    for(auto itGOBS = itGSYS->second.begin(); itGOBS != itGSYS->second.end(); ++itGOBS )
    {
      GOBS gobs = itGOBS->first;

#ifdef DEBUG
      obs_list << " " << gobs2str(gobs);
      obs_sats << " " << _m_satHav[gsys][gobs];
      obs_expt << " " << _m_obsExp[gsys][X];
#endif       
      string expZ(XML_NAN);
      string expU(XML_NAN);
      if( _m_okElev[gsys][gobs] > 0 ){
        if( _m_woElev[gsys][gobs] < _m_obsHav[gsys][gobs] ) expZ = int2str(_m_obsExp[gsys][gobs]); // [gsys][X] for version < v2.0.0:
        if( _m_woElev[gsys][gobs] < _m_cutHav[gsys][gobs] ) expU = int2str(_m_cutExp[gsys][gobs]); // [gsys][X] for version < v2.0.0:
      }
      
      xml_node OBS = GNODE.append_child("obs");
      _default_attr(OBS, "type",     trim(gobs2str(gobs)), true);
      _default_node(OBS, "nsat",     trim(int2str(_m_satHav[gsys][gobs])).c_str(), true);
      _default_node(OBS, "expt",     trim(expZ).c_str(),                           true);
      _default_node(OBS, "have",     trim(int2str(_m_obsHav[gsys][gobs])).c_str(), true);
      _default_node(OBS, "expt_usr", trim(expU).c_str(),                           true);
      _default_node(OBS, "have_usr", trim(int2str(_m_cutHav[gsys][gobs])).c_str(), true);
       
      if( gobs_code(gobs) ){
        string mpth(XML_NAN);
        if( _m_mp[gsys].find(gobs) != _m_mp[gsys].end() ){
          mpth = dbl2str(  _m_mp[gsys][gobs][epo_mxsum][""].first
                         / _m_mp[gsys][gobs][epo_mxsum][""].second, 1);
        }
        _default_node(OBS, "mpth", trim(mpth).c_str(), true);
      }

      if( gobs_phase(gobs) ){
        _default_node(OBS, "slps", trim(int2str(_m_totSlp[gsys][gobs])).c_str(), true);
      }
    }
  }
  
  // the POSITION sequency should be after all SYSTEMs (above)
  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS )
  {
    GSYS gsys = itGSYS->first;

    if( double_eq(_m_xyz_est[gsys][0],0.0) && 
        double_eq(_m_xyz_est[gsys][1],0.0) && 
        double_eq(_m_xyz_est[gsys][2],0.0)  ) continue;
       
    xml_node POS = _XDAT.append_child("position");
    xml_node CRD =   POS.append_child("coordinate");
    xml_node DEV =   POS.append_child("deviation");
    _default_attr(POS, "type", trim(t_gsys::gsys2str(gsys)), true);
    _default_attr(POS, "source", string("analysis"), true);
    string crdlabel[] = { "X", "Y", "Z" };

    xml_node XYZ = CRD.append_child("gml:Point");
    _default_attr(XYZ, "gml:id",  string(t_gsys::gsys2str(gsys)+"_analysis"), true);
    _default_node(XYZ, "gml:pos", trim(dbl2str(_m_xyz_est[gsys][0])
                  +dbl2str(_m_xyz_est[gsys][1])
                  +dbl2str(_m_xyz_est[gsys][2])).c_str(), true);
    for( int i = 0; i < 3; ++i ){
     xml_node XXX = DEV.append_child("axis");
     XXX.append_child(pugi::node_pcdata).set_value(trim(dbl2str(_m_xyz_rep[gsys][i])).c_str());
      _default_attr(XXX, "name", crdlabel[i], true);
    }
  }   
}


// Section: XML Output 2 skyplot (qc observations)
// ---------------------------
void t_gxtrqc::_xml_full()
{
  if( !_xml_ext ) return;
  
  _XEXT = _default_node(_ROOT, "data_epo");
  set<GSYS> sys = _gobs->sys(_site);
  
  /// TODO: If _elev switched off and the _band number is switched on, than it is still blank data.
  /// TODO: check the gap
  if( _elev < 1 ) return;

  for(auto itSYS  = _qcdata->qc_tim.dat.begin();
           itSYS != _qcdata->qc_tim.dat.end(); ++itSYS)
  {
    GSYS gsys = itSYS->first;
    xml_node sysNode = _XEXT.append_child("system");  _default_attr(sysNode, "type", t_gsys::gsys2str( gsys ), true);
    
    for(auto itEPO  = itSYS->second.begin();
             itEPO != itSYS->second.end(); ++itEPO)
    {
      xml_node epochNode = sysNode.append_child("epoch");
      t_gtime tt = itEPO->first;
      string  ts = tt.str_ymdhms(); // timeData
              ts.replace(10,1,"T");

      int nsat = 0;
      for(auto itPRN = _sat_view.begin(); itPRN != _sat_view.end(); ++itPRN )
      {
        string prn = itPRN->first;
//      if( _stat_obs.find(gsys) == _stat_obs.end() ) continue;
        if( t_gsys::gsys2char(gsys) != prn[0] ) continue;  // DATA MISSING
//      if( !_qcdata->qc_tim.dat[gsys][tt][prn].health && _useHealth >= ALL_HEALTH ){ continue; } // UNHEALTHY

        for(auto itBEG = _sat_view[prn].begin(); itBEG != _sat_view[prn].end(); ++itBEG ){
          if( itBEG->first > tt || itBEG->second < tt ) continue;
          nsat++;          
        }
      }

//      int sec = (itTime.gwk()*604800) + (itTime.dow()*86400)+ (itTime.sod()+round(itTime.dsec()));

      _default_attr(epochNode, "time", ts,   true);
      _default_attr(epochNode, "nsat", nsat, true);

      for(auto itSAT  = itEPO->second.begin();
               itSAT != itEPO->second.end(); ++itSAT )
      {
        string sat = itSAT->first;

        if( _stat_obs.find(gsys) == _stat_obs.end() ) continue; // DATA MISSING
        if( !_qcdata->qc_tim.dat[gsys][tt][sat].health && _useHealth >= ALL_HEALTH ){ continue; } // UNHEALTHY

        xml_node obsNode = epochNode.append_child("sat");
		    _default_attr(obsNode, "id", sat, true);

        if( _elev  ){
//          _default_attr(obsNode, "ele", itSAT->second.ele, true);
//          _default_attr(obsNode, "azi", itSAT->second.azi, true);
          _default_node(obsNode, "ele", trim(dbl2str(itSAT->second.ele)).c_str(), true);
          _default_node(obsNode, "azi", trim(dbl2str(itSAT->second.azi)).c_str(), true);
        }
        
        if( _band ){
//          _default_attr(obsNode, "cbn", itSAT->second.cbn, true);
//          _default_attr(obsNode, "lbn", itSAT->second.lbn, true);
          _default_node(obsNode, "cbn", trim(int2str(itSAT->second.cbn)).c_str(), true);
          _default_node(obsNode, "lbn", trim(int2str(itSAT->second.lbn)).c_str(), true);
        }

        // Multipath
        if( _mult ){
//        xml_node mpt_node = obsNode.append_child("multipath");
          for(auto itOBS  = itSAT->second.mpt.begin();
                   itOBS != itSAT->second.mpt.end(); ++itOBS )
          {
//            double val = itOBS->second.first;
//
//            if( val > 0.0 ){
              xml_node qc_obs = obsNode.append_child("obs");              
//              xml_node qc_obs = mpt_node.append_child("obs");
//              xml_node qc_obs = _default_node(mpt_node, "obs", trim(dbl2str(val)).c_str(), true);

              _default_attr(qc_obs, "type",  trim(gobs2str(itOBS->first)),  true);
              _default_node(qc_obs, "mpth",  trim(dbl2str(itOBS->second.first,2)).c_str(), true);
//            _default_attr(qc_obs, "value", trim(dbl2str(val,2)), true);
//            }
          }
        }
        // Signal-to-noise
        if( _stnr ){
//          xml_node snr_node = obsNode.append_child("signalnoise");
          for(auto itOBS  = itSAT->second.snr.begin();
                   itOBS != itSAT->second.snr.end(); ++itOBS )
          {
//            double val = itOBS->second;
//
//            if( val > 0.0 ){
              xml_node qc_obs = obsNode.append_child("obs");
//              xml_node qc_obs = snr_node.append_child("obs");

              _default_attr(qc_obs, "type",  trim(gobs2str(itOBS->first)),  true);
              _default_node(qc_obs, "mpth",  trim(dbl2str(itOBS->second,2)).c_str(), true);
//              _default_attr(qc_obs, "value", trim(dbl2str(val,2)), true);
//            }
          }
        }
      }
    }
  }
}



// XML NAVIGATION DATA store
// -----------------------------
void t_gxtrqc::_xml_navi()
{
  gtrace("t_gxtrqc::_xml_navi");

  for(auto itGSYS = _stat_obs.begin(); itGSYS != _stat_obs.end(); ++itGSYS)
  {
    GSYS gsys = itGSYS->first;
    string gs = t_gsys::gsys2str(gsys);
     
    xml_node SYS = _XNAV.append_child( "system" );
    _default_attr(SYS, "type", trim(gs), true );
     
    int nsat = _gnav->nsat(gsys);
    int have = _gnav->have(gsys, _beg, _end);

    _default_attr(SYS, "nsat", nsat,                  true ); // mandatory
    _default_node(SYS, "have", int2str(have).c_str(), true ); // mandatory
  }
}


// Set log
// -----------------------------
void t_gxtrqc::_setOut()
{
  gtrace("t_gxtrqc::_setOut");

  string tmp;
  int toff;
  t_gtime::t_tsys tsys;
  t_gxtr::_setOut();

  if( _xtr ){ if( _xtr->is_open() ) _xtr->close(); delete _xtr; _xtr = 0; }  
  tmp  = dynamic_cast<t_gsetout*>(_set)->outputs("xtr");
  toff = dynamic_cast<t_gsetout*>(_set)->file_toff("xtr");
  tsys = dynamic_cast<t_gsetout*>(_set)->file_tsys("xtr");
  if( !tmp.empty() ){
    substitute(tmp,"$(rec)",_site, false);
    
    _xtr = new t_glog;
    _xtr->tsys(tsys);
    _xtr->toff(toff);
    _xtr->mask(tmp);
    _xtr->append(false);
    _xtr->verb(dynamic_cast<t_gsetout*>(_set)->verb());
  }

  tmp  = dynamic_cast<t_gsetout*>(_set)->outputs("xml");
  toff = dynamic_cast<t_gsetout*>(_set)->file_toff("xml");
  tsys = dynamic_cast<t_gsetout*>(_set)->file_tsys("xml");
  if( !tmp.empty() ){
    substitute(tmp,"$(rec)",_site, false);
    substitute(tmp,GFILE_PREFIX,"");

    t_gtime file_tm = t_gtime::current_time(tsys);
            file_tm.add_secs(toff*60.0);
    _name = file_tm.str(tmp);

    _xml = true;
    xml_node root = _doc;
     
    _ROOT = this->_default_node( root, _root.c_str() );
    _XFMT = this->_default_node(_ROOT, "meta");
//  _XHDR = this->_default_node(_ROOT, "head"); // to enable multi-header output
    _XNAV = this->_default_node(_ROOT, "navi");
    _XDAT = this->_default_node(_ROOT, "data");

    // GML ID
    _default_attr(_ROOT, "gml:id", _site+"_"+_obs_beg.str("%Y%m%dT%H%M%S")+"_"
                                            +_obs_end.str("%Y%m%dT%H%M%S"), true );
    _default_attr(_ROOT, "xmlns",              string("https://software.pecny.cz/anubis"), true );
    _default_attr(_ROOT, "xmlns:gml",          string("http://www.opengis.net/gml/3.2"), true );
    _default_attr(_ROOT, "xmlns:xsi",          string("https://www.w3.org/2001/XMLSchema/instance"), true );
//  _default_attr(_ROOT, "xsi:schemaLocation", string("https://software.pecny.cz/anubis qc_gnss.xsd"), true );
    _default_attr(_ROOT, "schemaVersion",      string("1.1"), true);        
  }
/*
  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("xqc");
  if( !tmp.empty() ){
    substitute(tmp,"$(rec)",_site, false);
    substitute(tmp,GFILE_PREFIX,"");
    _xml2 = new t_gxml;                    //////// !!!!!!!!!!!
//  _name = tmp;                           //////// !!!!!!!!!!!

    xml_node root = _doc;                  //////// !!!!!!!!!!!
     
    _ROOT2 = this->_default_node( root, _root.c_str() ); //////// !!!!!!!!!!!
    _XFMT2 = this->_default_node(_ROOT2, "meta");
//  _XHDR2 = this->_default_node(_ROOT2, "head"); // to enable multi-header output
    _XNAV2 = this->_default_node(_ROOT2, "navi");
    _XDAT2 = this->_default_node(_ROOT2, "data");
    _XEXT2 = this->_default_node(_ROOT2, "full");

    // GML ID
    _default_attr(_ROOT, "gml:id", _site+"_"+_obs_beg.str("%Y%m%dT%H%M%S")+"_"
                                            +_obs_end.str("%Y%m%dT%H%M%S"), true );
    _default_attr(_ROOT, "xmlns",              string("https://software.pecny.cz/anubis"), true );
    _default_attr(_ROOT, "xmlns:gml",          string("http://www.opengis.net/gml/3.2"), true );
    _default_attr(_ROOT, "xmlns:xsi",          string("https://www.w3.org/2001/XMLSchema/instance"), true );
//  _default_attr(_ROOT, "xsi:schemaLocation", string("https://software.pecny.cz/anubis qc_gnss.xsd"), true );
    _default_attr(_ROOT, "schemaVersion",      string("1.1"), true);
  }
*/
}

// Estimation of position based on code data and few first epochs
// -------------------------------
t_gtriple t_gxtrqc::_xyz_approx()
{
   t_gtriple crd(0,0,0);   
   
   if(_gnav == 0) return crd;

   t_gsppflt gspp(_site, _set);

   gspp.setDAT(_gobs, _gnav);
   t_gallprod prod;
   gspp.setOUT(&prod);
   gspp.setOBJ(_gobj);
   gspp.glog(_log);

   t_gtime epo = _obs_sync;
   
   while(epo < _obs_end || epo == _obs_end )
   {  
      prod.clear();

      gspp.processBatch(epo,epo);

      shared_ptr<t_gprodcrd> p = static_pointer_cast<t_gprodcrd>( prod.get(_site, t_gdata::POS, epo) );

      if( p && ! p->xyz().zero() ){
        crd = p->xyz(); break;
      }else{ epo = epo + 30; continue; }
   }
   
   return crd;
}        
     
// Clear all containers for summary
// -----------------------------
void t_gxtrqc::_clear()
{
  gtrace("t_gxtrqc::_clear");
   
   _sat_view.clear();
   
   _m_totAll.clear();
   _m_totSlp.clear();
   _m_totEpo.clear();
   _m_totSat.clear();
   _m_totSig.clear();
   
   _m_obsHav.clear();
   _m_obsExp.clear();
   _m_cutHav.clear();
   _m_cutExp.clear();
   
   _m_okElev.clear();
   _m_woElev.clear();

   _m_FineEpo.clear();
   _m_UnusEpo.clear(); 
   _m_UnusSat.clear();
   
   _m_nSatEpo.clear();
   _m_Breaks.clear();
   _m_Slips.clear();
   _m_SlipsGap.clear();

   _stat_obs.clear();
   _stat_ele.clear();
   _stat_smp.clear();
   
   _m_mp.clear();
   
}
   

// Clear all containers for summary
// -----------------------------
int t_gxtrqc::_get_nsat(GSYS gsys)
{
  gtrace("t_gxtrqc::_get_nsat");
   
  int numbSat = 0;

  if( _stat_obs.find(gsys) == _stat_obs.end() ) return numbSat;

  for(auto itGOBS = _stat_obs[gsys].begin(); itGOBS != _stat_obs[gsys].end(); ++itGOBS )
    if( numbSat < _m_satHav[gsys][itGOBS->first] )
        numbSat = _m_satHav[gsys][itGOBS->first];

  return numbSat;
}


// Get header information
// ----------------------
int t_gxtrqc::_get_rnxhdr()
{
  shared_ptr<t_gobj> obj = _gobj->obj(_site);

  t_rnxhdr rnxhdr;

  if( obj && obj->isrec() ){
    _m_rnxHdr = dynamic_pointer_cast<t_grec>(obj)->gethdr();

    if( _m_rnxHdr.size() == 0 ){
      if( _log ) _log->comment(0,"gxtrqc","Warning: Receiver " + obj->id() + " has no header information!");
      else               cerr << "gxtrqc - Warning: Receiver "<< obj->id()<< " has no header information!.\n";
      return -1;
    }

  }else{
    if( _log ) _log->comment(0,"gxtrqc","Warning: no object found.");
    else               cerr << "gxtrqc - Warning: no object found.\n";
    return -1;
  }
  return 0;
}

} // namespace
