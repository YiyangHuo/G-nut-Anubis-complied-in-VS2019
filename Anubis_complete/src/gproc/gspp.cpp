
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
#include "gproc/gspp.h"
#include "gmodels/gsppmodel.h"

namespace gnut {

// Constructor
t_gspp::t_gspp(string mark, t_gsetbase* set)
:
  _grec(nullobj),
  _gallobj(0),
  _gallbias(0),   
  _weight(SINEL),
  _observ(IONO_FREE),
  _valid_crd_xml(false),
  _valid_ztd_xml(false),
  _site(mark),
  _set(set),
  _log(0),
  _glog(0),
  _res(0),  
  _gobs(0),
  _gnav(0),
  _gmet(0),
  _gion(0),
  _allprod(0),
  _phase(false),
  _gnss(GNS),
  _initialized(false),
  _use_ecl(false)       
{
   _setLog();       // set log  (before using it later on!)

   _crd_begStat = _crd_endStat = FIRST_TIME;
   _ztd_begStat = _ztd_endStat = FIRST_TIME;

   _valid_ztd_xml = false;

   _get_settings();
   
   _gModel   = new t_gsppmodel(_site,_log,_set);   
}

// Destructor
t_gspp::~t_gspp()
{
   if( _gModel   ) delete _gModel;
   
   if( _log ){ if( _log->is_open() ){ _log->close(); }; delete _log; }
   if( _res ){ if( _res->is_open() ){ _res->close(); }; delete _res; }   
}

// Set obs, nav
// --------------------------------
void t_gspp::setDAT(t_gallobs* gobs, t_gallnav* gnav)
{      
  _gobs = gobs;
  _gnav = gnav;
   
  //this->fixObsTypes();

}

// set output for products
// ----------------------------
void t_gspp::setOUT(t_gallprod* products)
{
	this->_setOut();
	_allprod = products;

}
   
// Set OBJ
// --------------------------------
void t_gspp::setOBJ(t_gallobj* gallobj)
{      
   _gallobj = gallobj;
   
   if(_gallobj)	_grec = _gallobj->obj(_site);
}   
 
// set DCB products
// ----------------------------
void t_gspp::setDCB(t_gallbias* bias)
{
  _gallbias = bias;
  _gModel->setBIAS(bias);
}
   
// set FCB products
// ----------------------------
void t_gspp::setFCB(t_gallbias* bias)
{
	_gallfcb = bias;
}

// Tropo is used or not
// -------------------------
void t_gspp::tropo(bool tropo)
{
   this->_tropo_est = tropo;
}


// Tropo slants are provided
// -------------------------
void t_gspp::tropo_slant(bool slant)
{
   this->_tropo_slant = slant;
}


// Tropo is used or not
// -------------------------
void t_gspp::phase(bool phase)
{
   this->_phase = phase;
}   

// set used GNSS
// --------------------------
void t_gspp::setgnss(GSYS sys)
{
   this->_success = true;
   this->_gnss = sys;
}

// set map of signals used in processing
void t_gspp::fixObsTypes()
{
  // find out frequency of individual signals from gallnav
  t_gallobs::t_map_frq mfrq = _gobs->frqobs(_site);

  // find out most frequent signals
  for (map<string, map<GOBSBAND, map<GOBS, int>>>::iterator itPRN = mfrq.begin(); itPRN != mfrq.end(); itPRN++) {
    string prn = itPRN->first;

    for (map<GOBSBAND, map<GOBS, int>>::iterator itBAND = itPRN->second.begin(); itBAND != itPRN->second.end(); itBAND++) {
      GOBSBAND band = itBAND->first;

      int max = 0;
      for (map<GOBS, int>::iterator itGOBS = itBAND->second.begin(); itGOBS != itBAND->second.end(); itGOBS++) {
        string str = gobs2str(itGOBS->first);
        if (str.compare(0, 1, "C") != 0 && str.compare(0, 1, "L") != 0) continue;
        if (itGOBS->second > max) {
          max = itGOBS->second;
          _signals[prn][band] = str2gobsattr(str);
          //	      cout << "Fixed Types: " << prn << " " << str << endl;
        }
      }
    }
  }
}
   
// Get settings
// -----------------------------
int t_gspp::_get_settings()
{
//   _sats          = dynamic_cast<t_gsetgnss*>(_set)->sat();
   _tropo_est     = dynamic_cast<t_gsetproc*>(_set)->tropo();
   _iono_est      = dynamic_cast<t_gsetproc*>(_set)->iono();
   _tropo_grad    = dynamic_cast<t_gsetproc*>(_set)->tropo_grad();
   _tropo_slant   = dynamic_cast<t_gsetproc*>(_set)->tropo_slant();
   _ztd_mf        = dynamic_cast<t_gsetproc*>(_set)->tropo_mf();
   _grd_mf        = dynamic_cast<t_gsetproc*>(_set)->grad_mf();   
   _crd_est       = dynamic_cast<t_gsetproc*>(_set)->crd_est();   
   _sampling      = dynamic_cast<t_gsetgen*>(_set)->sampling();
   _scale         = dynamic_cast<t_gsetgen*>(_set)->sampling_scalefc();
   _minElev       = dynamic_cast<t_gsetproc*>(_set)->minimum_elev();
   _sig_init_crd  = dynamic_cast<t_gsetproc*>(_set)->sig_init_crd();
   _sig_init_ztd  = dynamic_cast<t_gsetproc*>(_set)->sig_init_ztd();
   _sig_init_vion = dynamic_cast<t_gsetproc*>(_set)->sig_init_vion();
   _sig_init_grd  = dynamic_cast<t_gsetproc*>(_set)->sig_init_grd();  
   _sig_init_glo  = dynamic_cast<t_gsetproc*>(_set)->sig_init_glo();
   _sig_init_gal  = dynamic_cast<t_gsetproc*>(_set)->sig_init_gal();
   _sig_init_bds  = dynamic_cast<t_gsetproc*>(_set)->sig_init_bds();
   _sig_init_qzs  = dynamic_cast<t_gsetproc*>(_set)->sig_init_qzs();   

   _sigCodeGPS    = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(GPS);
   _sigCodeGLO    = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(GLO);
   _sigCodeGAL    = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(GAL);
   _sigCodeBDS    = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(BDS);
   _sigCodeQZS    = dynamic_cast<t_gsetgnss*>(_set)->sigma_C(QZS);
   _sigPhaseGPS   = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(GPS);
   _sigPhaseGLO   = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(GLO);
   _sigPhaseGAL   = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(GAL);
   _sigPhaseBDS   = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(BDS);
   _sigPhaseQZS   = dynamic_cast<t_gsetgnss*>(_set)->sigma_L(QZS);

   _pos_kin       = dynamic_cast<t_gsetproc*>(_set)->pos_kin();
   _weight        = dynamic_cast<t_gsetproc*>(_set)->weighting();
   _observ        = dynamic_cast<t_gsetproc*>(_set)->obs_combin();
   _use_ecl       = dynamic_cast<t_gsetproc*>(_set)->use_eclipsed();     

   if( _sampling == 0 ) _sampling = dynamic_cast<t_gsetgen*>(_set)->sampling_default();

  return 1;
}

// Set Out
// -----------------------------
void t_gspp::_setOut()
{
  string tmp( dynamic_cast<t_gsetout*>(_set)->outputs("res") );
  if( !tmp.empty() ){
    substitute(tmp,"$(rec)",_site, false);
    _res = new t_glog;
    _res->tsys(t_gtime::GPS);
    _res->mask(tmp);
    _res->append( dynamic_cast<t_gsetout*>(_set)->append() );
  }
}

// Set log
// -----------------------------
void t_gspp::_setLog()
{  
  string tmp( dynamic_cast<t_gsetout*>(_set)->outputs("ppp") );
  if( !tmp.empty() ){
    substitute(tmp,"$(rec)",_site, false);

	t_gtime beg = dynamic_cast<t_gsetgen*>(_set)->beg();
	if( beg==FIRST_TIME){ beg = t_gtime::current_time(t_gtime::GPS); } // real-time model

	substitute(tmp, "$(rec)", _site, false);
	substitute(tmp, "$(doy)", int2str(beg.doy()), false);

    _log = new t_glog;
    _log->tsys(t_gtime::GPS);
    _log->mask(tmp);
    _log->append( dynamic_cast<t_gsetout*>(_set)->append() );
    _log->verb(   dynamic_cast<t_gsetout*>(_set)->verb() );
  }
}   
   
   
} // namespace
