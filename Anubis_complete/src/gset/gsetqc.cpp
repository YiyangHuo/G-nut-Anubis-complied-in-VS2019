
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
#include <sstream>
#include <algorithm>

#include "gset/gsetqc.h"
#include "gio/gxtrqc.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetqc::t_gsetqc() 
 : t_gsetbase()
{
  _set.insert(XMLKEY_QC);
  _summ           = 1;
  _head           = 1;
  _stat           = 1;
  _gaps           = 1;
  _band           = 1;
  _prep           = 1;
  _elev           = 1;
  _mult           = 1;
  _stnr           = 1;
  _sinf           = 1;
  _gkpi           = 1;
  _calc           = 1;
  _step           = STP;         // [s]
  _tgap           = GAP;         // [s]
  _tpcs           = PCS;         // [s]
  _nsat           = NSAT;        // number of columns for sat-specific reporting
  _mp_nepochs     = MP_NEPOCHS;  // # of epochs
  _mp_limit       = MP_LIMIT;    // sigma-multiplicator for MP cycle-slip & outlier detection
  _mp_all         = false;
  _ele_cut        = 10.0;
  _pos_cut        = 5.0;
  _pos_int        = 900;         // [s]
  _pos_kin        = false;
  _ele_app        = false;
  _ele_new        = true;
  _sat_rec        = false;
  
  _use_health     = ALL_HEALTH;

  _dV_lim         = 5;
  _dH_lim         = 5;
  _dG_lim         = 6;  
}


// Destructor
// ----------
t_gsetqc::~t_gsetqc()
{}


// Return value
// ----------
int   t_gsetqc::summ()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_sum").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int   t_gsetqc::head()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_hdr").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int   t_gsetqc::stat()
{
  _gmutex.lock();

  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_obs").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int   t_gsetqc::gaps()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_gap").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int   t_gsetqc::band()
{
  _gmutex.lock();
  
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_bnd").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int   t_gsetqc::prep()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_pre").as_int();

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
int   t_gsetqc::elev()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_ele").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int   t_gsetqc::calc()
{
  _gmutex.lock();
  
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_est").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int   t_gsetqc::mult()
{
  _gmutex.lock();
  
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_mpx").as_int();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
int   t_gsetqc::stnr()
{
  _gmutex.lock();
  
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_snr").as_int();

  _gmutex.unlock(); return tmp;
}   

// Return value
// ----------
int t_gsetqc::sinf()
{
  _gmutex.lock();
  
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_sat").as_int();

  _gmutex.unlock(); return tmp;
}   
   
// Return value
// ----------
int t_gsetqc::gkpi()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sec_kpi").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double    t_gsetqc::step()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("int_stp").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int    t_gsetqc::tgap()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("int_gap").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int    t_gsetqc::tpcs()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("int_pcs").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int    t_gsetqc::nsat()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("col_sat").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int    t_gsetqc::mp_nepochs()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("mpx_nep").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
bool t_gsetqc::mp_all()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("mpx_all").as_bool();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetqc::mp_limit()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("mpx_lim").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetqc::ele_cut()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("ele_cut").as_double();

  if(tmp <  0.0) tmp =  0.0;
  if(tmp > 25.0) tmp = 25.0;

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetqc::pos_cut()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("ele_cut").as_double();

  if(tmp <  0.0) tmp =  0.0;
  if(tmp > 75.0) tmp = 75.0;

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
bool t_gsetqc::pos_kin()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("pos_kin").as_bool();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int t_gsetqc::pos_int()
{
  _gmutex.lock();
   
  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("pos_int").as_int();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
bool t_gsetqc::ele_app()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("ele_app").as_bool();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
bool t_gsetqc::ele_new()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("ele_new").as_bool();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
double t_gsetqc::dV_lim()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("dV_lim").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetqc::dH_lim()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("dH_lim").as_double();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
double t_gsetqc::dG_lim()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("dG_lim").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
bool t_gsetqc::sat_rec()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("sat_rec").as_bool();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
USE_HEALTH t_gsetqc::useHealth()
{
   _gmutex.lock();

   string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_QC).attribute("use_health").value();

   transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
   USE_HEALTH UH = str2useHealth(tmp);
   
   if(UH == DEF_HEALTH){
      xml_node parent = _doc.child(XMLKEY_ROOT);
      xml_node node = _default_node(parent, XMLKEY_QC);
      tmp = useHealth2str(_use_health);
      _default_attr(node,"use_health", tmp);
      
      UH = _use_health;
  }   
   
  _gmutex.unlock(); return UH;
}

// settings check
// ----------
void t_gsetqc::check()
{
  _gmutex.lock();
     
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_QC);

  // check existence of attributes
  _default_attr(node,"sec_sum",     _summ);
  _default_attr(node,"sec_hdr",     _head);
  _default_attr(node,"sec_est",     _calc);
  _default_attr(node,"sec_obs",     _stat);
  _default_attr(node,"sec_gap",     _gaps);
  _default_attr(node,"sec_bnd",     _band);
  _default_attr(node,"sec_pre",     _prep);
  _default_attr(node,"sec_ele",     _elev);
  _default_attr(node,"sec_mpx",     _mult);
  _default_attr(node,"sec_sat",     _sinf);
  
  _default_attr(node,"int_stp",     _step);
  _default_attr(node,"int_gap",     _tgap);
  _default_attr(node,"int_pcs",     _tpcs);
  _default_attr(node,"col_sat",     _nsat);
  _default_attr(node,"mpx_nep",     _mp_nepochs);
  _default_attr(node,"mpx_lim",     _mp_limit);
  _default_attr(node,"mpx_all",     _mp_all);
  _default_attr(node,"pos_kin",     _pos_kin);
  _default_attr(node,"pos_int",     _pos_int);
  _default_attr(node,"pos_cut",     _pos_cut);
  _default_attr(node,"ele_cut",     _ele_cut);
  _default_attr(node,"ele_new",     _ele_new);
  _default_attr(node,"ele_app",     _ele_app);
  _default_attr(node,"sat_rec",     _sat_rec);
  _default_attr(node,"use_health",  _use_health);

  _default_attr(node,"dV_lim",      _dV_lim);
  _default_attr(node,"dH_lim",      _dH_lim);
  _default_attr(node,"dG_lim",      _dG_lim);  
  
  _gmutex.unlock(); return;
}


// help body
// ----------
void t_gsetqc::help()
{
  _gmutex.lock();

  cerr << " <qc \n"
       << "   sec_sum=\""     <<  _summ          << "\" \n"
       << "   sec_hdr=\""     <<  _head          << "\" \n"
       << "   sec_obs=\""     <<  _stat          << "\" \n"
       << "   sec_est=\""     <<  _calc          << "\" \n"
       << "   sec_gap=\""     <<  _gaps          << "\" \n"
       << "   sec_bnd=\""     <<  _band          << "\" \n"
       << "   sec_pre=\""     <<  _prep          << "\" \n"
       << "   sec_ele=\""     <<  _elev          << "\" \n"
       << "   sec_mpx=\""     <<  _mult          << "\" \n"
       << "   sec_snr=\""     <<  _stnr          << "\" \n"
       << "   sec_sat=\""     <<  _sinf          << "\" \n"       

       << "   int_stp=\""     <<  _step          << "\" \n"
       << "   int_gap=\""     <<  _tgap          << "\" \n"
       << "   int_pcs=\""     <<  _tpcs          << "\" \n"
       << "   col_sat=\""     <<  _nsat          << "\" \n"
       << "   mpx_nep=\""     <<  _mp_nepochs    << "\" \n"
       << "   mpx_lim=\""     <<  _mp_limit      << "\" \n"
//     << "   mpx_all=\""     <<  _mp_all        << "\" \n"
       << "   pos_kin=\""     <<  _pos_kin       << "\" \n"
       << "   pos_int=\""     <<  _pos_int       << "\" \n"
       << "   pos_cut=\""     <<  _pos_cut       << "\" \n"
       << "   ele_cut=\""     <<  _ele_cut       << "\" \n"
//     << "   ele_new=\""     <<  _ele_new       << "\" \n"
//     << "   ele_app=\""     <<  _ele_app       << "\" \n"
       << "   sat_rec=\""     <<  _sat_rec       << "\" \n"
       << "   use_health=\""  <<  _use_health    << "\" \n"

       << "   dV_lim=\""      <<  _dV_lim        << "\" \n"
       << "   dH_lim=\""      <<  _dH_lim        << "\" \n"
       << "   dG_lim=\""      <<  _dG_lim        << "\" \n"  
       << " />\n";

  cerr << "\t<!-- quality check description:\n"
       << "\t sec_sum     [0-9]    .. summary statistics\n"
       << "\t sec_hdr     [0-9]    .. header metadata check\n"
       << "\t sec_est     [0-9]    .. estimated values\n"
       << "\t sec_obs     [0-9]    .. observation statistics\n"
       << "\t sec_gap     [0-9]    .. gap & pieces\n"
       << "\t sec_bnd     [0-9]    .. observation bands\n"
       << "\t sec_pre     [0-9]    .. cycle-slip, clock-jumps\n"
       << "\t sec_ele     [0-9]    .. azimuth/elevation (if navigation)\n"
       << "\t sec_mpx     [0-9]    .. multipath calculation\n"
       << "\n sec_snr     [0-9]    .. signal-to-noise ratio\n"
       << "\n sec_sat     [0-9]    .. satellite information\n"

       << "\t int_stp     int[s]   .. interval for time-spacing\n"
       << "\t int_gap     int[s]   .. interval for gap identification\n"
       << "\t int_pcs     int[s]   .. interval for smalle pieces identification\n"
       << "\t col_sat     int[#]   .. number of columns for sat-specific reporting\n"
       << "\t mpx_nep     int[#]   .. number of epochs for multipath calculation\n"
//     << "\t mpx_all     bool     .. high-resolution multipath: estimate all or interpolate\n"
       << "\t mpx_lim     double   .. sigma-multiplicator for MP cycle-slip & outlier detection\n"
       << "\t pos_kin     bool     .. kinematic receiver (true = kinematic)\n"
       << "\t pos_int     int      .. positioning interval\n"
       << "\t pos_cut     double   .. positioning elevation angle cut-off (degrees)\n"
       << "\t ele_cut     double   .. user elevation cut-off (only for expt/have, degrees)\n"
//     << "\t ele_new     bool     .. new strategy for cut-off + horizon estimates\n"
//     << "\t ele_app     bool     .. approximated elevations (true = less precise)\n"
       << "\t sat_rec     bool     .. expected satellite (all | receiver tracking)\n"
       << "\t use_health  enum     .. using method of satellite healhy status (position|statistics|all)\n"
  
       << "\t dV_lim      double   .. vertical limit position \n"
       << "\t dH_lim      double   .. horizontal limit position \n"
       << "\t dG_lim      doubel   .. GDOP limit position \n"
       << "\t -->\n\n";

  _gmutex.unlock(); return;
}


// convert str to USE_HEALTH enum
// ----------   
USE_HEALTH t_gsetqc::str2useHealth(string s)
{
  USE_HEALTH UH;
   
       if(s.compare("POSITION")   == 0) UH = POS_HEALTH;
  else if(s.compare("STATISTICS") == 0) UH = STT_HEALTH;
  else if(s.compare("ALL")        == 0) UH = ALL_HEALTH;
  else {
    UH = DEF_HEALTH;
    stringstream ostr;
    ostr << "Unsupported combination (" << s <<  ")! Used default value (" << useHealth2str(_use_health) << ")";
    _add_log("gsetqc", ostr.str());
  }
  return UH;
}

// convert USE_HEALTH enum to string
// ---------
string t_gsetqc::useHealth2str(USE_HEALTH UH)
{
  switch (UH) {
    case POS_HEALTH:   return "POSITION";
    case STT_HEALTH:   return "STATISTICS";
    case ALL_HEALTH:   return "ALL";
    case DEF_HEALTH:   return "NOT DEFINED";
  }
  return "";
}   

} // namespace
