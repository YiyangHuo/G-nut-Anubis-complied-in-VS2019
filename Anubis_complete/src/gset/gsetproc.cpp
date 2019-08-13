
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

#include "gset/gsetproc.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetproc::t_gsetproc() 
 : t_gsetbase()
{
   _set.insert(XMLKEY_PROC);
  _phase          = true;
  _tropo          = true;
  _iono           = true;  
  _tropo_grad     = false;
  _tropo_slant    = false;
  _tropo_model    = SAASTAMOINEN;
  _tropo_mf       = GMF;
  _iono_mf        = ICOSZ;  
  _grad_mf        = TILTING;
  _obs_weight     = SINEL;
  _res_type       = RES_NORM;
  _obs_combin     = IONO_FREE;

  _sig_init_ztd   = 0.1;
  _sig_init_vion  = 10;
  _sig_init_grd   = 0.0005;  
  _sig_init_amb   = 1000.0;
  _sig_init_crd   = 100.0;
  _sig_init_glo   = 1000.0;
  _sig_init_gal   = 1000.0;   
  _sig_init_bds   = 1000.0;
  _sig_init_qzs   = 1000.0;
  _minimum_elev   = 10;
  _max_res_norm   = 3; 
  _crd_est        = "EST";
  _pos_kin        = false;
  _use_eclipsed   = false;
  _auto_band      = false;
}


// Destructor
// ----------
t_gsetproc::~t_gsetproc()
{}


// Return value
// ----------
bool   t_gsetproc::tropo()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo").as_bool();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
bool   t_gsetproc::iono()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("iono").as_bool();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
bool   t_gsetproc::tropo_slant()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_slant").as_bool();

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
bool   t_gsetproc::tropo_grad()
{
  _gmutex.lock();

  // temporary for legacy mode with keyword "gradient"
  bool tmp = false;
  bool tmp1 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_grad").as_bool(); 
  bool tmp2 = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("gradient").as_bool();   // legacy mode
  if( tmp1 || tmp2 ) tmp = true;

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
bool   t_gsetproc::phase()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("phase").as_bool();

  _gmutex.unlock(); return tmp;
}   
   
// Return value
// ----------
bool   t_gsetproc::pos_kin()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("pos_kin").as_bool();

  _gmutex.unlock(); return tmp;
} 
   
// Return value
// ----------
bool   t_gsetproc::use_eclipsed()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("use_eclipsed").as_bool();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
bool   t_gsetproc::auto_band()
{
  _gmutex.lock();
   
  bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("auto_band").as_bool();

  _gmutex.unlock(); return tmp;
}   
   
// Return value
// ----------
double t_gsetproc::sig_init_ztd()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_ztd").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_vion()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_vion").as_double();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
double t_gsetproc::sig_init_grd()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_grd").as_double();

  _gmutex.unlock(); return tmp;
}   

// Return value
// ----------
double t_gsetproc::sig_init_crd()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_crd").as_double();

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
double t_gsetproc::sig_init_amb()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_amb").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_glo()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_glo").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetproc::sig_init_gal()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_gal").as_double();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
double t_gsetproc::sig_init_bds()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_bds").as_double();

  _gmutex.unlock(); return tmp;
}   

// Return value
// ----------
double t_gsetproc::sig_init_qzs()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("sig_init_qzs").as_double();

  _gmutex.unlock(); return tmp;
}   
   
// Return value
// ----------
double t_gsetproc::minimum_elev()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("minimum_elev").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetproc::max_res_norm()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("max_res_norm").as_double();

  _gmutex.unlock(); return tmp;
}   

// Return value
// ----------
TROPMODEL t_gsetproc::tropo_model()
{
   _gmutex.lock();
   
   string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_model").value();

   transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
   TROPMODEL TM = str2tropmodel(tmp);
   
   if(TM == DEF_TROPMODEL){
      xml_node parent = _doc.child(XMLKEY_ROOT);
      xml_node node = _default_node(parent, XMLKEY_PROC);
      tmp = tropmodel2str(_tropo_model);
      _default_attr(node,"tropo_model", tmp);
      
      TM = _tropo_model;
  }   
   
  _gmutex.unlock(); return TM;
}

// Return value
// ----------
GRDMPFUNC t_gsetproc::grad_mf()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("grad_mf").value();

  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  GRDMPFUNC MF = str2grdmpfunc(tmp);

   if(MF == DEF_GRDMPFUNC){
      xml_node parent = _doc.child(XMLKEY_ROOT);
      xml_node node = _default_node(parent, XMLKEY_PROC);
      tmp = grdmpfunc2str(_grad_mf);
      _default_attr(node,"grad_mf", tmp);
      
      MF = _grad_mf;
  }
   
  _gmutex.unlock(); return MF;
}

// Return value
// ----------   
OBSWEIGHT t_gsetproc::weighting()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("obs_weight").value();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  OBSWEIGHT WG = str2obsweight(tmp);

   if(WG == DEF_OBSWEIGHT){
      xml_node parent = _doc.child(XMLKEY_ROOT);
      xml_node node = _default_node(parent, XMLKEY_PROC);
      tmp = obsweight2str(_obs_weight);
      _default_attr(node,"obs_weight", tmp);     
      
      WG = _obs_weight;
  }
   
  _gmutex.unlock(); return WG;
}
 
// Return value
// ----------   
RESIDTYPE t_gsetproc::residuals()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("residuals").value();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  RESIDTYPE RS = str2residtype(tmp);

   if(RS == DEF_RESIDTYPE){
      xml_node parent = _doc.child(XMLKEY_ROOT);
      xml_node node = _default_node(parent, XMLKEY_PROC);
      tmp = residtype2str(_res_type);
      _default_attr(node,"residuals", tmp);     
      
      RS = _res_type;
  }
   
  _gmutex.unlock(); return RS;
}

// Return value
// ----------   
OBSCOMBIN t_gsetproc::obs_combin()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("obs_combination").value();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  OBSCOMBIN OC = str2obscombin(tmp);

   if(OC == DEF_OBSCOMBIN){
      xml_node parent = _doc.child(XMLKEY_ROOT);
      xml_node node = _default_node(parent, XMLKEY_PROC);
      tmp = obscombin2str(_obs_combin);
      _default_attr(node,"obs_combination", tmp);     
      
      OC = _obs_combin;
  }
   
  _gmutex.unlock(); return OC;
}   
   
 // Return value
// ----------
ZTDMPFUNC t_gsetproc::tropo_mf()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("tropo_mf").value();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  ZTDMPFUNC MF = str2ztdmpfunc(tmp);

  if(MF == DEF_ZTDMPFUNC){
     xml_node parent = _doc.child(XMLKEY_ROOT);
     xml_node node = _default_node(parent, XMLKEY_PROC);     
     tmp = ztdmpfunc2str(_tropo_mf);
     _default_attr(node,"tropo_mf", tmp);
     
     MF = _tropo_mf;
  }
   
  _gmutex.unlock(); return MF;
}

// Return value
// ----------
IONMPFUNC t_gsetproc::iono_mf()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("iono_mf").value();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  IONMPFUNC MF = str2ionmpfunc(tmp);

  if(MF == DEF_IONMPFUNC){
     xml_node parent = _doc.child(XMLKEY_ROOT);
     xml_node node = _default_node(parent, XMLKEY_PROC);     
     tmp = ionmpfunc2str(_iono_mf);
     _default_attr(node,"iono_mf", tmp);
     
     MF = _iono_mf;
  }
   
  _gmutex.unlock(); return MF;
}
   
// Return value
// ----------
CONSTRPAR t_gsetproc::crd_est()
{
  _gmutex.lock();
   
  CONSTRPAR constr = EST;
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_PROC).attribute("crd_constr").value();

  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  if( tmp.empty() ) tmp = _crd_est; // default

  if(      tmp.compare("FIX") == 0 ) constr = FIX;
  else if( tmp.compare("EST") == 0 ) constr = EST;
   
  _gmutex.unlock(); return constr;         
}   

// settings check
// ----------
void t_gsetproc::check()
{
  _gmutex.lock();
   
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_PROC);

  // check existence of attributes
  _default_attr(node,"phase", _phase);
  _default_attr(node,"tropo",       _tropo);
  _default_attr(node,"iono",        _iono);  
  _default_attr(node,"tropo_grad",  _tropo_grad);
  _default_attr(node,"tropo_model", _tropo_model);
  _default_attr(node,"tropo_slant", _tropo_slant);
  _default_attr(node,"tropo_mf",    ztdmpfunc2str(_tropo_mf));
  _default_attr(node,"iono_mf",     ionmpfunc2str(_iono_mf));  
  _default_attr(node,"grad_mf",     grdmpfunc2str(_grad_mf));   
  _default_attr(node,"obs_weight",  obsweight2str(_obs_weight));   
  _default_attr(node,"sig_init_crd", _sig_init_crd);
  _default_attr(node,"sig_init_ztd", _sig_init_ztd);
  _default_attr(node,"sig_init_grd", _sig_init_grd);  
  _default_attr(node,"sig_init_vion", _sig_init_vion);
  _default_attr(node,"sig_init_amb", _sig_init_amb);
  _default_attr(node,"sig_init_glo", _sig_init_glo);
  _default_attr(node,"sig_init_gal", _sig_init_gal);
  _default_attr(node,"sig_init_bds", _sig_init_bds);
  _default_attr(node,"sig_init_qzs", _sig_init_qzs);   
  _default_attr(node,"minimum_elev", _minimum_elev);
  _default_attr(node,"max_res_norm", _max_res_norm);   
  _default_attr(node,"crd_constr", _crd_est);
  _default_attr(node,"pos_kin", _pos_kin);
  _default_attr(node,"use_eclipsed", _use_eclipsed);
  _default_attr(node,"auto_band", _auto_band);  
  
  _gmutex.unlock(); return; 
}


// help body
// ----------
void t_gsetproc::help()
{
  _gmutex.lock();
   
  cerr << " <process \n"
       << "   phase=\""          <<  _phase          << "\" \n"
       << "   tropo=\""          <<  _tropo          << "\" \n"
       << "   iono=\""           <<  _iono           << "\" \n"       
       << "   tropo_grad=\""     <<  _tropo_grad     << "\" \n"     
       << "   tropo_model=\""    <<  tropmodel2str(_tropo_model)    << "\" \n"
       << "   tropo_slant=\""    <<  _tropo_slant    << "\" \n"
       << "   tropo_mf=\""       <<  ztdmpfunc2str(_tropo_mf)   << "\" \n"
       << "   iono_mf=\""        <<  ionmpfunc2str(_iono_mf)   << "\" \n"       
       << "   grad_mf=\""        <<  grdmpfunc2str(_grad_mf)    << "\" \n"
       << "   obs_weight=\""     <<  obsweight2str(_obs_weight) << "\" \n"
       << "   residuals=\""      <<  residtype2str(_res_type) << "\" \n"
       << "   obs_combination=\""<<  obscombin2str(_obs_combin) << "\" \n"	 
       << "   sig_init_crd=\""   <<  _sig_init_crd   << "\" \n"
       << "   sig_init_ztd=\""   <<  _sig_init_ztd   << "\" \n"
       << "   sig_init_vion=\""  <<  _sig_init_vion  << "\" \n"  
       << "   minimum_elev=\""   <<  _minimum_elev   << "\" \n"
       << "   max_res_norm=\""   <<  _max_res_norm   << "\" \n"     
       << " />\n";
   
  cerr << "\t<!-- process description:\n"
       << "\t phase [0,1]   .. use carrier phase data\n"
       << "\t tropo [0,1]   .. estimate troposphere\n"
       << "\t tropo_grad    .. tropospheric horizontal gradient models\n"
       << "\t tropo_model   .. tropospheric model (SAAS, GPT, ...)\n"
       << "\t tropo_slant   .. tropospheric slant delays produced\n"
       << "\t tropo_mf      .. tropospheric mapping function (COSZ, GMF, ...)\n"
       << "\t grad_mf       .. tropo gradient mapping function (TILTING,CHEN_HERRING, BAR_SEVER ...)\n"
       << "\t obs_weight    .. observation elevation dependant weighting (EQUAL, SINEL, SINEL2, SINEL4, CODPHA, MLTPTH)\n"
       << "\t residuals     .. type of residuals (RES_ORIG, RES_NORM, RES_ALL)\n"
       << "\t sig_init_crd  .. accuracy of initial coordinates [m]\n"
       << "\t sig_init_ztd  .. accuracy of initial zenith path delay [m]\n"
       << "\t sig_init_vion .. accuracy of initial vertical iono path delay [m]\n"  
       << "\t minimum_elev  .. elevation angle cut-off [degree]\n"
       << "\t max_res_norm  .. maximal normalized residuals\n"
       << "\t -->\n\n";

   _gmutex.unlock(); return;
}

// convert str to GRDMPFUNC enum
// ----------   
GRDMPFUNC t_gsetproc::str2grdmpfunc(string mf)
{
   GRDMPFUNC MF; 
   
        if(mf.compare("TILTING")      == 0) MF = TILTING;
   else if(mf.compare("CHEN_HERRING") == 0) MF = CHEN_HERRING;
   else if(mf.compare("BAR_SEVER")    == 0) MF = BAR_SEVER;
   else {
      MF = DEF_GRDMPFUNC;
      stringstream ostr;
      ostr << "Unsupported GRD mapping function (" << mf <<  ")! Used default value (" << grdmpfunc2str(_grad_mf) << ")";
      _add_log("gsetproc", ostr.str());
   }
      
   return MF;
}

// convert str to ZTDMPFUNC enum
// ----------   
ZTDMPFUNC t_gsetproc::str2ztdmpfunc(string mf)
{
   ZTDMPFUNC MF; 
   
        if(mf.compare("COSZ") == 0) MF = COSZ;
   else if(mf.compare("GMF")  == 0) MF = GMF;
   else {
      MF = DEF_ZTDMPFUNC;
      stringstream ostr;
      ostr << "Unsupported ZTD mapping function (" << mf <<  ")! Used default value (" << ztdmpfunc2str(_tropo_mf) << ")";
      _add_log("gsetproc", ostr.str());      
   }
   

   return MF;
}

// convert str to IONMPFUNC enum
// ----------   
IONMPFUNC t_gsetproc::str2ionmpfunc(string mf)
{
   IONMPFUNC MF; 
   
        if(mf.compare("COSZ") == 0) MF = ICOSZ;
   else if(mf.compare("QFAC") == 0) MF = QFAC;
   else if(mf.compare("NONE") == 0) MF = NONE;   
   else {
      MF = DEF_IONMPFUNC;
      stringstream ostr;
      ostr << "Unsupported ION mapping function (" << mf <<  ")! Used default value (" << ionmpfunc2str(_iono_mf) << ")";
      _add_log("gsetproc", ostr.str());      
   }
   

   return MF;
}

// convert str to OBSWEIGHT enum
// ----------   
OBSWEIGHT t_gsetproc::str2obsweight(string wg)
{
   OBSWEIGHT WG;
   
        if(wg.compare("EQUAL")  == 0) WG = EQUAL;
   else if(wg.compare("SINEL")  == 0) WG = SINEL;
   else if(wg.compare("SINEL2") == 0) WG = SINEL2;
   else if(wg.compare("SINEL4") == 0) WG = SINEL4;
   else if(wg.compare("CODPHA") == 0) WG = CODPHA;   
   else if(wg.compare("MLTPTH") == 0) WG = MLTPTH;
   else {
      WG = DEF_OBSWEIGHT;
      stringstream ostr;
      ostr << "Unsupported observation weighting model (" << wg <<  ")! Used default value (" << obsweight2str(_obs_weight) << ")";
      _add_log("gsetproc", ostr.str());      
   }
   
   
   return WG;
}
 
// convert str to TROPMODEL enum
// ----------   
TROPMODEL t_gsetproc::str2tropmodel(string tm)
{
   TROPMODEL TM;
   
        if(tm.compare("SAASTAMOINEN") == 0) TM = SAASTAMOINEN;
   else if(tm.compare("DAVIS")        == 0) TM = DAVIS;
   else if(tm.compare("HOPFIELD")     == 0) TM = HOPFIELD;
   else if(tm.compare("MOPS")         == 0) TM = MOPS;
   else if(tm.compare("GPTW")         == 0) TM = GPTW;   
   else if(tm.compare("GPT2W")        == 0) TM = GPT2W;
   else if(tm.compare("GAL27")        == 0) TM = GAL27;
   else if(tm.compare("GALTROPO27")   == 0) TM = GALTROPO27;
   else if(tm.compare("EXTERN")       == 0) TM = EXTERN;   
   else {
      TM = DEF_TROPMODEL;
      stringstream ostr;
      ostr << "Unsupported tropospheric model (" << tm <<  ")! Used default value (" << tropmodel2str(_tropo_model) << ")";
      _add_log("gsetproc", ostr.str());      
   }
   
   
   return TM;
}   
 
// convert str to RESIDTYPE enum
// ----------   
RESIDTYPE t_gsetproc::str2residtype(string rs)
{
   RESIDTYPE RS;
   
        if(rs.compare("RES_ORIG") == 0) RS = RES_ORIG;
   else if(rs.compare("RES_NORM") == 0) RS = RES_NORM;
   else if(rs.compare("RES_ALL")  == 0) RS = RES_ALL;
   else {
      RS = DEF_RESIDTYPE;
      stringstream ostr;
      ostr << "Unsupported type of residuals (" << rs <<  ")! Used default value (" << residtype2str(_res_type) << ")";
      _add_log("gsetproc", ostr.str());      
   }
   
   
   return RS;
}
 
// convert str to OBSCOMB enum
// ----------   
OBSCOMBIN t_gsetproc::str2obscombin(string oc)
{
   OBSCOMBIN OC;
   
        if(oc.compare("IONO_FREE")  == 0) OC = IONO_FREE;
   else if(oc.compare("RAW_SINGLE") == 0) OC = RAW_SINGLE;
   else if(oc.compare("RAW_DOUBLE") == 0) OC = RAW_DOUBLE;
   else if(oc.compare("RAW_ALL")    == 0) OC = RAW_ALL;
   else {
      OC = DEF_OBSCOMBIN;
      stringstream ostr;
      ostr << "Unsupported observations combination (" << oc <<  ")! Used default value (" << obscombin2str(_obs_combin) << ")";
      _add_log("gsetproc", ostr.str());      
   }
   
   
   return OC;
}
   
// convert GRDMPFUNC enum to str
// ---------
string t_gsetproc::grdmpfunc2str(GRDMPFUNC MF)
{
   switch (MF) {
    case TILTING:       return "TILTING";
    case CHEN_HERRING:  return "CHEN_HERRING";
    case BAR_SEVER:     return "BAR_SEVER";
    case DEF_GRDMPFUNC: return "NOT DEFINED";
   }
  return "";
}
   
// convert ZTDMPFUNC enum to str
// ---------
string t_gsetproc::ztdmpfunc2str(ZTDMPFUNC MF)
{
   switch (MF) {
    case COSZ:          return "COSZ";
    case GMF:           return "GMF";
    case DEF_ZTDMPFUNC: return "NOT DEFINED";
   }   
   return "";
}

// convert IONMPFUNC enum to str
// ---------
string t_gsetproc::ionmpfunc2str(IONMPFUNC MF)
{
   switch (MF) {
    case COSZ:          return "COSZ";
    case QFAC:          return "QFAC";
    case NONE:          return "NONE";    
    case DEF_IONMPFUNC: return "NOT DEFINED";
   }   
   return "";
}

// convert OBSWEIGHT enum to str
// ---------   
string t_gsetproc::obsweight2str(OBSWEIGHT WG)
{
   switch (WG) {
    case EQUAL:         return "EQUAL";
    case SINEL:         return "SINEL";
    case SINEL2:        return "SINEL2";
    case SINEL4:        return "SINEL4";
    case CODPHA:        return "CODPHA";
    case MLTPTH:        return "MLTPTH";      
    case DEF_OBSWEIGHT: return "NOT DEFINED";
   }
   return "";
}
   
// convert TROPMODEL enum to str
// ---------   
string t_gsetproc::tropmodel2str(TROPMODEL TM)
{
   switch (TM) {
    case SAASTAMOINEN:  return "SAASTAMOINEN";
    case DAVIS:         return "DAVIS";       
    case HOPFIELD:      return "HOPFIELD";
    case MOPS:          return "MOPS";         
    case GPTW:          return "GPTW";        
    case GPT2W:         return "GPT2W";         
    case GAL27:         return "GAL27";       
    case GALTROPO27:    return "GALTROPO27";
    case EXTERN:        return "EXTERN";
    case DEF_TROPMODEL: return "NOT DEFINED";	
   }
   return "";
}   

// convert RESIDTYPE enum to str
// ---------
string t_gsetproc::residtype2str(RESIDTYPE RS)
{
   switch (RS) {
    case RES_ORIG:      return "RES_ORIG";
    case RES_NORM:      return "RES_NORM";
    case RES_ALL:       return "RES_ALL";
    case DEF_RESIDTYPE: return "NOT DEFINED";
   }
  return "";
}
   
// convert OBSCOMB enum to str
// ---------
string t_gsetproc::obscombin2str(OBSCOMBIN OC)
{
   switch (OC) {
    case IONO_FREE:      return "IONO_FREE";
    case RAW_SINGLE:     return "RAW_SINGLE";
    case RAW_DOUBLE:     return "RAW_DOUBLE";
  	case RAW_ALL:        return "RAW_ALL";
    case DEF_OBSCOMBIN:  return "NOT DEFINED";
   }
  return "";
}   
   
} // namespace
