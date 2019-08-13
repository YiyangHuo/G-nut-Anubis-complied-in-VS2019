
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

#include <iomanip>
#include <sstream>

#include "gset/gsetflt.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetflt::t_gsetflt() 
 : t_gsetbase()
{
   _set.insert(XMLKEY_FLT);
  _noise_clk      = 1000.0;
  _noise_crd      = 100.0;   
  _rndwk_glo      = 20.0;
  _rndwk_gal      = 20.0;
  _rndwk_bds      = 20.0;
  _rndwk_qzs      = 20.0;   
  _rndwk_ztd      = 3.0;
  _noise_vion     = 0.1;
  _rndwk_grd      = 0.3;   
  _method_flt     = "kalman";
  _method_smt     = "RTSSVD";
  _reset_amb      = 0;
  _reset_par      = 0;
  _delay          = 0;
  _smooth         = false;
}


// Destructor
// ----------
t_gsetflt::~t_gsetflt()
{}


// Return value
// ----------
string t_gsetflt::method_flt()
{
   _gmutex.lock();
   
   string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("method_flt").value();

   _gmutex.unlock(); return tmp;
}

// Return value
// ----------
string t_gsetflt::method_smt()
{
   _gmutex.lock();
   
   string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("method_smt").value();
   
   _gmutex.unlock(); return tmp;
}


// Return value
// ----------
double t_gsetflt::noise_clk()
{
   _gmutex.lock();
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("noise_clk").as_double();
   
   _gmutex.unlock(); return tmp;
}


// Return value
// ----------
double t_gsetflt::noise_crd()
{
   _gmutex.lock();
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("noise_crd").as_double();
   
   _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetflt::rndwk_glo()
{
   _gmutex.lock();
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("rndwk_glo").as_double();

   _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetflt::rndwk_gal()
{
  _gmutex.lock();
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("rndwk_gal").as_double();

   _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetflt::rndwk_bds()
{
   _gmutex.lock();
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("rndwk_bds").as_double();
   
   _gmutex.unlock(); return tmp;
}   
   
// Return value
// ----------
double t_gsetflt::rndwk_qzs()
{
   _gmutex.lock();
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("rndwk_qzs").as_double();

   _gmutex.unlock(); return tmp;
}   
   
// Return value
// ----------
double t_gsetflt::rndwk_ztd()
{
   _gmutex.lock(); 
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("rndwk_ztd").as_double();

   _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
double t_gsetflt::noise_vion()
{
   _gmutex.lock(); 
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("noise_vion").as_double();

   _gmutex.unlock(); return tmp;
}   

// Return value
// ----------
double t_gsetflt::rndwk_grd()
{
   _gmutex.lock(); 
   
   double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("rndwk_grd").as_double();

   _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int t_gsetflt::reset_amb()
{
   _gmutex.lock(); 
   
   int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("reset_amb").as_int();

   _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int t_gsetflt::reset_par()
{
   _gmutex.lock(); 
   
   int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("reset_par").as_int();

   _gmutex.unlock(); return tmp;
}

// Return value
// ----------
int t_gsetflt::smt_delay()
{
   _gmutex.lock(); 
   
   int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("smt_delay").as_int();

   _gmutex.unlock(); return tmp;
}

   
   // Return value
// ----------
bool t_gsetflt::smooth()
{
   _gmutex.lock(); 
   
   bool tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_FLT).attribute("smooth").as_bool();

   _gmutex.unlock(); return tmp;
}
   
// settings check
// ----------
void t_gsetflt::check()
{
  _gmutex.lock();

  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_FLT);

  // check existence of attributes
  _default_attr(node,"method_flt",_method_flt);
  _default_attr(node,"method_smt",_method_smt);   
  _default_attr(node,"noise_clk", _noise_clk );
  _default_attr(node,"noise_crd", _noise_crd );
  _default_attr(node,"rndwk_glo", _rndwk_glo );
  _default_attr(node,"rndwk_gal", _rndwk_gal );
  _default_attr(node,"rndwk_bds", _rndwk_bds );
  _default_attr(node,"rndwk_qzs", _rndwk_qzs );   
  _default_attr(node,"rndwk_ztd", _rndwk_ztd );
  _default_attr(node,"noise_vion",_noise_vion );
  _default_attr(node,"rndwk_grd", _rndwk_grd );   
  _default_attr(node,"reset_amb", _reset_amb );
  _default_attr(node,"reset_par", _reset_par );
  _default_attr(node,"smt_delay", _delay );
  _default_attr(node,"smooth", _smooth );   

   _gmutex.unlock(); return;
}


// help body
// ----------
void t_gsetflt::help()
{
  _gmutex.lock();

  cerr << " <filter \n"
       << "   method_flt=\""     <<  _method_flt  << "\" \n"
       << "   method_smt=\""     <<  _method_smt  << "\" \n"     
       << "   noise_clk=\""      <<  _noise_clk   << "\" \n"
       << "   noise_crd=\""      <<  _noise_crd   << "\" \n"
       << "   rndwk_ztd=\""      <<  _rndwk_ztd   << "\" \n"
       << "   noise_vion=\""     <<  _noise_vion  << "\" \n"  
       << "   rndwk_grd=\""      <<  _rndwk_grd   << "\" \n"
       << "   smt_delay=\""      <<  _delay   << "\" \n"
       << "   smooth=\""         <<  _smooth   << "\" \n"     
       << "  />\n";

  cerr << "\t<!-- filter description:\n"
       << "\t method_flt    .. type of filtering method (SRCF, Kalman)\n"
       << "\t method_smt    .. type of smoothing method (RTS)\n"     
       << "\t noise_clk     .. white noise for clocks \n"
       << "\t noise_crd     .. white noise for coordinates \n"
       << "\t rndwk_ztd     .. random walk process for ZTD [mm/sqrt(hour)] \n"
       << "\t noise_vion    .. white noise process for VION [mm/sqrt(hour)] \n"  
       << "\t -->\n\n";
   
   _gmutex.unlock(); return;   
}

} // namespace
