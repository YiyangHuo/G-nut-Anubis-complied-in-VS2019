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

#include "gset/gsetnwm.h"
#include "gall/gallsurf.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetnwm::t_gsetnwm() 
 : t_gsetbase()
{
  _set.insert(XMLKEY_NWM);
  _min_rad      =       0;         // switched off (currently not used)
  _min_lat      = MIN_LAT;
  _max_lat      = MAX_LAT;
  _min_lon      = MIN_LON;
  _max_lon      = MAX_LON;
  _grid_mult    =     1.0;         // 4 grid points around requested point
  _decay_min    = DECAY_MIN;       // minimum allowed value for decay parameters
  _decay_max    = DECAY_MAX;       // maximum allowed value for decay parameters
  _vert_mult    =       0;         // no high-densification
  _vert_temp    = "linear";        // TEMP_LINEAR
  _vert_zwd     = "de";            // W_DE
  _surf_zwd     = "integral";      // W_INTEGR
  _vert_adj     = "lmm";           // LMM
  _vert_fit     = "pow_surf";      // POW_SURF
  _interp_time  = "linear";        // LINEAR
  _interp_vert  = "interpolate";   // INTEGRATE
  _interp_space = "ver2hor";       // VER2HOR
  _interp_plane = "bilinear";      // BILINEAR

  _m_refr = REFR_COEF();
  t_map_refr::iterator it;
  for( it = _m_refr.begin(); it != _m_refr.end(); ++it ) _refr_list  += " " + it->first;
  transform(_refr_list.begin(), _refr_list.end(), _refr_list.begin(), ::toupper);
  _refr_name    = "BEVIS";
}


// Destructor
// ----------
t_gsetnwm::~t_gsetnwm()
{}


// Return set
// ----------
set<string> t_gsetnwm::param()
{
  _gmutex.lock();
   
  set<string> tmp = t_gsetbase::_setval( XMLKEY_NWM, "param" );

  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetnwm::point()
{
  _gmutex.lock();
   
  set<string> tmp = t_gsetbase::_setval( XMLKEY_NWM, "point" );

  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetnwm::surface()
{
  _gmutex.lock();
   
  set<string> tmp = t_gsetbase::_setval( XMLKEY_NWM, "surface" );

  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetnwm::profile()
{
  _gmutex.lock();
   
  set<string> tmp = t_gsetbase::_setval( XMLKEY_NWM, "profile" );

  _gmutex.unlock(); return tmp;
}

// Return set
// ----------
set<string> t_gsetnwm::assess()
{
  _gmutex.lock();
   
  set<string> tmp = t_gsetbase::_setval( XMLKEY_NWM, "assess" );

  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetnwm::clean()
{
  _gmutex.lock();
   
  set<string> tmp = t_gsetbase::_setval( XMLKEY_NWM, "clean" );

  _gmutex.unlock(); return tmp;
}

/*
// Get level request
// ----------
set<string> t_gsetnwm::surface()
{
  set<string> tmp;
  string str;
  string id = "level";
  string ID = id; transform(ID.begin(), ID.end(), ID.begin(), ::toupper);

  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).first_child(); node; node = node.next_sibling() ){
    if( node.name() == id || node.name() == ID ){
      istringstream is( node.child_value() );
      while( is >> str && !is.fail() ){
	transform(str.begin(), str.end(), str.begin(), ::toupper);
	tmp.insert( str );
//	cout << "LEVEL = |" << str << "|\n";
      }
    }
  }
  return tmp;
}


// Get level request
// ----------
set<string> t_gsetnwm::param()
{
  set<string> tmp;
  string str;
  string id = "param";
  string ID = id; transform(ID.begin(), ID.end(), ID.begin(), ::toupper);

  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).first_child(); node; node = node.next_sibling() ){
    if( node.name() == id || node.name() == ID ){
      istringstream is( node.child_value() );
      while( is >> str && !is.fail() ){
	transform(str.begin(), str.end(), str.begin(), ::toupper);
	tmp.insert( str );
//	cout << "SURFACE = |" << str << "|\n";
      }
    }
  }
  return tmp;
}
*/

// Return value
// ----------
double t_gsetnwm::grid_mult()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("grid_mult").as_double();

  _gmutex.unlock(); return tmp;
}
   
// Return value
// ----------
double t_gsetnwm::min_rad()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("min_rad").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetnwm::min_lat()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("min_lat").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetnwm::max_lat()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("max_lat").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetnwm::min_lon()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("min_lon").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetnwm::max_lon()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("max_lon").as_double();

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
double t_gsetnwm::decay_min()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("decay_min").as_double();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
double t_gsetnwm::decay_max()
{
  _gmutex.lock();
   
  double tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("decay_max").as_double();

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
int t_gsetnwm::vert_mult()
{
  _gmutex.lock();

  int tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("vert_mult").as_int();

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
NWMFIT t_gsetnwm::vert_fit()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("vert_fit").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  if( tmp.empty() ) tmp = _vert_fit; // default

  NWMFIT fit = POW_SURF;
  if(      tmp.compare("LOG")      == 0 ) fit = LOG;
  else if( tmp.compare("POWSURF")  == 0 ) fit = POW_SURF;
  else if( tmp.compare("POW_SURF") == 0 ) fit = POW_SURF;
  else if( tmp.compare("POW")      == 0 ) fit = POW;
  else if( _log ){ _log->comment(0,"gsetnwm","Warning not valid option - vert_fit: "+tmp+" [used:"+_vert_fit+"]" ); }
       else{               cerr << "gsetnwm - Warning not valid option - vert_fit: "+tmp+" [used:"+_vert_fit+"]\n"; }
   
  _gmutex.unlock(); return fit;
}


// Return value
// ----------
NWMADJ t_gsetnwm::vert_adj()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("vert_adj").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  if( tmp.empty() ) tmp = _vert_adj; // default

  NWMADJ adj = LMM;
  if(      tmp.compare("LMM") == 0 ) adj = LMM;
  else if( tmp.compare("LSQ") == 0 ) adj = LSQ;
  else if( _log ){ _log->comment(0,"gsetnwm","Warning not valid option - vert_adj: "+tmp+" [used:"+_vert_adj+"]" ); }
       else{               cerr << "gsetnwm - Warning not valid option - vert_adj: "+tmp+" [used:"+_vert_adj+"]\n"; }
		
  _gmutex.unlock(); return adj;
}


// Return value
// ----------
PROFTEMP t_gsetnwm::vert_temp()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("vert_temp").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  if( tmp.empty() ) tmp = _vert_temp; // default

  PROFTEMP t = dT_LINEAR;
  if(      tmp.compare("LINEAR")   == 0 ) t = dT_LINEAR;
  else if( tmp.compare("CONST")    == 0 ) t = dT_CONST;
  else if( tmp.compare("ZERO")     == 0 ) t = dT_ZERO;
  else if( _log ){ _log->comment(0,"gsetnwm","Warning not valid option - vert_temp: "+tmp+" [used:"+_vert_temp+"]" ); }
       else{               cerr << "gsetnwm - Warning not valid option - vert_temp: "+tmp+" [used:"+_vert_temp+"]\n"; }
		
  _gmutex.unlock(); return t;
}


// Return value
// ----------
PROFZWD t_gsetnwm::vert_zwd()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("vert_zwd").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  if( tmp.empty() ) tmp = _vert_zwd; // default

  PROFZWD zwd = dW_DE;
  if(      tmp.compare("DE_MOPS") == 0 ) zwd = dW_DE_MOPS;
  else if( tmp.compare("DE")      == 0 ) zwd = dW_DE;
  else if( tmp.compare("MOPS")    == 0 ) zwd = dW_MOPS;
  else if( _log ){ _log->comment(0,"gsetnwm","Warning not valid option - vert_zwd: "+tmp+" [used:"+_vert_zwd+"]" ); }
       else{               cerr << "gsetnwm - Warning not valid option - vert_zwd: "+tmp+" [used:"+_vert_zwd+"]\n"; }
		
  _gmutex.unlock(); return zwd;
}


// Return value
// ----------
SURFZWD t_gsetnwm::surf_zwd()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("surf_zwd").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  if( tmp.empty() ) tmp = _surf_zwd; // default

  SURFZWD zwd = W_INTEGR;
  if(      tmp.compare("INTEGRAL") == 0 ||
           tmp.compare("INTEGR")   == 0 ||
           tmp.compare("ZWD")      == 0 ) zwd = W_INTEGR;
  else if( tmp.compare("MOPS")     == 0 ) zwd = W_MOPS;
  else if( tmp.compare("DE_G")     == 0 ) zwd = W_DE_GAMMA;
  else if( tmp.compare("DE_L")     == 0 ) zwd = W_DE_LAMBDA;
  else if( tmp.compare("DE_A")     == 0 ) zwd = W_DE_AUTO;
  else if( tmp.compare("DE")       == 0 ) zwd = W_DE;
  else if( tmp.compare("ESA_L")    == 0 ) zwd = W_ESA_L;
  else if( tmp.compare("ESA")      == 0 ) zwd = W_ESA;
  else if( _log ){ _log->comment(0,"gsetnwm","Warning not valid option - surf_zwd: "+tmp+" [used:"+_surf_zwd+"]" ); }
       else{               cerr << "gsetnwm - Warning not valid option - surf_zwd: "+tmp+" [used:"+_surf_zwd+"]\n"; }

  _gmutex.unlock(); return zwd;
}


// Return value
// ----------
string t_gsetnwm::refr_coef()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("refr_coef").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  if( tmp.empty() )
     tmp = _refr_name;
  else if( _m_refr.find(tmp) == _m_refr.end() ){
     tmp = _refr_name;
     cerr << "Refractivity coefficients not supported: "+tmp+". Used Bevis_1994.\n";
  }

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetnwm::interp_time()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("interp_time").as_string();

  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetnwm::interp_vert()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("interp_vert").as_string();

  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);

  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetnwm::interp_space()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("interp_space").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetnwm::interp_plane()
{
  _gmutex.lock();
     
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_NWM).attribute("interp_plane").as_string();
   
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
   
  _gmutex.unlock(); return tmp;
}


// settings check
// ----------
void t_gsetnwm::check()
{
  _gmutex.lock();
     
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_NWM);

  // check existence of attributes
  _default_attr(node,"min_rad",      _min_rad);
  _default_attr(node,"min_lat",      _min_lat);
  _default_attr(node,"max_lat",      _max_lat);
  _default_attr(node,"min_lon",      _min_lon);
  _default_attr(node,"max_lon",      _max_lon);
  _default_attr(node,"grid_mult",    _grid_mult);
  _default_attr(node,"decay_min",    _decay_min);
  _default_attr(node,"decay_max",    _decay_max);
  _default_attr(node,"refr_coef",    _refr_name);
  _default_attr(node,"surf_zwd",     _surf_zwd);
  _default_attr(node,"vert_mult",    _vert_mult);
  _default_attr(node,"vert_temp",    _vert_temp);
  _default_attr(node,"vert_zwd",     _vert_zwd);
  _default_attr(node,"vert_adj",     _vert_adj);
  _default_attr(node,"vert_fit",     _vert_fit);
  _default_attr(node,"interp_time",  _interp_time);
  _default_attr(node,"interp_vert",  _interp_vert);
  _default_attr(node,"interp_space", _interp_space);
  _default_attr(node,"interp_plane", _interp_plane);

   _gmutex.unlock(); return;
}


// help body
// ----------
void t_gsetnwm::help()
{
  _gmutex.lock();
   
  cerr << " <nwm"
//     << "   min_rad=\""      <<  _min_rad       << "\""     
       << "   min_lat=\""      <<  _min_lat       << "\""
       << "   max_lat=\""      <<  _max_lat       << "\""
       << "   min_lon=\""      <<  _min_lon       << "\""
       << "   max_lon=\""      <<  _max_lon       << "\""
       << "   grid_mult=\""    <<  _grid_mult     << "\""
       << "   refr_coef=\""    <<  _refr_name     << "\""
       << "   decay_min=\""    <<  _decay_min     << "\""
       << "   decay_max=\""    <<  _decay_max     << "\""
       << "   vert_mult=\""    <<  _vert_mult     << "\""
       << "   vert_temp=\""    <<  _vert_temp     << "\""
       << "   vert_zwd=\""     <<  _vert_zwd      << "\""
       << "   vert_adj=\""     <<  _vert_adj      << "\""     
       << "   vert_fit=\""     <<  _vert_fit      << "\""
       << "   surf_zwd=\""     <<  _surf_zwd      << "\""
       << "   interp_time\""   <<  _interp_time   << "\""
       << "   interp_vert\""   <<  _interp_vert   << "\""
       << "   interp_space\""  <<  _interp_space  << "\""
//     << "   interp_plane\""  <<  _interp_plane  << "\""
       << " />\n";

  cerr << "   <param> ZHD ZWD IWV PRES TEMP TM DE DW DT   </param>    \t\t <!-- output parameters (etc) -->\n"
       << "   <profile> surf full half high               </profile>  \t\t <!-- output profiles -->\n"
//     << "   <surface> geoid elipsoid surface model      </surface>  \t\t <!-- output reference surface -->\n"
//     << "   <assess> TM_SURF ZWD_SURF ZWD_VERT TRO_TIES </assess>   \t\t <!-- output assessment functions -->\n"
//     << "   <point>  grid point                         </point>    \t\t <!-- output points -->\n"
       << " </nwm>";

  cerr << "\t<!-- NWM procesing description:\n"
//     << "\t min_rad        double                       .. radius factor     (select nearest grid points)\n"
       << "\t min_lat        double                       .. minimum latitude  (region)\n"
       << "\t max_lat        double                       .. maximum latitude  (region)\n"
       << "\t min_lon        double                       .. minimum longitude (region)\n"
       << "\t max_lon        double                       .. maximum longitude (region)\n"
       << "\t grid_mult      double                       .. multiplicator for getting grid dLat/dLon around points\n"
       << "\t decay_min      double                       .. minimum value allowed for vertical decay parameters\n"
       << "\t decay_max      double                       .. maximum value allowed for vertical decay parameters\n"
       << "\t refr_coef      string                       .. refractivity coefficients: "+_refr_list+"\n"
       << "\t vert_mult      [0|1,2,..]                   .. multiplicator for increase of vertical resolution\n"
       << "\t vert_temp      [linear|const|zero]          .. method of vertical T   representation\n"
       << "\t vert_zwd       [de|de_mops|mops]            .. method of vertical ZWD representation\n"
       << "\t vert_adj       [lmm|lsq]                    .. method of vertical adjustement\n"
       << "\t vert_fit       [pow_surf|pow|log]           .. method of vertical approximation\n"
       << "\t surf_zwd       [zwd|mops|de|de_G|de_L|de_X] .. method of surface  ZWD calculation\n"
       << "\t interp_time    [linear|spline]              .. method of temporal approximation\n"
       << "\t interp_vert    [interpolate|scale]          .. method of vertical approximation\n"
       << "\t interp_space   [ver2hor|hor2ver]            .. method of spatial  approximation\n"
//     << "\t interp_plane   [bilinear|bispline]          .. method of planar   approximation\n"
       << "\t -->\n\n";

  _gmutex.unlock(); return;
}

} // namespace
