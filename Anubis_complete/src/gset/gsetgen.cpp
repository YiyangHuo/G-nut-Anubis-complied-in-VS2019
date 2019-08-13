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
#include <algorithm>

#include "gutils/gnss.h"
#include "gutils/gsys.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetgen.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetgen::t_gsetgen() 
 : t_gsetbase(),
   _dec(0)
{ 
  _set.insert(XMLKEY_GEN);
   
  // initiate STRINGS
  set<string>::const_iterator itSAT;

  t_map_sats gnss_sats = GNSS_SATS();
  t_map_sats::const_iterator itGNS;
  for( itGNS = gnss_sats.begin(); itGNS != gnss_sats.end(); ++itGNS ){
   _sys += " " + t_gsys::gsys2str( itGNS->first );
 }
}

// Destructor
// ----------
t_gsetgen::~t_gsetgen()
{}


// Return gtime
// ----------
t_gtime t_gsetgen::beg()
{
  _gmutex.lock();

  string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value( "beg" );
  substitute(str,"\n","");
  substitute(str,"\"","");

  t_gtime gt(t_gtime::GPS);
   
  if( str.empty() ){
    _gmutex.unlock();
    t_gtime tmp(FIRST_TIME);
    return tmp;
  }
   
  gt.from_str( "%Y-%m-%d %H:%M:%S", trim(str) );
  _gmutex.unlock(); return gt;
}


// Return gtime
// ----------
t_gtime t_gsetgen::end()
{
  _gmutex.lock();
   
  string str = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value( "end" );
  substitute(str,"\n","");   
  substitute(str,"\"","");   

  t_gtime gt(t_gtime::GPS);
   
  if( str.empty() ){
    _gmutex.unlock();
    t_gtime tmp(LAST_TIME);
    return tmp;
  }
      
  gt.from_str( "%Y-%m-%d %H:%M:%S", trim(str) );
  _gmutex.unlock(); return gt;
}


// Return sampling
// ----------
double t_gsetgen::sampling()
{
  _gmutex.lock();
   
  string str  = _doc.child(XMLKEY_ROOT).child(XMLKEY_GEN).child_value( "int" );

  // delete spaces
  str.erase(remove(str.begin(), str.end(), ' '), str.end());
   
  double tmp = str2dbl(str);  

  if( str.find(".") != string::npos ){
    _dec = str.substr(str.find(".")+1).length(); // decimal digits resolution
  }
   
  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetgen::sys()
{
  _gmutex.lock();
   
  set<string> xcl, tmp = t_gsetbase::_setval( XMLKEY_GEN, "sys" );

  // exclude starting with '-'
  set<string>::iterator itSYS, itTMP;
  for( itTMP = tmp.begin(); itTMP != tmp.end(); ){
   if( (*itTMP)[0] == '-' ){ xcl.insert(*itTMP); itSYS=itTMP; ++itTMP; tmp.erase(itSYS); }
   else ++itTMP;
  }

  // if empty, complete, i.e. if only exclusions listed
  if( tmp.size() == 0 ){     
    t_map_sats gnss_sats = GNSS_SATS();
    t_map_sats::const_iterator itGNS;
    // loop over all systems
    for( itGNS = gnss_sats.begin(); itGNS != gnss_sats.end(); ++itGNS ){
      string gs = t_gsys::gsys2str( itGNS->first );
      if( xcl.find("-"+gs) == xcl.end() ) tmp.insert(gs);
    }
  }
  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetgen::rec()
{
  _gmutex.lock();
   
  set<string> tmp = t_gsetbase::_setval( XMLKEY_GEN, "rec" );
  _gmutex.unlock(); return tmp;
}


// settings check
// ----------
void t_gsetgen::check()
{
 _gmutex.lock();

  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_GEN);

  _default_node(node, "BEG", "");                             // all!
  _default_node(node, "END", "");	                      // all!
  _default_node(node, "SYS", "");	                      // all!
  _default_node(node, "REC", "");	                      // none 
  _default_node(node, "INT", int2str(DEF_SAMPLING).c_str());  // default
   
//  xml_node attr_int = _default_node(node, "INT", int2str(DEF_SAMPLING).c_str());  // default
//  _default_attr(attr_int, "unit", bool(DEF_SAMPUNIT));             // default

  // TO CHECK USER-SETUP CONSISTENCY
  _gmutex.unlock();
  if( floor( sampling() ) < 1 || _dec > 0 ){

    if( sampling() < 0.0 ){
      _default_node(node, "INT", "0.0", true); // reset!
      if( _log ) _log->comment(0,"gsetgen","Warning: sampling rate settings negative: reset to 0");
                         cerr << "gsetgen - Warning: sampling rate settings negative: reset to 0\n";
    }
    else{
      if( _log ) _log->comment(1,"gsetgen: sampling rate settings above 1Hz recognized");
    }
  }

  return;
}


// help body
// ----------
void t_gsetgen::help()
{
  _gmutex.lock();
  t_gtime beg(t_gtime::GPS); beg = beg - beg.sod();
  t_gtime end(t_gtime::GPS); end = beg + 86399;

  cerr << "\n <gen>\n"
       << "   <beg> \"" << beg.str_ymdhms() <<"\" </beg>\n"  // FIRST_TIME.str("\"%Y-%m-%d %H:%M:%S\"")
       << "   <end> \"" << end.str_ymdhms() <<"\" </end>\n"  // LAST_TIME.str("\"%Y-%m-%d %H:%M:%S\"")
       << "   <sys> "   << _sys             <<  " </sys>\n"  // GNSS systems
       << "   <rec> GOPE WTZR POTS                </rec>\n"  // list of site identificators
       << "   <int>"+int2str(DEF_SAMPLING)+"</int>\n"
       << " </gen>\n";

  cerr << "\t<!-- general description:\n"
       << "\t beg    .. beg time          (default: all)\n"
       << "\t end    .. end time          (default: all)\n"
       << "\t int    .. data sampling     (default: 30s)\n"
       << "\t sys    .. GNSS system(s)    (default: all)\n"
       << "\t rec    .. GNSS receiver(s)  (rec active list, e.g.: GOPE ONSA WTZR ... )\n"
       << "\t -->\n\n";

   _gmutex.unlock(); return;
}

} // namespace
