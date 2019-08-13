
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

#include "gset/gsetinp.h"
#include "gutils/gfileconv.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Convertor for INP formats
// ----------
IFMT t_gsetinp::str2ifmt(const string& s)
{ 
  string tmp = s;
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
  if( tmp == "RINEXC"  ) return RINEXC_INP;
  if( tmp == "RINEXC_PRD" ) return RINEXC_PRD;
  if( tmp == "RINEXC_REF" ) return RINEXC_REF;   
  if( tmp == "RINEXO"  ) return RINEXO_INP;
  if( tmp == "RINEXN"  ) return RINEXN_INP;
  if( tmp == "RINEXN_PRD" ) return RINEXN_PRD;
  if( tmp == "RINEXN_REF" ) return RINEXN_REF;   
  if( tmp == "RINEXM"  ) return RINEXM_INP;
  if( tmp == "BNCOBS"  ) return BNCOBS_INP;
  if( tmp == "BNCRTCM" ) return BNCRTCM_INP;
  if( tmp == "BNCRTCM_PRD" ) return BNCRTCM_PRD;
  if( tmp == "BNCRTCM_REF" ) return BNCRTCM_REF;   
  if( tmp == "TROSINEX") return TROSINEX_INP;
  if( tmp == "TROSINEX0") return TROSINEX0_INP;
  if( tmp == "TROGRID" ) return TROGRID_INP;
  if( tmp == "RTCM"    ) return RTCM_INP;
  if( tmp == "TSMET"   ) return TSMET_INP;
  if( tmp == "SP3"     ) return SP3_INP;
  if( tmp == "SP3_PRD" ) return SP3_PRD;
  if( tmp == "SP3_REF" ) return SP3_REF;
   
  if( tmp == "ATX"     ) return ATX_INP;
  if (tmp == "BLQ") return BLQ_INP;
  if (tmp == "FCB") return FCB_INP;

  if( tmp == "GRIB_CHMI" ) return GRIB_CHMI_INP;
  if( tmp == "GRIB_HARM" ) return GRIB_HARM_INP;
  if( tmp == "GRIB_ERA"  ) return GRIB_ERA_INP;
  if( tmp == "NCD3_ICS"  ) return NCD3_ICS_INP;   

  if( tmp == "RAO_TEXT"  ) return RAO_TEXT_INP;
  if( tmp == "RAO_IGRA"  ) return RAO_IGRA_INP;
  if( tmp == "RAO_BADC"  ) return RAO_BADC_INP;

  if( tmp == "BIASINEX"  ) return BIASINEX_INP;
  if( tmp == "BIABERN"  ) return BIABERN_INP;  
    
  if( tmp == "IONEX"   ) return IONEX_INP;    

  return IFMT(-1);
}

// Convertor for INP formats
// ----------
string t_gsetinp::ifmt2str(const IFMT& f)
{ 
  switch ( f ){        
   case RINEXO_INP:  return "RINEXO";
   case RINEXC_INP:  return "RINEXC";
   case RINEXC_PRD:  return "RINEXC_PRD";
   case RINEXC_REF:  return "RINEXC_REF";     
   case RINEXN_INP:  return "RINEXN";
   case RINEXN_PRD:  return "RINEXN_PRD";
   case RINEXN_REF:  return "RINEXN_REF";     
   case RINEXM_INP:  return "RINEXM";
   case BNCOBS_INP:  return "BNCOBS";
   case BNCRTCM_INP: return "BNCRTCM";
   case BNCRTCM_PRD: return "BNCRTCM_PRD";
   case BNCRTCM_REF: return "BNCRTCM_REF";     
   case TROSINEX_INP:return "TROSINEX";
   case TROSINEX0_INP:return "TROSINEX0";
   case TROGRID_INP: return "TROGRID";
   case RTCM_INP:    return "RTCM";
   case TSMET_INP:   return "TSMET";
   case SP3_INP:     return "SP3";
   case SP3_PRD:     return "SP3_PRD";
   case SP3_REF:     return "SP3_REF";     
   case ATX_INP:     return "ATX";
   case BLQ_INP:     return "BLQ";
   case FCB_INP:     return "FCB";
   case GRIB_CHMI_INP:  return "GRIB_CHMI";
   case GRIB_HARM_INP:  return "GRIB_HARM";
   case GRIB_ERA_INP:   return "GRIB_ERA";
   case NCD3_ICS_INP:   return "NCD3_ICS";
   case RAO_TEXT_INP:   return "RAO_TEXT";
   case RAO_IGRA_INP:   return "RAO_IGRA";
   case RAO_BADC_INP:   return "RAO_BADC";
   case BIASINEX_INP:   return "BIASINEX";
   case BIABERN_INP:    return "BIABERN";    
   case IONEX_INP:      return "IONEX";    

   default:             return "UNDEF";
  }
  return "UNDEF";
}

// Constructor
// ----------
t_gsetinp::t_gsetinp()
  : t_gsetbase()
{
  _set.insert(XMLKEY_INP);
  _chkNavig  = true;
  _chkHealth = true;
  _corrStream = "";
}


// Destructor
// ----------
t_gsetinp::~t_gsetinp()
{}


// Get formats input size
// ----------
int t_gsetinp::input_size(const string& fmt)
{
  _gmutex.lock();
   
  int tmp = _inputs(fmt).size();

  _gmutex.unlock(); return tmp;
}


// Get formats inputs (all in multimap)
// ----------
multimap<IFMT,string> t_gsetinp::inputs_all()
{
  _gmutex.lock();
   
  multimap<IFMT,string> map;

  set<string> ifmt = _iformats();
  set<string>::const_iterator itFMT = ifmt.begin();
  while( itFMT != ifmt.end() )
  {
    string fmt = *itFMT;
    IFMT  ifmt = str2ifmt( fmt );
    vector<string> inputs = _inputs( fmt );
    vector<string>::const_iterator itINP = inputs.begin();
    while( itINP != inputs.end() ){
      map.insert(map.end(), pair<IFMT,string>( ifmt, *itINP ));
      itINP++;
    }
    itFMT++;
  }
  _gmutex.unlock(); return map;
}


// Get formats inputs
// ----------
vector<string> t_gsetinp::inputs(const string& fmt)
{
   return _inputs(fmt);
}


// Get input formats
// ----------
set<string> t_gsetinp::iformats()
{
   return _iformats();  
}

// Get correction stream
// ----------
string t_gsetinp::corrStream()
{
  _gmutex.lock();

  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).child("bncrtcm").attribute("stream").value();

  _gmutex.unlock(); return tmp;
}

  
// Get formats inputs
// ----------
vector<string> t_gsetinp::_inputs(const string& fmt)
{
  vector<string> tmp;
  set<string> list;
  string str;

  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling() ){
    if( node.name() == fmt ){       
      istringstream is( node.child_value() );
      while( is >> str && !is.fail() ){
        if( str.find("://") == string::npos ) str = GFILE_PREFIX + str;
        if( list.find(str) == list.end() ){
          tmp.push_back( str );
          list.insert( str );
        }else{
          if( _log ) _log->comment(1,"gsetinp","READ: "+str+" multiple request ignored");
        }
      }
    }
  }
  return tmp;
}

   
// Get input formats
// ----------
set<string> t_gsetinp::_iformats()
{
  set<string> tmp;
  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).first_child(); node; node = node.next_sibling() ){
    tmp.insert( node.name() );
  }
  return tmp;
}

// Checking navigation Rinex
// ----------
bool t_gsetinp::chkNavig()
{
  bool tmp;

  tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).attribute("chk_nav").as_bool();

  return tmp;
}
	 
// Checking satellite healthy status
// ----------
bool t_gsetinp::chkHealth()
{
   bool tmp;

   tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).attribute("chk_health").as_bool();

	 return tmp;
}	 

// settings check
// ----------
void t_gsetinp::check()
{
  _gmutex.lock();
   
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_INP);

  // check supported input formats (see IFMT enum !)
  set<string> ifmt = _iformats();
  set<string>::const_iterator itFMT = ifmt.begin();
  while( itFMT != ifmt.end() ){
    string fmt = *itFMT;
    IFMT  ifmt = str2ifmt( fmt );
    if( ifmt < 0 ){
      _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child( node.child(fmt.c_str()) );
      if( _log ) _log->comment(0, "Warning: " + fmt + " inp format not implemented [gsetinp::check()]!");
      itFMT++;
      continue;
    }
  
    // check application-specific output format
    if( _IFMT_supported.find( ifmt ) == _IFMT_supported.end() ){
      _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).remove_child( node.child(fmt.c_str()) );
      if( _log ) _log->comment(0, "Warning: " + fmt + " inp format not supported by this application!");
      else                cerr << "Warning: " + fmt + " inp format not supported by this application!\n";       
    }
    itFMT++;            
  }
   
   _default_attr(node,"chk_nav",    _chkNavig);
   _default_attr(node,"chk_health", _chkHealth);	 

   xml_node nodeBNCRTCM = _doc.child(XMLKEY_ROOT).child(XMLKEY_INP).child("bncrtcm");
   _default_attr(nodeBNCRTCM,"_corrStream", _corrStream);
     
   _gmutex.unlock(); return;
}


// settings help
// ----------
void t_gsetinp::help()
{
  _gmutex.lock();
   
  
  cerr << " <inputs>\n"
       << "   <rinexo> file://dir/name </rinexo> \t\t <!-- obs RINEX decoder -->\n"
       << "   <rinexn> file://dir/name </rinexn> \t\t <!-- nav RINEX decoder -->\n"
       << " </inputs>\n";

  cerr << "\t<!-- inputs description:\n"
       << "\t <decoder> path1 path2 path3  </decoder>\n"
       << "\t ... \n"
       << "\t where path(i) contains [file,tcp,ntrip]:// depending on the application\n"
       << "\t -->\n\n";

   _gmutex.unlock(); return;  
}

} // namespace
