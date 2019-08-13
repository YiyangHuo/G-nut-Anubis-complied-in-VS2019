
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

#include "gset/gsetgnss.h"
#include "gutils/gnss.h"
#include "gutils/gobs.h"
#include "gutils/gsys.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetgnss::t_gsetgnss() 
 : t_gsetbase()
{
  _set.insert(XMLKEY_GNSS);
   
//  t_map_gnss::const_iterator itGNSS;
//  t_map_type::const_iterator itTYPE;
//  t_map_attr::const_iterator itBAND;
//  t_vec_attr::const_iterator itATTR;
  stringstream os;

  _sigma_def[GPS] = t_gpair(2.0, 0.02);
  _sigma_def[GLO] = t_gpair(4.0, 0.04);
  _sigma_def[GAL] = t_gpair(3.0, 0.03);
  _sigma_def[BDS] = t_gpair(5.0, 0.03);
  _sigma_def[QZS] = t_gpair(5.0, 0.03);
  _sigma_def[IRN] = t_gpair(5.0, 0.03);

  _maxres_def[GPS] = t_gpair(10.0, 0.08);
  _maxres_def[GLO] = t_gpair(15.0, 0.08);
  _maxres_def[GAL] = t_gpair(15.0, 0.08);
  _maxres_def[BDS] = t_gpair(15.0, 0.08);
  _maxres_def[QZS] = t_gpair(15.0, 0.08);
  _maxres_def[IRN] = t_gpair(15.0, 0.08);
   
  t_map_gnss gnss_data = GNSS_DATA_PRIORITY();
  for( auto itGNSS = gnss_data.begin(); itGNSS != gnss_data.end(); ++itGNSS ) // using C++11 initializer
  {
    GSYS gsys = itGNSS->first;
    string gs = t_gsys::gsys2str( gsys );      
     
    for( auto itBAND  = gnss_data[gsys].begin();
              itBAND != gnss_data[gsys].end();
            ++itBAND )
    {
      GOBSBAND gobsband = itBAND->first;
      string band = gobsband2str( gobsband );
      _band_str[gsys].push_back(     band );
//    _band_obs[gsys].push_back( gobsband );
       
      for( auto itTYPE  = gnss_data[gsys][gobsband].begin();
                itTYPE != gnss_data[gsys][gobsband].end(); 
              ++itTYPE )
      {
        GOBSTYPE gobstype = itTYPE->first;
        string type = gobstype2str( gobstype );
        os << gs << "  band: " << band << "  type: " << type << "  attr:";

        // ONLY THOSE NOT YET INCLUDED
        set<string> type_search( _type_str[gsys].begin(), _type_str[gsys].end() );
        if( type_search.find( type ) == type_search.end() )
          {     
            
            _type_str[gsys].push_back(     type );
            //        _type_obs[gsys].push_back( gobstype );
          }
        
        for( auto itATTR  = gnss_data[gsys][gobsband][gobstype].begin();
             itATTR != gnss_data[gsys][gobsband][gobstype].end(); 
             ++itATTR )
          {       
            GOBSATTR gobsattr = *itATTR;
            string attr = gobsattr2str( gobsattr );
            os << " " << attr;
            
            // ONLY THOSE NOT YET INCLUDED
            set<string> attr_search( _attr_str[gsys].begin(), _attr_str[gsys].end() );
            if( attr_search.find( attr ) == attr_search.end() )
              {      
                _attr_str[gsys].push_back(     attr );
                //          _attr_obs[gsys].push_back( gobsattr );
              }
          }
        os << endl;
      }
    }
    os << endl;
  }

#ifdef DEBUG
  cout << endl << "GNSS DEFAULT SETTINGS:\n" << os.str();
#endif
}


// Destructor
// ----------
t_gsetgnss::~t_gsetgnss()
{}


// Return set (to replace gsetgen general functions, used in processing modules)
// ---------- 
set<string> t_gsetgnss::sat()
{   
  _gmutex.lock();

  set<string> tmp;
  set<string>::const_iterator it;

  t_map_sats gnss_sats = GNSS_SATS();
  t_map_sats::const_iterator itGNS;
  for( itGNS = gnss_sats.begin(); itGNS != gnss_sats.end(); ++itGNS ){
    GSYS gsys = itGNS->first;
    string gs = t_gsys::gsys2str(gsys);
    set<string> sats = t_gsetbase::_setval( gs, "sat" );
   
    // AUTO SET
    if( sats.size() == 0 ) sats = gnss_sats[gsys];
     
    for( it = sats.begin(); it != sats.end(); ++it ) tmp.insert( *it );
     
  }
   
  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetgnss::sat(GSYS gsys, bool def) // default=true enforce to fill 
{
  _gmutex.lock();
  set<string> tmp = t_gsetbase::_setval( _gsys(gsys), "sat" );

  // AUTO SET
  if( def && tmp.size() == 0 ) tmp = GNSS_SATS()[gsys];
   
  _gmutex.unlock(); return tmp;
}


// Return set
// ----------
set<string> t_gsetgnss::nav(GSYS gsys, bool def) // default=true enforce to fill 
{
  _gmutex.lock();  
  set<string> tmp = t_gsetbase::_setval( _gsys(gsys), "nav" );
   
  // AUTO SET
  if( def && tmp.size() == 0 ) tmp = GNSS_GNAV()[gsys];

  _gmutex.unlock(); return tmp;
}



// Return set
// ----------
set<string> t_gsetgnss::obs(GSYS gsys, bool def)
{   
  _gmutex.lock();
  
  string str;
  set<string> tmp;

  vector<string> type = t_gsetbase::_vecval( _gsys(gsys), "type" );
  vector<string> band = t_gsetbase::_vecval( _gsys(gsys), "band" );
  vector<string> attr = t_gsetbase::_vecval( _gsys(gsys), "attr" );
   
  vector<string>::const_iterator itT;
  vector<string>::const_iterator itB;
  vector<string>::const_iterator itA;

  string gs = t_gsys::gsys2str(gsys);
  ostringstream os;
  os << gs << "  band: " << band.size() << "  type: " << type.size() << "  attr: " << attr.size() << endl;

  int bset = band.size();
  int tset = type.size();
  int aset = attr.size();
   
  // AUTO SET
  if( bset == 0 ) band.assign( _band_str[gsys].begin(), _band_str[gsys].end() );
  if( tset == 0 ) type.assign( _type_str[gsys].begin(), _type_str[gsys].end() );
  if( aset == 0 ) attr.assign( _attr_str[gsys].begin(), _attr_str[gsys].end() );
   
  for( itB = band.begin(); itB != band.end(); ++itB ){
    GOBSBAND gband = char2gobsband( (*itB)[0] );
    string       b = gobsband2str( gband );
    if( gband == 999 ) continue;

    for( itT = type.begin(); itT != type.end(); ++itT ){
      GOBSTYPE gtype = char2gobstype( (*itT)[0] );
      string       t = gobstype2str( gtype );
      if( gtype == 999 ) continue;

      os << gs << "  band: " << b << "  type: " << t << "  attr:";

      for( itA = attr.begin(); itA != attr.end(); ++itA ){
        GOBSATTR gattr = char2gobsattr( (*itA)[0] );
        string       a = gobsattr2str( gattr );
        if( gattr == 999 ) continue;
     
        str  = gobstype2str( gtype );
        str += gobsband2str( gband );
        str += gobsattr2str( gattr );

        if( def || (bset || tset || aset) ) tmp.insert( str );
// if( def || (bset && tset && aset) ) tmp.insert( str2gobs( str ) );
        os << " " << a;
    
      }
      os << endl;
    }
  }
#ifdef DEBUG
  cout << endl << "GNSS USER SETTINGS:\n" << os.str();
#endif

  _gmutex.unlock(); return tmp;      
}
   
// Return set
// ----------
set<string> t_gsetgnss::gobs(GSYS gsys)
{   
  _gmutex.lock();
  
  string str;
  set<string> tmp;

  vector<string> vgobs = t_gsetbase::_vecval( _gsys(gsys), "gobs" );
  
  for(auto it = vgobs.begin(); it != vgobs.end(); it++){
    tmp.insert(*it);
  }

#ifdef DEBUG
  cout << "Set GOBS: ";
  for(auto it = tmp.begin(); it != tmp.end(); it++) cout << *it << " ";
  cout << endl;
#endif

  _gmutex.unlock(); return tmp;      
}   


// Return set
// ----------
vector<GOBSTYPE> t_gsetgnss::type(GSYS gsys)
{
  _gmutex.lock();

  vector<GOBSTYPE> v_tmp = _type( gsys );
   
  _gmutex.unlock(); return v_tmp;
}


// Return set
// ----------
vector<GOBSBAND> t_gsetgnss::band(GSYS gsys)
{
  _gmutex.lock();

  vector<GOBSBAND> v_tmp = _band( gsys );

  _gmutex.unlock(); return v_tmp;
}


// Return set
// ----------
vector<GOBSATTR> t_gsetgnss::attr(GSYS gsys)
{
  _gmutex.lock();   

  vector<GOBSATTR> v_tmp = _attr( gsys );

  _gmutex.unlock(); return v_tmp;
}

// Return set
// ----------
double t_gsetgnss::sigma_C(GSYS gsys)
{
  _gmutex.lock();   

  double tmp = _sigma_C( gsys );

  _gmutex.unlock(); return tmp;
}
   
// Return set
// ----------
double t_gsetgnss::sigma_L(GSYS gsys)
{
  _gmutex.lock();   

  double tmp = _sigma_L( gsys );

  _gmutex.unlock(); return tmp;
}
   
// Return set
// ----------
double t_gsetgnss::maxres_C(GSYS gsys)
{
  _gmutex.lock();   

  double tmp = _maxres_C( gsys );

  _gmutex.unlock(); return tmp;
}
   
// Return set
// ----------
double t_gsetgnss::maxres_L(GSYS gsys)
{
  _gmutex.lock();   

  double tmp = _maxres_L( gsys );

  _gmutex.unlock(); return tmp;
}   

// settings check
// ----------
void t_gsetgnss::check()
{
  _gmutex.lock();

/*
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
   
  for( itGNSS = _gnss.begin(); itGNSS != _gnss.end(); ++itGNSS ){

    xml_node node  = _default_node( parent, itGNSS->c_str() );
    
    _default_node( node, "sat"  ); 
    _default_node( node, "type" ); 
    _default_node( node, "band" ); 
    _default_node( node, "attr" );     
  }
*/
  _gmutex.unlock(); return;
}


// convert GNSS
// ----------
string t_gsetgnss::_gsys(GSYS gsys)
{ 
  string tmp = t_gsys::gsys2str(gsys);
  transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
  return tmp;
}



// Return set
// ----------
vector<GOBSTYPE> t_gsetgnss::_type(GSYS gsys)
{ 
  vector<string> v_str = t_gsetbase::_vecval( _gsys(gsys), "type" );
  vector<GOBSTYPE> v_tmp;
  for( vector<string>::const_iterator it = v_str.begin(); it != v_str.end(); ++it ){
    string str = *it; 
    transform(str.begin(), str.end(), str.begin(), ::toupper);
    GOBSTYPE gobstype = str2gobstype( str );
    if( gobstype != TYPE ) v_tmp.push_back(gobstype);
  };
  return v_tmp;
}


// Return set
// ----------
vector<GOBSBAND> t_gsetgnss::_band(GSYS gsys)
{   
  vector<string> v_str = t_gsetbase::_vecval( _gsys(gsys), "band" );
  vector<GOBSBAND> v_tmp;
  for( vector<string>::const_iterator it = v_str.begin(); it != v_str.end(); ++it ){
    string str = *it;
    transform(str.begin(), str.end(), str.begin(), ::toupper);
    GOBSBAND gobsband = str2gobsband( str );
    if( gobsband != BAND ) v_tmp.push_back(gobsband);
  };
  return v_tmp;
}


// Return set
// ----------
vector<GOBSATTR> t_gsetgnss::_attr(GSYS gsys)
{
  vector<string> v_str = t_gsetbase::_vecval( _gsys(gsys), "attr" );
  vector<GOBSATTR> v_tmp;
  for( vector<string>::const_iterator it = v_str.begin(); it != v_str.end(); ++it ){
    string str = *it; 
    transform(str.begin(), str.end(), str.begin(), ::toupper);
    GOBSATTR gobsattr = str2gobsattr( str );
    if( gobsattr != ATTR ) v_tmp.push_back(gobsattr);
  };
  return v_tmp;
}


// Return set
// ----------
double t_gsetgnss::_sigma_C(GSYS gsys)
{
  double sig = 0.0;

  sig = t_gsetbase::_dblatt( _gsys(gsys), "sigma_C" );

  // default value
  if( double_eq(sig, 0.0) ) sig = _sigma_def[gsys][0];  // t_gpair(sig_code, sig_phase)
   
  return sig;
}

// Return set
// ----------
double t_gsetgnss::_sigma_L(GSYS gsys)
{
  double sig = 0.0;
  sig = t_gsetbase::_dblatt( _gsys(gsys), "sigma_L" );

  // default value
  if( double_eq(sig, 0.0) ) sig = _sigma_def[gsys][1];   // t_gpair(sig_code, sig_phase)
   
  return sig;
}
   
// Return set
// ----------
double t_gsetgnss::_maxres_C(GSYS gsys)
{
  double res = 0.0;

  res = t_gsetbase::_dblatt( _gsys(gsys), "maxres_C" );

  // default value
  if( double_eq(res, 0.0) ) res = _maxres_def[gsys][0];  // t_gpair(res_code, res_phase)
   
  return res;
}

// Return set
// ----------
double t_gsetgnss::_maxres_L(GSYS gsys)
{
  double res = 0.0;
   
  res = t_gsetbase::_dblatt( _gsys(gsys), "maxres_L" );

  // default value
  if( double_eq(res, 0.0) ) res = _maxres_def[gsys][1];   // t_gpair(res_code, res_phase)
   
  return res;
}   
   
// help body
// ----------
void t_gsetgnss::help()
{
  _gmutex.lock();

  cerr << " <gps>            \t\t  <!-- any GNSS constellation: GPS GLO GAL BDS SBS QZS -->\n"
       << "   <sat>  </sat>  \t\t  <!-- list of GPS satellites: G01 G02 .. or empty(ALL) -->\n"
       << "   <type> </type> \t\t  <!-- list of GPS  obs types: C L D S P or empty(ALL)  -->\n"
       << "   <band> </band> \t\t  <!-- list of GPS  obs bands: 1 2 5 or empty(ALL) -->\n"
       << "   <attr> </attr> \t\t  <!-- list of PGS attributes: A B C D I L M N P Q S W X Y Z or empty(ALL) -->\n"
       << " </gps>\n\n";

  cerr << " <glo>            \t\t  <!-- any GNSS constellation: GPS GLO GAL BDS SBS QZS -->\n"
       << "   <sat>  </sat>  \t\t  <!-- list of GPS satellites: R01 R02 .. or empty(ALL) -->\n"
       << "   <type> </type> \t\t  <!-- list of GPS  obs types: C L D S P or empty(ALL)  -->\n"
       << "   <band> </band> \t\t  <!-- list of GPS  obs bands: 1 2 3 or empty(ALL) -->\n"
       << "   <attr> </attr> \t\t  <!-- list of PGS attributes: A B C D I L M N P Q S W X Y Z or empty(ALL) -->\n"
       << " </glo>\n\n";

  cerr << "\t<!-- gnss observation definition:\n"
       << "\t      .. ID name\n"
       << "\t name   .. full name\n"
       << "\t desc   .. description\n"
       << "\t -->\n\n";

  _gmutex.unlock(); return;   
}

} // namespace
