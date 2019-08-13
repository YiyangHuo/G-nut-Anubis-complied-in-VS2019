
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

#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include "gdata/gdata.h" 
 
using namespace std;

namespace gnut {

// constructor
// ----------
t_gdata::t_gdata()
  : t_gmonit("gdata"),
    _type(NONE),
    _group(GRP_NONE)
{ 
  _log = 0;
}



// copy constructor
// ----------
t_gdata::t_gdata( const t_gdata& data )
  : t_gmonit( data )
{ 
  _log   = data.glog();
  _type  = data.id_type();
  _group = data.id_group();
}


// destructor
// ----------
t_gdata::~t_gdata(){

}


// assignment operator
// ----------
t_gdata& t_gdata::operator=( const t_gdata& data )
{ 
  _log   = data.glog();
  _type  = data.id_type();
  _group = data.id_group();
  return *this;
}


// set data type
// ----------
int t_gdata::id_type( ID_TYPE t )
{
  unsigned int last = LAST;
  for( unsigned int i = 0; i < last ; ++i ){
    if( t == ID_TYPE(i) ){
      return _type = t;
    }
  }
   
  if( _log ){
    string tmp("Warning: Unknown data id. Reset to 0!\n"); // t !
    _log->comment(1,"gdata",tmp);
  }
  _type = t_gdata::NONE;
  _moni_id = str_type();
  return 0;
}


// set group type
// ----------
int t_gdata::id_group( ID_GROUP g )
{
  unsigned int last = GRP_LAST;
  for( unsigned int i = 0; i < last ; ++i ){
    if( g == ID_GROUP(i) ){
      return _group = g;
    }
  }
  if( _log ){
    string tmp("Warning: Unknown group type. Reset to 0!\n"); // g !
    _log->comment(1,"gdata",tmp);
  }
  _group = t_gdata::GRP_NONE;
  return 0;
}

// ID_TYPE to string
// -----------------     
string t_gdata::type2str(ID_TYPE type)
{
   string str = "";
   
  switch( type ){
   case  NONE    :  str = "NONE";    break;
   case  OBS     :  str = "OBS";     break;
   case  OBSGNSS :  str = "OBSGNSS"; break;

   case  QCDATA  :  str = "QCDATA";  break;
   case  QCPROD  :  str = "QCPROD";  break;

   case  EPH     :  str = "EPH";     break;
// case  EPHALL  :  str = "EPHALL";  break; 
   case  EPHGPS  :  str = "EPHGPS";  break;
   case  EPHGLO  :  str = "EPHGLO";  break;
   case  EPHGAL  :  str = "EPHGAL";  break;
   case  EPHQZS  :  str = "EPHQZS";  break;
   case  EPHBDS  :  str = "EPHBDS";  break;
   case  EPHIRN  :  str = "EPHIRN";  break;     
   case  EPHSBS  :  str = "EPHSBS";  break;
   case  EPHPREC :  str = "EPHPREC"; break;
   case  EPHRTCM :  str = "EPHRTCM"; break;
	    
   case  ALLGIO  :  str = "ALLGIO" ; break;
   case  ALLNAV  :  str = "ALLNAV" ; break;
   case  ALLPREC :  str = "ALLPREC"; break;
   case  ALLRTCM :  str = "ALLRTCM"; break;
   case  ALLOBS  :  str = "ALLOBS" ; break;
   case  ALLOBJ  :  str = "ALLOBJ" ; break;
   case  ALLPCV  :  str = "ALLPCV" ; break;
   case  ALLOTL  :  str = "ALLOTL" ; break;
     
   case  STRBUFF :  str = "STRBUFF"; break;
   case  POS     :  str = "POS";     break;
   case  POST    :  str = "POST";    break;    
   case  MET     :  str = "MET";     break;
   case  ION     :  str = "ION";     break;
   case  TRP     :  str = "TRP";     break;
   case  TRPSLT  :  str = "TRPSLT";  break;
	    
   case  PCV     :  str = "PCV";     break;
   case  OTL     :  str = "OTL";     break;
   case  BIAS    :  str = "BIAS";    break;
   case  ERP     :  str = "ERP";     break;
	    
   case  SOL     :  str = "SOL";     break;
   case  LAST    :  str = "UNDEF";   break;
   default       :  str = "UNDEF";
  }
   return str;
}
   
   
// get data type
// ----------
string t_gdata::str_type() const
{
   string type = type2str(_type);   

   return type;
}


// get data type
// ---------- 
string t_gdata::str_group() const
{
  string group;
  switch( _group ){
   case GRP_NONE    : group = "GRP_NONE";    break;
   case GRP_OBSERV  : group = "GRP_OBSERV";  break;
   case GRP_QC      : group = "GRP_QC";      break;
   case GRP_EPHEM   : group = "GRP_EPHEM";   break;
   case GRP_PRODUCT : group = "GRP_PRODUCT"; break;
   case GRP_MODEL   : group = "GRP_MODEL";   break; 
   case GRP_SOLUT   : group = "GRP_SOLUT";   break; 
   case GRP_GRID    : group = "GRP_GRID";    break;
   case GRP_GIO     : group = "GRP_GIO";     break;
   case GRP_LAST    : group = "GRP_UNDEF";   break;
   default          : group = "GRP_UNDEF";
  }
return group;
}



// show function
// ----------
string t_gdata::show(int verb)
{
  gtrace("t_gdata::show");
  ostringstream os("");
   
//  if( (!_log) || ( _log  && _log->verb() >= verb )){
    os << "gdata: nothing to show \n";
//  }
   
//  cout << _moni_id << ":" << endl; 
//  char buff[200 + 1] = "";  
//  sprintf( buff, "%i\n", 2222222222 );
//  sprintf( buff, "%s: %i\n", _moni_id, _data ); // NEFUNGUJE
//  buff[200+1] = '\0';

  return os.str();
}

} // namespace
