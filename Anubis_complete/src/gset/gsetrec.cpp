
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

#include "gset/gsetrec.h"
#include "gutils/gsysconv.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gsetrec::t_gsetrec() 
 : t_gsetbase()
{
  _set.insert(XMLKEY_REC);
  _id         = "";
  _name       = "";
  _desc       = "";
  _domes      = "";
  _X          = 0.0;
  _Y          = 0.0;
  _Z          = 0.0;
  _dX         = 0.0;
  _dY         = 0.0;
  _dZ         = 0.0;
  _dE         = 0.0;
  _dN         = 0.0;
  _dU         = 0.0;   
  _ZTD        = 0.0;
  _ZTD_sig    = 0.0;
  _CRD_sig    = "";
  _rec        = "";
  _ant        = "";
  _beg        = FIRST_TIME;
  _end        = LAST_TIME;
  _overwrite  = false;
}


// Destructor
// ----------
t_gsetrec::~t_gsetrec()
{}


// Return value
// ----------
string t_gsetrec::id()
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("id").value();
   
  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetrec::id(string s)
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str()).attribute("id").value();
   
  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetrec::name(string s)
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str()).attribute("name").value(); 
   
  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetrec::desc(string s)
{
  _gmutex.lock();
   
  string tmp = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str()).attribute("desc").value();
   
  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetrec::rec(string s, string t)
{
  _gmutex.lock();
         
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id",  s.c_str());
  string tmp = site.find_child_by_attribute("set", "beg", t.c_str()).attribute("rec").value(); 
   
  _gmutex.unlock(); return tmp;
}


// Return value
// ----------
string t_gsetrec::ant(string s, string t)
{
  _gmutex.lock();
   
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id",  s.c_str());
  string tmp = site.find_child_by_attribute("set", "beg", t.c_str()).attribute("ant").value();

  _gmutex.unlock(); return tmp;
}

// Return value
// ----------
shared_ptr<t_grec> t_gsetrec::grec(string s, t_glog* glog)
{
  _gmutex.lock();

  string  str;
  shared_ptr<t_grec>  tmp = make_shared<t_grec>();
   
  // return if no rec id found
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str() );
  if( ! site ) {_gmutex.unlock(); return 0;}   
   
  t_gtime begdef(FIRST_TIME);
  t_gtime enddef(LAST_TIME);
  t_gtime beg(begdef), end(enddef);
  t_gtriple xyz = _get_crd_xyz(s);
  t_gtriple blh = _get_crd_blh(s);
  t_gtriple eccneu = _get_ecc_neu(s);
   
#ifdef DEBUG
  t_gtriple XYZ;  ell2xyz(blh.crd_array(),XYZ.crd_array(),true);
  t_gtriple ELL;  xyz2ell(XYZ.crd_array(),ELL.crd_array(),true);

cerr << "REC: " << fixed << setprecision(3)
     << setw(14) << xyz.crd(0)
     << setw(14) << xyz.crd(1)
     << setw(14) << xyz.crd(2) << setprecision(5)
     << setw(14) << blh.crd(0)
     << setw(14) << blh.crd(1)
     << setw(14) << blh.crd(2) 
     << endl;
cerr << "REC: " << fixed << setprecision(3)
     << setw(14) << XYZ.crd(0)
     << setw(14) << XYZ.crd(1)
     << setw(14) << XYZ.crd(2) << setprecision(5)
     << setw(14) << ELL.crd(0)
     << setw(14) << ELL.crd(1)
     << setw(14) << ELL.crd(2)
     << endl;
#endif

  if (double_eq(xyz.crd(0), 0.0) &&
      double_eq(xyz.crd(1), 0.0) &&
      double_eq(xyz.crd(2), 0.0) 
  ){
     if( double_eq(blh.crd(2), HSL_UNKNOWN) ){
       if( _log ) _log->comment(0,"gsetrec","warning - invalid coordinates for rec: " + s );
       else               cerr << "gsetrec:  warning - invalid coordinates for rec: " + s << endl;
     }else{	     
       ell2xyz(blh.crd_array(),xyz.crd_array(),true);
     }
  }
   
  // find first node for rec/site s
  string str_beg_gen = _doc.child(XMLKEY_ROOT).child("gen").child_value( "beg" );
  substitute(str_beg_gen,"\"","");
  t_gtime gen_beg(t_gtime::GPS);   
  if( !str_beg_gen.empty() ) gen_beg.from_str( "%Y-%m-%d %H:%M:%S", trim(str_beg_gen) ); 
  else gen_beg = FIRST_TIME;

  string str_end_gen = _doc.child(XMLKEY_ROOT).child("gen").child_value( "end" );
  substitute(str_end_gen,"\"","");
  t_gtime gen_end(t_gtime::GPS);   
  if( !str_end_gen.empty() ) gen_end.from_str( "%Y-%m-%d %H:%M:%S", trim(str_end_gen) );
  else gen_end = FIRST_TIME;   

  tmp->glog(glog);
  tmp->overwrite( site.attribute("overwrite").as_bool() );
  tmp->id( site.attribute("id").value() );
  tmp->desc( site.attribute("desc").value() );
  tmp->name( site.attribute("name").value() );
  tmp->domes( site.attribute("domes").value() );
  tmp->crd( xyz,                           beg, end); // FULL TIME
  tmp->eccneu( eccneu,                     beg, end); // FULL TIME
  tmp->rec( site.attribute("rec").value(), beg, end); // FULL TIME
  tmp->ant( site.attribute("ant").value(), beg, end); // FULL TIME

#ifdef DEBUG   
  cout << site << " " << gen_beg.str_ymdhms() << " "
                      << gen_end.str_ymdhms() << " " 
                  " " << beg.str_ymdhms() << " "
                      << end.str_ymdhms() << " " 
                      << fixed << setprecision(3)
                      << xyz[0] << " " << xyz[1] << endl;
#endif

  for( xml_node set = site.child("set"); set; set = set.next_sibling("set") ){ // RETURNS REDUNDANT --> problem? !!!!

    str = set.attribute("beg").value(); if(!str.empty()) beg.from_str("%Y-%m-%d %H:%M:%S",str); else beg = begdef;
    str = set.attribute("end").value(); if(!str.empty()) end.from_str("%Y-%m-%d %H:%M:%S",str); else end = enddef;
    
    // default if values are not correct
    if( enddef < end ){ end = enddef; cout << "END VALUE OUT OF INTERVAL !" << endl; }
    if( beg < begdef ){ beg = begdef; cout << "BEG VALUE OUT OF INTERVAL !" << endl; }

    // check BEG+END, gap + WARNINGS ?!  No, simply use the BEG..
    str = set.attribute("rec").value();  if(!str.empty()){ tmp->rec(str,beg,end);}// cout << "SETTING REC: " << str << endl;}
    str = set.attribute("ant").value();  if(!str.empty()){ tmp->ant(str,beg,end);}// cout << "SETTING ANT: " << str << endl;}     

    string x = set.attribute("X").value();
    string y = set.attribute("Y").value();
    string z = set.attribute("Z").value();   
    if( !x.empty() && !y.empty() && !z.empty()) { 
       xyz[0] = str2dbl(x);
       xyz[1] = str2dbl(y);
       xyz[2] = str2dbl(z);
       tmp->crd(xyz,beg,end);
    }

    string dx = set.attribute("DX").value(); 
    string dy = set.attribute("DY").value();
    string dz = set.attribute("DZ").value();         
    if( !dx.empty() && !dy.empty() && !dz.empty()) { 
       xyz[0] = str2dbl(dx);
       xyz[1] = str2dbl(dy);
       xyz[2] = str2dbl(dz);
       tmp->eccxyz(xyz,beg,end);       
    }
     
    string dn = set.attribute("DN").value();
    string de = set.attribute("DE").value();
    string du = set.attribute("DU").value();
    if( !dx.empty() && !dy.empty() && !dz.empty()) { 
       xyz[0] = str2dbl(dn);
       xyz[1] = str2dbl(de);
       xyz[2] = str2dbl(du);
       tmp->eccneu(xyz,beg,end);       
    }     
     
#ifdef DEBUG
    cout << " beg="   << set.attribute("beg").value()
         << " end="   << set.attribute("end").value()
         << " rec="   << set.attribute("rec").value()
         << " ant="   << set.attribute("ant").value()
         << endl;
#endif
  }

#ifdef DEBUG
    t_gtime tt(t_gtime::GPS);  // today time
    tt.from_ymdhms(2010,11,2,0,0,0.0);
    cout << "REC: "   << tmp->id() << setw(12) << "[" + tmp->name() + "]"
         << " tim="   << tt.str("%Y-%m-%d %H:%M:%S[%T]")
         << " rec="   << tmp->rec(tt)
         << " ant="   << tmp->ant(tt)
         << " crd="   << tmp->crd(tt)[0]
         << endl;
#endif

   _gmutex.unlock(); return tmp;
}


// Return value
// ----------
t_gtime t_gsetrec::beg(string s, string t)
{
  _gmutex.lock();
   
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id",  s.c_str());
  string val = site.find_child_by_attribute("set", "beg", t.c_str()).attribute("beg").value();
  if( val.empty() ) { _gmutex.unlock(); return _beg; } // default value

  t_gtime tt(t_gtime::GPS); tt.from_str("%Y-%m-%d %H:%M:%S", val);
  _gmutex.unlock(); return tt;
}


// Return value
// ----------
t_gtime t_gsetrec::end(string s, string t)
{
  _gmutex.lock();
     
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id",  s.c_str());
  string val = site.find_child_by_attribute("set", "beg", t.c_str()).attribute("end").value();
  if( val.empty() ) { _gmutex.unlock(); return _end; }// default value

  t_gtime tt(t_gtime::GPS); tt.from_str("%Y-%m-%d %H:%M:%S", val);
  _gmutex.unlock(); return tt;
}


// Get formats inputs
// ----------
set<string> t_gsetrec::objects()
{
   return _objects();
}

// Get formats inputs
// ----------
set<string> t_gsetrec::_objects()
{
  set<string> tmp;
  string str;

  for( xml_node node = _doc.child(XMLKEY_ROOT).first_child(); node; node = node.next_sibling() ){
    string name = node.name();
    if( name.compare(XMLKEY_REC) == 0 ){
      string str = node.attribute("id").as_string();
      if( ! str.empty() ) tmp.insert( str );
    }
  }
  return tmp;
}

// Return value
// ----------
double t_gsetrec::_aprDX()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("DX").as_double();
}

// Return value
// ----------
double t_gsetrec::_aprDY()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("DY").as_double();
}

// Return value
// ----------
double t_gsetrec::_aprDZ()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("DZ").as_double();
}


// Return value
// ----------
double t_gsetrec::_aprDE()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("DE").as_double();
}

// Return value
// ----------
double t_gsetrec::_aprDN()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("DN").as_double();
}

// Return value
// ----------
double t_gsetrec::_aprDU()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("DU").as_double();
}

// Return value
// ----------
double t_gsetrec::_aprZTD()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("ZTD").as_double();
}

// Return value
// ----------
double t_gsetrec::_sigZTD()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("ZTD_sig").as_double();
}

// Return value
// ----------
string t_gsetrec::_sigCRD()
{
  return _doc.child(XMLKEY_ROOT).child(XMLKEY_REC).attribute("CRD_sig").as_string();
}



// Return value --> IN FUTURE OBSOLETE, BUT STILL USED IN GPPPFILTER (CHANGE TO _get_crd() below )
// ----------
int t_gsetrec::get_crd_xyz(t_gtriple& xyz, string s)
{
  _gmutex.lock();

  xyz = _get_crd_xyz(s);

  if (double_eq(xyz.crd(0), 0.0) &&
      double_eq(xyz.crd(1), 0.0) &&
      double_eq(xyz.crd(2), 0.0)){
     
    _gmutex.unlock(); return 0;
  }
 
  _gmutex.unlock(); return 1;
}

/*
// Return value --> IN FUTURE OBSOLETE, BUT STILL USED IN GPPPFILTER (CHANGE TO _get_crd() below )
// ----------
int t_gsetrec::get_ecc_neu(t_gtriple& neu, string s)
{
  _gmutex.lock();

  neu = _get_ecc_neu(s);

  //if (double_eq(neu.crd(0), 0.0) &&
  //    double_eq(neu.crd(1), 0.0) &&
  //    double_eq(neu.crd(2), 0.0)){
  //   
  //  _gmutex.unlock(); return 0;
  //}
 
  _gmutex.unlock(); return 1;
}

// Return value --> IN FUTURE OBSOLETE, BUT STILL USED IN GPPPFILTER (CHANGE TO _get_crd() below )
// ----------
int t_gsetrec::get_crd_blh(t_gtriple& blh, string s)
{
  _gmutex.lock();

  blh = _get_crd_blh(s);

  if (double_eq(blh.crd(0), 0.0) &&
      double_eq(blh.crd(1), 0.0) &&
      double_eq(blh.crd(2), HSL_UNKNOWN)){

    _gmutex.unlock(); return 0;     
  }

  _gmutex.unlock(); return 1;
}
*/


// Return value
// ----------
t_gtriple t_gsetrec::_get_crd_xyz(string s)
{
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str() );
   
  t_gtriple xyz(site.attribute("X").as_double(),
	        site.attribute("Y").as_double(), 
	 	site.attribute("Z").as_double());
   
  return xyz;
}


// Return value
// ----------
t_gtriple t_gsetrec::_get_ecc_neu(string s)
{
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str() );
   
  t_gtriple neu(site.attribute("DN").as_double(),
	        site.attribute("DE").as_double(), 
	 	site.attribute("DU").as_double());
  return neu;
}
   
// Return value
// ----------
t_gtriple t_gsetrec::_get_ecc_xyz(string s)
{
  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str() );
   
  t_gtriple xyz(site.attribute("DX").as_double(),
                site.attribute("DY").as_double(), 
                site.attribute("DZ").as_double());
  return xyz;
}   

// Return value
// ----------
t_gtriple t_gsetrec::_get_crd_blh(string s)
{
  t_gtriple zero(0.0,0.0,HSL_UNKNOWN);

  xml_node site = _doc.child(XMLKEY_ROOT).find_child_by_attribute(XMLKEY_REC, "id", s.c_str() );
  xml_attribute attr;

  if( (attr = site.attribute("LAT")) &&
      (attr = site.attribute("LON")) &&
      (attr = site.attribute("HSL")) 
  ){     
    t_gtriple blh(site.attribute("LAT").as_double(),
 	          site.attribute("LON").as_double(), 
	   	  site.attribute("HSL").as_double());
    double pres, temp, undu = 0.0;
    _ggpt.gpt_v1(51544.0, blh[0]*D2R, blh[1]*D2R, blh[2], pres, temp, undu); // RAD
    blh[2] += undu; // ADD UNDULATION BECAUSE ONLY HSL (ABOVE MEAN SEE LEVEL) SUPPORTED
    return blh;
  }

  t_gtriple xyz = _get_crd_xyz(s);

  if (!double_eq(xyz.crd(0), 0.0) &&
      !double_eq(xyz.crd(1), 0.0) &&
      !double_eq(xyz.crd(2), 0.0)
  ){
    t_gtriple blh;
    if( xyz2ell(xyz,blh,true) == 1 ) return blh;
    else                             return zero;
  }

  return zero;
}

// settings check
// ----------
void t_gsetrec::check()
{
  _gmutex.lock();
   
  // check existence of nodes/attributes
  xml_node parent = _doc.child(XMLKEY_ROOT);
  xml_node node = _default_node(parent, XMLKEY_REC);

  // check existence of attributes
  _default_attr(node,"X", _X);
  _default_attr(node,"Y", _Y);
  _default_attr(node,"Z", _Z);
  _default_attr(node,"DX", _dX);
  _default_attr(node,"DY", _dY);
  _default_attr(node,"DZ", _dZ);
  _default_attr(node,"DN", _dN);
  _default_attr(node,"DE", _dE);
  _default_attr(node,"DU", _dU);
  _default_attr(node,"overwrite", _overwrite);
   
  _default_attr(node,"ZTD", _ZTD);
  _default_attr(node,"ZTD_sig", _ZTD_sig);
  _default_attr(node,"CRD_sig", _CRD_sig);

  // CHECK USER-SETUP CONSISTENCY

  // check duplicities
  map<string, xml_node> mchild;
  for( xml_node node = _doc.child(XMLKEY_ROOT).child(XMLKEY_REC); node;   ){
    xml_node rec = node;
    node = node.next_sibling(XMLKEY_REC); // get next immediately
     
    // check ID & NAME & DOMES
    string id = rec.attribute("id").value();
    string name = rec.attribute("name").value();
    string domes = rec.attribute("domes").value();

    if( name.empty() ){  // complete NAME
      rec.insert_attribute_after("name", rec.attribute("id")) = id.c_str();
    }
    if( domes.empty() ){  // complete NAME
      if( name.length() < 14  ) { domes = "XXXXXXXXX"; }
      if( name.length() >= 14 ) { domes = tail(name,9); name = name.substr(0, name.size() - 10);}
      rec.insert_attribute_after("domes", rec.attribute("name")) = domes.c_str();
      rec.insert_attribute_after("name", rec.attribute("id")) = name.c_str();
    }

    if( mchild.find(id) == mchild.end() ){
      mchild[id] = rec;
    }else{
      _doc.child(XMLKEY_ROOT).remove_child(rec);
      if( _log ) _log->comment(0,"gsetrec","warning - removed duplicated record [id] for rec:" + id );
      else               cerr << "gsetrec:  warning - removed duplicated record [id] for rec:" + id << endl;
    }

    // check 'set' duplicity
    map<string, xml_node> mattr;
    for( xml_node set = rec.child("set"); set;    ){
      xml_node tmp = set;
      set = set.next_sibling("set"); // get next immediately       

      string beg = tmp.attribute("beg").value();       
      if( mattr.find(beg) == mattr.end() ){
        mattr[beg] = set;
      }else{
        rec.remove_child(tmp);
        if( _log ) _log->comment(0,"gsetrec","warning - removed duplicated record [set] for rec:" + id + ", beg:" + beg );
        else               cerr << "gsetrec:  warning - removed duplicated record [set] for rec:" + id + ", beg:" + beg << endl;
      }
    }
  }

  _gmutex.unlock(); return;
}


// help body
// ----------
void t_gsetrec::help()
{
  _gmutex.lock();

  cerr << " <rec"
       <<        " id=\"GOPE\""
       <<        " name=\"GOPE\""
       <<        " domes=\"11502M002\""
       <<        " desc=\"Geodetic Observatory Pecny, Czech Republic\" >\n"
       << "  <set"
       <<        " beg=\"1995 05 13 00 00 00\""
       <<        " end=\"1997 06 11 00 00 00\""
       <<        " rec=\"TRIMBLE 4000SSE\""
       <<        " ant=\"TRM14532.00     NONE\""
       << "\n\t"
       <<        " X=\"3979316.0\""
       <<        " Y=\"1050312.0\""
       <<        " Z=\"4857067.0\""
       <<        " dX=\"0.0\""
       <<        " dY=\"0.0\""
       <<        " dZ=\"0.0\""
       <<        " dN=\"0.0\""
       <<        " dE=\"0.0\""
       <<        " dU=\"0.0\""
       <<  "  />\n"
       << "  <set"
       <<        " beg=\"1997 06 11 00 00 00\""
       <<        " end=\"1997 06 20 00 00 00\""
       <<        " rec=\"SPP GEOTRACER100\""
       <<        " ant=\"TRM14532.00     NONE\""
       <<  "  />\n"
       << "  <set"
       <<        " beg=\"1997 06 20 00 00 00\""
       <<        " end=\"1999 11 04 00 00 00\""
       <<        " rec=\"TRIMBLE 4000SSE\""
       <<        " ant=\"TRM14532.00     NONE\""
       <<  "  />\n"
       << "  <set"
       <<        " beg=\"1999 11 04 00 00 00\""
       <<        " end=\"2000 07 24 00 00 00\""
       <<        " rec=\"ASHTECH Z18\""
       <<        " ant=\"ASH701073.3     SNOW\""
       <<  "  />\n"
       << "  <set"
       <<        " beg=\"2000 07 24 00 00 00\""
       <<        " end=\"2000 10 04 00 00 00\""
       <<        " rec=\"TRIMBLE 4000SSE\""
       <<        " ant=\"TRM14532.00     NONE\""
       <<  "  />\n"
       << "  <set"
       <<        " beg=\"2000 10 04 00 00 00\""
       <<        " end=\"2001 07 18 00 00 00\""
       <<        " rec=\"ASHTECH Z18\""
       <<        " ant=\"ASH701073.3     SNOW\""
       <<  "  />\n"
       << "  <set"
       <<        " beg=\"2006 07 14 00 00 00\""
       <<        " end=\"2009 12 14 00 00 00\""
       <<        " rec=\"ASHTECH Z18\""
       <<        " ant=\"TPSCR3_GGD      CONE\""
       <<  "  />\n"
       << "  <set"
       <<        " beg=\"2009 12 14 00 00 00\""
       <<        " end=\"2013 03 19 00 00 00\""
       <<        " rec=\"TPS NETG3\""
       <<        " ant=\"TPSCR.G3        TPSH\""
       <<  "  />\n"
       <<  " </rec>\n";
   
  cerr << "\t<!-- receiver description:\n"
       << "\t id     .. ID name\n"
       << "\t name   .. site name\n"
       << "\t domes  .. site domes\n"
       << "\t desc   .. description\n"
       << "\t X[YZ]  .. X,Y,Z-coordinate [m]\n"
       << "\t dX[YZ] .. X,Y,Z-eccentricity [m]\n"
       << "\t -->\n\n";

  _gmutex.unlock(); return;   
}

} // namespace
