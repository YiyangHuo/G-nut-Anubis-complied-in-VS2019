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

#include <set>
#include <iostream>

#include "gio/gxtr.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructors
// ---------
 t_gxtr::t_gxtr()
 : t_gxml("QC_GNSS")
{
}
   
t_gxtr::t_gxtr(t_gsetbase* set,
             const string& pgm, const t_gtime& dt)
 : t_gxml("QC_GNSS"),
   _set(set),
   _gobj(0),
   _gnav(0),
   _nav(0),
   _log(0),
   _pgm(pgm),
   _dt(dt),
   _site(""),
   _smp_req(0.0),
   _nsat(NSAT)     
{
  gtrace("t_gxtr::constructor");
  _pos[0] = _pos[1] = _pos[2] = 0.0;
  gtrace("t_gxtr::t_gxtr");

  _get_settings();      
}

// Destructor
// ---------
t_gxtr::~t_gxtr()
{
  gtrace("t_gxtr::desctructor");
  if( _nav ){ if( _nav->is_open() ) _nav->close(); delete _nav; }
}

// public interface for navig
// ---------------------------
void t_gxtr::navig(string site)
{
  gtrace("t_gxtr::navig");
  if( site.empty() && _site.empty() ) return;
  if( ! site.empty() ) _site = site;

   _setOut();
   
  if( _log ) _log->comment(0, "gxtr", "SITE: " + _site + " " + _nav->mask());

   t_gtime t1,t2;
   t1 = t2 = t_gtime::current_time();

   ostringstream os;
   _navig(os);  t1=t2; t2=t_gtime::current_time(); if(_log) _log->comment(0,"gxtr",_site+": " + dbl2str(t2-t1) );
   
  if( _nav ){
    _nav->write("# " + _pgm + "\n\n"); // pgm version
    _nav->write(os.str());
  }
   
}   

// set data & gnav
// ---------
void t_gxtr::setNAV(t_gallnav* gnav)
{   
  gtrace("t_gxtr::setNAV");
  _gnav = gnav;
}
   
// set object
// ---------
void t_gxtr::setOBJ(t_gallobj* gobj)
{   
  gtrace("t_gxtr::setOBJ");
  _gobj = gobj;
}

// set site
// ---------
void t_gxtr::setSIT(string site)
{   
  gtrace("t_gxtr::setSIT");
  _site = site;
}

// help
// ---------
void t_gxtr::help(){}


// check
// ---------
void t_gxtr::check(){}

   
// PROTECTED
// ------------

// Set log
// -----------------------------
void t_gxtr::_setOut()
{
  gtrace("t_gxtr::_setOut");

  string tmp;
  if( _nav ){ if( _nav->is_open() ) _nav->close(); delete _nav; _nav = 0; }   
  tmp = dynamic_cast<t_gsetout*>(_set)->outputs("nav");
  if( !tmp.empty() ){
    substitute(tmp,"$(rec)",_site, false);
    _nav = new t_glog;
    _nav->tsys(t_gtime::GPS);
    _nav->mask(tmp);
    _nav->append(false);
    _nav->verb(dynamic_cast<t_gsetout*>(_set)->verb());
  }
}   

// generate columns: SITE X Y Z ELE AZI PRN
// -------------------------------
void t_gxtr::_navig(ostringstream& os)
{
  gtrace("t_gxtr::_navig");
   
  shared_ptr<t_grec> setrec = _rec->grec(_site);
  if( setrec == 0 ){
    if( _log ) _log->comment(0,"gxtr","warning[gxtr]: no receiver settings.");
    else       cerr << "*** warning[gxtr]: no receiver in settings.\n";
    return;
  }   

   os << W(20) << "#yyyy-mm-dd hr:mm:ss"
      << W(6)  << "site"
      << W(13) << "rec_X[m]"   << W(13) << "rec_Y[m]"   << W(13) << "rec_Z[m]"
      << setprecision(5)
      << W(10) << "ele[deg]"   << W(10) << "azi[deg]"   << W(6)  << "prn"
      << setprecision(3)
      << W(15) << "sat_X[m]"   << W(15) << "sat_Y[m]"   << W(15) << "sat_Z[m]"
      << endl;
   
  set<string> sats = _gnav->satellites();
  set<string>::const_iterator itPRN;
     
  t_gtime tmp_beg = _set_beg;
  t_gtime tmp_end = _set_end;

   // CHECK IF BEG/END DEFINED IN SETTING
   // -- validity time in classes include gallnav/gallprec (i.e. NAV x SP3 etc)
   // -- IMPROVE beg_time/end_time in gallnav/gallprec classes (i.e. NAV x PREC)
  if( tmp_beg == FIRST_TIME || tmp_end == LAST_TIME ){
    for( itPRN = sats.begin(); itPRN != sats.end(); itPRN++ ){
      string prn = *itPRN;
      if( tmp_beg == FIRST_TIME ) tmp_beg = _gnav->beg_time(prn);
      if( tmp_end == LAST_TIME  ) tmp_end = _gnav->end_time(prn);
    }
  }

  // LOOP OVER SAMPLING RATE
  for(t_gtime epo = tmp_beg; epo < tmp_end; epo = epo + _smp_req){
    t_gtriple crd = setrec->crd(epo);
    
    for(set<string>::iterator itPRN = sats.begin(); itPRN != sats.end(); itPRN++){
      string prn = *itPRN;
      t_gtime nav_beg, nav_end;
      
      /* SHOULD NOT BE NECESSARY SINCE BEING HANDLED IN GALLNAV/GALLPREC
       if(_gnav->id_type() == t_gdata::ALL_NAV){
       nav_beg = _gnav->beg_gnav(prn) - 2*3600;;
       nav_end = _gnav->end_gnav(prn) + 2*3600;
       }else if(_gnav->id_type() == t_gdata::ALLPREC){
       nav_beg = dynamic_cast<t_gallprec*>(_gnav)->beg_data(prn) - 15*60;
       nav_end = dynamic_cast<t_gallprec*>(_gnav)->end_data(prn) + 15*60;
       }
       if(epo < nav_beg || epo > nav_end) continue;
       */
#ifdef DEBUG    
    cout << "prn:" << prn << " epo: " << epo.str_ymdhms() 
                          << " beg: " << nav_beg.str_ymdhms() 
                          << " end: " << nav_end.str_ymdhms() << endl;
#endif

      double xyz[3]  = {0.0, 0.0, 0.0};
      double vel[3]  = {0.0, 0.0, 0.0};
      double var[3]  = {0.0, 0.0, 0.0};
    
      if(_gnav->pos(prn, epo, xyz, var, vel) < 0) continue;
      
      t_gtriple xyz_sat(xyz);
      t_gtriple ell_sit;
      xyz2ell(crd, ell_sit, false);

      t_gtriple xyz_rho = xyz_sat - crd;
      
      double rho = (xyz_sat - crd).norm();
      
      t_gtriple neu_sat;
      xyz2neu(ell_sit, xyz_rho, neu_sat);
      
      double NE2 = neu_sat[0]*neu_sat[0] + neu_sat[1]*neu_sat[1];
      double ele = acos(sqrt(NE2)/rho);
      if( neu_sat[2]<0.0 ) continue;

      ele *= R2D;

      double azi = atan2(neu_sat[1], neu_sat[0]);
      if (azi < 0) azi += 2*G_PI;
      azi *= R2D;
    
    os << W(20) << epo.str_ymdhms(" ")
       << W(6)  << _site    << fixed
       << setprecision(3) 
       << W(13) << crd[0]   << W(13) << crd[1]   << W(13) << crd[2]
       << setprecision(5)  
       << W(10) << ele      << W(10) << azi      << W(6)  << prn
       << setprecision(3) 
       << W(15) << xyz_sat[0] << W(15) << xyz_sat[1] << W(15) << xyz_sat[2]
       << endl;
    }
  }
}
       
// Get settings
// -----------------------------
void t_gxtr::_get_settings()
{
  gtrace("t_gxtr::_get_settings");

  _smp_req          = dynamic_cast<t_gsetgen*>(_set)->sampling();
  _set_beg          = dynamic_cast<t_gsetgen*>(_set)->beg();
  _set_end          = dynamic_cast<t_gsetgen*>(_set)->end();
  _rec              = dynamic_cast<t_gsetrec*>(_set);   

}

   // section title
// ---------
void t_gxtr::_section(ostringstream& os, string tit, int verb)
{
  os << endl 
     << "#====== " << tit << " (v." << verb << ")"  << endl;
}


// subsection
// ---------
void t_gxtr::_subsec(ostringstream& os, string tit)
{
#ifdef DEBUG
   os << "#" << tit << endl;
#endif
}


// legend
// ---------
void t_gxtr::_legend(ostringstream& os, string tit, string tit2, int w)
{
  _setkey(os, tit, '#'); os << W(8) << tit2;

  for(int i = 1; i <= _nsat; ++i )
    os << W(w) << 'x' << W(2) << setfill('0') << i << setfill(' ');

  os << endl;
}


// section title
// ---------
void t_gxtr::_setkey(ostringstream& os, string key, char vrb, string epo)
{
  if( epo == "XXXX-XX-XX XX:XX:XX" ) epo = _ref.str_ymdhms();
   
  os << W(1)          << vrb
     << W(6)  << left << key << " "
     << W(19) << left << epo << right;
}
   
}// namespace
