
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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "gdata/gsatdata.h"
#include "gmodels/gephplan.h"
#include "gutils/gtypeconv.h"
#include "gutils/gsysconv.h"
#include "gutils/gmatrixconv.h"

using namespace std;

namespace gnut {

// -----------
t_gsatdata::t_gsatdata(string site, string sat, const t_gtime& t)
  : t_gobsgnss(site,sat,t),
    _satcrd(0.0,0.0,0.0),
  	_satpco(0.0,0.0,0.0),
    _clk(0.0),
    _ele(0.0),
    _azi(0.0),
    _rho(0.0),
    _eclipse(false),
    _mfH(0.0),
    _mfW(0.0),
    _mfG(0.0),
    _wind(0.0),     
    _low_prec(false),
    _slipf(false),
    _beta_val(999),
    _orb_angle_val(999), 
    _yaw(999)

{ 
  id_type(t_gdata::SATDATA);
  id_group(t_gdata::GRP_OBSERV);
}


// -----------
t_gsatdata::t_gsatdata(const t_gobsgnss& obs)
  : t_gobsgnss(obs), // t_gobsgnss(obs.site() ,obs.sat(), obs.epoch() ),
    _satcrd(0.0,0.0,0.0),
	  _satpco(0.0,0.0,0.0),
    _clk(0.0),
    _ele(0.0),
    _azi(0.0),
    _rho(0.0),
    _eclipse(false),
    _mfH(0.0),
    _mfW(0.0),
    _mfG(0.0),
    _wind(0.0),
    _low_prec(false),
    _slipf(false),
    _beta_val(999),
    _orb_angle_val(999),
		_yaw(999)
{        
  id_type(t_gdata::SATDATA);
  id_group(t_gdata::GRP_OBSERV);
}


// -----------
t_gsatdata::~t_gsatdata()
{
#ifdef DEBUG
  cout << "GSATDATA - destruct POCAT: " << site() << " " << sat() << " "
       << epoch().str("  %Y-%m-%d %H:%M:%S[%T] ") << fixed << setprecision(3) << endl;
   
  vector<GOBS> v_obs = this->obs();
  vector<GOBS>::iterator itOBS = v_obs.begin();
  for( ; itOBS != v_obs.end(); ++itOBS ) 
      cout << " " << t_gobs::gobs2str( *itOBS ) << ":" << this->getobs( *itOBS );
  cout << endl;

  cout << "GSATDATA - destruct KONEC: " << site() << " " << sat() << " "
       << epoch().str("  %Y-%m-%d %H:%M:%S[%T] ");
  cout.flush();   
#endif   
}

// -----------
void t_gsatdata::addpco( const t_gtriple& pco )
{
  gtrace("t_gsatdat::addpco");
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _satpco = pco;
  _gmutex.unlock();  
  return;
}

// -----------
void t_gsatdata::addcrd( const t_gtriple& crd )
{
  gtrace("t_gsatdat::addcrd");
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _satcrd = crd;
  _gmutex.unlock();  
  return;
}

// -----------
void t_gsatdata::addvel( const t_gtriple& vel )
{
  gtrace("t_gsatdat::addvel");
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _satvel = vel;
  _gmutex.unlock();  
  return;
}   
   
//----------
int t_gsatdata::addprd( t_gallnav* gnav, bool corrTOT, bool msk_health )
{
   gtrace("t_gsatdat::addprd");

   _gmutex.lock();
   
   _low_prec = false;
  
   int irc = this->_addprd(gnav, corrTOT, msk_health);

   _gmutex.unlock(); return irc;
}
   
//----------
int t_gsatdata::addprd_nav( t_gallnav* gnav, bool corrTOT, bool msk_health )
{
  gtrace("t_gsatdat::addprd_nav");

   _gmutex.lock();
   
   _low_prec = true;
  
   int irc = this->_addprd(gnav, corrTOT, msk_health);
  
   _gmutex.unlock(); return irc;
}
   

// Compute rho, ele, azi ...
// -----------------------
int t_gsatdata::cmpVal(t_gtriple& xyz)
{
   gtrace("t_gsatdata::cmpVal"); 
   
   if( double_eq(_satcrd[0], 0.0) ||
       double_eq(_satcrd[0], 0.0) ||
       double_eq(_satcrd[0], 0.0) ) {
     cerr << "Satellite position has not been calculated" << endl;
     return -1;
   }
   
   
   if( double_eq(xyz[0], 0.0) ||
       double_eq(xyz[0], 0.0) ||
       double_eq(xyz[0], 0.0) ) {
     cerr << "Station position has not been calculated" << endl;   
     return -1;
   }

   t_gtriple neu_sat;
   t_gtriple ell_sit;
   xyz2ell(xyz, ell_sit, false);
   t_gtriple xyz_rho = _satcrd - xyz;

   xyz2neu(ell_sit, xyz_rho, neu_sat);
      
   // Correct Earth rotation
   ColumnVector xRec(3);
   double rho0 = (_satcrd.crd_cvect() - xyz.crd_cvect()).norm_Frobenius();
   double dPhi = OMEGA * rho0 / CLIGHT;

   xRec(1) = xyz[0] * cos(dPhi) - xyz[1] * sin(dPhi);
   xRec(2) = xyz[1] * cos(dPhi) + xyz[0] * sin(dPhi);
   xRec(3) = xyz[2];  
		 
   double tmp = (_satcrd.crd_cvect() - xRec).norm_Frobenius();
   _rho = tmp;

   double NE2 = neu_sat[0]*neu_sat[0] + neu_sat[1]*neu_sat[1];
   double ele = acos(sqrt(NE2)/_rho);
   if( neu_sat[2]<0.0 ) {
      ele *= - 1.0;
   }

   if( sqrt(NE2)/_rho > 1.0 )
        _ele = 0.0;
   else _ele = ele;

   double azi = atan2(neu_sat[1], neu_sat[0]);
   if (azi < 0) azi += 2*G_PI;
      _azi = azi;

#ifdef DEBUG
   cout << site() << " " << xyz << endl
        << sat() << " " << _satcrd << endl
        << "ele: " << ele*R2D << " azi: " << azi*R2D << endl << endl;
#endif
   
   return 1;
}



// -----------
void t_gsatdata::addclk( double clk )
{
  gtrace("t_gsatdat::addclk");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _clk = clk;
  _gmutex.unlock();
}


// -----------
void t_gsatdata::addele( double ele )
{
  gtrace("t_gsatdat::addele");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _ele = ele;
  _gmutex.unlock();  
}


// -----------
void t_gsatdata::addazi( double azi )
{
  gtrace("t_gsatdat::addazi");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _azi = azi;
  _gmutex.unlock();  
}


// -----------
void t_gsatdata::addrho( double rho )
{
  gtrace("t_gsatdat::addrho");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _rho = rho;
  _gmutex.unlock();  
}
   
// -----------
void t_gsatdata::addmfH( const double& mfH )
{
  gtrace("t_gsatdat::addmfH");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _mfH = mfH;
  _gmutex.unlock();  
}
   
// -----------
void t_gsatdata::addmfW( const double& mfW )
{
  gtrace("t_gsatdat::addmfW");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _mfW = mfW;
  _gmutex.unlock();  
}   

// -----------
void t_gsatdata::addmfG( const double& mfG )
{
  gtrace("t_gsatdat::addmfG");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  _mfG = mfG;
  _gmutex.unlock();  
}   
   

// -----------
t_gtriple t_gsatdata::satcrd()
{
  gtrace("t_gsatdat::satcrd");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
//  t_gtriple crd(_satcrd[0], _satcrd[1], _satcrd[2] );
 
   _gmutex.unlock();  
  return _satcrd;
}

// -----------
t_gtriple t_gsatdata::satpco()
{
  gtrace("t_gsatdat::satcrd");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
 
   _gmutex.unlock();  
  return _satpco;
}

// -----------
t_gtriple t_gsatdata::satvel()
{
  gtrace("t_gsatdat::satvel");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
 
   _gmutex.unlock();  
  return _satvel;
}   

// -----------
double t_gsatdata::clk()
{
  gtrace("t_gsatdat::clk");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _clk;
  _gmutex.unlock();  
  return tmp;
}


// -----------
double t_gsatdata::ele()
{
  gtrace("t_gsatdat::ele");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _ele;
  _gmutex.unlock(); 
  return tmp;
}

// -----------
double t_gsatdata::ele_deg()
{
  gtrace("t_gsatdat::ele_deg");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _ele*180.0/G_PI;
  _gmutex.unlock();
  return tmp;
}


// -----------
double t_gsatdata::azi()
{
  gtrace("t_gsatdat::azi");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _azi;
  _gmutex.unlock();
  return tmp;
}


// -----------
double t_gsatdata::rho()
{
  gtrace("t_gsatdat::rho");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  double tmp = _rho;
  _gmutex.unlock();
  return tmp;
}


// valid
// ----------
bool t_gsatdata::valid()
{
  gtrace("t_gsatdat::valid");
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  bool tmp = this->_valid();
  _gmutex.unlock();
  return tmp;
}

// -----------
void t_gsatdata::addslip(const bool& flag)
{
	gtrace("t_gsatdat::addslip");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();
	_slipf = flag;
	_gmutex.unlock();
}

// -----------
bool t_gsatdata::islip()
{
	gtrace("t_gsatdat::isSlip");

#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif
	_gmutex.lock();
	bool tmp = _slipf;
	_gmutex.unlock();
	return tmp;
}

double t_gsatdata::beta()
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif

   _gmutex.lock();
   double tmp = _beta();
   _gmutex.unlock();
   return tmp;   
}

double t_gsatdata::orb_angle()
{
#ifdef BMUTEX
	boost::mutex::scoped_lock lock(_mutex);
#endif

   _gmutex.lock();
   double tmp = _orb_angle();
   _gmutex.unlock();
   return tmp;   
}   
   
// clean data
// ----------
void t_gsatdata::clear()
{ 
  gtrace("t_gsatdat::clear");   
   
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();  
  this->_clear();
  _gmutex.unlock();
   return;
}

//----------
int t_gsatdata::_addprd( t_gallnav* gnav, bool corrTOT, bool msk_health )
{
  gtrace("t_gsatdat::_addprd");

  string satname( _satid );   
  
  GSYS gs = this->gsys();   
    
  GOBSBAND b1, b2;
  b1 = b2 = BAND;
  
  // automatic selection of two bands for IF LC
  set<GOBSBAND> bands = _band_avail_code();
  auto itBAND = bands.begin();
  if(corrTOT){
    if(bands.size() < 2) {
      cout << "t_gsatdata bands.size() < 2" << endl;
      if(_log) _log->comment(2, "gsatdata", "At least two bands are necessary for TOT correction in sat pos/clk calculation!");
      return -1;
    }
    b1 = *itBAND;
    itBAND++;
    b2 = *itBAND;
  }
  
  _gmutex.unlock();
  double P3 = this->P3(b1, b2);
  _gmutex.lock();   
  
  //test for observations availability
  if( gnav == 0 ) {
    if( _log )
    _log->comment(2, "gsatdata"," satellite " + satname 
                  + _epoch.str_ymdhms("  t_gallnav pointer is not available "));
    //      cout << "Add prd " << epoch().str_hms() << " " << sat() << " t_gallnav pointer is not available" << endl;
    return -1;
  }        
  
  //test for observations availability
  if( double_eq(P3, 0.0) && corrTOT) {
    if( _log )
    _log->comment(2, "gsatdata"," satellite " + satname 
                  + _epoch.str_ymdhms(" P3 = 0!"));
    //                   cout << "Add prd " << epoch().str_hms() << " " << sat() << " P3 = 0;" << endl;
    return -1;      
  }     
  
  double xyz[3]  = {0.0, 0.0, 0.0};
  double vel[3]  = {0.0, 0.0, 0.0};
  double var[3]  = {0.0, 0.0, 0.0};
  double clk =    0.0;
  double dclk =   0.0;
  double clkrms = 0.0;   
  
  if( satname.substr(0,1) != "G" &&
      satname.substr(0,1) != "R" &&
      satname.substr(0,1) != "E" &&
      satname.substr(0,1) != "J" &&
//       satname.substr(0,1) != "S" &&       
      satname.substr(0,1) != "C" )
    {
      if( _log ) {
        _log->comment(2, "gsatdata"," satelite " + satname 
                      + _epoch.str_ymdhms(" Undefined satellite system! "));
      }
      //      cout << "Undefined satellite system" << endl;
      return -1; 
    }
  
  t_gtime epoT(t_gtime::GPS);
  double satclk  = 0.0;
  double satclk2 = 1.0;
  int cnt = 0;
   
  if(corrTOT){
    while( fabs(satclk-satclk2) > 1.e-3/CLIGHT ){
      satclk2 = satclk;
      epoT = _epoch - P3/CLIGHT - satclk;
      
      int irc = gnav->clk( satname, epoT, &clk, &clkrms, &dclk, msk_health );

#ifdef DEBUG
   cout << "gsatdata: " << _epoch.str_ymdhms() << " " << satname << fixed << setprecision(3)
        << " " << satclk << " " << satclk2 << " " << (satclk - satclk2)*CLIGHT
        << " " <<    clk << " " << clkrms << " " << dclk << " " << P3  << " " << irc << endl;
#endif	 
      if( irc < 0 || cnt++ > 25 ){
        if( _log ) _log->comment(2, "gsatdata"," satelite " + satname 
                                 + _epoch.str_ymdhms(" clocks not calculated (irc|iter) for epoch: "));
        //        cerr << satname + _epoch.str_ymdhms(" clocks not calculated (irc|iter) for epoch: ") << endl;
        return -1;
      }
      satclk = clk;
    }
  }else{
    epoT = _epoch;
    int irc = gnav->clk( satname, epoT, &satclk, &clkrms, &dclk, msk_health );
    if( irc < 0 ){	    
      if( _log )
      _log->comment(2, "gsatdata"," satelite " + satname 
                    + _epoch.str_ymdhms(" clocks not calculated for epoch "));
      //				 cout << satname + _epoch.str_ymdhms(" clocks not calculated for epoch ") << endl;
      return -1;
    }
  }

  int irc = 0;
  irc = gnav->pos( satname, epoT, xyz, var, vel, msk_health );

//   if (_low_prec) irc = gnav->nav( satname, epoT, xyz, var, vel ); // sometimes problem !?
//   else           irc = gnav->pos( satname, epoT, xyz, var, vel ); // OK!
   // =======================================================================
   // TADY JE OBCAS PROBLEM S NAV() PRO VYPOCET ELEVACE/AZIMUTU V ANUBIS!!!!
   // =======================================================================

  if( irc < 0 ) {
    if( _log )
    _log->comment(2, "gsatdata"," satelite " + satname
                  + _epoch.str_ymdhms(" coordinates not calculated for epoch "));
    //    cout << satname + _epoch.str_ymdhms(" coordinates not calculated for epoch ") << endl;
    return -1;
  }  
  
  t_gtriple txyz(xyz);
  t_gtriple tvel(vel);   
  
  // relativistic correction
  // WARNING: GLONASS clk already include the correction if broadcast eph are used !!!!!
  if( gs != GLO ||  
      (gs == GLO && gnav->id_type() == t_gdata::ALLPREC) ){
	  double rel = 2.0 * (txyz[0] * vel[0] + txyz[1] * vel[1] + txyz[2] * vel[2]) / CLIGHT / CLIGHT; //default
	  shared_ptr<t_geph> eph= gnav->find(satname, epoT);
	  if (gs == BDS && eph && gnav->id_type() == t_gdata::ALLRTCM) {
		  shared_ptr <t_gnavbds> gnavb = dynamic_pointer_cast <t_gnavbds>(eph);
		  double clk0 = 0.0, dclk0 = 0.0, clkrms0 = 0.0;
		  gnavb->clk(epoT, &clk0, &clkrms0, &dclk0, msk_health);
		  rel = gnavb->rel();
	  }
	  else if (gs == GAL&&eph && gnav->id_type() == t_gdata::ALLRTCM) {
		  shared_ptr <t_gnavgal> gnave = dynamic_pointer_cast <t_gnavgal>(eph);
		  double clk0 = 0.0, dclk0 = 0.0, clkrms0 = 0.0;
		  gnave->clk(epoT, &clk0, &clkrms0, &dclk0, msk_health);
		  rel = gnave->rel();
	  }

	  if (rel == 0.0) {
		  if (_log)
			  _log->comment(0, "gsatdata", " satelite " + satname
				  + _epoch.str_ymdhms(" relativity correction not calculated for epoch "));		  
		  return -1;
	  }

	  satclk -= rel;
  }     
   
  // filling gsatdata
  _satcrd = txyz;
  _satvel = tvel;
  _clk    = satclk*CLIGHT;
 
#ifdef DEBUG   
  ostringstream os;
  os << "gsatdata "<< satname
        << " CRD " << fixed << setprecision(3)
        <<  "  "   << epoT.str_ymdhms()
        << " X "   << setw(14) << txyz[0]
        << " Y "   << setw(14) << txyz[1]
        << " Z "   << setw(14) << txyz[2]
        << " T "   << setw(14) << satclk*1000.0*1000.0;
  //if (_log) _log->comment(2, "gsatdata", os.str());
  cout<<os.str()<<endl;
#endif
  return 1;
}   
   

// clean internal function
// ----------
void t_gsatdata::_clear()
{
  gtrace("t_gsatdat::_clear");   
   
  t_gobsgnss::_clear();
//  _satcrd[0] = 0.0;
//  _satcrd[1] = 0.0;
//  _satcrd[2] = 0.0;
  _ele = 0.0;
  _azi = 0.0;
  _rho = 0.0;
  _clk = 0.0;
}



// clean internal function
// ----------
bool t_gsatdata::_valid() const
{   
  gtrace("t_gsatdat::_valid");   
   
   if( _rho   == 0.0  ||   // single validity identification for gsatdata!
       t_gobsgnss::_valid() ) return false;

  return true; 
}

// reset eclipsing flag
// ---------------------
void  t_gsatdata::setecl(bool ecl)
{
	_eclipse = ecl;
}

// get eclipsing
// ---------------------
bool t_gsatdata::ecl()
{
  gtrace("t_gsatdat::ecl");   
   
   return _eclipse;
}
   
// determi ne wheather eclipsed or not
// -------------------------------------
void t_gsatdata::addecl(map<string, t_gtime>& lastEcl)
{
   gtrace("t_gsatdat::addecl");   

   if(fabs(_beta()) < EPS0_GPS && fabs(_orb_angle()) < EPS0_GPS){
      _eclipse = true;
      lastEcl[_satid] = _epoch;
      return;
   }else{
			auto itLast = lastEcl.find(_satid);
			if(itLast != lastEcl.end()){
				 double tdiff = _epoch.diff(itLast->second);
				 if(abs(tdiff) <= POST_SHADOW) _eclipse = true;
				 else _eclipse = false;
			}
   }   
}
   
// add postfit residuals
//------------------
void t_gsatdata::addres(RESIDTYPE restype, GOBSTYPE type, double res)
{
   _gmutex.lock();
   
   if(restype == RES_ORIG){
	  if(type == TYPE_C) _code_res_orig.push_back(res);
	  if(type == TYPE_L) _phase_res_orig.push_back(res);
   }

   if(restype == RES_NORM){
	  if(type == TYPE_C) _code_res_norm.push_back(res);
	  if(type == TYPE_L) _phase_res_norm.push_back(res);
   }   
   
   _gmutex.unlock();
   return;
}
      

// get postfit residuals
//------------------
vector<double> t_gsatdata::residuals(RESIDTYPE restype, GOBSTYPE type)
{
   _gmutex.lock();
   
   vector<double> res;
   
   if(restype == RES_ORIG){
	  if(     type == TYPE_C) res = _code_res_orig;
	  else if(type == TYPE_L) res = _phase_res_orig;
   }
   
   if(restype == RES_NORM){
	  if(     type == TYPE_C) res = _code_res_norm;
	  else if(type == TYPE_L) res = _phase_res_norm;
   }   
   
   _gmutex.unlock(); return res;
}

// clean residuals
// ---------------   
void t_gsatdata::clear_res(RESIDTYPE restype)
{
   if(restype == RES_ORIG){
	  _code_res_orig.clear();
	  _phase_res_orig.clear();
   }

   if(restype == RES_NORM){
	  _code_res_norm.clear();
	  _phase_res_norm.clear();   
   }  
}
   
// add postfit residuals
//------------------
void t_gsatdata::addwind(const double& wind)
{
   _gmutex.lock();
   
   _wind = wind;

   _gmutex.unlock();
   return;
}
      

// get stored wind up
//------------------
double t_gsatdata::wind()
{
   _gmutex.lock();
   
   _gmutex.unlock(); return _wind;
}   
      
// get hydrostatic mapping factor
//------------------
double t_gsatdata::mfH()
{
   _gmutex.lock();
   
   _gmutex.unlock(); return _mfH;
}
   
// get wet mapping factor
//------------------
double t_gsatdata::mfW()
{
   _gmutex.lock();
   
   _gmutex.unlock(); return _mfW;
}   
   
// get tropo gradient mapping factor
//------------------
double t_gsatdata::mfG()
{
   _gmutex.lock();
   
   _gmutex.unlock(); return _mfG;
}   
   
// Sun elevation relative to orbital plane
double t_gsatdata::_beta()
{
  gtrace("t_gsatdat::beta");   
   
   // test if already calculated
   if(!double_eq(_beta_val, 999)) return _beta_val;
   
   if(_satcrd.zero()) return 999;
   
   double beta = 0.0;
   double dt   = 300;
   
   double dmjd = _epoch.dmjd();
   t_gephplan eph;
   ColumnVector Sun = eph.sunPos(dmjd, false).crd_cvect();   //ICRF
   double gmt    = eph.gmst(dmjd);
   double gmt_dt = eph.gmst(dmjd + dt/86400.0);
   
   ColumnVector Satcrd = _satcrd.crd_cvect();
   ColumnVector Satvel = _satvel.crd_cvect();
   ColumnVector Satcrd_dt = Satcrd + Satvel*dt;   // 300s for extrapolation
   
   // prec. and nut. matrix should not change significantly in dt
   t_geop80 eop;
   Matrix prec = eop.precMatrix(dmjd);
   Matrix nut  = eop.nutMatrix(dmjd);
      
   // ITRF -> ICRF
   Satcrd     = prec.i() * nut.i() * rotZ(-gmt)*Satcrd;
   Satcrd_dt  = prec.i() * nut.i() * rotZ(-gmt_dt)*Satcrd_dt;

   ColumnVector n = crossproduct(Satcrd, Satcrd_dt);
   
   n /= n.NormFrobenius(); 
   
   ColumnVector nSun = Sun / Sun.NormFrobenius();

   double cosa = DotProduct(nSun,n);
//   if(cosa < 0) cosa *= -1;
   
   beta = G_PI/2.0 - acos(cosa);

   _beta_val = beta;
   
#ifdef DEBUG
   cout << "Beta angle " << sat() << " " << _epoch.str_hms() << endl
     << " Sat pos: " << _satcrd
     << " Sat vel: " << _satvel
     << " Beta: "    << beta*R2D << endl;
#endif   
   
  return beta; 
}

// Orbit angle
double t_gsatdata::_orb_angle()
{
  gtrace("t_gsatdat::orb_angle");   

   // test if already calculated
   if(!double_eq(_orb_angle_val, 999)) return _orb_angle_val;

   if(_satcrd.zero()) return 999;

   double mi = 0.0;
   double dt = 30;
   
   double dmjd = _epoch.dmjd();
   t_gephplan eph;
   ColumnVector Sun = eph.sunPos(dmjd, false).crd_cvect();   //ICRF
   double gmt    = eph.gmst(dmjd);
   double gmt_dt = eph.gmst(dmjd + dt/86400.0);         

   ColumnVector Satcrd = _satcrd.crd_cvect();
   ColumnVector Satvel = _satvel.crd_cvect();
   ColumnVector Satcrd_dt = Satcrd + Satvel*dt;   // 30s for extrapolation   

   // prec. and nut. matrix should not change significantly in dt
   t_geop80 eop;
   Matrix prec = eop.precMatrix(dmjd);
   Matrix nut  = eop.nutMatrix(dmjd);
      
   // ITRF -> ICRF
   Satcrd     = prec.i() * nut.i() * rotZ(-gmt)*Satcrd;
   Satcrd_dt  = prec.i() * nut.i() * rotZ(-gmt_dt)*Satcrd_dt;   
   
   ColumnVector n = crossproduct(Satcrd, Satcrd_dt);
   
   ColumnVector es = Satcrd / Satcrd.NormFrobenius();
      
   ColumnVector eSun = Sun / Sun.NormFrobenius();

   ColumnVector p = crossproduct(Sun,n);
   p /= p.NormFrobenius();

   double E = acos(DotProduct(es, p));
   /* mi = G_PI / 2.0 + (DotProduct(es, eSun) <= 0 ? -E : E); */
   /* if (mi < -G_PI / 2.0) mi += 2.0 * G_PI; */
   /* else if (mi >= G_PI / 2.0) mi -= 2.0*G_PI; */

	 double SunSat = acos(DotProduct(es, eSun));

	 if( SunSat > G_PI/2){
      if(E <= G_PI/2) mi = G_PI/2 - E;
      else            mi = G_PI/2 - E;
	 }else{
      if(E <= G_PI/2) mi = G_PI/2 + E;
      else            mi = E - G_PI - G_PI/2;
	 }

   _orb_angle_val = mi;
     
#ifdef DEBUG
   cout << "Orbit angle " << sat() << " " << _epoch.sod() << " " //str_hms() << " "
     << " Sat pos: " << _satcrd
     << " Sat vel: " << _satvel
     << " Angle: "    << mi*R2D << endl;
#endif   
   
  return mi; 
}   
   
} // namespace
