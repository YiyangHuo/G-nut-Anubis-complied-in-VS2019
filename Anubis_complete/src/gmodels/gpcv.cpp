
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
#include <iostream>
#include <iomanip>

#include "gmodels/gpcv.h"
#include "gutils/gsysconv.h"
#include "gutils/gtypeconv.h"
#include "gutils/gmatrixconv.h"
#include "gmodels/ginterp.h"
#include "gmodels/gpppmodel.h"

using namespace std;

namespace gnut {   

// -------------------------------
// T_GPCV - single antenn patterns
// -------------------------------

// constructor
// ----------
t_gpcv::t_gpcv()
 : _trans(true),    // transmitter (default yes)
   _anten(""),      // antenna type
   _ident(""),      // antenna identification
   _svcod(""),      // SVN code
   _method(""),     // calibartion method
   _source(""),     // source of calibration
   _snxcod("")      // sinex code 
{
  _beg.tsys(t_gtime::GPS);
  _end.tsys(t_gtime::GPS);
  id_type(t_gdata::PCV);   
}

// destructor
// ----------
t_gpcv::~t_gpcv(){}
   
   
// return correction ARP/CoM -> phase center
// zen [rad] - for receiver antenna, it is zenith angle to the satellite
//           - for satellite antenna, it is nadir angle (receiver zenith minus
//             geocentric (space) angle btw XYZsat and XYZrec)
// azi [rad] - for receiver antenna -> need to be implemented !
//           - for satellite antenna not interesting (only if horizontal eccentricity exists!)
// ----------
double t_gpcv::pco( double zen, double azi, GFRQ frq )
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();

  double corr = 0.0;
  if( zen > G_PI ){
    ostringstream lg;
    lg << "not valid zenith angle:" << zen*R2D << endl;
    if( _log ) _log->comment(0,"gpcv",lg.str());
     _gmutex.unlock(); return corr;
  }

  // satellite only PCO (Z-offset) mapped to rec-sat direction (approximated)
  if( _mappco.find(frq) != _mappco.end() ){
     corr = _mappco[frq][2] * cos(zen);  // for satellite zen should to be zen-alfa
  }
   
  // SATELLITE HORIZONTAL ECCENTRICITIES NOT YET IMPLEMENTED !!!!
  // RECEIVER PCO NOT YET IMPLEMENTED (azimut/zenith dependent)
   
  _gmutex.unlock(); return corr;
}

// Update satellite coordinates according to pco
// -----------------------------------------------
int t_gpcv::pcoS( t_gsatdata& satdata, t_gtriple& pco, GOBS_LC lc)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();

  t_gtime epo         = satdata.epoch();
  t_gtriple Satcrd_t  = satdata.satcrd();
  t_gtriple Satvel_t  = satdata.satvel();      // TEMPORARY
  GSYS sys            = satdata.gsys();
  string sat          = satdata.sat();  
  ColumnVector Satcrd = Satcrd_t.crd_cvect();
  ColumnVector Satvel = Satvel_t.crd_cvect();     // TEMPORARY
   
  t_gtriple apcf1, apcf2, apcf3, apcf4, apcf5, apcLC;

// JD: New flexible way of defining L3 frequency for multi-GNSS
  GOBSBAND b1 = t_gsys::band_priority(sys, FREQ_1);
  GOBSBAND b2 = t_gsys::band_priority(sys, FREQ_2);
  GOBSBAND b3 = t_gsys::band_priority(sys, FREQ_3);
  GOBSBAND b4 = t_gsys::band_priority(sys, FREQ_4);
  GOBSBAND b5 = t_gsys::band_priority(sys, FREQ_5);

     GFRQ f1 = t_gsys::freq_priority(sys, FREQ_1);
     GFRQ f2 = t_gsys::freq_priority(sys, FREQ_2);
     GFRQ f3 = t_gsys::freq_priority(sys, FREQ_3);
     GFRQ f4 = t_gsys::freq_priority(sys, FREQ_4);
     GFRQ f5 = t_gsys::freq_priority(sys, FREQ_5);

   if( _mappco.find(f1) != _mappco.end() && _mappco.find(f2) != _mappco.end() ){
      apcf1 = _mappco[f1];
      apcf2 = _mappco[f2];
      apcf3 = _mappco[f3];
      apcf4 = _mappco[f4];
      apcf5 = _mappco[f5];
//    cout <<  "JD: pcoS maps " << satdata.sat()
//                       << " " << t_gfreq::gfreq2str(f1) << " " << apcf1
//                       << " " << t_gfreq::gfreq2str(f2) << " " << apcf2 << endl;
   }else{
      t_map_offs offs = GNSS_PCO_OFFSETS();      
      t_map_offs::iterator itSYS = offs.find(sys);
      if(itSYS != offs.end()){
				 t_map_pcos::iterator itb1 = itSYS->second.find(b1);
				 t_map_pcos::iterator itb2 = itSYS->second.find(b2);
				 t_map_pcos::iterator itb3 = itSYS->second.find(b3);
				 t_map_pcos::iterator itb4 = itSYS->second.find(b4);
				 t_map_pcos::iterator itb5 = itSYS->second.find(b5);
	 
				 if(itb1 != itSYS->second.end()) {
						apcf1[0] = itb1->second[0];
						apcf1[1] = itb1->second[1];
						apcf1[2] = itb1->second[2];	    
				 }
				 
				 if(itb2 != itSYS->second.end()) {
						apcf2[0] = itb2->second[0];
						apcf2[1] = itb2->second[1];
						apcf2[2] = itb2->second[2];
				 }
	 
				 if(itb3 != itSYS->second.end()) {
						apcf3[0] = itb3->second[0];
						apcf3[1] = itb3->second[1];
						apcf3[2] = itb3->second[2];
				 }	 
				 
				 if(itb4 != itSYS->second.end()) {
						apcf4[0] = itb4->second[0];
						apcf4[1] = itb4->second[1];
						apcf4[2] = itb4->second[2];
				 } 

				 if(itb5 != itSYS->second.end()) {
						apcf5[0] = itb5->second[0];
						apcf5[1] = itb5->second[1];
						apcf5[2] = itb5->second[2];
				 } 
				 
      }
	 }

	 for (int i=0; i<=2; i++){
			apcf1[i] /= 1000.0;
			apcf2[i] /= 1000.0;
			apcf3[i] /= 1000.0;
			apcf4[i] /= 1000.0;
			apcf5[i] /= 1000.0;     
	 }      

   // PCO for linear combination
	 if (lc == LC_IF){
     double koef1 = 0.0;
     double koef2 = 0.0;       
     satdata.coef_ionofree(BAND_1, koef1, BAND_2, koef2);
     apcLC = apcf1*koef1 + apcf2*koef2;
   }else if (lc == LC_L1){ apcLC = apcf1;
   }else if (lc == LC_L2){ apcLC = apcf2;
   }else {apcLC = apcf1;}
   		 
	 pco = apcLC;	
   
#ifdef DEBUG
   cout << "PCO calculation for " << satdata.sat() << " " << anten() << " " << epo.str_hms() << endl;
   cout << fixed << setprecision(4);
   cout << "pco f1: " << apcf1 << " m" << endl
        << "pco f2: " << apcf2 << " m" << endl      
        << "pco lc: " << apcLC << " m" << endl;
//   int ooo; cin >> ooo; 
#endif

   _gmutex.unlock(); return 1;
}

// Receiver pco
// -----------------------------------------------
int t_gpcv::pcoR( t_gsatdata& satdata, t_gtriple& dx, t_gtriple& site, GOBS_LC lc)
{
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
  _gmutex.lock();
   
  t_gtriple apcf1, apcf2, apcLC;
  
  GSYS sys = satdata.gsys();
  
  // Temporary - GAL site antenna calibration not available
  if(sys == GAL) sys = GPS;
  
// JD: New flexible way of defining L3 frequency for multi-GNSS
  GFRQ f1 = t_gsys::freq_priority(sys, FREQ_1);
  GFRQ f2 = t_gsys::freq_priority(sys, FREQ_2);

  if( _mappco.find(f1) == _mappco.end() || _mappco.find(f2) == _mappco.end() ){
    _gmutex.unlock(); return -1;
  }

  apcf1 = _mappco[f1];
  apcf2 = _mappco[f2];
   
//  cout <<  "JD: pcoR maps " << anten()
//                     << " " << t_gfreq::gfreq2str(f1) << " " << apcf1
//                     << " " << t_gfreq::gfreq2str(f2) << " " << apcf2 << endl;

  for (int i=0; i<=2; i++){
     apcf1[i] /= 1000.0;
     apcf2[i] /= 1000.0;     
  }
   
   // PCO for linear combination
    if (lc == LC_IF){
      double koef1 = 0.0;
      double koef2 = 0.0;
      satdata.coef_ionofree(BAND_1, koef1, BAND_2, koef2);              
      apcLC = apcf1*koef1 + apcf2*koef2;
    }else if(lc == LC_L1) {apcLC = apcf1;
    }else if(lc == LC_L2) {apcLC = apcf2;
    }else apcLC = apcf2;
   
   t_gtriple ell(0.0, 0.0, 0.0);
   xyz2ell(site, ell, false);
   neu2xyz(ell, apcLC, dx);
   
#ifdef DEBUG
   cout << "PCO calculation for " << anten() << endl;
   cout << fixed << setprecision(4);
   cout << "pco f1: " << apcf1 << " m" << endl
        << "pco f2: " << apcf2 << " m" << endl      
        << "pco lc: " << apcLC << " m" << endl
        << fixed << setprecision(5)
        << "Sat dX:\n" << dx << endl;
#endif
   
   _gmutex.unlock(); return 1;   
}

// return type
// ---------------
int t_gpcv::pco_proj(double& corr, t_gsatdata& satdata, t_gtriple& site, t_gtriple& dx, GOBS_LC lc)
{
   gtrace("t_gpcv::pco_proj");
   
   t_gtriple satcrd = satdata.satcrd();
   t_gtriple SatRec(0.0, 0.0, 0.0);
   SatRec = satcrd - site;
   ColumnVector e = SatRec.unitary();
   
   corr = DotProduct(dx.crd_cvect(),e);
   
   return 1;
}


// return type
// ----------
int t_gpcv::pcvS( double& corrLC, t_gsatdata& satdata, t_gtriple& site, GOBS_LC lc)
{   
  
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();

  GSYS sys = satdata.gsys();

// JD: New flexible way of defining L3 frequency for multi-GNSS
   GFRQ f1 = t_gsys::freq_priority(sys, FREQ_1);
   GFRQ f2 = t_gsys::freq_priority(sys, FREQ_2);
// GFRQ f3 = t_gsys::freq_priority(sys, FREQ_3);
// GFRQ f4 = t_gsys::freq_priority(sys, FREQ_4);
// GFRQ f5 = t_gsys::freq_priority(sys, FREQ_5);

  if( _mapzen.find(f1) == _mapzen.end() || _mapzen.find(f2) == _mapzen.end() ){
    _gmutex.unlock(); return -1;
  }
   
   map<double, double> mapDataf1;
   mapDataf1 = _mapzen[f1];
   
   map<double, double> mapDataf2;
   mapDataf2 = _mapzen[f2];   
     
   double eleS = satdata.ele();     
   double rS = satdata.satcrd().crd_cvect().NormFrobenius();
   double rR =         site.crd_cvect().NormFrobenius();
   
   double sinz = (rR/rS)*cos(eleS);
   double zen = asin(sinz) * R2D;

   double corrf1 = 0.0;
   double corrf2 = 0.0;
   
   t_ginterp interp;
   if ( interp.linear(mapDataf1, zen, corrf1) < 0 || interp.linear(mapDataf2, zen, corrf2) < 0 ){
      _gmutex.unlock(); return -1;
   }

//  cout <<  "JD: pcvS maps " << satdata.sat()
//                     << " " << t_gfreq::gfreq2str(f1) << " " << corrf1
//                     << " " << t_gfreq::gfreq2str(f2) << " " << corrf2 << endl;

   // [mm] -> [m]
   corrf1 /= 1000.0;
   corrf2 /= 1000.0;   
   
   // PCV for linear combination
    if (lc == LC_IF){
      double koef1 = 0.0;
      double koef2 = 0.0;       
      satdata.coef_ionofree(BAND_1, koef1, BAND_2, koef2);              
      corrLC = corrf1*koef1 + corrf2*koef2;
    }else if(lc == LC_L1) {corrLC = corrf1;
    }else if(lc == LC_L2) {corrLC = corrf2;
    }else corrLC = corrf1;   
   
#ifdef DEBUG
   cout << "Sat PCV interpolation" << endl;
   cout << "Nadir = "    << zen 
	<< " corrf1 =  " << corrf1*1000 << " mm, "
	<< " corrf2 =  " << corrf2*1000 << " mm, "
	<< " corrLC =  " << corrLC*1000 << " mm, "	
	<< " PRN: " << anten() << endl;
#endif
   
   _gmutex.unlock(); return 1;
}

// return type
// ----------
int t_gpcv::pcvR( double& corrLC, t_gsatdata& satdata, GOBS_LC lc)
{   
#ifdef BMUTEX   
   boost::mutex::scoped_lock lock(_mutex);
#endif 
   _gmutex.lock();
  
   GSYS sys = satdata.gsys();
  
   // Temporary - GAL site antenna calibration not available
   if(sys == GAL) sys = GPS;
   
   // JD: New flexible way of defining L3 frequency for multi-GNSS
   GFRQ f1 = t_gsys::freq_priority(sys, FREQ_1);
   GFRQ f2 = t_gsys::freq_priority(sys, FREQ_2);
// GFRQ f3 = t_gsys::freq_priority(sys, FREQ_3);
// GFRQ f4 = t_gsys::freq_priority(sys, FREQ_4);
// GFRQ f5 = t_gsys::freq_priority(sys, FREQ_5);

   double zen = G_PI/2.0 - satdata.ele();        
   zen *= R2D;
   
   double azi = satdata.azi();
   azi *= R2D;
   
   double corrf1 = 0.0;
   double corrf2 = 0.0;

   if( _azi_dependent(f1) &&  _azi_dependent(f2)) {  // AZI-dependant calibration available
     t_gpair p_az(azi, zen);
   
     map<t_gpair, double> mapDataf1;   
     map<t_gpair, double> mapDataf2;

     for(map<GFRQ, t_map_A>::iterator itGFRQ = _mapazi.begin(); itGFRQ != _mapazi.end(); itGFRQ++){
       GFRQ f = itGFRQ->first;      
       if(f != f1 && f != f2) continue;
	 
       map<double, t_map_Z>::iterator itA1  = itGFRQ->second.lower_bound(azi);
       if(itA1 == itGFRQ->second.end() || itA1 == itGFRQ->second.begin()) { _gmutex.unlock(); return -1; }
       map<double, t_map_Z>::iterator itA2  = itA1; itA2--;
       
       map<double, double>::iterator itZ1  = itA1->second.lower_bound(zen);
       if(itZ1 == itA1->second.end() || itZ1 == itA1->second.begin()) { _gmutex.unlock(); return -1; }
       map<double, double>::iterator itZ2  = itZ1; itZ2--;

       map<double, double>::iterator itZ3  = itA2->second.lower_bound(zen);
       if(itZ3 == itA2->second.end() || itZ3 == itA2->second.begin()) { _gmutex.unlock(); return -1; }      
       map<double, double>::iterator itZ4  = itZ3; itZ4--;      
      
       t_gpair p1(itA1->first, itZ1->first);
       t_gpair p2(itA1->first, itZ2->first);            
       t_gpair p3(itA2->first, itZ3->first);
       t_gpair p4(itA2->first, itZ4->first);
	 
       if (f == f1) {
         mapDataf1[p1] = itZ1->second; //cout << p1[0] << " : " << p1[1] << "  " << itZ1->second << endl;
         mapDataf1[p2] = itZ2->second; //cout << p2[0] << " : " << p2[1] << "  " << itZ2->second << endl;
         mapDataf1[p3] = itZ3->second; //cout << p3[0] << " : " << p3[1] << "  " << itZ3->second << endl;
         mapDataf1[p4] = itZ4->second; //cout << p4[0] << " : " << p4[1] << "  " << itZ4->second << endl;
       }else if (f == f2) {
         mapDataf2[p1] = itZ1->second;
         mapDataf2[p2] = itZ2->second;
         mapDataf2[p3] = itZ3->second;
         mapDataf2[p4] = itZ4->second;	 
       }
     }     
   
     t_ginterp interp;
     if ( interp.bilinear(mapDataf1, p_az, corrf1) < 0 || interp.bilinear(mapDataf2, p_az, corrf2) < 0 ){
       _gmutex.unlock(); return -1;
     }
   }else{  // AZI-dependant calibration NOT available (only NOAZI)
     if( _mapzen.find(f1) == _mapzen.end() || _mapzen.find(f2) == _mapzen.end() ) {      	   
       _gmutex.unlock(); return -1;
     }
     map<double, double> mapDataf1;
     mapDataf1 = _mapzen[f1];
     
     map<double, double> mapDataf2;
     mapDataf2 = _mapzen[f2];
     
     t_ginterp interp;
     if ( interp.linear(mapDataf1, zen, corrf1) < 0 || interp.linear(mapDataf2, zen, corrf2) < 0 ) {	   
       _gmutex.unlock(); return -1;
     }            
   }
   
//  cout <<  "JD: pcvR maps " << anten()
//                     << " " << t_gfreq::gfreq2str(f1) << " " << corrf1
//                     << " " << t_gfreq::gfreq2str(f2) << " " << corrf2 << endl;

   // [mm] -> [m]
   corrf1 /= 1000.0;
   corrf2 /= 1000.0;
   
   // PCV for linear combination
   if (lc == LC_IF){
     double koef1 = 0.0;
     double koef2 = 0.0;
     satdata.coef_ionofree(BAND_1, koef1, BAND_2, koef2);              
     corrLC = corrf1*koef1 + corrf2*koef2;
   }else if (lc == LC_L1) {corrLC = corrf1;
   }else if (lc == LC_L2) {corrLC = corrf2;    
   }else corrLC = corrf1;   
   
#ifdef DEBUG
   cout << "Rec PCV interpolation" << endl;
   cout << "Zenith = "    << zen 
	<< " corrf1 =  " << corrf1*1000 << " mm, "
	<< " corrf2 =  " << corrf2*1000 << " mm, "
	<< " corrLC =  " << corrLC*1000 << " mm, "	
	<< " Ant: " << anten() << " " 
        << " Sat: " << satdata.sat() << endl;
//   int ooo; cin >> ooo;     
#endif   
   
   _gmutex.unlock(); return 1;   
}

// does the calibration contain azi-depenedant data?   
bool t_gpcv::_azi_dependent(GFRQ f)
{     
   t_map_azi::iterator it = _mapazi.find(f);
   
   bool ret = false; 
   
   if(it == _mapazi.end()){
      ret = false; 
   }else{
      
      int size = _mapazi[f].size();
      if(size > 0) {ret = true;}
      else {ret = false;}
      
   }
     
   return ret;   
}
   
} // namespace
