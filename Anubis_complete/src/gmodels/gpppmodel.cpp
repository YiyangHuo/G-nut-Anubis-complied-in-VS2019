
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

#include "gmodels/gpppmodel.h"
#include "gdata/gpppdata.h"

namespace gnut {

   // Constructors
   // -------------------
   t_gpppmodel::t_gpppmodel(string site, t_glog* glog, t_gsetbase* settings)
      : t_gsppmodel(site, glog, settings),
      _gallobj(0)
   {
      gtrace("t_gpppmodel::t_gpppmodel");

      /* _site = site; */
      /* _settings = settings; */
      /* _log = glog; */

      if (_trpModStr == EXTERN)  _tropoModel = make_shared<t_gtropo>();

      //   if(      _trpModStr.compare("met") == 0 ){ _tropoModel = make_shared<t_gtropo>(site); } // first option for site
      //   else if( _trpModStr.compare("nwm") == 0 ){ _tropoModel = make_shared<t_gtropo>();     } // second option for site

      _phase = dynamic_cast<t_gsetproc*>(_settings)->phase();
      _grad_mf = dynamic_cast<t_gsetproc*>(_settings)->grad_mf();
   }

   // Destructor
   // --------------------------
   t_gpppmodel::~t_gpppmodel()
   {
      gtrace("t_gpppmodel::~t_gpppmodel");
   }


   // wind-up correction
   // ----------
   double t_gpppmodel::windUp(t_gsatdata& satdata, const ColumnVector& rRec)
   {
     gtrace("t_gpppmodel::windUp");

     t_gtime epoch = satdata.epoch();
     string prn = satdata.sat();
     ColumnVector rSat = satdata.satcrd().crd_cvect();
     //GSYS sys          = satdata.gsys();

     double Mjd = epoch.mjd(true) + epoch.sod(true) / 86400.0;

     // First time - initialize to zero
     // -------------------------------
     map<string, double>::iterator it = _windUpTime.find(prn);
     if (it == _windUpTime.end()){
       _windUpTime[prn] = Mjd;
       _windUpSum[prn] = 0.0;
     }
     
     // Compute the correction for new time
     // -----------------------------------
     else if (_windUpTime[prn] != Mjd) {
       
       _windUpTime[prn] = Mjd;
       
       ColumnVector rho = rRec - rSat;
       rho /= rho.norm_Frobenius();
       
       ColumnVector i(3);
       ColumnVector j(3);
       ColumnVector k(3);

       // attitude model
       string antype = "";
       if (_gallobj != 0){
         shared_ptr<t_gobj>  sat_obj = _gallobj->obj(satdata.sat());
         shared_ptr<t_gpcv>  sat_pcv;
         if (sat_obj != 0)  sat_pcv = sat_obj->pcv(satdata.epoch());
         if (sat_pcv != 0)  antype = sat_pcv->anten();
       }
       attitude(satdata, antype, i, j, k);

       // Effective Dipole of the GPS Satellite Antenna
       // ---------------------------------------------
       ColumnVector dipSat = i - rho * DotProduct(rho, i) - crossproduct(rho, j);
       
       // Receiver unit Vectors rx, ry
       // ----------------------------
       ColumnVector rx(3);
       ColumnVector ry(3);
       
       double recEll[3]; xyz2ell(rRec.data(), recEll, false);
       double neu[3];
       
       neu[0] = 1.0;
       neu[1] = 0.0;
       neu[2] = 0.0;
       neu2xyz(recEll, neu, rx.data());
       
       neu[0] = 0.0;
       neu[1] = -1.0;
       neu[2] = 0.0;
       neu2xyz(recEll, neu, ry.data());

       // Effective Dipole of the Receiver Antenna
       // ----------------------------------------
       ColumnVector dipRec = rx - rho * DotProduct(rho, rx) + crossproduct(rho, ry);
       
       // Resulting Effect
       // ----------------
       double alpha = DotProduct(dipSat, dipRec) / (dipSat.norm_Frobenius() * dipRec.norm_Frobenius());
       
       if (alpha > 1.0) alpha = 1.0;
       if (alpha < -1.0) alpha = -1.0;
       
       double dphi = acos(alpha) / 2.0 / G_PI;  // in cycles

       if (DotProduct(rho, crossproduct(dipSat, dipRec)) < 0.0) dphi = -dphi;
       
       _windUpSum[prn] = floor(_windUpSum[prn] - dphi + 0.5) + dphi;
     }
     
#ifdef DEBUG
      cout << "Wind up " << epoch.str_ymdhms() << " " << prn << " " << _windUpSum[prn] << " " << satdata.orb_angle()*R2D << endl;
//    int ooo; cin >> ooo;
#endif   
     
     satdata.addwind(_windUpSum[prn]);
     return _windUpSum[prn];
   }
   

   // set gallobj pointer
   // --------------------
   void t_gpppmodel::setOBJ(t_gallobj* obj)
   {
     _gallobj = obj;
     if (obj){
       _grec = obj->obj(_site);         
     }
   }

  // model computed range value (phase/code)
  // ----------
  double t_gpppmodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& gsatdata, t_gobs& gobs, bool com)
  {
    gtrace("t_gpppmodel::cmpObs");

    double spp_model = t_gsppmodel::cmpObs(epoch, param, gsatdata, gobs);
    if (spp_model < 0) return -1;

    // Cartesian coordinates to ellipsodial coordinates
    t_gtriple xyz, ell;
    ColumnVector cRec(3);
    
    if (param.getCrdParam(_site, xyz) > 0){
      cRec = xyz.crd_cvect();
    }
    else{
      xyz = _grec->crd_arp(epoch);
      cRec = xyz.crd_cvect();
    }
    xyz2ell(xyz, ell, false);
      
    t_gtriple satcrd = gsatdata.satcrd();
    ColumnVector cSat = satcrd.crd_cvect();

    // Wind up correction 
    double wind = 0.0;
    if (gobs.is_phase()) {        
      double wavelength = 0.0;
      if(_observ == IONO_FREE) {
            
        GSYS gs = gsatdata.gsys();
        GOBSBAND b1 = t_gsys::band_priority(gs, FREQ_1);
        GOBSBAND b2 = t_gsys::band_priority(gs, FREQ_2);
        
        wavelength = gsatdata.wavelength_L3(b1, b2);
      }else wavelength = gsatdata.wavelength(gobs.band());
         
      if (fabs(gsatdata.wind()) > 0) wind = gsatdata.wind() * wavelength;
      else wind = windUp(gsatdata, cRec) * wavelength;
    }     
    
    // Phase center variation correction
    double pcv_R = 0.0;  
    double pcv_S = 0.0;  
    double pco_R = 0.0;  
    double pco_S = 0.0;  
    if (_gallobj != 0){
      shared_ptr<t_gobj>  sat_obj = _gallobj->obj(gsatdata.sat());
      shared_ptr<t_gobj> site_obj = _gallobj->obj(_site);
      
      shared_ptr<t_gpcv>  sat_pcv;
      shared_ptr<t_gpcv> site_pcv;

      if (sat_obj != 0)  sat_pcv  = sat_obj->pcv(epoch);
      if (site_obj != 0) site_pcv = site_obj->pcv(epoch);

      GOBS_LC lc = LC_IF;
      if(_observ != IONO_FREE){
        if(t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_1) lc = LC_L1;
        if(t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_2) lc = LC_L2;
        if(t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_3) lc = LC_L3;
        if(t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_4) lc = LC_L4;
        if(t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_5) lc = LC_L5;
      }
      
      if(sat_pcv != 0 && com){
        // Satellite phase center variation
        sat_pcv->pcvS(pcv_S, gsatdata, xyz, lc);
        // Satellite phase center offset
        t_gtriple pco(0,0,0);
        if(sat_pcv->pcoS(gsatdata, pco, lc) > 0) {
          string antenna = sat_pcv->anten();
          t_gtriple dx(0,0,0);
          ColumnVector i(3);
          ColumnVector j(3);
          ColumnVector k(3);        
          this->attitude(gsatdata, antenna, i, j, k);
          dx[0] = pco[0]*i(1) + pco[1]*j(1) + pco[2]*k(1);
          dx[1] = pco[0]*i(2) + pco[1]*j(2) + pco[2]*k(2);
          dx[2] = pco[0]*i(3) + pco[1]*j(3) + pco[2]*k(3);
          sat_pcv->pco_proj(pco_S, gsatdata, xyz, dx, lc);
		  gsatdata.addpco(dx);
        }
        
      }
      
      if(site_pcv != 0){
        // Receiver phase center variation
        site_pcv->pcvR(pcv_R, gsatdata, lc);
        // Receiver phase center offset
        t_gtriple dx(0.0, 0.0, 0.0);          
        if(site_pcv->pcoR(gsatdata, dx, xyz, lc) > 0) {
          site_pcv->pco_proj(pco_R, gsatdata, xyz, dx, lc);
        }
        pco_R *= -1; 
      }
    }
    
#ifdef DEBUG
      cout << gsatdata.epoch().str_hms() <<  " " 
         << gsatdata.sat() << " " << gobs2str(gobs.gobs())
          << fixed << setprecision(3)
         << "  SPP_model "  << setw(15)<< spp_model
         << "  wind: "      << setw(6) << wind
         << "  Rec PCO: "   << setw(6) << pco_R          
         << "  Rec PCV: "   << setw(6) << pcv_R
         << "  Sat PCV: "   << setw(6) << pcv_S
         << "  Sat PCO: "   << setw(6) << pco_S    
         << endl;    
//         int ooo; cin >> ooo;
#endif          

      // Return value
      return spp_model +
        wind  +
        pco_R +
        pco_S +        
        pcv_R +
        pcv_S;
  }
  

   // Compute troposperic delay
   // -----------
   double t_gpppmodel::tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple ell, t_gsatdata& satdata)
   {
      gtrace("t_gpppmodel::tropoDelay");

      if (_tropoModel == 0){
         if (_log) _log->comment(0, "gppp", "Tropo Model setting is not correct. Default used! Check config.");
         else                cerr << "gppp - Tropo Model setting is not correct. Default used! Check config.\n";
         _tropoModel = make_shared<t_saast>();
      }

      double ele = satdata.ele();
      double azi = satdata.azi();

      double delay = 0.0;
      double zwd = 0.0;
      double zhd = 0.0;

      int i, j, k;
      i = param.getParam(_site, t_gpar::TRP, "");
      j = param.getParam(_site, t_gpar::GRD_N, "");
      k = param.getParam(_site, t_gpar::GRD_E, "");

      if (i >= 0){
         zwd = param[i].value();
         zhd = param[i].apriori();
      }
      else{
         if (_tropoModel != 0){
            zwd = _tropoModel->getZWD(ell, epoch);
            zhd = _tropoModel->getZHD(ell, epoch);
         }
      }

      double mfh, mfw, dmfh, dmfw;
      if (_tropo_mf == GMF){
         t_gmf mf;
         mf.gmf(epoch.mjd(), ell[0], ell[1], ell[2], G_PI / 2.0 - ele,
            mfh, mfw, dmfh, dmfw);
      }
      else if (_tropo_mf == COSZ){
         mfh = mfw = 1 / sin(ele);
         dmfh = dmfw = -(cos(ele)) / (sin(ele)*sin(ele));
      }

      satdata.addmfH(mfh);
      satdata.addmfW(mfw);

      delay = mfh * zhd + mfw * zwd;

      double grdN, grdE;
      grdN = grdE = 0.0;

      if (j >= 0 && k >= 0){
         if (_grad_mf == TILTING){
            grdN = param[j].value() * dmfw * cos(azi);
            grdE = param[k].value() * dmfw * sin(azi);
            satdata.addmfG(dmfw);
         }
         else if (_grad_mf == CHEN_HERRING){
            double mfg = 1.0 / (sin(ele)*tan(ele) + 0.0032);
            grdN = param[j].value() * 1000.0 * mfg * cos(azi);  grdN /= 1000.0;
            grdE = param[k].value() * 1000.0 * mfg * sin(azi);  grdE /= 1000.0;
            satdata.addmfG(mfg);
         }
         else if (_grad_mf == BAR_SEVER){
            double mfg = mfw * (1 / tan(ele));
            grdN = param[j].value() * mfg * cos(azi);
            grdE = param[k].value() * mfg * sin(azi);
            satdata.addmfG(mfg);
         }

         delay += grdN + grdE;
      }


#ifdef DEBUG
      cout << epoch.str("EPOCH: %H:%M:%S") << endl << fixed << setprecision(3);
      cout << satdata.sat() << " Ell:" << ell[0] << " " << ell[1] << " " << ell[2]
         << " Hydrostatic part: " << zhd //_tropoModel->getZHD(ell, epoch)
         << " Wet part: "        << zwd
         << " mfh: " << mfh
         << " mfw: " << mfw
         << " Delay: " << delay << endl << endl;
      //      int ooo; cin >> ooo;
#endif

      return delay;
   }

// satellite attitude model
// --------------------------------   
int t_gpppmodel::attitude_old(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    if(satdata.gsys() == GAL){      
         _ysm(satdata, i, j, k);
    }else{
      int irc = _yaw(satdata,antype, i, j, k);
      if (irc == 0)return 0;
   }
    
   return 1;
}
    
// satellite attitude model
// --------------------------------   
int t_gpppmodel::attitude(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
  if(satdata.gsys() == GPS){
    if(antype == "BLOCK II"){
      _attitude_GPSIIA(satdata, antype, i, j, k);
    }else if(antype == "BLOCK IIA"){
      _attitude_GPSIIA(satdata, antype, i, j, k);
    }else if(antype.find("BLOCK IIR") != string::npos){
      _attitude_GPSIIR(satdata, antype, i, j, k);
    }else if(antype == "BLOCK IIF"){
      _attitude_GPSIIF(satdata, antype, i, j, k);
    }else _ysm(satdata, i, j, k);
     
  }else if(satdata.gsys() == GLO){    
    _attitude_GLO(satdata, antype, i, j, k);
         
  }else if(satdata.gsys() == GAL){
    if(antype == "GALILEO-1"){
      _attitude_GAL1(satdata, antype, i, j, k);
    }else if(antype == "GALILEO-2"){
      _attitude_GAL2(satdata, antype, i, j, k);
    }else _ysm(satdata, i, j, k);
    
  }else if(satdata.gsys() == BDS){    
    _attitude_BDS(satdata, antype, i, j, k);
    
  }else if(satdata.gsys() == QZS){    
    _attitude_QZS(satdata, antype, i, j, k);        
  }
  
  if(i.NormFrobenius() == 0 || j.NormFrobenius() == 0 || k.NormFrobenius() == 0) return 0;
   
#ifdef DEBUG
   double angl = satdata.orb_angle();
   if(angl > G_PI/2) angl -= 2.0*G_PI;
   cout << satdata.sat() << " " << satdata.ele_deg() << " " << angl*R2D << " " << satdata.yaw()*R2D << " " << satdata.beta()*R2D << endl;  
#endif
    
  return 1;
}   
   
// Yaw-steering mode attitude model
// --------------------------------   
void t_gpppmodel::_ysm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
   double MJD = satdata.epoch().dmjd();
   
   t_gtriple satcrd_t  = satdata.satcrd();   
   ColumnVector satcrd = satcrd_t.crd_cvect();

   // Satelite-Earth unit vector
   k = -satcrd;

   // Along solar panel unit vector  
   ColumnVector sun = _ephplan.sunPos(MJD).crd_cvect();
    j = crossproduct(k,sun);

   // complete to satelite fixed right-hand coord system
   i = crossproduct(j,k);
    
    i = i / i.NormFrobenius();
    j = j / j.NormFrobenius();
    k = k / k.NormFrobenius();    

    _last_beta[satdata.sat()] = satdata.beta();
    
    double yaw = atan2(-satdata.beta(), sin(satdata.orb_angle()));
    satdata.yaw(yaw);
}
   
// orbit normal mode (i toward the velocity)
// --------------------------------   
void t_gpppmodel::_onm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{   
    double yaw = 0;
    _yaw2ijk(satdata, yaw, i, j, k);    
}   

// yaw model of different satellite
// Nominal attitude (same as GPS BLOCK II/IIA) without yaw maneuver for MEO and IGSO satellites ();
// Yaw - fixed attitude mode used for GEO satellites
// --------------------------------   
int t_gpppmodel::_yaw(t_gsatdata& satdata,string antype, ColumnVector& xs, ColumnVector& ys, ColumnVector& zs)
{
#if 0
   double exs[3], eys[3], ezs[3], r[3];
   int i;

   int opt = 2; //2 :precise yaw
   
   bool yawFix = t_gsys::bds_geo(satdata.sat()) || (satdata.beta()*R2D <= 4 && satdata.beta()*R2D >= -4);
   if (satdata.gsys() == GSYS::BDS && yawFix){  //Only for BDS, opt=3
      opt = 3;
   }

   double rs[6];
   for (int i = 0; i < 3; i++)  rs[i] = satdata.satcrd()[i];
   for (int i = 3; i < 6; i++)  rs[i] = satdata.satvel()[i - 3];

   string prn = satdata.sat();
   t_gtime epoch = satdata.epoch();

   for (i = 0; i < 3; i++) r[i] = -rs[i];
   if (!normv3(r, ezs)) return 0;
   /* satellite yaw attitude model */
   double  ep[6];
   ep[0] = epoch.year(); ep[1] = epoch.mon(); ep[2] = epoch.day();
   ep[3] = epoch.hour(); ep[4] = epoch.mins(); ep[5] = epoch.secs();
   gtime_t time = epoch2time(ep);

   if (!normv3(r, ezs)) return 0;
   if (!sat_yaw(time, prn.c_str(), antype.c_str(), opt, rs, exs, eys)) return 0;
   for (int i = 0; i < 3; i++)
   {
      xs(i + 1) = exs[i];
      ys(i + 1) = eys[i];
      zs(i + 1) = ezs[i];
   }
#endif
   return 1;
}


// Attitude modelling for GPS Block IIA
void t_gpppmodel::_attitude_GPSIIA(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R_GPSIIA[] = {
         0.1046,0.1230,0.1255,0.1249,0.1003,0.1230,0.1136,0.1169,0.1253,0.0999,
         0.1230,0.1230,0.1230,0.1230,0.1092,0.1230,0.1230,0.1230,0.1230,0.1230,
         0.1230,0.1230,0.1230,0.0960,0.0838,0.1284,0.1183,0.1230,0.1024,0.1042,
         0.1230,0.1100,0.1230
    };
    int sat = str2int(satdata.sat().substr(1,2));
    double R = R_GPSIIA[sat - 1] * D2R;

    double beta0 = atan2(MUDOT_GPS, R);
    
    if(fabs(satdata.orb_angle()) > G_PI/2 && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn(satdata, R, i, j, k);
    }else if(fabs(satdata.orb_angle()) < G_PI/2 && fabs(satdata.beta()) < beta0){
         // midnight maneuver
         _midnight_turn_GPSIIA(satdata, R, i, j, k);
    }else{
         _ysm(satdata, i, j, k); 
    }  

}
    
// Attitude modelling for GPS Block IIR    
void t_gpppmodel::_attitude_GPSIIR(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R = 0.2*D2R;             // maximal yaw hardware rate
    double beta0 = atan2(MUDOT_GPS, R);
    
    if(fabs(satdata.orb_angle()) > G_PI/2 && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn(satdata, R, i, j, k);
    }else if(fabs(satdata.orb_angle()) < G_PI/2 && fabs(satdata.beta()) < beta0){
         // midnight maneuver
         _midnight_turn(satdata, R, i, j, k);
    }else{              
         _ysm(satdata, i, j, k);
    }  

    i = -1 * i;      // X away from the Sum
    j = -1 * j;
    satdata.yaw(satdata.yaw() + G_PI);  
    
}
    
// noon maneuver
// -------------------------------------------
void t_gpppmodel::_noon_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta0 = atan2(MUDOT_GPS, R);
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    
    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }

    if (beta >= 0) R *= -1;
    
    double mi_s = G_PI - sqrt( beta0 * fabs(beta) - beta*beta );
    double yaw;
    double yaw_nom = atan2(-tan(beta), sin(mi));

    if(mi >= mi_s || mi <= 0){
         if(mi < 0) mi += 2.0*G_PI;
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
         if ((beta >= 0 && yaw < yaw_nom) || (beta < 0 && yaw > yaw_nom)) yaw = yaw_nom;

    }else yaw = yaw_nom;
    
    _yaw2ijk(satdata, yaw, i, j, k);
}
    
// midnight maneuver
// -------------------------------------------
void t_gpppmodel::_midnight_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta0 = atan2(MUDOT_GPS, R);
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    
    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }
    
    if (beta < 0) R *= -1;
    
    double mi_s = - sqrt( beta0 * fabs(beta) - beta*beta );
    double yaw;
    double yaw_nom = atan2(-tan(beta), sin(mi));
    
    if(mi >= mi_s){
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
         if ((beta >= 0 && yaw > yaw_nom) || (beta < 0 && yaw < yaw_nom)) yaw = yaw_nom;
    }else yaw = yaw_nom;
    
    _yaw2ijk(satdata, yaw, i, j, k);       
}

// midnight maneuver
// -------------------------------------------
void t_gpppmodel::_midnight_turn_GPSIIA(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    double mi_s = - sqrt(EPS0_GPS*EPS0_GPS - beta*beta);
    double mi_e = - mi_s;         

    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }
    
    if (beta < 0) R *= -1;
    
    double yaw;
    
    if (mi_s <= mi && mi < mi_e) {
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GPS;
    }else if(mi_e <= mi && mi < mi_e + POST_SHADOW*MUDOT_GPS) satdata.setecl(true);
    else yaw = atan2(-tan(beta), sin(mi));
    
    _yaw2ijk(satdata, yaw, i, j, k);       
}   
    
// midnight maneuver for GPS Block IIF
// -------------------------------------------
void t_gpppmodel::_midnight_turn_GPSIIF(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double  tan_beta = tan(beta);
    double mi   = satdata.orb_angle();
    
    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) beta *= -1;
    }
    
    if (beta < 0) R *= -1;
    
    double mi_s = -acos(cos(EPS0_GPS) / cos(beta)); 
    double mi_e = - mi_s;
    double sin_mi_s = sin(mi_s);
    double mi_f = MUDOT_GPS*(atan2(-tan_beta, -sin_mi_s) - atan2(-tan_beta, sin_mi_s)) / R + mi_s;

    double yaw;
    
    if (mi_s <= mi && mi < mi_f) {
         yaw = atan2(-tan_beta, sin_mi_s) + R*(mi - mi_s) / MUDOT_GPS;
         _yaw2ijk(satdata, yaw, i, j, k);     
      }
      else if (mi_f <= mi && mi < mi_e) {
         yaw = atan2(-tan_beta, -sin_mi_s);
          _yaw2ijk(satdata, yaw, i, j, k);       
      }else _ysm(satdata, i, j, k);
    

}
    
// midnight maneuver for GLONASS-M
// -------------------------------------------
void t_gpppmodel::_midnight_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double  tan_beta = tan(beta);
    double mi   = satdata.orb_angle();
    
    if (beta < 0) R *= -1;
    
    double mi_s = -acos(cos(EPS0_GLO) / cos(beta)); 
    double mi_e = - mi_s;
    double sin_mi_s = sin(mi_s);
    double mi_f = MUDOT_GLO*(atan2(-tan_beta, -sin_mi_s) - atan2(-tan_beta, sin_mi_s)) / R + mi_s;

    double yaw;
    
    if (mi_s <= mi && mi < mi_f) {
         yaw = atan2(-tan_beta, sin_mi_s) + R*(mi - mi_s) / MUDOT_GLO;
      }
      else if (mi_f <= mi && mi < mi_e) {
         yaw = atan2(-tan_beta, -sin_mi_s);
      }      
    
    _yaw2ijk(satdata, yaw, i, j, k);
}   

// noon maneuver
// -------------------------------------------
void t_gpppmodel::_noon_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    
    if (beta >= 0) R = -R;
    
    double mi_s, mi_e, sin_mi_s, B, yaw;
    int c = 0;

    for(mi_s = 178.6*D2R; c < 4; c++) {
         sin_mi_s = sin(mi_s);
         B = -beta*cos(mi_s) / (beta*beta + sin_mi_s*sin_mi_s);
         mi_s = (atan(-beta / sin_mi_s) + mi_s*B + G_PI*R / MUDOT_GLO - G_PI / 2.0) / (R / MUDOT_GLO + B);
    }
    if(beta >= 0) mi_s = 2.0*G_PI - mi_s;
      mi_e = 2.0*G_PI - mi_s;

      if (mi_s <= mi && mi < mi_e) {
         yaw = atan2(-tan(beta), sin(mi_s)) + R*(mi - mi_s) / MUDOT_GLO;
      }
    
    _yaw2ijk(satdata, yaw, i, j, k);    
}
    
// Attitude modelling for GPS Block IIR-M
void t_gpppmodel::_attitude_GPSIIRM(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    _attitude_GPSIIR(satdata, antype, i, j, k);
}
    
// Attitude modelling for GPS Block IIF
void t_gpppmodel::_attitude_GPSIIF(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R_noon = 0.11*D2R;             // maximal yaw hardware rate during noon turn
    const double R_midn = 0.06*D2R;             // maximal yaw hardware rate during midnight turn
    double beta0 = atan2(MUDOT_GPS, R_noon);

    if(fabs(satdata.orb_angle()) > G_PI/2 && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn(satdata, R_noon, i, j, k);
    }else if(fabs(satdata.orb_angle()) < EPS0_GPS && fabs(satdata.beta()) < EPS0_GPS){
         // midnight maneuver
         _midnight_turn_GPSIIF(satdata, R_midn, i, j, k);
    }else{
         _ysm(satdata, i, j, k); 
    }  
}
    
// Attitude modelling for GLONASS
void t_gpppmodel::_attitude_GLO(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    const double R = 0.25*D2R;             // maximal yaw hardware rate
    double beta0 = 2.0*D2R;
    double mi = satdata.orb_angle();
    double mi_s = 176.8*D2R;
    
    if(mi > mi_s && mi < 2.0*G_PI - mi_s && fabs(satdata.beta()) < beta0){
         // noon maneuver
         _noon_turn_GLOM(satdata, R, i, j, k);
    }else if(fabs(satdata.orb_angle()) < EPS0_GLO && fabs(satdata.beta()) < EPS0_GLO){
         // midnight maneuver
         _midnight_turn_GLOM(satdata, R, i, j, k);
    }else{
         _ysm(satdata, i, j, k); 
    }  
}   
    
// noon maneuver for Galileo IOV
// -------------------------------------------
void t_gpppmodel::_noon_turn_GAL1(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();
    
    double beta_y = 2*D2R;
    double beta_x = 15*D2R;
    
    double Sx =  sin(mi)*cos(beta);
    double Sy = -sin(beta);

    double sinby = sin(beta_y);
    double sinbx = sin(beta_x);   
    if(Sy<0) sinby *= -1;
    
    double Shy = (sinby + Sy)/2 + cos(G_PI*fabs(Sx)/sinbx) * (sinby - Sy)/2;
          
    double yaw = atan2(Shy, Sx);
    
    // test of beta sign change
    auto itSAT = _last_beta.find(satdata.sat());
    if(itSAT != _last_beta.end()){
         if(itSAT->second * beta < 0) yaw *= -1;
    }
    
    _yaw2ijk(satdata, yaw, i, j, k);
}
    
// Attitude modelling for Galileo IOV
void t_gpppmodel::_attitude_GAL1(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{   
   double beta0 = 2.0*D2R;
    double mi_n0 = (180-15)*D2R;
    double mi_m0 = 15*D2R;
    
   if(fabs(satdata.orb_angle()) > mi_n0 && fabs(satdata.beta()) < beta0){
     // noon maneuver
      _noon_turn_GAL1(satdata, i, j, k);
   }else if(fabs(satdata.orb_angle()) < mi_m0 && fabs(satdata.beta()) < beta0){
      // midnight maneuver
      _noon_turn_GAL1(satdata, i, j, k);
   }else{
         _ysm(satdata, i, j, k);
    }

//   i = -1 * i;      // X away from the Sum
//   j = -1 * j;  
//   satdata.yaw(satdata.yaw() + G_PI);
}

// noon maneuver for Galileo FOC
// -------------------------------------------
void t_gpppmodel::_noon_turn_GAL2(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{   
    double beta = satdata.beta();
    double mi   = satdata.orb_angle();

    // test of beta sign change
    auto it = _last_beta.find(satdata.sat());
    if(it != _last_beta.end()){
         if(it->second * beta < 0) beta *= -1;
    } 
    
    double init = G_PI/2;
    if(beta > 0) init *= -1;
    
    double yaw_s;
    t_gtime epo_s;
    t_gtime epo = satdata.epoch();
    
    auto itYAW = _last_yaw.find(satdata.sat());
    if(itYAW != _last_yaw.end()){
         yaw_s = itYAW->second;
    }else{
         yaw_s = atan2(-tan(beta), sin(mi));
         _last_yaw[satdata.sat()] = yaw_s;
    }
    
    auto itEPO = _last_epo.find(satdata.sat());
    if(itEPO != _last_epo.end()){
         epo_s = itEPO->second;
    }else{
         epo_s = epo;
         _last_epo[satdata.sat()] = epo_s;
    }  
    
   double yaw = init + (yaw_s - init) * cos((2.0*G_PI/5656)*(epo - epo_s));
       
    _yaw2ijk(satdata, yaw, i, j, k);
}   

// Attitude modelling for Galileo FOC
void t_gpppmodel::_attitude_GAL2(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{   
   double beta0 = 4.1*D2R;
    double mi_n0   = (180-10)*D2R;
    double mi_m0   = 10*D2R;
    
   if(fabs(satdata.orb_angle()) > mi_n0 && fabs(satdata.beta()) < beta0){
     // noon maneuver
      _noon_turn_GAL2(satdata, i, j, k);
   }else if(fabs(satdata.orb_angle()) < mi_m0 && fabs(satdata.beta()) < beta0){
      // midnight maneuver
      _noon_turn_GAL2(satdata, i, j, k);
   }else{
      _ysm(satdata, i, j, k);
         _last_yaw[satdata.sat()] = satdata.yaw();
         _last_epo[satdata.sat()] = satdata.epoch();        
   }

//   i = -1 * i;      // X away from the Sum
//   j = -1 * j;  
//   satdata.yaw(satdata.yaw() + G_PI);
}
       
// Attitude modelling for BeiDou
void t_gpppmodel::_attitude_BDS(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{

    if( t_gsys::bds_geo(satdata.sat()) ){
         _onm(satdata, i, j, k);
    }else{
         if( fabs(satdata.beta()) <= 4.0*D2R ){
             _onm(satdata, i, j, k);       
         }else _ysm(satdata, i, j, k);
    }  
}
    
// Attitude modelling for QZSS
void t_gpppmodel::_attitude_QZS(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{   
    if( fabs(satdata.beta()) <= 20*D2R ){
         _onm(satdata, i, j, k);  // i oposite velocity direction
         i = -1 * i;
         j = -1 * j;
    }else _ysm(satdata, i, j, k);    
}
       
// Calculate satellite-fixed vectors from yaw angle    
void t_gpppmodel::_yaw2ijk(t_gsatdata& satdata, double& yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k)
{   
  satdata.yaw(yaw);   // store yaw angle

  ColumnVector satcrd = satdata.satcrd().crd_cvect();
  ColumnVector satvel = satdata.satvel().crd_cvect();       

  if(satcrd.NormFrobenius() == 0 || satcrd.NormFrobenius() == 0) return;
    
  // ITRF -> ICRF velocity
  satvel(1) -= OMEGA*satcrd(2);
  satvel(2) += OMEGA*satcrd(1);
   
  ColumnVector n = crossproduct(satcrd, satvel);   
    
  // Satelite-Earth unit vector
  k = -satcrd;
  k /= k.NormFrobenius();    

  ColumnVector ex = crossproduct(n, satcrd);
  
  ex /= ex.NormFrobenius();
  n  /=  n.NormFrobenius();
    
       
  double cosy = cos(yaw);
  double siny = sin(yaw);
  for (int r = 1; r <= 3; r++) {
    i(r) = -siny*n(r) + cosy*ex(r);
    j(r) = -cosy*n(r) - siny*ex(r);
  }     
}
    
} // namespace
