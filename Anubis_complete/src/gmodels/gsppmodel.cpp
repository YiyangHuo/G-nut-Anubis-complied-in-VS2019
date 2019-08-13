
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


#include "gmodels/gsppmodel.h"
#include "gmodels/gtropoblind.h"
#include "gdata/gpppdata.h"
#include "gprod/gprodion.h"

namespace gnut {      
   
// Constructors
// -------------------
t_gsppmodel::t_gsppmodel()
{
  gtrace("t_gsppmodel::t_gsppmodel");
   
  _tropoModel = make_shared<t_gtropo>();
  _gallbias = 0;
}
   
t_gsppmodel::t_gsppmodel(string site, t_glog* glog, t_gsetbase* settings)
: _observ(IONO_FREE)
{
  gtrace("t_gsppmodel::t_gsppmodel");
   
  _tropoModel = make_shared<t_gtropo>();
  _gallbias = 0;
  
  _settings = settings;
  _site = site;
  _log = glog;
  _phase = false;
  
  set<GSYS> systems = GNSS_SUPPORTED();
  for(set<GSYS>::iterator it = systems.begin(); it != systems.end(); it++){
    _maxres_C[*it] = dynamic_cast<t_gsetgnss*>(_settings)->maxres_C(*it);
    _maxres_L[*it] = dynamic_cast<t_gsetgnss*>(_settings)->maxres_L(*it);
  }

  _maxres_norm      = dynamic_cast<t_gsetproc*>(_settings)->max_res_norm();
  _tropo_mf         = dynamic_cast<t_gsetproc*>(_settings)->tropo_mf();      
  _trpModStr        = dynamic_cast<t_gsetproc*>(_settings)->tropo_model();
  
  _resid_type       = dynamic_cast<t_gsetproc*>(_settings)->residuals();
  _observ           = dynamic_cast<t_gsetproc*>(_settings)->obs_combin();
  
  if      (_trpModStr == SAASTAMOINEN)  _tropoModel = make_shared<t_saast>();
  else if (_trpModStr == DAVIS)         _tropoModel = make_shared<t_davis>();
  else if (_trpModStr == HOPFIELD)      _tropoModel = make_shared<t_hopf>();
  else if (_trpModStr == MOPS)          _tropoModel = make_shared<t_blindmops>();
// else if (_trpModStr == GPTW)          _tropoModel = make_shared<t_blindgpt2w>();
// else if (_trpModStr == GPT2W)         _tropoModel = make_shared<t_blindgpt2w>();
// else if (_trpModStr == GAL27)         _tropoModel = make_shared<t_blindgal27>();
// else if (_trpModStr == GALTROPO27)    _tropoModel = make_shared<t_blindgal27>();
}

// Destructor
// ---------------------
t_gsppmodel::~t_gsppmodel()
{
  gtrace("t_gsppmodel::~t_gsppmodel");
}


// Outliers detection based on chi2 testing of normalized residuals
// ----------
int t_gsppmodel::outlierDetect_chi(      vector<t_gsatdata>& data,
                          SymmetricMatrix& Qx,
                    const SymmetricMatrix& Qsav, 
                    const ColumnVector& v)
{
  gtrace("t_gsppmodel::outlierDetect_chi");
   
  vector<t_gsatdata>::iterator it;
  vector<t_gsatdata>::iterator itMaxV;
   
  double maxV = 0.0;
   
  int ii = 0;
  for( it = data.begin(); it != data.end(); it++ ){ 
    ++ii;

    if( maxV == 0.0 || v(ii)*v(ii) > maxV ){
      
      maxV = v(ii)*v(ii);
      itMaxV = it;
    }
      
    if( _phase ){ 
      ++ii;
      if( maxV == 0.0 || v(ii)*v(ii) > maxV ){
        
        maxV = v(ii)*v(ii);
        itMaxV = it;
      }   
    }      
  }
  
  if( maxV > 5.024 ){ // chi2(100) = 2.706; chi2(050) = 3.841; chi2(025) = 5.024; chi2(010) = 6.635; chi2(005) = 7.879;

    if( _log)
    _log->comment(2, "gpppmodel", _site + " outlier " + itMaxV->sat()
                  + " size:" + int2str(data.size())
                  + " v_norm: "   + dbl2str(maxV)
                  + " " + itMaxV->epoch().str_hms() );
    data.erase(itMaxV);
    Qx = Qsav;
    return 1;
  }
  
   return 0;
}
   
// Outliers detection
// ----------
int t_gsppmodel::outlierDetect(      vector<t_gsatdata>& data,
                      SymmetricMatrix& Qx,
                const SymmetricMatrix& Qsav, 
                const ColumnVector& v)
{
   gtrace("t_gsppmodel::outlierDetect");   
      
   vector<t_gsatdata>::iterator itMaxVcodeGPS;
   vector<t_gsatdata>::iterator itMaxVcodeGLO;
   vector<t_gsatdata>::iterator itMaxVcodeGAL;
   vector<t_gsatdata>::iterator itMaxVcodeBDS;
   vector<t_gsatdata>::iterator itMaxVcodeQZS;   
   
   vector<t_gsatdata>::iterator itMaxVphaseGPS;
   vector<t_gsatdata>::iterator itMaxVphaseGLO;
   vector<t_gsatdata>::iterator itMaxVphaseGAL;
   vector<t_gsatdata>::iterator itMaxVphaseBDS;
   vector<t_gsatdata>::iterator itMaxVphaseQZS;

   double maxVcodeGPS, maxVcodeGLO, maxVcodeGAL, maxVcodeBDS, maxVcodeQZS;
   maxVcodeGPS = maxVcodeGLO = maxVcodeGAL = maxVcodeBDS = maxVcodeQZS = 0.0;
   
   double maxVphaseGPS, maxVphaseGLO, maxVphaseGAL, maxVphaseBDS, maxVphaseQZS;
   maxVphaseGPS = maxVphaseGLO = maxVphaseGAL = maxVphaseBDS = maxVphaseQZS = 0.0;

   // deviding multi-freq residual vector into single-freq vectors
   vector<ColumnVector> v_frqs = _devide_res(v);
   
   // find maximal code/phase residuals for individual GNSS   
   for(unsigned int i = 0; i < v_frqs.size(); i++){   
     double maxres = _maxres(v_frqs[i], GPS, false, data, itMaxVcodeGPS);
     if(maxres > maxVcodeGPS) maxVcodeGPS = maxres;
            maxres = _maxres(v_frqs[i], GLO, false, data, itMaxVcodeGLO);
     if(maxres > maxVcodeGLO) maxVcodeGLO = maxres;
            maxres = _maxres(v_frqs[i], GAL, false, data, itMaxVcodeGAL);
     if(maxres > maxVcodeGAL) maxVcodeGAL = maxres;
            maxres = _maxres(v_frqs[i], BDS, false, data, itMaxVcodeBDS);
     if(maxres > maxVcodeBDS) maxVcodeBDS = maxres;
            maxres = _maxres(v_frqs[i], QZS, false, data, itMaxVcodeQZS);   
     if(maxres > maxVcodeQZS) maxVcodeQZS = maxres;
   
            maxres = _maxres(v_frqs[i], GPS, true,  data, itMaxVphaseGPS);
     if(maxres > maxVphaseGPS) maxVphaseGPS = maxres;
             maxres = _maxres(v_frqs[i], GLO, true,  data, itMaxVphaseGLO);
     if(maxres > maxVphaseGLO) maxVphaseGLO = maxres;
             maxres = _maxres(v_frqs[i], GAL, true,  data, itMaxVphaseGAL);
     if(maxres > maxVphaseGAL) maxVphaseGAL = maxres;
            maxres = _maxres(v_frqs[i], BDS, true,  data, itMaxVphaseBDS);
     if(maxres > maxVphaseBDS) maxVphaseBDS = maxres;
            maxres = _maxres(v_frqs[i], QZS, true,  data, itMaxVphaseQZS);
     if(maxres > maxVphaseQZS) maxVphaseQZS = maxres;
   }

#ifdef DEBUG
    cout << "Max res range: " << maxVcodeGPS << " " << itMaxVcodeGPS->sat() << endl;
    cout << "Max res phase: " << maxVphaseGPS << " "<< itMaxVphaseGPS->sat() << endl;
#endif
   //Only detect the outliers for the constellations with maximum outliers
   //The GLONASS is set as the basic reference
   double maxvc = maxVcodeGLO, maxvp = maxVphaseGLO;
   GSYS   maxsys = GLO;
   if (maxVcodeGPS > maxvc || maxVphaseGPS > maxvp) { maxvc = maxVcodeGPS; maxvp = maxVphaseGPS; maxsys = GPS; }
   if (maxVcodeGAL > maxvc || maxVcodeGAL > maxvp)  { maxvc = maxVcodeGAL; maxvp = maxVcodeGAL;  maxsys = GAL; }
   if (maxVcodeBDS > maxvc || maxVphaseBDS > maxvp) { maxvc = maxVcodeBDS; maxvp = maxVphaseBDS; maxsys = BDS; }
   if (maxVcodeQZS > maxvc || maxVphaseQZS > maxvp) { maxvc = maxVcodeQZS; maxvp = maxVphaseQZS; maxsys = QZS; }


   if(maxsys==GLO && _check_outl(false, maxVcodeGLO,  itMaxVcodeGLO, data)) { data.erase(itMaxVcodeGLO);  Qx = Qsav; return 1; }
   if(maxsys==GLO && _check_outl(true, maxVphaseGLO, itMaxVphaseGLO, data) ) { data.erase(itMaxVphaseGLO); Qx = Qsav; return 1; }
   
   if(maxsys==GPS && _check_outl(false, maxVcodeGPS,  itMaxVcodeGPS, data) ) { data.erase(itMaxVcodeGPS);  Qx = Qsav; return 1; }
   if(maxsys==GPS && _check_outl(true, maxVphaseGPS, itMaxVphaseGPS, data) ) { data.erase(itMaxVphaseGPS); Qx = Qsav; return 1; }
   
   if(maxsys==GAL &&  _check_outl(false, maxVcodeGAL,  itMaxVcodeGAL, data) ) { data.erase(itMaxVcodeGAL);  Qx = Qsav; return 1; }
   if(maxsys==GAL &&  _check_outl(true, maxVphaseGAL, itMaxVphaseGAL, data) ) { data.erase(itMaxVphaseGAL); Qx = Qsav; return 1; }
   
   if(maxsys==BDS &&  _check_outl(false, maxVcodeBDS,  itMaxVcodeBDS, data) ) { data.erase(itMaxVcodeBDS);  Qx = Qsav; return 1; }
   if(maxsys==BDS &&  _check_outl(true, maxVphaseBDS, itMaxVphaseBDS, data) ) { data.erase(itMaxVphaseBDS); Qx = Qsav; return 1; }
   
   if(maxsys==QZS &&  _check_outl(false, maxVcodeQZS,  itMaxVcodeQZS, data) ) { data.erase(itMaxVcodeQZS);  Qx = Qsav; return 1; }
   if(maxsys==QZS &&  _check_outl(true, maxVphaseQZS, itMaxVphaseQZS, data) ) { data.erase(itMaxVphaseQZS); Qx = Qsav; return 1; }      
   
   return 0;
}

   
// Outliers detection
// ----------
int t_gsppmodel::outlierDetect( vector<t_gsatdata>& data,
                             SymmetricMatrix& Qx,
                       const SymmetricMatrix& Qsav)
{
  gtrace("t_gsppmodel::outlierDetect");   
      
  vector<t_gsatdata>::iterator itMaxVcodeNORM  = data.end();
  vector<t_gsatdata>::iterator itMaxVphaseNORM = data.end();
   
  vector<t_gsatdata>::iterator itMaxVcodeORIG  = data.end();
  vector<t_gsatdata>::iterator itMaxVphaseORIG = data.end();

  vector<t_gsatdata>::iterator itDataErase = data.end();

  double maxVcodeNORM  = 0.0;   
  double maxVphaseNORM = 0.0;
   
  double maxVcodeORIG  = 0.0;   
  double maxVphaseORIG = 0.0;   
  
  // find maximal code/phase residuals
  maxVcodeNORM  = _maxres(false, data, itMaxVcodeNORM,  RES_NORM);
  maxVphaseNORM = _maxres(true,  data, itMaxVphaseNORM, RES_NORM);

  maxVcodeORIG  = _maxres(false, data, itMaxVcodeORIG,  RES_ORIG);
  maxVphaseORIG = _maxres(true,  data, itMaxVphaseORIG, RES_ORIG);   

#ifdef DEBUG
  if(itMaxVcodeNORM  != data.end()) cout << "Max res range norm: " << fixed << setprecision(3) << maxVcodeNORM  << " " << itMaxVcodeNORM->sat()  << endl;
  if(itMaxVphaseNORM != data.end()) cout << "Max res phase norm: " << fixed << setprecision(3) << maxVphaseNORM << " " << itMaxVphaseNORM->sat() << endl;
  if(itMaxVcodeORIG  != data.end()) cout << "Max res range orig: " << fixed << setprecision(3) << maxVcodeORIG  << " " << itMaxVcodeORIG->sat()  << endl;
  if(itMaxVphaseORIG != data.end()) cout << "Max res phase orig: " << fixed << setprecision(3) << maxVphaseORIG << " " << itMaxVphaseORIG->sat() << endl;
  int ooo; cin >> ooo;
#endif

  if( _check_outl(true,  maxVphaseNORM, itMaxVphaseNORM, maxVphaseORIG, itMaxVphaseORIG, itDataErase, data) ) { data.erase(itDataErase); Qx = Qsav; return 1; }   
  if( _check_outl(false, maxVcodeNORM,  itMaxVcodeNORM,  maxVcodeORIG,  itMaxVcodeORIG,  itDataErase, data) ) { data.erase(itDataErase); Qx = Qsav; return 1; }
   
  return 0;
}   
   
// model for computed range value 
// ----------
double t_gsppmodel::cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata& gsatdata, t_gobs& gobs, bool com)
{   
  gtrace("t_gsppmodel::cmpObs");   
   
  // Cartesian coordinates to ellipsodial coordinates
  t_gtriple xyz, ell;
  ColumnVector cRec(3);

  if ( param.getCrdParam(_site, xyz) > 0 ){
    cRec = xyz.crd_cvect();
  }else{
    xyz  = _grec->crd_arp(epoch);
    cRec = xyz.crd_cvect();
  }      
  xyz2ell(xyz, ell, false);
  
  t_gtriple satcrd = gsatdata.satcrd();
  ColumnVector cSat = satcrd.crd_cvect();          
  
  string sat  = gsatdata.sat(); 
  t_gtime epo = gsatdata.epoch();
  
  // Tropospheric wet delay correction
  double trpDelay = 0;
  trpDelay = tropoDelay(epoch, param, ell, gsatdata);   
  if(fabs(trpDelay) > 50) return -1;

  // Receiver clock correction 
  double clkRec = 0.0;
  int i;
  i = param.getParam(_site, t_gpar::CLK, "");
  if ( i >= 0 ) clkRec = param[i].value();
  else
  if( _log && _log->verb() >= 1) 
    _log->comment(1, "gppp", _site + " ! warning:  Receiver Clock is not included in parameters!");   
  
  // Inter frequency clocks bias
  double ifcb = 0.0;
  i = param.getParam(_site, t_gpar::IFCB_F3, gsatdata.sat());   
  if ( i >= 0 && t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_3 ) {
    ifcb = -param[i].value();
  }
  i = param.getParam(_site, t_gpar::IFCB_F4, gsatdata.sat());   
  if ( i >= 0 && t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_4 ) {
    ifcb = -param[i].value();
  }
  i = param.getParam(_site, t_gpar::IFCB_F5, gsatdata.sat());   
  if ( i >= 0 && t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_5 ) {
    ifcb = -param[i].value();
  }  
  
  // Inter frequency code bias FREQ_3
  double ifb = 0.0;
  i = param.getParam(_site, t_gpar::IFB_C3, "");
  if ( i >= 0 && gobs.is_code() && t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_3) {
    ifb = param[i].value();
    //cout << epoch.str_ymdhms() << " " << gsatdata.sat() << " " << fixed << setprecision(3) << setw(9) << ifb_c3 << endl;
  }
  
  // Inter frequency code bias FREQ_4
  i = param.getParam(_site, t_gpar::IFB_C4, "");   
  if ( i >= 0 && gobs.is_code() && t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_4) {
    ifb = param[i].value();
    //cout << epoch.str_ymdhms() << " " << gsatdata.sat() << " " << fixed << setprecision(3) << setw(9) << ifb_c4 << endl;
  }
  
  // Inter frequency code bias FREQ_5
  i = param.getParam(_site, t_gpar::IFB_C5, "");   
  if ( i >= 0 && gobs.is_code() && t_gsys::band2freq(gsatdata.gsys(), gobs.band()) == FREQ_5) {
    ifb = param[i].value();
    //cout << epoch.str_ymdhms() << " " << gsatdata.sat() << " " << fixed << setprecision(3) << setw(9) << ifb_c5 << endl;
  }  

  // GLONASS system time offset
  double glonass_offset = 0;
  if (gsatdata.gsys() == GLO){
    int i;
    i = param.getParam(_site, t_gpar::GLO_ISB, "");
    if ( i >= 0 ) {
      glonass_offset = param[i].value();
      //cout << "GLO Offset: " << sat << " " << glonass_offset << " " << gsatdata.epoch().str_hms() << endl;
    }     
    
  }
  
  // Galileo system time offset
  double galileo_offset = 0;
  if (gsatdata.gsys() == GAL){
    int i;
    i = param.getParam(_site, t_gpar::GAL_ISB, "");
    if ( i >= 0 ) { galileo_offset = param[i].value();
      //cout << "GAL Offset: " << sat << " " << galileo_offset << " " << gsatdata.epoch().str_hms() << endl;       
    }
     
  }
   
  // BaiDou system time offset
  double beidou_offset = 0;
  if (gsatdata.gsys() == BDS){
    int i;
    i = param.getParam(_site, t_gpar::BDS_ISB, "");
    if ( i >= 0 ) {beidou_offset = param[i].value();
      //cout << "BDS ISB: " << sat << " " << param[i].value() << " " << gsatdata.epoch().str_hms() << endl;
    }
  }
  
  // QZSS system time offset
  double qzss_offset = 0;
  if (gsatdata.gsys() == QZS){
    int i;
    i = param.getParam(_site, t_gpar::QZS_ISB, "");
    if ( i >= 0 ) qzss_offset = param[i].value();
  }
  
  // ionosphere delay
  double ionDelay = 0.0;
  if(_observ != IONO_FREE) {
    ionDelay = ionoDelay(epoch, param, ell, gsatdata, gobs);
  }
  
  // Code bias
  double cbias = 0.0;
  GSYS sys   = gsatdata.gsys();
  double fk = gsatdata.frequency(gobs.band());
  FREQ_SEQ freq = t_gsys::band2freq(gsatdata.gsys(), gobs.band());  
  
  GOBS g;
  if(gobs.attr() == ATTR) g = gsatdata.id_range(gobs.band());  // automatic GOBS selection
  else g = gobs.gobs();                                       // specific GOBS    
  
//cout << _site << " DCB: " << sat << " " << epo.str_ymdhms() << " C1W-C2W " << corrDCB/CLIGHT*1e9 << endl;
//cout << "DCB: " << sat << " " << epo.str_ymdhms() << " C1C-C5Q " << corrDCB_15 << endl;  
  if(gobs.is_code() && _observ != IONO_FREE){
    // find DCB
    double corrDCB = 0.0;
//    double corrDCB_15 = 0.0;
    if(_gallbias) {
      if(sys == GPS) {
        if(g < 1000) {
          corrDCB = _gallbias->get(epo, sat, C1W, C2W);
          //corrDCB_15 = _gallbias->get(epo, sat, C1W, C5Q);
        }else  {          
          corrDCB = _gallbias->get(epo, sat, P1, P2);
        }
      }else if(sys == GAL) {
        if(g < 1000) corrDCB = _gallbias->get(epo, sat, C1X, C5X);
        else         corrDCB = _gallbias->get(epo, sat, C1, C2);
      }
    }

    // Satellite GPS code bias P1-P2
    if(sys == GPS && (freq == FREQ_1 || freq == FREQ_2)) {
      double f1 = G01_F;      
      double alfa = (f1*f1) / (fk*fk);
      double beta = (G02_F*G02_F)/(G01_F*G01_F - G02_F*G02_F);      
      cbias = -alfa * beta * corrDCB; // P1-P2
    }
    
    // // Satellite GPS code bias P1-P5
    // if(sys == GPS && freq == FREQ_3) {
    //   double beta = (G01_F*G01_F)/(G01_F*G01_F - G05_F*G05_F);
    //   cbias = - beta * corrDCB_15; // P1-P5
    // }

    // Satellite GAL code bias E1-E5a    
    if(sys == GAL && (freq == FREQ_1 || freq == FREQ_2)) {
      double f1 = E01_F;      
      double alfa = (f1*f1) / (fk*fk);
      double beta = (E05_F*E05_F)/(E01_F*E01_F - E05_F*E05_F);
      cbias = alfa * beta * corrDCB; // E1-E5      
    }    

    // Receiver GPS code bias P1-P2
    int idcb;
    idcb = param.getParam(_site, t_gpar::P1P2G_REC, "");
    if ( idcb >= 0 && sys == GPS && (freq == FREQ_1 || freq == FREQ_2) ){
      double f1 = G01_F;      
      double alfa = (f1*f1) / (fk*fk);
      double beta = (G02_F*G02_F)/(G01_F*G01_F - G02_F*G02_F);
      cbias -= alfa * beta * param[idcb].value();
    }
    
    // Receiver GAL code bias E1-E5a
    idcb = param.getParam(_site, t_gpar::P1P2E_REC, "");
    if ( idcb >= 0 && sys == GAL && (freq == FREQ_1 || freq == FREQ_2) ){
      double f1 = E01_F;
      double alfa = (f1*f1) / (fk*fk);
      double beta = (E05_F*E05_F)/(E01_F*E01_F - E05_F*E05_F);
      cbias -= alfa * beta * param[idcb].value();
    }    
  }
   

#ifdef DEBUG
   cout << gsatdata.epoch().str_hms() <<  " " 
        << gsatdata.sat() << " " << gobs2str(g)
       << fixed << setprecision(3)
//        << "  rho: "       << gsatdata.rho()
        << "  ele: "       << setw(6) << gsatdata.ele()*180.0/G_PI
//        << "  H:   "       << setw(6) << ell[2]
        << "  rec clk: "   << setw(6) << clkRec
        << "  sat clk: "   << setw(6) << gsatdata.clk()
//        << "  ifcb:  "   << setw(6) << ifcb
        << "  Tropdelay: " << setw(6) << trpDelay
        << "  Ionodelay: " << setw(6) << ionDelay
//        << "  GLO_ISB " << setw(6) << glonass_offset
//        << "  GAL_ISB " << setw(6) << galileo_offset
//        << "  BDS_ISB " << setw(6) << beidou_offset
//        << "  QZS_ISB " << setw(6) << qzss_offset
        << endl;
//   int ooo; cin >> ooo;
#endif          

   
   // Return value
  return gsatdata.rho() + 
         clkRec         - 
         gsatdata.clk() +
         trpDelay       +
         ifcb           +
         ifb            +
         glonass_offset +
         galileo_offset +
         beidou_offset  +
         qzss_offset    +
         ionDelay      +
         cbias;
}

// Compute troposperic delay
// -----------
double t_gsppmodel::tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple ell, t_gsatdata& satdata)
{        
  gtrace("t_gsppmodel::tropoDelay");   

  if (_tropoModel == 0){
    if( _log ) _log->comment(0, "gppp","Tropo Model setting is not correct. Default used! Check config.");
    else                cerr << "gppp - Tropo Model setting is not correct. Default used! Check config.\n";
    _tropoModel = make_shared<t_saast>();
  }
   
  double ele = satdata.ele();
  
  double delay = 0.0;
  double zwd   = 0.0;
  double zhd   = 0.0;
   
  int i;
  i = param.getParam(_site, t_gpar::TRP, "");
  if ( i >= 0 ){
    zwd = param[i].value();
    zhd = param[i].apriori();
  }else{
    if(_tropoModel != 0){
      zwd = _tropoModel->getZWD(ell, epoch);
      zhd = _tropoModel->getZHD(ell, epoch);
    }
  }
   
  if ( _tropo_mf == GMF ){
    double gmfh, gmfw, dgmfh, dgmfw;
    t_gmf mf;
    mf.gmf( epoch.mjd(), ell[0], ell[1], ell[2], G_PI/2.0-ele,
            gmfh, gmfw, dgmfh, dgmfw );
    //      delay = gmfh * _tropoModel->getZHD(ell, epoch) + gmfw * zwd;
    delay = gmfh * zhd + gmfw * zwd;
      
#ifdef DEBUG
    cout << epoch.str("EPOCH: %H:%M:%S") << endl << fixed << setprecision(3);
    cout << "Ell:" << ell[0] << " " << ell[1] << " " << ell[2]
    << " Hydrostatic part: " << zhd // _tropoModel->getZHD(ell, epoch)
    << " Wet part: "        << zwd
    << " gmfh: " << gmfh
    << " gmfw: " << gmfw
    << " Delay: " << delay << endl << endl;
    int ooo; cin >> ooo;
#endif
  }else if ( _tropo_mf == COSZ ){
    double mf = 1/sin(ele);
    //      delay = 1/sin(ele) * _tropoModel->getZHD(ell, epoch) +1/sin(ele) * zwd;
    delay =  mf * zhd + mf * zwd;
  }
  
  return delay;
}   

// Compute ionospheic delay
// -----------
double t_gsppmodel::ionoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& gsatdata, t_gobs& gobs)
{        
  gtrace("t_gsppmodel::ionoDelay");

  double iono_param = 0.0;
  double iono_model = 0.0;

  double mf = 1 / sqrt(1.0 - pow(R_SPHERE / (R_SPHERE + 450000.0) * sin(G_PI / 2.0 - gsatdata.ele()), 2));
  double f1 = G01_F;
  double fk = gsatdata.frequency(gobs.band());
  double alfa = 0.0;
  if (gobs.is_phase()) {
    alfa = -(f1 * f1) / (fk * fk);
  }
  if (gobs.is_code()) {
    alfa = (f1 * f1) / (fk * fk);
  }
  
  // ionosphere slant delay parameter
  int i = param.getParam(_site, t_gpar::SION, gsatdata.sat());
  if ( i >= 0 ) {        
    iono_param = alfa * param[i].value();
  }
  
  // ionosphere vertical delay paremeter
  i = param.getParam(_site, t_gpar::VION, gsatdata.sat());
  if ( i >= 0 ) {        
    iono_param = alfa * mf * param[i].value();
  }

    // ionosphere model
  if (_gallion && i < 0) {
    t_gtriple ipp_ell(0.0, 0.0, 0.0);
    ell2ipp(gsatdata, site_ell, ipp_ell);
    iono_model = _gallion->iono(ipp_ell[0] * R2D, ipp_ell[1] * R2D, epoch);
    iono_model = ( iono_model * 40.28 * 1e16 ) / (G01_F*G01_F);
    
    iono_model *= alfa * mf;    
  }
  
  #ifdef DEBUG
    cout << "IONO delay: " << gsatdata.sat() << " " << epoch.str_ymdhms() << " model = " << iono_model << " iono_par = "  << iono_param << endl;
  //  int ooo; cin >> ooo;
  #endif

  return iono_model + iono_param;
}

// Find maximal residual
double t_gsppmodel::_maxres(const ColumnVector& v, GSYS gs, bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA)
{      
  unsigned int inc  = 2;   
   
  if(v.Nrows() == data.size()) {      // code data only 
    if(phase) return 0.0;
    inc = 1;
  }   
  
  vector<t_gsatdata>::iterator it;
  int ii = 1;
  if(phase) ii = 2;
  double maxres = 0.0;
   
  for( it = data.begin(); it != data.end(); it++){
    if(it->gsys() != gs) {
      ii += inc;
      continue;
    }      
    
    if( maxres == 0.0 || abs(v(ii)) > maxres ){
      maxres = abs(v(ii));
      itDATA = it;
    }
    
    ii += inc;
  }
  
  return maxres;      
}
   
   
// Find maximal residual
double t_gsppmodel::_maxres(bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA, RESIDTYPE res_type, GSYS gs)
{      
   
  vector<t_gsatdata>::iterator it;
   
  double maxres = 0.0;
  
  for( it = data.begin(); it != data.end(); it++){
    if(it->gsys() != gs && gs != GNS) continue;
      
    vector<double> res;
    if(phase) res = it->residuals(res_type, TYPE_L);
    else      res = it->residuals(res_type, TYPE_C);
    
    for(auto itRES = res.begin(); itRES != res.end(); itRES++){
      if( maxres == 0.0 || fabs(*itRES) > maxres ){
        maxres = fabs(*itRES);
        itDATA = it;
      }
     }        
    
  }
  
  return maxres;      
}   

// check maximal residual   
bool t_gsppmodel::_check_outl(bool phase, double& maxres, vector<t_gsatdata>::iterator& itData, vector<t_gsatdata>& data)
{
  map<GSYS, double> map_res;
  if(phase) map_res = _maxres_L;
  else      map_res = _maxres_C;
   
  GSYS gs;
  if(itData != data.end()) gs = itData->gsys();
  else return false;
                      
  if( (maxres > map_res[gs]  && _resid_type == RES_ORIG) ||
      (maxres > _maxres_norm && _resid_type == RES_NORM) ){
    
    if(phase) _logOutl(true,  itData->sat(), data.size(),  maxres, itData->ele_deg(), itData->epoch(), _resid_type);
    else      _logOutl(true, itData->sat(), data.size(), maxres, itData->ele_deg(), itData->epoch(), _resid_type);
    
    return true;
  }
  return false;
}   
   
// check maximal residual   
bool t_gsppmodel::_check_outl(bool phase, double& maxresNORM, vector<t_gsatdata>::iterator& itDataNORM, 
                              double& maxresORIG, vector<t_gsatdata>::iterator& itDataORIG,
                              vector<t_gsatdata>::iterator& itDataErase, vector<t_gsatdata>& data)
{  
  map<GSYS, double> map_res;
  if(phase) map_res = _maxres_L;
  else      map_res = _maxres_C;
   
  GSYS gs;
  if(itDataORIG != data.end()) gs = itDataORIG->gsys();
  else return false;
   
  if(_resid_type == RES_ORIG){
    if(maxresORIG > map_res[gs]){
      itDataErase = itDataORIG;
      if(phase) _logOutl(true,  itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), _resid_type);
      else      _logOutl(false, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), _resid_type);
      return true;
    }
  }else if(_resid_type == RES_NORM){
    if(maxresNORM > _maxres_norm){
      itDataErase = itDataNORM;
      if(phase) _logOutl(true,  itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), _resid_type);
      else      _logOutl(false, itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), _resid_type);
      return true;
    }
   }else if(_resid_type == RES_ALL){
     if(maxresORIG > map_res[gs] || maxresNORM > _maxres_norm){
       if(itDataORIG == itDataNORM){
         itDataErase = itDataORIG;
         if(phase) _logOutl(true,  itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RES_ORIG);
         else      _logOutl(false, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RES_ORIG);
         return true;     
       }else{
         double ratioORIG = maxresORIG /  map_res[gs];
         double ratioNORM = maxresNORM / _maxres_norm;
         if(ratioNORM >= ratioORIG){
           itDataErase = itDataNORM;
           if(phase) _logOutl(true,  itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), RES_NORM);
           else      _logOutl(false, itDataNORM->sat(), data.size(), maxresNORM, itDataNORM->ele_deg(), itDataNORM->epoch(), RES_NORM);
           return true;
         }else{
           itDataErase = itDataORIG;
           if(phase) _logOutl(true,  itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RES_ORIG);
           else      _logOutl(false, itDataORIG->sat(), data.size(), maxresORIG, itDataORIG->ele_deg(), itDataORIG->epoch(), RES_ORIG);
           return true;
         }
       }
     }
   }
  
   return false;
}
   
// logging outlier   
void t_gsppmodel::_logOutl(bool phase, string prn,int data_size, double maxres, double ele, t_gtime epo, RESIDTYPE resid_type)
{
  string obsType = "";
  string resType = "";
  if(phase) obsType = "phase";
  else      obsType = "range";
  
  if(resid_type == RES_NORM) resType = "Norm residual";
  if(resid_type == RES_ORIG) resType = "Orig residual";
  if(resid_type == RES_ALL)  resType = "All residual";
   
  ostringstream os;
  os << _site << " outlier (" << resType << ": "  << obsType << ") " << prn
     << " size:" << fixed << setw(2) << data_size
     << " v: "   << fixed << setw(16) << right << setprecision(3) << maxres
     << " ele: " << fixed << setw(6) << setprecision(2) << ele
     << " " << epo.str_hms();
  if( _log )
    _log->comment(1, "gppp", os.str());
      
}

// devide multi-freq residuals vector into single-freq
vector<ColumnVector> t_gsppmodel::_devide_res(const ColumnVector& v_orig)
{
  vector<ColumnVector> vec;
   
  unsigned int k = 1;
  if(_observ == RAW_DOUBLE) k = 2;
  
  if(k == 1){
    vec.push_back(v_orig);
    return vec;
  }
  
  ColumnVector v_L1(v_orig.Nrows() / k);
  ColumnVector v_L2(v_orig.Nrows() / k);
  
  int i = 1;
  int j = 1;
  while(i <= v_orig.Nrows() - 1){
    v_L1(j) = v_orig(i  );
    v_L2(j) = v_orig(i+1);
    j+=1;
    i+=k;
  }
  
  vec.push_back(v_L1);
  vec.push_back(v_L2);
  
  return vec;
}
   
   
} // namespace
