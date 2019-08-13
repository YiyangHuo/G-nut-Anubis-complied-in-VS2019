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

#include "gproc/gsppflt.h"
#include "gutils/gmatrixconv.h"
#include "gutils/gtimesync.h"
#include "gmodels/gmltpth.h"
#include "gprod/gprodclk.h"

namespace gnut {

// constructor
// ----------
t_gsppflt::t_gsppflt(string mark, t_gsetbase* set)
: t_gspp(mark, set),
   _minsat(static_cast<size_t>(SPP_MINSAT)),
   _filter(0),
   _smooth(false),
   _n_NPD_flt(0),
   _n_ALL_flt(0),
   _n_NPD_smt(0),
   _n_ALL_smt(0)
{
  gtrace("t_gsppflt::t_gsppflt");

  t_gspp::_get_settings();
  
  _smooth =  dynamic_cast<t_gsetflt*>(_set)->smooth();
   
  string fltModStr(dynamic_cast<t_gsetflt*>(_set)->method_flt() );
  if      (fltModStr.compare("kalman") == 0) _filter = new t_kalman();
  else if (fltModStr.compare("srcf")   == 0) _filter = new t_SRF();
  else
    if( _log ) _log->comment(0, "gsppflt", _site + " not correct filter setting - check XML-config");
    else                           cerr << _site + " not correct filter setting - check XML-config\n";   

  _trpStoModel = new t_randomwalk();
  _trpStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_ztd());
  
  _ionStoModel = new t_whitenoise(dynamic_cast<t_gsetflt*>(_set)->noise_vion());
  //_ionStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_vion()); 

  _gloStoModel = new t_randomwalk();
  _gloStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_glo());
   
  _galStoModel = new t_randomwalk();
  _galStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_gal());
   
  _bdsStoModel = new t_randomwalk();
  _bdsStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_bds());

  _qzsStoModel = new t_randomwalk();
  _qzsStoModel->setq(dynamic_cast<t_gsetflt*>(_set)->rndwk_qzs());
   
  _clkStoModel = new t_whitenoise(dynamic_cast<t_gsetflt*>(_set)->noise_clk() );
  _crdStoModel = new t_whitenoise(dynamic_cast<t_gsetflt*>(_set)->noise_crd() );

  _vBanc.ReSize(4); _vBanc = 0.0;   
   
  int ipar = 0;   
   
  // Add coordinates parameters
  if( _crd_est != FIX ){
     _param.addParam( t_gpar(_site, t_gpar::CRD_X, ++ipar, "") );
     _param.addParam( t_gpar(_site, t_gpar::CRD_Y, ++ipar, "") );
     _param.addParam( t_gpar(_site, t_gpar::CRD_Z, ++ipar, "") );      
    }         

  // Add receiver clock parameter
  _param.addParam( t_gpar(_site, t_gpar::CLK, ++ipar, "") );
   
  // Add tropospheric wet delay parameter
  if( _tropo_est ){
     t_gpar par_trp(_site, t_gpar::TRP, ++ipar, "");
     par_trp.setMF( _ztd_mf );
     _param.addParam(par_trp);
  }   
   
  // Filling init parameter covariance matrix
  _Qx.ReSize( _param.parNumber() ); _Qx = 0.0;     
  double crdInit = _sig_init_crd;
  double ztdInit = _sig_init_ztd; 
  for (unsigned int i=1; i<=_param.parNumber(); i++) {
    if (_param[i-1].parType == t_gpar::CRD_X )   
    _Qx(i,i) = crdInit * crdInit;
    else if (_param[i-1].parType == t_gpar::CRD_Y)
    _Qx(i,i) = crdInit * crdInit;
    else if (_param[i-1].parType == t_gpar::CRD_Z)
    _Qx(i,i) = crdInit * crdInit;
    else if (_param[i-1].parType == t_gpar::CLK)
    _Qx(i,i) = _clkStoModel->getQ() * _clkStoModel->getQ();
    else if (_param[i-1].parType == t_gpar::TRP)
    _Qx(i,i) = ztdInit * ztdInit;    
  }
   
  _resid_type = dynamic_cast<t_gsetproc*>(_set)->residuals();
  _auto_band  = dynamic_cast<t_gsetproc*>(_set)->auto_band();
  
  _success = true;
  _ifb3_init = false;
  _ifb4_init = false;
  _ifb5_init = false;  

} // end constructor


// destructor
// ----------
t_gsppflt::~t_gsppflt()
{ 
  gtrace("t_gsppflt::~t_gsppflt");
   
  if( _trpStoModel ) delete _trpStoModel;   
  if( _ionStoModel ) delete _ionStoModel;   
  if( _filter      ) delete _filter;   
  if( _gloStoModel ) delete _gloStoModel;
  if( _galStoModel ) delete _galStoModel;
  if( _bdsStoModel ) delete _bdsStoModel;
  if( _qzsStoModel ) delete _qzsStoModel;   
  if( _clkStoModel ) delete _clkStoModel;
  if( _crdStoModel ) delete _crdStoModel;
}

// Process observation batch
// ----------
int t_gsppflt::processBatch(const t_gtime& beg, const t_gtime& end)
{

  gtrace("t_gsppflt::processBatch");

  if( beg > end ){ cerr << _site << " - processing not started [beg > end]\n"; return -1; }

#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif 
  _gmutex.lock();

  _dop.set_data(_gnav, _gobs, _site);   
   
  t_gtime begT( beg );
  t_gtime endT( end );

  _param.setSite(_site);
  t_gtriple xyz, neu, ell;

  double subint = 0.1;
  if(_scale    > 0) subint = 1.0 / _scale;
  if(_sampling > 1) subint = pow(10,floor(log10(_sampling)));

  int sign  = int(fabs(_sampling)/_sampling);  
     
   t_gtime now( begT );

#ifdef DEBUG
   if( _log ) _log->comment(1,"t_gsppflt", begT.str_ymdhms( "Filtering: " ) + endT.str_ymdhms(" "));
#endif

   t_gtime timeT(now); 
     
//   cout << _site << ": Start filtering: " << _site << " " << now.str_ymdhms() << endl;
   while(now < end || now == end)
   {
     // synchronization of now != end
      if (now != end){
         if (!time_sync(now, _sampling, _scale, _log) ){ // now.sod()%(int)_sampling != 0 ){ 
            now.add_dsec(sign * subint / 100); // add_dsec used for synchronization!
            continue;
         }
         if( _sampling > 1 ) now.reset_dsec();
      }

     // clean/collect/filter epoch data
     _data.erase( _data.begin(), _data.end() );
     _data = _gobs->obs( _site, now );
      
     if( _data.size() == 0){
       if( _log ) {_log->comment(1,"gsppflt", _site + now.str_ymdhms(" no observation found at epoch: "));}
      if (_sampling > 1) now.add_secs(int(sign*_sampling));  // =<1Hz data
      else               now.add_dsec(sign*_sampling);       //  >1Hz data
       continue;
     }

      t_gtime obsEpo = _data.begin()->epoch();       
     _timeUpdate(obsEpo); 

     int irc_epo = _processEpoch(obsEpo);

     if( irc_epo < 0 ) {
       _success = false;
       if( _log ) {_log->comment(1,"gsppflt", _site + now.str_ymdhms(" epoch ") + " was not calculated" );}
       if (_sampling > 1) now.add_secs(int(sign*_sampling));  // =<1Hz data
       else               now.add_dsec(sign*_sampling);       //  >1Hz data
       continue;
     }else _success = true;

     if( _log && _log->verb() >= 1 )
    _log->comment(1,"gsppflt", _site + now.str_ymdhms(" processing epoch: "));

     if( _param.getCrdParam(_site, xyz) > 0 ){
     }else{
       if( _log ) _log->comment(1,"gsppflt", _site + now.str_ymdhms(" No coordinates included in params: "));
     }

     double clk = 0.0;
     double clk_std = 0.0;
     int iclk = _param.getParam(_site, t_gpar::CLK, "");
     if(iclk >= 0) {
       clk_std = _Qx(_param[iclk].index, _param[iclk].index);
       clk = _param[iclk].value();
     }
      
     xyz2ell(xyz, ell, false);

     // CRD using eccentricities
     t_gtriple xyz_ecc = xyz  - _grec->eccxyz(now); // MARKER + ECC = ARP


     if (_allprod != 0){
       shared_ptr<t_gprodcrd> prd_crd = make_shared<t_gprodcrd>(now); // now is rounded value !!
       shared_ptr<t_gprodclk> prd_clk = make_shared<t_gprodclk>(now); // now is rounded value !!
       prd_crd->xyz(xyz_ecc);
       prd_clk->clk(clk, clk_std);
       prd_crd->nSat(_nSat);
       prd_crd->nSat_excl(_nSat_excl);

       // store DOP
       set<string> sat_list;
       for(auto it = _data.begin(); it != _data.end(); it++) sat_list.insert(it->sat());
       _dop.set_sats(sat_list);
       _dop.calculate(now, xyz);
       double gdop = _dop.gdop();
       double pdop = _dop.pdop();
       double hdop = _dop.hdop();
       double vdop = _dop.vdop();

       prd_crd->set_val("GDOP",gdop);
       prd_crd->set_val("PDOP",pdop);
       prd_crd->set_val("HDOP",hdop);
       prd_crd->set_val("VDOP",vdop);

       _allprod->add(prd_crd, _site);
       _allprod->add(prd_clk, _site);       
     }                  

     _initialized = true;

    if (_sampling > 1) now.add_secs(int(sign*_sampling));  // =<1Hz data
    else               now.add_dsec(sign*_sampling);       //  >1Hz data            
   }        
   
   _gmutex.unlock();

  return 1;
}

// Process one epoch
// -----------------------------
int t_gsppflt::_processEpoch(const t_gtime& runEpoch)
{
  gtrace("t_gsppflt::_processEpoch");

#ifdef DEBUG
  cout << "gsppflt: Epoch - full data size: " << _data.size() << " "
       << runEpoch.str_hms() << endl;
#endif

  if(_grec == nullptr) {
    if(_log) _log->comment(1, "gsppflt", "No receiver settings available!!!");
    return -1;
  }
  
  _epoch = runEpoch;
  t_gtriple crdapr = _grec->crd_arp(_epoch);

  if (double_eq(crdapr[0], 0.0) && double_eq(crdapr[1], 0.0) &&
      double_eq(crdapr[2], 0.0)) {
    _valid_crd_xml = false;
  } else {
    _valid_crd_xml = true;
  }

  if (!_valid_crd_xml) _sig_init_crd = 100.0;

  Matrix A;
  DiagonalMatrix P;
  ColumnVector l, dx, dx_el;
  ColumnVector v, v_orig, v_norm;
  SymmetricMatrix Qsav, QsavBP;
  vector<t_gsatdata>::iterator it;
  vector<GOBSTYPE> obstypes;
  t_gallpar XsavBP;
   
  _cntrep = 0;    // number of iterations caused by outliers

  do{
    if( _prepareData() < 0 )  {
      if(_log) _log->comment(2, "gsppflt", "Preparing data failed!");
      return -1;
    }         

#ifdef DEBUG
    cout << "gsppflt: Epoch - proc data size: " << _data.size() << " " << runEpoch.str_hms() << endl;cout.flush();
#endif               

    QsavBP = _Qx; XsavBP = _param;
    _predict();

    if( _data.size() < _minsat ){
      if(_log) _log->comment(2, "gsppflt", "Not enough visible satellites!");
      _restore(QsavBP, XsavBP);
      return -1;
    }

    t_gtriple groundEll, groundXYZ; 

    if ( _param.getCrdParam(_site, groundXYZ) > 0){
    }else if( _valid_crd_xml ) { groundXYZ = _grec->crd_arp(_epoch); }

    xyz2ell(groundXYZ.crd_array(), groundEll.crd_array(), false);

    // define a number of measurements 
    unsigned int nObs = _data.size();
    unsigned int mult = 1;
    if(_observ == RAW_DOUBLE) {mult = 2; nObs *= 2;}
    if(_observ == RAW_ALL)    {mult = 2; nObs *= 5;} // reservation for 5 freq - not used raws will be removed
    if(_phase) {mult *= 2; nObs *= 2;}  // code + phase
     
#ifdef DEBUG
  cout << "gsppflt: nObs: " << _data.size() << " "<< nObs << " " << runEpoch.str_hms() << endl;cout.flush();
#endif               

    unsigned int nPar = _param.parNumber();
    A.ReSize(nObs, nPar); A = 0.0;
    l.ReSize(nObs);       l = 0.0;
    P.ReSize(nObs);       P = 0.0;
    dx.ReSize(nPar);     dx = 0.0;
    unsigned int iobs = 1;

    _frqNum.clear();

    // Create matrices and vectors for estimation
    // ----------------------------------------------------------------------------------
    // loop over all measurement in current epoch

    for( it = _data.begin(); it != _data.end(); ){

      if( _addObsP(*it, iobs, groundEll, A, l, P) > 0 ) obstypes.push_back(TYPE_C);
      else{ it = _data.erase(it); continue; } // just return next iterator

      if( _addObsL(*it, iobs, groundEll, A, l, P) > 0 ) obstypes.push_back(TYPE_L);
      else{ it = _data.erase(it); continue; } // just return next iterator

      ++it; // increase iterator if not erased
    }

//    cout << _epoch.str_hms() << endl;
//    cout << "Predict cov: \n" << fixed << setprecision(0) << _Qx << endl;

    if( iobs <= _minsat*mult) {
      if(_log) _log->comment(2, "gsppflt", "Not enough processable observations!");
      _restore(QsavBP, XsavBP);     
      return -1;
    }

     // delete zero rows/cols
     iobs--;
     A = A.Rows(1,iobs);
     P = P.SymSubMatrix(1,iobs);
     l = l.Rows(1,iobs);
     
    
     Qsav =  _Qx;

     _filter->update(A, P, l, dx, _Qx);     

     // increasing variance after update in case of introducing new ambiguity
    for( size_t iPar=0; iPar<_param.parNumber(); iPar++ ){
      if( _param[iPar].parType == t_gpar::AMB_IF ){
        string sat = _param[iPar].prn;
        if ( _newAMB.find(sat) != _newAMB.end() && _cntrep == 1){
          if(_newAMB[sat] == 1) _Qx(iPar+1, iPar+1) += 10;
          if(_newAMB[sat] == 2 && _Qx(iPar+1, iPar+1) > 0.01) _Qx(iPar+1, iPar+1) += 1;
          _newAMB[sat]++;
        }          
      }     
    }               

     // post-fit residuals
     v_orig = l - A*dx;
     
     // normalized post-fit residuals
     Matrix Qv = A*_Qx*A.t() + P.i();
     v_norm.ReSize(v_orig.Nrows());
     for(int i = 1; i <= v_norm.Nrows(); i++){
            v_norm(i) = sqrt(1/Qv(i,i)) * v_orig(i);
     }          

     int freedom = A.Nrows() - A.Ncols();
     if (freedom < 1) {
      if (_log) _log->comment(1,"t_gsppflt", "No redundant observations!");
//    _restore(QsavBP, XsavBP);     
//    return -1;
      freedom = 1;
     }
    
     ColumnVector vtPv = v_orig.t() * P * v_orig;
     _sig_unit = vtPv(1) / freedom;        
    
    /* Matrix N = A.t()*P*A; */
    /* Matrix Ql = A*N.i()*A.t(); */
    /* cout << fixed << setprecision(0) << Ql << endl; */
    /* int ooo; cin >> ooo; */
    
#ifdef DEBUG
//    SymmetricMatrix Corr = cov2corr(_Qx);
    cout << "Used satellites: ";
    for (unsigned int i = 0; i < _data.size(); i++) cout << _data[i].sat() << " ";
    cout << endl;
    cout << "Used parameters: ";
    for (unsigned int i = 0; i < _param.parNumber(); i++) cout << _param[i].str_type() << "_" << _param[i].prn << " ";
    cout << endl;
    cout << fixed   << setprecision(3)
         << _epoch.str_ymdhms() << endl
//         << " m0: " << setw(7) << _sig_unit << endl
//         << " dx: " << setw(7) << dx
         << " l.Nrows(): " << l.Nrows() << endl
//         << setw(7) << l 
//         << " v_orig.Nrows(): " << v_orig.Nrows() 
//         << " v_orig: "  << setw(7) << v_orig
//         << " v_norm: "  << setw(7) << v_norm
         << " A(" << A.Nrows() << "," << A.Ncols() << "):" << setw(7) << endl << A
//         << " P: "  << setw(7) << P
         << " Q: "  << setw(7) << _Qx
//         << " C: "  << setw(7) << Corr    
//         << " Param " << setw(7) << _param
         << endl; cout.flush();
    
    // for(int i = 1; i <= _param.parNumber(); i++) cout << _param[_param.getParam(i)].str_type() << " "; cout << endl;
    // for(int i = 1; i <= _param.parNumber(); i++) cout << _param[_param.getParam(i)].str_type() << " " << Corr.row(i) << " "; cout << endl;

//int ooo; cin >> ooo;
#endif

    // save post-fit residuals
    _save_residuals(v_norm, _data, RES_NORM);
    _save_residuals(v_orig, _data, RES_ORIG);          

  }while( _gModel->outlierDetect(_data, _Qx, Qsav) != 0 );

  // test purpose - print L1, L2, L3 residuals
  /* for(auto it = _data.begin(); it != _data.end(); it++){ */
  /*   vector<double> vres = it->residuals(RES_ORIG, TYPE_L); */
  /*   cout << it->epoch().str_ymdhms() << " " << it->sat(); */
  /*   for(auto itRES = vres.begin(); itRES != vres.end(); itRES++){ */
  /*      cout << fixed << setprecision(3) << setw(10) << *itRES; */
  /*   } */
  /*   if(vres.size() < 3) cout << setw(10) << 0.0; */
  /*   cout << endl; */
  /* } */
   
  _nSat_excl = _nSat - _data.size();   // number of excluded satellites due to different reasons

  try{
     _n_ALL_flt++;     
     Cholesky(_Qx);
  }catch(NPDException) {
//     cout << _site << " " << _epoch.str_ymdhms() << " NPDException " << endl;
     _n_NPD_flt++;
  }
 // end of covariance matrix testing     

  if( _data.size() < _minsat ){
    if( _log ) _log->comment(2, "gsppflt", _site + _epoch.str_ymdhms(" epoch ")
                        + " skipped: "+int2str(_data.size())+" < _minsat)");
//    cout << "Epoch " << _epoch.str_hms() << " not calculated due to _minsat." << endl;     
    _restore(QsavBP, XsavBP);
    return -1;
  }

  // edit and save ZHD - for post-fit residual
  // (Must be before updating crd due to maintain consistency with a priory)
  t_gtriple Ell, XYZ;    
  if (_param.getCrdParam(_site, XYZ) > 0){
  }else if ( _valid_crd_xml ) { XYZ = _grec->crd_arp(_epoch); }
  xyz2ell( XYZ, Ell, false);   
  int itrp = _param.getParam(_site, t_gpar::TRP, "");
  if(itrp >= 0) {
    if(_gModel->tropoModel() != 0){
     _param[itrp].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch) );
    }
  }
   
  // Edit and save parematres
  for( unsigned int iPar = 0; iPar < _param.parNumber(); iPar++ ) {
    _param[iPar].value( _param[iPar].value() + dx(_param[iPar].index));
  }

  return 1;   
   
}

// weight coef for observations
// ---------------------------
double t_gsppflt::_weightObs(t_gsatdata& satdata, t_gobs& go)
{
   double weight_coef = 0;
   
   switch (_weight) {
    case DEF_OBSWEIGHT: cerr << "gsppflt: WeightObs (default) should not happened!\n"; break;
    case EQUAL:  weight_coef = 1; break;
    case SINEL:  //weight_coef = sin(satdata.ele()); break;
      if (satdata.ele()*180.0 / G_PI > 40) weight_coef = 1.0;
      else weight_coef = sin(satdata.ele())*sin(satdata.ele());
    case SINEL2: weight_coef = pow( sin(satdata.ele()), 2 ); break;
    case SINEL4: weight_coef = pow( sin(satdata.ele()), 4 ); break;
    case CODPHA: 
      if( go.is_code()  ) weight_coef = pow( sin(satdata.ele()), 4 );
      if( go.is_phase() ) weight_coef = pow( sin(satdata.ele()), 2 );      
      break;            
    case MLTPTH: weight_coef = sin(satdata.ele());   
      if( go.is_code() ){
    t_gtime tt = satdata.epoch();
    t_gmltpth gmltpth(_gobs, _set, _site);
    t_gmltpth::t_map_mpobs map;
    vector<GOBS> obs = satdata.obs();
    pair<double, int> p;
    for (vector<GOBS>::iterator it = obs.begin(); it != obs.end(); it++) {
       if( gobs_phase(*it) ){
          map[*it][FIRST_TIME][""] = p;
       }  
    }   
    GOBS g = go.gobs();
    map[g][FIRST_TIME][""] = p;
    
    gmltpth.multipath(satdata.sat(), satdata.gsys(), satdata.epoch(), 5, map);
    double mltpth = map[g][satdata.epoch()][satdata.sat()].first;
    double mltpth_weight = 0;
    if(mltpth > 0) mltpth_weight = 1/(mltpth*1e-2*mltpth*1e-2);
    //      if (satdata.sat() == "G30") cout << satdata.epoch().str_ymdhms() << " " << fixed << setprecision(3) << weight_coef << " " << mltpth << " " << weight_coef + tanh(mltpth_weight/100) << endl;
    weight_coef += tanh(mltpth_weight/100);   
      }   
   }
   
   return weight_coef;
}

// add pseudo observations as a constrains
// ---------------------------
int t_gsppflt::_addPseudoZTD(unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P)
{
   gtrace("t_gsppflt::_addPseudoZTD");

   // Design Matrix
   if( _epoch > _ztd_begStat && _epoch < _ztd_endStat )
     {
#ifdef DEBUG
   cout << _epoch.str_ymdhms(" epoch: ") << endl;
   cout << _ztd_begStat.str_ymdhms("_ztd_begStat: ") << endl;
   cout << _ztd_endStat.str_ymdhms("_ztd_endStat: ") << endl;
#endif
   
   Matrix_addRC(A, A.Nrows()+1, 0);
   int i = _param.getParam(_site, t_gpar::TRP, "");
   A(A.Nrows(), i+1) = 1;
   
   // tropo (wet part) is constrained to this value 
   double TRP_fix = 0;
   if( _gModel->tropoModel() != 0 ) TRP_fix = _aprox_ztd_xml - _gModel->tropoModel()->getZHD(ell, _epoch);       
       
   // Reduced measurement
   Vector_add(l, l.Nrows()+1);
   l(l.Nrows()) = _param[i].value() - TRP_fix;
   
   // weight matrix
   Matrix_addRC(P, P.Nrows()+1);
   P(P.Nrows(), P.Ncols()) = 99999;

   iobs++;
    }
   
   
   return 1;
}
           
// add data to A, l, P - code
// ---------------------------
int t_gsppflt::_addObsP(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P)
{
   gtrace("t_gsppflt::_addObsP");
   
   t_gsys gsys(satdata.gsys());
   GSYS gs = gsys.gsys();
  
  GOBSBAND b1, b2, b3, b4, b5;
  
  if(_auto_band){    // automatic dual band selection -> for Anubis purpose
    set<GOBSBAND> bands = satdata.band_avail(_phase);
    auto itBAND = bands.begin();
    if(bands.size() < 2) return -1;
    b1 = *itBAND;
    itBAND++;
    b2 = *itBAND;
  }else{            // select fix defined band according the table
    b1 = t_gsys::band_priority(gs, FREQ_1);
    b2 = t_gsys::band_priority(gs, FREQ_2);   
    b3 = t_gsys::band_priority(gs, FREQ_3);
    b4 = t_gsys::band_priority(gs, FREQ_4);
    b5 = t_gsys::band_priority(gs, FREQ_5);
  }
  
   t_gobs gobs1, gobs2, gobs3, gobs4, gobs5;                       
   
   _getgobs(satdata.sat(), TYPE_C, b1, gobs1);
   _getgobs(satdata.sat(), TYPE_C, b2, gobs2);
   _getgobs(satdata.sat(), TYPE_C, b3, gobs3);
   _getgobs(satdata.sat(), TYPE_C, b4, gobs4);
   _getgobs(satdata.sat(), TYPE_C, b5, gobs5);      
//cout << satdata.epoch().str_hms() << " " << satdata.sat() << " " << gobs2str(satdata.id_range(gobs1.band())) << " " 
//                                                                 << gobs2str(satdata.id_range(gobs2.band())) << " "
//                                                                 << gobs2str(satdata.id_range(gobs3.band())) << " "
//                                                                 << gobs2str(satdata.id_range(gobs4.band())) << " "  
//                                                                 << gobs2str(satdata.id_range(gobs5.band())) << endl;

   double PIF, P1, P2, P3, P4, P5;
   PIF = P1 = P2 = P3 = P4 = P5 = 0.0;

   PIF = satdata.P3(gobs1, gobs2);
   P1  = satdata.obs_C(gobs1);
   P2  = satdata.obs_C(gobs2);    
   if(b3 != BAND) P3 = satdata.obs_C(gobs3);
   if(b4 != BAND) P4 = satdata.obs_C(gobs4);
   if(b5 != BAND) P5 = satdata.obs_C(gobs5);      

   if(_observ == IONO_FREE  &&  double_eq(PIF, 0.0)) return -1;
   if(_observ == RAW_SINGLE &&  double_eq(P1, 0.0)) return -1;
   if((_observ == RAW_DOUBLE || _observ == RAW_ALL) && (double_eq(P1, 0.0) || double_eq(P2, 0.0))) return -1;   

   double weight_coef_b1 = _weightObs(satdata, gobs1);
   double weight_coef_b2 = _weightObs(satdata, gobs2);
   double weight_coef_b3 = _weightObs(satdata, gobs3);
   double weight_coef_b4 = _weightObs(satdata, gobs4);
   double weight_coef_b5 = _weightObs(satdata, gobs5);   
   double weight_coef = 0;
   
   if(_observ == IONO_FREE){
     double koef1 = 0.0;
     double koef2 = 0.0;
     satdata.coef_ionofree(b1, koef1, b2, koef2);   
     weight_coef = (weight_coef_b1 + weight_coef_b2)/2;
   }else weight_coef = weight_coef_b1;

   int addIobs = 1;
   if(P3 > 0) addIobs++;
   if(P4 > 0) addIobs++;
   if(P5 > 0) addIobs++;     

   int usedFrq = addIobs + 1;
   if(_observ == IONO_FREE || _observ == RAW_SINGLE) _frqNum[satdata.sat()] = 1;
   if(_observ == RAW_DOUBLE) _frqNum[satdata.sat()] = 2;
   if(_observ == RAW_ALL)    _frqNum[satdata.sat()] = usedFrq;
   
   double sigCode = 0.0;
   switch( gs ){
   case GPS : sigCode = _sigCodeGPS; break;
   case GLO : sigCode = _sigCodeGLO; break;
   case GAL : sigCode = _sigCodeGAL; break;
   case BDS : sigCode = _sigCodeBDS; break;
   case QZS : sigCode = _sigCodeQZS; break;
   default  : sigCode = 0.0;
   }    

  // Create weight matrix
  P(iobs,  iobs  ) = weight_coef * (1/(sigCode * sigCode));
  if(_observ == RAW_DOUBLE ||
     _observ == RAW_ALL)      P(iobs+1,iobs+1) = weight_coef_b2 * (1/(sigCode * sigCode));
  if(_observ == RAW_ALL)
  for(int i = 2; i <= addIobs; i++) {
    if(i == 2) P(iobs+i,iobs+i) = weight_coef_b3 * (1/(sigCode * sigCode));
    if(i == 3) P(iobs+i,iobs+i) = weight_coef_b4 * (1/(sigCode * sigCode));
    if(i == 4) P(iobs+i,iobs+i) = weight_coef_b5 * (1/(sigCode * sigCode));      
  }
  
  bool com = true;
  if( _gnav ) com = _gnav->com();
  
   //Create reduced measurements (prefit residuals)
   double modObsP = 0.0;
   if(_observ == IONO_FREE ) {
     modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs1, com); if (modObsP < 0) return -1;
     _applyDCB(satdata, PIF, &gobs1, &gobs2);
     l(iobs  ) = PIF - modObsP;

   }else if(_observ == RAW_SINGLE) {
     modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs1, com); if (modObsP < 0) return -1;
     _applyDCB(satdata, P1, &gobs1); 
     l(iobs  ) = P1 - modObsP;
   }else if(_observ == RAW_DOUBLE){
     modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs1, com); if (modObsP < 0) return -1;
     _applyDCB(satdata, P1, &gobs1);
     l(iobs  ) = P1 - modObsP;
     modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs2, com); if (modObsP < 0) return -1;    
     _applyDCB(satdata, P2, &gobs2);
     l(iobs+1) = P2 - modObsP;
   }else if(_observ == RAW_ALL){
     for(int i = 0; i <= addIobs; i++) {
       double Pi = 0.0;
       if(i == 0) {modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs1, com); _applyDCB(satdata, P1, &gobs1); Pi = P1;}
       if(i == 1) {modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs2, com); _applyDCB(satdata, P2, &gobs2); Pi = P2;}
       if(i == 2) {modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs3, com);                                 Pi = P3;}
       if(i == 3) {modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs4, com);                                 Pi = P4;}
       if(i == 4) {modObsP = _gModel->cmpObs(_epoch, _param, satdata, gobs5, com);                                 Pi = P5;}
       if(double_eq(modObsP, 0.0) || double_eq(Pi, 0.0)) continue;
       l(iobs+i) = Pi - modObsP;
     }
     
   }

   // Create first design matrix
   if ( gs == GLO || 
      gs == GPS ||
      gs == GAL ||
      gs == BDS ||
      gs == QZS ){
      
      for( unsigned int ipar=1; ipar <= _param.parNumber(); ipar++ ){
                                    A(iobs,  ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs1);
       if(_observ == RAW_DOUBLE || 
          _observ == RAW_ALL)       A(iobs+1,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs2);
       if(_observ == RAW_ALL){
         for(int i = 2; i <= addIobs; i++) {         
           if(i == 2) A(iobs+i,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs3);
           if(i == 3) A(iobs+i,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs4);
           if(i == 4) A(iobs+i,ipar) = _param[ipar-1].partial( satdata, _epoch, ell, gobs5);
         }
       }
      }    
   } else return -1;
   
#ifdef DEBUG       
      cout << satdata.sat()     << fixed << setprecision(10)
     << "  PIF: (" << BAND_1 << "," << BAND_2 << ") :"  << setw(12) << satdata.P3(b1,b2)
           << "  cmpObs: "  << setw(12) << modObsP
           << endl;
#endif   
   
   if(_observ == RAW_DOUBLE) iobs += 2;
   else if(_observ == RAW_ALL) iobs = iobs + addIobs + 1;   
   else iobs++;
   
   return 1;
}

// add data to A, l, P - carrier phase
// ---------------------------
int t_gsppflt::_addObsL(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P)
{
   gtrace("t_gsppflt::_addObsL");   
   return 1;
}


// Prepare data: 
// filter, bancroft, members in gsatdata
// ----------
int t_gsppflt::_prepareData()
{   
  gtrace("t_gsppflt::_prepareData");
       
  vector<t_gsatdata>::iterator it = _data.begin();
  while (it != _data.end()) {       
    switch (_gnss) {
      case GNS: break;
      case GPS: if (it->gsys() != GPS) { it = _data.erase(it); continue; } break;
      case GLO: if (it->gsys() != GLO) { it = _data.erase(it); continue; } break;
      case GAL: if (it->gsys() != GAL) { it = _data.erase(it); continue; } break;
      case BDS: if (it->gsys() != BDS) { it = _data.erase(it); continue; } break;   
      case QZS: if (it->gsys() != QZS) { it = _data.erase(it); continue; } break;
      case IRN: if (it->gsys() != IRN) { it = _data.erase(it); continue; } break;   
      case SBS: if (it->gsys() != SBS) { it = _data.erase(it); continue; } break;
    }
    it++;
   }
  
  Matrix BB;

  BB.ReSize(_data.size(), 4);
  BB = 0.0;   
  int iobs = 0;
  vector<t_gsatdata>::iterator iter = _data.begin();

  _nSat = _data.size();  // configured GNSS is considered (see erase in above swich)

  while( iter != _data.end() ){
          
    GSYS gs = iter->gsys();

    GOBSBAND b1, b2;

    if(_auto_band){
      set<GOBSBAND> bands;
      bands = iter->band_avail(_phase);
      auto itBAND = bands.begin();
      if(bands.size() < 2) {iter++; continue;}
      b1 = *itBAND;
      itBAND++;
      b2 = *itBAND;
    }else{
      b1 = t_gsys::band_priority(gs, FREQ_1);
      b2 = t_gsys::band_priority(gs, FREQ_2);
    }

    iter->glog(_log);

    // check data availability
    double P3 = iter->P3(b1, b2);
    double L3 = iter->L3(b1, b2);          

#ifdef DEBUG
    cerr << "iter->epoch/_epoch: " << iter->epoch().str_hms() << _epoch.str_hms(" / ")
         << " " << t_gsys::gsys2str(gs) << " " << iter->sat() << " " << iter->channel() 
    << " " << fixed << setprecision(3) << P3 << " " << L3 << endl;
#endif

    if( double_eq(L3, 0.0) && _phase) {
      if( _log ) {
        ostringstream str;
        str << "prepareData: erasing data due to no phase double bands observation, " <<
        "epo: " << _epoch.str_hms() << ", "        <<
        "prn: " << iter->sat() << " (" << iter->channel() << ")" ;
        _log->comment(2, "gsppflt", str.str());
      }
      iter = _data.erase(iter); // !!!! OK u vektoru funguje bez ++ (erase zvysuje pointer !)      
      continue;
    }
            
    if( double_eq(P3, 0.0) ) {
      if( _log ) {
        ostringstream str; 
        str << "prepareData: erasing data due to no code double bands observation, " <<
        "epo: " << _epoch.str_hms() << ", "        <<
        "prn: " << iter->sat();
        _log->comment(2, "gsppflt", str.str());
      }
      iter = _data.erase(iter); // !!!! OK u vektoru funguje bez ++ (erase zvysuje pointer !)      
      continue;       
    }

    if( _satPos( _epoch, *iter ) < 0 ){
      if( _log ) {
        ostringstream str; 
        str << "prepareData: erasing data since _satPos failed, "
            << "epo: " << _epoch.str_hms() << ", "        <<
            "prn: " << iter->sat();
        _log->comment(2, "gsppflt", str.str());
      }       
#ifdef DEBUG
    cerr << "prepareData: sat pos not calculated "
    << _epoch.str_ymdhms(_site+" "+iter->sat()+" "+int2str(_data.size())+" ")
    << endl; cerr.flush();
#endif       
      iter = _data.erase(iter); // !!!! OK u vektoru funguje bez ++ (erase zvysuje pointer !)

    }else{
      iobs++;
      BB(iobs,1) = iter->satcrd().crd(0);
      BB(iobs,2) = iter->satcrd().crd(1);
      BB(iobs,3) = iter->satcrd().crd(2);
      BB(iobs,4) = iter->P3(b1,b2) + iter->clk();
      iter++;
      continue;
    }
  }

  if( _data.size() < _minsat ) {
    if( _log )
    _log->comment(2, "gsppflt", _epoch.str_ymdhms( _site + " epoch ")
                  + " skipped (Bancroft not calculated: "+int2str(_data.size())+" < _minsat)");
//      cout << t_gsys::gsys2str(_gnss) << " " << _epoch.str_ymdhms( _site + " epoch ") << " " << int2str(_data.size()) <<  " skipped (Bancroft not calculated: _data.size < _minsat)" << endl;
    return -1;
  }

  BB = BB.Rows(1,iobs);    // delete zero rows
  
  if( BB.Nrows() < static_cast<int>(_minsat) ){
    if( _log && _log->verb() >= 2 )
    _log->comment(2, "gsppflt", _epoch.str_ymdhms(_site + " epoch ")
                  + " skipped (Bancroft not calculated: BB.Nrows < _minsat)");
//      cout << _epoch.str_ymdhms( _site + " epoch ") << int2str(_data.size()) <<  " skipped (Bancroft not calculated: _data.size < _minsat)" << endl;
    return -1;
  }

  _vBanc = 0.0;
  gbancroft(BB, _vBanc);   
   
  if( !_initialized && (_vBanc.Rows(1,3) - _grec->crd_arp(_epoch).crd_cvect()).NormFrobenius() > 1000) _valid_crd_xml = false;
  
  if ( _valid_crd_xml && !_initialized ) {
    _vBanc(1) = _grec->crd_arp(_epoch)[0];
    _vBanc(2) = _grec->crd_arp(_epoch)[1];
    _vBanc(3) = _grec->crd_arp(_epoch)[2];
  }
  
  t_gtriple test_xyz(_vBanc);
  t_gtriple test_ell;
  xyz2ell(test_xyz, test_ell, true);
  double radius = test_xyz.norm(); 
  if(radius < B_WGS - 500) {
    string warning = "WARNING: Unexpected site (" + _site + ") coordinates from Bancroft. Orbits/clocks or code observations should be checked.";
    if(_log)  _log->comment(1, "gsppflt", warning + _epoch.str_ymdhms(" Epoch "));
    if(_glog) _glog->comment(2, "gsppflt", warning + _epoch.str_ymdhms(" Epoch "));      
    if(!_phase) {return -1;}      
  }
   
   
#ifdef DEBUG
   cout << _epoch.str_ymdhms() << endl; cout.flush();
//   cout << "\nBB\n" << fixed << setprecision(10) << BB << endl;
   cout << "vBanc\n"  << fixed << setprecision(10)<< _vBanc << endl; cout.flush();
//   cout << "Clk [mili sec]" << fixed << setprecision(10) << _vBanc(4)*1000/CLIGHT << endl;
//   cout << "_Qx\n"  << fixed << setprecision(3)<< _Qx << endl;    cout.flush();
//   int ooo; cin >> ooo;
#endif

  //Compute sat elevation and rho
  iter = _data.begin();
  while( iter != _data.end() ){
    t_gtriple xyz_r, xyz_s, neu_s;
    xyz_s = iter->satcrd();
    if ( _crd_est == FIX && _valid_crd_xml ){
      xyz_r = _grec->crd_arp(_epoch);
    }else if ( _crd_est == EST ){
      if( !_initialized ){
        xyz_r.set(_vBanc.Rows(1,3));
      }else{
        if(_pos_kin){ xyz_r.set(_vBanc.Rows(1,3));
        }else _param.getCrdParam(_site, xyz_r);
      }       
    }else {
      cerr << "Undefined CRD constraining" << endl;
      return -1;
    }
    
#ifdef DEBUG
    cout << "Coordinates: " << iter->sat() << " " << iter->epoch().str_hms()
         << fixed << setprecision(10) << endl
         << " X_R " << setw(12) << xyz_r[0]
         << " Y_R " << setw(12) << xyz_r[1]
         << " Z_R " << setw(12) << xyz_r[2] << endl
         << " X_S " << setw(12) << xyz_s[0]
         << " Y_S " << setw(12) << xyz_s[1]
         << " Z_S " << setw(12) << xyz_s[2] << endl; cout.flush();
     int ooo; cin >> ooo;
#endif
    
    t_gtriple xyz_rho = xyz_s - xyz_r;
    t_gtriple ell_r;
    xyz2ell(xyz_r, ell_r, false);
    
    xyz2neu(ell_r, xyz_rho, neu_s);

    // Earth rotation correction
    t_gtriple xRec;
    double rho0 = sqrt(pow(xyz_r[0]-xyz_s[0],2) + pow(xyz_r[1]-xyz_s[1],2) + pow(xyz_r[2]-xyz_s[2],2));
    double dPhi = OMEGA * rho0 / CLIGHT;
    xRec[0] = xyz_r[0] * cos(dPhi) - xyz_r[1] * sin(dPhi);
    xRec[1] = xyz_r[1] * cos(dPhi) + xyz_r[0] * sin(dPhi);
    xRec[2] = xyz_r[2];

    // Apply tides 
    _apply_tides(_epoch, xRec);

    double tmp = (iter->satcrd().crd_cvect() - xRec.crd_cvect()).norm_Frobenius();
    
    iter->addrho( tmp );    
    //cout << iter->sat() << " " << iter->epoch().str_ymdhms() << fixed << setprecision(3) << " " << neu_s[0] << " " << neu_s[1] << " " << neu_s[2] << endl;
    double NE2 = neu_s[0]*neu_s[0] + neu_s[1]*neu_s[1];
    double ele = acos(sqrt(NE2)/iter->rho());
    if( neu_s[2]<0.0 ) {
      ele *= - 1.0;
    }
    
    if( sqrt(NE2)/iter->rho() > 1.0 )
    iter->addele(0.0);
    else iter->addele(ele);
    
    double azi = atan2(neu_s[1], neu_s[0]);
    if (azi < 0) azi += 2*G_PI;
    iter->addazi(azi);
    if(!_use_ecl) iter->addecl(_lastEcl);
    
    // check elevation cut-off 
    if ( iter->ele_deg() < _minElev ){
      ostringstream os;        
      os << "Erasing " << iter->sat() << " data due to low elevation angle (ele = " << fixed << setprecision(1) << iter->ele_deg()
         << ") " << iter->epoch().str_ymdhms();
      if( _log ) _log->comment(2, "gsppflt", os.str()); 
      iter = _data.erase(iter); // !!!! zveda iterator !!!!
      if( _data.size() < _minsat ){
        if( _log ) _log->comment(2, "", _site + _epoch.str_ymdhms(" epoch ") + int2str(_data.size()) + " skipped (data.size < _minsat)");             
        return -1;
      }
      continue;
    }
    
    // Printing beta and orbit angles for deep verbosity
    ostringstream os;
    os << iter->sat() << " " << iter->epoch().str_ymdhms() << " " << fixed << setprecision(1) << iter->beta() << " " << iter->orb_angle() << endl;
    if(_log) _log->comment(4, "gsppflt", os.str()); 
    
    // check satellite eclipsing
    if ( iter->ecl() ){
      ostringstream os;
      os << "Erasing " << iter->sat() << " data due to satellite eclipsing (beta = " << fixed << setprecision(1) << iter->beta()
         << " ,orbit angle = " << iter->orb_angle() << ") " << iter->epoch().str_ymdhms();
      if( _log ) _log->comment(2, "gsppflt", os.str());
      iter = _data.erase(iter); // !!!! zveda iterator !!!!
      if( _data.size() < _minsat ){
        if( _log ) _log->comment(2, "", _site + _epoch.str_ymdhms(" epoch ") + int2str(_data.size()) + " skipped (data.size < _minsat)");             
        return -1;
      }    
    } else iter++;       
    
  }
    return 1;
}


// Tides are not applied in SPP
// ---------------------------------
int t_gsppflt::_apply_tides(t_gtime& _epoch, t_gtriple& xRec)
{
  return 1;
}


// Predict state vector and covariance matrix
// =======================================================
void t_gsppflt::_predict()
{
  gtrace("t_gsppflt::_predict");

  _cntrep++;
  
  _Noise.ReSize(_Qx.Nrows()); _Noise = 0;
   
  // Predict coordinates, clock and troposphere
  // state vector is the same, covariance matrix is predicted with noise   
  int i;
   
  double crdInit       = _sig_init_crd;
  double ztdInit       = _sig_init_ztd;
   
  i = _param.getParam(_site, t_gpar::CRD_X, "");
  if ( i >= 0 )  {
    if ( !_initialized || _Qx(i + 1, i + 1) == 0.0) {
      _param[i].value( _vBanc(1) );
      _Qx(i+1,i+1) = crdInit * crdInit;
    }else{
      if( _epoch > _crd_begStat && _epoch < _crd_endStat){
        if(_smooth) _Noise(i+1,i+1) = 0.0;
      }else{ 
        if(_pos_kin) _param[i].value( _vBanc(1) );
        if (_cntrep == 1 && _success) _Qx(i+1,i+1) += _crdStoModel->getQ() * _crdStoModel->getQ();
        if(_smooth) _Noise(i+1,i+1) = _crdStoModel->getQ() * _crdStoModel->getQ();  
      }
    }
  }

  i = _param.getParam(_site, t_gpar::CRD_Y, "");
  if ( i >= 0 )  {
    if ( !_initialized || _Qx(i + 1, i + 1) == 0.0) {
      _param[i].value( _vBanc(2) );
      _Qx(i+1,i+1) = crdInit * crdInit;
    }else{
      if( _epoch > _crd_begStat && _epoch < _crd_endStat){
        if(_smooth) _Noise(i+1,i+1) = 0.0;
      }else{        
        if(_pos_kin) _param[i].value( _vBanc(2) );     
        if (_cntrep == 1 && _success) _Qx(i+1,i+1) += _crdStoModel->getQ() * _crdStoModel->getQ();
        if(_smooth) _Noise(i+1,i+1) = _crdStoModel->getQ() * _crdStoModel->getQ();  
      }
    }      
  }   
   
  i = _param.getParam(_site, t_gpar::CRD_Z, "");
  if ( i >= 0 )  {
    if ( !_initialized || _Qx(i + 1, i + 1)==0.0) {
      _param[i].value( _vBanc(3) );
      _Qx(i+1,i+1) = crdInit * crdInit;
    }else{
      if( _epoch > _crd_begStat && _epoch < _crd_endStat){
        if(_smooth) _Noise(i+1,i+1) = 0.0;
      }else{
        if(_pos_kin) _param[i].value( _vBanc(3) );     
        if (_cntrep == 1 && _success) _Qx(i+1,i+1) += _crdStoModel->getQ() * _crdStoModel->getQ();
        if(_smooth) _Noise(i+1,i+1) = _crdStoModel->getQ() * _crdStoModel->getQ();  
      }     
    }      
  }
  
  i = _param.getParam(_site, t_gpar::CLK, "");
  if ( i >= 0 )  {
    _param[i].value( _vBanc(4) );
    for( unsigned int jj = 1; jj <= _param.parNumber(); jj++ ) _Qx(i+1, jj) = 0.0;
    _Qx(i+1,i+1) = _clkStoModel->getQ() * _clkStoModel->getQ();
    if(_smooth) _Noise(i+1,i+1) = _clkStoModel->getQ() * _clkStoModel->getQ();      
  }   
  
  i = _param.getParam(_site, t_gpar::GLO_ISB, "");
  if ( i >= 0 )  {
    if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
      _Qx(i+1,i+1) = _sig_init_glo * _sig_init_glo;
    } else {
      if (_cntrep == 1) _Qx(i+1,i+1) += _gloStoModel->getQ();
      if(_smooth) _Noise(i+1,i+1)  = _gloStoModel->getQ();
    }      
  }
  
  i = _param.getParam(_site, t_gpar::GLO_IFCB, "");
  if ( i >= 0 )  {
    if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
      _Qx(i+1,i+1) = _sig_init_glo * _sig_init_glo;
    } else {
      if (_cntrep == 1) _Qx(i+1,i+1) += _gloStoModel->getQ();
      if(_smooth) _Noise(i+1,i+1)  = _gloStoModel->getQ();
    }      
  }   
   
  i = _param.getParam(_site, t_gpar::GAL_ISB, "");
  if ( i >= 0 )  {
    if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
      _Qx(i+1,i+1) = _sig_init_gal * _sig_init_gal;
    } else {
      if (_cntrep == 1) _Qx(i+1,i+1) += _galStoModel->getQ();
      if(_smooth) _Noise(i+1,i+1)  = _galStoModel->getQ();
    }      
  }
   
  i = _param.getParam(_site, t_gpar::BDS_ISB, "");
  if ( i >= 0 )  {
    if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
      _Qx(i+1,i+1) = _sig_init_bds* _sig_init_bds;
    } else {
      if (_cntrep == 1) _Qx(i+1,i+1) += _bdsStoModel->getQ();
      if(_smooth) _Noise(i+1,i+1)  = _bdsStoModel->getQ();
    }      
  }   
   
  i = _param.getParam(_site, t_gpar::QZS_ISB, "");
  if ( i >= 0 )  {
    if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
      _Qx(i+1,i+1) = _sig_init_qzs * _sig_init_qzs;
    } else {
      if (_cntrep == 1) _Qx(i+1,i+1) += _qzsStoModel->getQ();
      if(_smooth) _Noise(i+1,i+1)  = _qzsStoModel->getQ();
    }      
  }   
   
  i = _param.getParam(_site, t_gpar::TRP, "");
  if ( i >= 0 )  {
    if (!_initialized || _Qx(i + 1, i + 1) == 0.0){
      // for determination elipsoidal crd of receiver
      t_gtriple Ell, XYZ;
      
      if (_param.getCrdParam(_site, XYZ) > 0){
      }else if ( _valid_crd_xml ) { XYZ = _grec->crd_arp(_epoch); }
      xyz2ell( XYZ, Ell, false);
      
      if(_gModel->tropoModel() != 0){
        _param[i].apriori(_gModel->tropoModel()->getZHD(Ell, _epoch) );  
      }
      
      if ( _valid_ztd_xml ) {
        _param[i].value( _aprox_ztd_xml - _gModel->tropoModel()->getZHD(Ell, _epoch) );
      }else{
        if(_gModel->tropoModel() != 0){
          _param[i].value(_gModel->tropoModel()->getZWD(Ell, _epoch) );  
        }
      }
      
      _Qx(i+1,i+1) = ztdInit * ztdInit;
    }else{
      if( _epoch > _ztd_begStat && _epoch < _ztd_endStat){
        if(_smooth) _Noise(i+1,i+1) = 0.0;
      }else{
        if(_cntrep == 1)  _Qx(i+1,i+1) += _trpStoModel->getQ();
        if(_smooth) _Noise(i+1,i+1) = _trpStoModel->getQ(); 
      }
    }    
  }

  // add/remove inter system bias
  _syncSys();

  // add/remove ionosphere delay
  _syncIono();
  
  // add/remove inter-frequency biases
  _syncIFB();  

  // predict ionosphere delay  
  vector<t_gsatdata>::iterator it;
  for( it = _data.begin(); it != _data.end(); it++){  
    i = _param.getParam(_site, t_gpar::SION, it->sat());
    if ( i >= 0 ) {
      if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
        _Qx(i + 1, i + 1) = _sig_init_vion * _sig_init_vion;
      } else {
        if (_cntrep == 1) _Qx(i + 1, i + 1) += 1;
        if (_smooth) _Noise(i + 1, i + 1) = 1;
      }
    }
    i = _param.getParam(_site, t_gpar::VION, it->sat());
    if ( i >= 0 ) {      
      if(_gion) {   // Introduce new apriory value from IONEX if available
        t_gtriple site_xyz(0.0, 0.0, 0.0);
        t_gtriple site_ell(0.0, 0.0, 0.0);
        t_gtriple ipp_ell(0.0, 0.0, 0.0);
          
        _param.getCrdParam(_site, site_xyz);
        xyz2ell(site_xyz, site_ell, false);          
        ell2ipp(*it, site_ell, ipp_ell);
        
        double newval = _gion->iono(ipp_ell[0] * R2D, ipp_ell[1] * R2D, it->epoch());
        newval = ( newval * 40.28 * 1e16 ) / (G01_F*G01_F);
        
        _param[i].value(newval);
      }
      
      if (!_initialized || _Qx(i + 1, i + 1) == 0.0) {
        _Qx(i + 1, i + 1) =  _sig_init_vion * _sig_init_vion;
      }
      else {                
        if (_cntrep == 1) _Qx(i + 1, i + 1) += _ionStoModel->getQ() * _ionStoModel->getQ(); 
        if (_smooth) _Noise(i + 1, i + 1) = _ionStoModel->getQ() * _ionStoModel->getQ();
      }
    }         
  }   
  
  // predict inter-frequency clock bias
  for( it = _data.begin(); it != _data.end(); it++){
    i = _param.getParam(_site, t_gpar::IFCB_F3, it->sat());
    if ( i >= 0 ) {
      if (_cntrep == 1) _Qx(i+1,i+1) += 0.001;
      if(_smooth) _Noise(i+1,i+1)  = 0.001;
    }
    i = _param.getParam(_site, t_gpar::IFCB_F4, it->sat());
    if ( i >= 0 ) {
      if (_cntrep == 1) _Qx(i+1,i+1) += 0.001;
      if(_smooth) _Noise(i+1,i+1)  = 0.001;
    }
    i = _param.getParam(_site, t_gpar::IFCB_F5, it->sat());
    if ( i >= 0 ) {
      if (_cntrep == 1) _Qx(i+1,i+1) += 0.001;
      if(_smooth) _Noise(i+1,i+1)  = 0.001;
    }    
  }
  
}
   
// Restore state and covariance matrix
// ---------------------
void t_gsppflt::_restore(const SymmetricMatrix& QsavBP, const t_gallpar& XsavBP)
{
  _Qx    = QsavBP;
  _param = XsavBP;
}   
       
// Satelite position
// ----------
int t_gsppflt::_satPos(t_gtime& epo, t_gsatdata& gsatdata)
{
  gtrace("t_gsppflt::_satPos");

  string satname  =  gsatdata.sat();

  int i;
  i = gsatdata.addprd(_gnav); //add sat crd and clk
  if( i < 0 ) return i;

  return 1;     
}

// Get satellite number in epoch for specific system
// -----------------------------
int t_gsppflt::_numSat(GSYS gsys)
{
  int num = 0;
   
  vector<t_gsatdata>::iterator it;
  for (it = _data.begin(); it != _data.end(); it++){
    GSYS satSys = it->gsys();
    if (satSys == gsys) num++;
    else continue;
  }
   
  return num;
}

// Update time for RDW processes
// ---------------------------
void t_gsppflt::_timeUpdate(const t_gtime& epo)
{
  // tropo delay  
  _trpStoModel->updateTime(epo);
  _trpStoModel->setTcurr(epo);

  // iono delay  
  //_ionStoModel->updateTime(epo);
  //_ionStoModel->setTcurr(epo);  
   
  // GLO ISB
  _gloStoModel->updateTime(epo);
  _gloStoModel->setTcurr(epo);
  
  // GAL ISB
  _galStoModel->updateTime(epo);
  _galStoModel->setTcurr(epo);
  
  // BDS ISB
  _bdsStoModel->updateTime(epo);
  _bdsStoModel->setTcurr(epo);
   
   // QZS ISB
  _qzsStoModel->updateTime(epo);
  _qzsStoModel->setTcurr(epo);   
} 

// Sync inter-GNSS systems bias
// -----------------------------------
void t_gsppflt::_syncSys()
{
  gtrace("t_gsppflt::_syncSys() " + _epoch.str_hms());
   
//  set<string> mapPRN; 
  bool obsGps = false;
  bool obsGlo = false;   
  bool obsGal = false;
  bool obsBds = false;
  bool obsQzs = false;   
  vector<t_gsatdata>::iterator it;
  for( it = _data.begin(); it != _data.end(); it++){       // loop over all observations

    if (it->gsys() == GPS) obsGps = true;          
    if (it->gsys() == GLO) obsGlo = true;
    if (it->gsys() == GAL) obsGal = true;
    if (it->gsys() == BDS) obsBds = true;
    if (it->gsys() == QZS) obsQzs = true;     
  }   
  
  bool onlyGlo = false;
  bool onlyGal = false;
  bool onlyBds = false;
  bool onlyQzs = false;
  
  if(obsGlo && !obsGps            && !obsGal && !obsBds && !obsQzs) onlyGlo = true;
  if(obsGal && !obsGps && !obsGlo            && !obsBds && !obsQzs) onlyGal = true;
  if(obsBds && !obsGps && !obsGlo && !obsGal            && !obsQzs) onlyBds = true;
  if(obsQzs && !obsGps && !obsGlo && !obsGal && !obsBds           ) onlyQzs = true;   
  
   
  bool parGlo = false;   
  bool parGal = false;
  bool parBds = false;
  bool parQzs = false; 
  for (unsigned int i = 0; i < _param.parNumber(); i++){
     
    if (_param[i].parType == t_gpar::GLO_ISB) parGlo = true;     
    if (_param[i].parType == t_gpar::GAL_ISB) parGal = true;
    if (_param[i].parType == t_gpar::BDS_ISB) parBds = true;
    if (_param[i].parType == t_gpar::QZS_ISB) parQzs = true;     
  }        

  // Add GLO ISB parameter
  if ( !parGlo && obsGlo && !onlyGlo){
    t_gpar newPar(_data.begin()->site(), t_gpar::GLO_ISB, _param.parNumber()+1, "");
    newPar.value( 0.0 );
    _param.addParam(newPar);
    Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
    _Qx(_param.parNumber(), _param.parNumber()) = _sig_init_glo * _sig_init_glo; 
#ifdef DEBUG   
       cout << "GLO_ISB was added!" << " Epoch: " 
            << _epoch.str_ymdhms() << endl;
#endif
  }      
   
  // Add GAL ISB parameter
  if ( !parGal && obsGal && !onlyGal){
    t_gpar newPar(_data.begin()->site(), t_gpar::GAL_ISB, _param.parNumber()+1, "");
    newPar.value( 0.0 );
    _param.addParam(newPar);
    Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
    _Qx(_param.parNumber(), _param.parNumber()) = _sig_init_gal * _sig_init_gal; 
#ifdef DEBUG   
    cout << "GAL_ISB was added!" << " Epoch: " 
         << _epoch.str_ymdhms() << endl;
#endif
  }

  // Add BDS ISB parameter
  if ( !parBds && obsBds && !onlyBds){
    t_gpar newPar(_data.begin()->site(), t_gpar::BDS_ISB, _param.parNumber()+1, "");
    newPar.value( 0.0 );
    _param.addParam(newPar);
    Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
    _Qx(_param.parNumber(), _param.parNumber()) = _sig_init_bds * _sig_init_bds; 
#ifdef DEBUG   
       cout << "BDS_ISB was added!" << " Epoch: " 
            << _epoch.str_ymdhms() << endl;
#endif
  }

  // Add QZS ISB parameter
  if ( ! parQzs && obsQzs && !onlyQzs){
    t_gpar newPar(_data.begin()->site(), t_gpar::QZS_ISB, _param.parNumber()+1, "");
    newPar.value( 0.0 );
    _param.addParam(newPar);
    Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
    _Qx(_param.parNumber(), _param.parNumber()) = _sig_init_qzs * _sig_init_qzs;
#ifdef DEBUG   
       cout << "QZS_ISB was added!" << " Epoch: " 
            << _epoch.str_ymdhms() << endl;
#endif
  }      

  // Remove GLO ISB paremeter
  if (parGlo && ! obsGlo){
    int i = _param.getParam(_data.begin()->site(), t_gpar::GLO_ISB, "");
    Matrix_remRC(_Qx, _param[i].index, _param[i].index);
    _param.delParam(i);
    _param.reIndex();
#ifdef DEBUG   
    cout << "GLO_ISB was removed!" << " Epoch: " 
         << _epoch.str_ymdhms() << endl;
#endif     
  }   
   
  // Remove GAL ISB paremeter
  if (parGal && ! obsGal){
    int i = _param.getParam(_data.begin()->site(), t_gpar::GAL_ISB, "");
    Matrix_remRC(_Qx, _param[i].index, _param[i].index);
    _param.delParam(i);
    _param.reIndex();
#ifdef DEBUG   
    cout << "GAL_ISB was removed!" << " Epoch: " 
         << _epoch.str_ymdhms() << endl;
#endif     
  }

  // Remove BDS ISB paremeter
  if (parBds && ! obsBds){
    int i = _param.getParam(_data.begin()->site(), t_gpar::BDS_ISB, "");
    Matrix_remRC(_Qx, _param[i].index, _param[i].index);
    _param.delParam(i);
    _param.reIndex();     
#ifdef DEBUG   
       cout << "BDS_ISB was removed!" << " Epoch: " 
            << _epoch.str_ymdhms() << endl;
#endif     
  }   

  // Remove QZS ISB paremeter
  if (parQzs && ! obsQzs){
    int i = _param.getParam(_data.begin()->site(), t_gpar::QZS_ISB, "");
    Matrix_remRC(_Qx, _param[i].index, _param[i].index);
    _param.delParam(i);
    _param.reIndex();
#ifdef DEBUG   
       cout << "QZS_ISB was removed!" << " Epoch: " 
            << _epoch.str_ymdhms() << endl;
#endif     
  }   

} // void     

// Add/Remove ionosphere delay
void t_gsppflt::_syncIono()
{
  gtrace("t_gpppflt::_syncIono");
  
  // Add ionosphere parameter and appropriate rows/columns covar. matrix   
  set<string> mapPRN; 
   
  vector<t_gsatdata>::iterator it;
  for( it = _data.begin(); it != _data.end(); it++){       // loop over all tracked satellites
    mapPRN.insert( it->sat() );
    
    if(_observ == IONO_FREE) return;
    
    if(_observ == RAW_DOUBLE || _observ == RAW_ALL){
/*      
      // add ionosphere slant delay - not in case of single frequency
      if( _param.getParam(_site, t_gpar::SION, it->sat()) < 0 ){
        t_gpar parSION(it->site(), t_gpar::SION, _param.parNumber()+1, it->sat());
        parSION.value(5.0);
        _param.addParam(parSION);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = 10 * 10;
      }
*/
      // add ionosphere vertical delay - not in case of single frequency
      if( _param.getParam(_site, t_gpar::VION, it->sat()) < 0 && _iono_est){
        t_gpar parVION(it->site(), t_gpar::VION, _param.parNumber()+1, it->sat());
        if(_gion) {
          t_gtriple site_xyz(0.0, 0.0, 0.0);
          t_gtriple site_ell(0.0, 0.0, 0.0);
          t_gtriple ipp_ell(0.0, 0.0, 0.0);
          
          _param.getCrdParam(_site, site_xyz);
          xyz2ell(site_xyz, site_ell, false);          
          ell2ipp(*it, site_ell, ipp_ell);
          
          double newval = _gion->iono(ipp_ell[0] * R2D, ipp_ell[1] * R2D, it->epoch());
          newval = ( newval * 40.28 * 1e16 ) / (G01_F*G01_F);
          
          parVION.value(newval);
        }else parVION.value(2);  // apriory VTEC delay = 2 m
        _param.addParam(parVION);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = _sig_init_vion * _sig_init_vion;
        
//cout << "Pridavam VION " << it->sat() << " " << it->epoch().str_hms() << endl;
      }
      
    }
    
  }  // end loop over all observations
  
  // Remove params and appropriate rows/columns covar. matrix
  for (unsigned int i=0; i<=_param.parNumber()-1; i++) {   
    if (_param[i].parType == t_gpar::SION   ||
        _param[i].parType == t_gpar::VION ) {
      
      string sat = _param[i].prn;
      
      set<string>::iterator prnITER = mapPRN.find(sat);
      if (prnITER == mapPRN.end()) {   

        Matrix_remRC(_Qx, _param[i].index, _param[i].index);
        _param.delParam(i);
        _param.reIndex();
        i--;
      }    
    } 
  }   
  
  return;
}
   
// Add/Remove inter-freq. biases
void t_gsppflt::_syncIFB()
{
  gtrace("t_gpppflt::_syncIFB");
  
  set<string> mapPRN; 
   
  vector<t_gsatdata>::iterator it;
  for( it = _data.begin(); it != _data.end(); it++){       // loop over all tracked satellites
    mapPRN.insert( it->sat() );
    
    if(_observ == IONO_FREE) return;
    
    if(_observ == RAW_ALL){

      // add inter-frequency clock bias for FREQ_3
      if(it->contain_freq(FREQ_3) && _param.getParam(_site, t_gpar::IFCB_F3, it->sat()) < 0 ){
        t_gpar parIFCB(it->site(), t_gpar::IFCB_F3, _param.parNumber()+1, it->sat());
        parIFCB.value(0.0);
        _param.addParam(parIFCB);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = 1 * 1;
      }
      
      // add inter-frequency clock bias for FREQ_4
      if(it->contain_freq(FREQ_4) && _param.getParam(_site, t_gpar::IFCB_F4, it->sat()) < 0 ){
        t_gpar parIFCB(it->site(), t_gpar::IFCB_F4, _param.parNumber()+1, it->sat());
        parIFCB.value(0.0);
        _param.addParam(parIFCB);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = 1 * 1;
      }
      
      // add inter-frequency clock bias for FREQ_5
      if(it->contain_freq(FREQ_5) && _param.getParam(_site, t_gpar::IFCB_F5, it->sat()) < 0 ){
        t_gpar parIFCB(it->site(), t_gpar::IFCB_F5, _param.parNumber()+1, it->sat());
        parIFCB.value(0.0);
        _param.addParam(parIFCB);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = 1 * 1;
      }      
      
      // add inter-frequency code bias for FREQ_3 (if not already)
      if(!_ifb3_init && it->contain_freq(FREQ_3) && _param.getParam(_site, t_gpar::IFB_C3, "") < 0 ){ 
        t_gpar parIFB(it->site(), t_gpar::IFB_C3, _param.parNumber()+1, "");
        parIFB.value(0.0);
        _param.addParam(parIFB);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = 100 * 100;
        _ifb3_init = true;
      }      
      
      // add inter-frequency code bias for FREQ_4 (if not already)
      if(!_ifb4_init && it->contain_freq(FREQ_4) && _param.getParam(_site, t_gpar::IFB_C4, "") < 0 ){
        t_gpar parIFB(it->site(), t_gpar::IFB_C4, _param.parNumber()+1, "");
        parIFB.value(0.0);
        _param.addParam(parIFB);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = 100 * 100;
        _ifb4_init = true;
      }
      
      // add inter-frequency code bias for FREQ_5 (if not already)
      if(!_ifb5_init && it->contain_freq(FREQ_5) && _param.getParam(_site, t_gpar::IFB_C5, "") < 0 ){
        t_gpar parIFB(it->site(), t_gpar::IFB_C5, _param.parNumber()+1, "");
        parIFB.value(0.0);
        _param.addParam(parIFB);
        Matrix_addRC(_Qx, _param.parNumber(), _param.parNumber());
        _Qx(_param.parNumber(), _param.parNumber()) = 100 * 100;
        _ifb5_init = true;
      }      
      
    }
    
  }  // end loop over all satellites
  
  // Remove params and appropriate rows/columns covar. matrix
  for (unsigned int i=0; i<=_param.parNumber()-1; i++) {   
    if (_param[i].parType == t_gpar::IFCB_F3 ||
        _param[i].parType == t_gpar::IFCB_F4 ||
        _param[i].parType == t_gpar::IFCB_F5  ) {
      
      string sat = _param[i].prn;
      
      set<string>::iterator prnITER = mapPRN.find(sat);
      if (prnITER == mapPRN.end()) {   

        Matrix_remRC(_Qx, _param[i].index, _param[i].index);
        _param.delParam(i);
        _param.reIndex();
        i--;
      }    
    } 
  }   
  
  return;
}   
      
// get t_gobs instance for particular prn, band, and attribute   
int t_gsppflt::_getgobs(string prn, GOBSTYPE type, GOBSBAND band, t_gobs& gobs)
{
  map<string, map<GOBSBAND, GOBSATTR>>::iterator itPRN = _signals.find(prn);
  if( itPRN == _signals.end() ) {
    _signals[prn][band] = ATTR;
    itPRN = _signals.find(prn);
  }   
  
  map<GOBSBAND, GOBSATTR>::iterator itBAND = itPRN->second.find(band);
  if( itBAND == itPRN->second.end() ) {
    _signals[prn][band] = ATTR;
  }
  
  GOBSATTR attr = _signals[prn][band];
        
  gobs.type(type);
  gobs.band(band);
  gobs.attr(attr);

  return 1;
}
   
// save observations residuals   
void t_gsppflt::_save_residuals(ColumnVector& v, vector<t_gsatdata>& data, RESIDTYPE restype)
{
  int mult = 1;
  if(_phase) mult = 2;
  int k = 0;   
  
  for(auto it = data.begin(); it != data.end(); it++){
    
    vector<double> res;
    
    auto itFrq = _frqNum.find(it->sat());
    if(itFrq == _frqNum.end()) continue;
    
    // separate frequencies     
    for(int i = 1; i <= itFrq->second * mult; i++) {      
      res.push_back( v(i+k) );
    }
    k += itFrq->second * mult;
    
    // clear all previously stored residuals
    it->clear_res(restype);
    
    // separate code and phase and storing
    unsigned int sz = res.size();
    if(_phase) {
      for(unsigned int i = 0; i < sz/2; i++)  it->addres(restype, TYPE_C, res[i]);
      for(unsigned int i = sz/2; i < sz; i++) it->addres(restype, TYPE_L, res[i]);
    }else for(unsigned int i = 0; i < sz; i++) it->addres(restype, TYPE_C, res[i]);
    
  }      
}
   
// Apply DCB correction
int t_gsppflt::_applyDCB(t_gsatdata& satdata, double& P, t_gobs* gobs1, t_gobs* gobs2)
{ 
  GSYS gsys = satdata.gsys();
  GOBS g1_R2, g1_R3, g1_ref;
  GOBS g2_R2, g2_R3, g2_ref;;
  g1_R2 = g1_R3 = g2_R2 = g2_R3 = X;
  if(gsys == GPS || gsys == GLO){
    g1_R2 = P1;  g1_R3 = C1W;
    g2_R2 = P2;  g2_R3 = C2W;
  }else if(gsys == GAL){
    g1_R2 = C1;  g1_R3 = C1X;
    g2_R2 = C5;  g2_R3 = C5X;
  }else if(gsys == BDS){
    g1_R2 = C1;  g1_R3 = C1I;
    g2_R2 = C7;  g2_R3 = C7I;
  }
  
  GOBS g1 = X; 
  GOBS g2 = X;
    
  if(gobs1->attr() == ATTR) g1 = satdata.id_range(gobs1->band());  // automatic GOBS selection
  else g1 = gobs1->gobs();                                         // specific GOBS

  if(gobs2 != 0){
    if(gobs2->attr() == ATTR) g2 = satdata.id_range(gobs2->band());  // automatic GOBS selection
    else g2 = gobs2->gobs();                                         // specific GOBS  
  }
  
  FREQ_SEQ freq1 = t_gsys::band2freq(satdata.gsys(), gobs1->band());
  
  double corr1 = 0.0;
  double corr2 = 0.0;
  
  t_gtime epoch = satdata.epoch();
  string sat    = satdata.sat();

  if(g1 < 1000) { g1_ref = g1_R3; g2_ref = g2_R3; }
  else          { g1_ref = g1_R2; g2_ref = g2_R2; }
  
  if(_gallbias){
    if(freq1 == FREQ_1 && g1 != g1_ref) {
      corr1 = _gallbias->get(epoch, sat, g1_ref, g1);      
    }
    if(freq1 == FREQ_2 && g1 != g2_ref) {
      corr1 =  _gallbias->get(epoch, sat, g2_ref, g1);    
    }

    if(gobs2){
      FREQ_SEQ freq2 = t_gsys::band2freq(satdata.gsys(), gobs2->band());
      if(freq2 == FREQ_1 && g2 != g1_ref) {
        corr2 = _gallbias->get(epoch, sat, g1_ref, g2);        
      }
      if(freq2 == FREQ_2 && g2 != g2_ref) {
        corr2 = _gallbias->get(epoch, sat, g2_ref, g2);        
      }
    }
  }else{
    if(_log) _log->comment(1, "gsppflt", "Warning: Code biases container not available!");
  }
      
  if(gobs2 != 0){
    double alfa, beta;
    alfa = beta = 0.0;
    satdata.coef_ionofree(gobs1->band(), alfa, gobs2->band(), beta);
    // cout << "Opravuji Codo obs: " << satdata.sat() << " " << satdata.epoch().str_ymdhms() << " " 
    //      << gobs2str(g1) << " " << corr1/CLIGHT*1e9 << " " << gobs2str(g2) << " " << corr2/CLIGHT*1e9 << endl;
    if(_observ == IONO_FREE) P += alfa * corr1 + beta * corr2;    
  }else {
    P += corr1;
  }
  
  return 1;
}
   
} // namespace
