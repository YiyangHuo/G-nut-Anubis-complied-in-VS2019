
#ifndef GSPPFLT_H
#define GSPPFLT_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements spp client
  Version: $ Rev: $

  2014-11-24 /PV: created

-*/

#include "gproc/gspp.h"
#include "gmodels/gdop.h"
#include "gproc/gflt.h"
#include "gall/gallpar.h"
#include "gall/gallbias.h"
#include "gall/gallprod.h"
#include "gutils/gsysconv.h"
#include "gmodels/gstochasticmodel.h"
#include "gset/gsetflt.h"

namespace gnut {

#define SPP_MINSAT  (6)    // minimum number of satellites

class t_gsppflt : public virtual t_gspp
{
 public:
   t_gsppflt(string mark, t_gsetbase* set);
   virtual ~t_gsppflt();
   
   virtual int processBatch(const t_gtime& beg, const t_gtime& end);

   void minsat(size_t minsat){ _minsat = (minsat<5)?5:minsat; }
   
 protected:
   virtual void   _predict();
   virtual void   _restore(const SymmetricMatrix& Qsav, const t_gallpar& Xsav);
   virtual int    _satPos(t_gtime&, t_gsatdata&);
   virtual int    _prepareData();
   virtual int    _apply_tides(t_gtime& _epoch, t_gtriple& xRec);
   virtual int    _processEpoch(const t_gtime& runEpoch);      
           int   _addObsP(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P);
   virtual int   _addObsL(t_gsatdata& satdata, unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P);
   virtual int   _addPseudoZTD(unsigned int& iobs, t_gtriple& ell, Matrix& A, ColumnVector& l, DiagonalMatrix& P);
   
         double   _weightObs(t_gsatdata& satdata, t_gobs& go);

   virtual void   _timeUpdate(const t_gtime& epo);
           void   _syncSys();
           void   _syncIono();
           void   _syncIFB();   
           void   _save_residuals(ColumnVector& v, vector<t_gsatdata>& satdata, RESIDTYPE restype);
		   double _applySICB(string prn, double elev, GOBSBAND freq);
           int    _applyDCB(t_gsatdata& satdata, double& P, t_gobs* gobs1, t_gobs* gobs2 = 0);
   
   vector<t_gsatdata> _data;
   map<string, int>   _newAMB;
   set<string>        _slips;
   unsigned int       _minsat;
   double             _sig_unit;   
   
    
// Models   
   t_randomwalk*       _trpStoModel;
   t_whitenoise*       _ionStoModel;
   //t_randomwalk*       _ionStoModel;
   t_randomwalk*       _gloStoModel;
   t_randomwalk*       _galStoModel;
   t_randomwalk*       _bdsStoModel;
   t_randomwalk*       _qzsStoModel;
   t_whitenoise*       _clkStoModel;
   t_whitenoise*       _crdStoModel;    
   
   
// Parameters and covariance matrix   
   t_gallpar           _param;
   SymmetricMatrix     _Qx;   
   
// Noise matrix
   DiagonalMatrix      _Noise;
   
// Epoch time
   t_gtime           _epoch;

// Estimation objects   
   t_gflt*           _filter;
   
   int               _numSat(GSYS gsys);   
   ColumnVector      _vBanc; 
   t_gdop            _dop; 
   int               _cntrep;      
   bool              _smooth;
   unsigned int      _n_NPD_flt;
   unsigned int      _n_ALL_flt;
   unsigned int      _n_NPD_smt;
   unsigned int      _n_ALL_smt;
   
   RESIDTYPE         _resid_type;   
   
   int                  _getgobs(string prn, GOBSTYPE type, GOBSBAND band, t_gobs& gobs);
   map<string, int>     _frqNum;
   map<string, t_gtime> _lastEcl;
   
   bool                 _ifb3_init;
   bool                 _ifb4_init;
   bool                 _ifb5_init;
   
   bool                 _auto_band;
};

} // namespace

#endif // GSPPFLT_H
