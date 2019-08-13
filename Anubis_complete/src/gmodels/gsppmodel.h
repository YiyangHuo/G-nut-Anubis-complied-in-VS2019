
#ifndef GSPPMODEL_H
#define GSPPMODEL_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: various SPP models
  Version: $ Rev: $

  2014-11-27 /PV: created

-*/

#include <string>
#include <map>
#include <cmath>

#include "gmodels/gmodel.h"

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetrec.h"
#include "gutils/gsysconv.h"
#include "gutils/gnss.h"
#include "gmodels/ggmf.h"
#include "gutils/gnss.h"

using namespace std;

namespace gnut {   
         
class t_gsppmodel : public t_gmodel
{
 public:
   t_gsppmodel(string site,t_glog* glog, t_gsetbase* settings);
   t_gsppmodel();
  virtual ~t_gsppmodel();

   virtual int outlierDetect(vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, const ColumnVector&);  // OLD
   virtual int outlierDetect(vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&);
   virtual int outlierDetect_chi(  vector<t_gsatdata>& data, SymmetricMatrix& Qx, const SymmetricMatrix&, const ColumnVector&);  // OLD 
   
   virtual double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata&, t_gobs& gobs, bool com = false);
   double tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& satdata);   
   double ionoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple site_ell, t_gsatdata& satdata, t_gobs& gobs);      
      
 protected:
   double               _maxres(const ColumnVector& v, GSYS gs, bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA);  // OLD
   double               _maxres(bool phase, vector<t_gsatdata>& data, vector<t_gsatdata>::iterator& itDATA, RESIDTYPE res_type, GSYS gs = GNS);
   bool                 _check_outl(bool phase, double& maxres, vector<t_gsatdata>::iterator& itData, vector<t_gsatdata>& data); // OLD
   bool                 _check_outl(bool phase, double& maxresNORM, vector<t_gsatdata>::iterator& itDataNORM,
 									                  double& maxresORIG, vector<t_gsatdata>::iterator& itDataORIG,
									                  vector<t_gsatdata>::iterator& itDataErase, vector<t_gsatdata>& data);
   void                 _logOutl(bool phase, string prn, int data_size, double maxres, double ele, t_gtime epo, RESIDTYPE resid_type);
   vector<ColumnVector> _devide_res(const ColumnVector& v_orig);
   
   map<GSYS, double>   _maxres_C;
   map<GSYS, double>   _maxres_L;
   double              _maxres_norm;
   shared_ptr<t_gobj>  _grec;  
   t_gtriple           _antennal_height;
   TROPMODEL           _trpModStr;
   ZTDMPFUNC           _tropo_mf;   
   RESIDTYPE           _resid_type;
   OBSCOMBIN           _observ;
   
};

} // namespace

#endif //  GSPPMODEL_H
