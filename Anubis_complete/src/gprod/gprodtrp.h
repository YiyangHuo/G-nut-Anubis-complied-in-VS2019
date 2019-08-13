
#ifndef GPRODTRP_H
#define GPRODTRP_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: class storing tropospheric (ZTD+GRD) products
  Version: $Rev:$

  2011-03-25 /JD: created

-*/

#include <iostream>

#include "gprod/gprod.h"
#include "gprod/gprodmet.h"
#include "gmodels/gnwmbase.h"
#include "gmodels/ggmf.h"
#include "gset/gsetproc.h"

using namespace std;

namespace gnut {

class t_gprodtrp : public t_gprodmet
{
 public:
   t_gprodtrp(const t_gtime& t, shared_ptr<t_gobj> pt = nullobj);
  ~t_gprodtrp();

//enum TRO_ENUM { ZTD, ZHD, ZWD, PRES, EPRE, TEMP, MTEMP, RELH, SHUM, IWV, DT, DE, DW }
//typedef map<TRO_ENUM, pair<double,double> t_map_met;
//typedef map<TRO_ENUM, double> t_map_rms;
//typedef map<string, pair<double,double> > t_map_met;

  virtual void   zhd( const double& val, const double& rms = 0.0);
  virtual double zhd();
  virtual double zhd_rms();
   
  virtual void   zwd( const double& val, const double& rms = 0.0);
  virtual double zwd();
  virtual double zwd_rms();

  virtual void   ztd( const double& val, const double& rms = 0.0);   
  virtual double ztd();
  virtual double ztd_rms();
   
  virtual void   grd( const double& grdN,       const double& grdE,
  	              const double& rmsN = 0.0, const double& rmsE = 0.0);
   
  virtual void   grd_N( const double& grd, const double& rms = 0.0);
  virtual double grd_N() const;
  virtual double grd_N_rms() const;

  virtual void   grd_E( const double& grd, const double& rms = 0.0);
  virtual double grd_E() const;
  virtual double grd_E_rms() const;   

  virtual void   ztd_mf(const ZTDMPFUNC& mf);
  virtual void   grd_mf(const GRDMPFUNC& mf);   
   
 protected:
  double    _grd_N;
  double    _grd_E;
  double    _grd_N_rms;
  double    _grd_E_rms;
  ZTDMPFUNC _ztd_mf;
  GRDMPFUNC _grd_mf;
   
  bool     _lock_ztd;
   
};

} // namespace

#endif
