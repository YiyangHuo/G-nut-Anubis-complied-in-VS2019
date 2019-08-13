
#ifndef GPRODMET_H
#define GPRODMET_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: class storing meteorological data
  Version: $Rev:$

  2016-04-14 /JD: created

-*/

#include <map>
#include <iostream>

#include "gprod/gprod.h"
#include "gmodels/gnwmbase.h"

using namespace std;

namespace gnut {

enum METEO_ID { METEO_PRES, METEO_EPRE, METEO_RHUM, METEO_SHUM, METEO_TEMP, METEO_MTEMP, METEO_TDEW,
                METEO_ZHD,  METEO_ZWD,  METEO_ZTD,  METEO_IWV,
                METEO_DT,   METEO_DM,   METEO_DE,   METEO_DW,   METEO_DI,
                METEO_UNKNOWN };

class t_gprodmet : public t_gprod
{
 public:
   t_gprodmet(const t_gtime& t, shared_ptr<t_gobj> pt = nullobj);
  ~t_gprodmet();

  virtual set<METEO_ID> get_params();
   
  virtual void   met( const METEO_ID key, const double& val, const double& rms = 0.0);
  virtual double met( const METEO_ID key);
  virtual double rms( const METEO_ID key);

  virtual void   pres( const double& val, const double& rms = 0.0);   // [hPa]
  virtual double pres();
  virtual double pres_rms();

  virtual void   epre( const double& val, const double& rms = 0.0);   // [hPa]
  virtual double epre();
  virtual double epre_rms();

  virtual void   rhum( const double& val, const double& rms = 0.0);   // [%] 
  virtual double rhum();
  virtual double rhum_rms();

  virtual void   shum( const double& val, const double& rms = 0.0);   // [g/kg]
  virtual double shum();
  virtual double shum_rms();

  virtual void   temp( const double& val, const double& rms = 0.0);   // [K]
  virtual double temp();
  virtual double temp_rms();

  virtual void   mtemp( const double& val, const double& rms = 0.0);   // [K]
  virtual double mtemp();
  virtual double mtemp_rms();

  virtual void   tdew( const double& val, const double& rms = 0.0);   // [K]
  virtual double tdew();
  virtual double tdew_rms();

  virtual void   ztd( const double& val, const double& rms = 0.0);   // [mm]
  virtual double ztd();
  virtual double ztd_rms();

  virtual void   zhd( const double& val, const double& rms = 0.0);   // [mm]
  virtual double zhd();
  virtual double zhd_rms();

  virtual void   zwd( const double& val, const double& rms = 0.0);   // [mm]
  virtual double zwd();
  virtual double zwd_rms();

  virtual void   iwv( const double& val, const double& rms = 0.0);   // [kg/m2]
  virtual double iwv();
  virtual double iwv_rms();

  virtual void   dtemp( const double& val, const double& rms = 0.0);  // [K/km]
  virtual double dtemp();
  virtual double dtemp_rms();

  virtual void   dmtemp( const double& val, const double& rms = 0.0);  // [K/km]
  virtual double dmtemp();
  virtual double dmtemp_rms();

  virtual void   dzwd( const double& val, const double& rms = 0.0);  // [-]
  virtual double dzwd();
  virtual double dzwd_rms();

  virtual void   depre( const double& val, const double& rms = 0.0);  // [-]
  virtual double depre();
  virtual double depre_rms();

  virtual void   diwv( const double& val, const double& rms = 0.0);  // [-]
  virtual double diwv();
  virtual double diwv_rms();

 protected:
	       
  map<METEO_ID,double>   _map_val;
  map<METEO_ID,double>   _map_rms;

};

} // namespace

#endif
