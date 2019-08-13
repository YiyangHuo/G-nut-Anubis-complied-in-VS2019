
#ifndef GSETFLT_H
#define GSETFLT_H

#define XMLKEY_FLT "filter"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements filter setting
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gsetbase.h"

using namespace std;

namespace gnut {

class t_gsetflt : public virtual t_gsetbase
{
 public:
   t_gsetflt();
  ~t_gsetflt();

   void check();                                  // settings check
   void help();                                   // settings help

   string method_flt();
   string method_smt();   
   double noise_clk();
   double noise_crd();
   
   double rndwk_glo();
   double rndwk_gal();
   double rndwk_bds();
   double rndwk_qzs();
   
   double rndwk_ztd();
   double noise_vion();
   double rndwk_grd();   
   int    reset_amb();
   int    reset_par();
   int    smt_delay();
   bool   smooth();

 protected:
   string _method_flt;                            // type of filtering method (kalman, SRCF)
   string _method_smt;                            // type of filtering method (kalman, SRCF)   
   double _noise_clk;                             // white noise for receiver clock [m]
   double _noise_crd;                             // white noise for coordinates [m]
   double _rndwk_glo;                             // random walk process for GLONASS system time offset
   double _rndwk_gal;                             // random walk process for Galileo system time offset
   double _rndwk_bds;                             // random walk process for BeiDou system time offset
   double _rndwk_qzs;                             // random walk process for ZQSS system time offset   
   double _rndwk_ztd;                             // random walk process for ZTD [mm/sqrt(hour)]
   double _noise_vion;                            // white noise process for VION [mm/sqrt(hour)]   
   double _rndwk_grd;                             // random walk process for tropo grad [mm/sqrt(hour)]   
   int    _reset_amb;                             // interval for reseting ambiguity [s]
   int    _reset_par;                             // interval for reseting CRD, ZTD, AMB [s]
   int    _delay;                                 // delay for backwart smoothing [s]
   bool   _smooth;                                // backward smoothing

 private:
};

} // namespace

#endif
