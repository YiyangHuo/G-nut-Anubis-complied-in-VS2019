
#ifndef GSETQC_H
#define GSETQC_H

#define XMLKEY_QC "qc"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements data extraction setting
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut {

// POS_HEALTH is a minimum if INP chk_health=true
// order gives a priority sequence!   
enum USE_HEALTH { DEF_HEALTH, POS_HEALTH, STT_HEALTH, ALL_HEALTH };

class t_gsetqc : public virtual t_gsetbase
{
 public:
   t_gsetqc();
  ~t_gsetqc();
  
   void check();                                 // settings check
   void help();                                  // settings help

   int    summ();
   int    head();
   int    stat();
   int    gaps();
   int    band();
   int    prep();
   int    elev();
   int    mult();
   int    stnr();
   int    sinf();
   double step();
   int    tgap();
   int    tpcs();
   int    calc();
   int    nsat();
   int    gkpi();
   int    mp_nepochs();
   bool   mp_all();
   double mp_limit();
   double ele_cut();
   double pos_cut();
   int    pos_int();
   bool   pos_kin();
   bool   ele_app();
   bool   ele_new();
   bool   sat_rec();

   double dV_lim();
   double dH_lim();
   double dG_lim();   

   USE_HEALTH useHealth();
   
   USE_HEALTH str2useHealth(string s);

   string useHealth2str(USE_HEALTH UH);

 protected:
   
   int    _summ;                                 // summary statistics (verbosity)
   int    _head;                                 // check header metadata (verbosity)
   int    _stat;                                 // observation statistics (verbosity)
   int    _calc;                                 // Calculation aprox position
   int    _gaps;                                 // gap & pieces statistics (verbosity)
   int    _band;                                 // gsys band statistics (verbosity)
   int    _prep;                                 // cycle-slips/clock-jumps (verbosity)
   int    _elev;                                 // elevation/azimuth (verbosity)
   int    _mult;                                 // multipath (verbosity)
   int    _stnr;                                 // SNR (verbosity)
   int    _sinf;                                 // satellite info
   int    _gkpi;                                 // selected kpi (dedicated service)
   double _step;                                 // step for time spacing [s]
   int    _tgap;                                 // interval for gaps [s]
   int    _tpcs;                                 // interval for pieces [s]
   int    _nsat;                                 // number of columns for sat-specific reporting
   int    _mp_nepochs;                           // number of epochs for multipath estimation
   bool   _mp_all;                               // high-resolution mulitipath: estimate all or interpolate
   double _mp_limit;                             // sigma-multiplicator for MP cycle-slip & outlier detection
   double _ele_cut;                              // elevation mask (cut-off angle)
   double _pos_cut;                              // elevation mask for positioning
   int    _pos_int;                              // sampling interval for positioning
   bool   _pos_kin;                              // static/kinematic receiver (true == kinematic)
   bool   _ele_app;                              // approximated elevations (true = less precise)
   bool   _ele_new;                              // new strategy for cut-off + horizon estimates

   double _dV_lim;                               // Vertical position limit
   double _dH_lim;                               // Horizontal position limit
   double _dG_lim;                               // GDOP limit

   bool   _sat_rec;                              // expected satellite (all | receiver tracking)

   USE_HEALTH _use_health;                       // method of using satellite healthy status

 private:
};

} // namespace

#endif
