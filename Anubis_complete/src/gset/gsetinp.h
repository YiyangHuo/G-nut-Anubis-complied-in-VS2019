
#ifndef GSETINP_H
#define GSETINP_H

#define XMLKEY_INP "inputs"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements input setting class
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <map>
#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut {

// the order is important here!
enum IFMT {
     FMT_INP,
     IONEX_INP,
     SP3_INP,
     SP3_REF,
     SP3_PRD,     
     ATX_INP,      // need to be before all OBSERVATIONS to create objects
     RINEXC_INP,
     RINEXC_REF,
     RINEXC_PRD,     
     RINEXM_INP,
     RINEXN_INP,
     RINEXN_REF,
     RINEXN_PRD,
     RINEXO_INP,   // need to be after RINEXN (for GLONASS chanels)!
     BNCOBS_INP,   // need to be after RINEXN (for GLONASS chanels)!
     TROSINEX_INP,
     TROSINEX0_INP,
     RTCM_INP,
     TSMET_INP,
     BLQ_INP,
     BNCRTCM_INP,  // need to be after RINEXN to fill the corrections !
     BNCRTCM_REF,
     BNCRTCM_PRD,     
     TROGRID_INP,
     FCB_INP,

     GRIB_CHMI_INP,
     GRIB_HARM_INP,
     GRIB_ERA_INP,
     NCD3_ICS_INP,

     RAO_TEXT_INP,
     RAO_IGRA_INP,
     RAO_BADC_INP,
  
     BIASINEX_INP,
     BIABERN_INP
     
};


class t_gsetinp : public virtual t_gsetbase
{
 public:
   t_gsetinp();
  ~t_gsetinp();
  
   void check();                                 // settings check
   void help();                                  // settings help

   static IFMT   str2ifmt( const string& s );
   static string ifmt2str( const IFMT&   f );
   
   int input_size(const string& fmt);            // get format input size
   multimap<IFMT,string> inputs_all();           // get format inputs (all in multimap)
   vector<string> inputs(const string& fmt);     // get format inputs (ordered)
   set   <string> iformats();                    // get formats       (ordered)
   string corrStream();                          // get rtcm correction stream

   bool chkNavig();                              // checking navigation Rinex
   bool chkHealth();                             // checking satellite healthy status
   
 protected:
   vector<string> _inputs(const string& fmt);
   set   <string> _iformats();

   set<IFMT> _IFMT_supported;                    // vector of supported IFMTs (app-specific)
   
   bool      _chkNavig;
   bool      _chkHealth;

   string    _corrStream;
 private:

};

} // namespace

#endif
