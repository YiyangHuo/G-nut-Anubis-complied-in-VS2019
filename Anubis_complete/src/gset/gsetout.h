
#ifndef GSETOUT_H
#define GSETOUT_H

#define XMLKEY_OUT "outputs"
#define DEFAULT_FILE_VER ""         // default version for file format
#define DEFAULT_FILE_UPD     0      // default auto update [min] for file saving
#define DEFAULT_FILE_LEN     0      // default auto length [min] for file saving
#define DEFAULT_FILE_SMP     0      // default auto sample [sec] for file saving
#define DEFAULT_FILE_OFF     0      // default file offset [min] for file saving
#define DEFAULT_FILE_SYS  "UTC"     // default file system       for file saving

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements output setting class
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut {

// the order is important here!
enum OFMT {
     XXX_OUT,
     TSA_OUT,
     MPT_OUT,     
     LOG_OUT,
     SUM_OUT,
     XTR_OUT,
     NAV_OUT,
     XML_OUT,
     XQC_OUT,
     PPP_OUT,
     FLT_OUT,
     SMT_OUT,
     LSQ_OUT,
     RES_OUT,
     GRD_OUT,
     SRF_OUT,
     FIT_OUT,
     CFG_OUT,
     RINEXN_OUT,
     RINEXN2_OUT,
     RINEXO_OUT,
     RINEXC_OUT,
     RTCM_OUT,
     BNCOBS_OUT,
     BNCRTCM_OUT,
     SINEX_OUT,
     TROSINEX_OUT,
     TROSINEX0_OUT,
     KML_OUT,
     PRDQC_OUT
};

class t_gsetout : public virtual t_gsetbase
{
 public:
   t_gsetout();
  ~t_gsetout();

   static OFMT   str2ofmt( const string& s );
   static string ofmt2str( const OFMT&   f );

   void check();                               // settings check
   void help();                                // settings help

   // attributes
   int verb();                                 // get verbosity attribute
   bool append();                              // get append request

   virtual void _upd_glog();                   // update glog (mask,verb)
     
   // elements
   int output_size(const string& fmt);         // get format output size
   string outputs(const string& fmt);          // get string outputs
   set<string> oformats();                     // get formats

   string version(const string& fmt);          // get string output version
   
   int             file_toff(const string& fmt);  // get time offset [min] for output filename
   t_gtime::t_tsys file_tsys(const string& fmt);  // get time system       for output filename

   int   out_update(const string& fmt);        // get update [min] for output file update
   int   out_length(const string& fmt);        // get length [min] for output file content
   float out_sample(const string& fmt);        // get sample [sec] for output file data

 protected:
   string _outputs(const string& fmt);
   set<string> _oformats();

   set<OFMT> _OFMT_supported;                  // vector of supported OFMTs (app-specific)
   bool  _append;                              // append mode
   int   _verb;                                // output verbosity   
   int   _upd;                                 // update [min] for output file update
   int   _len;                                 // length [min] for output file content
   float _smp;                                 // sample [sec] for output file data
   
 private:

};

} // namespace

#endif
