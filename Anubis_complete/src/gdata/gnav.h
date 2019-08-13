
#ifndef GNAV_H
#define GNAV_H 
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements ephemerides and navigation classes
  Version: $ Rev: $

  2011-01-10 /JD: created

  TODO: overwritting mode
        return first/last
        limit for number of data inclusions per satellites
-*/

#include "gdata/geph.h"
#include "gutils/gsys.h"
#include "gutils/gtime.h"

#define MAX_RINEXN_REC     29   // maximum number of RINEXN records for any system !!

#define MAX_NAV_TIMEDIFF  3600*2    // NAV GNS valitity interval [s]
#define MAX_GPS_TIMEDIFF  3600*2    // NAV GPS validity interval [s]
#define MAX_GLO_TIMEDIFF    60*17   // NAV GLO validity interval [s]
#define MAX_GAL_TIMEDIFF  3600*3    // NAV GAL validity interval [s]
#define MAX_BDS_TIMEDIFF  3600      // NAV BDS validity interval [s]
#define MAX_SBS_TIMEDIFF   360      // NAV SBS validity interval [s]
#define MAX_IRN_TIMEDIFF  3600*2    // NAV IRN validity interval [s]
#define MAX_QZS_TIMEDIFF  3600      // NAV QZS validity interval [s]

using namespace std;

namespace gnut {

typedef double t_gnavdata[MAX_RINEXN_REC];

class t_gnav : public t_geph {

 public:
  t_gnav();
  virtual ~t_gnav();
  
  static  int nav_validity( GSYS gs );  // get GNSS NAV validity (half-interval) [s]

  virtual int data2nav( string sat, const t_gtime& ep, const t_gnavdata& data ){ return -1; }
  virtual int nav2data(                                      t_gnavdata& data ){ return -1; }

  virtual int iod() const { return -1; }
  virtual int rec() const { return MAX_RINEXN_REC; }

//  virtual int pos( t_gtime t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL ){ return 1; }
//  virtual int clk( t_gtime t, double*    clk, double*   var = NULL, double*  dclk = NULL ){ return 1; }
  virtual bool healthy() const;
  virtual string health_str() const { return _health_str(); }
  virtual int chk(set<string>& msg){ return 1; }

 protected:
  virtual bool _healthy() const { return true; }
  virtual string _health_str() const { return ""; }

 private:

};

} // namespace

#endif
