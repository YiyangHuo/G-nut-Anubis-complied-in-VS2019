
#ifndef GNAVSBS_H
#define GNAVSBS_H 
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $ Rev: $

  2011-02-14 /JD: created
  2018-08-13 /JD: updated

-*/

#include <vector>

#include "../newmat/newmat.h"
#include "gdata/gnav.h"
#include "gutils/gtime.h"
#include "gutils/gconst.h"
#include "gutils/gcommon.h"
#include "gutils/gtriple.h"

#define MAX_RINEXN_REC_SBS 15

//#define MIN_SBS_RADIUS 25400       // km
//#define MAX_SBS_RADIUS 25600       // km
//#define MAX_SBS_RADDIF 25          // km
//#define MAX_SBS_CLKDIF 15          // ns
//#define RAD_SBS_FACTOR 0.001       // m->km
//#define CLK_SBS_FACTOR 1000000000  // sec->ns

using namespace std;

namespace gnut {

/* ura values (ref [3] 20.3.3.3.1.3) [meters] */
//static const double ura_eph[]=
//  {2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,3072.0,6144.0,0.0};

class t_gnavsbs : public t_gnav {

 public:
  t_gnavsbs();
  virtual ~t_gnavsbs();

  // pointers to support NULL if not requested!
  virtual int pos( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true );
  virtual int clk( const t_gtime& t, double*    clk, double*   var = NULL, double*  dclk = NULL, bool chk_health = true );
  int chk(set<string>& msg);
  int ura( double acc ) const;
  int data2nav( string sat, const t_gtime& ep, const t_gnavdata& data );   
  int nav2data( t_gnavdata& data );

  int iod() const { return _iod; }
  int rec() const { return MAX_RINEXN_REC_SBS; }
   
  t_timdbl param( const NAVDATA& n );
  int      param( const NAVDATA& n, double val );

  string line() const;
  string linefmt() const;
   
 protected:
   bool _healthy() const;

 private:   
   
   double   _x;       // position X [km]
   double   _x_d;     // velocity X [km/s]
   double   _x_dd;    // acceleration X [km/s^2]
   double   _y;       // position Y [km]
   double   _y_d;     // velocity Y [km/s]
   double   _y_dd;    // acceleration Y [km/s^2]
   double   _z;       // position Z [km]
   double   _z_d;     // velocity Z [km/s]
   double   _z_dd;    // acceleration Z [km/s^2]
   double   _f0;      // SV clock bias [s]
   double   _tki;     // Transmission time of message in GPS seconds of weak
   double   _health;  // health 0 = OK
   double   _C_rms;   // Accuracy code [m]
   int      _iod;     // Issue of Data Navigation
   double   _f1;      // SV relative frequency bias []
   t_gtime  _toc;     // Epoch of ephemerides   
   
};

} // namespace

#endif
