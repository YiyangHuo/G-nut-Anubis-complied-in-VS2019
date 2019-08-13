
#ifndef GNAVBDS_H
#define GNAVBDS_H 
 
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

#include <set>

#include "gdata/gnav.h"
#include "gutils/gtime.h"
#include "gutils/gconst.h"
#include "gutils/gcommon.h"
#include "gutils/gtriple.h"

#define MAX_RINEXN_REC_BDS  29

using namespace std;

namespace gnut {

// ICD BeiDou v2
#define SC2R        3.1415926535898 // semi-circle to radian

// RTCM3 BDS EPH decoding consistency with BNC
//////////////////////////////////////////////////////////
#define BDSTOINT(type, value) static_cast<type>(round(value))
#define BDSADDBITS(a, b) {bitbuffer = (bitbuffer<<(a)) \
                       |(BDSTOINT(long long,b)&((1ULL<<a)-1)); \
                       numbits += (a); \
                       while(numbits >= 8) { \
                       buffer[size++] = bitbuffer>>(numbits-8);numbits -= 8;}}
#define BDSADDBITSFLOAT(a,b,c) {long long i = BDSTOINT(long long,(b)/(c)); \
                             BDSADDBITS(a,i)};

class t_gnavbds : public t_gnav {

 public:
  t_gnavbds();
  virtual ~t_gnavbds();

  // pointers to support NULL if not requested!
  int pos( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true ); //[m]
  int clk( const t_gtime& t, double*    clk, double*   var = NULL, double*  dclk = NULL, bool chk_health = true ); //[s]
  int chk(set<string>& msg);
  int ura( double acc ) const;
  int data2nav( string, const t_gtime& ep, const t_gnavdata& data );   
  int nav2data( t_gnavdata& data );

  int iod() const;
  //int iod() const { return _iode; }
  int rec() const { return MAX_RINEXN_REC_BDS; }
  double tgd() { return (C02_F*C02_F) / (C02_F*C02_F - C07_F * C07_F)*_tgd[0] - (C07_F*C07_F) / (C02_F*C02_F - C07_F * C07_F)*_tgd[1]; }
  double rel() { return _rel;  }
  virtual bool chktot(const t_gtime& t);
   
  t_timdbl param( const NAVDATA& n );
  int      param( const NAVDATA& n, double val );
   
  string line() const;
  string linefmt() const;

 private:
  bool    _healthy() const;
  void    _ecc_anomaly( double dt, double& Ek, double& dEk );
  double  _mean_motion();
  const int	 _getIODC() const;

  int     _iode;        // issue of ephemeris
  int     _iodc;        // issue of clocks
  int     _health;      // sv health (0:ok)
  t_gtime _toe;         // time of ephemerides 
  t_gtime _toc;         // time of clocks
  t_gtime _tot;         // time of transmission
  double  _a;           // major semiaxis [m] 
  double  _e;           // eccentricity 
  double  _m;           // mean anomaly at t0 [rad] 
  double  _i;           // inclination [rad] 
  double  _idot;        // inclination change [rad/sec]
  double  _omega;       // argument of perigee [rad]
  double  _OMG;         // longit. of ascend. node at weekly epoch [rad]
  double  _OMGDOT;      // longit. of ascend. node at weekly epoch's change [rad/sec]
  double  _dn;          // mean motion difference, delta n [rad/sec] 
  double  _crc;         // amplit. of cos harm. terms - radius [m]
  double  _cic;         // amplit. of cos harm. terms - inclination [rad]
  double  _cuc;         // amplit. of cos harm. terms - arg-of-latitude [rad]
  double  _crs;         // amplit. of sin harm. terms - radius [m]
  double  _cis;         // amplit. of sin harm. terms - inclination [rad]
  double  _cus;         // amplit. of sin harm. terms - arg-of-latitude [rad]
  double  _f0;          // sv clock parameters [sec]
  double  _f1;          // sv clock parameters [sec/sec]
  double  _f2;          // sv clock parameters [sec/sec^2]
  double  _acc;         // sv accuracy index (URA) 
  double  _tgd[2];      // group delay parameters [sec]  0: B1/B3,  1: B2/B3
  double  _rel;			//relativity correction calculated with ICD
};

} // namespace

#endif
