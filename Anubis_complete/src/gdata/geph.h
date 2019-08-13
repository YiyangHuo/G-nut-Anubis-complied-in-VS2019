
#ifndef GEPH_H
#define GEPH_H 
 
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

#include <memory>

#include "gio/gio.h"
#include "gdata/gdata.h"
#include "gutils/gsys.h"
#include "gutils/gtime.h"

using namespace std;
    
namespace gnut {

enum NAVDATA {
        NAV_UNDEF,
        NAV_A, NAV_E, NAV_M, NAV_I, NAV_IDOT, NAV_OMEGA, NAV_OMG, NAV_OMGDOT, NAV_DN,
        NAV_CRC, NAV_CIC, NAV_CUC, NAV_CRS, NAV_CIS, NAV_CUS, NAV_F0, NAV_F1, NAV_F2,
        NAV_X, NAV_XD, NAV_XDD, NAV_Y, NAV_YD, NAV_YDD, NAV_Z, NAV_ZD, NAV_ZDD,
        NAV_IOD, NAV_HEALTH,
        NAV_TGD0, NAV_TGD1, NAV_TGD2, NAV_TGD3
      };
   
typedef pair<t_gtime, double>  t_timdbl;
   

class t_geph : public t_gdata {

 public:
  t_geph();
  virtual ~t_geph();
   
  // pointers to support NULL if not requested!
  virtual int clk( const t_gtime& t, double*    clk, double*   var = NULL, double*  dclk = NULL, bool chk_health = true ){ return -1; }  // [s]
  virtual int pos( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true ){ return -1; }  // [m]
  virtual int nav( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true ){ 
              return this->pos(t, xyz, var, vel, chk_health); }  // [m]
   
  virtual bool healthy() const { return true; }
  virtual string health_str() const { return ""; }
  virtual int chk() const { return -1; }

  virtual GNAVTYPE gnavtype(bool full = true) const { return NAV; }
  virtual int src(bool full = true) const { return  0; }
   
  virtual string linefmt() const { return ""; }
  virtual string line() const { return ""; }
  virtual void print() const { cout << linefmt(); }

  virtual t_timdbl param( const NAVDATA& n );
  virtual int      param( const NAVDATA& n, double val );
  virtual bool     param_cyclic( const NAVDATA& n );

  virtual void gio(shared_ptr<t_gio> p){ _gio_ptr = p; }
  shared_ptr<t_gio> gio(){ return _gio_ptr; }

  void clear();
  bool valid();
  void valid(bool validity) {_validity = validity;}   

  GSYS    gsys() const;                                          // GNSS system
  string  gsat() const;                                          // satellite number

  // POZDEJI JEN GSAT a vse pres MUTEX !!!
  string  sat()      const{ return _sat; }                       // satellite number  
  double  interval() const{ return _interval; }                  // get validity interval
  t_gtime epoch()    const{ return _epoch; }                     // reference epoch
  t_gtime begin()    const{ return _epoch - _interval/2; }       // beg of validity
  t_gtime end()      const{ return _epoch + _interval/2; }       // end of validity

  virtual bool chktot(const t_gtime& t) { return true; }
   
 protected:
  virtual void _clear();
  virtual bool _valid() const;

  string            _sat;         // satellite number
  t_gtime           _epoch;       // reference epoch
  double            _interval;    // validity interval
  bool              _validity;       // validity
   
  shared_ptr<t_gio> _gio_ptr;

 private:

};

} // namespace

#endif
