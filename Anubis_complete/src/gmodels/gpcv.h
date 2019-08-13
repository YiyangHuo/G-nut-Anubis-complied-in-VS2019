
#ifndef GPCV_H
#define GPCV_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements antenna phase offsets and variations
  Version: $ Rev: $

  2011-12-05 /JD: created

-*/

#include <vector>

#include "../newmat/newmat.h"
#include "gdata/gdata.h"
#include "gutils/gconst.h"
#include "gutils/gsys.h"
#include "gutils/gtime.h"
#include "gmodels/gephplan.h"
#include "gdata/gsatdata.h"

using namespace std;

namespace gnut {

// ----------
class t_gpcv : public t_gdata {

 public:
  t_gpcv();
  virtual ~t_gpcv();

  typedef map<double, double>    t_map_Z;      // map of ZEN-dependence
  typedef map<double, t_map_Z>   t_map_A;      // map of AZI-dependence & ZEN-dependence
   
  typedef map<GFRQ, t_gtriple>   t_map_pco;    // North/East/Up phase centre offsets     (N-freqs)
  typedef map<GFRQ, t_map_Z  >   t_map_zen;    // map of AZI-dependence                  (N-freqs)
  typedef map<GFRQ, t_map_A  >   t_map_azi;    // map of AZI-dependence & ZEN-dependence (N-freqs)

  bool   is_transmitter()const{return _trans;} // model is for transmitter (satellite)
  bool   is_receiver()const{  return !_trans;} // model is for receiver (station or LEO satellite)

  void   anten(string s){ _anten = s; }        // set/get antenna type definition (ATX field 0-19)
  string anten()const{ return _anten; }
  void   ident(string s){ _ident = s; }        // set/get antenna identification (ATX field 20-39)
  string ident()const{ return _ident; }
   
  void   svcod(string s){ _svcod = s;          // set/get satellite (svn) code (ATX field 40-59)
               if( s.empty() ) _trans = false;
                          else _trans = true;} // _trans=false/true (receiver/transmitter)
  string svcod()const{ return _svcod; }

  string pcvkey()const{ if(is_transmitter())   // pcvall first map key (antenna/satellite prn)
                             return ident();
                        else return anten(); }

  string pcvtyp()const{ if(is_transmitter())   // pcvall second map key (individual/type callibration)
                             return ""; 
                        else return ident(); }

  void   method(string s){ _method = s; }      // set/get method of callibration
  string method()const{ return _method; }
  void   source(string s){ _source = s; }      // set/get source of callibration
  string source()const{ return _source; }
  void   snxcod(string s){ _snxcod = s; }      // set/get sinex code
  string snxcod()const{ return _snxcod; }

  void    beg(t_gtime t){ _beg = t; }          // set/get valid from
  t_gtime beg()const{ return _beg;  }
  void    end(t_gtime t){ _end = t; }          // set/get valid until
  t_gtime end()const{ return _end;  }

  void   dazi(double d){ _dazi = d; }          // set/get azimuth sampling
  double dazi()const{ return _dazi; }
  void   dzen(double d){ _dzen = d; }          // set/get zenith sampling
  double dzen()const{ return _dzen; }
  void   zen1(double d){ _zen1 = d; }          // set/get zenith1
  double zen1()const{ return _zen1; }
  void   zen2(double d){ _zen2 = d; }          // set/get zenith2
  double zen2()const{ return _zen2; }

  // correct satellite coord
  double  pco( double zen, double azi, GFRQ f ); // pco correction - !!!! TEMP !!!!

  int     pcoS( t_gsatdata& satdata, t_gtriple& pco, GOBS_LC lc );                                // pco correction - satellite   
  int     pcoR( t_gsatdata& satdata, t_gtriple& dx, t_gtriple& site, GOBS_LC lc);                 // pco correction - receiver
  int     pco_proj(double& corr, t_gsatdata& satdata, t_gtriple& site, t_gtriple& dx, GOBS_LC lc);// pco projection into LoS
  int     pcvS(    double& corr, t_gsatdata& sat, t_gtriple& site, GOBS_LC lc );                  // pcb correction - satellite (nadir)
  int     pcvR(    double& corr, t_gsatdata& sat, GOBS_LC lc );                                   // pcv correction - site (zenith)
   
  void      pco(GFRQ f,   const t_gtriple& t ){ _mappco[f] = t; }  // set/get PCO
  t_gtriple pco(GFRQ f) { return _mappco[f]; }      // CANNOT BE CONST ?
  t_map_pco pco() const { return _mappco; }
   
  void      pcvzen(GFRQ f,const t_map_Z& t){ _mapzen[f] = t; }  // set/get PCO
  t_map_Z   pcvzen(GFRQ f){           return _mapzen[f]; }      // CANNOT BE CONST ?
  t_map_zen pcvzen() const{           return _mapzen; }
   
  void      pcvazi(GFRQ f,const t_map_A& t){ _mapazi[f] = t; }  // set/get PCO
  t_map_A   pcvazi(GFRQ f){           return _mapazi[f]; }      // CANNOT BE CONST ?
  t_map_azi pcvazi() const{           return _mapazi; }

 private:
  bool   _azi_dependent(GFRQ f);                    // does the calibration contain azi-depenedant data?
   
  t_gephplan     _ephplan;
   
  bool           _trans;      // transmitter[true], receiver[false]
  string         _anten;      // antenna type
  string         _ident;      // antenna identification
  string         _svcod;      // svn code
   
  string         _method;     // method of calibration
  string         _source;     // source of calibration
  string         _snxcod;     // sinex code
  t_gtime        _beg;        // valid from
  t_gtime        _end;        // valid until
  double         _dazi;       // azimuth sampling
  double         _dzen;       // zenith sampling
  double         _zen1;       // zenith1
  double         _zen2;       // zenith2

  t_map_pco      _mappco;     // map of PCOs (all frequencies)
  t_map_zen      _mapzen;     // map of NOAZI values (all frequencies)
  t_map_azi      _mapazi;     // map of AZI-dep values (all frequencies)

};

} // namespace

#endif
