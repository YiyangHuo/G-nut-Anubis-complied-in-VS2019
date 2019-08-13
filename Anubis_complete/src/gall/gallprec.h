
#ifndef GALLPREC_H
#define GALLPREC_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implementation of SP3 ephemeris container
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include "gall/gallnav.h"
#include "gdata/gephprec.h"
#include "gutils/gtriple.h"
#include "gmodels/gpoly.h"

using namespace std;

namespace gnut {  

typedef map<string,     double>       t_map_dat;     // single data for a single satellite
typedef map<t_gtime, t_map_dat>       t_map_epo;     // all data for a single satellite
typedef map<string,  t_map_epo>       t_map_prn;     // all data for all satellites

class t_gallprec : public t_gallnav {
  
 public:
  t_gallprec();
  virtual ~t_gallprec();

  typedef map<string, shared_ptr<t_gephprec> > t_map_sp3; 

  virtual bool health( string sat, const t_gtime& t );  // inherited from gallnav to fix healthy
   
  virtual int pos(    string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true ); // [m]
  virtual int nav(    string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true ); // [m] GNAV quality
  virtual int clk(    string sat, const t_gtime& t, double* clk,    double* var    = NULL, double* dclk   = NULL, bool chk_mask = true ); // [s]
  virtual int clk_sp3(string sat, const t_gtime& t, double* clk,    double* var    = NULL, double* dclk   = NULL); // [s]

  // ALTERNATIVES for direct interpolations - using gephprec!
  virtual int pos_int(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL); // [m]
  virtual int clk_int(string sat, const t_gtime& t, double* clk,    double* var    = NULL, double* dclk   = NULL); // [s]
  // ALTERNATIVES for direct interpolations !
  virtual int pos_alt(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL); // [m]
  virtual int clk_alt(string sat, const t_gtime& t, double* clk,    double* var    = NULL, double* dclk   = NULL); // [s]

//  int addpos( string sat, const t_gtime& ep, double xyzt[4], double var[4] );  // [m], [m]
  int addvel( string sat, const t_gtime& ep, double xyzt[4], double var[4] );  // [m], [m]
  int addclk( string sat, const t_gtime& ep, double  clk[3], double var[3] );  // [s]

  int addpos( string sat, const t_gtime& ep, t_gtriple xyz, double t, t_gtriple dxyz, double dt);  // [m], [usec]
//  int addvel( string sat, const t_gtime& ep, double xyzt[4], double var[4] );  // [m], [m]
//  int addclk( string sat, const t_gtime& ep, double  clk[3], double var[3] );  // [s]

  void use_clkrnx(bool b){ _clkrnx = b; }   
  void use_clksp3(bool b){ _clksp3 = b; }
  void use_clknav(bool b){ _clknav = b; }
  void use_posnav(bool b){ _posnav = b; }   
  void clean_outer( const t_gtime& beg = FIRST_TIME, const t_gtime& end = LAST_TIME );

  virtual set<string> satellites();                                      // get all satellites
  virtual unsigned int nepochs( const string& prn );                     // get number of epochs
  virtual t_gtime beg_data( string prn = "");                            // get first t_gephprec epoch
  virtual t_gtime end_data( string prn = "");                            // get  last t_gephprec epoch
  // IMPROVE beg_time/end_time to distinguish GALLNAV/GALLPREC - t_satview !
  virtual t_gtime beg_time( string prn = ""){ return beg_gnav(prn); }    // get first t_gephprec epoch
  virtual t_gtime end_time( string prn = ""){ return end_gnav(prn); }    // get first t_gephprec epoch
  // =========================================================
  virtual t_gtime beg_clk(  string prn = "");                            // get first precise clocks epoch
  virtual t_gtime end_clk(  string prn = "");                            // get last precise clocks epoch
  virtual t_gtime beg_prec( string prn = "");                            // get first precise polynomials epoch
  virtual t_gtime end_prec( string prn = "");                            // get last precise polynomials epoch

 protected:      
  virtual shared_ptr<t_geph> _find( string sat, const t_gtime& t );  // find appropriate t_geph element
  virtual int         _get_crddata( string sat, const t_gtime& t );  // fill PT,X,Y,Z vectors
  virtual int         _get_clkdata( string sat, const t_gtime& t );  // fill CT,C vectors
   
  // also gallnav member can be used here !!!
  t_map_sat         _mapprec;     // map of sp3 polynomials
  t_map_prn         _mapsp3;      // precise orbits&clocks (SP3) - full discrete data sets
  t_map_prn         _mapclk;      // precise clocks (CLOCK-RINEX) - full discrete data sets
   
 private:
  t_map_sp3         _prec;        // CACHE: single SP3 precise ephemeris for all satellites
  unsigned int      _degree_sp3;  // polynom degree for satellite sp3 position and clocks
  double            _sec;         // default polynomial units
  t_gtime           _ref;         // selected reference epoch for crd data/polynomials
  t_gtime           _clkref;      // selected reference epoch for clk data/polynomials
  bool              _clkrnx;      // true: use               clk from Rinex Clocks
  bool              _clksp3;      // true: use alternatively clk from sp3 (~15min!)
  bool              _clknav;      // true: use alternatively nav (low-precise clocks)
  bool              _posnav;      // true: use alternatively nav (low-precise orbits)

// CACHE for approximative solutions
  map<string,t_gtime>   _poly_beg;
  map<string,t_gtime>   _poly_end;
  map<string,t_gpoly>   _poly_x;
  map<string,t_gpoly>   _poly_y;
  map<string,t_gpoly>   _poly_z;

// BEGIN OF TEMPORARY (ALTERNATIVE for direct interpolation)
  vector<double>    _PT;          // vector of time-difference (X for polynomials)
  vector<t_gtime>   _T;           // vector of full time       (X for polynomials)
  vector<double>    _X;           // vector of x-coordinate    (Y for polynomials)
  vector<double>    _Y;           // vector of y-coordinate    (Y for polynomials)
  vector<double>    _Z;           // vector of z-coordinate    (Y for polynomials)

  vector<double>    _CT;          // vector of time-difference (X for polynomials)
  vector<double>    _C;           // vector of clk correction  (Y for polynomials)
// END OF TEMPORARY (ALTERNATIVE)

};

} // namespace

#endif
