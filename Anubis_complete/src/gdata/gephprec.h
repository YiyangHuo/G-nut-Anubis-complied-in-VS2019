
#ifndef GEPHPREC_H
#define GEPHPREC_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: derive eph to implement precise ephemerides class
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include "geph.h"
#include "gmodels/gpoly.h"

#define UNDEFVAL_CLK 999999999.999
#define UNDEFVAL_POS         0.000

#define MAXDIFF_CLK 300.0
#define MAXDIFF_EPH 900.0

using namespace std;

namespace gnut {

class t_gephprec : public t_geph {
  
 public:
  t_gephprec();
  virtual ~t_gephprec();

  // pointers to support NULL if not requested!
  int pos(     const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL, bool chk_health = true ); // [m]
  int clk(     const t_gtime& t, double*    clk, double*   var = NULL, double*  dclk = NULL, bool chk_health = true ); // [s]
   
  int pos_int( const t_gtime& t, double  xyz[3], double var[3] = NULL, double vel[3] = NULL ); // [m])
  int clk_int( const t_gtime& t, double*    clk, double*   var = NULL, double*  dclk = NULL ); // [s]
   
  int add(string sat, vector<t_gtime> t, 
	              const vector<double>& x, 
	              const vector<double>& y,
	              const vector<double>& z, 
	              const vector<double>& c);

  int chk()const{ return 1; }
  string str() const { return ""; }
  int print(){ cout << str(); return 0; }
   
  void degree( int d){ _clear(); _degree = d; }          // set/get degree of polynomials
  unsigned int degree()const{ return _degree; }
   
  unsigned int ndata()const{ return _degree+1; }         // get number of needed data
  unsigned int interval()const{ return (size_t)_poly_x.span(); } // get validity span
  unsigned int xref()const{ return (size_t)_poly_x.xref(); }     // get xreference value
  bool valid(const t_gtime& t) const;                    // check validity (incl.data span)
   
 protected:
  void _clear();
  bool _valid_clk() const;
  bool _valid_crd() const;   

 private:
  unsigned int      _degree;      // polynomial degree
  t_gpoly           _poly_x;      // X polynomials
  t_gpoly           _poly_y;      // Y polynomials
  t_gpoly           _poly_z;      // Z polynomials
  t_gpoly           _poly_c;      // C polynomials
 
  // _epoch is reference epoch
  vector<double>    _dt;          // time difference to _epoch
  vector<double>    _xcrd;        // x-coordinate
  vector<double>    _ycrd;        // y-coordinate
  vector<double>    _zcrd;        // z-coordinate
  vector<double>    _clkc;        // clock correction
};

} // namespace

#endif
