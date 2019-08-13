
#ifndef GPRODCLK_H
#define GPRODCLK_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2017-11-10 /JD: created

-*/

#include <iostream>

#include "gprod/gprod.h"
#include "gutils/gtriple.h"

using namespace std;

namespace gnut {

class t_gprodclk : public t_gprod {

 public:
  t_gprodclk(const t_gtime& t, shared_ptr<t_gobj> pt = nullobj);
  virtual ~t_gprodclk();
   
  void   clk( const double& val, const double& rms = 0.0);
  double clk();
  double clk_rms();

  void   icb( const double& val, const double& rms = 0.0);
  double icb();
  double icb_rms();

 protected:
  double      _clk;
  double      _clk_rms;   
  double      _icb;
  double      _icb_rms;

 private:

};

} // namespace

#endif
