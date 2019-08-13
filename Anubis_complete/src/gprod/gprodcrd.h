
#ifndef GPRODCRD_H
#define GPRODCRD_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2011-03-25 /JD: created

-*/

#include <iostream>

#include "gprod/gprod.h"
#include "gutils/gtriple.h"

using namespace std;

namespace gnut {

 enum COV_TYPE { COV_XY, COV_XZ, COV_YZ };
       
class t_gprodcrd : public t_gprod {

 public:
  t_gprodcrd(const t_gtime& t, shared_ptr<t_gobj> pt = nullobj);
  virtual ~t_gprodcrd();
   
  void       xyz(const t_gtriple& xyz);
  t_gtriple  xyz() const;
   
  void       xyz_rms(const t_gtriple& xyz_rms);
  t_gtriple  xyz_rms() const;
  t_gtriple  xyz_var() const;

  void       apr(const t_gtriple& apr);
  t_gtriple  apr() const;   

  void       apr_rms(const t_gtriple& apr_rms);
  t_gtriple  apr_rms() const;
  t_gtriple  apr_var() const;
   
  void       cov(COV_TYPE type, double& cov);
  double     cov(COV_TYPE type) const;
   
 protected:
  t_gtriple   _xyz;
  t_gtriple   _xyz_rms;   
   
  t_gtriple   _apr;   
  t_gtriple   _apr_rms;

  double      _xy_cov;
  double      _xz_cov;
  double      _yz_cov;

 private:

};

} // namespace

#endif
