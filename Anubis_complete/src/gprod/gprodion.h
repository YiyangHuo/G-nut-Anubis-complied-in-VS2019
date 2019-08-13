
#ifndef GPRODION_H
#define GPRODION_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2018-04-10 /PV: created

-*/

#include <iostream>

#include "gprod/gprod.h"

using namespace std;

namespace gnut {

class t_gprodion : public t_gprod {

 public:
  t_gprodion(const t_gtime& t, shared_ptr<t_gobj> pt = nullobj);
  virtual ~t_gprodion();
   
  typedef map<double, map<double, double> > t_map_grd;

  double const tec(const double& lat, const double& lon);                     // get TEC for particular location
  void         tec(const double& lat, const double& lon, const double& val);  // set TEC for a grid point

  void         clear();
 protected:
  t_map_grd _tec_map;

 private:

};

} // namespace

#endif
