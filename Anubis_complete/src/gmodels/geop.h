
#ifndef GEOP_H
#define GEOP_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: earth orientation parameters, base class
           Neither precession nor nutation apllied -> return identity matrix
  Version: $ Rev: $

  2013-03-27 /PV: created
  
-*/

#include "gdata/gdata.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"

using namespace std;

namespace gnut {   

// ----------
class t_geop {

 public:
  t_geop();
  virtual ~t_geop();
 
  Matrix nutMatrix(double mjd);
  Matrix precMatrix(double mjd);
   
// protected:
   double _normangle(double x);   
};

} // namespace

#endif
