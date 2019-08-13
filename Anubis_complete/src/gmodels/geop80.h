
#ifndef GEOP80_H
#define GEOP80_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: earth orientation parameters (IAU Conventions 1980)
  Version: $ Rev: $

  2013-03-27 /PV: created

 Reference: 
      Explanatory Supplement to the Astronomical Almanac,
      P. Kenneth Seidelmann (ed), University Science Books (1992),
      Section 3.222 (p111).
 
      Derived from SOFA software (www.iausofa.org)
  
-*/

#include "gdata/gdata.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gmodels/geop.h"

using namespace std;

namespace gnut {   

// ----------
class t_geop80 : public t_geop {

 public:
  t_geop80();
  virtual ~t_geop80();
 
  Matrix nutMatrix(double mjd);
  Matrix precMatrix(double mjd_1);
   
 protected:
   double _frac(double x);      
};

} // namespace

#endif
