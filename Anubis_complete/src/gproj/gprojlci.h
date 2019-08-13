
#ifndef GPROJLCI_H
#define GPROJLCI_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: Spherical to Lambert Conformal Conic projection (index-based)
  Version: $ Rev: $

  2015-02-19 /JD: created

-*/

#include <string>
#include <iomanip>
#include <iostream>

#include "gutils/gpair.h"
#include "gproj/gproj.h"

using namespace std;

namespace gnut {
 
class t_gprojlci : public t_gproj {
	
 public:
  t_gprojlci( double lon0, double lat0, double mesh,
	      double xoff, double yoff, double radi );
  virtual ~t_gprojlci(){};

  virtual t_gpair ll2proj(const t_gpair& latlon, int dec = 0); // from spherical LL[deg] to LCI XY[idx]  (idx ge 0!)
  virtual t_gpair proj2ll(const t_gpair& xy_idx);              // from LCI XY[idx] to spherical LL[deg]  (idx ge 0!)

  bool operator==(const t_gprojlci& prj) const;
   
  double lon0()const{ return _lon0; }
  double lat0()const{ return _lat0; }
  double mesh()const{ return _mesh; }
  double radi()const{ return _radi; }
  double xoff()const{ return _xoff; }
  double yoff()const{ return _yoff; }

 protected:

  double  _lon0;
  double  _lat0;
  double  _mesh;
  double  _radi;
  int     _xoff;
  int     _yoff;
   
  double  _K;
  double  _R0;
};

} // namespace

#endif
