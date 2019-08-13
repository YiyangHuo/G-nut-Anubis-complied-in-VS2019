
#ifndef GSYSCONV_H
#define GSYSCONV_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: reference system conversion utilities
  Version: $ Rev: $

  2011-04-26 /PV: created
  2012-04-06 /JD: extracted coordinate system conversion utilities only

-*/

#include <string>

#include "gutils/gtriple.h"
#include "gdata/gsatdata.h"

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

using namespace std;

namespace gnut {

int ell2xyz(const double* Ell, double* XYZ, bool degrees);           // True/False:  LatLon in Degrees/Radians!
int xyz2ell(const double* XYZ, double* Ell, bool degrees);           // True/False:  LatLon in Degrees/Radians!
int xyz2neu(const double* XYZ, const double* XYZ_Ref, double* neu);
int neu2xyz(const double* Ell, const double* neu,     double* xyz);  // ! Radians only: Ell[0], Ell[1]

int rao2xyz_rot(const ColumnVector& pos, const ColumnVector& vel, Matrix& R);
   
int rao2xyz(const ColumnVector& pos, const ColumnVector& vel,
            const ColumnVector& rao,       ColumnVector& xyz);
int xyz2rao(const ColumnVector& pos, const ColumnVector& vel,
                  ColumnVector& xyz,       ColumnVector& rao);   
   
int  ell2xyz(const t_gtriple &ell, t_gtriple &xyz, bool degrees);
int  xyz2ell(const t_gtriple &crd, t_gtriple &ell, bool degrees);    // True/False:  LatLon in Degrees/Radians!
int  xyz2ell2(t_gtriple &crd, t_gtriple &ell, bool degrees);         // True/False:  LatLon in Degrees/Radians!
void xyz2neu(t_gtriple &ell, t_gtriple &xyz, t_gtriple &neu);        // Radians only: Ell[0], Ell[1]
void neu2xyz(t_gtriple &ell, t_gtriple &neu, t_gtriple &xyz);        // Radians only: Ell[0], Ell[1]
   
int xyz2neu(t_gtriple& xyz, SymmetricMatrix& Q_xyz, SymmetricMatrix& Q_neu);

int ell2ipp(t_gsatdata& satdata, t_gtriple& ell_site, t_gtriple& ell_ipp);
   
} // namespace
   
#endif // GSYSCONV_H
