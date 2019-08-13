
#ifndef GBLINDMOPS_H
#define GBLINDMOPS_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: subroutine focused on RTCA MOPS blind model
  Version: $ Rev: $

  2014-03-21 /ME: created
 
 
  Reference:

    RTCA (2006), Minimum operational performance standards for Global Positioning 
                 System/Wide Area Augmentation System airborne equipment,
    Radio Technical Commission for Aeronautics, SC-159, publication DO-229D, Washington, D. C.

-*/

#include <map>
#include <string>
#include <set>
#include <cmath>

#include "../newmat/newmat.h"

#include "gmodels/ginterp.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {   

// ----------
class t_mops { 

 public:
   t_mops(){};
   virtual ~t_mops(){};
   
   void mops(  double& lat, double& H, double& elev, double& DOY, double& ZHD, double& ZWD);
   int delay(  double& H, double& elev, vector<double>& EVAL, double& ZHD0, double& ZHD, double& ZWD0, double& ZWD, double& mf );
     
 protected:
   int _tabval( double& lat, Matrix AVG, Matrix VAR, vector<double>& I_AVG, vector<double>& I_VAR );
   int _extrap( double& DOY, double& DOYmin, vector<double>& AVG, vector<double>& VAR, vector<double>& EVAL );
 
 private:
   static double a[30];
   static double b[30];

   
};

} // namespace

#endif