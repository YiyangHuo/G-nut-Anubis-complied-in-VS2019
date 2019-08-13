
#ifndef GINTERP_H
#define GINTERP_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implementation various interpolation tecnique
  Version: $ Rev: $

  2013-06-05 /PV: created

-*/

#include <vector>
#include <map>

#include "../newmat/newmat.h"

#include "gutils/gtime.h"
#include "gutils/gpair.h"
#include "gdata/gdata.h"

using namespace std;

namespace gnut {   

// ----------
class t_ginterp : public t_gdata 
{
 public:   
            t_ginterp();
   virtual ~t_ginterp();

   enum INTERP_1D { LINEAR, SPLINE };
   enum INTERP_2D { BILINEAR, IDW };
   enum INTERP_3D { VER2HOR, HOR2VER };
   enum INTERP_HT { INTERPOLATE, SCALE };

   INTERP_1D str_to_interp_1d(const string& str);
   INTERP_2D str_to_interp_2d(const string& str);
   INTERP_3D str_to_interp_3d(const string& str);
   INTERP_HT str_to_interp_ht(const string& str);

   string    interp_1d_to_str(const INTERP_1D& typ);
   string    interp_2d_to_str(const INTERP_2D& typ);
   string    interp_3d_to_str(const INTERP_3D& typ);
   string    interp_ht_to_str(const INTERP_HT& typ);
   
   int linear(        map<double, double>& data, double val, double& fval);
   int spline(        map<double, double>& data, double val, double& fval);
   int linear(  const map<t_gtime, double>& data, const t_gtime& epo, double& fval);
   int spline(  const map<t_gtime, double>& data, const t_gtime& epo, double& fval);
   int bilinear(const map<t_gpair, double>& data, const t_gpair& req_pos, double& fval);
   int idw(     const map<t_gpair, double>& data, const t_gpair& req_pos, double& fval);   
   
 protected:
   INTERP_1D  _interp_1d;
   INTERP_2D  _interp_2d;
   INTERP_3D  _interp_3d;
   INTERP_HT  _interp_ht;

};

} // namespace

#endif
