
#ifndef GSETNWM_H
#define GSETNWM_H

#define XMLKEY_NWM "nwm"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements data extraction setting
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gnwmsurf.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut {

class t_gsetnwm : public virtual t_gsetbase
{
 public:
   t_gsetnwm();
  ~t_gsetnwm();
  
   void check();                                 // settings check
   void help();                                  // settings help

   double   min_rad();
   double   min_lat();
   double   max_lat();
   double   min_lon();
   double   max_lon();
   double   grid_mult();
   double   decay_min();
   double   decay_max();
   int      vert_mult();
   NWMFIT   vert_fit();
   NWMADJ   vert_adj();
   PROFTEMP vert_temp();
   PROFZWD  vert_zwd();
   SURFZWD  surf_zwd();
   string   refr_coef();
   string   interp_time();
   string   interp_vert();
   string   interp_space();
   string   interp_plane();

   set<string> param();                          // get param request
   set<string> point();                          // get point request
   set<string> surface();                        // get surface request
   set<string> profile();                        // get profile request
   set<string> assess();                         // get assess request
   set<string> clean();                          // get clean request

 protected:

   set<string> _param;
   set<string> _point;
   set<string> _surface;
   set<string> _profile;
   set<string> _assess;

   double      _min_rad;                        // radius factor     (select nearest grid points)
   double      _min_lat;	                // minimum latitude  (region)
   double      _max_lat;	                // maximum latitude  (region)
   double      _min_lon;	                // minimum longitude (region)
   double      _max_lon;	                // maximum longitude (region)
   double      _grid_mult;	                // multiplicator for getting grid dLat/dLon around points
   double      _decay_min;                      // min limit for decay parameters
   double      _decay_max;                      // max limit for decay parameters
   string      _refr_name;	                // refractivity coefficients
   string      _refr_list;	                //
   int         _vert_mult;	                // multiplicator for increase of vertical resolution
   string      _vert_temp;	                // method of vertical T   representation
   string      _vert_zwd;	                // method of vertical ZWD representation
   string      _vert_adj;	                // method of vertical adjustement
   string      _vert_fit;	                // method of vertical approximation
   string      _surf_zwd;	                // method of vertical approximation
   string      _interp_time;	                // method of temporal approximation
   string      _interp_vert;	                // method of vertical approximation
   string      _interp_space;	                // method of spatial  approximation
   string      _interp_plane;	                // method of planar   approximation
   t_map_refr  _m_refr;

 private:
};

} // namespace

#endif
