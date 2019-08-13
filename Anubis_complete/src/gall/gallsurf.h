
#ifndef GALLSURF_H
#define GALLSURF_H

// TO SUPPORT BOTH LAT/LON  Y/X coordinate systems
#define MIN_LAT -9999  //   -90
#define MAX_LAT +9999  //    90
#define MIN_LON -9999  //     0
#define MAX_LON +9999  //   360

#define DEFAULT_SAMPLING 30
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: container for meteo data
  Version: $ Rev: $

  2013-09-25 /JD: created

-*/

#include <map>
#include <set>
#include <string>
#include <memory>
#include <string>
#ifdef BMUTEX
  #include <boost/thread/mutex.hpp>
#endif

#include "gio/giof.h"
#include "gdata/gdata.h"
#include "gutils/gpair.h"
#include "gutils/gtime.h"
#include "gutils/gconst.h"
#include "gutils/gtriple.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gnwmsurf.h"
#include "gmodels/ginterp.h"
#include "gmodels/gfuncplane.h"
#include "gproj/gproj.h"
#include "gset/gsetbase.h"
#include "../newmat/newmat.h"

using namespace std;

namespace gnut {  

// ======================================================
//                interpolation       mesh_in_KM
//                 (grid_intv)        (grid_dist)
// ======================================================
//                 N/lat  E/lon     N/lat[km]   E/lon[km]
// ======================================================
//  meter          loc_N  loc_E     loc_N       loc_E
//  index          idx_N  idx_E     mesh_N      mesh_E
//  degree         deg_N  deg_E     sphere_N    sphere_E
// =======================================================
//
class t_gallsurf : public t_ginterp  // already inherit g_data
{
 public:
  t_gallsurf();
  t_gallsurf(t_gsetbase* gset);
  t_gallsurf(double minLat, double maxLat, double minLon, double maxLon);

  virtual ~t_gallsurf();
   
  typedef map<t_gpair, shared_ptr<t_gnwmsurf>> t_map_dat;// map of latitude/longitude data (single epo)
  typedef map<t_gtime, t_map_dat>              t_map_met;// map of all meteo data (all epochs)

  virtual int add(shared_ptr<t_gnwmsurf> met,bool proj=false);  // add single element

  virtual set<t_gtime> epochs();                         // get all NWM epochs
  virtual set<t_gpair> pairs(const t_gtime& epo);        // get all longitude/latitudes pairs (for epoch)

  void   grid_type(string s){ _grid_type = s; }          // grid type
  string grid_type(){ return _grid_type; }               // grid type
   
  void   grid_idx(bool ll){ _grid_ll = ll; }             // grid latlon index
  bool   grid_idx(){ return _grid_ll; }                  // grid latlon index

  void   grid_intv(double iLat, double iLon){ _grid_lat = fabs(iLat);             // grid interval [various units]
                                              _grid_lon = fabs(iLon); }
  void   grid_dist(double iLat, double iLon){ _dist_lat = fabs(iLat);             // grid distance [km]
                                              _dist_lon = fabs(iLon); }

  double vert_dist(){ return  _vert_dst; }                                        // [m] distance (maximum value)
  double grid_dist(){ return((_dist_lat>_dist_lon)?_dist_lat:_dist_lon)*1000.0;}  // [m] distance (maximum value)

  bool   reg_auto(){ return _reg_auto; }
  double min_lat(){ return _min_lat; }
  double max_lat(){ return _max_lat; }
  double min_lon(){ return _min_lon; }
  double max_lon(){ return _max_lon; }
   
  int                 gproj(const shared_ptr<t_gproj>& prj); // set projection (only if not set before!)
  shared_ptr<t_gproj> gproj(){ return _proj; }               // get projection pointer (null if not set)
   
  // meteo parameter interpolation
  int interpolate_point(const t_gtriple& ell,
                        const t_gtime&   beg,
	                const t_gtime&   end,
 	                map<t_gtime, t_map_metdata>          *m_data = 0,
 	 	        map<t_gtime, t_map_metsurf>          *m_surf = 0,
 	 	        map<t_gtime, shared_ptr<t_gnwmsurf>> *m_prof = 0,
		        double sampling = 0 );

   // analyse gradients
   int analyse_gradient( const t_gtriple& ell,                            // ELL:BLH of point of interest
			 const double& cut_off,                           // [deg] elevation cut-off
                         const double sampling,                           // [sec] sampling interval
			 map<t_gtime, map<MET_DATA, double> >& rf_data);  // ref.point interpolated data (interpolate_point)

   int add_param(MET_DATA par);
   int add_param(MET_SURF par);

//int grid_data(const& MET_DATA, t_gpair&, LatLon );     // print grid
//int grid_data(const& MET_SURF, t_gpair&, LatLon );     // print surf
     
//  virtual int bilinear(const t_gpair& LatLon, const t_gtime& epo); // spatial interpolation - bilinear
//  virtual int kriging(const t_gpair& LatLon, const t_gtime& epo);  // spatial interpolation - kriging
//  virtual int spline(const t_gpair& LatLon, const t_gtime& epo);   // spatial interpolation - spline

 protected:
   virtual int _new_element(shared_ptr<t_gnwmsurf> met);      // add element in container (in heap)

   virtual int _set_region(const double minLat,
	   	           const double maxLat, 
	   	           const double minLon,
		           const double maxLon);              // set region (should be done at the constructor only)

   virtual void _set_output();
   virtual void _set_param( const set<string>& param );
   virtual int _get_grid_step();                              // get Lat/Lon sampling step (equidistant!)
   virtual int _get_node_epochs(const t_gtime& epo,           // get node epochs (single epoch)
			              set<t_gtime>& set_dt);  
   virtual int _get_node_epochs(const t_gtime& beg,           // get node epochs (interval)
			        const t_gtime& end,
			              set<t_gtime>& set_dt);  
   virtual int _get_grid_points(const t_gpair& pt,            // get data points
			        const set<t_gtime>& set_dt,
			  	      set<t_gpair>& set_pt);
   virtual int _get_grid_interp(const t_gtriple& pt,          // get maps of GRID/REQ.PT for all NWM epochs
				const set<t_gtime>& set_epo,
				const set<t_gpair>& set_pt,
                                const vector<MET_SURF> vec_surf,
				t_map_met& map_grid,
				t_map_met& map_prof);
   virtual int _interp_temporal(const t_gtime& beg,
			        const t_gtime& end,
		                const double& sampling,
	  		        const map<t_gtime,double>& inp,
	 	                      map<t_gtime,double>& out);

   t_map_met    _map_met;                                     
   string       _ver;
   string       _grid_type;
   bool         _overwrite;
   bool         _reg_auto;
   double       _min_lat;
   double       _max_lat;
   double       _min_lon;
   double       _max_lon;
   
   bool         _grid_ll;
   double       _grid_lat;
   double       _grid_lon;
   double       _dist_lat;
   double       _dist_lon;
   double       _vert_dst;

   PROFTEMP     _vert_temp;
   PROFZWD      _vert_zwd;
   SURFZWD      _surf_zwd;
   NWMFIT       _vert_fit;
   NWMADJ       _vert_adj;

   vector<MET_DATA> _vec_par;
   vector<MET_SURF> _vec_srf;

   shared_ptr<t_gproj> _proj;

   t_gsetbase*  _set;
   t_giof*      _out;
   t_giof*      _sum;
   t_giof*      _grd;
   t_giof*      _srf;
   t_giof*      _fit;
};

} // namespace

#endif
