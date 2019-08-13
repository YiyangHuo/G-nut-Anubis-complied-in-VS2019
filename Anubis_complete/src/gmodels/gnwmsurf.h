
#ifndef GNWMSURF_H
#define GNWMSURF_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implementation of NWM surface data
  Version: $ Rev: $

  2013-06-25 /JD: created

-*/

#include <vector>
#include <string>
#include <map>
#include <set>

#include "gmodels/gnwmbase.h"
#include "gutils/gtime.h"
#include "gio/glog.h"
#include "gset/gsetnwm.h"
#include "gset/gsetout.h"

#define RATIO_LIM 0.25      // limit for fabs(dE-dW) difference for estimating ratio
#define RATIO_EW  64.3      // ratio offset (roughly at surface)
#define RATIO_DH  0.066     // ratio vertical approximation [per 1m]
#define RATIO_XH  900.0     // ratio vertical approximation maximum height (above surface)
#define DECAY_MIN   0.0     // minimum value for vertical decay parameters
#define DECAY_MAX  50.0     // minimum value for vertical decay parameters

using namespace std;

namespace gnut {   

typedef map<MET_SURF, double> t_map_metsurf;
typedef map<MET_DATA, double> t_map_metdata;
     
// ----------
class t_gnwmsurf : public t_gnwmbase {

 public:
   t_gnwmsurf();
   t_gnwmsurf(set<MET_DATA> data, set<MET_SURF> surf);
   virtual ~t_gnwmsurf();
     
   virtual void   epoch(const t_gtime& t){_epoch = t; }   // set/get epoch
   virtual t_gtime epoch()const{       return _epoch; }

   virtual void glog(t_glog* l){            _log = l; }   // set/get glog pointer
   virtual t_glog* glog()const{          return _log; }

   virtual void decay_min(double d){  _decay_min = d; }   // min value allow for decay parameters
   virtual void decay_max(double d){  _decay_max = d; }   // max value allow for decay parameters
   virtual void vert_temp(PROFTEMP s){_vert_temp = s; }   // T vertical scaling model
   virtual void vert_zwd(PROFZWD s){  _vert_zwd  = s; }   // ZWD profile scaling model
   virtual void surf_zwd(SURFZWD s){  _surf_zwd  = s; }   // ZWD surface model

   virtual double   decay_min()const{ return _decay_min;} // min value allow for decay parameters
   virtual double   decay_max()const{ return _decay_max;} // max value allow for decay parameters
   virtual PROFTEMP vert_temp()const{ return _vert_temp;} // T vertical scaling model
   virtual PROFZWD  vert_zwd()const{  return  _vert_zwd;} // ZWD profile scaling model
   virtual SURFZWD  surf_zwd()const{  return  _surf_zwd;} // ZWD surface model

   virtual int    add_surf(MET_SURF typ, double val, bool reset = false);// set value for selected surf parameter
   virtual int    add_surf(MET_DATA typ, double val, bool reset = false);// set value for selected surf parameter
   virtual int    add_surfdata(MET_DATA typ, double val); // set value for DATA, but surface (e.g. TEMP = TSURF, ..)
   virtual double get_surf(MET_SURF typ) const;           // get value for selected surf parameter
   virtual double get_surf(MET_DATA typ) const;           // get value for selected surf parameter
   virtual double get_surfdata(MET_DATA typ) const;       // get value for DATA, but surface (e.g. TEMP = TSURF, ..)

   virtual double get_in_hgp(MET_DATA typ, double hgp);   // get value converted to geopotential height [m]
   virtual double get_in_hms(MET_DATA typ, double hms);   // get value converted to mean sea level height [m]
   virtual double get_in_hel(MET_DATA typ, double hel);   // get value converted to elipsoidal height [m]
   virtual double get_interp(MET_DATA typ, double hel);   // get value intepolate in profile (not applicable here)

   virtual int    add_rms(MET_SURF typ, double val);      // set value rms
   virtual int    add_rms(MET_DATA typ, double val);      // set value rms
   virtual double get_rms(MET_SURF typ) const;            // get value rms 
   virtual double get_rms(MET_DATA typ) const;            // get value rms

   virtual int    add_itr(MET_SURF typ, int val);         // set iterations
   virtual double get_itr(MET_SURF typ) const;            // get iterations

   virtual double lat()const{return get_surf(LAT);}       // get latitude
   virtual double lon()const{return get_surf(LON);}       // get longitude
   virtual double hgt()const{return get_surf(HEL);}       // get geometric height (elipsoidal)
   
   double vert_scaling(MET_DATA type,                     // vertical scaling function
                         double geopX,                    // required geopotential height
                         double geop,                     // reference geopotential height
                         double data,                     // reference data at geop
                         double temp = NWM_UNKNOWN,       // reference temp at geop (optional)
                         double pres = NWM_UNKNOWN );     // zwd vertical dependence

   int    calc_surface();                                 // calculate surface heights,gravity,geoid
   double radius_wgs84(const double& lat) const;          // WGS84 Earth radius(lat), Mahoney (2001) 
   double gravity_wgs84(const double& lat) const;         // WGS84 gravity(lat), Mahoney (2001)
   
   double geop2geom(const double& lat,                    // geop [gpm] -> geom [m], Mahoney (2001)
		    const double& geom) const;
   double geom2geop(const double& lat,                    // geom [m] -> geop [gpm], Mahoney (2001)
		    const double& geom) const;

   void print_surface() const;                            // print surface values

 protected:
   t_glog*        _log;
   t_gtime        _epoch;
   t_map_metsurf  _map_surf;      // surface values
   t_map_metsurf  _map_rms;       // surface value rms
   t_map_metsurf  _map_iter;      // surface iterations
   t_map_metdata  _map_data;      // surface values
   t_map_metdata  _map_drms;      // surface value rms

   set<MET_SURF>  _set_sf;        // settings
   set<MET_DATA>  _set_dt;        // settings

   double         _decay_min;     // min value for decay parameters
   double         _decay_max;     // max value for decay parameters
   PROFTEMP       _vert_temp;     // TEMP vertical scaling model 
   PROFZWD        _vert_zwd;      // ZWD vertical scaling model
   SURFZWD        _surf_zwd;      // ZWD surface model

};

} // namespace

#endif
