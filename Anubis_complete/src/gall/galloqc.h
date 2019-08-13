
#ifndef GALLOQC_H
#define GALLOQC_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: observation container for QC
  Version: $ Rev: $

  2016-03-10 /JD: created

-*/

#include "gall/gallobs.h"
#include "gset/gsetqc.h"

using namespace std;

namespace gnut {  

static const set<double> ELEVATIONS = { 0, 5, 10, 15, 20, 30, 50, 70, 90 };

class t_galloqc : public t_gallobs {

 public:
  t_galloqc();
  virtual ~t_galloqc();

  virtual void gset(t_gsetbase*);
   
  typedef map<string, pair<int,int> >          t_map_bsat;
  typedef map<t_gtime, t_map_bsat>             t_map_bepo;
  typedef map<GSYS,    t_map_bepo>             t_map_bsys;
   
  typedef map<string, pair<int,int> >          t_map_stt_sat;  // data statistics (have,total)
  typedef map<GOBS,   t_map_stt_sat>           t_map_stt_obs;
  typedef map<GSYS,   t_map_stt_obs>           t_map_stt_sys;
  typedef map<double, int>                     t_map_stt_smp;

  typedef map<int,             int >           t_map_sky_ele;  // elevation statistics
  typedef map<GOBS,   t_map_sky_ele>           t_map_sky_obs;
  typedef map<GSYS,   t_map_sky_obs>           t_map_sky_sys;

  static  GOBSATTR select_attr(const t_map_stt_sys& stt,
                         GSYS gs, GOBSTYPE gt, GOBSBAND gb );  // attr selection

  virtual void beg(const t_gtime& t){ _beg = t; }              // set first effective t_gobs epoch (for all functions)
  virtual void end(const t_gtime& t){ _end = t; }              // set last  effective t_gobs epoch (for all functions)

  virtual t_gtime beg_obs(const string& site, double smpl = 0.0 );  // get first t_gobs epoch for site
  virtual t_gtime end_obs(const string& site);                      // get last  t_gobs epoch for site

  virtual t_map_stt_sys stat_obs(string site,                  // get obs statistics
				 double* smp_est,                                      // set sampling (get if zero!)
				 double* min_ele,                                      // get minimum elevation
				 double  cut_off,                                      // req mask for elevation cut-off
//       const t_gtime& beg = FIRST_TIME,                      // req start period
   			 const t_gtime& end = LAST_TIME,                       // req end period
		 	   t_map_sky_sat* sky_sat = 0,                           // based on skyview
		 	   t_map_stt_smp* smp_stt = 0,                           // get sample histogram
		 	   t_map_sky_sys* sky_stt = 0);                          // get obs statistics based on skyview
			                      
  virtual int gaps(string site, int g, map<t_gtime,int>& l, double& min, double& max); // gaps number
  virtual int block(string site, int b, int g, map<t_gtime,int>& l);                   // short block data
  virtual int nbands(string site, t_map_bsys& nb, double sampl = 0.0);                 // phase/code band number
  virtual int add_elevations(string site, t_gtriple xyz_rec, t_gallnav* gnav);         // add elevations for site
  virtual int apr_elevations(string site, t_gtriple xyz_rec, t_gallnav* gnav);         // add elevations for site
  virtual int est_elevations(string site, t_gtriple xyz_rec, t_gallnav* gnav,
                             const t_gtime& beg, const t_gtime& end,                   // should be setup reasonably !
                             const double& mask,
                             t_map_sky_sat& sat_zero,                                  // out: mask for zero cut-off
                             t_map_sky_sat& sat_mask );                                // out: mask for user cut-off

 protected:

   t_gtime      _beg;           // effective data begin (for all functions)
   t_gtime      _end;           // effective data end   (for all functions)
   double       _cutoff;
   USE_HEALTH   _useHealth;
   
};

} // namespace

#endif
