
#ifndef GALLNAV_H
#define GALLNAV_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $ Rev: $

  2011-02-14 /JD: created

-*/

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <memory>

#include "gdata/gnav.h"
#include "gdata/geph.h"
#include "gdata/gdata.h"
#include "gdata/gnavgps.h"
#include "gdata/gnavglo.h"
#include "gdata/gnavgal.h"
#include "gdata/gnavqzs.h"
#include "gdata/gnavsbs.h"
#include "gdata/gnavbds.h"
#include "gdata/gnavirn.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gtetrad.h"

#define MAX_GPS_PRN 32
#define MAX_GLO_PRN 24
#define MAX_GAL_PRN 30
#define NAV_BUF   1024


using namespace std;

namespace gnut {  
   
class t_gallnav : public t_gdata {

 public:
  t_gallnav();
  virtual ~t_gallnav();

  typedef multimap<t_gtime, shared_ptr<t_geph> > t_map_ref; // all data for a single satellite
  typedef map<string, t_map_ref>                 t_map_sat; // all data for all satellites
   
//  typedef map<t_gtime,double>               t_map_sky_epo; // sat visibility of a single sat/site
//  typedef map<string, t_map_sky_epo>        t_map_sky_sat; // sat visibility of a single satellite

//  virtual int visibility( t_map_sky_sat& skysat, const t_gtriple& xyz,
//		          const t_gtime& beg = FIRST_TIME,
//		          const t_gtime& end = LAST_TIME );

          void chk_health(bool b){ _chk_health = b; }        // set nav healthy status control
          void chk_navig(bool b){ _chk_navig = b; }          // set nav internal quality control
          void chk_tot(bool b) { _chk_tot = b; }             // test if tot < t
   
  virtual bool health( string sat, const t_gtime& t );       // get satellite health

  virtual int nav(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true ); // [m]
  virtual int pos(string sat, const t_gtime& t, double  xyz[3], double  var[3] = NULL, double  vel[3] = NULL, bool chk_mask = true ); // [m]
  virtual int clk(string sat, const t_gtime& t, double* clk,    double* var    = NULL, double* dclk   = NULL, bool chk_mask = true ); // [s]  

  virtual void print(string sat, const t_gtime& t);

  virtual void clean_invalid();                        // clean invalid messages
  virtual void clean_duplicit();                       // clean duplicit messages
  virtual void clean_outer( const t_gtime& beg = FIRST_TIME,
                            const t_gtime& end = LAST_TIME );

  virtual t_gtime beg_gnav( string prn = "");                            // get first t_gnav epoch
  virtual t_gtime end_gnav( string prn = "");                            // get last  t_gnav epoch
  // IMPROVE beg_time/end_time to distinguish GALLNAV/GALLPREC - t_satview !
  virtual t_gtime beg_time( string prn = ""){ return beg_gnav(prn); }    // get first t_gnav epoch
  virtual t_gtime end_time( string prn = ""){ return end_gnav(prn); }    // get first t_gnav epoch
  // =========================================================

  virtual set<string> satellites();                                  // get all satellites
  virtual set<GSYS>   systems();                                     // get all systems
  virtual unsigned int nepochs( const string& prn );                 // get number of epochs
  virtual int add( shared_ptr<t_gnav> nav );                         // add single navigation message

  virtual vector<t_gtime>            vec_epo( string prn = "" );     // get list of nav epochs
  virtual set<t_gtime>            vec_epo( GSYS gsys );           // get list of nav epochs  
  virtual vector<shared_ptr<t_geph>> vec_nav( string prn = "" );     // get list of nav messages

  virtual map<string, t_gtriple>  map_xyz(set<string> prns, const t_gtime& beg);       // get list of calculated crd   
  virtual map<string, double>     map_clk(set<string> prns, const t_gtime& beg);       // get list of calculated clk
   
   map<string, map<shared_ptr<t_geph>, t_gtriple> > multi_xyz(set<string> prns, const t_gtime& epo);
   map<string, map<shared_ptr<t_geph>, double> >    multi_clk(set<string> prns, const t_gtime& epo);
   
  virtual void multi(bool b){ _multimap = b; }                       // set/get multimap mode
  virtual bool multi(){ return _multimap; }

  virtual void overwrite(bool b){ _overwrite = b; }                  // set/get overwrite mode
  virtual bool overwrite(){ return _overwrite; }                     // .. not for multimap, but derived classe

  virtual void com(bool b){ _com = b; }                              // position/clock reference point
  virtual bool com(){ return _com; }                                 // position/clock reference point
   
  virtual void offset( int i ){ _offset = i; }                       // set/get offset for satdata
  virtual int  offset(){ return _offset; }

  virtual int nsat(GSYS gs);                                         // get number of satellites
  virtual int intv(GSYS gs);                                         // get interval between messages
  virtual int have(GSYS gs, const t_gtime& beg, const t_gtime& end); // get existing number of messages
  virtual int expt(GSYS gs, const t_gtime& beg, const t_gtime& end); // get expected number of messages
  virtual int excl(GSYS gs, const t_gtime& beg, const t_gtime& end); // get excluded number of messages

  virtual int consolidate( double cfdi = 0.0);                       // consolidate (confident interval, 0 = auto)

  virtual shared_ptr<t_geph> find( string sat, const t_gtime& t, bool chk_mask = true );   // find appropriate t_geph element (interface only)
  virtual vector<shared_ptr<t_geph>> find_mult( string sat, const t_gtime& t );   // find appropriate t_geph elements (interface only)
                                                                                  // satellite healthy status is not checked 
   
 protected:   
//  virtual int _visibility( t_map_sky_sat& skysat, const t_gtriple& xyz,
//                           const t_gtime& beg = FIRST_TIME,
//		                       const t_gtime& end = LAST_TIME );

  virtual shared_ptr<t_geph> _find( string sat,          const t_gtime& t, bool chk_mask = true );  // find appropriate t_geph element
  vector<shared_ptr<t_geph>> _find_mult( string sat,     const t_gtime& t );  // find vector of appropriate t_geph elements
                                                                              // satellite healthy status is not checked 
  virtual shared_ptr<t_geph> _find( string sat, int iod, const t_gtime& t );  // find appropriate t_gnav element

  bool               _com;         // position/clock reference point (com = true; apc = false);
  int                _offset;      // offset for RTCM corrections
  int                _nepoch;      // maximum number of epochs (0 = keep all)
  t_map_sat          _mapsat;      // map over all satellites (positions,..?)
  bool               _multimap;    // use multimap for redundant records
  bool               _overwrite;   // overwrite mode (for derived classes with MAP)
  bool               _chk_health;  // check satellite health (navigation)
  bool               _chk_navig;   // check navigation messages (internal)
  shared_ptr<t_geph> _null;        // null pointer

  bool _chk_tot;
 //private:
};

} // namespace

#endif
