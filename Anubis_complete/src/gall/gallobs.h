
#ifndef GALLOBS_H
#define GALLOBS_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: observation container
  Version: $ Rev: $

  2012-05-02 /JD: created

-*/

#include <iostream>
#include <string.h>
#include <map>
#include <set>
#include <memory>

#include "gdata/gdata.h"
#include "gdata/gobsgnss.h"
#include "gdata/gsatgnss.h"
//#include "gdata/gsdbs.h"
#include "gdata/gsatdata.h"
#include "gutils/gsatview.h"
#include "gutils/gtypeconv.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gutils/gnss.h"
#include "gutils/gsys.h"
#include "gutils/gobs.h"
#include "gset/gsetgen.h"

#define DIFF_SEC_NOMINAL 0.905  // [sec] returns observations within +- DIFF_SEC for 1Hz

// normalize for 1Hz !
#define DIFF_SEC(a) (((a)>(0.0)&&(a)<(1.0))?(DIFF_SEC_NOMINAL*a):(DIFF_SEC_NOMINAL))

using namespace std;

namespace gnut {  

typedef shared_ptr<t_gobsgnss> t_spt_gobs;
   
class t_gallobs : public t_gdata {

 public:
  t_gallobs();
  virtual ~t_gallobs();

  enum XDATA { XDATA_BEG, XDATA_END, XDATA_SMP, XDATA_SYS };    // QC

  typedef map<XDATA, int>                      t_map_xdat;      // map of site filtered data (QC)
  struct  t_xfilter { t_map_xdat xdat; t_gtime beg, end; };     // filtered data (QC)
  typedef map<string, t_xfilter>               t_map_xfil;      // map of file filtered data (QC)
  typedef map<string, t_map_xfil>              t_map_xflt;      // map of  all filtered data (QC)


  typedef map<string, t_spt_gobs>              t_map_osat;      // all data-types/single epoch
  typedef map<t_gtime,t_map_osat>              t_map_oref;      // all data-types/all epochs/single object
  typedef map<string, t_map_oref>              t_map_oobj;      // all data-types/all epochs/all objects

  typedef map<string, map<GOBSBAND, map<GOBS, int>>> t_map_frq; // signals occurance
   
  virtual void gset(t_gsetbase*);
   
  virtual set<string> stations();                                                // get all station
  virtual set<GSYS> sys(string site);
  virtual set<GSYS> sys(const string& site, const t_gtime& t);                   // get all systems for epoch

  virtual bool isSite(string site);
   
  virtual set<string> sats(const string& site, const t_gtime& t, GSYS gnss);     // get all satellites for epoch t and system

//  virtual t_gsdbs sdbs(const string& site, const t_gtime& t);                  // single differ. btw sats in epoch t

  virtual vector<t_gsatdata> obs(const string& site, const t_gtime& t);          // get all t_gsatdata for epoch t
  virtual vector<t_spt_gobs> obs_pt(const string& site, const t_gtime& t);       // get all t_gobsgnss pointers for epoch t
  virtual vector<t_spt_gobs> obs_prn_pt(const string& site, const string& prn, 
                                        const t_gtime& beg, const t_gtime& end); // get all t_gobsgnss pointers for prn in interval
  virtual vector<t_gtime> epochs(const string& site);                            // get all t_gepochs for site   
   
  virtual void clean_outer(const string& site = "",
                           const t_gtime& beg = FIRST_TIME,
	                         const t_gtime& end = LAST_TIME);

  virtual t_gtime beg_obs(const string& site, double smpl = 0.0 );  // get first t_gobs epoch for site
  virtual t_gtime end_obs(const string& site);                      // get last  t_gobs epoch for site

  void      xdata(string site, string file, t_xfilter xflt);        // add site-specific filtered data/epochs
  t_xfilter xdata(string site, string file);                        // get site-specific filtered data/epochs

  int addobs(t_spt_gobs obs);                                       // add single station observation (P and L in meters !)

  void overwrite(bool b){ _overwrite = b; }                         // set/get overwrite mode
  bool overwrite(){ return _overwrite; } 

  void maxepoch( unsigned int i ){ _nepoch = i; }                   // set/get maximum number of epochs stored
  unsigned int maxepoch(){ return _nepoch; }

  unsigned int nepochs(const string& site);                         // get number of epochs for station
  unsigned int nepochs(const string& site, const t_gtime& beg, const t_gtime& end, double sampl, map<GSYS, pair<int,int> >& n);  // get number of epochs for station expect/have according to sampl

  virtual t_map_osat find( string site, const t_gtime& t);          // find appropriate t_gobsgnss element for site/epoch
  virtual double     find( string site, const t_gtime& t, const string& prn, const GOBS& gobs);
   
  t_map_frq frqobs(string site);                                    // get number of occurance of individual signals
   
 protected:
  virtual set<string>        _sats(const string& site, const t_gtime& t, GSYS gnss);
  virtual set<GSYS>          _gsys(const string& site, const t_gtime& t);
  virtual vector<t_gsatdata> _gobs(const string& site, const t_gtime& t);

  int _find_epo(const string& site, const t_gtime& epo, t_gtime& tt);  // find epoch from the map w.r.t. DIFF_SEC

  t_gsetbase*       _set;
  unsigned int      _nepoch;      // maximum number of epochs (0 = keep all)
  t_map_oobj        _mapobj;      // map over all objects (receivers)
  t_map_xflt        _filter;      // structure of stations/files filtered data (QC)
  bool              _overwrite;   // rewrite/add only mode
  set<string>       _sys;         // systems settings
  double            _smp;         // sampling settings
  double            _scl;         // sampling scale-factor

 private:
};

} // namespace

#endif
