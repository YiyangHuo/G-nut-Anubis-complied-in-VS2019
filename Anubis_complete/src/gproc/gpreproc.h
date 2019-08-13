
#ifndef PREPROCESSING_H
#define PREPROCESSING_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: preprocessing, header
  Version: $ Rev: $

  2012-09-03 /PV: created

-*/

#include <sstream>

#include "../newmat/newmat.h"

#include "gall/gallobs.h"
#include "gall/gallnav.h"
#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"
#include "gset/gsetproc.h"

using namespace gnut;

namespace gnut {

class t_gobs_pair;

class t_gpreproc : public t_gmonit
{
 public:
    t_gpreproc(t_gallobs* obs, t_gsetbase* settings);
   ~t_gpreproc();

   void glog(t_glog* l){ _log = l; }                 // set/get glog pointer
   t_glog* glog()const{ return _log; }

   int preprocess(string site, const t_gtime& beg, const t_gtime& end, double sampl, bool sync, bool save = false);
   void setNav(t_gallnav* nav);

   void setSite(string site);
   
   map<t_gtime, map<string, map<GOBS, double> > > getSlips();
   map<t_gtime, map<string, map<GOBS, int> > > getSlipsGap();   
   map<t_gtime, int> getClkjump();

   typedef map<int, map<t_gobs_pair, double> > t_map_slp;
   typedef map<int, vector<t_gobs_pair> > t_vec_slp;   
   
 protected:
   t_gallobs*         _obs;
   t_gallnav*         _nav;   
   string             _site;
   t_glog*            _log;
   t_gsetbase*        _set;
   double             _sigCode, _sigPhase;
   double             _sigCode_GLO, _sigPhase_GLO;
   double             _sigCode_GAL, _sigPhase_GAL;
   double             _sigCode_BDS, _sigPhase_BDS;
   double             _sigCode_QZS, _sigPhase_QZS;
   double             _sigCode_IRN, _sigPhase_IRN;   
   double             _sumS;   
   double             _scl;
   t_gtriple          _pos;
   set<string>        _sys;   
   set<string>        _sat;
   
   bool               _beg_end;

   map<string, double>  _dI;                                   // ionospheric delay   
   map<string, int>     _msoffset;
   map<t_gtime, int>    _mbreaks;                              // for logging
     
   t_map_slp _m_lcslp;
   t_vec_slp _v_lcslp;   
   
   map<t_gtime, map<string, map<GOBS, double> > > _mslips;    // for logging - slip due to true CS
   map<t_gtime, map<string, map<GOBS, int> > > _mslipsGap;    // for logging - slip due to data gap => 1 (epoch gap), 2(sat gap), 3(GOBS gap)

   int      _slip(t_gobsgnss* gobs1, t_spt_gobs gobs2);
   double   _slipThreshold(string LC, t_spt_gobs obs, t_gband& band1, t_gband& band2, double sampl = 30.0);
   int      _jumps(t_gobsgnss* gobs1, t_spt_gobs gobs2);
   void     _repair(vector<t_spt_gobs> epoData, double dL);   

   int      _transform(t_spt_gobs gobs,  bool save);   
   
   void     _save(t_spt_gobs gobs, const map<GOBS,double>& slips);
   void     _remove_slip(vector<t_spt_gobs> gobs);

   void _common(set<GOBSBAND>& set1, set<GOBSBAND>& set2);
   void _common(set<GOBS>& set1, set<GOBS>& set2);
   double _disf(t_gobsgnss* gobs1, t_spt_gobs gobs2, t_gobs& s1, t_gobs& s2);
   void   _iono(t_gobsgnss* gobs1, t_spt_gobs gobs2, t_gobs& s1, t_gobs& s2);
   double _findSlp(int& i, t_gobs_pair& gpair);
   void   _gapReport(vector<shared_ptr<t_gobsgnss>> epoData);
   void   _gapReport(shared_ptr<t_gobsgnss> data);   
   void   _compare(vector<t_gobsgnss> data1, vector<shared_ptr<t_gobsgnss>> data2);
};

class t_gobs_pair
{
 public:
   t_gobs_pair(t_gobs& gobs1, t_gobs& gobs2);
   t_gobs obs1;
   t_gobs obs2;
   double val;
   bool operator<(const t_gobs_pair& t) const;
};

} // namespace

#endif
