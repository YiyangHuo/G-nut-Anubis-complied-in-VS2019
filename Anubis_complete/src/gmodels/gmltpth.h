
#ifndef MLTPTH_H
#define MLTPTH_H

#define MP_UNITS   (100) // multipath [m]->[cm]
#define W(x)  setw(x)    // print field width

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: multipath estimation
           the class is derived from t_gxtr
  Version: $ Rev: $

  2012-11-05 /PV: created
 
-*/

#include "gall/gallobs.h"
#include "gset/gsetbase.h"
#include "gset/gsetgen.h"

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

using namespace std;

namespace gnut {   

// ----------
class t_gmltpth {

 public:
   t_gmltpth(t_gallobs* gobs, t_gsetbase* gset, string site);
   virtual ~t_gmltpth();

   typedef map<string, pair<double,int> > t_map_mpsat;
   typedef map<t_gtime, t_map_mpsat >     t_map_mpepo;
   typedef map<GOBS, t_map_mpepo >        t_map_mpobs;        // used only for re-ordering
   typedef map<GSYS, t_map_mpobs >        t_map_mpsys;        // used only for re-ordering   
   
   void multipath(const string& prn, GSYS gsys, const t_gtime& epo,
		  const unsigned int& nepo, t_map_mpobs& m_obs);   

 protected:
   double       _sampling;
   double       _mp_limit;
   t_gallobs*   _gobs;
   t_gsetbase*  _set;
   string       _site;
   
   
#ifdef BMUTEX
   boost::mutex _mutex;
#endif   
};

} // namespace

#endif
