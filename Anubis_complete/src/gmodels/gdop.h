
#ifndef GDOP_H
#define GDOP_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: dilution of precesion - gdop, pdop, hdop, vdop, tdop

  Version: $Rev:$

  2014-01-21 /PV: created

-*/

#include <iostream>
#include <set>

#include "../newmat/newmat.h"
#include "gall/gallnav.h"
#include "gall/gallobs.h"
#include "gutils/gtriple.h"

using namespace std;

namespace gnut {   

class t_gdop {
 public:
   t_gdop();
   t_gdop(t_gallnav* gnav, t_gallobs* gobs, string site);
   t_gdop(t_gallnav* gnav, set<string> sats);
  ~t_gdop();
      
   int    calculate(const t_gtime& epoch, t_gtriple& rec, GSYS gnss = GNS);
   double pdop();
   double gdop();   
   double tdop();
   double hdop();
   double vdop();
   void   set_data(t_gallnav* gnav, t_gallobs* gobs, string site);
   void   set_log(t_glog* glog);
   void   set_sats(set<string>& sats);
   
 private:       
   string          _site;   // site name
   t_gallnav*      _gnav;   // ephemerides
   t_gallobs*      _gobs;   // observations
   set<string>     _sats;   // set of visible satellites
   t_gtriple       _rec;    // receiver position
   SymmetricMatrix _Qx;     // variance-covariance matrix
   t_glog*         _log;
};

} // namespace

#endif
