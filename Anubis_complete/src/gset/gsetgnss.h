
#ifndef GSETGNSS_H
#define GSETGNSS_H

#define XMLKEY_GNSS "gnss"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements data extraction setting
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <map>
#include <vector>
#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gsys.h"
#include "gutils/gobs.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut {

class t_gsetgnss : public virtual t_gsetbase
{
 public:
   t_gsetgnss();
  ~t_gsetgnss();
  
   void check();                                 // settings check
   void help();                                  // settings help

   set<string>  sat();                            // all systems
   set<string>  sat(GSYS gsys, bool def=true);    // single system def=true (give default values at least)
   set<string>  obs(GSYS gsys, bool def=true);    // single system def=true (give default values at least)
   set<string>  nav(GSYS gsys, bool def=true);    // single system def=true (give default values at least)

   set<string> gobs(GSYS gsys);                   // extending gobs list with complete singals
   
   vector<GNAVTYPE> gnav(GSYS gsys);
   vector<GOBSTYPE> type(GSYS gsys);
   vector<GOBSBAND> band(GSYS gsys);
   vector<GOBSATTR> attr(GSYS gsys);
   
   double           sigma_L(GSYS gsys);
   double           sigma_C(GSYS gsys);
   
   double           maxres_L(GSYS gsys);
   double           maxres_C(GSYS gsys);   

 protected:
   vector<GOBSTYPE> _type(GSYS gsys);
   vector<GOBSBAND> _band(GSYS gsys);
   vector<GOBSATTR> _attr(GSYS gsys);
   string           _gsys(GSYS gsys);
   double           _sigma_L(GSYS gsys);
   double           _sigma_C(GSYS gsys);    
   double           _maxres_L(GSYS gsys);
   double           _maxres_C(GSYS gsys);
   
   map<GSYS,vector<string> > _band_str;           // default set
   map<GSYS,vector<string> > _type_str;           // default set
   map<GSYS,vector<string> > _attr_str;           // default set

   map<GSYS, t_gpair>        _sigma_def;          // default set
   map<GSYS, t_gpair>        _maxres_def;         // default set

// map<GSYS,vector<GOBSBAND> > _band_obs;         // default set
// map<GSYS,vector<GOBSTYPE> > _type_obs;         // default set
// map<GSYS,vector<GOBSATTR> > _attr_obs;         // default set

 private:
};

} // namespace

#endif
