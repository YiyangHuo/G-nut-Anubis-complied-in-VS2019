
#ifndef GSETREC_H
#define GSETREC_H

#define XMLKEY_REC "rec"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements receiver object settings
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <map>
#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gdata/grec.h"
#include "gset/gsetbase.h"
#include "gutils/gtypeconv.h"
#include "gutils/gtriple.h"
#include "gmodels/ggpt.h"

#define HSL_UNKNOWN -9999  // unknown height ebove mean see level

using namespace std;

namespace gnut {

class t_gsetrec : public virtual t_gsetbase
{
 public:
   t_gsetrec();
  ~t_gsetrec();

   void check();                                  // settings check
   void help();                                   // settings help

   string id();

   int get_crd_xyz(t_gtriple& xyz, string s);

   set<string> objects();                         // get all objects IDs
   shared_ptr<t_grec> grec(string s, t_glog* glog = 0);  // ID: return grec object for selected ID
   
   // all below can be removed in future
   string domes(string s);                        // ID: domes
   string name(string s);                         // ID: name [id]
   string desc(string s);                         // ID: full name description
   string   id(string s);                         // ID: internal name
   t_gtime beg(string s, string t = "");          // ID: start time
   t_gtime end(string s, string t = "");          // ID: end time
   string  rec(string s, string t = "");          // ID: receiver name
   string  ant(string s, string t = "");          // ID: antenna name

 protected:
   t_gtriple _get_crd_xyz(string s);
   t_gtriple _get_ecc_neu(string s);
   t_gtriple _get_ecc_xyz(string s);   
   t_gtriple _get_crd_blh(string s);              // lat, lon, hmsl(above see level)

   t_gpt  _ggpt;                                  // to know undulation (if HSL only provided)
   double _aprX(string s);
   double _aprY();
   double _aprZ();
   double _aprDX();
   double _aprDY();
   double _aprDZ();
   double _aprDE();
   double _aprDN();
   double _aprDU();
   double _aprZTD();
   double _sigZTD();
   string _sigCRD();
   set<string> _objects();                        // get all objects names

   string _rec;                                   // default receiver name
   string _ant;                                   // default antenna name
   string _id;                                    // receiver id
   string _name;                                  // receiver name
   string _desc;                                  // receiver description
   string _domes;                                 // receiver monumentation domes
   t_gtime _beg;                                  // default begin time
   t_gtime _end;                                  // default end time
   double _X;                                     // receiver X-coordinate [m] 
   double _Y;                                     // receiver Y-coordinate [m] 
   double _Z;                                     // receiver Z-coordinate [m]
   double _dX;                                    // receiver X-eccentricity [m]
   double _dY;                                    // receiver Y-eccentricity [m] 
   double _dZ;                                    // receiver Z-eccentricity [m]
   double _dE;                                    // receiver E-eccentricity [m]
   double _dN;                                    // receiver N-eccentricity [m]
   double _dU;                                    // receiver U-eccentricity [m]   
   double _ZTD;
   double _ZTD_sig;
   string _CRD_sig;
   bool   _overwrite;
   
 private:
};

} // namespace

#endif
