
#ifndef GXTR_H
#define GXTR_H

#define W(x)  setw(x)    // print field width
#define NSAT        (36) // number of svns for GNSS (considered to be reported)

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements extration class
  Version: $ Rev: $

  2013-01-10 /JD: created

-*/

#include <sstream>

#include "gio/gxml.h"
#include "gset/gsetrec.h"
#include "gset/gsetgen.h"
#include "gset/gsetout.h"
#include "gall/gallobj.h"
#include "gall/gallprec.h"

using namespace std;
using namespace pugi;

namespace gnut {

class t_gxtr : public t_gxml                
{
 public:
   t_gxtr(t_gsetbase* set, const string& pgm, const t_gtime& dt);
   t_gxtr();   
   virtual ~t_gxtr();

   virtual void navig(string site);
   
   virtual void setNAV(t_gallnav* gnav);
   virtual void setOBJ(t_gallobj* gobj);
   virtual void setSIT(string site);
   
   void setBeg(t_gtime& beg) {_set_beg = beg;}
   void setEnd(t_gtime& end) {_set_end = end;}
   
   void glog(t_glog* l){ _log = l; }                          // set/get glog pointer
   t_glog* glog(){ return _log; }
   
   virtual void check();                                      // check
   virtual void help();                                       // help
   
 protected:
   virtual void _get_settings();   
   virtual void _navig(ostringstream& os);
   virtual void _setOut();

   void _section(ostringstream& os, string tit, int verb = 0);
   void _subsec(ostringstream& os, string tit);
   void _legend(ostringstream& os, string tit, string tit2 = "", int w = 2);
   void _setkey(ostringstream& os, string key, char vrb = ' ', string epo = "XXXX-XX-XX XX:XX:XX");   
   
   t_gsetbase* _set;
   t_gallobj*  _gobj;
   t_gallnav*  _gnav;
   t_glog*     _nav;
   t_glog*     _log;
   string      _pgm; // app name
   t_gtime     _dt;  // app modification
   t_gtriple   _pos;
   string      _site;
   double      _smp_req;
   t_gtime     _ref;
   int         _nsat;   
      
   t_gtime     _set_beg, _set_end;
   t_gsetrec*  _rec;   
};
   
} // namespace
   
#endif