
#ifndef GNOTE_H
#define GNOTE_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements motes (messages,warning,errors)
  Version: $ Rev: $

  2017-08-04 /JD: created

-*/

#include <stdio.h>
#include <fstream>
#include <string>

#include "gio/glog.h"
#include "gutils/gmutex.h"

#define BUF_SIZE 1024

using namespace std;

namespace gnut {

enum t_note{ GERROR, GWARNING, GMESSAGE };

class t_gnote {
   
 public:
   t_gnote(t_note n, string s);
   virtual ~t_gnote();
   
   string str()const{ return _str(); }
   string text()const{ return _text; }
   t_note status()const{ return _stat; }

   bool            operator<(const t_gnote& n) const;
   bool            operator==(const t_gnote& n) const;
   friend ostream& operator<<(ostream& os, const t_gnote& n);

 protected:
   virtual string _str()const;
   
   string  _text;          // note text
   t_note  _stat;          // note status

 private:
     
};

} // namespace

#endif
