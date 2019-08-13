
#ifndef GSETVAR_H
#define GSETVAR_H

#define XMLKEY_VAR "var"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements processing setting
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <string>
#include <iostream>

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "gset/gsetbase.h"

using namespace std;
using namespace pugi;

namespace gnut {
   
class t_gsetvar : public virtual t_gsetbase
{
 public:
   t_gsetvar();
  ~t_gsetvar();

  void   set_attr(string key, string val);
  string get_attr(string key);

 protected:
   
 private:
};

} // namespace

#endif
