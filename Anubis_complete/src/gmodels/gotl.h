
#ifndef  GOTL_H
#define  GOTL_H

/* ----------------------------------------------------------------------
   (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
       Research Institute of Geodesy, Topography and Cartography
       Ondrejov 244, 251 65, Czech Republic
   
   Purpose: implementation of ocean tide loading
   Version: $Rev:$
       
   2013-03-29 /PV: created
   
-*/

#include <string>

#include "../newmat/newmat.h"

#include "gdata/gdata.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {   

class t_gotl : public t_gdata {
   
 public:   
   t_gotl();
   ~t_gotl();
   
   string site();
   double lat();
   double lon();
   Matrix data();
   void setdata(const string& site, const double& lon, const double& lat, const Matrix& data);
   
 private:
   string _site;
   Matrix _data;
   double _lat;
   double _lon;
};

} // namespace

#endif
