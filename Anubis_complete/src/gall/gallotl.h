
#ifndef  GALLOTL_H
#define  GALLOTL_H

/* ----------------------------------------------------------------------
   (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
       Research Institute of Geodesy, Topography and Cartography
       Ondrejov 244, 251 65, Czech Republic
   
   Purpose: container along sites for BLQ data
   Version: $Rev:$
       
   2013-03-29 /PV: created
   
-*/

#include <string>
#include <map>

#include "../newmat/newmat.h"
#include "../newmat//newmatio.h"

#include "gdata/gdata.h"
#include "gutils/gtime.h"
#include "gmodels/gotl.h"

using namespace std;

namespace gnut {  

class t_gallotl : public t_gdata {
   
 public:   
   t_gallotl();
   ~t_gallotl();
   
   int data(Matrix& data, const string& site);
   double lat(const string& site);
   double lon(const string& site);
   void add(t_gotl& otl);
   void print();
 private:
   map<string, t_gotl> _mapotl;
};

} // namespace

#endif
