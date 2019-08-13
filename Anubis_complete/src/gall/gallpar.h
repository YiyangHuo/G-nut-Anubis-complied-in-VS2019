
#ifndef GALLPAR_H
#define GALLPAR_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: parametres container
  Version: $ Rev: $

  2011-04-18 /PV: created

-*/

#include <string>
#include <set>

#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gmodels/gpar.h"

using namespace std;

namespace gnut {  

class t_gallpar
{
 public:
   void addParam(t_gpar);
   void delParam(int i);
   void delAllParam();
   set<int> delEmpty();
   vector<int> delAmb();
   void resetAllParam();
   int getParam(string site, t_gpar::t_type, string prn,
       t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;
   int getParam(int index);
   unsigned int parNumber() const;
   unsigned int ambNumber() const;   
   int getCrdParam(string site, t_gtriple& crd,
		   t_gtime beg = FIRST_TIME, t_gtime end = LAST_TIME) const;

   ColumnVector get_cvect(t_gallpar& par);
   t_gpar&   operator[](const size_t idx);
   t_gallpar operator-(t_gallpar& par);
   t_gallpar operator+(t_gallpar& par);   
   static int mult(const Matrix& K, t_gallpar& par, t_gallpar& res);
   void   reIndex();
   void   decIndex(int i);
   void   incIndex(int i);   
   void   setSite(string site);
   set<string> ambs();
   int sum(t_gallpar& X1, t_gallpar& X2);
   
   friend ostream& operator<<(ostream& os, t_gallpar& x);
   
 private:
   vector<t_gpar> _vParam;
};

} // namespace

#endif
