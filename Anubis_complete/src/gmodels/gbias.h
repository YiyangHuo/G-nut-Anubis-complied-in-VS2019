
#ifndef GBIAS_H
#define GBIAS_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements GNSS code biases
  Version: $ Rev: $

  2012-11-05 /JD: created

-*/

#include <vector>

#include "../newmat/newmat.h"
#include "gdata/gdata.h"
#include "gdata/gobsgnss.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {   

// ----------
class t_gbias : public t_gdata {

   public:
   t_gbias();
   virtual ~t_gbias();
   
   void  set(t_gtime beg, t_gtime end, double d, GOBS obs1, GOBS obs2 = X);    // add single differential bias
   void  set(double d, GOBS obs1, GOBS obs2 = X);                                // add single differential bias
   void  ref(GOBS ref);
   
   double bias(bool meter = true);                                          // get signgle differential bias
   
   GOBS gobs()const{ return _gobs; }
   GOBS ref()  const{ return _ref;   }      
   
   void    beg(t_gtime t){ _beg = t; }          // set/get valid from
   t_gtime beg()const{ return _beg;  }
   void    end(t_gtime t){ _end = t; }          // set/get valid until
   t_gtime end()const{ return _end;  }
   
   bool valid(const t_gtime& epo);
   
   private:

   t_gtime        _beg;        // valid from
   t_gtime        _end;        // valid until
   
   GOBS           _gobs;
   GOBS           _ref;
   double         _val;
 

};

} // namespace

#endif
