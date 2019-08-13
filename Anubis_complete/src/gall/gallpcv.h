
#ifndef GALLPCV_H
#define GALLPCV_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $ Rev: $

  2011-02-14 /JD: created

-*/

#include <iostream>
#include <string.h>

#include "gdata/gdata.h"
#include "gutils/gsys.h"
#include "gutils/gconst.h"
#include "gutils/gtime.h"
#include "gmodels/gpcv.h"

using namespace std;

namespace gnut {  

// ----------
class t_gallpcv : public t_gdata {
 
  public:
   t_gallpcv();
   virtual ~t_gallpcv();

   typedef map<t_gtime, shared_ptr<t_gpcv>> t_map_epo;     // map of 1-ant/1-sn/N-epochs
   typedef map<string,  t_map_epo>          t_map_num;     // map of 1-ant/N-sn/N-epochs
   typedef map<string,  t_map_num>          t_map_pcv;     // map of N-ant/N-sn/N-epochs

   shared_ptr<t_gpcv> find( string ant, string ser, 
			    const t_gtime& t );                  // find appropriate t_gpcv element
     
   int addpcv( shared_ptr<t_gpcv> pcv );                         // add single antenn pattern (PCV)
   shared_ptr<t_gpcv> gpcv( string ant, string num, 
			    const t_gtime& t );                  // get single antenn pattern (PCV)

   vector<string> antennas();                                    // returns vector of all antennas

   void overwrite(bool b){  _overwrite = b; }                    // set/get overwrite mode
   bool overwrite(){ return _overwrite; }

 protected:
   virtual shared_ptr<t_gpcv>
           _find( string ant, string ser, const t_gtime& t );    // find appropriate t_gpcv element
   
 private:
   t_map_pcv      _mappcv;     // complete PCV-map
   bool           _overwrite;  // rewrite/add only mode

   shared_ptr<t_gpcv> _pcvnull;

};

} // namespace

#endif
