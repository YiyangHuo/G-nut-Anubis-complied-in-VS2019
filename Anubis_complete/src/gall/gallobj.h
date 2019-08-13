
#ifndef GALLOBJ_H
#define GALLOBJ_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: container for all objects
  Version: $ Rev: $

  2012-12-04 /JD: created

-*/

#include <vector>
#include <memory>

#ifdef BMUTEX
  #include <boost/thread/mutex.hpp>
#endif

#include "gdata/grec.h"
#include "gdata/gtrn.h"
#include "gutils/gsysconv.h"
#include "gutils/gtypeconv.h"
#include "gall/gallotl.h"

using namespace std;

namespace gnut {  

// ----------
class t_gallobj : public t_gdata {

 public:      
  t_gallobj();
  t_gallobj(t_gallpcv* pcv, t_gallotl* otl);
  virtual ~t_gallobj();

  typedef map<string, shared_ptr<t_gobj>> t_map_obj;          // allocated data of a single object
   
  void setPCV(t_gallpcv*  pcv);
  void setOTL(t_gallotl*  otl);   

  int add(shared_ptr<t_gobj> obj);                             // set/get single obj element
  shared_ptr<t_gobj> obj(string s);

  void sync_pcvs();                                            // synchronize PCVs
  void read_satinfo(t_gtime& epo);
   
//virtual vector<string>             objects_id(t_gdata::ID_TYPE id=NONE);  // get all object names
//virtual vector<shared_ptr<t_gobj>> objects(t_gdata::ID_TYPE id=NONE);     // get all object elements
  virtual map<string,shared_ptr<t_gobj>> objects(t_gdata::ID_TYPE id=NONE); // get all object elements
//virtual vector<string> obj_sat();                            // get all satellites
//virtual vector<string> obj_rec();                            // get all receivers
 
 virtual int obj_num();                                        // get object numbera
   
 private:
  t_map_obj      _mapobj;                                      // map of all objects
  t_gallpcv*     _gpcv;
  t_gallotl*     _gotl;
  void           _aloctrn();
};

} // namespace

#endif
