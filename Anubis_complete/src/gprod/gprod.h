
#ifndef GPROD_H
#define GPROD_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2011-03-25 /JD: created

-*/

#include <set>
#include <iostream>

#include "gdata/gdata.h"
#include "gdata/gobj.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

static shared_ptr<t_gobj> nullobj;

class t_gprod : public t_gdata {   

 public:

  t_gprod(const t_gtime& t, shared_ptr<t_gobj> pt = nullobj);
  virtual ~t_gprod();

  typedef map<string, pair<double,double> > t_map_prod;

  t_gtime epoch()const{ return _epo; }
  void    epoch(const t_gtime& epo) { _epo = epo; }

  shared_ptr<t_gobj> obj()const{ return _obj; }
  string  obj_id()const{ if(_obj) return _obj->id(); else return ""; }

  int set_val(const string& str, const double& val, const double& rms = 0.0);
  int get_val(const string& str,       double& val,       double& rms);
  int get_val(const string& str,       double& val);  

  set<string> list_id();

  int  nSat()              {return _nSat;}
  void nSat(const int& n)  {_nSat = n;}

  int  nSat_excl()             {return _nSat_excl;}
  void nSat_excl(const int& n) {_nSat_excl = n;}
   
 protected:

  t_gtime                    _epo;
  shared_ptr<t_gobj>         _obj;
  t_map_prod                 _prod;
  t_map_prod::const_iterator itPROD;

  int                        _nSat;
  int                        _nSat_excl;

 private:

};

} // namespace

#endif
