
#ifndef GALLPROD_H
#define GALLPROD_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: container for all products
  Version: $ Rev: $

  2014-05-14 /PV created

-*/

#include <map>
#include <set>
#include <string>
#include <memory>
#include <string>

#include "gprod/gprod.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {  

// ----------
class t_gallprod : public t_gdata {

 public:
  t_gallprod( );
  virtual ~t_gallprod();
  
  typedef map<t_gtime, shared_ptr<t_gprod>>  t_map_epo;
  typedef map<ID_TYPE, t_map_epo> t_map_id;
  typedef map<string, t_map_id>  t_map_prd;
   
  
  int add(shared_ptr<t_gprod> prod, string site = "");
  shared_ptr<t_gprod> get(const string& site, ID_TYPE type, const t_gtime& t);
  void rem(const string& site, ID_TYPE type, const t_gtime& t);
   
  set<string>   prod_sites();
  set<ID_TYPE>  prod_types(const string& site);
  set<t_gtime>  prod_epochs(const string& site, ID_TYPE type);

  void clear();
  void clean_outer( const t_gtime& beg = FIRST_TIME, const t_gtime& end = LAST_TIME );

  double iono(const double& lat, const double& lon, const t_gtime& epo, bool interp = true);

 protected:
  shared_ptr<t_gprod> _find(const string& site, ID_TYPE type, const t_gtime& t); // find element

  t_map_prd  _map_prod;
};

} // namespace

#endif
