
#ifndef GPROJ_H
#define GPROJ_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: Spherical to projection transformations
  Version: $ Rev: $

  2015-02-19 /JD: created

-*/

#include <string>
#include <iomanip>
#include <iostream>

#include "gutils/gpair.h"
#include "gutils/gtriple.h"
#include "gutils/gconst.h"

using namespace std;

namespace gnut {

enum PROJ { PROJ_LCC, PROJ_LCI, PROJ_UNI, PROJ_ZERO };
   
class t_gproj {
      
 public:
           t_gproj(){ _name = ""; _id = PROJ_ZERO; }
  virtual ~t_gproj(){};

  virtual t_gpair   ll2proj(const t_gpair& latlon, int dec = 0 ) = 0;
  virtual t_gpair   proj2ll(const t_gpair& xy_idx              ) = 0;

  virtual t_gtriple ll2proj(const t_gtriple& latlon, int dec = 0 );
  virtual t_gtriple proj2ll(const t_gtriple& xy_idx              );
   
  virtual PROJ id(){ return _id; }
  virtual string name(){ return _name; }

  bool operator==(const t_gproj& prj){ return true; }
   
 protected:
  string    _name;
  PROJ      _id;
   
};
   
} // namespace

#endif
