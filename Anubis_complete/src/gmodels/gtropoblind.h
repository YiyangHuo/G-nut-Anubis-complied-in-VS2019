
#ifndef GTROPOBLIND_H
#define GTROPOBLIND_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements troposphere model class (blind models)
  Version: $ Rev: $

  2016-01-10 /JD: created

-*/

#include "gmodels/gtropo.h"
#include "gmodels/gblindmops.h"
#include "gmodels/ggpt.h"

using namespace std;

namespace gnut {   

// ----------------------------- BLIND MODELS ------------------------------------
class t_blindmops : public t_gtropo {
 public:
   t_blindmops(){ _model = new t_mops(); }
  ~t_blindmops(){}

   virtual double getZHD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
   
 private:
   t_mops* _model;
};

} // namespace

#endif
