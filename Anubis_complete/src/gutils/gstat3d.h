
#ifndef GSTAT3D_H
#define GSTAT3D_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: various statistics function
  Version: $ Rev: $

  2014-04-18 /PV: created

-*/

#include "gutils/gstat.h"
#include "gutils/gtriple.h"

#include <vector>

using namespace std;

namespace gnut {

class t_gstat3d : public t_gstat
{
 public:
   t_gstat3d(vector<t_gtriple>& data, double conf_interv = CONF_INTERV);

   int calc_median();
   int calc_stat();

   t_gtriple  get_std3d();
   t_gtriple  get_mean3d();
   t_gtriple  get_median3d();
   t_gtriple  calc_mean();
   t_gtriple  calc_std();
   
 protected:
   void       _add_data(vector<t_gtriple>& data);
   
   vector<t_gtriple> _data3d;
   t_gtriple         _std3d;
   t_gtriple         _mean3d;
   t_gtriple         _median3d;
};

} // namespace

#endif
