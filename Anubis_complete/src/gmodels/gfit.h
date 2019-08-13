
#ifndef GFIT_H
#define GFIT_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements vector fitting (pure virtual function)
  Version: $ Rev: $

  2014-03-21 /JD: created

-*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <string>
#include <cmath>
#include <vector>
#include <string>
#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

#include "gio/glog.h"
#include "gmodels/gfunc.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"

#define MAX_FIT_ITER 100

using namespace std;

namespace gnut {   

class t_gfit
{   
 public:
           t_gfit( t_gfunc* model, int max_iter = MAX_FIT_ITER );
  virtual ~t_gfit(){};
      
  void glog(t_glog* glog){ _log = glog; }
  virtual int fit( const ColumnVector&   X, const ColumnVector&   Y,
                         ColumnVector& par,       ColumnVector& err, int& i, string& msg) = 0;

 protected:
  t_glog*       _log;
  t_gfunc*      _func;
  double        _lim_par;       // convergence criterium limit for parameter
  int           _max_iter;      // max iterations

};

} // namespace

#endif
