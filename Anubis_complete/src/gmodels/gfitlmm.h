
#ifndef GFITLMM_H
#define GFITLMM_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements vector fitting via Levenberg-Marquardt Method
  Version: $ Rev: $

  2014-03-21 /JD: created

-*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <map>
#include <cmath>
#include <vector>

#include "gmodels/gfit.h"

namespace gnut {   
     
class t_gfitlmm : public t_gfit {
   
 public:
           t_gfitlmm( t_gfunc* model );
  virtual ~t_gfitlmm(){};
      
  virtual int fit( const ColumnVector&   X, const ColumnVector&   Y,
                         ColumnVector& par,       ColumnVector& err, int& i, string& msg);

  virtual int fit( const ColumnVector&   X, const ColumnVector&   Y, const DiagonalMatrix&  P,
                         ColumnVector& par,       ColumnVector& err, int& i, string& msg);

 protected:

  double        _lambda;      // LMM - lambda parameter
  double        _limit;
  double        _lim_grd;     // convergence criterium limit for gradient
};

} // namespace

#endif
