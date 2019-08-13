
#ifndef FLT_H
#define FLT_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements Kalman filter class
  Version: $ Rev: $

  2011-01-10 /PV: created

-*/

#include "../newmat/newmat.h"
#include "../newmat/newmatap.h"

#include "gall/gallpar.h"

using namespace std;

namespace gnut {

// Base class for all filtering technique
// -------------------------------------------------------

class t_gflt 
{
 public:
   t_gflt();
   virtual ~t_gflt();
   
   virtual void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q) {};
   virtual void update_AmbEl(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, const SymmetricMatrix& Q, SymmetricMatrix& Q_el) {};
   virtual void smooth(t_gallpar& Xp, t_gallpar& Xu, SymmetricMatrix& Qp, SymmetricMatrix& Qu, t_gallpar& Xsm, SymmetricMatrix& Qsm, DiagonalMatrix& Noise) {};
 protected:
};

// --------------------------------------------------------
 

// Classical formule for Kalman filter
// -------------------------------------------------------
class t_kalman : public t_gflt
{
   
 public:

  t_kalman();
 ~t_kalman();
   
 void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);
};
//-------------------------------------------------------

// Square root covariance filter
// ------------------------------------------------------
class t_SRF : public t_gflt
{
 public:
   t_SRF();
   ~t_SRF();
   
  void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);
};
// ------------------------------------------------------
 
// Square root information filter
// ------------------------------------------------------
class t_SRIF : public t_gflt
{
 public:
   t_SRIF();
   ~t_SRIF();
   
   void update(const Matrix& A, const DiagonalMatrix& P, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Q);
};

} // namespace

#endif

