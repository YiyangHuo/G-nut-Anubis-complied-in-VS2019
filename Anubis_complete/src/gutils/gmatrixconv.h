
#ifndef GMATRIXCONV_H
#define GMATRIXCONV_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: matrix conversion utilities
  Version: $ Rev: $

  2011-04-26 /PV: created
  2012-04-06 /JD: extracted matrix conversion utilities from old utils.h

-*/

#include <ostream>   // need for NEWMATIO !
#include <vector>

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

using namespace std;

namespace gnut {

// probably will be obsolete using gtriple everywhere !!!!!
// -----
ColumnVector array2colVec(double*);                     // convert array to columnVector

void Matrix_remRC(SymmetricMatrix&, int row, int col);  // remove   r_th row and c_th column in SymMatrix
void Matrix_remRC(DiagonalMatrix&,  int row);  // remove   r_th row and c_th column in DiagonalMatrix
void Matrix_rem(SymmetricMatrix&, vector<int>&);  // remove rows and columns stored in set
void Matrix_addRC(SymmetricMatrix&, int row, int col);  // add zero r_th row and c_th column in SymMatrix
void Matrix_addRC(Matrix&,          int row, int col);  // add zero r_th row and c_th column in Matrix
void Matrix_addRC(DiagonalMatrix&,  int row);           // add zero r_th row and r_th column in DiagMatrix
void Matrix_swap(SymmetricMatrix&,  int a, int b);      // swap a_th and b_th row and column in SymMatrix
void Matrix_cpRC(SymmetricMatrix, SymmetricMatrix&, int r, int c); // copy row and col
void Vector_add(ColumnVector&, int row);           // add zero r_th row in Column Vector

void indexing(const ColumnVector& v, const ColumnVector& v_sorted, vector<int>& index);
   
DiagonalMatrix SR(DiagonalMatrix& D);
Matrix rotX(double Angle);
Matrix rotY(double Angle);
Matrix rotZ(double Angle);

SymmetricMatrix cov2corr(SymmetricMatrix& Q);
} // namespace

#endif