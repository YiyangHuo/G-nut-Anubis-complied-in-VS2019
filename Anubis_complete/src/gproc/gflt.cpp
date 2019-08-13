
/* ----------------------------------------------------------------------
* G-Nut - GNSS software development library
* 
(c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
Research Institute of Geodesy, Topography and Cartography
Ondrejov 244, 251 65, Czech Republic

This file is part of the G-Nut C++ library.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <chrono>


#include "gproc/gflt.h"
#include "gall/gallfltmat.h"
#include "gutils/gmatrixconv.h"

using namespace std;
using namespace std::chrono;

namespace gnut {


// Constructors
t_gflt::t_gflt()
{
}

t_kalman::t_kalman()
{   
}

t_SRF::t_SRF()
{
}
   
t_SRIF::t_SRIF()
{
}
   
// Destructors
t_gflt::~t_gflt()
{
}

t_kalman::~t_kalman()
{   
}

t_SRF::~t_SRF()
{
}
   
t_SRIF::~t_SRIF()
{
}
   
// Methods
// 
// Update filter
// ---------------
// Input arguments:
// A  ... first design matrix
// Ql ... variance-covariance matrix of measurements
// l  ... reduced measurements
// dx ... state vector
// Qx ... variance-covariance matrix of state
// 
// Output arguments:
// dx ... state vector
// Qx ... variance-covariance matrix of state
// 
void t_kalman::update(const Matrix& A, const DiagonalMatrix& Pl, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Qx)
{
  gtrace("t_kalman::update");
//  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  Matrix K;            // Kalman gain
  Matrix NN;
  
  NN = Pl.i() + A * Qx * A.t(); 

  K = Qx * A.t() * NN.i();
  
//  cout  << "Qx" << endl << Qx << endl << "KAQx" << endl << K*A*Qx << endl << "K" << K << endl;
  
  dx = K * (l-A*dx);               // update state vector
  Qx << Qx - K * A * Qx;               // update variance-covariance matrix of state
//  cout << "Qx_upd" << endl << Qx << endl;
//  high_resolution_clock::time_point t2 = high_resolution_clock::now();
//  auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//  cout << duration << endl;
}

void t_SRF::update(const Matrix& A, const DiagonalMatrix& Pl, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Qx)
{
  gtrace("t_SRF::update");   

//  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  int nObs = A.Nrows();
  int nPar = A.Ncols();
  
  UpperTriangularMatrix SS;
      
  try{
     SS = Cholesky(Qx).t();
  }
  catch (NPDException) {
     return;
  }
  
  Matrix SA = SS*A.t();
  Matrix SRF(nObs+nPar, nObs+nPar); SRF = 0;
  for (int ii = 1; ii <= nObs; ++ii) {       
    SRF(ii,ii) = 1.0 / sqrt(Pl(ii,ii));
  }
  
  SRF.SubMatrix   (nObs+1, nObs+nPar, 1, nObs) = SA;
  SRF.SymSubMatrix(nObs+1, nObs+nPar)          = SS;
  
  UpperTriangularMatrix UU;
  QRZ(SRF, UU);
  
  SS = UU.SymSubMatrix(nObs+1, nObs+nPar);
  UpperTriangularMatrix SH_rt = UU.SymSubMatrix(1, nObs);
  Matrix YY  = UU.SubMatrix(1, nObs, nObs+1, nObs+nPar);
  
  UpperTriangularMatrix SHi = SH_rt.i();
  
  Matrix KT  = SHi * YY;
  SymmetricMatrix Hi; Hi << SHi * SHi.t();
  
  dx = KT.t() * l;
  Qx << (SS.t() * SS);  
  
//  high_resolution_clock::time_point t2 = high_resolution_clock::now();
//  auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//  cout << duration << endl; 
  
}

// needs debugging
void t_SRIF::update(const Matrix& A, const DiagonalMatrix& Pl, const ColumnVector& l, ColumnVector& dx, SymmetricMatrix& Qx)
{
  int nObs = A.Nrows();
  int nPar = A.Ncols();
  SymmetricMatrix INF(nPar); INF = 0;
  INF = Qx.i();
  UpperTriangularMatrix R = Cholesky(INF).t();
  
  ColumnVector z(nPar); z=0;
  
  Matrix SRIF(nObs+nPar, nPar+1); SRIF=0;

  SRIF.SymSubMatrix(1,nPar) = R;
  SRIF.SubMatrix(nPar+1, nPar+nObs, 1, nPar) = A;
  SRIF.SubMatrix(1, nPar, nPar+1, nPar+1) = z;
  SRIF.SubMatrix(nPar+1, nPar+nObs, nPar+1, nPar+1) = l;
  
  /* SRIF.SubMatrix(1, nObs, 1, nPar) = A; */
  /* SRIF.SubMatrix(nObs+1, nPar+nObs, 1, nPar) = R;    */
  /* SRIF.SubMatrix(1, nObs, nPar+1, nPar+1) = l; */
  /* SRIF.SubMatrix(nObs+1, nPar+nObs, nPar+1, nPar+1) = z; */
  
  UpperTriangularMatrix MU;
  QRZ(SRIF, MU);
  
  R = 0; z = 0;
  R = MU.SymSubMatrix(1,nPar);
  z = MU.SubMatrix(1, nPar, nPar+1, nPar+1);
  
  INF << R * R.t();
  
  dx = R.i() * z;
  Qx = INF.i();
  //   cout << "nPar = " << nPar << endl << "nObs = " << nObs << endl;
  //   cout << "MU dim = (" << MU.Nrows() << ", " << MU.Ncols() << " )" << endl;
  //  cout << Qx << endl;
  //  cout << dx << endl;
  //  int ooo; cin >> ooo;       
}
  
} // namespace
