
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

#include <iostream>
#include <iomanip>
#include <cmath>

#include "gutils/gmatrixconv.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// Convert array to column vector
// ----------
ColumnVector array2colVec(double arr[3])
{
//int dim = sizeof(arr)/sizeof(arr[0]);
  ColumnVector vec(3);
  for( int i=0; i<=2; i++ ){
    vec(i+1) = arr[i];
  }
  return vec;
}


// remove r_th row and c_th column in SymMatrix
// ----------
void Matrix_remRC(SymmetricMatrix& Q, int row, int col)
{   
  SymmetricMatrix Qt(Q.Nrows());
  Qt = Q;
  Q.ReSize(Q.Nrows()-1);
  int rr = 1;
    for( int r = 1; r <= Qt.Nrows(); r++ ){
		   
      if( r == row ) continue;
      int cc = 1;
	       
      for( int c = 1; c <= Qt.Ncols(); c++ ){
	 
        if( c == col ) continue;
        Q(rr, cc) = Qt(r,c);
	cc++;
      }	
    rr++;
  }
}

void Matrix_remRC(DiagonalMatrix& D, int row)
{
  SymmetricMatrix Dt(D.Nrows());
  Dt = D;
  D.ReSize(D.Nrows()-1);
  int rr = 1;
    for( int r=1; r<=Dt.Nrows(); r++ ){
	if (r == row) continue;
	D(rr, rr) = Dt(r,r);
        rr++;  
    }	
}


void Matrix_rem(SymmetricMatrix& Q, vector<int>& ind)
{
  vector<int>::iterator it;
  vector<int>::iterator it2;   
  for (it = ind.begin(); it != ind.end(); it++){
     Matrix_remRC(Q, *it, *it);
     for (it2 = it; it2 != ind.end(); it2++) (*it2)--;
  }
}


// add zero r_th row and c_th column in SymMatrix
// ----------
void Matrix_addRC(SymmetricMatrix& Q, int row, int col)
{
   SymmetricMatrix Qt(Q.Nrows());
   Qt = Q;   
   Q.ReSize(Q.Nrows()+1); Q = 0.0;
   
  bool bc = false;
  bool br = false;
  int  rr = 1;
   
  for( int r=1; r<=Q.Nrows(); r++ ){
    if( r == row ){
       br = true; 
       continue;
    }
       
    int cc = 1;
    for( int c=1; c<=Q.Ncols(); c++ ){
      if( c == col ){
        bc = true; 
        continue;
      }
       
      if(  br &&  bc ) Q(r,c)  = Qt(rr,cc);
      if( !br &&  bc ) Q(r,c)  = Qt(rr,cc);
      if(  br && !bc ) Q(r,c)  = Qt(rr,cc);
      if( !br && !bc ) Q(r,c)  = Qt(rr,cc);
	
      cc++;
    }
    rr++;
  }  
}


// add zero r_th row and c_th column in Matrix
// ----------
void Matrix_addRC(Matrix& Q, int row, int col)
{
  Matrix Qt(Q.Nrows(), Q.Ncols());
  Qt = Q;
  if(       row == 0 ){ Q.ReSize(Q.Nrows()  , Q.Ncols()+1); Q = 0;
  }else if( col == 0 ){ Q.ReSize(Q.Nrows()+1, Q.Ncols()  ); Q = 0;
  }else{                Q.ReSize(Q.Nrows()+1, Q.Ncols()+1); Q = 0;
  }

  
  bool bc = false;
  bool br = false;
  int  rr = 1;
  for( int r=1; r<=Q.Nrows(); r++ ){
    if( r == row ){
      br = true; 
      continue;
    }
      
    int cc = 1;
    for( int c=1; c<=Q.Ncols(); c++ ){
      if( c == col ){
        bc = true; 
        continue;
      }
	  
      if(  br &&  bc ) Q(r,c)  = Qt(rr,cc);
      if( !br &&  bc ) Q(r,c)  = Qt(rr,cc);
      if(  br && !bc ) Q(r,c)  = Qt(rr,cc);
      if( !br && !bc ) Q(r,c)  = Qt(rr,cc);
      cc++;
    }       
  rr++;
  }
}


// add zero r_th row and r_th column in DiagMatrix
// ----------
void Matrix_addRC(DiagonalMatrix& Q, int row){
   
  DiagonalMatrix Qt(Q.Nrows());
  Qt = Q;
  Q.ReSize(Q.Nrows()+1); Q = 0;
  
  int rr = 1;
  for( int r=1; r<=Q.Nrows(); r++){
    if( r == row) continue;
    else{
      Q(r,r) = Qt(rr,rr); 
      rr++;
    }      
  } 
}


// add zero r_th row and c_th column in SymMatrix
// ----------
void Vector_add(ColumnVector& V, int row)
{   
  ColumnVector Vt(V.Nrows());
  Vt = V;
  V.ReSize(V.Nrows()+1); V = 0;
   
  int rr = 1;
  for( int r=1; r<=V.Nrows(); r++ ){
    if( r == row )  continue;
    else{     
      V(r) = Vt(rr);
      rr++;
    }	
  } 
}


// Rotation Matrix
// -----------------------------------
Matrix rotX(double Angle) {
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix UU(3,3);
    UU[0][0] = 1.0;  UU[0][1] = 0.0;  UU[0][2] = 0.0;
    UU[1][0] = 0.0;  UU[1][1] =  +C;  UU[1][2] =  +S;
    UU[2][0] = 0.0;  UU[2][1] =  -S;  UU[2][2] =  +C;
    return UU;
}
  
Matrix rotY(double Angle) {
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix UU(3,3);
    UU[0][0] =  +C;  UU[0][1] = 0.0;  UU[0][2] =  -S;
    UU[1][0] = 0.0;  UU[1][1] = 1.0;  UU[1][2] = 0.0;
    UU[2][0] =  +S;  UU[2][1] = 0.0;  UU[2][2] =  +C;
    return UU;
}
  
Matrix rotZ(double Angle) {
    const double C = cos(Angle);
    const double S = sin(Angle);
    Matrix UU(3,3);
    UU[0][0] =  +C;  UU[0][1] =  +S;  UU[0][2] = 0.0;
    UU[1][0] =  -S;  UU[1][1] =  +C;  UU[1][2] = 0.0;
    UU[2][0] = 0.0;  UU[2][1] = 0.0;  UU[2][2] = 1.0;
    return UU;
}

// swap a_th and b_th row and column in SymMatrix
// -------------
void Matrix_swap(SymmetricMatrix& Q, int a, int b)
{
   // Copy symmetric matrix to normal matrix
   Matrix T(Q.Nrows(), Q.Ncols());
   for (int r = 1; r <= Q.Nrows(); r++){
      for (int c = 1; c <= Q.Ncols(); c++){
	 T(r,c) = Q(r,c);
      }      
   }
      
   // Swap rows   
   for (int i = 1; i <= Q.Ncols(); i++){
      double temp = T(a,i);      
      T(a,i) = T(b,i);
      T(b,i) = temp;
   }
     
   // Swap columns
   for (int i = 1; i <= Q.Nrows(); i++){
      double temp = T(i,a);      
      T(i,a) = T(i,b);
      T(i,b) = temp;
   }

   // Back to symmetric matrix
   for (int r = 1; r <= T.Nrows(); r++){
      for (int c = 1; c <= T.Ncols(); c++){
	 Q(r,c) = T(r,c);
      }      
   }   
}

// Copy row and col
// --------------
void Matrix_cpRC(SymmetricMatrix Q1, SymmetricMatrix& Q2, int r, int c)
{
   if (Q1.Nrows() != Q2.Nrows()){
      cerr << "Matrix_cpRC: not compatible dimension" << endl; 
   }
   
   for (int i = 1; i <= Q1.Nrows(); i++){
     for (int j = 1; j <= Q1.Ncols(); j++){
	if (i == r || j == c) Q2(i,j) = Q1(i,j);
     }
   }
   
}

// sqrt of diagonal matrix
// ---------------
DiagonalMatrix SR(DiagonalMatrix& D)
{
   DiagonalMatrix SR( D.Nrows() );
   
   for (int i = 1; i <= D.Nrows(); i++){
      SR(i,i) = sqrt( D(i,i) );
   } 
   
   return SR;
}

// find indexes of sorted elements in a original vector
// ---------------------   
void indexing(const ColumnVector& v, const ColumnVector& v_sorted, vector<int>& index)
{     
   for(int i = 1; i <= v_sorted.Nrows(); i++){
      for(int j = 1; j <= v.Nrows(); j++){
	 if( double_eq(v_sorted(i),v(j)) ) {
	    index.push_back(j);
	    break;
	 }	 
      }      
   }
}

// Variance-covarriance matrix to correlation matirx
// -----------------------------
SymmetricMatrix cov2corr(SymmetricMatrix& Q)
{
  SymmetricMatrix C(Q.Nrows());
  
  for(int i = 1; i <= Q.Nrows(); i++){
    for(int j = 1; j <= i; j++){    
      if(i == j) C(i,j) = 1;
      C(i,j) = Q(i,j) / ( sqrt(Q(i,i)) * sqrt(Q(j,j)) );
    }
  }
  
  return C;
}
   
   
} // namespace
