
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

#include "gmodels/gfuncplane.h"

using namespace std;

namespace gnut {   
   // LSQ (plane)
   // ---------
   //   int t_gfunc_plane::LSQ_plane(const ColumnVector& X, const ColumnVector& Y, const ColumnVector& Z )                             // !!! MOZNO NESKOR MOZEM FUNKCIU PRETAZIT A NA VSTUPE VOLAT TRI VEKTORY DAT
   int t_gfunc_plane::LSQ_plane(map<t_gpair, double>& DATA, map<t_gpair, double>& OUT, ColumnVector& fit)
     {
	int n = DATA.size();
	
	vector<double> vXorig, vYorig;
	vector<double> vX; ColumnVector X(n);                        /*for example: X = lat */
	vector<double> vY; ColumnVector Y(n);                        /*for example: Y = lon */
	vector<double> vZ; ColumnVector Z(n);                        /*for example: Z = ZTD */

	for(map<t_gpair, double>::iterator I = DATA.begin(); I!= DATA.end(); ++I)
	  {

	     vXorig.push_back(I->first.crd(0));
	     vYorig.push_back(I->first.crd(1));

	     vX.push_back(I->first.crd(0));
	     vY.push_back(I->first.crd(1));
//	     vX.push_back(6378*sin(I->first.crd(0)/(180/G_PI)));
//	     vY.push_back(6378*sin(I->first.crd(1)/(180/G_PI)));
//	     vX.push_back(2*6378*sin( (I->first.crd(0)/2)/(180/G_PI) ));
//	     vY.push_back(2*6378*sin( (I->first.crd(1)/2)/(180/G_PI) ));
	     vZ.push_back(I->second);
	  }
	
	for( int i = 1; i<=n; ++i)
	  {
	     X(i) = vX[i-1];
	     Y(i) = vY[i-1];
	     Z(i) = vZ[i-1];
	  }
#ifdef DEBUG
	for(int i = 0; i <n; i++)
	  {
	     cout << vX[i] << "  " << vY[i] << "  " << vZ[i] << endl; 
 	  }
#endif
	if(X.size()==0 || Y.size() == 0 || Z.size() == 0){ cerr << "Warning: Some of the input vectors is empty." << endl; return -1;}
	
	if(X.size() != Y.size() || X.size() != Z.size() || Y.size() != Z.size()){ cerr << "Warning: Input vector sizes are supposet to be equal" << endl; return -1;}
	
	// Settings
	 
	double eps = 1e-6;      // convergence criterium for gradient
	
	int Npar   = 3;         // number of fitted parameters
	int Niter  = 100;       // max iterations
	
	int it     = 0;         // ini. iteration
	int stop   = 0;         // stop criterium
	
	// Procedure
//	ColumnVector fit(Npar);
	ColumnVector rms(Npar);
//	for(int i = 1; i <=Npar; ++i){fit(i)=0.0;} // LEPE INICIOVAT VNE - rozumnymi hodnotami! FUNGUJE ALE DOBRE
     // cout << "fit " << fit << endl;
	ColumnVector val(Z.size());
	
//	double err, err_LM;
	
	Matrix A(X.size(), Npar);            // design matrix
	Matrix M(Npar, Npar);                // covariance matrix
	ColumnVector N(Npar);                // right side vector
	ColumnVector dp(Npar);               // fitted parameters   
	ColumnVector Con(Npar);

	ColumnVector DIF(n);
	
	while(stop == 0 && it <= Niter)
	  {
	     it++;
	     value(X, Y, fit, val);
	     
	     for(int i = 1; i<=n; ++i) DIF(i) = Z(i)-val(i);

	     //err    = norm_sq(Differ);
	     
	     // A: Designe matrix
	     deriv(X, Y, A);

	     M  = A.t() * A;
	     N  = A.t()*DIF;
	     dp = M.i()*N;
	
	     for(int j = 1; j<=fit.size(); ++j){fit(j) = fit(j)+dp(j);}
//cout << "fit DP "  << dp  << endl;
//cout << "fit A  "  << A << endl;
	     
	     // end of procedure
		      
	     // tests of convergence
	     Real mv = N.MaximumAbsoluteValue();
	     if(mv <= eps && it > 2)
	       {
#ifdef DEBUG
		  cout << " Z = aX + bY + c                                        "           << endl;
		  cout << "----------------                                        "           << endl;
		  cout << fixed << setprecision(6) << "** Convergence in AtdY   ** "           << endl;
		  cout << fixed << setprecision(6) << "** Iterations            ** " << it     << endl;
		  cout << fixed << setprecision(6) << "** Fitted parameters (a) ** " << fit(1) << endl;
		  cout << fixed << setprecision(6) << "                     (b)    " << fit(2) << endl;
		  cout << fixed << setprecision(6) << "                     (c)    " << fit(3) << endl;
		  cout << fixed << setprecision(6) << "** Epsilon               ** " << eps    << endl;
#endif
		  stop = 1;
	       }
	     
	     if(it == Niter)
	       {
#ifdef DEBUG
		  cout << "** Maximum number of iterations reached without convergence! ** " << endl;
		  cout << fixed << setprecision(6) << "** Fitted parameters (a) ** " << fit(1) << endl;
		  cout << fixed << setprecision(6) << "                     (b)    " << fit(2) << endl;
		  cout << fixed << setprecision(6) << "                     (c)    " << fit(3) << endl;
#endif
		  stop = 1;
	       }
	  }

	// output
	vector<double>out;
	for( int k = 1; k <=n; ++k){out.push_back(val(k));}
	for( int l = 0; l <n; ++l)
	  {
	     t_gpair temp00(vXorig[l], vYorig[l]);	
	     OUT[temp00]=out[l];
	  }
	
	return 0;
     }
   
   
   // function (plane)
   // ---------
   int t_gfunc_plane::value( const ColumnVector& X, const ColumnVector& Y, const ColumnVector& fit, ColumnVector& val )
     {
	/*
	 if( fit.size() != 2 ){
	 cout << "# warning: not properly initialized function (plane)\n"; return -1;
	 }
	 */
	for( int i=1; i<=X.size(); ++i ){
	   val(i) = fit(1)*X(i) + fit(2)*Y(i) + fit(3);
	}
	return 0;
     }
   // ---------
   int t_gfunc_plane::deriv( const ColumnVector& X, const ColumnVector& Y, Matrix& val )
     {
	/*
	 if( fit.size() != 2 ){
	 cout << "# warning: not properly initialized function (plane)\n"; return -1;
	 }
	 */ 
	
	for( int i=1; i<=X.size(); ++i ){
	   val(i,1) = X(i);
	   val(i,2) = Y(i);
	   val(i,3) = 1.0;
	}
	return 0;
     }
   
} // namespace
