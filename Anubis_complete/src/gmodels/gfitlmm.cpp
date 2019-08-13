
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

#include "gmodels/gfitlmm.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_gfitlmm::t_gfitlmm( t_gfunc* model )
 : t_gfit( model )
{
  _limit    =   0.0;
  _lambda   =  1e-1;   // lambda parameter                          // 1e-1 ok!
  _lim_par  =  1e-3;   // convergence criterium limit for parameter // 1e-3 ok!
  _lim_grd  =  1e-4;   // convergence criterium limit for gradient  // 1e-4 ok!
}


// fitting  multi-dimensional
// ---------
int t_gfitlmm::fit( const ColumnVector&   X, const ColumnVector& Y,
	 	          ColumnVector& par,       ColumnVector& err,
		    int& niter, string& msg)
{
  int nobs = X.size();  // number of observations

  DiagonalMatrix P(nobs);
  for(int i = 1; i<=nobs; ++i) P(i,i) = 1.0;

  return fit(X,Y,P,par,err,niter,msg);
}


// fitting  multi-dimensional
// --------------------------
int t_gfitlmm::fit( const ColumnVector&   X, const ColumnVector& Y, const DiagonalMatrix& P,
	 	          ColumnVector& par,       ColumnVector& err,
		    int& niter, string& msg)
{
   gtrace("t_fitlmm::fit");

   bool updJacob   = true;        // update Jacobi Matrix
   int nobs        =   X.size();  // number of observations
   int npar        = par.size();  // number of fitted parameters
   int iter        = 0;           // iteration
   int stop        = 0;           // stop criterium   
   double limit    = _limit;
   double limit_LM = _limit;
   double lambda   = _lambda;

   Matrix         J    ( nobs, npar );
   Matrix         H    ( npar, npar );
   Matrix         H_lm ( npar, npar );
   ColumnVector   G    ( npar );
   ColumnVector   Con  ( npar );

   ColumnVector   fit  ( npar );
   ColumnVector   fit0 ( npar );  fit0 = par; // fit = fit0 + dfit
   ColumnVector   dfit ( npar );

   ColumnVector AprRes( nobs ), FitRes( nobs );
   ColumnVector AprEst( nobs ), FitEst( nobs );

   while(stop == 0 && iter <= _max_iter)
   {
     ++iter;

     if( updJacob )
     {
   
#ifdef DEBUG	
       for(int i = 1; i <= nobs; i++)
         cerr << " MODEL FIT2 " << i << " "  << X(i) << " " << Y(i) << " " << nobs << endl;
#endif	

       if( _func->value(X, fit0, AprEst) < 0 ){
	 msg = "gfitlmm error - evaluation not possible";
	 return -1;
       }

       AprRes = Y - AprEst;
       limit  = pow(AprRes.norm_Frobenius(),2.0);

       if( _func->deriv(X, fit0, J) < 0 ){
  	 msg = "gfitlmm error - derivation not possible";
	 return -1;
       }

       J = -J;
       H =  J.t() * P * J;
     }
     
     H_lm = H;
     for(int i = 1; i<=npar; ++i) H_lm(i,i) += lambda*H(i,i);

     G    = J.t() * P * AprRes;
     dfit = - H_lm.i() * G;
     fit  = fit0 + dfit;

     if( _func->value(X, fit, FitEst) < 0 ){
       msg = "gfitlmm error - evaluation not possible";
       return -1;
     }
      
     FitRes = Y - FitEst;

     limit_LM = pow(FitRes.norm_Frobenius(),2.0);
	
     if( limit_LM < limit ){
       fit0     = fit;
       limit    = limit_LM;
       lambda   = lambda/10.0;
       updJacob = true;
     }else{
       lambda   = lambda*10.0;
       updJacob = false;
     }
      
     for(int i = 1; i<=npar; ++i) Con(i) = dfit(i) / fit0(i);

     // tests of convergence
     Real mvG = G.MaximumAbsoluteValue();
     Real mvP = Con.MaximumAbsoluteValue();

     if( iter > 2  &&  mvG <= _lim_grd ) stop = 1;
     if( iter > 2  &&  mvP <= _lim_par ) stop = 2;
     if( iter == _max_iter             ) stop = 9;
      
     if( _log && _log->verb() >= 2 ){
       ostringstream os;
       os << fixed << setprecision(6)
	  << "converge itr[" <<     iter << "]"
	         << "  par[" << _lim_par << "]:" << setw(12) << mvP
	         << "  grd[" << _lim_grd << "]:" << setw(12) << mvG;
       _log->comment(4,"gfitlmm",os.str());

       os.str(""); os.clear(); os << fixed << setprecision(6) << "stop:" << setw(3) << stop;
       switch( stop ){
         case 1: os << " itr[" << iter << "]  par[" << _lim_par << "]:" << setw(12) << mvP;
 	         _log->comment(4,"gfitlmm",os.str()); break;
         case 2: os << " itr[" << iter << "]  grd[" << _lim_grd << "]:" << setw(12) << mvG; 
	         _log->comment(4,"gfitlmm",os.str()); break;
         case 9: os << " itr[" << iter << "]  max[" << _max_iter << "]"; 
	         _log->comment(4,"gfitlmm",os.str()); break;
       }
     }
    }

   // SSR: Sum of square of residuals
   double SSR = 0.0;
   for(int i = 1; i <= nobs; i++)
     SSR = SSR + pow(FitRes(i),2);
   
   // RMS of residuals
   double RMS = sqrt(SSR / nobs);

   // Covariance Matrix   
   Matrix Cov(npar, npar);  Cov = H.i();

   // standard error of parameters
   ColumnVector DCovSQRT(npar);
   for(int i = 1; i<=npar; ++i)  DCovSQRT(i) = sqrt(Cov(i,i));

#ifdef DEBUG   
   cerr << fixed << setprecision(3)
        << " FIT: " << fit0(1)
        << " RMS: " << RMS * DCovSQRT(1)
        << " ITR: " << iter
        << " MAX: " << _max_iter
        << endl;
#endif

 par    = fit0;
 err    = RMS * DCovSQRT;
 niter  = iter-1;
 return 0;
}

} // namespace
