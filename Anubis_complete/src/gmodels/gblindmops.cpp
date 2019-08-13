
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

 
  Reference:

    RTCA (2006), Minimum operational performance standards for Global Positioning 
                 System/Wide Area Augmentation System airborne equipment,
    Radio Technical Commission for Aeronautics, SC-159, publication DO-229D, Washington, D. C.

-*/

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "gutils/gconst.h"
#include "gmodels/gblindmops.h"

using namespace std;

namespace gnut {   

     	   // MEAN:    LAT      P0      T0     e0   beta0  lambda0
double t_mops::a[30]={ 15, 1013.25, 299.65, 26.31, 0.00630, 2.77,     //    15 or <                            
                       30, 1017.25, 294.15, 21.79, 0.00605, 3.15,     //    30                                 
                       45, 1015.75, 283.15, 11.66, 0.00558, 2.57,     //    45                                 
                       60, 1011.75, 272.15,  6.78, 0.00539, 1.81,     //    60                                 
                       75, 1013.00, 263.65,  4.11, 0.00453, 1.55 };   //    75 or >                            

      // VARIATION:    LAT    P0    T0     e0  beta0  lambda0
double t_mops::b[30]={ 15,  0.00,  0.00, 0.00, 0.0000, 0.00,          //    15 or <                        
                       30, -3.75,  7.00, 8.85, 0.0025, 0.33,          //    30                             
                       45, -2.25, 11.00, 7.24, 0.0032, 0.46,          //    45
                       60, -1.75, 14.00, 5.36, 0.0081, 0.74,          //    60                             
                       75, -0.50, 14.50, 3.39, 0.0062, 0.30 };        //    75 or >                        

// procedures/functions
// --------------------

// ---------
void t_mops::mops(double& lat, double& H, double& elev, double& DOY, double& ZHD, double& ZWD)
{
   
   int n = 5; // number of meteo parameters
   
   Matrix AVG(n,n+1);
   AVG << a;

   Matrix VAR(n,n+1);
   VAR << b;
   
   double DOYmin;
   if(lat >= 0.0) DOYmin = 28.0;
   else           DOYmin = 211.0;

   vector<double> I_AVG, I_VAR;
   vector<double> EVAL;


   _tabval(lat, AVG, VAR, I_AVG, I_VAR);   
   _extrap(DOY, DOYmin, I_AVG, I_VAR, EVAL);
   
#ifdef DEBUG   
   for(int i = 0; i<5; i++) cout << EVAL[i] << endl;
#endif
   
   // tropo delay calculation
   double ZHD0, ZWD0, mf;
   delay(H, elev, EVAL, ZHD0, ZHD, ZWD0, ZWD, mf);
   
#ifdef DEBUG   
   cout << fixed << setprecision(4) << "  DOY:    " << DOY
                                    << "  LAT:    " << lat
                                    << "  H:      " << H 
                                    << "  elev:   " << elev
                                    << "  DOYmin: " << DOYmin
                                    << "  ZHD0:   " << ZHD0 
                                    << "  ZHD:    " << ZHD 
                                    << "  ZWD0:   " << ZWD0 
                                    << "  ZWD:    " << ZWD 
                                    << "  ZTD:    " << ZHD+ZWD 
                                    << "  mf:     " << mf 
                                    << "  RTE:    " << 0.12*mf 
                                    << "  STD:    " << mf*(ZHD+ZWD) << endl;
#endif
}

// Delay calculation   
// ---------
int t_mops::delay( double &H, double &elev, vector<double> &EVAL, double &ZHD0, double &ZHD, double &ZWD0, double &ZWD, double &mf )
{
   /*settings*/
   double k1 = 77.604;  /* [K/mbar]   */
   double k2 = 382000;  /* [K^2/mbar] */
   double rd = 287.054; /* [J/(kg*K)] */
   double gm = 9.784;   /* [m/s^2]    */
   double g  = 9.80665; /* [m/s^2]    */
   ZHD0 = 1e-6*k1*rd*EVAL[0]/gm;
   ZWD0 = (1e-6*k2*rd)/(gm*(EVAL[4]+1)-EVAL[3]*rd) * (EVAL[2]/EVAL[1]);
   ZHD  = pow( (1-(EVAL[3]*H/EVAL[1])), (g/(Rd*EVAL[3])) ) * ZHD0;
   ZWD  = pow( (1-(EVAL[3]*H/EVAL[1])), ( ((EVAL[4]+1)*g)/(rd*EVAL[3]) -1) ) * ZWD0;
   mf   = 1.001 / sqrt(0.002001 + pow(sin(elev),2) );
   
   return 0;
}

   

// tabular reading   
// ---------
int t_mops::_tabval( double &lat, Matrix AVG, Matrix VAR, vector<double> &I_AVG, vector<double> &I_VAR )
{
  double nlat = lat;
  if( lat < 0 ) nlat = abs(lat);
  
  int n = 5;
  if( nlat < 15 ){
    for(int i = 1; i <= n; ++i){
      I_AVG.push_back( AVG(1,i+1) );
      I_VAR.push_back( VAR(1,i+1) );
    }
  }

  else if( nlat > 75 ){
    for(int i = 1; i <= n; ++i){
      I_AVG.push_back( AVG(5,i+1) );
      I_VAR.push_back( VAR(5,i+1) );
    }
  }

  else if( nlat == 15.0 || nlat == 30.0 || nlat == 45.0 || nlat == 60.0 || nlat == 75.0 ){
    for(int i = 1; i <= n; ++i){
      if( nlat == VAR(i,1) ){
        for(int j = 2; j <= n+1; ++j){
          I_AVG.push_back( AVG(i,j) );
          I_VAR.push_back( VAR(i,j) );
        }
      }
    }
  }

  else{
    t_ginterp interp;
    double fvalavg, fvalvar;

    for(int j = 1; j <= n; j++){
      map<double, double> I;
      map<double, double> J;
      for(int i = 1; i < n; i++)
      {
        if(nlat > AVG(i,1) && nlat < AVG(i+1,1)){
          I[AVG(i,  1)] = AVG(i,  j+1);
          I[AVG(i+1,1)] = AVG(i+1,j+1);
          J[VAR(i,  1)] = VAR(i,  j+1);
          J[VAR(i+1,1)] = VAR(i+1,j+1);
	}
      }

      interp.linear(I, nlat, fvalavg);
      interp.linear(J, nlat, fvalvar);

      I_AVG.push_back( fvalavg );
      I_VAR.push_back( fvalvar );

      I.clear();
      J.clear();
    }
  }
   
#ifdef DEBUG   
 for(int i = 0; i<n; ++i) cout << fixed << setprecision(4) << I_AVG[i] << "  :  " << I_VAR[i] << endl;
#endif

 return 0;
}



// ---------
int t_mops::_extrap( double &DOY, double &DOYmin, vector<double> &AVG,
	  	                                  vector<double> &VAR,
		                                  vector<double> &EVAL )
{
   for(unsigned int i = 0; i<5; ++i) EVAL.push_back( AVG[i] - VAR[i]*cos( 2*G_PI*(DOY-DOYmin)/365.25 ) );
   return 0;
}

} // namespace
