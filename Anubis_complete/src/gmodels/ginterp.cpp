
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
#include <iostream>
#include <iomanip>

#include "gmodels/ginterp.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_ginterp::t_ginterp()
 : _interp_1d(LINEAR),
   _interp_2d(BILINEAR),
   _interp_3d(VER2HOR),
   _interp_ht(INTERPOLATE)
{}

// destructor
// ----------
t_ginterp::~t_ginterp(){}


// ----------
t_ginterp::INTERP_1D t_ginterp::str_to_interp_1d(const string& str)
{
  if(      str.compare("LINEAR") == 0 ) return LINEAR;
  else if( str.compare("SPLINE") == 0 ) return SPLINE;
  else cerr << "# warning - unknown 1D interpolation method: " << str << endl;
	
  return LINEAR;
}

// ----------
t_ginterp::INTERP_2D t_ginterp::str_to_interp_2d(const string& str)
{
  if(      str.compare("BILINEAR") == 0 ) return BILINEAR;
  else cerr << "# warning - unknown 2D interpolation method: " << str << endl;
	
  return BILINEAR;
}

// ----------
t_ginterp::INTERP_3D t_ginterp::str_to_interp_3d(const string& str)
{
  if(      str.compare("HOR2VER") == 0 ) return HOR2VER;
  else if( str.compare("VER2HOR") == 0 ) return VER2HOR;
  else cerr << "# warning - unknown 3D interpolation method: " << str << endl;
	
  return VER2HOR;
}

// ----------
t_ginterp::INTERP_HT t_ginterp::str_to_interp_ht(const string& str)
{
  if(      str.compare("INTERPOLATE") == 0 ) return INTERPOLATE;
  else if( str.compare("SCALE")       == 0 ) return SCALE;
  else cerr << "# warning - unknown HT interpolation method: " << str << endl;
	
  return INTERPOLATE;
}

// ----------
string t_ginterp::interp_1d_to_str(const INTERP_1D& typ)
{
  switch( typ )
  {
   case LINEAR : return "LINEAR"; break;
   case SPLINE : return "SPLINE"; break;
   default     : cerr << "# warning - unknown 1D interpolation method\n";
  }
  return "LINEAR";
}

// ----------
string t_ginterp::interp_2d_to_str(const INTERP_2D& typ)
{
  switch( typ )
  {
   case BILINEAR : return "BILINEAR"; break;
   default     : cerr << "# warning - unknown 2D interpolation method\n";
  }
  return "BILINEAR";
}

// ----------
string t_ginterp::interp_3d_to_str(const INTERP_3D& typ)
{
  switch( typ )
  {
   case HOR2VER : return "HOR2VER"; break;
   case VER2HOR : return "VER2HOR"; break;
   default     : cerr << "# warning - unknown 3D interpolation method\n";
  }
  return "VER2HOR";
}

// ----------
string t_ginterp::interp_ht_to_str(const INTERP_HT& typ)
{
  switch( typ )
  {
   case INTERPOLATE : return "INTERPOLATE"; break;
   case SCALE       : return "SCALE";       break;
   default          : cerr << "# warning - unknown HT interpolation method\n";
  }
  return "INTERPOLATE";
}


// ----------------------------
// Linear Interpolation (double)
// ----------------------------
int t_ginterp::linear(map<double, double>& data, double val, double& fval)
{
   map<double, double>::iterator it2, it1;
   it2 = data.lower_bound(val);

   if(it2!=data.end() && double_eq(it2->first, val) ){
      fval = it2->second;
      return 0;
   }
   if( it2 == data.begin() || it2 == data.end() ){
#ifdef DEBUG      
      cerr << "# warning: extrapolation is not allowed " << it2->first << " " << val << "!\n";
#endif      
      return -1;
   }
   
   it1 = it2;
   it1--;
   
   Matrix x(2,2), y(2,1);
   x << it1->first << 1
     << it2->first << 1;
   y << it1->second 
     << it2->second;
   
   Matrix coef(2,1);
   coef = x.i()*y;
   
   fval = coef(1,1)*val + coef(2,1);

#ifdef DEBUG   
   cerr << "t_ginterp::linear: " << "x = " << val << " fx = " << fval << endl;
#endif
  return 0;
}


// ----------------------------
// Spline Interpolation (double)
// ----------------------------
int t_ginterp::spline(map<double, double>& data, double val, double& fval)
{
   if(data.size() == 1)
     {
	cerr << "Warning: For using of Cubic Spline Interpolation Method, the dimension of input data vector equal to (at least): dimm = 3, is requested " << endl; return -1;
     }
   else if(data.size() == 2)
     {
	t_ginterp CallLin;
	CallLin.linear(data, val, fval);
     }
   else if(data.size() > 2)
     {
//	map<double, double>::const_iterator I = data.upper_bound( val );
	map<double, double>::const_iterator II;

	vector<double> X, Y;
	for(II = data.begin(); II != data.end(); ++II)
	  {
	     X.push_back(II->first);
	     Y.push_back(II->second);
	  }
	
#ifdef DEBUG	
	for(II = data.begin(); II != data.end(); ++II) cout << II->first << "  :  " << II->second << endl;
	
	cout << "interpolovana hodnota " << val << endl;
#endif

	double    e = val;

	
	vector<double>  h, lambda, eta, B, C, D, g, CoefA, CoefB, CoefC, CoefD, AA, BB, FI;

	// vector h: OK
	for(unsigned int i = 0; i < X.size()-1; i++)
	  {
	     h.push_back( X[i+1]-X[i] );
	  }
	// lambda: OK, eta: OK, g
	for(unsigned int i = 0; i < X.size()-2; i++)
	  {
	     lambda.push_back( h[i] / (h[i] + h[i+1]) );
	     eta.push_back( 1.0 - lambda[i] );
	     B.push_back(2);     // (4.6.12.) diag = 2;
	  }
	for(unsigned int j = 1; j < X.size()-2; j++)
	  {
	     C.push_back( lambda[j] );
	  }
	for(unsigned int j = 0; j < X.size()-3; j++)
	  {
	     D.push_back(eta[j]);
	  }
	
	Matrix A(X.size()-2,X.size()-2), M(X.size()-2, X.size()-2); A = 0.0;
	for(unsigned int i = 1; i<=X.size()-2; i++)
	  for(unsigned int j = 1; j<=X.size()-2; j++)
	    {
	       A(i,i)=B[j-1];     // ok
	    }
	for(unsigned int i = 1; i<X.size()-2; i++)
	  for(unsigned int j = 1; j<X.size()-2; j++)
	    {
	       if( j==i)
		 {
		    A(i,j+1)=D[j-1];
		    A(i+1,j)=C[j-1];
		 }
	    }
	
	ColumnVector G(X.size()-2), MM(X.size());
	for(unsigned int i=0; i<X.size()-2; i++) // ok
	  {
	     g.push_back( ( 6.0/(h[i] + h[i+1]) ) * ((( Y[i+2]-Y[i+1])/h[i+1]) - (( Y[i+1]-Y[i])/h[i]) ) );
	  }
	
	for(unsigned int i=1; i<=g.size(); i++)
	  {
	     G(i) = g[i-1];
	  }
	
	// Solving of linear equations
	M = A.i()*G;
	double M0 = 0.0;
	double Mn = 0.0;
   
	MM.row(1) << M0;
	for(unsigned int i = 2; i<X.size(); i++) // ok
	  {
	     MM.row(i) << M[i-2];
	  }
	MM.row(X.size()) << Mn;
	
	// Coeficients a, b, c, d
	for(unsigned int i = 0; i<X.size()-1; i++)  // ok
	  {
	     CoefA.push_back( Y[i] );
	     CoefB.push_back( ((Y[i+1]-Y[i])/h[i]) - (2*MM[i]+MM[i+1])*h[i]/6 );
	     CoefC.push_back( 0.5*MM[i] );
	     CoefD.push_back( (MM[i+1] - MM[i])/(6*h[i]) );
	     
	     // Interpotion
	                                                                                  
	     AA.push_back( ((Y[i+1] - Y[i])/h[i]) - (h[i]/6)*(MM[i+1]-MM[i]) );
	     BB.push_back( Y[i] - MM[i]*(pow(h[i],2)/6) );
	     FI.push_back( MM[i] * (pow((X[i+1]-e),3)/(6*h[i])) + MM[i+1]*((pow((e - X[i]),3))/(6*h[i])) + AA[i]*(e - X[i]) + BB[i] );
	     if(e == X[i])
	       {
		  fval = Y[i];
	       }
	     else if(X[i] < e && e < X[i+1])
	       {
		  fval = FI[i];
	       }
	  }
     }
   else
     {
	cerr << "Warning (ginterp-spline): something wrong with input data!" << endl;
     }
   
   
#ifdef DEBUG   
   cerr << "t_ginterp::cubic spline: " << "x = " << val << " fx = " << fval << endl;
#endif
  return 0;
}


// ----------------------------
// Linear Interpolation (time)
// ----------------------------
int t_ginterp::linear(const map<t_gtime, double>& data, const t_gtime& epo, double& fval)
{
  if(       data.size() == 0 ){ 
    return -1;
     
  }else if( epo == data.begin()->first ){ // data.size() == 1 && ){

    fval = data.begin()->second;
    return 0;
     
  }else if( epo == data.rbegin()->first ){ // data.size() == 1 && ){
      
    fval = data.rbegin()->second;
     return 0;
     
  }else if( data.size() == 2 ){
     map<t_gtime, double>::const_iterator itDAT;
     map<double,double> tmp;

    for( itDAT = data.begin(); itDAT != data.end(); ++itDAT )
    {
      t_gtime dt(itDAT->first);
      double gt = ( dt.mjd()  - epo.mjd() )*86400
	          + dt.sod()  - epo.sod()
	          + dt.dsec() - epo.dsec();               // all in sec

      tmp[ gt ] = itDAT->second;
//      cout << " seconds = " << itDAT->first.str_ymdhms() << " " << gt << endl;
    }

    int irc = linear( tmp, 0.0, fval );
    return irc;

  }else{
    cerr << "# warning: a problem in time linear interpolation! [#node:" 
         << data.size() << "] " << epo.str_ymdhms() << endl ;
    return -1;
  }
   
  return 0;
}

   // ----------------------------
   // Spline Interpolation (time)
   // ----------------------------
   int t_ginterp::spline(const map<t_gtime, double>& data, const
			 t_gtime& epo, double& fval)
     {
	map<t_gtime, double>::const_iterator I;
	map<double, double> tmp;
	for( I = data.begin(); I != data.end(); ++I )
	  {
	     t_gtime dt(I->first);
	     double gt = ( dt.mjd()  - epo.mjd() )*86400
	       + dt.sod()  - epo.sod()
		 + dt.dsec() - epo.dsec();               // all in sec
	     
	     tmp[ gt ] = I->second;
	     
	  }
	
	int irc = spline( tmp, 0.0, fval );
	return irc;
	
	return 0;
     }
   


// ----------------------------
// bilinear interpolation
// 
// 0 .. 11 (bottom-left)         21 *---------* 22
// 1 .. 12 (bottom-right)           |         |
// 2 .. 21 (top-left)               |         |
// 3 .. 22 (top-right)           11 *---------* 12
// ----------------------------
int t_ginterp::bilinear(const map<t_gpair, double>& data, const t_gpair& req_pos, double& fval)
{
   vector<t_gpair> v_pt;
   map<t_gpair,double>::const_iterator itDAT = data.begin();
   for( ; itDAT != data.end(); ++itDAT){ 
     v_pt.push_back( itDAT->first );      
#ifdef DEBUG      
     cout << "# interp:" << fixed << setprecision(0)
	  << " " << itDAT->first[0]
	  << ":" << itDAT->first[1] 
	  << " " << setprecision(3)
	  << " " << data.at(v_pt.at(v_pt.size()-1))
	  << endl;
#endif
   }

   // no interpolation (1 point)
   if( v_pt.size() == 1 ){  fval = data.at(v_pt[0]);
      
   // linear interpolation (2 points)
   }else if( v_pt.size() == 2 ){
     double x = 0.0;
     double a = 0.0;
     double b = 0.0;
      
     if( v_pt[1][0] - v_pt[0][0] == 0 )
     {
       x = req_pos[1];  // lat=const
       a = v_pt[0][1];
       b = v_pt[1][1];	
     }else{
       x = req_pos[0];  // lon=const
       a = v_pt[0][0];
       b = v_pt[1][0];
     }
       
     if( a - b == 0 ){
       cout << "# warning: Bilinear interpolation failed (2 points) \n";
       fval = 0.0;
       return -1;
	 
     }else{
       fval = data.at(v_pt[0]) + (data.at(v_pt[1]) - data.at(v_pt[0])) * (x - a) / (a - b);
     }
      
   // bilinear interpolation (4 points)
   }else if( v_pt.size() == 4 ){
         
#ifdef DEBUG
     cout << "# Linear interpolation (4): " << fixed << setprecision(0)
	  << "  " << v_pt[0][0] << ":" << v_pt[0][1]
	  << "  " << v_pt[1][0] << ":" << v_pt[1][1]
	  << "  " << v_pt[2][0] << ":" << v_pt[2][1]
	  << "  " << v_pt[3][0] << ":" << v_pt[3][1]
	  << " ." << setprecision(3)
	  << "  " << req_pos[0] << "-" << req_pos[1] 
	  << "  " << data.at(v_pt[0])
	  << "  " << data.at(v_pt[1])
	  << "  " << data.at(v_pt[2])
	  << "  " << data.at(v_pt[3])
	  << endl;
#endif

     double a = data.at(v_pt[0]) * (v_pt[3][0] - req_pos[0]) * (v_pt[3][1] - req_pos[1]);
     double b = data.at(v_pt[1]) * (req_pos[0] - v_pt[2][0]) * (v_pt[2][1] - req_pos[1]);
     double c = data.at(v_pt[2]) * (v_pt[1][0] - req_pos[0]) * (req_pos[1] - v_pt[1][1]);
     double d = data.at(v_pt[3]) * (req_pos[0] - v_pt[0][0]) * (req_pos[1] - v_pt[0][1]);
        
     if( (v_pt[2][0]-v_pt[1][0]) * (v_pt[1][1]-v_pt[0][1]) == 0 )
     {
       cout << "# warning: Bilinear interpolation failed (4 points)\n";
       fval = 0.0;
       return -1;
     }else{
       fval = ( a+b+c+d ) / ( (v_pt[2][0]-v_pt[1][0]) * (v_pt[1][1]-v_pt[0][1]) );
     }      
   }else{
     fval = 0.0;
     cout << "# warning: Bilinear interpolation failed (X points ?)\n";
     return -1;
   }

   return 0;
}

   
int t_ginterp::idw(const map<t_gpair, double>& data, const t_gpair& req_pos, double& fval)
{
if(data.size()<3){
  fval = 0.0;
  cout << "warning: Input data's size is not appropriate for IDW spatial interpolation technique!" << endl;
   return -1;
}
// weights & idw
double wght = 0.0;
double vals = 0.0;   
for(map<t_gpair,double>::const_iterator I=data.begin(); I!=data.end(); ++I){
#ifdef DEBUG
   cout << I->first[0] << "  " 
        << req_pos[0]  << "  "
        << I->first[1] << "  " 
        << req_pos[1]  << "  "
        <<  sqrt((pow((I->first.crd(0)-req_pos.crd(0)), 2.0))+(pow((I->first.crd(1)-req_pos.crd(1)), 2.0)))
        << endl;
#endif
   wght += 1.0/( sqrt((pow((I->first.crd(0)-req_pos.crd(0)), 2.0))+(pow((I->first.crd(1)-req_pos.crd(1)), 2.0))));
   vals += I->second/( sqrt((pow((I->first.crd(0)-req_pos.crd(0)), 2.0))+(pow((I->first.crd(1)-req_pos.crd(1)), 2.0))));
}

fval = vals/wght;

return 0;
}
   
} // namespace
