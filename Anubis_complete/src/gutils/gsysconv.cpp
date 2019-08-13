
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

#include <cmath>
#include <iostream>
#include <iomanip>

#include "gutils/gconst.h"
#include "gutils/gsysconv.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// Ellipsoidal to XYZ coordinates
// ----------
int ell2xyz(const double* Ell, double* XYZ, bool degrees )
{  
  const double bell = Aell*(1.0-1.0/Finv) ;
  const double e2   = (Aell*Aell-bell*bell)/(Aell*Aell) ;

  double nn;
  double hsl = Ell[2]; // [m] above mean sea level
   
  if( degrees ){
    nn = Aell/sqrt( 1.0 - e2 * sin(Ell[0]*D2R) * sin(Ell[0]*D2R));
    XYZ[0] =        (nn+hsl) * cos(Ell[0]*D2R) * cos(Ell[1]*D2R);
    XYZ[1] =        (nn+hsl) * cos(Ell[0]*D2R) * sin(Ell[1]*D2R);
    XYZ[2] = ((1-e2)*nn+hsl) * sin(Ell[0]*D2R);
  }else{
    nn = Aell/sqrt( 1.0 - e2 * sin(Ell[0]) * sin(Ell[0]));
    XYZ[0] =        (nn+hsl) * cos(Ell[0]) * cos(Ell[1]);
    XYZ[1] =        (nn+hsl) * cos(Ell[0]) * sin(Ell[1]);
    XYZ[2] = ((1-e2)*nn+hsl) * sin(Ell[0]);
  }

  return 1;
}
   

// XYZ to Ellipsoidal coordinates
// ----------
int xyz2ell(const double* XYZ, double* Ell, bool degrees )
{    
  const double bell = Aell*(1.0-1.0/Finv) ;
  const double e2   = (Aell*Aell-bell*bell)/(Aell*Aell) ;
  const double e2c  = (Aell*Aell-bell*bell)/(bell*bell) ;
                            
  double nn, ss, zps, hOld, phiOld, theta, sin3, cos3;
   
  ss    = sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]);   
   
  if( double_eq(ss, 0.0) ){ Ell[0] = -999; Ell[1] = -999; Ell[2] = -999; return 1; }
   
  zps   = XYZ[2]/ss;
  theta = atan( (XYZ[2]*Aell) / (ss*bell) );
  sin3  = sin(theta) * sin(theta) * sin(theta);
  cos3  = cos(theta) * cos(theta) * cos(theta);
   
  // Closed formula
  Ell[0] = atan( (XYZ[2] + e2c * bell * sin3) / (ss - e2 * Aell * cos3) );
  Ell[1] = atan2(XYZ[1],XYZ[0]);
  nn = Aell/sqrt(1.0-e2*sin(Ell[0])*sin(Ell[0]));
  Ell[2] = ss / cos(Ell[0]) - nn;

  const int MAXITER = 100;
  for( int ii = 1; ii <= MAXITER; ii++ ){
    nn     = Aell/sqrt(1.0-e2*sin(Ell[0])*sin(Ell[0]));
    hOld   = Ell[2];
    phiOld = Ell[0];
    Ell[2] = ss/cos(Ell[0])-nn;
    Ell[0] = atan(zps/(1.0-e2*nn/(nn+Ell[2])));

    if( fabs(phiOld-Ell[0]) <= 1.0e-11 && fabs(hOld-Ell[2]) <= 1.0e-5 ){

      // always convert longitude to 0-360
      if( Ell[1] < 0.0 ) Ell[1] += 2 * G_PI;
 
      if( degrees ){ Ell[0] *= R2D; Ell[1] *= R2D; }
       
      return 0;
    }
  }

  return 1;
}


// XYZ to Topocentric/Local (North, East, Up) coordinates
// ----------
int xyz2neu(const double* XYZ, const double* XYZ_Ref, double* neu)
{
  double ele[3];
  xyz2ell(XYZ_Ref, ele, false);
   
  ColumnVector r(3);
   
  r << XYZ[0]-XYZ_Ref[0] << XYZ[1]-XYZ_Ref[1] << XYZ[2]-XYZ_Ref[2];
   
  double sinPhi = sin(ele[0]);
  double cosPhi = cos(ele[0]);
  double sinLam = sin(ele[1]);
  double cosLam = cos(ele[1]);
  
  neu[0] = - sinPhi*cosLam * r(1)
           - sinPhi*sinLam * r(2)
           + cosPhi        * r(3);
  
  neu[1] = - sinLam * r(1)
           + cosLam * r(2);
  
  neu[2] = + cosPhi*cosLam * r(1)
           + cosPhi*sinLam * r(2)
           + sinPhi        * r(3);
 
  return 1;
}


// Topocentric/Local (NEU) to XYZ coordinates
// ----------
int neu2xyz(const double* Ell, const double* neu,     double* xyz)
{
   double sinPhi = sin(Ell[0]);
   double cosPhi = cos(Ell[0]);
   double sinLam = sin(Ell[1]);
   double cosLam = cos(Ell[1]);
   
   xyz[0] = - sinPhi*cosLam * neu[0]
            - sinLam        * neu[1]
            + cosPhi*cosLam * neu[2];
   
   xyz[1] = - sinPhi*sinLam * neu[0]
            + cosLam        * neu[1]
            + cosPhi*sinLam * neu[2];
   
   xyz[2] = + cosPhi        * neu[0]
            + sinPhi        * neu[2];

  return 1;
}

// rotation matrix for rao2xyz   
int rao2xyz_rot(const ColumnVector& pos, const ColumnVector& vel, Matrix& R)
{   
   ColumnVector along  = vel / vel.norm_Frobenius();
   ColumnVector cross  = crossproduct(pos, vel); cross /= cross.norm_Frobenius();
   ColumnVector radial = crossproduct(along, cross);
   
   R.Column(1) = radial;
   R.Column(2) = along;
   R.Column(3) = cross;  
   
   return 1;
}
  
   
// Radial/Along-track/Out-of-plane (Satell-system) to XYZ coordinates
// ---------------
int rao2xyz(const ColumnVector& pos, const ColumnVector& vel,
            const ColumnVector& rao, ColumnVector& xyz)
{          
  ColumnVector along  = vel / vel.norm_Frobenius();
  ColumnVector cross  = crossproduct(pos, vel); cross /= cross.norm_Frobenius();
  ColumnVector radial = crossproduct(along, cross);
                      
  Matrix RR(3,3);
  RR.Column(1) = radial;
  RR.Column(2) = along;
  RR.Column(3) = cross;
  
  xyz = RR * rao;

  Matrix R(3,3); R = 0;
  rao2xyz_rot(pos, vel, R);   
  
  return 1;   
   
/*   
  Matrix R(3,3);
  rao2xyz_rot(pos, vel, R);
                        
  xyz = R * rao;

  return 1;
*/
}
   
//  XYZ coordinates to Radial/Along-track/Out-of-plane (Satell-system)
// ---------------
int xyz2rao(const ColumnVector& pos, const ColumnVector& vel,
                  ColumnVector& xyz,       ColumnVector& rao)
{                                      
  Matrix R(3,3);    
  rao2xyz_rot(pos, vel, R);
   
  rao = R.t() * xyz;

  return 1;
}   

// ELL coordinates to XYZ coordinates
int ell2xyz(const t_gtriple &ell, t_gtriple &xyz, bool degrees)
{
   double ELL[3] = {ell[0], ell[1], ell[2]};   
   double XYZ[3] = {0.0, 0.0, 0.0};
   int irc = ell2xyz(ELL, XYZ, degrees);
   
   if(irc > 0){
      xyz[0] = XYZ[0];
      xyz[1] = XYZ[1];
      xyz[2] = XYZ[2];
      return 1;
   }else return -1;
   
}
   
   
// ok, 1-2 iterations only!
int xyz2ell(const t_gtriple &crd, t_gtriple &ell, bool degrees )
{
    const double bell = Aell*(1.0-1.0/Finv) ;
    const double e2   = (Aell*Aell-bell*bell)/(Aell*Aell) ;
    const double e2c  = (Aell*Aell-bell*bell)/(bell*bell) ;

    double nn, ss, zps, hOld, phiOld, theta, sin3, cos3;

    ss    = sqrt(crd[0]*crd[0]+crd[1]*crd[1]) ;

    if( double_eq(ss, 0.0) ){ ell[0] = -999; ell[1] = -999; ell[2] = -999; return 1; }

    zps   = crd[2]/ss ;
    theta = atan( (crd[2]*Aell) / (ss*bell) );
    sin3  = sin(theta) * sin(theta) * sin(theta);
    cos3  = cos(theta) * cos(theta) * cos(theta);

    ell.set(0, atan( (crd[2] + e2c * bell * sin3) / (ss - e2 * Aell * cos3) ));
    ell.set(1, atan2(crd[1],crd[0])) ;
    nn     = Aell/sqrt(1.0-e2*sin(ell[0])*sin(ell[0])) ;
    ell.set(2, ss / cos(ell[0]) - nn);
#ifdef DEBUG
    cout << "XYZ2ELL: " << 0 << " " << fixed << setprecision(6) << crd[0] << " " << crd[1] << " " << crd[2]
       <<     "  " <<          setprecision(6) << ell[0] << " " << ell[1] << " " << ell[2]
       <<  " nn: " << nn << " ss: " << ss  << endl;
#endif
    const int MAXITER = 100;
    for (int ii = 1; ii <= MAXITER; ii++)
    {
        nn     = Aell/sqrt(1.0-e2*sin(ell[0])*sin(ell[0])) ;
        hOld   = ell[2] ;
        phiOld = ell[0] ;
        ell.set(2, ss/cos(ell[0])-nn);
        ell.set(0, atan(zps/(1.0-e2*nn/(nn+ell[2])))) ;
#ifdef DEBUG
   cout << "XYZ2ELL: " << ii << " " << fixed << setprecision(6) << crd[0] << " " << crd[1] << " " << crd[2]
       <<     "  " <<          setprecision(6) << ell[0] << " " << ell[1] << " " << ell[2]
       <<  " nn: " << nn << " ss: " << ss << " hOld: " << hOld << "  phiOld: " << phiOld << endl;
#endif
        if ( fabs(phiOld-ell[0]) <= 1.0e-11 && fabs(hOld-ell[2]) <= 1.0e-5 )
        {
            // always convert longitude to 0-360
            if( ell[1] < 0.0 ) ell[1] += 2*G_PI;
 
            if( degrees ){ ell[0] *= R2D; ell[1] *= R2D; }

            return 0;
        }
    }
    
    return 1;
}

// ok, only more iterations
int xyz2ell2(t_gtriple &crd, t_gtriple &ell, bool degrees )
{
    const double A  = Aell;
    const double B  = A*(1.0-1.0/Finv);
    const double E2 = (A*A-B*B)/(A*A);

    double P = sqrt(crd[0]*crd[0]+crd[1]*crd[1]);
    if( double_eq(P, 0.0) ){ ell[0] = -999; ell[1] = -999; ell[2] = -999; return 1; }

    ell[0] = atan2( crd[2], P*(1-E2) );
    ell[1] = atan2( crd[1], crd[0]   );

    double V = A/sqrt(1.0-E2*sin(ell[0])*sin(ell[0]));
    ell[2] = P / cos(ell[0]) - V;

#ifdef DEBUG
    cout << "XYZ2ELL: " << 0 << " " << fixed << setprecision(6) << crd[0] << " " << crd[1] << " " << crd[2]
       <<     "  " <<          setprecision(6) << ell[0] << " " << ell[1] << " " << ell[2]
       <<  "  v: " << V << "  p: " << P
       << endl;
#endif
    const int MAXITER = 100;
    for (int ii = 1; ii <= MAXITER; ii++)
    {
        double helOld = ell[2] ;
        double latOld = ell[0] ;

        V      = A/sqrt(1.0 - E2*sin(ell[0])*sin(ell[0]));
        ell[2] = P/cos(ell[0])- V;
        ell[0] = atan2( crd[2]/P , 1.0 - E2*V/(V + ell[2]));
#ifdef DEBUG
      cout << "XYZ2ELL " << ii << " " << fixed << setprecision(6) << crd[0] << " " << crd[1] << " " << crd[2]
       <<     "  " <<          setprecision(6) << ell[0] << " " << ell[1] << " " << ell[2]
       <<  "  v: " << V << "  p: " << P << " helOld: " << helOld << "  latOld: " << latOld*R2D << endl;
#endif
        if ( fabs(latOld-ell[0]) <= 1.0e-11 && fabs(helOld-ell[2]) <= 1.0e-5 )
        {
            // always convert longitude to 0-360
            if( ell[1] < 0.0 ) ell[1] += 2*G_PI;
 
            if( degrees ){ ell[0] *= R2D; ell[1] *= R2D; }

            return 0;
        }
    }
    
    return 1;
}

//void xyz2neu(t_gtriple &crd, t_gtriple &crd_site, t_gtriple &neu)
void xyz2neu(t_gtriple &ell, t_gtriple &xyz, t_gtriple &neu)   
{
    double sinPhi = sin(ell[0]);
    double cosPhi = cos(ell[0]);
    double sinLam = sin(ell[1]);
    double cosLam = cos(ell[1]);

    neu.set(0, - sinPhi*cosLam * xyz[0]
               - sinPhi*sinLam * xyz[1]
               + cosPhi        * xyz[2]);

    neu.set(1, - sinLam * xyz[0]
               + cosLam * xyz[1]);

    neu.set(2, + cosPhi*cosLam * xyz[0]
               + cosPhi*sinLam * xyz[1]
               + sinPhi        * xyz[2]);
}

void neu2xyz(t_gtriple &ell, t_gtriple &neu, t_gtriple &xyz )
{
    double sinPhi = sin(ell[0]);
    double cosPhi = cos(ell[0]);
    double sinLam = sin(ell[1]);
    double cosLam = cos(ell[1]);

    xyz.set(0, - sinPhi*cosLam * neu[0]
               - sinLam        * neu[1]
               + cosPhi*cosLam * neu[2]);

    xyz.set(1, - sinPhi*sinLam * neu[0]
               + cosLam        * neu[1]
               + cosPhi*sinLam * neu[2]);

    xyz.set(2, + cosPhi        * neu[0]
               + sinPhi        * neu[2]);
}

int xyz2neu(t_gtriple& xyz, SymmetricMatrix& Q_xyz, SymmetricMatrix& Q_neu) {
  
  t_gtriple ele(0,0,0);
  xyz2ell(xyz, ele, false);  
  
  double sinPhi = sin(ele[0]);
  double cosPhi = cos(ele[0]);
  double sinLam = sin(ele[1]);
  double cosLam = cos(ele[1]);
   
  int n_par = Q_xyz.Ncols();
  Matrix R(n_par, n_par); R = 0;
  for(int i = 1; i <= n_par; i++){
    R(i,i) = 1; 
  }

  R(1,1) = -sinPhi*cosLam; R(1,2) = -sinPhi*sinLam; R(1,3) = cosPhi;
  R(2,1) = -sinLam;        R(2,2) =  cosLam;        R(2,3) = 0;
  R(3,1) =  cosPhi*cosLam; R(3,2) =  cosPhi*sinLam; R(3,3) = sinPhi;

  Q_neu << R*Q_xyz*R.t();
  
  return 1;
}
   
int ell2ipp(t_gsatdata& satdata, t_gtriple& ell_site, t_gtriple& ell_ipp)
{
  double radius = R_SPHERE;

  double ippE = acos((radius/(radius+450000.0)) * cos(satdata.ele()));
  double eps  = ippE - satdata.ele();
  double ipp_lat = ell_site[0] + eps*cos(satdata.azi());
  if(fabs(ipp_lat*R2D) > 80) return -1;
  double ipp_lon = ell_site[1] + eps*sin(satdata.azi())/cos(ipp_lat);
  if(ipp_lon < 0) ipp_lon += 2.0*G_PI;
  if(ipp_lon > 2.0*G_PI) ipp_lon -= 2.0*G_PI;

  ell_ipp[0] = ipp_lat;
  ell_ipp[1] = ipp_lon;

  // GPStk formula
  double p = G_PI/2.0 - satdata.ele() - asin((radius*cos(satdata.ele())) / (radius+450000.0));
  t_gtriple ell_ipp_g;
  ell_ipp_g[0] = asin(sin(ell_site[0])*cos(p) + cos(ell_site[0])*sin(p)*cos(satdata.azi()));
  ell_ipp_g[1] = ell_site[1] + asin(sin(p)*sin(satdata.azi())/cos(ell_site[1]));

//  ell_ipp = ell_ipp_g;
  // -------------
  
  return 1;
}

} // namespace
