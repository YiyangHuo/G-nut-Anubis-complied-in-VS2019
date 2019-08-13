
#ifndef GFUNCLIN_H
#define GFUNCLIN_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements different functions (and derivatives) for 2D data fitting
  Version: $ Rev: $

  2014-03-21 /JD: created

-*/

#include <map>

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

#include "gmodels/gfunc.h"

using namespace std;

namespace gnut {   


// y = fit(0) - fit(1) * x
// ----------------------------
class t_gfunc_linH0 : public t_gfunc
{   
 public:
            t_gfunc_linH0( const map<string, double>& cf ) : t_gfunc(cf){};
   virtual ~t_gfunc_linH0(){};

   virtual int value( const ColumnVector& dat,
		      const ColumnVector& fit, ColumnVector& val );
   virtual int deriv( const ColumnVector& dat,  
		      const ColumnVector& fit, Matrix& val );
};


// y = fit(2) - fit(1) * (x-H0)
// ----------------------------
class t_gfunc_linHS : public t_gfunc
{   
 public:
            t_gfunc_linHS( const map<string, double>& cf ) : t_gfunc(cf){};
   virtual ~t_gfunc_linHS(){};

   virtual int value( const ColumnVector& dat,
		      const ColumnVector& fit, ColumnVector& val );
   virtual int deriv( const ColumnVector& dat,  
		      const ColumnVector& fit, Matrix& val );
};


// y = a0 - fit*(x - H0)
// ---------------------
class t_gfunc_linH1 : public t_gfunc
{   
 public:
            t_gfunc_linH1( const map<string, double>& cf ) : t_gfunc(cf){};
   virtual ~t_gfunc_linH1(){};

   virtual int value( const ColumnVector& dat,
		      const ColumnVector& fit, ColumnVector& val );
   virtual int deriv( const ColumnVector& dat,  
		      const ColumnVector& fit, Matrix& val );
};

// y = fit(0) + x*fit(1)
// ---------------------
class t_gfunc_linT : public t_gfunc
{   
 public:
            t_gfunc_linT( const map<string, double>& cf ) : t_gfunc(cf){};
   virtual ~t_gfunc_linT(){};

   virtual int value( const ColumnVector& dat,
		      const ColumnVector& fit, ColumnVector& val );
   virtual int deriv( const ColumnVector& dat,  
		      const ColumnVector& fit, Matrix& val );
};

   
// y = fit1*cos(x) + fit2*sin(x)
// -----------------------------
class t_gfunc_linCS : public t_gfunc
{   
 public:
            t_gfunc_linCS( const map<string, double>& cf ) : t_gfunc(cf){};
   virtual ~t_gfunc_linCS(){};

   virtual int value( const ColumnVector& dat,
		      const ColumnVector& fit, ColumnVector& val );
   virtual int deriv( const ColumnVector& dat,  
		      const ColumnVector& fit, Matrix& val );
};
   

// y = [ fit1*cos(x) + fit2*sin(x) ] / H_tropo [km]
// ------------------------------------------------
class t_gfunc_linCSH : public t_gfunc
{   
 public:
            t_gfunc_linCSH( const map<string, double>& cf ) : t_gfunc(cf){};
   virtual ~t_gfunc_linCSH(){};

   virtual int value( const ColumnVector& dat,
		      const ColumnVector& fit, ColumnVector& val );
   virtual int deriv( const ColumnVector& dat,  
		      const ColumnVector& fit, Matrix& val );
};

} // namespace

#endif