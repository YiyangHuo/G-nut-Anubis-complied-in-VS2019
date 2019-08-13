
#ifndef GFUNCPLANE_H
#define GFUNCPLANE_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements different functions (and derivatives) for 2D data fitting
  Version: $ Rev: $

  2014-03-21 /JD: created

-*/

#include <map>
#include <vector>

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

#include "gutils/gpair.h"

using namespace std;

namespace gnut {   

/*
// ax + by + cz + d = 0 
// --------------------
class t_gfunc_plane : public t_gfunc
{   
 public:
            t_gfunc_plane( const map<string, double>& cf ) : t_gfunc(cf){};
   virtual ~t_gfunc_plane(){};

   virtual int value( const ColumnVector& dat,
		      const ColumnVector& fit, ColumnVector& val );
   virtual int deriv( const ColumnVector& dat,  
		      const ColumnVector& fit, Matrix& val );
};
*/
   // ax + by + cz + d = 0 
   // --------------------
   class t_gfunc_plane
     {   
      public:
	t_gfunc_plane(){};
	virtual ~t_gfunc_plane(){};

//	virtual	int LSQ_plane(const ColumnVector& X, 
//			      const ColumnVector& Y, 
//			      const ColumnVector& Z );
//			      
	virtual int LSQ_plane(map<t_gpair, double>& DATA, map<t_gpair, double>& OUT, ColumnVector& fit);
	
	virtual int value( const ColumnVector& Xdat, 
			   const ColumnVector& Ydat,
			   const ColumnVector& fit, ColumnVector& val );
	
	virtual int deriv( const ColumnVector& Xdat,  
			   const ColumnVector& Ydat,  
			   Matrix& val );
};
   
} // namespace

#endif