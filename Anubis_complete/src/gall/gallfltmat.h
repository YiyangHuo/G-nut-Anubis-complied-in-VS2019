
#ifndef GALLFLTMAT_H
#define GALLFLTMAT_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $ Rev: $

  2013-08-14 /PV: created

-*/

#include <vector>
#include <set>

#ifdef BMUTEX
  #include <boost/thread/mutex.hpp>
#endif

#include "gdata/gdata.h"
#include "gproc/gfltmat.h"

using namespace std;

namespace gnut {  

// ----------
class t_gallfltmat : public t_gdata {
   
 public:
  t_gallfltmat();
  virtual ~t_gallfltmat();
  vector<t_gfltmat> data;

  set<int> addAmb;
   
  void add(const t_gfltmat& mat);
  int common(t_gallpar& par1, SymmetricMatrix& Q1, t_gallpar& par2, SymmetricMatrix& Q2, DiagonalMatrix* Noise = NULL);
  int syncSMT(t_gallpar& Xu, SymmetricMatrix& Qu, t_gallpar& Xsm, SymmetricMatrix& Qsm, DiagonalMatrix* Noise = NULL);
  int checkSlp(t_gallpar& Xu, SymmetricMatrix& Qu, t_gallpar& Xp, SymmetricMatrix& Qp, t_gallpar& Xsm, SymmetricMatrix& Qsm, set<string>& slips);
  int checkOutl(t_gallpar& Xp, SymmetricMatrix& Qp, t_gallpar& Xsm, SymmetricMatrix& Qu);
  void clear();
  void clear_addAmb();   
   
  // Not necessary
  int reorder(t_gallpar& X1, SymmetricMatrix& Q1, t_gallpar& X2, SymmetricMatrix& Q2);
  int reorder2(t_gallpar& X1, SymmetricMatrix& Q1, t_gallpar& X2, SymmetricMatrix& Q2);
   
 protected:
};

} // namespace

#endif
