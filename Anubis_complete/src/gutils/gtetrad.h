
#ifndef GTETRAD_H
#define GTETRAD_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements 4D representation (e.g. coordinates + time)
  Version: $ Rev: $

  2012-05-11 /JD: created

-*/

#include <iostream>
#include <string.h>

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

using namespace std;

namespace gnut {

class t_gtetrad {

 public:
  t_gtetrad();
  t_gtetrad(double x, double y, double z, double t);
  t_gtetrad(double crd[4]);
  t_gtetrad(const ColumnVector& crd);
  virtual ~t_gtetrad();
  
  t_gtetrad& operator=(const t_gtetrad& other);     // assignment operator
  t_gtetrad  operator+(const t_gtetrad& other) const;     //
  bool       operator==(const t_gtetrad& tr) const; // equal operator
  double&    operator[](const size_t idx);          // get a reference of element
  double     operator[](const size_t idx) const;    // get value of element
  friend ostream& operator<<(ostream& os, const t_gtetrad& x);
     
  double        crd(int idx) const;                 // get single element
  void          set(int idx, double newValue);      // set single element
  void          set(const ColumnVector&);           // set array by ColumnVector
  void          set(double crd[4]);                 // set array by array
  double*       crd_array();                        // get array
  ColumnVector  crd_cvect();                        // get ColumnVector
  t_gtetrad&    crd_tetrad();                       // get tetrad
  ColumnVector  unitary();                          // get unit ColumnVector
  bool          zero();                             // true: zero elements, false: not zero elements

 protected:

 private:
  double         _crd[4];

};

} // namespace

#endif