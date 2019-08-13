
#ifndef GTRIPLE_H
#define GTRIPLE_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements vector of 3D (e.g. coordinates)
  Version: $ Rev: $

  2012-05-11 /JD: created

-*/

#include <iostream>
#include <string.h>

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

#include "gutils/gpair.h"

using namespace std;

namespace gnut {

class t_gtriple {

 public:
  t_gtriple();
  t_gtriple(double x, double y, double z);
  t_gtriple(double crd[3]);
  t_gtriple(const ColumnVector& crd);
  virtual ~t_gtriple();
  
  t_gtriple& operator=(const t_gtriple& other);     // assignment operator
  t_gtriple& operator+=(const t_gtriple& other);
  t_gtriple& operator-=(const t_gtriple& other);
  t_gtriple& operator*=(const double& x);
  t_gtriple& operator/=(const double& x);

  t_gtriple  operator+(const t_gtriple& other) const;
  t_gtriple  operator-(const t_gtriple& other) const;
  t_gtriple  operator*(const double& x) const;
  t_gtriple  operator/(const double& x) const;
   
  bool       operator==(const t_gtriple& tr) const; // equal operator
  bool       operator!=(const t_gtriple& tr) const; // not equal operator
  bool       operator<(const t_gtriple& tr)  const; // equal for sorting
  double&    operator[](const size_t idx);          // get a reference of element
  double     operator[](const size_t idx) const;    // get value of element
  friend ostream& operator<<(ostream& os, const t_gtriple& x);
     
  double        crd(int idx) const;                 // get single element
  double        norm() const;                       // Eukl. norm
  double        sum() const;                        // Sum of items
  void          set(int idx, double newValue);      // set single element
  void          set(const ColumnVector&);           // set array by ColumnVector
  void          set(double crd[3]);                 // set array by array
  double*       crd_array();                        // get array
  ColumnVector  crd_cvect();                        // get ColumnVector
  t_gtriple&    crd_tripl();                        // get triple
  ColumnVector  unitary();                          // get unit ColumnVector
  t_gpair       gpair() const;
  bool          zero();                             // true: zero elements, false: not zero elements
   
 protected:

 private:
  double         _crd[3];

};

} // namespace

#endif
