
#ifndef GPAIR_H
#define GPAIR_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements 2D coordinates representation (e.g. horizontal coordinates)
  Version: $ Rev: $

  2012-05-11 /JD: created

-*/

#include <iostream>
#include <string.h>

#include "../newmat/newmat.h"
#include "../newmat/newmatio.h"

using namespace std;

namespace gnut {

class t_gpair {

 public:
  t_gpair();
  t_gpair(double x, double y);
  t_gpair(double crd[2]);
  t_gpair(const ColumnVector& crd);
  virtual ~t_gpair();
  
  t_gpair&   operator=(const t_gpair& other);       // assignment operator
  t_gpair    operator+(const t_gpair& other) const; //
  bool       operator==(const t_gpair& tr) const;   // equal operator
  bool       operator<( const t_gpair& tr) const;
  double&    operator[](const size_t idx);          // get a reference of element
  double     operator[](const size_t idx) const;    // get value of element

  friend ostream& operator<<(ostream& os, const t_gpair& x);
   
  double        crd(int idx) const;                 // get single element
  void          set(int idx, double newValue);      // set single element
  void          set(const ColumnVector&);           // set array by ColumnVector
  void          set(double crd[2]);                 // set array by array
  double*       crd_array();                        // get array
  ColumnVector  crd_cvect();                        // get ColumnVector
  t_gpair&      crd_pair();                         // get pair
  ColumnVector  unitary();                          // get unit ColumnVector
  bool          zero();                             // true: zero elements, false: not zero elements

 protected:

 private:
  double         _crd[2];

};

} // namespace

#endif
