
#ifndef GOBS_H
#define GOBS_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: definition of GNSS observation types
  Version: $ Rev: $

  2012-09-26 /JD: created

-*/

#include <map>
#include <list>
#include <vector>
#include <string>
#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gutils/gnss.h"

using namespace std;

namespace gnut {

// class
// ----------
class t_gattr {
   
 public:
   t_gattr(){ _gattr = ATTR; };
   t_gattr(GOBSATTR a ){ _gattr = a; };
  ~t_gattr(){};

  virtual void attr(const GOBSATTR& a);              // set attr
  virtual GOBSATTR attr()const;                      // get attr
   
  virtual bool operator==(const t_gattr& g)const;
  virtual bool valid()const;
   
 protected:
  GOBSATTR _gattr;
};


// class
// ----------
class t_gband : public t_gattr {
   
 public:
   t_gband():t_gattr(){ _gband = BAND; };
   t_gband(GOBSBAND b, GOBSATTR a ):t_gattr(a){ _gband = b; };
  virtual ~t_gband(){};

  virtual void band(const GOBSBAND& g);              // set band
  virtual GOBSBAND band()const;                      // get band

  virtual void gattr(const t_gattr& g);              // set t_gattr
  virtual t_gattr gattr()const;                      // get t_gattr

  virtual bool operator==(const t_gband& g)const;
  virtual bool valid()const;
   
 protected:
  GOBSBAND _gband;
};


// class
// ----------
class t_gobs : public t_gband {

 public:
   t_gobs():t_gband(){ _gtype = TYPE; };
   t_gobs(GOBSTYPE t, GOBSBAND b, GOBSATTR a ):t_gband(b,a){ _gtype = t; };
   t_gobs(GOBS g){ gobs(g); };
  virtual ~t_gobs(){};
   
  virtual void type(const GOBSTYPE& t);               // set type (only, inherit band&attr!)  
  virtual GOBSTYPE type()const;                       // get type

  virtual void gband(const t_gband& g);               // set attr
  virtual t_gband gband()const;                       // get attr 

  int gobs(const GOBS& g);                            // set type (only! inherit)
  int gobs(const string& s);
   
  GOBS gobs()const;                                   // get gobs enum

  bool operator==(const t_gobs& g)const;

  bool valid()const;
  bool is_code()const;
  bool is_phase()const;
   
 protected:

  GOBSTYPE _gtype;
};

// ------------------------------------------------------------------------------------------------------
// STATIC MAPS
// ------------------------------------------------------------------------------------------------------

/*
// get system 1st band (preferred b)
// ----------
int freq_1st(GSYS gs, GOBSTYPE b=BAND_1)
{
  if( b < 1 ||  b > 8 ) b = 1;
  if( b > 3 &&  b < 5 ) b = 1;

  // exceptions
  if( gs == GAL  && b == 2 ) return 1;
  if( gs == SBS  && b == 5 ) return 1;
  if( gs == BDS  && b == 2 ) return 1;
   
  return b;
}
*/

} // namespace

#endif // GOBS_H

