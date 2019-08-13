
#ifndef GTRN_H
#define GTRN_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements transmitter class
  Version: $ Rev: $

  2013-06-20 /PV: created

-*/

#include <stdio.h>
#include <string>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gdata/gobj.h"
#include "gdata/gdata.h"

using namespace std;

namespace gnut {

class t_gtrn : public t_gobj {

 public:
   t_gtrn();
//   t_gtrn(const t_gtrn& obj);
   virtual ~t_gtrn();

//   typedef map<t_gtime, int>        t_mapchk;  // for changing channel canal of Glonass
   
   virtual bool isrec() {return true;}
   virtual bool istrn() {return false;}
   
//   void chk(int chk, const t_gtime& beg, const t_gtime& end);
//   int  chk(const t_gtime& t) const;       // set/get receiver
  
   virtual void channel(int chk);
   virtual int  channel() const;
   
 protected:
//   t_mapchk        _mapchk;
  int _channel;   // temporaryly only one value. Must be enhance via _mapchk 
 private:
};

} // namespace

#endif
