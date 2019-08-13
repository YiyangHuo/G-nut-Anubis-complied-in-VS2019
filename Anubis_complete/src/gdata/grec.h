
#ifndef GREC_H
#define GREC_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements receiver class
  Version: $ Rev: $

  2011-01-10 /JD: created

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

class t_grec : public t_gobj {

 public:
   t_grec();
//   t_grec(const t_grec& obj);
   virtual ~t_grec();

   typedef map<t_gtime,string>      t_maprec;
   typedef map<t_gtime,t_rnxhdr>    t_maphdr;
   
   virtual void addhdr(const t_rnxhdr& hdr, const t_gtime& epo, string path);
   t_maphdr gethdr();
   t_rnxhdr gethdr(const t_gtime& epo);
   
   t_maprec get_maprec() {return _maprec;}

   void rec(string rec, const t_gtime& beg, const t_gtime& end = LAST_TIME );
   string rec(const t_gtime& t) const;       // set/get receiver  
   
   virtual bool isrec() {return true;}
   virtual bool istrn() {return false;}
   
   virtual void compare(shared_ptr<t_grec> grec, const t_gtime& tt);

   virtual vector<t_gtime> rec_id() const;

   void fill_rnxhdr(const t_rnxhdr& rnxhdr);
   
 protected:
   void _fill_rnxhdr(const t_rnxhdr& rnxhdr);
   t_rnxhdr _gethdr(const t_gtime& epo);   
   
   void   _rec(string rec, const t_gtime& beg, const t_gtime& end = LAST_TIME );
   string _rec(const t_gtime& t) const;  
   
   t_maprec        _maprec;                  // map of receviers
   t_maphdr        _maphdr;                  // map of rinex header information
   
 private:
};

} // namespace

#endif
