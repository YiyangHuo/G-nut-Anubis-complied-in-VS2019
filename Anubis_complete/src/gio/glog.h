
#ifndef GLOG_H
#define GLOG_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements glog class derived from giof
  Version: $ Rev: $

  2011-10-10 /JD: created

-*/

#include <string>
#include <vector>

#include "gio/giof.h"

#define CACHE_LINES 30

using namespace std;

namespace gnut {

class t_glog : public t_giof {

 public:
   t_glog( string mask = "" );
   virtual ~t_glog();
   
   void comment(int l, const string& str);
   void comment(int l, const string& ide, const string& str);

   void cache_size(int i){ _size = i; }
   int  cache_size(){ return _size; } const
	
   void verb(int i){ _verb = i; }
   int  verb(){ return _verb; } const

   void clear(){ _cache.clear(); }

 protected:

   int             _verb;             // verbosity
   int             _size;             // cache size
   vector<string>  _cache;            // cache for messages
   t_gmutex        _log_gmutex;       // special mutex for comments
   
#ifdef BMUTEX
   boost::mutex    _log_mutex; 
#endif
   
 private:
     
};

} // namespace

#endif
