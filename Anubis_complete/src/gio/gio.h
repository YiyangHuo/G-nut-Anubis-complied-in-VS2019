
#ifndef GIO_H
#define GIO_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements gtime class (day and precise time)
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include <stdio.h>
#include <fstream>
#include <string>

#include "gio/glog.h"
#include "gutils/gmutex.h"
#include "gcoders/gcoder.h"

#define BUF_SIZE 1024

using namespace std;

namespace gnut {

class t_gcoder;
class t_gio {

 public:
   t_gio();
   virtual ~t_gio();

   virtual void run_read();
//   virtual void run_read(void*){ run_read(); }
   virtual void run_write();

   virtual void coder( t_gcoder* coder ){ _coder = coder; }
   
   virtual int init_write(){ return _opened = _init_common(); }
   virtual int init_read(){  return _opened = _init_common(); }
   
   virtual int stop_write(){ _stop_common(); return 0; }
   virtual int stop_read(){  _stop_common(); return 0; }

   virtual int    path( string str );                      // set/get path
   virtual string path()const{ return _path; }

   void   file( const char* f ){ _giof.mask(f); }          // set/get local i/o file
   string file(){ return _giof.mask(); }

   void   verb( int i ){ _verb = i; }                      // set/get verbosity
   int    verb(){ return _verb; }

   void   stop(){ _stop = 1; }                             // stop
   size_t size(){ return _size; }                          // get size
   
   int   running(){  return _running; }
   int    opened(){  return _opened;  }
   int connected(){  return _opened;  }
   
   void glog(t_glog* l){ _log = l; }                       // set/get glog pointer
   t_glog* glog(){ return _log; }
   
   bool            operator<(const t_gio& n) const;
   bool            operator==(const t_gio& n) const;
   friend ostream& operator<<(ostream& os, const t_gio& n);
   
 protected:
   virtual int _gio_write( const char* buff, int size ) = 0;
   virtual int _gio_read(  char* buff, int size ) = 0;
   
   virtual int _locf_write( const char* buff, int size );
   virtual int _locf_read(  char* buff, int size );
   
   virtual int _init_common();
   virtual int _stop_common();

   t_glog*         _log;              // log pointer
   int             _fd;               // file descriptor
   size_t          _size;             // buffer size
   string          _path;             // URL-like path (e.g. file:///home/honza/file)
   t_giof          _giof;             // local file
//   char*           _loc_buff;         // local buffer
   int             _count;            // record counter
   int             _verb;             // verbosity
   int             _stop;             // require a stop at run() loop
   int             _opened;           // 1: opened/connected
   int             _running;          // running
   t_gcoder*       _coder;            // decoder/encoder
   t_gmutex        _gmutex;
   
#ifdef BMUTEX
   boost::mutex    _mutex;            // mutual exlusion
#endif

 private:
     
};

} // namespace

#endif
