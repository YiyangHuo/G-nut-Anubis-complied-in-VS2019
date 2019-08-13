
#ifndef GFILE_H
#define GFILE_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements file io based on gio class
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include <stdio.h>
#include <fstream>
#include <string>

#include "gio/gio.h"
#include "gio/giof.h"
#include "gio/giofz.h"

// special buffer size for file reading
// --> must be bellow gcoder maximum limit !
#define FILEBUF_SIZE 20480
#define FILEHDR_SIZE    48

using namespace std;

namespace gnut {

class t_gfile : public t_gio {

 public:
   t_gfile();
   virtual ~t_gfile();
   
   virtual int init_write();
   virtual int init_read();
   
   virtual int irc()const;                   // get irc status
   virtual bool eof();                       // integrate gzip/ascii
   virtual string md5sum();                  // get md5sum
   virtual string mask();                    // integrate gzip/ascii

   virtual bool compressed(){ return _gzip; }

   virtual int path( string str );           // set/get file://dir/name
   virtual string path();
   virtual string name();
   virtual void reset();                     // reset path

 protected:
   virtual int _gio_write( const char* buff, int size );
   virtual int _gio_read(  char* buff, int size );
   virtual int _stop_common();               // reimplementaion
   virtual void _set_gzip( string name );    // compressed ?
   
   virtual int _read( char* b, int s);       // integrate gzip/ascii
   virtual int _write( const char* b, int s);// integrate gzip/ascii

   int         _irc;
   bool        _gzip;                        // compressed
   t_giof*     _file;                        // ascii file
   t_giofz*   _zfile;                        // gzip  file

 private:
};

} // namespace
 
#endif
