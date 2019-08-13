
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.
 
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

/* paths:
 *           file://dir/file
 *           http://user:pasw@host:port
 *            ftp://user:pasw@host:port
 *            tcp://host:port
 *          ntrip://mount/user:pasw@host:port
-*/

#include <cstring>
#include <iostream>

#include "gio/gio.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// constructor
// ---------
t_gio::t_gio()
 : _log(0),
   _fd(-1),
   _size(BUF_SIZE),
   _path(""),
   _giof(),
   _count(0),
   _verb(0),
   _stop(0),
   _opened(0),
   _running(0) 
{
  gtrace("t_gio::construct");

  _coder = 0;
//  char* loc_buff = new char[size()];
}


// destructor
// ---------
t_gio::~t_gio()
{
  gtrace("t_gio::destruct");
}


// set path
// ---------
int t_gio::path( string str )
{
  gtrace("t_gio::path(str)");

  if( str == "" ) return -1;
   
  _path = str;
  return 1;
}


// overloading << operator
// ---------
ostream& operator<<(ostream& os, const t_gio& n)
{
  os << n.path();
  return os;
}


// overloading == operator
// ---------
bool t_gio::operator==(const t_gio& n) const
{
  return ( n.path() == _path );
}


// overloading < operator
// ---------
bool t_gio::operator<(const t_gio& n) const
{
  return ( n.path() < _path );
}


// start reading (could be run in a separate thread)
// ---------
void t_gio::run_read()
{ 
  gtrace("t_gio::run_read");

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  if( _coder ){ _coder->clear(); } // be sure always clean decoder
  vector<string> errmsg;
   
  char* loc_buff = new char[_size];

  if( init_read() < 0 ){
    if( _log ) _log->comment(0,"gio"," warning - initialization failed");
    t_gtime::gmsleep(3000); // [ms]
//  sleep(3);               // [s] linux only
  }
  int nbytes  = 0;
  int decoded = 0;
  _stop = 0;
  _running = 1;

//  cout << "gio_read:" << _gio_read( loc_buff, _size ) << endl;

  while( (( nbytes = _gio_read( loc_buff, _size ) ) > 0 )
	  && _stop != 1 ){
     
//    for(int i = 0; i<nbytes ; i++){ cout << loc_buff[i]; }

    // archive the stream
    _locf_write(loc_buff,nbytes);

    if( _coder && nbytes > 0 ){
       decoded = _coder->decode_data( loc_buff, nbytes, _count, errmsg );
    }
    else {
       if( _log ) _log->comment(2,"gio","0 data decoded");       
       decoded = 0;
    }

    if( _log && _log->verb() >= 9) _log->comment(9,"gio","read: " + int2str(nbytes)
			                           + " decoded: " + int2str(decoded));
  }

  _stop_common(); 
  delete [] loc_buff;

  _gmutex.unlock(); return;
}


// start reading (could be run in a separate thread)
// ---------
void t_gio::run_write()
{
  gtrace("t_gio::run_write");

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( _coder ){ _coder->clear(); } // be sure always clean decoder
	
  vector<string>  errmsg;
  _stop = 0;

  char* loc_buff = new char[_size];
   
  if( _opened != 1 && init_write( ) < 0 ){
    if( _log ) _log->comment(0,"gio"," warning - initialization failed");
    _gmutex.unlock(); return;
  }

  int nbytes = 0;
  _running = 1;
  do{
     nbytes = 0;
     
     // 1. try to read from local file source
     nbytes = _locf_read(loc_buff,_size);

     // 2. try reading in a loop
//    while( nbytes < 0 ){
//       nbytes = _locf_read(loc_buff,_size);
//     }

     // 3. read from encoder
     if( nbytes < 0 && _coder ){
       nbytes = _coder->encode_data( loc_buff, _size, _count, errmsg );
     }

  // if( nbytes > 0 )   cout << " nbytes = encoded : " << nbytes << " \n";
  // t_gtime::gmsleep(100);
  // for(int i = 0; i<nbytes ; i++){ cout << loc_buff[i]; }
  }while( _stop == 0 && (_gio_write( loc_buff, nbytes ) > 0 ) );

  _stop_common();
  if( _log ) _log->comment(2,"gio","end of read [stop or gio_write problem]");
  delete [] loc_buff;
  _gmutex.unlock(); return;
}


// local log file source
// ---------
int t_gio::_locf_read( char* buff, int size )
{   
  gtrace("t_gio::_locf_read");
  
  if( _giof.mask() == "" ) return -1;

  return _giof.read(buff,size);
}


 // local log file archive
// ---------
int t_gio::_locf_write( const char* buff, int size )
{ 
  gtrace("t_gio::_locf_write");

  if( _giof.mask() == "" ) return -1;

  int irc = _giof.write( buff, size );
//  if( fprintf(_locf_ptr, buff, size) == 0 ) return 0;
   
//  for(int i = 0; i<size ; i++){ _locf_stream << buff[i]; }

//  if( irc ) return -1;
//  return size;
  return irc;
}


// common function for initialization
// ---------
int t_gio::_init_common()
{ 
  gtrace("t_gio::_init_common");

  return 1; 
}



// common function for socket/file close
// ---------
int t_gio::_stop_common()
{ 
  gtrace("t_gio::_stop_common");

  _opened  = 0;  // for TCP/NTRIP now, need to be checked for gfile! 
 _running = 0;
  return _running;
}
   
} // namespace