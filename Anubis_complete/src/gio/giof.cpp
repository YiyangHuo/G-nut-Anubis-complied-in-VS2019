
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

#include <iostream>

#include "md5/md5.h"
#include "gio/giof.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"
#include "gutils/gfileconv.h"

using namespace std;

namespace gnut {

// constructor
// ---------
t_giof::t_giof( string mask )
  : _irc(0),
    _mask(mask),
    _name(""),
    _omode(ios::trunc),
    _repl(false),
    _toff(0),
    _loop(false),
    _tsys(t_gtime::UTC)
{
  substitute(mask,GFILE_PREFIX,"");

  this->exceptions( fstream::failbit | fstream::badbit );
  if( _mask.find('%') != string::npos ) _repl = true;
  _name = _replace();

}


// destructor
// ---------
t_giof::~t_giof()
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  if( this->is_open() ){
    try{
      this->close();
    }catch( exception e ){
      cerr << "Error: Exception closing file " << _name << ": " << e.what() << endl;
    }
  }
  _gmutex.unlock();
}


// set mask
// ---------
int t_giof::mask( string mask )
{
  gtrace("t_giof::mask "+_name);
   
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();

  _irc = 0;
  substitute(mask,GFILE_PREFIX,"");

  this->clear();
  if( this->is_open() ){
    try{
#ifdef DEBUG
	  cout << "mask: closing log-file: " << _name << endl;
#endif  
      this->close();
    }catch( exception e ){
      cerr << "Error: Exception closing file " << _name << ": " << e.what() << endl;
      _irc++; _gmutex.unlock(); return -1;
    }
  }
  if( mask.find('%') != string::npos ) _repl = true;
  _mask = mask;
  _name = _replace();
  _gmutex.unlock(); return 1;
}


// local log file archive
// ---------
int t_giof::write( const string& s )
{
  gtrace("t_giof::write(string) "+_name);

  return this->write( s.c_str(), s.size() );
}
   

// local log file archive
// ---------
int t_giof::write( const char* buff, int size )
{     
  gtrace("t_giof::write(buff,size) "+_name);
   
  if( _mask == "" || size < 1 ) return -1;

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  string name = _replace();

  this->clear();
  // new name || not opened
  if( ! this->is_open() || name !=  _name ){
     
    if( this->is_open() ){
      try{ 
        this->close();
#ifdef DEBUG
        cout << "closing log-file: " << _name << " [" << name << "]\n";
#endif  
      }catch( exception e ){
        cerr << "Error: Exception closing file " << _name << ": " << e.what() << endl;
        _irc++; _gmutex.unlock(); return -1;
      }
    }
     
    try{ 
      _name = name; make_path( _name );
      this->open( _name.c_str(), ios::binary | ios::out | _omode );
#ifdef DEBUG
     cout << " opening file " << _name << endl;
#endif
    }catch( exception e ){
      cerr << "Error: Exception opening file " << _name << ": " << e.what() << endl;
      _irc++; _gmutex.unlock(); return -1;
    }
  }

  try{
    fstream::write( buff, size );
    fstream::flush();
  }catch( exception e ){
    cerr << "Error: Exception writting to file " << _name << ": " << e.what() << endl;
    _irc++; _gmutex.unlock(); return -1;
  }

  _gmutex.unlock(); return size;
}


// local log file source
// ---------
int t_giof::read( char* buff, int size )
{
  gtrace("t_giof::read(buff,size) "+_name);
   
  if( _mask == ""   ) return -1;

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _name = _replace();

  int nbytes;
  /* we might need to read twice because the last read read exactly until eof,
   * thus nbytes/gcount() would be 0 this time , which we don't want */
  for(int i = 0; i < 2; i++) {
      
  this->clear();
  if( ! this->is_open() ){
    try{
     make_path( _name );
     this->open( _name.c_str(), ios::binary | ios::in );
    }catch( exception e ){
      cerr << "Error: Exception opening file " << _name << ": " << e.what() << endl;
      _irc++; _gmutex.unlock(); return -1;
    }
  }

  nbytes = size;
  try{	
    fstream::read( buff, size );
  }catch( exception e ){
    if( this->eof() ){
      nbytes = (int)this->gcount();

      if( _loop ) this->close();

      if( nbytes < 1 ){
          // there was no data left before eof, let's try opening the file again
          continue; 
      }
    }else{	  
      cerr << "Error: Exception reading from file " << _name << ": " << e.what() << endl;
      _irc++; _gmutex.unlock(); return -1;
    }
  }
  break; // nbytes is more than 0, we can break
  }
#ifdef DEBUG
  cout << "giof_read  [" << nbytes << "] " << buff << "\n name "  << _name << "\n"; cout.flush();
#endif

  _gmutex.unlock(); return nbytes;
}


// append mode
// ---------
void t_giof::append(const bool& b)
{ 
  gtrace("t_giof::append "+_name);

  if( b ) _omode = ios::app;
  else    _omode = ios::trunc;
}


// set time system for replacement
// ---------
void t_giof::tsys(t_gtime::t_tsys ts)
{ 
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _tsys = ts;
  _gmutex.unlock();
   return;
}


// get time system for replacement
// ---------
t_gtime::t_tsys t_giof::tsys()
{ 
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   t_gtime::t_tsys tsys = _tsys;
  _gmutex.unlock();
  return tsys;
}


// get md5sum
// ----------
string t_giof::md5sum()
{
  gtrace("t_giof::md5sum "+_name);
   
  stringstream buffer;
  ifstream tmp( mask() );
  buffer << tmp.rdbuf();

  return md5( buffer.str() );
}


// evaluate name
// ---------
string t_giof::_replace()
{ 
  if( ! _repl ) return _mask;
  
  t_gtime file_tm = t_gtime::current_time(_tsys);
          file_tm.add_secs(_toff*60.0);
  return  file_tm.str(_mask);
}


} // namespace
