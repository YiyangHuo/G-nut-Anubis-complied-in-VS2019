
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

#include <cstring>
#include <sstream>

#include "gio/gfile.h"
#include "gutils/gcommon.h"
#include "gutils/gfileconv.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gfile::t_gfile()
 : t_gio(),
   _irc(0),
   _gzip(false)
{
  gtrace("t_gfile::construct");
  
  _file  = 0;
  _zfile = 0;
  _size = FILEBUF_SIZE;
}


// destructor
// ----------
t_gfile::~t_gfile()
{
  gtrace("t_gfile::destruct");

  reset();
}


// set path
// ----------
int t_gfile::path( string str )
{
  gtrace("t_gfile::path(str)");

  if( str == "" ) return -1;

  size_t idx1 = 0, idx2 = 0;
  string prefix = GFILE_PREFIX;

  // check if file
  idx2 = str.find(prefix);
  if( idx2 == string::npos ){ str = prefix+str; idx2 = 0; }

  // file path
  idx1 = idx2+prefix.length();
  idx2 = str.length();
  if( idx2 != string::npos && idx2 > idx1  ){
      
    string name = str.substr(idx1,idx2-idx1);
    _set_gzip(name);

    if( _gzip && _zfile == 0 ) _zfile = new t_giofz(name.c_str());
    else if(      _file == 0 ) _file  = new t_giof( name.c_str());
    else return -1;

    if( _gzip && _zfile ) _zfile->mask(name);
    else if(      _file )  _file->mask(name);
    else{ cerr << "Error: cannot assign file to path!\n";  _irc++; return -1;  }

    ostringstream ltmp; 
    ltmp << "File: " << int(idx1) << ":" << int(idx2) << " = " << name;
    if( _log ) _log->comment(3, "gfile", ltmp.str());
    
  }else{
    if( _log ) _log->comment(0, "gfile", "warning: path does not contain file://dir/name  [check file://]");
    return -1; 
  }

  t_gio::path(str);
   
  return 1;
}


// get path
// ----------
string t_gfile::path(){
  return GFILE_PREFIX + mask();
}


// get name
// ----------
string t_gfile::name(){
   
  return mask();
}


// irc
// ----------
int t_gfile::irc()const
{
  if( _gzip && _zfile ) return (_zfile->irc() + _irc);
  else if(     _file  ) return ( _file->irc() + _irc);
  
  return _irc;
}


// eof
// ----------
bool t_gfile::eof()
{
  if( _gzip && _zfile ) return _zfile->eof();
  else if(     _file  ) return  _file->eof();
  
  return true;
}


// mask
// ----------
string t_gfile::mask()
{
  if( _gzip && _zfile ) return _zfile->mask();
  else if(      _file ) return  _file->mask();
  
  return "";
}


// get md5
// ----------
string t_gfile::md5sum()
{
  gtrace("t_gfile::md5sum");
  
  if( _gzip && _zfile ) return _zfile->md5sum();
  else if(      _file ) return  _file->md5sum();
  
  return "";
}

   
// init read function (read header)
// ----------
int t_gfile::init_read()
{
  gtrace("t_gfile::init_read");
  
  if( ! _coder ){ return 0; }

  char* loc_buff = new char[FILEHDR_SIZE];

#ifdef DEBUG
  cout << "HEADER START " << mask() << "\n"; cout.flush();
#endif

  int nbytes = 0;
  vector<string>  errmsg;
//  while( ( nbytes = _gio_read( loc_buff, FILEHDR_SIZE ) ) >= 0 && _stop != 1 ){ // JD 2014-10-25
  while( ( nbytes = _gio_read( loc_buff, FILEHDR_SIZE ) ) > 0 && _stop != 1 ){

#ifdef DEBUG
      cout << "GIO_HEAD_READ:"; 
      for(int i = 0; i<nbytes ; i++){ cout << loc_buff[i]; } 
      cout << ":END\n";
#endif
     
    if( _coder->decode_head( loc_buff, nbytes, errmsg ) < 0 ) break;
  }

  if( nbytes < 0 ){ // header not completed properly ?
    ++_irc; if( _coder ) _coder->mesg(GERROR, "Error: Incomplete header identified."); 
  }

#ifdef DEBUG
  cout << "HEADER END " << mask() << "\n";
#endif

  delete [] loc_buff;
  return 1; 
}


// init write function (write header)
// ----------
int t_gfile::init_write()
{ 
  gtrace("t_gfile::init_write");  

  if( ! _coder ){ return 0; }
   
  char* loc_buff = new char[FILEHDR_SIZE];

#ifdef DEBUG
  cout << "HEADER START " << mask() << endl; cout.flush();
#endif
  vector<string> errmsg;
   
  int nbytes = 0;
  do{
    if( ( nbytes = _coder->encode_head( loc_buff, FILEHDR_SIZE, errmsg ) ) < 0 ) break;
  }while( (nbytes > 0) && (_gio_write( loc_buff, nbytes ) > 0) && _stop != 1  );

#ifdef DEBUG
  cout << "HEADER END " << mask() << endl; cout.flush();
#endif
   
  delete [] loc_buff;
      
  return 1; 
}

// reset
// ----------
void t_gfile::reset()
{     
  gtrace("t_gfile::reset");

  if(  _file ){ delete _file;   _file = 0; }
  if( _zfile ){ delete _zfile; _zfile = 0; }
}


// read data
// ----------
int t_gfile::_gio_read( char* buff, int size )
{     
  gtrace("t_gfile::_gio_read");

  if( mask() == "" ) return -1;

  int nbytes =  this->_read(buff,size);
  if( nbytes == 0 && this->eof() ) return -1;
  return nbytes;

//  int nbytes = 0
//  nbytes = this->_read(buff,size);
//  if( this->eof() || nbytes == -1 ){
//    if( nbytes == 0 ) nbytes = -1;
//    }
//  }
//  return nbytes;
}


// send data
// ----------
int t_gfile::_gio_write( const char* buff, int size )
{
  gtrace("t_gfile::_gio_write");

  if( mask() == "" ) return 0;

  this->_write( buff, size );

  return size;
}


// common function for file close
// ----------
int t_gfile::_stop_common()
{ 
  gtrace("t_gfile::_stop_common");

  return t_gio::_stop_common();
}

   
// compressed or not ?
// ----------
void t_gfile::_set_gzip( string name )
{ 
  gtrace("t_gfile::_set_gzip");

  if( // name.compare(name.length()-2,name.length(),".Z")  == 0 || ///  NOT SUPPORTED
      name.compare(name.length()-3,name.length(),".gz") == 0 )

        _gzip = true;
   else _gzip = false;

}


// ----------
int t_gfile::_read( char* b, int s )
{
  gtrace("t_gfile::_read");

  if( _gzip && _zfile ) return _zfile->read( b, s );
  else if(     _file  ) return  _file->read( b, s );
  
  return -1;
}
   

// ----------
int t_gfile::_write( const char* b, int s )
{
  gtrace("t_gfile::_write");

  if( _gzip && _zfile ) return _zfile->write( b, s );
  else if(     _file  ) return  _file->write( b, s );
  
  return -1;
}


} // namespace
