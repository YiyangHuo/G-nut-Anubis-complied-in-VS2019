
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

#include <vector>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <algorithm>

#include "gcoders/gcoder.h" 
#include "gutils/gtypeconv.h"
#include "gutils/gfileconv.h"
#include "gutils/gtimesync.h"
 
using namespace std;

namespace gnut {

// constructor
// ----------
t_gcoder::t_gcoder(t_gsetbase* s, string version, int sz, string id )
 : _class(id),
   _initialized(false),
   _gnss(true),
   _initver(version),
   _close_with_warning(true)
{
  _init();
   
  // after init!
//  _version = version; 
  _buffsz  = sz;
  _log     = 0;
  _set     = 0;
  _irc     = 0;
  _hdr     = false;
  
  _ss_position = 0;   
   
  _out_smp = 0;
  _out_len = 0;
  _out_epo = LAST_TIME;
   
// use malloc instead of new due to realocate function!
// ====================================================
  _buffer = (char*)malloc((_buffsz+1)*sizeof(char));  // due to realocate function!
   
  if( s ) _gset(s); // READ SETTINGS IN CONSTRUCT ONLY!

//  if( data != 0 ) _data[ data_id ] = data;

   _pgm = "Anubis";
}


/* ----------
 * destructor
 */
t_gcoder::~t_gcoder()
{
  // use free instead of delete due to realocate function!
  free(_buffer);

  if( _close_with_warning )
  {	
    std::sort( _notes.begin(), _notes.end() );
    for( auto it = _notes.begin(); it != _notes.end(); ++it ){
      if( _log ){ _log->comment( 0, "gcoder", it->str()+" .. "+base_name(_fname)); }
//    if( _log ){ _log->comment(-1, "gcoder", _fname+"  "+it->str()); }
      else{ cerr << it->status() << " " << it->str() << endl; }
    }
  }
}


/* ----------
 * read settings
 */
int t_gcoder::_gset(t_gsetbase* s)
{
  gtrace("t_gcoder::gset");

  if( ! s ) return -1;

  _set = s;
  _beg = dynamic_cast<t_gsetgen*>(_set)->beg();
  _end = dynamic_cast<t_gsetgen*>(_set)->end();
   
  _hdr = dynamic_cast<t_gsetgen*>(_set)->thin();
  
  if( _gnss )
  _sys = dynamic_cast<t_gsetgen*>(_set)->sys();
   
  _rec = dynamic_cast<t_gsetgen*>(_set)->rec();   
  _int = dynamic_cast<t_gsetgen*>(_set)->sampling();
  _scl = dynamic_cast<t_gsetgen*>(_set)->sampling_scalefc(); // scaling 10^decimal-digits

  // individually for each GNSS
//  t_map_gnss::const_iterator itGNSS;
//  for( itGNSS = GNSS_DATA_PRIORITY.begin(); itGNSS != GNSS_DATA_PRIORITY.end(); ++itGNSS )
  set<string>::const_iterator itGNSS;
  for( itGNSS = _sys.begin(); itGNSS != _sys.end(); ++itGNSS )
  {	
//    GSYS gsys = itGNSS->first;
    string gs = *itGNSS;
    GSYS gsys = t_gsys::str2gsys( gs );

    _sat[gsys] = dynamic_cast<t_gsetgnss*>(_set)->sat( gsys, false ); // empty if not set (to speed up!)
    _obs[gsys] = dynamic_cast<t_gsetgnss*>(_set)->obs( gsys, false ); // empty if not set (to speed up!)
    _nav[gsys] = dynamic_cast<t_gsetgnss*>(_set)->nav( gsys, false ); // empty if not set (to speed up!)

    // extend gobs list with completely defined signals
    set<string> sgobs = dynamic_cast<t_gsetgnss*>(_set)->gobs( gsys);
    for(auto it = sgobs.begin(); it != sgobs.end(); it++) _obs[gsys].insert(*it);
    
#ifdef DEBUG
    set<string>::const_iterator it;
    ostringstream os;
    os << gs << " NAV:";
    for( it = _nav[gsys].begin(); it != _nav[gsys].end(); ++it ) os << " " << *it;
    os << endl;

    os << gs << " OBS:";
    for( it = _obs[gsys].begin(); it != _obs[gsys].end(); ++it ) os << " " << *it;
    os << endl;

    os << gs << " SAT:";
    for( it = _sat[gsys].begin(); it != _sat[gsys].end(); ++it ) os << " " << *it;
    os << endl;

    cout << "GNSS filtering:\n" + os.str(); cout.flush();
#endif
  }

#ifdef DEBUG   
  cout << "GCODER settings:"
       << " _beg " << _beg.str("%Y-%m-%d %H:%M:%S")
       << " _end " << _end.str("%Y-%m-%d %H:%M:%S")
       << " _int " << _int
       << " _sys " << _sys.size() 
       << " _rec " << _rec.size()
       << " _sat " << _sat.size()
       << " _nav " << _nav.size()
       << " _obs " << _obs.size() << endl; cout.flush();
#endif
  return 0;
}


// set path
// ----------
void t_gcoder::path(string s)
{
   substitute(s,GFILE_PREFIX,"");

  _fname = s; 
}
   

/* ----------
 * add message
 */
void t_gcoder::mesg(t_note n, string s)
{
  t_gnote m(n,s);

  //if( find(_notes.begin(), _notes.end(), m) == _notes.end() ){
  //     _notes.push_back(m);

  // eliminate repeating messages
  bool exist = false;
  vector<t_gnote>::const_iterator it;
  for( it = _notes.begin(); it != _notes.end(); ++it ){
    if( *it == m ){ exist = true; }
  }
   
  if( !exist ) _notes.push_back(m);
}

/* ----------
 * get message
 */
vector<t_gnote> t_gcoder::mesg()
{
  return _notes;
}

   
/* ----------
 * initialize
 */
void t_gcoder::_init()
{
  gtrace("t_gcoder::_init");

  _version = _initver;
  _endpos  = 0;
}


/* ----------
 * clear
 */
void t_gcoder::clear()
{
  gtrace("t_gcoder::clear");

  // use free instead of delete due to realocate function!
  free(_buffer);  
   
  _init();
   
// use malloc instead of new due to realocate function!
// ====================================================
  _buffer = (char*)malloc((_buffsz+1)*sizeof(char));  // due to realocate function!
}


/* ----------
 * get single line from the buffer
 */
int t_gcoder::_getline( string& str, int from_pos )
{
  gtrace("t_gcoder::_getline");

  str = "";
#ifdef DEBUG
  cout << " from:" << from_pos << " to(end):"  << _endpos << endl;
#endif
  if( _endpos == 0 || from_pos >= _endpos) return -1;

  string tmp(_buffer); 
  size_t ifirst = 0;

  if( ( ifirst = tmp.find(crlf,from_pos)) == string::npos // &&
//    (ifirst = tmp.find(eof,tmp.size())) == string::npos        
//       tmp.c_str()[tmp.size()] != eof                         // end what about end of line ???
    ) return -1;
   
  str = tmp.substr( from_pos, ifirst+1 - from_pos); // + CRLF

  return str.length();
}


/* ----------
 * cummulate buffer
 */
int t_gcoder::_add2buffer( char* buff, int sz )
{   
  gtrace("t_gcoder::_add2buffer");

  if( sz == 0 && endpos() > 0 ) return endpos();
#ifdef DEBUG
  for(int i=0; i<sz; i++ ){
     cout << " IN " << i << " [" << *(buff+i) << "] ";
  }
#endif

  if( _endpos + sz > _buffsz ){
    int old = _buffsz;
     
    // INCREASE BUFFER USING REALLOC !
	_buffsz = (int)((_endpos + sz + 1)*BUFFER_INCREASE_FAC); //  (int)BUFFER_INCREASE_FAC == 1
    if(! (_buffer = (char *)realloc(_buffer, _buffsz*sizeof(char))) ){
      if( _log )  _log->comment(0,"gcoder","Error - reallocation not succesful!");
      exit(1);
    }

    // DO NOT EXCEED MAXIMUM BUFFER SIZE !
    if( _buffsz > MAXIMUM_BUFFER_SIZE ){
      if( _log )  _log->comment(0,"gcoder","Error - a decoding problem and buffer limited to " + int2str(MAXIMUM_BUFFER_SIZE));
      _consume(sz);

      // --> solved in individual gcoders for individual files!
      return -1; 
      // rather to kill the Anubis here !?
//    exit(1);

    }else{	       
      if( _log )  _log->comment(2,"gcoder","buffer encreased from " + int2str(old) + " to " + int2str(_buffsz));
    }
  }
   
  for(int i=0; i<sz; i++, _endpos++ ){
     *(_buffer + _endpos) = *(buff+i);
  }

  _buffer[_endpos] = '\0'; // temporary set end of string (will be occasionally replace by add2buff)

#ifdef DEBUG             // HERE WE CAN PRINT _BUFFER SINCE IT ENDS WITH '\0'
  cout << "add2buffer : " << sz << " -> " << _endpos << " BUFF:" << _buffer << ":END\n"; cout.flush();
#endif

  return sz;
}



/* ----------
 * remove from buffer
 */
int t_gcoder::_consume( int bytes_to_eat )
{
  gtrace("t_gcoder::_consume");

  if( bytes_to_eat <= 0        ){ return 0; }
  if( bytes_to_eat >  _endpos ||
      bytes_to_eat >  _buffsz  ){
     
    if( _log )  _log->comment(0,"gcoder","error: too many bytes to eat: " + int2str(bytes_to_eat)
                                                                    + " " + int2str(_buffsz)
                                                                    + " " + int2str(_endpos));
    else                cerr << "gcoder - error: too many bytes to eat: " + int2str(bytes_to_eat)
                                                                    + " " + int2str(_buffsz)
                                                                    + " " + int2str(_endpos) << endl;
    exit(1);
  }
   
  for( int eat = 0; eat + bytes_to_eat < _endpos; eat++ ){
    *(_buffer + eat) = *(_buffer + eat + bytes_to_eat );
  }
  _endpos -= bytes_to_eat;
  _buffer[_endpos] = '\0'; // temporary set end of string (will be occasionally replace by add2buff)

#ifdef DEBUG
  cout << "consumed: " << bytes_to_eat << " ... endpos = " << _endpos << "\n";
#endif  
  return bytes_to_eat;
}


// add specific data structure to the coder
// ----------
int t_gcoder::add_data( string data_id, t_gdata* data )
{
  gtrace("t_gcoder::add_data");

  if( data_id.empty() || data == 0 ) return -1;

  map<string,t_gdata*>::iterator it = _data.find(data_id);
  if( it != _data.end() ){
    if( _log ) _log->comment(1,"gcoder","warning: structure " + data_id + " already exists !" );
    return -1;
  }

#ifdef DEBUG
  cout << " gcoder: adding data structure: " << data_id << endl; cout.flush();
#endif
  _data[data_id] = data;

  // can be used to directly setup data containers in individual gcoders
  _add_data(data_id,data);

  return 1;
}

/*
// remove specific data structure from the coder
// ----------
int t_gcoder::rem_data( string data_id )
{
  if( data_id.empty() ) return -1;

  map<string,t_gdata*>::iterator it = _data.find(data_id);
  if( it != _data.end() ) return -1;

#ifdef DEBUG
  cout << " gcoder: removing data structure: " << data_id << endl;
#endif

  _data.erase(it);
   
  // directly remove data in
//  _rem_data();
   
  return 1;
}
*/

// remove specific data structure from the coder
// ----------
int t_gcoder::_fill_buffer( char* buff, int sz )
{
  gtrace("t_gcoder::_filter_buffer");

  _ss.seekg(0, _ss.end);
  long len = (long)_ss.tellg();
  len = len - _ss_position;
  _ss.seekg(_ss_position, ios_base::beg);

  int size = ( len < sz ) ? len : sz;

  _ss.clear();
  _ss.read(buff,size);

   
  if( _ss.fail() ){ cout << "HEAD: any problem ?\n";

  }else if( _ss.gcount() == 0 ){
//  cout << "HEAD: finished " << sz << " " << len << " " << _ss_position << " " << ss.gcount() << endl;
    _ss_position = 0;
    _ss.str(""); 
    _ss.clear();
  }else{
//  cout << "HEAD: read     " << sz << " " << len << " " << _ss_position << " "  << _ss.gcount()<< endl;
    _ss_position += size; // _ss.gcount();
  }

  return size;
}


// sampling filter for epochs (return true if the epoch fits sampling)
// ----------
bool t_gcoder::_filter_epoch(const t_gtime& epo)
{
  gtrace("t_gcoder::_filter_epoch");

  if( time_sync(epo,_int,_scl,_log) ){ 
//   cerr << epo.str_ymdhms("fitting: ")+dbl2str(epo.dsec()) << endl;
    return true;
  }
  return false;
}

   
// GNSS/sat filter (return true if the epoch fits gnss & sat)
// ----------
bool t_gcoder::_filter_gnss(const string& prn)
{
  gtrace("t_gcoder::_filter_gnss");
   
  string gs = t_gsys::char2str( prn[0] ); // DONT USE SYS FROM ARG WHICH MAY BE EMPTY
  GSYS gsys = t_gsys::str2gsys( gs );     // DONT USE SYS FROM ARG WHICH MAY BE EMPTY
  if(( _sys.size()       == 0 ||       _sys.find(  gs ) !=  _sys.end()      ) &&
     ( _sat[gsys].size() == 0 || _sat[gsys].find( prn ) != _sat[gsys].end() ))
  {
    return true;
  }
  return false;
}


} // namespace
