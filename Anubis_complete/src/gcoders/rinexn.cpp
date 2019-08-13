
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <memory>

#include "gdata/gnav.h"
#include "gdata/geph.h"
#include "gall/gallobj.h" 
#include "gall/gallnav.h"
#include "gutils/gtime.h"
#include "gcoders/rinexn.h"
#include "gutils/gtypeconv.h"
#include "gutils/gsys.h"

using namespace std;

namespace gnut {

/* ----------
 * constructor
 */
t_rinexn::t_rinexn( t_gsetbase* s, string version, int sz )
: t_gcoder( s, version, sz )
{
  _gnsssys  = 'G';
//_check    = true;
  _check_dt = FIRST_TIME;
  
  _gset(s); // HAVE TO BE EXPLICITLY CALLED HERE (AT THE END OF CONSTRUCTOR)
}


/* ----------
 * NAV-RINEX header
 */
int t_rinexn::decode_head( char* buff, int sz, vector<string>& errmsg)
{ 

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();

  if( t_gcoder::_add2buffer( buff, sz) == 0 ){ _mutex.unlock(); return 0; };  

  string line;
  int consume = 0;
  int tmpsize = 0;
  while( ( tmpsize = t_gcoder::_getline( line )) >= 0 ){
    
    consume += tmpsize;
    
    if( line.find("RINEX VERSION",60) != string::npos ){          // first line
      switch ( line[20] ){
        case 'N': _gnsssys = 'G'; break; // Navigation data - according to Rinex specification
        case 'G': _gnsssys = 'R'; break; // GLONASS NAVIGATION - occures sometimes in brdc
        default : { string lg("warning: not rinex navigation data file");
          if( _log )  _log->comment(0,"rinexn", lg);
          else  cerr << lg;
          _mutex.unlock(); return -1;
        }
      }
       
      switch ( line[40] ){
        case 'G': _gnsssys = 'G'; break; // GPS
        case 'R': _gnsssys = 'R'; break; // GLONASS
        case 'E': _gnsssys = 'E'; break; // GALILEO
        case 'C': _gnsssys = 'C'; break; // BEIDOU
        case 'J': _gnsssys = 'J'; break; // QZSS
        case 'S': _gnsssys = 'S'; break; // SBAS
        case 'I': _gnsssys = 'I'; break; // IRNSS
        case 'M': _gnsssys = 'M'; break; // MIXED
        case ' ': { 
          if( line[20] == 'N' ) _gnsssys = 'G';
          if( line[20] == 'G' ) _gnsssys = 'R';
          if( _log ) _log->comment(0,"rinexn","warning - RINEXN system not defined, used "+t_gsys::char2str(_gnsssys));
          break; 
        }
        default : { 
          string lg("warning: not supported satellite system"+line.substr(40,1));
          if( _log ) _log->comment(0,"rinexn", lg); 
          else       cerr << lg << endl;
          // _mutex.unlock(); return -1;
        }
      }
       
//      if( _gnsssys == ' ' ){
//        if( _log ) _log->comment(0,"rinexn","warning: RINEX system not defined, GPS set as default" );
//        _gnsssys = 'G';
//      }

      _version = trim(line.substr(0,9));
      
      if( _log && substitute(_version, " ", "") > 0 ){
        _log->comment(2,"rinexn", "reading VER: " + _version + " SYS: " + string(1,_gnsssys) );
      }

    }else if( line.find("DELTA-UTC",60) != string::npos ){
#ifdef DEBUG
      double a0 = strSci2dbl(line.substr( 3,  19 ));
      double a1 = strSci2dbl(line.substr( 22, 19 ));
      int    w  = static_cast<int>(strSci2dbl(line.substr( 42,  8 )));
      int    t  = static_cast<int>(strSci2dbl(line.substr( 53,  6 )));
      if( _log ) _log->comment(2,"rinexn","reading DELTA-UTC");
      cout << "--- set DELTA-UTC: " << a0 << " : " << a1 << " : " << t << " : " << w << endl;
#endif

    }else if( line.find("PGM / RUN BY / DATE",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexn","reading PGM / RUN BY / DATE");
       
    }else if( line.find("IONOSPHERIC CORR",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexn","reading IONOSPHERIC CORR");
       
    }else if( line.find("TIME SYSTEM CORR",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexn","reading TIME SYSTEM CORR");

    }else if( line.find("CORR TO SYSTEM",60) != string::npos ){      // <= RINEX v2.10
      if( _log ) _log->comment(2,"rinexn","reading CORR TO SYSTEM");

    }else if( line.find("D-UTC A0,A1,T,W,S,U",60) != string::npos ){ // == RINEX v2.11
      if( _log ) _log->comment(2,"rinexn","reading D-UTC A0,A1,T,W,S,U");

    }else if( line.find("ION ALPHA",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexn","reading ION_ALPHA");

    }else if( line.find("ION BETA",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexn","reading ION_BETA");
       
    }else if( line.find("LEAP SECONDS",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexn","reading LEAP SECONDS");
       
    }else if( line.find("COMMENT",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexn","reading COMMENT");
       
    }else if( line.find("END OF HEADER",60) != string::npos ){

      if( _log ) _log->comment(2,"rinexn","reading END OF HEADER ");
      t_gcoder::_consume(tmpsize);
      _mutex.unlock(); return -1;
    }
    t_gcoder::_consume(tmpsize);
  }

  _mutex.unlock(); return consume;
}

   
/* ----------
 * NAV-RINEX body
 */
int t_rinexn::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();
   
  if( t_gcoder::_add2buffer(buff, sz) == 0 ){ _mutex.unlock(); return 0; };

#ifdef DEBUG
  cout << " BUFFER : \n" << _buffer << "\n size = " << sz  << " END OF BUFFER \n\n"; cout.flush();
#endif
  
  t_gtime epoch;
  int b, e, s, l, i;
  int maxrec = MAX_RINEXN_REC;  // implicite
  string timstr;
  t_gnavdata data = { 0 };

  // RINEX v2.xx   
  if( _version[0] == '2' ){        b = 0; e = 22; l = 19; s = 3; // timstr = "%2s %2d %02d %02d %02d %02d %02d";
    // RINEX v3.xx
  }else if( _version[0] == '3' ){  b = 0; e = 23; l = 19; s = 4; // timstr = "%3s %4d %02d %02d %02d %02d %02d";
    // RINEX ???
  }else{                           b = 0; e = 23; l = 19; s = 4; // timstr = "%3s %4d %02d %02d %02d %02d %02d";
  }

  string line;
  int consume = 0;
  int tmpsize = 0;
  int recsize = 0;

  while( ( tmpsize = t_gcoder::_getline( line, 0 ) ) >= 0 ){
    
#ifdef DEBUG
    cout << "0: " << line << endl; cout.flush();
#endif

    consume += tmpsize;
    recsize += tmpsize;
    string epostr = line.substr(b,e);
    
    istringstream istr(line.substr(b,e) );
    istr.clear();

    string prn;
    int yr,mn,dd,hr,mi;
    int svn = 0;
    int irc = 0;
    float sec=0.0;
    char tmpbuff[82];   // RINEX 2 (80) + RINEX 3 (81)
    int min_sz = 23; // minimum size to be decoded (timestamp)
     
    switch( _version[0] ){
      case 2  : min_sz = 23; break;  // RINEX 2
      case 3  : min_sz = 24; break;  // RINEX 3
    }
       
    if( line.size() > 82 || size() <= min_sz ){ // avoid decoding such case
//      cout << "skip [" << size() << ":" << tmpsize <<  ":" << recsize <<  ":" << consume << "]:" << line; cout.flush();
      t_gcoder::_consume(tmpsize);
      recsize = consume = 0; break; // read buffer 
    }

    strncpy(tmpbuff, line.c_str(), min_sz );
    tmpbuff[line.size()] = '\0';
    
    if( _version[0] == '2' ){          // don't use '%2i' in scan, but '%2d' instead !
      irc = sscanf(tmpbuff, "%2i%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %5f", &svn,&yr,&mn,&dd,&hr,&mi,&sec);

      if( irc < 7 ){ // not success - remove this data
//        cout << "fail [" << irc << ":" << size() << ":" << tmpsize <<  ":" << recsize <<  ":" << consume << "]:[" << tmpbuff << "]\n"; cout.flush();
        t_gcoder::_consume(tmpsize);
        recsize = consume = 0; continue;
      }
      prn = t_gsys::eval_sat(svn, t_gsys::char2gsys(_gnsssys));
//        cout << "LINE [" << irc << "]:[" << tmpbuff << "]\n"; cout.flush();

    }else{ 
      int isec = 0; char sat[3 + 1]; sat[3] = '\0';  // don't use '%2i' in scan, but '%2d' instead !
      irc = sscanf(tmpbuff, "%c%2d%*[ ] %4d%*[ ] %2d%*[ ] %2d%*[ ]%2d%*[ ]%2d%*[ ]%2d", &sat[0],&svn,&yr,&mn,&dd,&hr,&mi,&isec);
      
      if( irc < 8 ){ // not success - remove this data
//        cout << "fail [" << irc << ":"  << size() << ":" << tmpsize <<  ":" << recsize <<  ":" << consume << "]:[" << tmpbuff << "]\n"; cout.flush();
        t_gcoder::_consume(tmpsize);
        recsize = consume = 0; continue;
      }
      prn = t_gsys::eval_sat(svn, t_gsys::char2gsys(sat[0]));

      sec = isec;
//      cout << "LINE [" << irc << "]:" << sat[0] << ":" << prn << "]:[" << tmpbuff << "]\n"; cout.flush();

    }
     
//    cout << "LINE IS OK:" << prn << ":" << line; cout.flush();
    shared_ptr<t_gnav> geph = make_shared<t_gnav>();
       
    // !!! SAT musi byt preveden na PRN (pro RINEX ver < 3) --> neni implementovano!
    if(       prn[0] == 'G' ){   maxrec = MAX_RINEXN_REC_GPS;  geph = make_shared<t_gnavgps>();  geph->glog(_log); epoch.tsys(t_gtime::GPS);
    }else if( prn[0] == 'R' ){   maxrec = MAX_RINEXN_REC_GLO;  geph = make_shared<t_gnavglo>();  geph->glog(_log); epoch.tsys(t_gtime::UTC);
    }else if( prn[0] == 'E' ){   maxrec = MAX_RINEXN_REC_GAL;  geph = make_shared<t_gnavgal>();  geph->glog(_log); epoch.tsys(t_gtime::GAL);
    }else if( prn[0] == 'J' ){   maxrec = MAX_RINEXN_REC_QZS;  geph = make_shared<t_gnavqzs>();  geph->glog(_log); epoch.tsys(t_gtime::GPS);   // QZSS is equal to GPS time
    }else if( prn[0] == 'S' ){   maxrec = MAX_RINEXN_REC_SBS;  geph = make_shared<t_gnavsbs>();  geph->glog(_log); epoch.tsys(t_gtime::GPS);   // SBAS is equal to GPS time
    }else if( prn[0] == 'C' ){   maxrec = MAX_RINEXN_REC_BDS;  geph = make_shared<t_gnavbds>();  geph->glog(_log); epoch.tsys(t_gtime::BDS);
    }else if( prn[0] == 'I' ){   maxrec = MAX_RINEXN_REC_IRN;  geph = make_shared<t_gnavirn>();  geph->glog(_log); epoch.tsys(t_gtime::GPS);   // ??
    }else{
      string lg("Warning: not supported satellite satellite system: "+prn); mesg(GWARNING,lg);
      if( _log )  _log->comment(0,"rinexn", lg); 
      else cerr << lg << endl;
    }

    epoch.from_ymdhms(yr, mn, dd, hr, mi, sec);

    // double check for corrupted files
//    if( fabs(epoch - _check_dt) < 7*86400 || _check_dt == FIRST_TIME || !_filter_gnss(prn) ){
//      cout << "dset: " << epoch.str_ymdhms() << " " << _check_dt.str_ymdhms() << " " << line; cout.flush();
//      _check_dt = epoch;
//    }else{
    if( fabs(epoch - _check_dt) > 7*86400 && _check_dt != FIRST_TIME ){
      string lg(prn+" strange epoch or corrupted file "+_fname); mesg(GWARNING,lg);
      if( _filter_gnss(prn) && _log ) _log->comment(0,"rinexn",lg);
//      cout << "date: " << epoch.str_ymdhms() << " " << _check_dt.str_ymdhms() << " " << line; cout.flush();
      t_gcoder::_consume(tmpsize);
      recsize = consume = 0; continue;
    }
     
#ifdef DEBUG       
    cout << "  line   = " << line 
         << "  satell = " << prn 
         << "  maxrec = " << maxrec 
         << epoch.str_ymdhms("  ep = ")
         << endl; cout.flush();
#endif
     
    if( tmpsize < 57+s ) break;

    data[0] = strSci2dbl(line.substr( 19+s, l ));
    data[1] = strSci2dbl(line.substr( 38+s, l ));
    data[2] = strSci2dbl(line.substr( 57+s, l ));

    i = 2;
    while( i < MAX_RINEXN_REC ){

      // incomplete record
      if( ( tmpsize = t_gcoder::_getline( line, recsize ) ) < 0 ){ break; }

      consume += tmpsize;
      recsize += tmpsize;
      if( ++i < maxrec ){ if( tmpsize>   s ) data[i] = strSci2dbl(line.substr(    s, l )); } //cout << "1 => " << data[i] << "\n"; }
      if( ++i < maxrec ){ if( tmpsize>19+s ) data[i] = strSci2dbl(line.substr( 19+s, l )); } //cout << "2 => " << data[i] << "\n"; }
      if( ++i < maxrec ){ if( tmpsize>38+s ) data[i] = strSci2dbl(line.substr( 38+s, l )); } //cout << "3 => " << data[i] << "\n"; }
      if( ++i < maxrec ){ if( tmpsize>57+s ) data[i] = strSci2dbl(line.substr( 57+s, l )); } //cout << "4 => " << data[i] << "\n"; }

#ifdef DEBUG
      cout << "prn = " << prn << " " << prn[0] << "  maxrec = " << maxrec
           << " " << i << epoch.str_ymdhms("  ep = ")
           << " " << data[i-3] << " " << data[i-2] << " " << data[i-1] << " " << data[i-0] << endl;
#endif
      
      // is record complete and filter-out GNSS systems
      if( geph && i+1 >= maxrec ){
        t_gcoder::_consume(recsize);
        recsize = 0;

        // filter GNSS and SAT
        if( !_filter_gnss(prn) ){
          if( _log ) _log->comment( 4, "rinexn", "skip "+prn);
          break;
        }

        if( epoch < _beg || epoch > _end ){
          if( _log ) _log->comment( 4, "rinexn", "skip "+prn);
          break;
	}

        geph->data2nav( prn, epoch, data );
        geph->gio(_gio_ptr.lock());
	 
	// reset check_dt (already filtered sat)
        if( geph->healthy() ){ _check_dt = geph->epoch(); };

        if( !_filter_gnav(geph, prn) ){
          if( _log ) _log->comment( 4, "rinexn", "skip navigation type: prn " + prn + " iod " + int2str(geph->iod()));
          break;
        }

        // collect 1-line messages
        // mesg(GWARNING,geph->linefmt());
	if( _log ) _log->comment(-3,geph->linefmt());
 
        // fill data
        map<string,t_gdata*>::iterator it = _data.begin();
        while( it != _data.end() ){
          
          if( it->second->id_type()  == t_gdata::ALLNAV ||
              it->second->id_group() == t_gdata::GRP_EPHEM){

            ((t_gallnav*)it->second)->add( geph );
            
          }else if( it->second->id_type() == t_gdata::ALLOBJ ){
            t_gallobj* all_obj = (t_gallobj*)it->second;
            shared_ptr<t_gobj>  one_obj = all_obj->obj(geph->sat());
            if( one_obj != 0 ) one_obj->glog(_log);	       	          
            if( geph->id_type() == t_gdata::EPHGLO ){
              int ch = dynamic_pointer_cast<t_gnavglo>(geph)->channel();
              if( one_obj != 0 ) one_obj->channel(ch);
            }
          }
          it++;
        }
        cnt++;
        break;
      }
    }
    if( recsize != 0 ) break; // break the initialization loop if read not finished correctly
  }
  _mutex.unlock(); return consume;
}


/* -------- *
 * ENCODER
 * -------- */

/* ----------
 * NAV-RINEX header
 */
int t_rinexn::encode_head( char* buff, int sz, vector<string>& errmsg)
{ 

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();

  t_gtime tt(t_gtime::UTC);
  t_gtime epoch(t_gtime::GPS);

  // check if data exists
  int count = 0;
  for( auto it = _data.begin(); it != _data.end(); ++it ){
    if( it->second->id_type()  == t_gdata::ALLNAV ||
        it->second->id_group() == t_gdata::GRP_EPHEM )
    {
      count += dynamic_cast<t_gallnav*>(it->second)->vec_nav().size();
      epoch  = dynamic_cast<t_gallnav*>(it->second)->beg_gnav("");
    }
  }

  if( _ss_position == 0 && count > 0 ){ // only if data

    if( _version < "3.00" ){
       
       // RINEX 2.x (GNSS-specific)
      if( _sys.size() != 1 ){
        if( _log ) _log->comment(0,"rinexn2","Warning: NAV no encoded, more than one GNSS system identified!");
        _mutex.unlock(); return -1;
      }
      set<string>::iterator itSYS = _sys.begin();
      switch( t_gsys::str2gsys(*itSYS) ){
        case GPS : _ss << "     2.11           N: GPS NAV DATA                         RINEX VERSION / TYPE" << endl; break;
        case GLO : _ss << "     2.11           G: GLO NAV DATA                         RINEX VERSION / TYPE" << endl; break;
        case GAL : _ss << "     2.12           N: GNSS NAV DATA    E: GAL NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case BDS : _ss << "     2.12           N: GNSS NAV DATA    C: BDS NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case QZS : _ss << "     2.12           N: GNSS NAV DATA    Q: QZS NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case SBS : _ss << "     2.12           N: GNSS NAV DATA    S: SBS NAV DATA     RINEX VERSION / TYPE" << endl; break;
        case IRN : _ss << "     2.12           N: GNSS NAV DATA    I: IRN NAV DATA     RINEX VERSION / TYPE" << endl; break;
        default :  _ss << "     2.12           N: GPS NAV DATA                         RINEX VERSION / TYPE" << endl; break;
      }
    }else{ // RINEX 3.x (GNSS-mixed)
                   _ss << "     3.02           N: GNSS NAV DATA    M (MIXED)           RINEX VERSION / TYPE" << endl;
  }

  //  << "G-Nut/Anubis        GOP/RIGTC           20120618 000000 UTC PGM / RUN BY / DATE"  << endl
  _ss << left << setw(20) << "G-Nut/"+ _pgm << setw(20) << "GOP/RIGTC"
      << tt.str("%Y%m%d %H%M%S UTC")                             << " PGM / RUN BY / DATE " << endl
      << right << setw(6) << epoch.leapsec() - TAI_GPS
      <<       "                                                      LEAP SECONDS        " << endl
      << "                                                            END OF HEADER       " << endl;
  }
   
  int size = _fill_buffer( buff, sz );

  _mutex.unlock(); return size;
}


/* ----------
 * NAV-RINEX data
 */
int t_rinexn::encode_data( char* buff, int sz, int& cnt, vector<string>& errmsg )
{ 

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _mutex.lock();

  // RINEX 2.x (GNSS-specific)
  if( _version < "3.00" && _sys.size() != 1 ){
    if( _log ) _log->comment(0,"rinexn2","Warning: NAV no encoded, more than one GNSS system identified!");
    _mutex.unlock(); return -1;
  }

  map<string,int> count_mesg;

  // fill data
  if( _ss_position == 0 ){

    map<string,t_gdata*>::iterator it = _data.begin();
    for( it = _data.begin(); it != _data.end(); ++it ){

      if( it->second->id_type()  == t_gdata::ALLNAV ||
          it->second->id_group() == t_gdata::GRP_EPHEM ){

        vector<shared_ptr<t_geph>> vec_nav = dynamic_cast<t_gallnav*>(it->second)->vec_nav();
        vector<shared_ptr<t_geph>>::iterator itNAV;
        for( itNAV = vec_nav.begin(); itNAV != vec_nav.end(); ++itNAV )
        {
          shared_ptr<t_gnav> gnav = dynamic_pointer_cast<t_gnav>(*itNAV);

          if( gnav->epoch() < _beg || gnav->epoch() > _end ){ continue; } // filter BEG/END

#ifdef DEBUG
  cout << "NAV: " << gnav->sat() << " " << gnav->epoch().str_ymdhms(" ")
	                         << " " << gnav->valid() << endl;
#endif         
          count_mesg[gnav->sat()] += 1;

          t_gnavdata data;
          if( gnav->nav2data( data ) == 0 ){

            if( _version < "3.00" ){
              _ss << setw(2) << str2int(gnav->sat().substr(1,2))
              << " " << gnav->epoch().str("%y %m %d %H %M") // y-digit year, float seconds
              << " " << fixed << setprecision(1) << setw(4) << gnav->epoch().secs() + gnav->epoch().dsec();
            }else{
              _ss << gnav->sat() + " " + gnav->epoch().str("%Y %m %d %H %M %S");  // 4-digit year !!!
            }
            _ss  << fixed << scientific << setprecision(12); // setfill('0')
	     
            for( int i=0; i<gnav->rec(); ++i ){
              if( (i-3)%4 == 0 ){
                
                if( _version < "3.00" )
                _ss << endl << setw(3) << ""; // RINEX 2: 3-CH space before first column first
                else _ss << endl << setw(4) << ""; // RINEX 3: 4-CH space before first column first
                
              }
              _ss << setw(19) << data[i];
            }
            _ss << endl;
          }
        }
      }
    }
  }

  int size = _fill_buffer( buff, sz );

  for( auto it = count_mesg.begin(); it != count_mesg.end(); ++it ){
    if( _log ) _log->comment(1,"rinexn","#NAV "+it->first+": "+int2str(it->second));
  }

  _mutex.unlock(); return size;
}

// settings
// ----------
//int t_rinexn::_gset(t_gsetbase* set)
//{ 
//
//   if( t_gcoder::_gset(set) != 0 ) return -1;
//
//  _check = dynamic_cast<t_gsetinp*>(_set)->chkNav();      // simple check
//   
//  return 0;
//}   
   
// filter out navigation mess. types
// ---------   
bool t_rinexn::_filter_gnav(shared_ptr<t_gnav> geph, const string& prn)
{
  gtrace("t_rinexn::_filter_gnav");
   
  bool ret = true; 
   
  GSYS gs = t_gsys::char2gsys(prn[0]);

   
  // currently only GALILEO can be filtered!
  if( gs == GAL ){
      
    GNAVTYPE gnav = dynamic_pointer_cast<t_gnavgal>(geph)->gnavtype();

    if( _nav[gs].size() == 0 ||
        _nav[gs].find(gnavtype2str(gnav)) != _nav[gs].end()){ ret = true;  }
    else{                                                      ret = false; }
  }
   
  return ret; // ALL OTHER SYSTEMS
}
   
} // namespace
