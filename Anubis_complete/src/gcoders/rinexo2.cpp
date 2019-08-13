
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
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
#include <algorithm>

#include "gcoders/rinexo2.h"
#include "gutils/gfileconv.h"
#include "gutils/gsysconv.h"
#include "gdata/grec.h"
#include "gdata/gobsgnss.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_rinexo2::t_rinexo2( t_gsetbase* s, string version, int sz )
  : t_gcoder( s, version, sz ),
  _site(""),
  _epo_beg(LAST_TIME),
  _epo_end(FIRST_TIME),
  _xbeg(0), _xend(0), _xsmp(0),
  _xsys(0) //,  _xsat(0), _xsig(0)
{
  this->clear();
}

   
// destructor
// ----------
t_rinexo2::~t_rinexo2(){

  // fill statistics of filtered data
  map<string,t_gdata*>::iterator itDAT = _data.begin();

  while( itDAT != _data.end() ){
    if( itDAT->second->id_type() != t_gdata::ALLOBS ){
      itDAT++; continue;
    }
     
    if( _log && _log->verb() >= 2 )
        _log->comment(2,"rinexo", _fname +
		        " [filtered data] xBEG:" + int2str(_xbeg) + " xEND:" + int2str(_xend)
	                              + " xSMP:" + int2str(_xsmp) + " xSYS:" + int2str(_xsys));

    t_gallobs::t_xfilter xflt;

    xflt.xdat[t_gallobs::XDATA_BEG] = _xbeg;
    xflt.xdat[t_gallobs::XDATA_END] = _xend;
    xflt.xdat[t_gallobs::XDATA_SMP] = _xsmp;
    xflt.xdat[t_gallobs::XDATA_SYS] = _xsys;

    xflt.beg = _epo_beg;
    xflt.end = _epo_end;
     
    ((t_gallobs*)itDAT->second)->xdata( _rnxhdr.markname(), _fname, xflt );
              
    itDAT++;  
  }
};


// clear
// ----------
void t_rinexo2::clear()
{ 
  t_gcoder::clear();

  _epoch    = t_gtime::current_time(t_gtime::GPS);
   
  _rnxhdr.clear();
   
  _line     = "";
  _flag     = '0';
  _nsat     = 0;
  _count    = 0;
  _tmpsize  = 0;
  _consume  = 0;
  _complete = true;
  _csys     = " ";
  _xbeg = _xend = _xsmp = _xsys = 0;
//_xsat = _xsig = 0;
   
  _pcosat.erase(_pcosat.begin(),_pcosat.end());
  _pcosys.erase(_pcosys.begin(),_pcosys.end());
  _pcoecc.erase(_pcoecc.begin(),_pcoecc.end());
  _mapobs.erase(_mapobs.begin(),_mapobs.end());

}

   
// ----------
// OBS-RINEX header
//   - read individual lines (tmpsize once used and consumed)
//   - read block of lines at once (tmpsize cummulated and then consumed)
// ----------
int t_rinexo2::_decode_head()
{ 
  gtrace("t_rinexo2::_decode_head");

  // -------- "PGM / RUN BY / DATE" --------
  if( _line.find("PGM / RUN BY / DATE",60) != string::npos ){
    _rnxhdr.program( trim(_line.substr( 0,20)) );
    _rnxhdr.runby  ( trim(_line.substr(20,20)) );
    t_gtime gtime(t_gtime::UTC);
    if( _line.substr(56,3) != "UTC" ) gtime.tsys(t_gtime::LOC);
       
    if(      gtime.from_str("%Y%m%d %H%M%S",    _line.substr(40,15)) == 0 ) ;
    else if( gtime.from_str("%Y-%m-%d %H-%M-%S",_line.substr(40,20)) == 0 ) ;
    else{    gtime = FIRST_TIME; }
    _rnxhdr.gtime(gtime); 
       
    if( _log && _log->verb() >= 2)
        _log->comment(2,"rinexo2","PGM / RUN BY / DATE: " + _rnxhdr.program()
		                                    + " " + _rnxhdr.runby() 
		                                    + " " + _rnxhdr.gtime().str_ymdhms() );


  // -------- "COMMENT" --------
  }else if( _line.substr(60,7).find("COMMENT") != string::npos ){
       
    _comment.push_back(_line.substr(0,60));
    if( _log ) _log->comment(2,"rinexo2","COMMENT");

     
  // -------- "MARKER NAME" --------
  }else if( _line.substr(60,11).find("MARKER NAME") != string::npos ){
    
    string guess, marker, src, typ, sit; // select   short (4-CHAR) or long (9-CHAR)

    // check sitename according to filename (if available)
    string filename = base_name( _fname );

    // first, check the MARKER NAME (RINEX HEADER)
    int i2; char c4[4],c3[3];
    const int size = 9; char tmpbuff[size+1];
    strncpy(tmpbuff, _line.substr(0,size).c_str(), size );  tmpbuff[size] = '\0';
    
    if( sscanf(tmpbuff, "%4c%2d%3c",c4,&i2,c3) == 3 ){ // guess the length (4-char vs 9-char)  %* ==> not robust!
          marker = guess = _line.substr(0,9); src="MARKER"; typ="long";  sit="(9-CH)"; } // LONG  SITE NAME
    else{ marker = guess = _line.substr(0,4); src="MARKER"; typ="short"; sit="(4-CH)"; } // SHORT SITE NAME
    transform(marker.begin(), marker.end(), marker.begin(), ::toupper);
    
    // second, check the FILE NAME
    size_t n = count(filename.begin(), filename.end(), '_');
    if( ! _fname.empty() ){
      if( _fname.length() > 8+4 && n >= 5 ){
            guess = filename.substr(0,9); src="FILE"; typ="long";  sit="(9-CH)"; } // LONG  FILE NAME
      else{ guess = filename.substr(0,4); src="FILE"; typ="short"; sit="(4-CH)"; } // SHORT FILE NAME
      transform(guess.begin(), guess.end(), guess.begin(), ::toupper);

      if( guess != marker ){
        if(_log ) _log->comment(1,"rinexo",src+": "+typ+" site name: "+guess+" "+sit );
        if( guess.substr(0,4) !=  marker.substr(0,4) ){
          mesg(GWARNING,"Warning: 4-CHAR: MARKER_NAME ["+marker+"] not equal to site from filename ["+guess+"]");
        }
      }
    }

    transform(guess.begin(),  guess.end(),  guess.begin(),  ::toupper);
    if(_log) _log->comment(1,"rinexo",src+": "+typ+" site name: "+guess+" "+sit+" identified" );

    // resolve situation of 4-CH MARKER with 9-CH in settings
    if( guess.length() == 4 ){
      for( auto itREC = _rec.begin(); itREC != _rec.end(); ++itREC ){
        if( itREC->length() > 4 && itREC->compare(0,4,guess,0,4) == 0 ){
          if(_log) _log->comment(1,"rinexo2","4-CH ["+guess+"] found in RINEX, but user requests 9-CH ["+(*itREC)+"]");
          guess = *itREC;
        }
      }
    }
    // resolve opposite situation (9-CH MARKER with 4-CH in settings)
    else if( guess.length() == 9 ){
      for( auto itREC = _rec.begin(); itREC != _rec.end(); ++itREC ){
        if( itREC->length() == 4 && itREC->compare(0,4,guess,0,4) == 0 ){
          if(_log) _log->comment(1,"rinexo2","9-CH ["+guess+"] found in RINEX, but user requests 4-CH ["+(*itREC)+"]");
          guess = *itREC;
        }
      }
    }
    _site = guess;
    _rnxhdr.markname( trim(marker) );             // capitals, cutted according to the station name guess
//  _rnxhdr.markname( trim(_line.substr(0,60)) ); // original text in MARKER_NAME
     

  // -------- "MARKER NUMBER" --------
  }else if( _line.substr(60,13).find("MARKER NUMBER") != string::npos ){
       
    _rnxhdr.marknumb( trim(_line.substr(0,60)) );
    if( _log ) _log->comment(2,"rinexo2","MARKER NUMB: " + _rnxhdr.marknumb());

     
  // -------- "MARKER TYPE" --------
  }else if( _line.substr(60,11).find("MARKER TYPE") != string::npos ){

    _rnxhdr.marktype( trim(_line.substr(0,20)) );
    if( _log ) _log->comment(2,"rinexo2","MARKER TYPE: " + _rnxhdr.marktype());

     
  // -------- "OBSERVER / AGENCY" --------
  }else if( _line.substr(60,17).find("OBSERVER / AGENCY") != string::npos ){

    _rnxhdr.observer( trim(_line.substr( 0,20)) );
    _rnxhdr.agency(   trim(_line.substr(20,40)) );
    if( _log ) _log->comment(2,"rinexo2","OBSERVER / AGENCY: " 
               + _rnxhdr.observer() + " " + _rnxhdr.agency());

     
  // -------- "REC # / TYPE / VERS" --------
  }else if( _line.substr(60,19).find("REC # / TYPE / VERS") != string::npos ){

    _rnxhdr.recnumb( trim(_line.substr(0, 20)) );
    _rnxhdr.rectype( trim(_line.substr(20,20)) );
    _rnxhdr.recvers( trim(_line.substr(40,20)) );
    if( _log ) _log->comment(2,"rinexo2","REC # / TYPE / VERS: "
  	       + _rnxhdr.recnumb() + " " + _rnxhdr.rectype() + " " + _rnxhdr.recvers());

     
  // -------- "ANT # / TYPE" --------
  }else if( _line.substr(60,12).find("ANT # / TYPE") != string::npos ){

    _rnxhdr.antnumb( trim(_line.substr(0, 20)) );
    _rnxhdr.anttype( trim(_line.substr(20,20)) );
    if( _log ) _log->comment(2,"rinexo2","ANT # / TYPE: "
	       + _rnxhdr.antnumb() + " " + _rnxhdr.anttype());

     
  // -------- "APPROX POSITION XYZ" --------
  }else if( _line.substr(60,19).find("APPROX POSITION XYZ") != string::npos ){

    t_gtriple apr; 
    apr[0] = str2dbl(_line.substr( 0,14));
    apr[1] = str2dbl(_line.substr(14,14));
    apr[2] = str2dbl(_line.substr(28,14));

    if( _log ) _log->comment(2,"rinexo2","APPROX POSITION XYZ: " 
               + dbl2str(apr[0]) + " " + dbl2str(apr[1]) + " " + dbl2str(apr[2]));
    _rnxhdr.aprxyz(apr);

     
  // -------- "ANTENNA: DELTA X/Y/Z --------
  }else if( _line.substr(60,20).find("ANTENNA: DELTA X/Y/Z") != string::npos ){

    t_gtriple ecc; 
    ecc[0] = str2dbl(_line.substr( 0,14));         
    ecc[1] = str2dbl(_line.substr(14,14)); 
    ecc[2] = str2dbl(_line.substr(28,14));
    if( _log ) _log->comment(2,"rinexo2","ANTENNA: DELTA X/Y/Z: "
               + dbl2str(ecc[0]) + " " + dbl2str(ecc[1]) + " " + dbl2str(ecc[2]));
    _rnxhdr.antxyz(ecc);
 
     
  // -------- "ANTENNA: DELTA H/E/N" --------
  }else if( _line.substr(60,20).find("ANTENNA: DELTA H/E/N") != string::npos ){
       
    t_gtriple ecc; 
    ecc[1] = str2dbl(_line.substr(14,14));  // East  // --> REVERSE with N !
    ecc[0] = str2dbl(_line.substr(28,14));  // North // --> REVERSE with E !
    ecc[2] = str2dbl(_line.substr( 0,14));  // Up
    if( _log ) _log->comment(2,"rinexo2","ANTENNA: DELTA H/E/N: "
               + dbl2str(ecc[0]) + " " + dbl2str(ecc[1]) + " " + dbl2str(ecc[2]));
    _rnxhdr.antneu(ecc);
 		if( !ecc.zero() &&  _rnxhdr.antxyz().zero() ){
			t_gtriple refell;  xyz2ell(_rnxhdr.aprxyz(), refell, false);
			t_gtriple antxyz;  neu2xyz(refell, ecc, antxyz);
      _rnxhdr.antxyz(antxyz);
		}

     
  // -------- "ANTENNA: PHASECENTER" --------
  }else if( _line.substr(60,20).find("ANTENNA: PHASECENTER") != string::npos ){

    t_gtriple ecc;
    ecc[0] = str2dbl(_line.substr( 5, 9));      // N or X
    ecc[1] = str2dbl(_line.substr(14,14));      // E or Y
    ecc[2] = str2dbl(_line.substr(28,14));      // U or Z
       
    _pcosat.push_back(_line[0]);                // GNSS (G/R/E/S/...)
    _pcosys.push_back(trim(_line.substr(2,3))); // coordinate system (NEU/XYZ)
    _pcoecc.push_back(ecc);                     // NEU or XYZ phase center vs. ARP eccentricities

    if( _log ) _log->comment(2,"rinexo2","ANTENNA: PHASECENTER: ");
//                     + _pcosat[_pcosat.size()-1]
//		 + " " + _pcosys[_pcsoys.size()-1]
//               + " " +  dbl2str(ecc[0])
//		 + " " +  dbl2str(ecc[1]) 
//		 + " " +  dbl2str(ecc[2]));

     
  // -------- "ANTENNA: B.SIGHT XYZ" --------
  }else if( _line.substr(60,20).find("ANTENNA: B.SIGHT XYZ") != string::npos ){

//    _recnumb  = trim(_line.substr(0, 19));
//    _rectype  = trim(_line.substr(20,19));
//    _recvers  = trim(_line.substr(40,19));
    if( _log ) _log->comment(2,"rinexo2","ANTENNA: B.SIGHT XYZ: " + string("not implemented"));
//             + _recnumb + " " + _rectype + " " + _recvers);

       
  // -------- "ANTENNA: ZERODIR AZI --------
  }else if( _line.substr(60,20).find("ANTENNA: ZERODIR AZI") != string::npos ){
    if( _log ) _log->comment(2,"rinexo2","ANTENNA: ZERODIR AZI: "   + string("not implemented"));


  // -------- "ANTENNA: ZERODIR XYZ" --------
  }else if( _line.substr(60,20).find("ANTENNA: ZERODIR XYZ") != string::npos ){
    if( _log ) _log->comment(2,"rinexo2","ANTENNA: ZERODIR XYZ: "   + string("not implemented"));


  // -------- "CENTER OF MASS: XYZ" --------
  }else if( _line.substr(60,19).find("CENTER OF MASS: XYZ") != string::npos ){
    if( _log ) _log->comment(2,"rinexo2","CENTER OF MASS: XYZ: "   + string("not implemented"));


  // -------- "# / TYPES OF OBSERV" --------
  }else if( _line.substr(60,19).find("# / TYPES OF OBSERV") != string::npos ){
       
    int fact = 1;
//  char sys( string(1,_rnxhdr.rnxsys());
    if( ! _csys.empty() && _csys[0] != ' ' ){
	 
      // guarantee that all relevant lines are available
      int num = str2int(_line.substr(0,6));
      int nlines = (int)ceil(num/9.0);         // upto 9 values/line
      int addsize = 0;
      string tmp;
      for( int nl = 0; nl < nlines;  ++nl ){
        int irc = t_gcoder::_getline(tmp, _tmpsize + addsize);
        if( irc >= 0 ){ 
          addsize += irc;
        }else{ 
          // cout << "INCOMPLETE RECORD(OBS TYPES): " << tmp << "\n";
          _complete = false;
          break;
        }	 
      }

      if( ! _complete ) return -1;

      // sufficient number of lines is already guaranteed here
      for(int ii=0; ii<num; ++ii){
        // make sure that the obs is not already in t_vobstypes
        string go = _line.substr(10 + 6*(ii%(9)), 2);
        GOBS obs = str2gobs( go );

        t_rnxhdr::t_vobstypes::iterator it = _mapobs[_csys].begin();

        // record exists
        for( it = _mapobs[_csys].begin(); it != _mapobs[_csys].end(); ++it ){
          if( it->first == obs ){ it->second = fact; break; }
        }
        
        // new record
        if( it == _mapobs[_csys].end() ) _mapobs[_csys].push_back( make_pair(obs, fact) );
        
        if( (ii+1)%9 == 0 && num > 9) // get newline and check (read next for old system)
        _tmpsize += t_gcoder::_getline(_line, _tmpsize );
      } // loop over # observations
      
      if( _log && _log->verb() >= 2 ){
        ostringstream lg;
        lg << "# / TYPES OF OBSERV: [" << _csys << "] ";
	
        t_rnxhdr::t_vobstypes::const_iterator it = _mapobs[_csys].begin();
        while( it != _mapobs[_csys].end() ){ 
          lg << " " << gobs2str( it->first );
          it++;	 
        }
        _log->comment(2,"rinexo2",lg.str());
      }
    }


  // -------- "INTERVAL" -------- 
  }else if( _line.substr(60,8).find("INTERVAL") != string::npos ){
       
    string sample = _line.substr(0,9);
    _rnxhdr.interval( str2dbl(sample) );
    if( _log ) _log->comment(2,"rinexo2","INTERVAL: " + dbl2str(_rnxhdr.interval()));

    if( _int > 1.0 && _rnxhdr.interval() < 1.0 && _scl == 1 ){

      int _dec = sample.substr(sample.find(".")+1).length()-1; // decimal digits resolution
      _scl = (int)pow(10, _dec);
      cerr << "#Notice [rinexo2]: use decimals in XML settings if high-rate file."
           << " Used <int> " << fixed << setprecision(_dec) << _int << " </int>\n";
    }

     
  // -------- "TIME OF FIRST OBS" -------- 
  }else if( _line.substr(60,17).find("TIME OF FIRST OBS") != string::npos ){

    string tsys = trim(_line.substr(48,3));
    t_gtime first; 
       
    if(      tsys == "GPS" )  first.tsys(t_gtime::GPS);
    else if( tsys == "GLO" )  first.tsys(t_gtime::GLO);
    else if( tsys == "GAL" )  first.tsys(t_gtime::GAL);
    else                      first.tsys(t_gtime::GPS);
       
    first.from_ymdhms(str2int(_line.substr( 0, 6)),        // year
                      str2int(_line.substr( 6, 6)),        // month
                      str2int(_line.substr(12, 6)),        // day
                      str2int(_line.substr(18, 6)),        // hours
                      str2int(_line.substr(24, 6)),        // minutes
                      str2dbl(_line.substr(30,13)) );      // double seconds

    if( _log  && _log->verb() >= 2 )
       _log->comment(2,"rinexo2",first.str_ymdhms("TIME OF FIRST OBS: ") );

    _rnxhdr.first(first); 


  // -------- "TIME OF LAST OBS" -------- 
  }else if( _line.substr(60,16).find("TIME OF LAST OBS") != string::npos ){

    string tsys = trim(_line.substr(48,3));
    t_gtime last;
       
    if(      tsys == "GPS" )  last.tsys(t_gtime::GPS);
    else if( tsys == "GLO" )  last.tsys(t_gtime::GLO);
    else if( tsys == "GAL" )  last.tsys(t_gtime::GAL);
    else                      last.tsys(t_gtime::GPS);

    last.from_ymdhms(str2int(_line.substr( 0, 6)),        // year
                     str2int(_line.substr( 6, 6)),        // month
                     str2int(_line.substr(12, 6)),        // day
                     str2int(_line.substr(18, 6)),        // hours
                     str2int(_line.substr(24, 6)),        // minutes
                     str2dbl(_line.substr(30,13)) );      // double seconds

    if( _log && _log->verb() >= 2 ) _log->comment(2,"rinexo2",last.str_ymdhms("TIME OF LAST OBS: ") );
    _rnxhdr.last(last);   
     
  // -------- "LEAP SECONDS" -------- 
  }else if( _line.substr(60,12).find("LEAP SECONDS") != string::npos ){  
       
    // implemented only 1st part - total LEAP SEC.
    _rnxhdr.leapsec( str2int(_line.substr(0,6)) );
    if( _log ) _log->comment(2,"rinexo2","LEAP SECONDS: " + int2str(_rnxhdr.leapsec()));


  // -------- "# OF SATELLITES" -------- 
  }else if( _line.substr(60,15).find("# OF SATELLITES") != string::npos ){
       
    _rnxhdr.numsats( str2int(_line.substr(0,6)) );
    if( _log ) _log->comment(2,"rinexo2","# OF SATELLITES: " + int2str(_rnxhdr.numsats()));


  // -------- "PRN / # OF OBS -------- 
  }else if( _line.substr(60,14).find("PRN / # OF OBS") != string::npos ){
    if( _log ) _log->comment(2,"rinexo2","PRN / # OF OBS: "         + string("not implemented"));

  }

  return 1;
}

   
// ----------
// OBS-RINEX body
//   - read block of lines for each epoch (tmpsize cummulated and then consumed)
// ----------
int t_rinexo2::_decode_data()
{  
  gtrace("t_rinexo2::_decode_data");

  vector<string> satellites;
  int irc = t_rinexo2::_read_epoch();

//if( _flag != '0' && _flag != '1' ){ // special event handling
//  if( _log ) _log->comment(2,"rinexo2",_epoch.str_ymdhms("RINEX 2.x special events not handled after "));
//  else                      cerr << _epoch.str_ymdhms("RINEX 2.x special events not handled after ") << "\n";

  // not ready to read
  if( irc < 0 ){ _complete = false; return -1; }

  // skip special lines
  if( irc > 0 ){

    int xlines = 1;           // special event,      skip only number of lines     
    if( irc >= 2 ){           // skip observations,  skip x-nsat lines
      if( _read_satvec( satellites ) < 0 ){ _complete = false; return -1; }
      xlines  = (int)floor((_mapobs[_csys].size()-1)/5)+1;
    }
//      cout << "skipping: " << xlines << " -> " << _nsat*xlines << " lines  at:" << _line;
    for( int i = 0; i < _nsat*xlines; ++i ){
      int sz = t_gcoder::_getline( _line, _tmpsize );
//      cout << " line:" << _line;
      if( sz <= 0 ){ _complete = false; return -1; }
      _tmpsize += sz;
    }
//  cout << "skip finished: " << _nsat*xlines << " lines  at:" << _line;
     
    t_gcoder::_consume(_tmpsize);
    _consume += _tmpsize;
    _tmpsize = 0;
    return 0;
  }


  // DATA READING 
  // ------------
  // read data for vector for satellites (RINEX v2.x specific)
  if( _read_satvec( satellites ) < 0 ){ _complete = false; return -1; }

  // loop over satellits records (x-lines)
  for( int i = 0; i < _nsat; i++ ){

    string sat = satellites[i];	   
    string key = sat;

    // read sat-line and define system/satellite key in _mapobs (RINEX v2.x specific)
    if( t_rinexo2::_read_syskey( sat, key ) < 0 ){ _complete = false; return -1; }

    // loop over all observation types for system/satellite key in _mapobs
    if( t_rinexo2::_read_obstypes( sat, key ) < 0 ){ _complete = false; return -1; }	 
  }

  // fill data - block (epoch-wise) is completed
  if( _complete ){

    if( _log && _log->verb() >= 5 ) _log->comment(5,"rinexo2",_epoch.str_ymdhms("reading finished at "));
    _count   += this->_fill_data();
     t_gcoder::_consume(_tmpsize);
    _consume += _tmpsize;
    _tmpsize = 0;
  }

  return 0;
}


// read sat-system key (G,E,R,.. M, or specific satellite for _mapobs)
// only RINEX 2.x version
// ----------
int t_rinexo2::_read_syskey( string& sat, string& key )
{
  gtrace("t_rinexo2::_read_syskey");

  // get newline and check the end of buffer	
  int addsize = t_gcoder::_getline( _line, _tmpsize );

  // getline succeed only if 'EOL' found
  if( _line.length() <= 0 ) return _stop_read();
   
  // select proper key (sat/sys) for observation list    
  if(      _mapobs.count(sat.substr(0,3)) != 0 ) key = sat.substr(0,3);
  else if( _mapobs.count(sat.substr(0,1)) != 0 ) key = sat.substr(0,1);
  else                                           key = _csys;          // string(1,_rnxhdr.rnxsys());

  if( key[0] != ' ' && key[0] != 'M' ) sat = t_gsys::eval_sat( sat, t_gsys::str2gsys(key) );

#ifdef DEBUG
  cout << " _READ_SYSKEY :"  << sat << "|" << sat.substr(0,1) 
       << "[" << key << "] " << " [" << _csys << "]\n";  // _rnxhdr.rnxsys()

  t_rnxhdr::t_obstypes::iterator it = _mapobs.begin();
  cout << "SYSKEY[" << _mapobs.size() << "]:";
  while( it != _mapobs.end() ){
    cout << ":" << it->first;
    it++;
  }
  cout << endl;
#endif

  // check if key exists
  if( _mapobs.find(key) == _mapobs.end() ){
    cerr << "rinexo2: Warning: mapobs syskey not found ["+key+"]\n";
       mesg(GWARNING,"Warning: mapobs syskey not found ["+key+"]");
    return _stop_read();
  }   

  // check if enough data (+1 EOL)
  if( addsize <= (int)_mapobs[key].size()*16 + 1 ){ _tmpsize += addsize;
  }else{ 
    if( _log ) _log->comment(0,"rinexo2", _epoch.str_ymdhms(sat + " data line probably incomplete: "));
                 mesg(GWARNING,"rinexo2:"+_epoch.str_ymdhms(sat + " data line probably incomplete: "));
    return _stop_read();
  }
  return 1;
}
   

// read epoch & number of satellites, return flag
// ----------
int t_rinexo2::_read_epoch()
{
  gtrace("t_rinexo2::_read_epoch");

  _nsat =  0;
  _flag = 'X';
   
  // cout << " epoch line  :" << _line ;
   
  // not enough data to recognize standard epoch or special event
  if( _line.length() < 32 ) return _stop_read();
  _flag = _line[28];
  _nsat = str2int( _line.substr(29,3) );
  if( _line.substr(0,6) == "      " ){
    switch ( _flag ){
     case '2' : 
        if( _log ) _log->comment(2,"rinexo2","Warning: start moving antenna identified, but not yet implemented!");
                               mesg(GWARNING,"Warning: start moving antenna identified, but not yet implemented!");
      //else              cerr <<  "rinexo2 - Warning: start moving antenna identified, but not yet implemented!\n";
        return 1;
        
     case '3' :
        if( _log ) _log->comment(2,"rinexo2","Warning: new site occupation identified, but not yet implemented!");
                               mesg(GWARNING,"Warning: new site occupation identified, but not yet implemented!");
      //else              cerr <<  "rinexo2 - Warning: new site occupation identified, but not yet implemented!\n";
        return 1;

     case '4' :
        if( _log ) _log->comment(2,"rinexo2","Warning: new header information identified, but not fully implemented!");
                               mesg(GWARNING,"Warning: new header information identified, but not fully implemented!");
      //else              cerr <<  "rinexo2 - Warning: new header information identified, but not fully implemented!\n"; 
        return 1;
 
     case '5' :
        if( _log ) _log->comment(2,"rinexo2","Warning: external event identified, but not yet implemented!");
                               mesg(GWARNING,"Warning: external event identified, but not yet implemented!");
       //else             cerr <<  "rinexo2 - Warning: external event identified, but not yet implemented!\n";
        return 1;
    }
  }	

  istringstream is(_line);
  is.clear();

  int yr,mn,dd,hr,mi;
  double sc;

  is >> yr >> mn >> dd >> hr >> mi >> sc >> _flag >> _nsat;

  // check success
  if( is.fail() ){ // cout << "rinexo: problem reading EPOCH:" << _line;    
    return _stop_read();
  }

  _epoch.from_ymdhms(yr, mn, dd, hr, mi, sc);
   
  // filter out data
  if( _flag != '0' &&
      _flag != '1'  ){                                                      return  1; }
  if( _epoch > _end ){ if( _epoch > _epo_end ){_epo_end = _epoch;} _xend++; return  2; }
  if( _epoch < _beg ){ if( _epoch < _epo_beg ){_epo_beg = _epoch;} _xbeg++; return  2; }

  if( !_filter_epoch(_epoch) ){ _xsmp++; return 2; }
   
#ifdef DEBUG	   
  if( _log ) _log->comment(0,"rinexo2","reading epoch: " + _epoch.str_ymd() + " " + _epoch.str_hms());
#endif

  return 0;
}


// read satellite list (vector)
// ----------
int t_rinexo2::_read_satvec( vector<string>& satellites )
{  
  gtrace("t_rinexo2::_read_satvec");
 
  int ii = 0;
  for( int i = 0; i < _nsat; i++, ii++ ){

    // getline succeed only if 'EOL' found
    if( (int)_line.length() <= 0 ) return _stop_read();

    // read extra line for each satellite (if complete)
    unsigned int idx = 32 + ii*3;
        
    if( ii > 11 ){
      ii = 0; idx = 32;
      unsigned int addsize = t_gcoder::_getline( _line, _tmpsize );
	 
      if( addsize < idx + 3 ) return _stop_read();
      else                    _tmpsize += addsize;
    }

    if( _line.length() >= idx + 3 ){

      string sat = t_gsys::eval_sat( _line.substr(idx,3) );
      satellites.push_back(sat);
    }
  }

  return _tmpsize;
}


// read satellite observation types
// ----------
int t_rinexo2::_read_obstypes( const string& sat, const string& sys )
{
  gtrace("t_rinexo2::_read_obstypes");

  int ii = 0;
  int addsize = 0;
  unsigned int idx = 0;

  // filter GNSS and SAT
  bool filter_sat = _filter_gnss(sat);
  if( ! filter_sat ){ _xsys++; if( _log ) _log->comment( 2, "rinexo2", "skip "+sat); }

  // general sys- or sat- specific group of observations
  string key = sat;
  t_rnxhdr::t_obstypes::const_iterator itMAP = _mapobs.find(sat);
  if( itMAP == _mapobs.end() ){ key = _csys; }

  // check if key exists
  if( _mapobs.find(key) == _mapobs.end() ){
    cerr << "rinexo2: Warning: _mapobs syskey not found [" << key << "]\n";
       mesg(GWARNING,"Warning: _mapobs syskey not found [" +  key +  "]");
    return _stop_read();
  }   
   
  t_rnxhdr::t_vobstypes::const_iterator it = _mapobs[key].begin();
   
//  t_spt_gobs obs = make_shared<t_gobsgnss>(_rnxhdr.markname(), sat, _epoch);
  t_spt_gobs obs = make_shared<t_gobsgnss>(_site, sat, _epoch);

  // loop over observation types
  while( it != _mapobs[key].end() ){

    idx = 16*ii;

    // read new line for satellite observations
    if( ii > 4 ){
      if( ( addsize = t_gcoder::_getline( _line, _tmpsize ) ) >= 0 ){
        _tmpsize += addsize;
        ii  = 0; 
        idx = 0;
      }else return _stop_read();
    }

    // last check (after new getline or for missing record)
    if( (_line.length() < (idx + 14) ) || trim(_line.substr(idx,14)).empty() ){
          if( filter_sat ){ _null_log( sat, gobs2str(it->first) ); }
    }else if( filter_sat ){
//        cerr << " JD fix: " << sat[0] << " go: " << gobs2str(it->first) << endl;
        string tmp = gobs2str(it->first);
        if( t_rinexo2::_fix_band( string(1, sat[0]), tmp ) ){
           _read_obs( idx, it, obs );
        }
    }

    ii++; it++;
       
  } // while (over observation types)

  if( filter_sat ) _vobs.push_back(obs);

#ifdef DEBUG
  cout << " RNX: "     << _version 
       << " sys: "     << sys
       << " sat: "     << sat
       << " key: "     << key
       << " epo: "     << _epoch.str_ymdhms()
//       << " code/phase ID"
//       << " - band1:"  << " " << gobs2str( obs.code(1)  )
//                       << " " << gobs2str( obs.phase(1) )
//       << " - band2:"  << " " << gobs2str( obs.code(2)  )
//                       << " " << gobs2str( obs.phase(2) )
       << endl;
#endif     

  return 1;
}


// fill single observation element
// ----------
int t_rinexo2::_read_obs(const unsigned int& idx,
                         const t_rnxhdr::t_vobstypes::const_iterator& it,
                         t_spt_gobs obs )
{
  gtrace("t_rinexo2::_read_obs");

  GOBS type = it->first;
  GSYS gsys = obs->gsys();

//  cout << obs->sat() << " _read_obs = " << setw(2) << idx << "  [" << type << "]" 
//       << " " << _epoch.str_ymdhms() << " " << _epoch.dsec() << endl;

  // filter OBS
  if( _obs[gsys].size() > 0 &&
      _obs[gsys].find( trim(gobs2str(type)) ) == _obs[gsys].end() ) return 0;

//string valstr = trim(_line.substr(idx,14)); if( valstr.empty() ) return 0; // eliminate empty string
//string valstr = _line.substr(idx,14).c_str();                              // trim not necessary, empty => 0.0
  string valstr = _line.substr(idx,14);                                      // trim not necessary, empty => 0.0
  double valdbl = str2dbl(valstr);    if( double_eq(valdbl, 0.0) ) return 0; // eliminate 0.000

  obs->addobs( type, valdbl * it->second );
  
#ifdef DEBUG
  cout << obs->sat()
       << " idx = " << setw(2) << idx << fixed << setprecision(3)
       << "  ["     << type
       << "] ["     << gobs2str(type)
//     << "] ["     << setw(14) << str2dbl(trim(_line.substr(idx,14))) * it->second
       << "] ["     << setw(14) << valdbl * it->second
       << "] "      << _epoch.str_ymdhms()
       << " "       << _epoch.dsec()
       << " dsize:" << obs->size()
       << endl; cout.flush();
#endif

  // only for C,L,D, but not SNR observations !!
  if( ( GOBS(type) >=  300 && GOBS(type)  < 400 ) ||  // signal-to-noise ratio
      ( GOBS(type) >= 1300 && GOBS(type) < 1400 ) ){  // signal-to-noise ratio

    // read LLI    
    if( _line.length()-1 >= (idx + 14+1) )
         obs->addlli( type, str2int(_line.substr(idx+14+2,1)));
    else obs->addlli( type, 0);

    // only for phase observations !
    if( GOBS(type) >= 100 && GOBS(type) < 200 ){
          
      // read SNR type instead of original OBS type (always 3-char in G-Nut = RNX2.x 2-char + space)
      GOBS snrtype = pha2snr(type);
      // ---------------------------------
      // slower alternative!
      // GOBS snrtype = str2gobs( "S" + gobs2str(type).substr(1,2));
      // ---------------------------------

      // read SNR observation and convert flag to 
      if( _line.length() -1 >= (idx + 14+2) ){
       
        string snr = trim(_line.substr(idx+14+1,1));
        int i = str2int(snr);

        if( snr == ""  )  obs->addobs( snrtype,  0.0);
        else if( i == 1 ) obs->addobs( snrtype,  6.0);
        else if( i == 2 ) obs->addobs( snrtype, 15.0);	 
        else if( i == 3 ) obs->addobs( snrtype, 20.0);
        else if( i == 4 ) obs->addobs( snrtype, 27.0);
        else if( i == 5 ) obs->addobs( snrtype, 33.0);
        else if( i == 6 ) obs->addobs( snrtype, 39.0);
        else if( i == 7 ) obs->addobs( snrtype, 45.0);
        else if( i == 8 ) obs->addobs( snrtype, 50.0);
        else if( i == 9 ) obs->addobs( snrtype, 60.0);
		 
      }else{ obs->addobs( snrtype,  0.0); }
    }
  }

  // verb() condition to reduce non-verbose mode time
  if( _log && _log->verb() >= 4 ){
    ostringstream lg;
    lg << obs->sat() + ": obs added "
       << _epoch.str_ymdhms()
       << fixed    << setprecision(3)
       << setw(4)  << gobs2str(type)
       << setw(14) << obs->getobs( type );	  
   _log->comment(4,"rinexo2",lg.str());   
 }
 return 1;
}
 

// fill header information
// ----------
int t_rinexo2::_fill_head()
{
  gtrace("t_rinexo2::_fill_head");
   
  int cnt = 0;

  map<string,t_gdata*>::iterator itDAT = _data.begin();
  while( itDAT != _data.end() ){

    if( itDAT->second->id_type() == t_gdata::ALLOBJ ){
       t_gallobj* all_obj = (t_gallobj*)itDAT->second;

       t_gtime epo  = _rnxhdr.first();

       shared_ptr<t_grec>  old_obj, new_obj;
       old_obj = dynamic_pointer_cast<t_grec>(all_obj->obj(_site));
//     old_obj = dynamic_pointer_cast<t_grec>(all_obj->obj(_rnxhdr.markname()));
       new_obj = make_shared<t_grec>();

#ifdef DEBUG
       if( _log ){
	 ostringstream ls;
	 ls << "OBJ set from RINEX: "
	    << _fname << " "
	    << fixed  << setprecision(4)
	    << setw(9)  << _site
	    << setw(9)  << _rnxhdr.markname()
	    << setw(15) << _rnxhdr.aprxyz()[0]
	    << setw(15) << _rnxhdr.aprxyz()[1]
	    << setw(15) << _rnxhdr.aprxyz()[2]
	    << setw(9)  << _rnxhdr.antneu()[0]
	    << setw(9)  << _rnxhdr.antneu()[1]
	    << setw(9)  << _rnxhdr.antneu()[2];
	  _log->comment(1,"rinexo", ls.str());
     //   cerr << ls.str() << endl;
       }
#endif
       new_obj->id(_site);                     // uniq ident
       new_obj->name(_site);                   // prefer before MARKER_NAME which might be incorrect
       new_obj->glog(_log);
       new_obj->addhdr(_rnxhdr, epo, _fname);  // add full rinex header

       // complete from RNX
       if( old_obj == 0 ){
         old_obj = make_shared<t_grec>();
         old_obj->id(_site);                   
         old_obj->name(_site);                 // prefer before MARKER_NAME which might be incorrect
         old_obj->glog(_log);
         old_obj->addhdr(_rnxhdr, epo, _fname);// suppose not necessary

         all_obj->add(old_obj);                // obj has been set -> add obj to gallobj	  
         if( _log ) _log->comment(1,"rinex", "Object created, using RINEX header: "+new_obj->id()+epo.str_ymdhms(" "));
       }
       old_obj->compare(new_obj, epo);
    }
     
    cnt++;
    itDAT++;  
  }

//cout << "End fill_head\n";

  return cnt;
}


// fill observation data structure
// ----------
int t_rinexo2::_fill_data()
{
  gtrace("t_rinexo2::_fill_data");

  int cnt = 0;
  map<string,t_gdata*>::iterator itDAT = _data.begin();
  // adding channel canal to t_gobsgnss from gallobj
  // !!!! must be before adding obs data to t_gallobs
  while( itDAT != _data.end() ){
    if( itDAT->second->id_type() != t_gdata::ALLOBJ ){
      itDAT++;
      continue;
    }
    vector<t_spt_gobs>::iterator itOBS = _vobs.begin();
    while( itOBS != _vobs.end() ){
      t_gallobj* all_obj = (t_gallobj*)itDAT->second;
      shared_ptr<t_gtrn> one_obj = dynamic_pointer_cast<t_gtrn>(all_obj->obj( (*itOBS)->sat() ));
      if( one_obj != 0 && one_obj->id_type() == t_gdata::TRN ){
        int ch = one_obj->channel();
        (*itOBS)->channel(ch);
      }
      itOBS++;
    }          
    itDAT++; 
  }
   
   
  // fill t_gobsgnss 
  itDAT = _data.begin();
  while( itDAT != _data.end() ){
     
    if( itDAT->second->id_type() != t_gdata::ALLOBS ){
      itDAT++;
      continue;
    }

    vector<t_spt_gobs>::const_iterator itOBS = _vobs.begin();	    
    while( itOBS != _vobs.end() ){
       
      if( _log && _log->verb() >= 4 )
          _log->comment(1,"rinexo2",(*itOBS)->sat() + " obs filled "
                                  + (*itOBS)->epoch().str_ymdhms()
                            + " " + (*itOBS)->site());

      ((t_gallobs*)itDAT->second)->addobs( *itOBS );
              
      itOBS++;
    }
    cnt++;
    itDAT++;  
  }

  return cnt;
}


// correct reading stop
// ----------
int t_rinexo2::_check_head()
{
  gtrace("t_rinexo2::_check_head");

  t_gtriple zero(0.0,0.0,0.0);

  if( _rnxhdr.rnxsys() == ' '    ){ mesg(GERROR,"Error: RINEX SYS not available!");   _irc++; } // mandatory
  if( _rnxhdr.rnxver().empty()   ){ mesg(GERROR,"Error: RINEX VERS not available!");  _irc++; } // mandatory
  if( _rnxhdr.program().empty()  &&
      _rnxhdr.program().empty()  ){ mesg(GERROR,"Error: PGM/RUNBY not available!");       _irc++; } // mandatory
  if( _rnxhdr.observer().empty() &&
      _rnxhdr.agency().empty()   ){ mesg(GERROR,"Error: OBSERVER/AGENCY not available!"); _irc++; } // mandatory
  if( _rnxhdr.markname().empty() ){ mesg(GERROR,"Error: MARKER NAME not available!"); _irc++; } // mandatory
  if( _rnxhdr.recnumb().empty()  ){ mesg(GERROR,"Error: REC NUMB not available!");    _irc++; } // mandatory
  if( _rnxhdr.rectype().empty()  ){ mesg(GERROR,"Error: REC TYPE not available!");    _irc++; } // mandatory
  if( _rnxhdr.recvers().empty()  ){ mesg(GERROR,"Error: REC VERS not available!");    _irc++; } // mandatory
  if( _rnxhdr.antnumb().empty()  ){ mesg(GERROR,"Error: ANT NUMB not available!");    _irc++; } // mandatory
  if( _rnxhdr.anttype().empty()  ){ mesg(GERROR,"Error: ANT TYPE not available!");    _irc++; } // mandatory
  if( _rnxhdr.mapobs().size()==0 ){ mesg(GERROR,"Error: OBS TYPE not available!");    _irc++; } // mandatory
  if( _rnxhdr.aprxyz() == zero   ){ mesg(GERROR,"Error: APR XYZ not available!");     _irc++; } // mandatory
  if( _rnxhdr.first() == FIRST_TIME
                                 ){ mesg(GERROR,"Error: FIRST TIME not available!");  _irc++; } // mandatory
  
//if( _rnxhdr.last()  == LAST_TIME
//                               ){ mesg(GWARNING,"Warning: LAST TIME not available!");       } // optional --> skip
//if( _rnxhdr.leapsec()  <= 0    ){ mesg(GWARNING,"Warning: LEAPSEC not available!");         } // optional --> skip
//if( _rnxhdr.numsats()  <= 0    ){ mesg(GWARNING,"Warning: NUMB SAT not available!");        } // optional --> skip
  if( _rnxhdr.interval() <= 0    ){ mesg(GWARNING,"Warning: INTERVAL not available!");        } // optional --> warning
//if( _rnxhdr.program().empty()  ){ mesg(GWARNING,"Warning: PGM empty field!");               } // optional --> skip
//if( _rnxhdr.runby().empty()    ){ mesg(GWARNING,"Warning: RUN BY empty field!");            } // optional --> skip
  if( _rnxhdr.marknumb().empty() ){ mesg(GWARNING,"Warning: MARKER NUMB empty field!");       } // optional --> warning
//if( _rnxhdr.observer().empty() ){ mesg(GWARNING,"Warning: OBSERVER empty field!");          } // optional --> skip
//if( _rnxhdr.marktype().empty() ){ mesg(GWARNING,"Warning: MARKER TYPE empty field!");       } // optional --> skip
//if( _rnxhdr.strength().empty() ){ mesg(GWARNING,"Warning: STRENGTH empty field!");          } // optional --> skip

  return _irc;
}
   
  
// correct reading stop
// ----------
int t_rinexo2::_null_log( const string& sat, const string& obstype )
{
  gtrace("t_rinexo2::_null_log");

  if( _log && _log->verb() >= 4 )
      _log->comment(4,"rinexo2",_epoch.str_ymdhms(sat + ": obs  null ") + " " + obstype );
  return 0;
}


// correct reading stop
// ----------
int t_rinexo2::_stop_read()
{   
  _tmpsize = 0;
  _complete = false;
  return -1;

}

   
// fix band
// ----------
int t_rinexo2::_fix_band(string sys, string& go)
{
  gtrace("t_rinexo2::_check_band");

  GSYS gsys = (sys.length() < 3 ) ? t_gsys::char2gsys( sys[0] ) : t_gsys::str2gsys(sys);

#ifdef DEBUG
  cerr << " sys : " << sys << " " << t_gsys::gsys2str(gsys)
       << " gobs: "  << go << "  GO : " << go.substr(1,1)
       << " "  << str2gobsband(go) 
       << " len:" << sys.length()
       << endl;
#endif

  // EXCEPTIONs FOR EXTENDED RINEX 2.11 (L7 etc)
  switch ( gsys ){
    case GLO : if( _version <= "2.12" && go[1] == '7' ){
                 if( _log ) _log->comment(0,"rinexo2","Warning: extended RINEX 2.11 format (GLO C7/L7 skipped)");
                 else               cerr << "rinexo2 - Warning: extended RINEX 2.11 format (GLO C7/L7 skipped)\n";
                                        mesg(GWARNING,"Warning: extended RINEX 2.11 format (GLO C7/L7 skipped)");
                 return 0;
              }; break;

    case GPS : if( _version <= "2.12" && go[1] == '7' ){
                 if( _log ) _log->comment(0,"rinexo2","Warning: extended RINEX 2.11 format (GPS C7/L7 skipped)");
                 else               cerr << "rinexo2 - Warning: extended RINEX 2.11 format (GPS C7/L7 skipped)\n";
                                        mesg(GWARNING,"Warning: extended RINEX 2.11 format (GPS C7/L7 skipped)");
                 return 0;
              }; break;
    case BDS : if( _version <= "2.12" ){
                 if( _log ) _log->comment(0,"rinexo2","Warning: extended RINEX 2.11 format (BDS skipped)");
                 else               cerr << "rinexo2 - Warning: extended RINEX 2.11 format (BDS skipped)\n";
                                        mesg(GWARNING,"Warning: extended RINEX 2.11 format (BDS skipped)");
                 return 0;
              }; break;
    case QZS : if( _version <= "2.12" ){
                 if( _log ) _log->comment(0,"rinexo2","Warning: extended RINEX 2.11 format (QZS skipped)");
                 else               cerr << "rinexo2 - Warning: extended RINEX 2.11 format (QZS skipped)\n";
                                        mesg(GWARNING,"Warning: extended RINEX 2.11 format (QZS skipped)");
                 return 0;
              }; break;
    case IRN : if( _version <= "2.12" ){
                 if( _log ) _log->comment(0,"rinexo2","Warning: extended RINEX 2.11 format (IRN skipped)");
                 else               cerr << "rinexo2 - Warning: extended RINEX 2.11 format (IRN skipped)\n";
                                        mesg(GWARNING,"Warning: extended RINEX 2.11 format (IRN skipped)");
                 return 0;
              }; break;
    
    default : return 1;
  }

  return 1;
}

} // namespace
