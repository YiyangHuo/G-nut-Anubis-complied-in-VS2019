
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

#include "gall/gallprec.h"
#include "gcoders/sp3.h"
#include "gutils/gtriple.h"
#include "gutils/gtypeconv.h"
 
using namespace std;

namespace gnut {

/* ----------
 * constructor
 */
t_sp3::t_sp3( t_gsetbase* s, string version, int sz )
  : t_gcoder( s, version, sz )
{
  gtrace("t_sp3::constructor");

  _start.tsys(t_gtime::GPS);
  _lastepo.tsys(t_gtime::GPS);
  _nrecord = -1;
  _nrecmax = -1;
}


/* ----------
 * SP3 header
 */
int t_sp3::decode_head( char* buff, int sz, vector<string>& errmsg)
{ 
  gtrace("t_sp3::decode_head");

  _mutex.lock();

  if( t_gcoder::_add2buffer( buff, sz) == 0 ){ _mutex.unlock(); return 0; };

  string tmp;
  int consume = 0;
  int tmpsize = 0;
  while( ( tmpsize = t_gcoder::_getline( tmp )) >= 0 ){

    consume += tmpsize;
     
    // first information
    if( tmp.substr(0,1) == "#" )
    {
       // first line
       if(  tmp.substr(1,1) != "a" ||  // 60-columns
            tmp.substr(1,1) != "c"     // 80-columns
         ){
	  
	 _version = tmp.substr(1,1);
         _start.from_str( "%Y %m %d  %H %M %S", tmp.substr(3,28) );
         _nepochs = str2int(tmp.substr(32,7));
         _orbrefs  = tmp.substr(46,5);
	 _orbtype  = tmp.substr(52,3);
     
       // second line
       }else if( tmp.substr(1,1) == "#" ){
         _orbintv  = (long)str2dbl(tmp.substr(24,14)); // [sec]
         if( _log ){
	   ostringstream ltmp;
           ltmp << "start time = " << _start.str(" %Y-%m-%d %H:%M:%S")
                <<     "  refs = " << _orbrefs
                <<     "  type = " << _orbtype
                <<     "  intv = " << _orbintv
                <<     "  nepo = " << _nepochs;
	       _log->comment(2,"sp3",ltmp.str());
	 }
       }else{
         if( _log ) _log->comment(1,"sp3"," unknown record: " + tmp);
       }
	  
    }else if( tmp.substr(0,2) == "+ " ){
      ostringstream ltmp("reading satellites:");
      for( int i = 9; i < 60 ; i=i+3 ){
         _prn.push_back( tmp.substr(i,3) );
 	 ltmp << " " << tmp.substr(i,3);
      }
      if( _log ) _log->comment(2,"sp3",ltmp.str());
       
    }else if( tmp.substr(0,2) == "++" ){
      ostringstream ltmp("reading accuracies:");
      for( int i = 9; i < 60 ; i=i+3 ){	    
         _acc.push_back( str2int(tmp.substr(i,3)) );
 	 ltmp << " " << tmp.substr(i,3);
      } 
      if( _log ) _log->comment(2,"sp3",ltmp.str());
       
    }else if( tmp.substr(0,2) == "%c" ){
      _timesys.push_back( tmp.substr( 3,2) );
      _timesys.push_back( tmp.substr( 6,2) );
      _timesys.push_back( tmp.substr( 9,3) );
      _timesys.push_back( tmp.substr(13,3) );
      _timesys.push_back( tmp.substr(17,4) );
      _timesys.push_back( tmp.substr(22,4) );
      _timesys.push_back( tmp.substr(27,4) );
      _timesys.push_back( tmp.substr(32,4) );
      _timesys.push_back( tmp.substr(37,5) );
      _timesys.push_back( tmp.substr(43,5) );
      _timesys.push_back( tmp.substr(49,5) );
      _timesys.push_back( tmp.substr(55,5) );
      if( _log ) _log->comment(2,"sp3","reading satellite systems");

    }else if( tmp.substr(0,2) == "%f" ){
      _accbase.push_back( str2int(tmp.substr( 3,9)) );
      _accbase.push_back( str2int(tmp.substr(14,9)) );
      if( _log ) _log->comment(2,"sp3","reading PV base");
       
    }else if( tmp.substr(0,2) == "%i" ){
      if( _log ) _log->comment(2,"sp3","additional info");
       
    }else if( tmp.substr(0,2) == "/*" ){
      if( _log ) _log->comment(2,"sp3","comments");

    }else if( tmp.substr(0,2) == "* " ){  // END OF HEADER !!!
      if( _log ) _log->comment(2,"sp3","END OF HEADER");
      _mutex.unlock(); return -1;

    }else{
      if( _log ) _log->comment(2,"sp3","END OF HEADER");
      else       cerr << "warning: unknown header message :" << tmp << endl;
    }
       
    t_gcoder::_consume(tmpsize);
  }

  _mutex.unlock(); return consume;
}

   
/* ----------
 * SP3 body
 */
int t_sp3::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
  gtrace("t_sp3::decode_data");
   
  _mutex.lock();

  if( t_gcoder::_add2buffer(buff, sz) == 0 ){ _mutex.unlock(); return 0; };
   
  string tmp;
//  int consume = 0;
//  int recsize = 0;
  int tmpsize = 0;
//  bool epoch_defined = false;

//  while( ( tmpsize = t_gcoder::_getline( tmp, recsize ) ) >= 0 ){
  while( ( tmpsize = t_gcoder::_getline( tmp, 0 ) ) >= 0 ){

#ifdef DEBUG
    cout << " 0: " << tmp ;
#endif

    // EPOCH record
    if( tmp.substr(0,1) == "*"   ||
	tmp.substr(0,3) == "EOF"
    ){

      if( tmp.substr(0,3) == "EOF" ){
	if( _log ) _log->comment(3,"sp3","EOF found");
        t_gcoder::_consume(tmpsize);
	_mutex.unlock(); return tmpsize;
      }

      if( _nrecord > 0 && _nrecord != _nrecmax ){
	ostringstream ltmp;
        ltmp << "warning: not equal number of satellites " << _nrecord <<  " " << _nrecmax << "!";
	if( _log ) _log->comment(0,"sp3",ltmp.str());
        else       cerr << ltmp.str();	 
      }

      char dummy;
      int yr,mn,dd,hr,mi;
      double sc;
      stringstream ss( tmp );
      ss >> dummy >> yr >> mn >> dd >> hr >> mi >> sc;
  
      if( ss.fail() ){
        printf("Warning: incorrect SP3 epoch record:%s \n\n", ss.str().c_str() );
        t_gcoder::_consume(tmpsize);
        _mutex.unlock(); return -1;
      }
       
      int sod = hr*3600 + mi*60 + (int)sc;
      _lastepo.from_ymd(yr,mn,dd,sod);
      ostringstream ltmp;
      ltmp << "reading EPOCH [" << _nrecord << "] - " << _lastepo.str(" %Y-%m-%d %H:%M:%S");
      if( _log ) _log->comment(3,"sp3",ltmp.str());
      _nrecord = -1;
    }
    
    // POSITION reccord
    if( tmp.substr(0,1) == "P" ){ // and epoch_defined ){

      t_gtriple  xyz(0.0, 0.0, 0.0);
      t_gtriple dxyz(0.0, 0.0, 0.0);
      double t = 0.0, dt = 0.9;
	 
      char sat[3+1]; sat[3] = '\0';
      char flg;
      double pos[4] = {0.0, 0.0, 0.0, 0.0};
//      double var[4] = {0.0, 0.0, 0.0, 0.0};
      stringstream ss( tmp );

      ss >> noskipws >> flg    >> sat[0] >> sat[1] >> sat[2] 
	 >>   skipws >> pos[0] >> pos[1] >> pos[2] >> pos[3];
       
      string prn = t_gsys::eval_sat( string(sat) );

      for(int i=0; i<3; i++)
        if( pos[i] == 0.0 )
	      xyz[i] = UNDEFVAL_POS; // gephprec
        else  xyz[i] = pos[i]*1000;  // km -> m
       
      if( pos[3] > 999999 )
	    t = UNDEFVAL_CLK;        // gephprec
      else  t = pos[3]/1000000;      // us -> s
	 
#ifdef DEBUG
     cout << " reading:" << _lastepo.str(" %Y-%m-%d %H:%M:%S ") << flg << " " << prn << " " << sat
          << " " << pos[0] << " " << pos[1] << " " << pos[2] << " " << pos[3]*1000000 << endl;
#endif

      // fill single data record
      if( !_filter_gnss(prn) ){
        if( _log ) _log->comment( 4, "sp3", "skip "+prn); 
      }else{	    
	map<string,t_gdata*>::iterator it = _data.begin();
        while( it != _data.end() ){
 	  if( it->second->id_type() == t_gdata::ALLPREC )
	     ((t_gallprec*)it->second)->addpos( prn, _lastepo, xyz, t, dxyz, dt );
	  it++;
        }
      }

      if( ss.fail() ){
        if( _log ) _log->comment(1,"sp3", "warning: incorrect SP3 data record: " + ss.str() );
        else       cerr << "warning: incorrect SP3 data record: " << ss.str() << endl;
        t_gcoder::_consume(tmpsize);
        _mutex.unlock(); return -1;
      }
      cnt++;
    }
     
    // VELOCITY reccord
    if( tmp.substr(0,1) == "V" ){ // and epoch_defined ){

      char sat[3+1]; sat[3] = '\0';
      char flg;
      double vel[4] = {0.0, 0.0, 0.0, 0.0};
      double var[4] = {0.0, 0.0, 0.0, 0.0};
       
      stringstream ss( tmp );
      ss >> noskipws >> flg    >> sat[0] >> sat[1] >> sat[2] 
	 >>   skipws >> vel[0] >> vel[1] >> vel[2] >> vel[3];
       
      string prn = t_gsys::eval_sat( string(sat) );
      
      // fill single data record
      map<string,t_gdata*>::iterator it = _data.begin();
      while( it != _data.end() ){

        // fill single data record
        if( !_filter_gnss(prn) ){
          if( _log ) _log->comment( 4, "sp3", "skip "+prn); 
        }else{
  	  map<string,t_gdata*>::iterator it = _data.begin();
          while( it != _data.end() ){
 	    if( it->second->id_type() == t_gdata::ALLPREC )
               ((t_gallprec*)it->second)->addvel( prn, _lastepo, vel, var );
	    it++;
          }
        }
      }

      if( ss.fail() ){
        if( _log ) _log->comment(1,"sp3", "warning: incorrect SP3 data record: " + ss.str() );
        else       cerr << "warning: incorrect SP3 data record: " << ss.str() << endl;
        t_gcoder::_consume(tmpsize);
        _mutex.unlock(); return -1;
      }
    }
     
    _nrecord++;
     
    if( _nrecord > _nrecmax ) _nrecmax = _nrecord;

    t_gcoder::_consume(tmpsize);
  }
  _mutex.unlock(); return tmpsize;
}

} // namespace
