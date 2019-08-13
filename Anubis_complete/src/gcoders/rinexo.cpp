
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
#include <algorithm>

#include "gcoders/rinexo.h"
#include "gdata/gobsgnss.h"
 
using namespace std;

namespace gnut {

// constructor
// ----------
t_rinexo::t_rinexo( t_gsetbase* s, string version, int sz )
  : t_rinexo3( s, version, sz )
{}

   
// OBS-RINEX header
//   - read individual lines (tmpsize once used and consumed)
//   - read block of lines at once (tmpsize cummulated and then consumed)
// ----------
int t_rinexo::decode_head( char* buff, int sz, vector<string>& errmsg)
{ 
  gtrace("t_rinexo::decode_head");

  _mutex.lock();

  int add = t_gcoder::_add2buffer(buff, sz);     
  if( add <  0 ){ _mutex.unlock(); return -1; };
  if( add == 0 ){ _mutex.unlock(); return  0; };

  _complete = true;
  _consume = 0;
  _tmpsize = 0;
  _line = "";

  this->_decode_head();

  _mutex.unlock(); return _consume;
}

   
// OBS-RINEX body
//   - read block of lines for each epoch (tmpsize cummulated and then consumed)
// ----------
int t_rinexo::decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{  
  gtrace("t_rinexo::decode_data");

  if( _hdr ) return -1; // thin execution, header only
  
//  // test empty marker name from header
//  if( _rnxhdr.markname().empty() ){
//     if( _log ) _log->comment(0,"rinexo", "Warning: empty marker name not allowed in RINEX header [skipped]." );
//                            mesg(GWARNING,"Warning: empty marker name not allowed in RINEX header!");
//     ++_irc;
//     return -1;
//  }   

  // filter-out unknown marker/receiver
  if( _rec.size() == 0 ){ _rec.insert(_site); }
  if( _rec.find( _site             ) == _rec.end() &&
      _rec.find( _site.substr(0,4) ) == _rec.end() ){
     if( _log ) _log->comment(2,"rinexo", _site + " marker name not found in configuration [skipped]." );
                           mesg(GMESSAGE, _site + " marker name not found in configuration!");
     ++_irc;
     return -1;
  }

  _flag     = '0';                   // special event flag
  _line     = "";                    // working line
  _nsat     = 0;                     // special event flag
  _tmpsize  = 0;                     // currently processing byte counter
  _consume  = 0;                     // total byte counter
  _complete = true;                  // flag for completed block
   
  _mutex.lock();

  // complete main buffer
  int add = t_gcoder::_add2buffer(buff, sz);
  if( add <  0 ){ _mutex.unlock(); return -1; };
  if( add == 0 ){ _mutex.unlock(); return  0; };
 
  this->_decode_data();   

  cnt = _count;
  _mutex.unlock(); return _consume;
}


// RINEX encoder
// ---------
int t_rinexo::encode_head(char* buff, int sz,           vector<string>& errmsg)
{
  gtrace("t_rinexo::encode_head");
  _mutex.lock();
  
  if( _ss_position == 0 ) {
    map<string,t_gdata*>::iterator it = _data.begin();
    for( it = _data.begin(); it != _data.end(); ++it ) {
      if( it->second->id_type()  == t_gdata::ALLOBJ ) {
          t_gallobj* all_obj = (t_gallobj*) it->second;
          map<string, shared_ptr<t_gobj>> allGrec = all_obj->objects(t_gdata::REC);
          //map<string, shared_ptr<t_gobj>>::const_iterator 
          if(allGrec.begin() != allGrec.end()) {
              //shared_ptr<t_gobj> oneGrec = allGrec.begin()->second;
              shared_ptr<t_grec> grecCasted = dynamic_pointer_cast<t_grec> (allGrec.begin()->second);
              t_grec::t_maphdr hdrFromGrec = grecCasted->gethdr();
              if(hdrFromGrec.begin() != hdrFromGrec.end()) {
                  t_rnxhdr rnxHdr = hdrFromGrec.begin()->second;
                  
                  // write RINEX VERSION / TYPE
                  const string rnxver = rnxHdr.rnxver();
                  if( rnxver != string("3.03") ){
                    if( _log ) _log->comment(0, "rinexo","Warning: not RINEX observation file version 3.03");
                    else                cerr << "rinexo - Warning: not RINEX observation file version 3.03\n";
//                                               _seterr("Warning: not RINEX observation file version!");
                                        
                    return -1;
                  }
                  _ss << string(9 - rnxver.size(), ' ') << rnxver
                          << "           " << "OBSERVATION DATA"
                          << "    " << rnxHdr.rnxsys() << "                   "
                          << "RINEX VERSION / TYPE" << endl;
                  
                  // write PGM / RUN BY / DATE
                  const string program = rnxHdr.program();
                  const string runby = rnxHdr.runby();
                  const t_gtime gtime = rnxHdr.gtime();
                  string theDate = gtime.str("%Y%m%d %H%M%S %T");
                  string::size_type locationOfLOC = theDate.find("LOC");
                  if(locationOfLOC != string::npos) theDate.replace(locationOfLOC, 3, "LCL");
                  
                  _ss << program << string(20 - program.size(), ' ')
                      << runby   << string(20 - runby.size(),   ' ')
                      << theDate << ' ' << "PGM / RUN BY / DATE" << endl;
                  
                  // write COMMENT (if present, all lines)
                  vector<string> comment = rnxHdr.comment();
                  for(vector<string>::iterator commIt = comment.begin(); commIt != comment.end(); ++commIt) {
                      _ss << *commIt << "COMMENT" << endl; 
                  }
                  
                  // write MARKER NAME 
                  const string markname = rnxHdr.markname();
                  _ss << markname << string(60 - markname.size(), ' ') << "MARKER NAME" << endl;
                  
                  // write MARKER NUMBER (if present)
                  const string marknumb = rnxHdr.marknumb();
                  if(!marknumb.empty())
                    _ss << marknumb << string(60 - marknumb.size(), ' ') << "MARKER NUMBER" << endl;
                  
                  // write MARKER TYPE
                  const string marktype = rnxHdr.marktype();
                  _ss << marktype << string(60 - marktype.size(), ' ') << "MARKER TYPE" << endl;
                  
                  // write OBSERVER / AGENCY
                  const string observer = rnxHdr.observer();
                  const string agency = rnxHdr.agency();
                  _ss << observer << string(20 - observer.size(), ' ')
                      << agency   << string(40 - agency.size(),   ' ')
                      << "OBSERVER / AGENCY" << endl;
                  
                  // write REC # / TYPE / VERS
                  const string recnumb = rnxHdr.recnumb();
                  const string rectype = rnxHdr.rectype();
                  const string recvers = rnxHdr.recvers();
                  _ss << recnumb << string(20 - recnumb.size(), ' ')
                      << rectype << string(20 - rectype.size(), ' ')
                      << recvers << string(20 - recvers.size(), ' ')  
                      << "REC # / TYPE / VERS" << endl;
                  
                  // write ANT # / TYPE
                  const string antnumb = rnxHdr.antnumb();
                  const string anttype = rnxHdr.anttype();
                  _ss << antnumb << string(20 - antnumb.size(), ' ')
                      << anttype << string(40 - anttype.size(), ' ')    
                      << "ANT # / TYPE" << endl;
                  
                  // write APPROX POSITION XYZ
                  const t_gtriple aprxyz = rnxHdr.aprxyz();
                  const double aprPosX = aprxyz[0];
                  const double aprPosY = aprxyz[1];
                  const double aprPosZ = aprxyz[2];
                  const string aprPosXstr = trim(dbl2str(aprPosX, 4));
                  const string aprPosYstr = trim(dbl2str(aprPosY, 4));
                  const string aprPosZstr = trim(dbl2str(aprPosZ, 4));
                  _ss << string(14 - aprPosXstr.size(), ' ') << aprPosXstr 
                      << string(14 - aprPosYstr.size(), ' ') << aprPosYstr
                      << string(14 - aprPosZstr.size(), ' ') << aprPosZstr
                      << "                  "
                      << "APPROX POSITION XYZ" << endl;
                  
                  // write ANTENNA: DELTA H/E/N
                  const t_gtriple antneu = rnxHdr.antneu();
                  const double antneuH = antneu[2]; // swapped with 0!! Up
                  const double antneuE = antneu[1]; // East
                  const double antneuN = antneu[0]; // swapped with 2!! North
                  const string antneuHstr = trim(dbl2str(antneuH, 4));
                  const string antneuEstr = trim(dbl2str(antneuE, 4));
                  const string antneuNstr = trim(dbl2str(antneuN, 4));
                  _ss << string(14 - antneuHstr.size(), ' ') << antneuHstr
                      << string(14 - antneuEstr.size(), ' ') << antneuEstr
                      << string(14 - antneuNstr.size(), ' ') << antneuNstr
                      << "                  "
                      << "ANTENNA: DELTA H/E/N" << endl;
                  
                  // write ANTENNA: DELTA X/Y/Z (if present)
                  const t_gtriple antxyz = rnxHdr.antxyz();
                  if (!double_eq(antxyz[0], 0) || !double_eq(antxyz[1], 0)
                          || !double_eq(antxyz[2], 0)) {
                    const double antDelX = antxyz[0];
                    const double antDelY = antxyz[1];
                    const double antDelZ = antxyz[2];
                    const string antDelXstr = trim(dbl2str(antDelX, 4));
                    const string antDelYstr = trim(dbl2str(antDelY, 4));
                    const string antDelZstr = trim(dbl2str(antDelZ, 4));
                    _ss << string(14 - antDelXstr.size(), ' ') << antDelXstr 
                        << string(14 - antDelYstr.size(), ' ') << antDelYstr
                        << string(14 - antDelZstr.size(), ' ') << antDelZstr
                        << "                  "
                        << "ANTENNA: DELTA X/Y/Z" << endl;
                  }
                  
                  // write ANTENNA:PHASECENTER (if present)
                  /* warning, all data get deleted in t_rinexo2::clear() which
                   * is called in t_gio::run_write() in the beginning, thus
                   * deleting any data you added through t_rinexo2::_decode_head().
                   Also after filling _pcosat,  _pcosys and _pcoecc in
                   * t_rinexo2::_decode_head(), it's not saved anywhere afterwards.
                   Also in t_rinexo2::_decode_head() it checks for "ANTENNA: PHASECENTER"
                   and not "ANTENNA:PHASECENTER" which is correct according to
                   * RINEX 3.03 standard */
                  for (unsigned int i = 0; i < _pcosat.size(); i++) {
                        const t_gtriple ecc = _pcoecc[i];
                        const string eccNX = trim(dbl2str(ecc[0], 4));
                        const string eccEY = trim(dbl2str(ecc[1], 4));
                        const string eccUZ = trim(dbl2str(ecc[2], 4));
                      
                        _ss << _pcosat[i] << ' '
                            << _pcosys[i] << string(3 - _pcosys[i].size(), ' ')
                            << string(9 - eccNX.size(), ' ') << eccNX 
                            << string(14 - eccEY.size(), ' ') << eccEY
                            << string(14 - eccUZ.size(), ' ') << eccUZ
                            << "                  "
                            << "ANTENNA:PHASECENTER" << endl;  
                  }
                  
                  // write ANTENNA: B.SIGHT XYZ (if present)
                  // not implemented
                  
                  // write ANTENNA: ZERODIR AZI (if present)
                  // not implemented
                  
                  // write ANTENNA: ZERODIR XYZ (if present)
                  // not implemented
                  
                  // write CENTER OF MASS: XYZ (if present)
                  // not implemented
                  
                  // write SYS / # / OBS TYPES
                  t_rnxhdr::t_obstypes mapobs = rnxHdr.mapobs();
                  t_rnxhdr::t_obstypes::iterator mapobsIt;
                  // for each satellite system
                  for(mapobsIt = mapobs.begin(); mapobsIt != mapobs.end(); ++mapobsIt) {
                      const string satSys = mapobsIt->first;
                      t_rnxhdr::t_vobstypes obsDescriptorFactorPairs = mapobsIt->second;
                      const int numObsTypes = obsDescriptorFactorPairs.size();
                      
                      _ss << satSys << "  " << setw(3) << numObsTypes;
                      
                      // satSys is 1 char, 2 spaces, numObsType with width 3 
                      int howManyOf60Written = 6;
                      int howManyWrittenOnLine = 0;
                      int howManyWrittenTotal = 0;
                      // for each obs descriptor associated with this sat sys
                      t_rnxhdr::t_vobstypes::iterator odfpIt;
                      for(odfpIt = obsDescriptorFactorPairs.begin(); odfpIt != obsDescriptorFactorPairs.end();
                              ++odfpIt) {
                          //GOBS gobs1 = (*odfpIt)->first;
                          GOBS gobs = (*odfpIt).first;
                          _ss << ' ' << gobs2str(gobs);
                          
                          howManyOf60Written += 4;
                          howManyWrittenOnLine++;
                          howManyWrittenTotal++;
                          if(  howManyWrittenOnLine == 13 || howManyWrittenTotal  == numObsTypes)
                          {
                            _ss << string(60 - howManyOf60Written, ' ') << "SYS / # / OBS TYPES" << endl;
                            if(howManyWrittenTotal != numObsTypes) _ss << "      ";
                            howManyWrittenOnLine = 0;
                            howManyOf60Written = 6;
                          }
                      }
                  }
                  
                  // write SIGNAL STRENGTH UNIT (if present)
                  const string strength = rnxHdr.strength();
                  if(!strength.empty())
                    _ss << strength << string(60 - strength.size(), ' ') << "SIGNAL STRENGTH UNIT" << endl;
                  
                  // write INTERVAL (if present)
                  const double interval = rnxHdr.interval();
                  if(!double_eq(interval, 0)) {
                      // 3 digits after decimal point
                      const string intervalStr = trim(dbl2str(interval, 3));
                      _ss << string(10 - intervalStr.size(), ' ') << intervalStr
                          << string(50, ' ') << "INTERVAL" << endl;
                  }
                  
                  // write TIME OF FIRST OBS
                  t_gtime first = rnxHdr.first();
                  string timeSys = first.sys();
                  int yearFO, monthFO, doMFO;
                  first.ymd(yearFO, monthFO, doMFO);
                  int hourFO, minuteFO, secondFO;
                  first.hms(hourFO, minuteFO, secondFO);
                  double dsec = first.dsec();
                  string secdsecStr = trim(dbl2str(((double)secondFO)+dsec, 7));
                  _ss << setw(6) << yearFO << setw(6) << monthFO
                          << setw(6) << doMFO
                          << setw(6) << hourFO << setw(6) << minuteFO
                          << string(13 - secdsecStr.size(), ' ') << secdsecStr
                          << "     " << string(3 - timeSys.size(), ' ')
                          << timeSys << "         " << "TIME OF FIRST OBS" << endl;
                  
                  // write TIME OF LAST OBS (if present)
                  t_gtime last = rnxHdr.last();
                  if (last != LAST_TIME) {
                    string timeSysLast = last.sys();
                    int yearLO, monthLO, doMLO;
                    last.ymd(yearLO, monthLO, doMLO);
                    int hourLO, minuteLO, secondLO;
                    last.hms(hourLO, minuteLO, secondLO);
                    double dsecLO = last.dsec();
                    string secdsecStrLO = trim(dbl2str(((double)secondLO)+dsecLO, 7));
                    _ss << setw(6) << yearLO << setw(6) << monthLO
                            << setw(6) << doMLO
                            << setw(6) << hourLO << setw(6) << minuteLO
                            << string(13 - secdsecStrLO.size(), ' ') << secdsecStrLO
                            << "     " << string(3 - timeSysLast.size(), ' ')
                            << timeSysLast << "         " << "TIME OF LAST OBS" << endl;
                  }
                  
                  // write RCV CLOCK OFFS APPL (if present)
                  // not implemented
                  
                  // write SYS / DCBS APPLIED (if present)
                  // not implemented
                  
                  // write SYS / PCVS APPLIED (if present)
                  // not implemented
                  
                  // write SYS / SCALE FACTOR (if present)
                  // for each satellite system
                  for(mapobsIt = mapobs.begin(); mapobsIt != mapobs.end(); ++mapobsIt) {
                      const string satSys = mapobsIt->first;
                      t_rnxhdr::t_vobstypes obsDescriptorFactorPairs = mapobsIt->second;
                      
                      set<int> extractedFactors;
                      t_rnxhdr::t_vobstypes::iterator odfpIt;
                      for(odfpIt = obsDescriptorFactorPairs.begin(); odfpIt != obsDescriptorFactorPairs.end(); ++odfpIt) {
                          extractedFactors.insert((*odfpIt).second);                          
                      }
                      if(extractedFactors.size() == 1
                              && extractedFactors.count(1) == 1) {
                          // contains just factors equal to 1, can be skipped
                          continue;
                      }
                      if (extractedFactors.size() == 1) {
                          _ss << satSys << " " << setw(4)
                              << (*(extractedFactors.begin()))
                              << string(54, ' ')
                              << "SYS / SCALE FACTOR" << endl;
                          continue;
                      }
                      
                      
                      // for each different factor associated with this sat sys
                      for(set<int>::iterator efIt = extractedFactors.begin(); efIt != extractedFactors.end(); ++efIt) {
                        const int currentFactor = (*efIt);
                        // first get num of obs types for this factor and sat sys 
                        int numObsTypes = 0;
                        vector<GOBS> neededObsDescrs;
                        for(odfpIt = obsDescriptorFactorPairs.begin(); odfpIt != obsDescriptorFactorPairs.end();
                                  ++odfpIt) {
                            if ((*odfpIt).second == currentFactor) {
                                numObsTypes++;
                                neededObsDescrs.push_back((*odfpIt).first);
                            }
                        }
                        _ss << satSys << " "  << setw(4) << currentFactor
                                      << "  " << setw(2) << numObsTypes;
                        
                        /* satSys is 1 char, 1 space, factor width 4, 2 spaces
                         * numObsType with width 2 */
                        int howManyOf60Written = 10;
                        int howManyWrittenOnLine = 0;
                        int howManyWrittenTotal = 0;
                        /* for each obs descriptor associated with this sat sys
                         * and this factor */
                        for(vector<GOBS>::iterator nodIT = neededObsDescrs.begin(); nodIT != neededObsDescrs.end(); ++nodIT) 
                        {
                            //GOBS gobs1 = (*odfpIt)->first;
                            
                            GOBS gobs = (*nodIT);
                            _ss << ' ' << gobs2str(gobs);

                            howManyOf60Written += 4;
                            howManyWrittenOnLine++;
                            howManyWrittenTotal++;
                            
                            if (howManyWrittenOnLine == 12 || howManyWrittenTotal == numObsTypes) {
                                _ss << string(60 - howManyOf60Written, ' ') << "SYS / SCALE FACTOR" << endl;
                                if(howManyWrittenTotal != numObsTypes) _ss << "          ";
                                howManyWrittenOnLine = 0;
                                howManyOf60Written = 10;
                            }
                        }
                      }   
                  }
                  
                  // write SYS / PHASE SHIFT
                  // not implemented (even though mandatory, at least in 3.03)
                  
                  // write GLONASS SLOT / FRQ #
                  // not implemented (even though mandatory, at least in 3.03)
                  
                  // write GLONASS COD/PHS/BIS
                  // not implemented (even though mandatory, at least in 3.03)
                  
                  // write LEAP SECONDS (if present)
                  // implemented only partially
                  const int leapSecs = rnxHdr.leapsec();
                  if(leapSecs != 0)
                    _ss << setw(6) << leapSecs << string(54, ' ') << "LEAP SECONDS" << endl;
                  
                  // write # OF SATELLITES (if present)
                  const int numSats = rnxHdr.numsats();
                  if(numSats != 0)
                    _ss << setw(6) << numSats << string(54, ' ') << "# OF SATELLITES" << endl;
                  
                  // write PRN / # OF OBS (if present)
                  // not implemented
                  
                  // write END OF HEADER
                  _ss << string(60, ' ') << "END OF HEADER" << endl;
              
                  break; // break out of for cycle going through _data
              }
          }
      }
    }
  }

  int size = _fill_buffer( buff, sz );
  
  _mutex.unlock();
  return size;  
}


// RINEX encode   
// ---------
int t_rinexo::encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg)
{
  gtrace("t_rinexo::encode_data");
  _mutex.lock();
  
  if( _ss_position == 0 ){
      map<string,t_gdata*>::iterator it = _data.begin();
      //string markname = "";
      t_rnxhdr::t_obstypes mapobs;
      
      for( it = _data.begin(); it != _data.end(); ++it ) {
        if( it->second->id_type()  == t_gdata::ALLOBJ ) {
            t_gallobj* all_obj = (t_gallobj*) it->second;
            map<string, shared_ptr<t_gobj>> allGrec = all_obj->objects(t_gdata::REC);
            //map<string, shared_ptr<t_gobj>>::const_iterator 
            if(allGrec.begin() != allGrec.end()) {
                //shared_ptr<t_gobj> oneGrec = allGrec.begin()->second;
                shared_ptr<t_grec> grecCasted = dynamic_pointer_cast<t_grec> (allGrec.begin()->second);
                t_grec::t_maphdr hdrFromGrec = grecCasted->gethdr();
                if(hdrFromGrec.begin() != hdrFromGrec.end()) {
                    t_rnxhdr rnxHdr = hdrFromGrec.begin()->second;
                    //markname = rnxHdr.markname();
                    mapobs = rnxHdr.mapobs();
                    break;
                }
            }
        }
      }
      
    //if (markname.size() > 0) {
    for( it = _data.begin(); it != _data.end(); ++it ) {
        if( it->second->id_type()  == t_gdata::ALLOBS ) {
          t_gallobs* all_obs = (t_gallobs*) it->second;
          set<string> allStations = all_obs->stations();
          if(allStations.empty()) {
              continue;
          }
          string markname = *(allStations.begin());

          vector<t_gtime> epochs = all_obs->epochs(markname);
          vector<t_gtime>::iterator epochsIT = epochs.begin();
          // for every pair "EPOCH record, OBSERVATION records"
          for( epochsIT = epochs.begin();
                  epochsIT != epochs.end(); ++epochsIT) {
            /* map from Satellite number from OBS record into a structure
            * that holds the obs, lli, sig strength trios */  
            map<string, shared_ptr<t_gobsgnss> > obsRecords
                    = all_obs->find(markname, *epochsIT);
            int yearEP, monthEP, doMEP;  
            (*epochsIT).ymd(yearEP, monthEP, doMEP);
            int hourEP, minuteEP, secondEP;
            (*epochsIT).hms(hourEP, minuteEP, secondEP);
            double dsec = (*epochsIT).dsec();
            string secdsecStr = trim(dbl2str(((double)secondEP)+dsec, 7));
            _ss << '>' <<  ' ' << yearEP
                    << ' ' << std::setfill('0') << setw(2) << monthEP
                    << ' ' << std::setfill('0') << setw(2) << doMEP
                    << ' ' << std::setfill('0') << setw(2) << hourEP
                    << ' ' << std::setfill('0') << setw(2) << minuteEP
                    << string(11 - secdsecStr.size(), ' ') << secdsecStr
                    << "  "
                    /* _flag was not saved in rinexo3::_read_epoch()!! 
                     Thus cannot distinguish between 0 and 1 and assuming 0
                     */
                    << '0'
                    << std::setfill(' ') << setw(3) << obsRecords.size() << endl;

            // for every OBSERVATION record (which means for every line)        
            map<string, shared_ptr<t_gobsgnss> >::iterator obsRecIT;
            for(obsRecIT =  obsRecords.begin(); obsRecIT != obsRecords.end(); ++obsRecIT) 
            {
                string satNum = obsRecIT->first;
                _ss << satNum;
                shared_ptr<t_gobsgnss> currObsRec = obsRecIT->second;
                t_rnxhdr::t_vobstypes obsDescrs = mapobs[satNum.substr(0,1)];
                t_rnxhdr::t_vobstypes::iterator obsDesIT;
                // for each trio obs, lli, sig strength
                for(obsDesIT = obsDescrs.begin(); obsDesIT != obsDescrs.end(); ++obsDesIT) 
                {
                    double fact = (double) (obsDesIT->second);
                    GOBS currObsDescr = obsDesIT->first;
                    double currObs = currObsRec->getobs(currObsDescr);
                    if(double_eq(fact, 0.0) || double_eq(currObs, NULL_GOBS)) {
                        _ss << string(16, ' ');
                        continue;
                    }
                    
                    currObs /= fact;
                    
                    string currObsStr = trim(dbl2str(currObs, 3));
                    _ss << string(14 - currObsStr.size(), ' ') << currObsStr;
                    
                    // now need lli and sig strength
                    // same condition as in t_rinexo2::_read_obs
                    if( ( GOBS(currObsDescr) >=  300 && GOBS(currObsDescr)  < 400 ) ||
                        ( GOBS(currObsDescr) >= 1300 && GOBS(currObsDescr) < 1400 ) ) {
                        int currLli = currObsRec->getlli(currObsDescr);
                        _ss << currLli;
                        
                        // again, same condition as in t_rinexo2::_read_obs
                        if( GOBS(currObsDescr) >= 100 && GOBS(currObsDescr) < 200 ){
                            GOBS snrtype = pha2snr(currObsDescr);
                            
                            double currSigStrengthConverted
                                = currObsRec->getobs(snrtype);
                            if(double_eq(NULL_GOBS, currSigStrengthConverted)) {
                                _ss << ' ';
                            } else {
                                if(double_eq(currSigStrengthConverted, 0.0)) {
                                    _ss << ' ';
                                } else if (double_eq(currSigStrengthConverted, 6.0)) {
                                    _ss << '1';
                                }else if (double_eq(currSigStrengthConverted, 15.0)) {
                                    _ss << '2';
                                }else if (double_eq(currSigStrengthConverted, 20.0)) {
                                    _ss << '3';
                                }else if (double_eq(currSigStrengthConverted, 27.0)) {
                                    _ss << '4';
                                }else if (double_eq(currSigStrengthConverted, 33.0)) {
                                    _ss << '5';
                                }else if (double_eq(currSigStrengthConverted, 39.0)) {
                                    _ss << '6';
                                }else if (double_eq(currSigStrengthConverted, 45.0)) {
                                    _ss << '7';
                                }else if (double_eq(currSigStrengthConverted, 50.0)) {
                                    _ss << '8';
                                }else if (double_eq(currSigStrengthConverted, 60.0)) {
                                    _ss << '9';
                                }else {
                                    _ss << ' ';
                                }
                            }
                        } else {
                            _ss << ' ';
                        }
                    } else {
                        _ss << string(2, ' ');
                    }
                }
                _ss << endl;
            }
          }
          break; // break out of for cycle going through _data
        }
    }
    //}  
  }

  int size = _fill_buffer( buff, sz );

  _mutex.unlock(); return size;
}


// OBS-RINEX header
//   - read individual lines (tmpsize once used and consumed)
//   - read block of lines at once (tmpsize cummulated and then consumed)
// ----------
int t_rinexo::_decode_head()
{ 
  gtrace("t_rinexo::_decode_head");
   
  while( _complete && ( ( _tmpsize = t_gcoder::_getline( _line )) >= 0 ) ){

#ifdef DEBUG
    cout << "LINE: " << _tmpsize << " " << _consume << " :" << _line; cout.flush();
#endif

    if( _tmpsize <= 61 ) break;
     
    _complete = true;

    // -------- "RINEX VERSION" --------
    if( _line.find("RINEX VERSION",60) != string::npos ){          // first line
      if( _line[20] != 'O' ){      // A1
        if( _log ){ _log->comment(0,"rinexo","Error: not RINEX observation file");  }
        else{               cerr << "rinexo - Error: not RINEX observation file\n"; }
                                 mesg(GERROR,"Error: not RINEX observation file!");  ++_irc;
        return (_consume = -1);
      }
      _version = trim(_line.substr(0,9)); // F9.2
      substitute(_version, " ", "");

      _rnxhdr.rnxver(_version);
      _rnxhdr.rnxsys( toupper(_line[40]) );     // A1 (G=GPS, R=GLO, E=GAL, S=SBAS, M=Mix)

      if( _log && substitute(_version, " ", "") > 0 ){ 
        _log->comment(2,"rinexo", "reading VER: " + _version + " SYS: " + string(1,_rnxhdr.rnxsys()) );
      }
      
      if( _rnxhdr.rnxsys() == ' ' ){
        if( _log )  _log->comment(0, "rinexo", "Warning: RINEX SYS not defined, GPS set as default" );
        if( _version < "3.00" ){ mesg(GWARNING,"Warning: RINEX SYS not defined!");          }
        else{                    mesg(GERROR,    "Error: RINEX SYS not defined!");  ++_irc; }
        _rnxhdr.rnxsys( 'G' );
      }
      _csys = string(1,_rnxhdr.rnxsys());

    // -------- "END OF HEADER" --------
    }else if( _line.find("END OF HEADER",60) != string::npos ){
      if( _log ) _log->comment(2,"rinexo","END OF HEADER ");
      t_gcoder::_consume(_tmpsize);
      _rnxhdr.mapobs(_mapobs);
      _rnxhdr.mapcyc(_mapcyc);
      _rnxhdr.glofrq(_glofrq);
      _rnxhdr.globia(_globia);
      _fill_head();
      if(      _version < "3.00" ) t_rinexo2::_check_head();
      else if( _version < "4.00" ) t_rinexo3::_check_head();
      return(_consume = -1);
    }

    else if( _version < "3.00" ) t_rinexo2::_decode_head();
    else if( _version < "4.00" ) t_rinexo3::_decode_head();


    // add _comment to _rnxhdr
    _rnxhdr.comment(_comment); 

#ifdef DEBUG
    if( _complete ) cout << "consume(tmp,rinex::head): " << fixed << setw(3) << _tmpsize
	                                         <<  " " << right << setw(3) << _endpos
	                                         << ": " << _line; cout.flush();
    else            cout << "not a full record read"     << _tmpsize
  	                                                 << " tmpsize: " << _tmpsize << endl
	                                 << " consume: " << _consume << endl; cout.flush();
#endif

    // -------- CONSUME --------
    if( _complete ) _consume += t_gcoder::_consume(_tmpsize);    
    else break;
  }

#ifdef DEBUG
  cout << "consume(tot,rinex::head): " << _consume << " [" << _endpos << "]\n"; cout.flush();
#endif

  return _consume;
}


   
// OBS-RINEX body
//   - read block of lines for each epoch (tmpsize cummulated and then consumed)
// ----------
int t_rinexo::_decode_data()
{  
  gtrace("t_rinexo::_decode_data");
   
  while( _complete && ( ( _tmpsize = t_gcoder::_getline( _line ) ) >= 0 ) ){
    
    _complete = true;
    _vobs.clear();

#ifdef DEBUG
    cout << "LINE: " << _tmpsize << " " << _consume << " :" << _line; cout.flush();
#endif

    // -------- DECODE DATA -------
    if(      _version < "3.00" ) t_rinexo2::_decode_data();
    else if( _version < "4.00" ) t_rinexo3::_decode_data();

  } // loop over lines
   
#ifdef DEBUG
  cout << "consume(tot,rinex::data): " << _consume << endl; cout.flush();
#endif

  return _consume;
}

} // namespace
