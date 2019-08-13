
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

#include <stdio.h>
#include <math.h>

#include "gdata/grec.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_grec::t_grec()
  : t_gobj()
    // t_gdata::id_type(REC) // NEFUNGUJE? 
{
//  cout << "CONSTRUCTOR t_grec \n"; cout.flush();
  id_type(REC);
}

/*
// copy constructor
// ----------
t_grec::t_grec(const t_grec& obj)
{
  _type = obj.id_type();
  _id   = obj.id();
  _name = obj.name();
       
  vector<t_gtime> vTIM = obj.rec_id();
  vector<t_gtime>::iterator itTIM = vTIM.begin();
  while( itTIM != vTIM.end() ){
    this->rec( *itTIM, obj.rec( *itTIM ));
    itTIM++;
  }
}
*/

// destructor
// ----------
t_grec::~t_grec()
{ 
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
  _maprec.clear();
  _gmutex.unlock();
//  cout << "DESTRUCTOR t_grec \n"; cout.flush();
}

// set receiver name
// ----------
void t_grec::rec(string rec, const t_gtime& beg, const t_gtime& end)
{  
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  _rec(rec, beg, end);
   
  _gmutex.unlock(); return;
}
   

// set receiver name
// ----------
void t_grec::_rec(string rec, const t_gtime& beg, const t_gtime& end)
{  
  t_maprec::iterator it =  _maprec.find(beg);

  if( !(beg < end) ){
    string lg = "Warning: " + _id + " not valid end time (end<beg) for antenna (beg:" + beg.str_ymdhms() + " -> end:" + end.str_ymdhms() + ")";     
    if( _log ) _log->comment(0,"grec",lg );
    else               cerr << "grec: " << lg << endl;
    return;
  }
   
  // begin record
  if( it == _maprec.end() ){     // not exists
      _maprec[beg] = rec;
     
  }else{                         // record exists
    if( it->first == LAST_TIME ||
        it->second.empty() ){

       _maprec[beg] = rec;
    }else{
       string lg = "Warning: " + _id + " valid object record cannot be changed for receiver";
       if( _log && _overwrite ) _log->comment(0,"grec",lg );
       else               cerr << "grec: " << lg << endl;
       return;
    }
  }

  // control end of record (with new beg search)
  it = _maprec.find(beg); it++;
   
  // beg was last in map (add final empty record)
  if( it == _maprec.end() ){
    _maprec[end] = "";

  }else{                                            // process end according to next value
    if( end < it->first ){                          // only if end is smaller then existing
      if( fabs(it->first - end) > 3600 ){           // significantly smaller!
	 if( it->second.empty() ) _maprec.erase(it); // remove obsolete empty record
	 _maprec[end] = "";

      }else{                                        // too close to next record
 	string lg = "Warning: " + _id + " 'rec' end tied to the existing value " + end.str("%Y-%m-%d %H:%M:%S -> ") + it->first.str("%Y-%m-%d %H:%M:%S");	 
        if( _log ) _log->comment(0,"grec",lg );
        else               cerr << "grec: " << lg << endl;
      }
    }
  }

  // remove duplicated empty records
  t_maprec::iterator itNEW = _maprec.begin();
  t_maprec::iterator itOLD = itNEW;
  while( itOLD != _maprec.end() ){    
    if( ++itNEW != _maprec.end() ){
      if( ( itNEW->second.empty() && itOLD->second.empty() ) ||
          ( itNEW->first == LAST_TIME ) )
      {
        _maprec.erase( itNEW++ );
      }
    }
    itOLD = itNEW;
  }
 
  return;
}

// get receiver name (>=t)
// ----------
string t_grec::rec(const t_gtime& t) const
{  
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  string tmp = "";
  tmp = _rec(t);
   
  _gmutex.unlock(); return tmp;
}
   
// get receiver name (>=t)
// ----------
string t_grec::_rec(const t_gtime& t) const
{  

  t_maprec::const_iterator it = _maprec.upper_bound(t);
  if( it == _maprec.begin() ){ return ""; } // not found (or not exists!)
  it--;
   
#ifdef DEBUG
	cout << " FOUND " << id() 
         << " time : " << it->first.str("%Y-%m-%d %H:%M:%S[%T]")
         << " < "      <<         t.str("%Y-%m-%d %H:%M:%S[%T]")
         << " "        << it->second
         << endl;
#endif   
   
  return it->second;

/*        
  string tmp("");
  map<t_gtime,string>::const_reverse_iterator it = _maprec.rbeg();
  while( it != _maprec.rend() ){
    if( t < it->first ) break;
    tmp = it->second;
//#ifdef DEBUG
	cout << " REC SEARCH " << id() 
         << " time : " << it->first.str("%Y-%m-%d %H:%M:%S[%T]")
         << " < "      <<         t.str("%Y-%m-%d %H:%M:%S[%T]")
         << " "        << tmp
         << endl;
//#endif
    it++;
  }
  return tmp;
*/
}


// get time tags
// ----------
vector<t_gtime> t_grec::rec_id() const
{
#ifdef BMUTEX   
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutex.lock();
   
  vector<t_gtime> tmp;
  t_maprec::const_iterator itMAP = _maprec.begin();
  while( itMAP != _maprec.end() ){

    tmp.push_back( itMAP->first );
    itMAP++;
  }
  _gmutex.unlock(); return tmp;
}

// add rinex header
// ----------------
void t_grec::addhdr(const t_rnxhdr& hdr, const t_gtime& epo, string path)
{
  _gmutex.lock();

  if( _maphdr.find(epo) == _maphdr.end() ){
     _maphdr[epo] = hdr;
     _maphdr[epo].path(path);
  }
   
  t_grec::t_maphdr::iterator itHDR;
  for( itHDR = _maphdr.begin(); itHDR != _maphdr.end(); ++itHDR )
  {
    t_rnxhdr rnxhdr = itHDR->second;
  }
  
  _gmutex.unlock();
}

// get all rinex headr
// -------------------
t_grec::t_maphdr t_grec::gethdr()
{   
  return _maphdr;   
}


// get one rinex headr
// -------------------
t_rnxhdr t_grec::gethdr(const t_gtime& epo)
{
   _gmutex.lock();
   
   t_rnxhdr  rnxhdr = _gethdr(epo);
   
   _gmutex.unlock();
   return rnxhdr;
}


// check consistency
// ----------
void t_grec::compare(shared_ptr<t_grec> grec, const t_gtime& tt)
{
   gtrace("t_grec::compare");

#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
 
  _gmutex.lock();

  t_rnxhdr rnxhdr = grec->gethdr(tt);
  grec->fill_rnxhdr(rnxhdr);

  _gmutex.unlock();
  t_gobj::compare( grec, tt );
  _gmutex.lock();

  string old, alt;
  old = trim(_rec(tt)); alt = trim(grec->rec(tt));
  if( old != alt ){      
    if ( old.empty() ) {
       _rec(alt, tt);
       if( _log ) _log->comment(0,"grec","Warning: object " + _id + " completed (Receiver): " + alt);
    }else if( _overwrite ) {	    
       _rec(alt, tt);
          if( _log ) _log->comment(0,"grec","Warning: object " + _id + " modified (Receiver): " + old + " -> " + alt);
    }else if( _log ) _log->comment(1,"gobj","Warning: object " + _id + " setting does not match Rinex header (Receiver): " + old + " -> " + alt);
  }

  // add hdr from grec if not exists at TT
  if( _maphdr.find(tt) == _maphdr.end() ){
    t_maphdr oth_head = grec->gethdr();
    if( oth_head.find(tt) != oth_head.end() ){
      t_rnxhdr head = oth_head[tt];
       _maphdr[tt] = head;
    }
     
  }
  
  _gmutex.unlock(); return;
}

// fill data members form rinex header 
// ---------------------   
void t_grec::fill_rnxhdr(const t_rnxhdr& rnxhdr)
{
   _gmutex.lock();

   _fill_rnxhdr(rnxhdr);

   _gmutex.unlock(); return;
}      
     
// get one rinex headr
// -------------------
t_rnxhdr t_grec::_gethdr(const t_gtime& epo)
{
  t_rnxhdr  rnxhdr;
  t_maphdr::iterator itREC = _maphdr.upper_bound(epo);
   
  if( itREC == _maphdr.begin() ) {
     return rnxhdr;
  }
   

  --itREC;  

  return itREC->second;
}
   
// fill data members form rinex header 
// ---------------------
void t_grec::_fill_rnxhdr(const t_rnxhdr& rnxhdr)
{
   t_gtime epo = rnxhdr.first();
     
   if( _name.empty() ) _name = rnxhdr.markname(); // NOT IF EXISTS! does not need to be the same

   _domes = rnxhdr.marknumb();
   _eccxyz(rnxhdr.antxyz(), epo);
   _eccneu(rnxhdr.antneu(), epo);
   _ant(rnxhdr.anttype(),   epo);
   _rec(rnxhdr.rectype(),   epo);
   _crd(rnxhdr.aprxyz(),    epo);   
}
   
   
} // namespace
