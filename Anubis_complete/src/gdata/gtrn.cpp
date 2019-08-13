
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

#include "gdata/gtrn.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gtrn::t_gtrn()
  : t_gobj(),
    _channel(255)
    // t_gdata::id_type(TRN) // NEFUNGUJE? 
{
//  cout << "CONSTRUCTOR t_gtrn \n"; cout.flush();
  id_type(TRN);
}

/*
// copy constructor
// ----------
t_gtrn::t_gtrn(const t_gtrn& obj)
{
  _type = obj.id_type();
  _id   = obj.id();
  _name = obj.name();
       
  vector<t_gtime> vTIM = obj.trn_id();
  vector<t_gtime>::iterator itTIM = vTIM.begin();
  while( itTIM != vTIM.end() ){
    this->trn( *itTIM, obj.trn( *itTIM ));
    itTIM++;
  }
}
*/

// destructor
// ----------
t_gtrn::~t_gtrn()
{ 
//  boost::mutex::scoped_lock lock(_mutex);

//  _mapchk.clear();
//  cout << "DESTRUCTOR t_gtrn \n"; cout.flush();
}

/*
void t_gtrn::chk(int chk, const t_gtime& beg, const t_gtime& end)
{
   _mapchk[beg] = chk;
}



int  t_gtrn::chk(const t_gtime& t) const
{

}
*/

void t_gtrn::channel(int chk)
{
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
   _gmutex.lock();
   _channel = chk;
   _gmutex.unlock();

   return;
}

int  t_gtrn::channel() const
{
#ifdef BMUTEX
   boost::mutex::scoped_lock lock(_mutex);
#endif
   _gmutex.lock();
   int tmp = _channel;
   _gmutex.unlock();
   
   return tmp;
}

} // namespace
