
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
#include <iomanip>

#include "gutils/gmutex.h"

namespace gnut {

// constructor
// ----------
t_gmutex::t_gmutex()
{
#if  defined _WIN32 || defined _WIN64
  InitializeCriticalSection( &_mutex );
#else
  pthread_mutex_init( &_mutex , 0 );
#endif
}


// destructor
// ----------
t_gmutex::~t_gmutex()
{
#if  defined _WIN32 || defined _WIN64
  DeleteCriticalSection( &_mutex );
#else
  pthread_mutex_destroy( &_mutex );
#endif
}


// ----------
void t_gmutex::lock()
{   
#if  defined _WIN32 || defined _WIN64
  EnterCriticalSection( &_mutex );
#else
  pthread_mutex_lock( &_mutex );
#endif
}
 

// ----------
void t_gmutex::unlock()
{
#if  defined _WIN32 || defined _WIN64
  LeaveCriticalSection( &_mutex );
#else
  pthread_mutex_unlock( &_mutex );
#endif
}


// ----------
#if  defined _WIN32 || defined _WIN64
CRITICAL_SECTION* t_gmutex::get_mutex()
#else
pthread_mutex_t* t_gmutex::get_mutex()
#endif
{
  return &_mutex;
}


/*
// TEMPORARY ?????
// ----------
t_glock::lock( mutex_type &mobj ) 
 : m_pMutex( &mobj ) 
{ 
  m_pMutex->Lock(); 
}


// ----------
t_glock::~glock() 
{ 
  m_pMutex->Unlock(); 
}
*/
   

} // namespace
