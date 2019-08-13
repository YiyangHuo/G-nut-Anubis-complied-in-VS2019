
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

#include "gprod/gprodion.h"
#include "gmodels/ginterp.h"

using namespace std;

namespace gnut {

// constructor
// ----------
t_gprodion::t_gprodion(const t_gtime& t, shared_ptr<t_gobj> pt)
 : t_gprod(t, pt)
{ 
   gtrace("t_gprodcrd::t_gprodion");
   
   id_type(ION);
   id_group(GRP_PRODUCT);   
}

// destructor
// --------------
t_gprodion::~t_gprodion()
{
}

// set TEC value of a grid point
void t_gprodion::tec(const double& lat, const double& lon, const double& val)
{
  _gmutex.lock();

  if(val < 0 && _log) {
    _log->comment(1,"gprodion", "Warning: storing negative value of TEC not allowed!");
    _gmutex.unlock();
    return;
  }

  if((lat > 90 || lat < -90 || lon > 0 || lon < 360) && _log){
    _log->comment(1, "gprodion", "Warning: TEC grid point out of range in storing!");
    _gmutex.unlock();
    return;
  }

  _tec_map[lat][lon] = val;
  
#ifdef DEBUG
  for(map<double, map<double, double>>::iterator itA = _tec_map.begin(); itA != _tec_map.end(); itA++)
  {
    for( map<double, double>::iterator itB = itA->second.begin(); itB != itA->second.end(); itB++)
    {
      cout << " Lat: " << itA->first 
           << " Lon: " << itB->first 
           << " Val: " << itB->second << endl;
    }
  }
#endif
  
  _gmutex.unlock();
  return;
}

// get TEC for particular location
double const t_gprodion::tec(const double& lat, const double& lon)
{
  _gmutex.lock();
  double val = 0.0;

  bool found_lat = false;
  bool found_lon = false;

  map<double, map<double, double>>::iterator itLAT1 = _tec_map.upper_bound(lat);
  map<double, map<double, double>>::iterator itLAT2;
  if (itLAT1 != _tec_map.end() && itLAT1 != _tec_map.begin()) {
    itLAT2 = itLAT1;
    itLAT2--;
    found_lat = true;
  }

  map<double, double>::iterator itLON1;
  if(found_lat) itLON1 = itLAT1->second.upper_bound(lon);
  map<double, double>::iterator itLON2;
  if (found_lat && itLON1 != itLAT1->second.end() && itLON1 != itLAT1->second.begin()) {
    itLON2 = itLON1;
    itLON2--;
    found_lon = true;
  }

  map<double, double>::iterator itLON3;
  map<double, double>::iterator itLON4;
  if (found_lat && itLON3 != itLAT2->second.end() && itLON3 != itLAT2->second.begin()) {
    itLON3 = itLAT2->second.upper_bound(lon);
    itLON4 = itLON3;
    itLON4--;
    found_lon = true;
  } else found_lon = false;

  map<t_gpair, double> grd_points;
//cout << "Found lat lon: " << found_lat << " " << found_lon << endl;
  if (found_lat && found_lon){

    t_gpair p1(itLAT1->first, itLON1->first);
    t_gpair p2(itLAT1->first, itLON2->first);
    t_gpair p3(itLAT2->first, itLON3->first);
    t_gpair p4(itLAT2->first, itLON4->first);
//cout << "MODEL epo: " << _epo.str_ymdhms() << endl;
    grd_points[p1] = itLON1->second; //cout << p1[0] << " : " << p1[1] << "  " << itLON1->second << endl;
    grd_points[p2] = itLON2->second; //cout << p2[0] << " : " << p2[1] << "  " << itLON2->second << endl;
    grd_points[p3] = itLON3->second; //cout << p3[0] << " : " << p3[1] << "  " << itLON3->second << endl;
    grd_points[p4] = itLON4->second; //cout << p4[0] << " : " << p4[1] << "  " << itLON4->second << endl;
  }else if(_log){
    ostringstream istr;
    istr << "Warning: TEC interpolation failed at the point (" << fixed << setprecision(5) << lat << "," << lon << ")!";
    _log->comment(1, "gprodion", istr.str());
  }

  t_ginterp interp;
  t_gpair point(lat, lon);
  if(interp.bilinear(grd_points, point, val) < 0 && _log) _log->comment(1, "t_gprodion", "Warning: TEC spatial interpolation failed!");

  _gmutex.unlock();
  return val;
}

// clear content
void t_gprodion::clear()
{
  _gmutex.lock();

  _tec_map.clear();

  _gmutex.unlock();
  return;
}


} // namespace
