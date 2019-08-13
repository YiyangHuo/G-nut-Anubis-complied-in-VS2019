
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

#include <cmath>
#include <iostream>
#include <iomanip>

#include "gutils/gstat3d.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

// Constructor
// ----------
t_gstat3d::t_gstat3d(vector<t_gtriple>& data, double conf_interv)
 : t_gstat(conf_interv)
{
  gtrace("t_gstat3d::constructor");

  _add_data(data);
}

// statistics
// ----------
int t_gstat3d::calc_stat()
{
  gtrace("t_gstat3d::calc_stat");
  if( _data3d.size() < 2 ) return -1;
   
  double max_res  = 0;
  double max_idx  = 0;
  double treshold = 0;
   
  this->calc_median(); 
  _mean3d = _median3d;

  do {      
    vector<t_gtriple> res;
    t_gtriple sum3(0.0,0.0,0.0);
    for( size_t i = 0; i < _data3d.size(); ++i ){
       t_gtriple diff(_data3d[i][0] - _mean3d[0],
        	      _data3d[i][1] - _mean3d[1], 
		      _data3d[i][2] - _mean3d[2]);
       res.push_back(diff);

#ifdef DEBUG       
       cout << fixed << "size: " << _data3d.size() << " iter/outl: " << _n_outl 
	    << " diff: " << diff[0]
	    << " data1: " << _data3d[i][0] 
	    << " mean1: " << _mean3d[0] << endl;
#endif
       sum3[0] += pow(diff[0], 2);
       sum3[1] += pow(diff[1], 2);
       sum3[2] += pow(diff[2], 2);
    }

    _std3d[0] = sqrt(sum3[0]/res.size());
    _std3d[1] = sqrt(sum3[1]/res.size());
    _std3d[2] = sqrt(sum3[2]/res.size());
      
    treshold = _cinterv * _std3d.norm();

    max_res = 0;
    max_idx = 0;

    for( size_t i = 0; i < _data3d.size(); ++i ){
       
#ifdef DEBUG       
         cout << fixed << setprecision(3) << i << " " << _data3d.size() 
	      << " " << _data3d[i][0] << " " << _data3d[i][1] << " " << _data3d[i][2] << endl;
#endif

      if( fabs(res[i].norm()) > max_res ){
#ifdef DEBUG
         cout << fixed << setprecision(3) << i << " res: " << res[i].norm() << " max: " << max_res << " !" << endl;
#endif	 
         max_res = fabs(res[i].norm());
         max_idx = i;
      }	 
    }
     
#ifdef DEBUF
    cout << fixed << "size: " << _data3d.size() << " iter/outl: " << _n_outl
         << " max_res = " << max_res << " " << treshold << " " << _std3d.norm() << endl;
#endif
    if( max_res > treshold ){
       _data3d.erase( _data3d.begin() + (long)max_idx );
       _data.erase(   _data.begin()   + (long)max_idx );
       _n_outl++;
    }      

    sum3[0] = sum3[1] = sum3[2] = 0.0;
    for( size_t i = 0; i < _data3d.size(); ++i ){
      sum3 = sum3 + _data3d[i];
    }
      
    this->calc_median();
    _mean3d = _median3d;

  }while(max_res > treshold && _data3d.size() > 1 );
   
  // mean after outl. detect.
  _mean3d = t_gstat3d::calc_mean();
  _std3d  = t_gstat3d::calc_std();

  _mean   = t_gstat::calc_mean();
  _std    = t_gstat::calc_std();
   
#ifdef DEBUG   
  cout << fixed << setprecision(3) 
       << " _mean3d: " << _mean3d  << " _std3d: " << _mean3d
       << " _mean: "   << _mean    << " _std: " << _std << endl;
#endif

  return 0;
}


// median
// ----------
int t_gstat3d::calc_median()
{
  gtrace("t_gstat3d::calc_median");
   
  // 1D
  t_gstat::calc_median();
     
  // 3D
  t_gstat stat1d;
  vector<double> a, b, c;
  for( size_t i = 0; i < _data3d.size(); ++i ){
      a.push_back(_data3d[i][0]);
      b.push_back(_data3d[i][1]);
      c.push_back(_data3d[i][2]);      
  }
  stat1d.add_data(a);
  if( stat1d.calc_median() == 0 ) _median3d[0] = stat1d.get_median();
  stat1d.add_data(b);
  if( stat1d.calc_median() == 0 ) _median3d[1] = stat1d.get_median();
  stat1d.add_data(c);
  if( stat1d.calc_median() == 0 ) _median3d[2] = stat1d.get_median();
   
  return 0;
}

// get values
// ----------
t_gtriple t_gstat3d::get_median3d()
{
  return _median3d;
}

// get values
// ----------
t_gtriple t_gstat3d::get_mean3d()
{
  return _mean3d;
}


// get values
// ----------
t_gtriple t_gstat3d::get_std3d()
{
  return _std3d;
}


// add data
// ----------
void t_gstat3d::_add_data(vector<t_gtriple>& data)
{
  gtrace("t_gstat3d::_add_data");
   
  // 1D initialization
  vector<double> data1d;
  for( unsigned int i = 0; i < data.size(); i++ )
     data1d.push_back(data[i].norm());

   t_gstat::_add_data(data1d);
   
  // 13 initialization
  _data3d.clear();
  _data3d = data;
   
  _mean3d[0] = 0.0;
  _mean3d[1] = 0.0;
  _mean3d[2] = 0.0;

  _std3d[0] = 0.0;
  _std3d[1] = 0.0;
  _std3d[2] = 0.0;   
   
  _median3d[0] = 0.0;
  _median3d[1] = 0.0;
  _median3d[2] = 0.0;   
}


// simple mean
// ----------
t_gtriple t_gstat3d::calc_mean()
{
   gtrace("t_gstat3d::calc_mean");
   t_gtriple sum(0.0,0.0,0.0);
   t_gtriple mean(0.0,0.0,0.0);
   
   for( size_t i = 0; i < _data3d.size(); ++i ){
      sum = sum + _data3d[i];
   }

   mean[0] = sum[0]/_data3d.size();
   mean[1] = sum[1]/_data3d.size();
   mean[2] = sum[2]/_data3d.size();
//   _mean = mean;
   return mean;
}


// std of mean
// ----------
t_gtriple t_gstat3d::calc_std()
{    
   gtrace("t_gstat3d::calc_std");
   t_gtriple std(0.0,0.0,0.0); 
   t_gtriple sum(0.0,0.0,0.0);
   for( size_t i = 0; i < _data.size(); ++i ){
      t_gtriple diff(_data3d[i][0] - _mean3d[0],
  	   	     _data3d[i][1] - _mean3d[1], 
  		     _data3d[i][2] - _mean3d[2] );
  
      sum[0] += pow(diff[0], 2);
      sum[1] += pow(diff[1], 2);
      sum[2] += pow(diff[2], 2);
   }
  
   std[0] = sqrt(sum[0]/_data3d.size());
   std[1] = sqrt(sum[1]/_data3d.size());
   std[2] = sqrt(sum[2]/_data3d.size());
   
  // _std = std;
   return std;
}

} // namespace
