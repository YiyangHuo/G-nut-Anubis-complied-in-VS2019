
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
#include <algorithm>

#include "gutils/gstat.h"
#include "gutils/gcommon.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {

// Constructor
// ----------
t_gstat::t_gstat(double conf_interv)
: _cinterv(conf_interv)
{
   gtrace("t_gstat::constructor");
   vector<double> tmp;
   _add_data(tmp);
   
   _rms = _std = _mean = _median = _min = _max = 0.0;
}

// Constructor
// ----------
t_gstat::t_gstat(vector<double>& data, double conf_interv)
: _cinterv(conf_interv)
{  
   gtrace("t_gstat::constructor");
   _add_data(data);
}
/*   
// Constructor
// ----------
 t_gstat::t_gstat(vector<double>& data, double conf_interv, double iqr_conf_lev)
 : _cinterv(conf_interv),
   _iqrsig(iqr_conf_lev)
{  
   gtrace("t_gstat::constructor");
   _add_data(data);
}   
*/
// add new data (and reset stats)
// ----------
void t_gstat:: add_data(vector<double>& data)
{
  _add_data(data);
}


// robust mean
// ----------
int t_gstat::robust_mean()
{
   gtrace("t_gstat::robust_mean");
   if( _data.size() < 2 ) return -1;
   
   int n = _data.size();   
   vector<double>::iterator it;

   // simple mean
   double count = 0;
   for( it = _data.begin(); it != _data.end(); ++it ) count += *it;
   double mean = count/n;
   
   double temp = 0;
   int c = 0;  // iteration number
   double std0, std;
   std0 = std = 0.0;
   while (fabs(temp - mean) > 0.01 && c <= 20){
     temp = mean;      
     // residuals and weighted mean
     count = 0;
     long double sump = 0;
     for( it = _data.begin(); it != _data.end(); ++it ){
        double v = fabs(*it - mean);
        sump += _p(v);
        count += _p(v)*(*it);
     }
     mean = count/sump;   
     // std
     count = 0;
     sump  = 0;
     for( it = _data.begin(); it != _data.end(); ++it ){
        double v = fabs(*it - mean);
        sump += _p(v);
        count += _p(v)*v*v;
     }      
     std0 = sqrt(count/(n-1));
     std = std0/sqrt(sump);
     c++;
   }
   
   _mean = mean;
   _std  = std;

   return 0;
}


// statistics
// ----------
int t_gstat::calc_stat(double std_ext)
{
   gtrace("t_gstat::calc_stat");
   if( _data.size() < 2 ) return -1;

   double max_res  = 0;
   double max_idx  = 0;
   double treshold = 0;
   
   this->calc_median();
   _mean = _median;

   do {
      vector<double> res;
      double sum = 0;
      for( size_t i = 0; i < _data.size(); ++i ){
         res.push_back(_data[i] - _mean);
         sum += pow(_data[i] - _mean, 2);
      }

      _std = sqrt(sum/res.size());
      
      if( !double_eq(std_ext, 0.0) ) treshold = _cinterv * std_ext;
      else                           treshold = _cinterv * _std;

      max_res = 0;
      max_idx = 0;

      for( size_t i = 0; i < _data.size(); ++i ){
	   
         if( fabs(res[i]) > max_res ){
           max_res = fabs(res[i]);
	   max_idx = i;
	 }	 
      }

      if ( max_res > treshold ) {
         _data.erase( _data.begin() + (long)max_idx );
	 _n_outl++;
      }      

      sum = 0;
      for( size_t i = 0; i < _data.size(); ++i ){
         sum += _data[i];
      }
      
      this->calc_median();
      _mean = _median;
      
   }while(max_res > treshold && _data.size() > 1 );

   _median = calc_median();
   _mean   = calc_mean();
   _std    = calc_std();
   _rms    = calc_rms();
   return 0;
}


// median
// ----------
int t_gstat::calc_median()
{
   gtrace("t_gstat::calc_median");
   _median = 0.0;

   vector<double> vec = _data;
   sort(vec.begin(), vec.end());
   
   int n = vec.size();
#ifdef DEBUG   
   for(int i = 0; i<n; i++)
     cout << i << "  " << vec[i] << endl;
#endif 

   if (n == 0) return -1;
   
   if (n%2 == 0) {
      int i = n/2;
      int j = i+1; 
   //   cout << i << "  " << j << "  " << vec[i-1] << "  " << vec[j-1] << endl;
      _median = (vec[i-1]+vec[j-1])/2;
   }else{
      int i = n/2 + 1; 
      _median = vec[i-1];
   }
   
   return 0;
}

// lower & upper quartile
int t_gstat::calc_lowuppq(double& lowq, double& uppq)
{     
  gtrace("t_stat::calc_lowuppq");
  _lowq = 0.0;
  _uppq = 0.0;
   
  vector<double> vec = _data;
  sort(vec.begin(), vec.end());
   
  int n = vec.size(); if(n==0) return -1;
   
  if(n%4 == 0){
    int idx_1 = n/4;
    int idx_2 = n/4 + 1;
    lowq = ( vec[idx_1] + vec[idx_2] )/2;
  }else if(n%4 != 0){
    double num = n/4.0;
    int idx = (int)ceil(num);
    lowq = vec[idx-1];
  }else{
     cout << "warning: Probem with lower quartile calculation!" << endl;
  }
   
  if(n%4 == 0){
    int idx_1 = 3*n/4;
    int idx_2 = 3*n/4 + 1;
    uppq = ( vec[idx_1] + vec[idx_2] )/2;
  }else if(n%4 != 0){
    double num = n/4.0;
    int idx = (int)ceil(3*num);
    uppq = vec[idx-1];
  }else{
    cout << "warning: Probem with upper quartile calculation!" << endl;
  }
   
  return 0;
}
   
 
// interquartile range
// -------------------
double t_gstat::calc_iqr()
{
  gtrace("t_gstat::calc_iqr");
   
  double iqr = 0.0;

  double min, max;
  this->calc_lowuppq(min, max);
   
  iqr = max - min;

  _iqr = iqr;
  return iqr;
}      

// interquartile limits
// --------------------
int t_gstat::calc_iqrlimits(double& lowb, double& uppb)
{
  gtrace("t_gstat::calc_iqrlimits");
   
  _lowb = 0.0;
  _uppb = 0.0;   

  this->calc_iqr();
  double min, max;
  this->calc_lowuppq(min, max);
   
  double condition = _cinterv*_iqr;
  lowb = min - condition;
  uppb = max + condition;
   
  _lowb = lowb;
  _uppb = uppb;   

  return 0;
}   

   
// simple mean
// ----------
double t_gstat::calc_mean()
{
   gtrace("t_gstat::calc_mean");

   double mean = 0.0;   
   double sum  = 0.0;
   
   if(_data.size() == 0) return mean;
   
   for( size_t i = 0; i < _data.size(); ++i ){
       sum += _data[i];
   }
   
   mean = sum/_data.size();

   _mean = mean;
   return mean;
}

// std of mean
// ----------
double t_gstat::calc_std()
{
   gtrace("t_gstat::calc_std");
   
   double std = 0.0;
   double sum = 0.0;
   
   this->calc_mean();
   
   for( size_t i = 0; i < _data.size(); ++i ){
      //_res.push_back(_data[i] - _mean);
      sum += pow(_data[i] - _mean, 2);
   }

   std = sqrt(sum/_data.size());
   
   _std = std;
   return std;
}

// mad of mean
// REF: Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median
// ----------
double t_gstat::calc_mad()
{
   gtrace("t_gstat::calc_mad");
   
//   double std = 0.0;
//   double sum = 0.0;
   
   this->calc_median();
   
   //the median is subtracted of each observation and becomes the series of absolute values
   vector<double> data0;
   for( size_t i = 0; i < _data.size(); ++i ){
	   data0.push_back(abs(_data[i] - _median));
   }   
   sort(data0.begin(), data0.end());

   //calculate the new median
   double median = 0.0;
   int n = data0.size();
   if (n%2 == 0) {
      int i = n/2;
      int j = i+1; 
      median = (data0[i-1]+data0[j-1])/2;
   }else{
      int i = n/2 + 1; 
      median = data0[i-1];
   }

   _mad = 1.4826*median;

   return _mad;
}
   
// variance
// --------
double t_gstat::calc_var()
{
   gtrace("t_gstat::calc_variance");
   
   double var = 0.0;

   this->calc_std();
   
   var = pow(_std,2.0);

   _var = var;
   return var;
}
   
   
// rms 
// ----------
double t_gstat::calc_rms()
{
   gtrace("t_gstat::calc_rms");
   
   double rms = 0.0;
   double sum = 0.0;   
   
   if(_data.size() <= 1) return rms;
   
   for( size_t i = 0; i < _data.size(); ++i ){
      sum += pow(_data[i], 2);
   }

   rms = sqrt(sum/_data.size());
   
   _rms = rms;
   return rms;
}
   
// get min & max val
// -----------------
int t_gstat::calc_minmax()
{
   gtrace("t_gstat::calc_minmax");
   _min = 0.0;
   _max = 0.0;
   
   if( _data.size() == 0 ) return -1;

   vector<double> vec = _data;
   sort(vec.begin(), vec.end());

   vector<double>::iterator i = vec.begin();
   vector<double>::reverse_iterator j = vec.rbegin();
   
//   _min = min = *i;
//   _max = max = *j;

   _min = *i;
   _max = *j;
   
   return 0;
}

// range of sample
// ---------------
double t_gstat::calc_ros()
{
  gtrace("t_gstat::calc_ros");
   
  double ros = 0.0;

  this->calc_minmax();
   
  ros = _max - _min;

  _ros = ros;
  return ros;
}   
   
// get values
// ----------
double t_gstat::get_median()
{
   return _median;
}

// get values
// ----------
double t_gstat::get_mean()
{
   return _mean;
}

// get values
// ----------
double t_gstat::get_std()
{
   return _std;
}
   
// get values
// ----------
double t_gstat::get_mad()
{
   return _mad;
}

// get values
// ----------
double t_gstat::get_var()
{
   return _var;
}   
  
// get values
// ----------
double t_gstat::get_rms()
{
   return _rms;
}   

// get values
// ----------

double t_gstat::get_min()
{
   return _min;
}

// get values
// ----------
double t_gstat::get_max()
{   
  return _max;
}   

// get values
// ----------
double t_gstat::get_ros()
{
   return _ros;
}
   
// get values
// ----------
double t_gstat::get_iqr()
{
   return _iqr;
}   

// get values
// ----------

int t_gstat::get_size()
{
   return _data.size();
}

// outliers
// ----------
int t_gstat::get_outl()
{
  return _n_outl; 
}

// add data
// ----------
void t_gstat::_add_data(vector<double>& data)
{
   gtrace("t_gstat::_add_data");

   _mean   = 0.0;
   _std    = 0.0;
   _median = 0.0;
   _n_outl = 0;
   _data.clear();
   _data   = data;   
}

// weights
// ----------
double t_gstat::_p(double v)
{
   return 1/sqrt(1+v*v/2);
}

   
} // namespace
