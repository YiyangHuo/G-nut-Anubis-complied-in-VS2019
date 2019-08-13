
#ifndef GSTAT_H
#define GSTAT_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: various statistics function
  Version: $ Rev: $

  2014-04-18 /PV: created

-*/

#include "gdata/gdata.h"

#include <vector>

#define CONF_INTERV 3       // confident interval factor
//#define IQR_SIG_LEV 2.5     // default setting of sig. level for box plot identification: 

using namespace std;

namespace gnut {

class t_gstat : public t_gdata 
{
 public:
   
   t_gstat(double conf_interv = CONF_INTERV);
   t_gstat(vector<double>& data, double conf_interv = CONF_INTERV);
//   t_gstat(vector<double>& data, double conf_interv = CONF_INTERV, double iqr_sig_lev = IQR_SIG_LEV);   

   void add_data(vector<double>& data);

   int robust_mean();
   int calc_median();
   int calc_minmax();
   int calc_stat(double std_ext = 0.0);
   int calc_lowuppq(double& lowq, double& uppq);
   int calc_iqrlimits(double& lowb, double& uppb);
   double calc_mean();
   double calc_std();
   double calc_mad();
   double calc_var();
   double calc_rms();
   double calc_ros();   
   double calc_iqr();      

   double get_min();
   double get_max();
   double get_std();
   double get_mad();
   double get_var();
   double get_rms();
   double get_mean();
   double get_median();
   double get_ros();
   double get_iqr();
   int    get_size();
   int    get_outl();

 protected:
   
   void   _add_data(vector<double>& data);
   double _p(double v);
   
   vector<double> _data;
   double _cinterv;
   double _iqrsig;
   double _std;
   double _mad;
   double _var;   
   double _rms;
   double _mean;
   double _median;
   double _min;
   double _max;
   double _lowq;
   double _uppq;
   double _ros;
   double _iqr;
   double _lowb;
   double _uppb;   
   int    _n_outl;
};
   
} // namespace

#endif
