
#ifndef GQCDATA_H
#define GQCDATA_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (services@gnutsoftware.com)

  Purpose: container for QC data
  Version: $Rev:$

  2018-05-20 /JD: created

-*/

#include <map>
#include <vector>

#include "gdata/gdata.h"
#include "gdata/gsatdata.h"
#include "gutils/gtime.h"
#include "gutils/gnss.h"

using namespace std;

namespace gnut {

struct t_qc_est {
  struct t_pos_line {
    t_gtriple xyz,blh;
    double gdop,pdop,hdop,vdop,clk;
    int nsat, xsat;
  };
  double                sampling;
  t_gtime               epo_beg, epo_end;
  map<GSYS,t_gtriple>   xyz, xyz_rep, blh, blh_rep;
  map<GSYS, int>        epo_used, epo_excl;

  map<GSYS,map<t_gtime,t_pos_line>> pos_line;
};

struct t_qc_kpi {
  map<GSYS,int>         nepo;
  map<GSYS,double>      dH_OK, dH_XX, dV_OK, dV_XX, GD_OK, GD_XX, DF_OK, DF_XX;
};

// time-dependent data
struct t_qc_tim {
  struct t_qctimv {
    double ele, azi;
    int    lbn, cbn;
    int    health;
    map<GOBS, double>           snr;
    map<GOBS, pair<double,int>> mpt;
  };

  map<GSYS, map<t_gtime, map<string, t_qctimv>>> dat;
};

struct t_qc_bnd {
  map<GSYS, map<t_gtime, map<string, pair<int,int>>>> data;
};


class t_gqcdata : public t_gdata 
{
 public:
  t_gqcdata();
 ~t_gqcdata();

  t_gtime    qc_epo;
  t_qc_est   qc_est;
  t_qc_tim   qc_tim;
  t_qc_bnd   qc_bnd;
  t_qc_kpi   qc_kpi;

};

} // namespace

#endif
