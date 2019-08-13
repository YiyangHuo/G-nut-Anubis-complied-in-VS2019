
#ifndef GSYNC_H
#define GSYNC_H 

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: sync/filter epochs
  Version: $ Rev: $

  2015-12-17 /JD: created

-*/

#include <cmath>

#include "gio/glog.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {

bool time_sync(const t_gtime& epo, double smp, double scl, t_glog* glog);

} // namespace

#endif // # GSYNC_H
