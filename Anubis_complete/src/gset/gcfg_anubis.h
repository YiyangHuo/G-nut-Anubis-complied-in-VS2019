
#ifndef GCFG_ANUBIS_H
#define GCFG_ANUBIS_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements settings for RINEX conversion, editing, QC ...
  Version: $ Rev: $

  2012-12-03 /JD: created

-*/

#include <string>
#include <sstream>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <signal.h>
#include <algorithm>
#include <memory>

#if defined __linux__
#include <pthread.h>
#endif

#include "gio/gfile.h"
#include "gio/gxtrqc.h"
#include "gio/gxtrgrc.h"
#include "gcoders/sp3.h"
#include "gcoders/rinexo2.h"
#include "gcoders/rinexo3.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexn.h"

#include "gall/gallprec.h"
#include "gall/gallobs.h"
#include "gall/galloqc.h"

#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetgnss.h"
#include "gset/gsetproc.h"
#include "gset/gsetflt.h"
#include "gset/gsetqc.h"


using namespace std;
using namespace pugi;

namespace gnut {

class t_gcfg_anubis : public t_gsetgen,
                      public t_gsetinp,
                      public t_gsetout,
                      public t_gsetrec,
                      public t_gsetgnss,
                      public t_gsetproc,
                      public t_gsetflt,
                      public t_gsetqc
{

 public:
   t_gcfg_anubis();
  ~t_gcfg_anubis();

   void check();                                 // settings check
   void help();                                  // settings help

 protected:

 private:

};

} // namespace

#endif