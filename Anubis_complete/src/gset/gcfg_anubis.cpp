
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

#include <iomanip>
#include <sstream>

#include "gset/gcfg_anubis.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gcfg_anubis::t_gcfg_anubis()
 : t_gsetbase(),
   t_gsetgen(),
   t_gsetinp(),
   t_gsetout(),
   t_gsetrec(),
   t_gsetgnss(),
   t_gsetproc(),
   t_gsetflt(),
   t_gsetqc()
{
  _IFMT_supported.insert(RINEXO_INP);
  _IFMT_supported.insert(RINEXN_INP);
  _IFMT_supported.insert(SP3_INP);
  _OFMT_supported.insert(LOG_OUT);
  _OFMT_supported.insert(PPP_OUT);
  _OFMT_supported.insert(XTR_OUT);
  _OFMT_supported.insert(XML_OUT);
  _OFMT_supported.insert(XQC_OUT);
  _OFMT_supported.insert(RINEXN_OUT);
  _OFMT_supported.insert(RINEXN2_OUT);
//_OFMT_supported.insert(RINEXO_OUT);
   
  _minimum_elev   = 0;
  _auto_band = true;

}


// Destructor
// ----------
t_gcfg_anubis::~t_gcfg_anubis()
{}


// settings check
// ----------
void t_gcfg_anubis::check()
{
  t_gsetgen::check();
  t_gsetinp::check();
  t_gsetout::check();
  t_gsetrec::check();
  t_gsetgnss::check();
  t_gsetproc::check();
  t_gsetflt::check();
  t_gsetqc::check();
}

// settings help
// ----------
void t_gcfg_anubis::help()
{
  t_gsetbase::help_header();
  t_gsetgen::help();
  t_gsetqc::help();
  t_gsetinp::help();
  t_gsetout::help();
  t_gsetgnss::help();
  t_gsetrec::help();
//t_gsetproc::help();
//t_gsetflt::help();
  t_gsetbase::help_footer();
}

} // namespace
