
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

#include "gutils/gtypeconv.h"
#include "gmodels/gstochasticmodel.h"

namespace gnut {

// Constant stochastic model

t_stochastic::t_stochastic()
{}

// Random walk stochastic model

t_randomwalk::t_randomwalk() : t_stochastic()
{
  _dSig = 9999.9;
   
  t_gtime Tfirst(1980, 01, 01, 00, 00, 00, t_gtime::GPS);
  setTprev(Tfirst);
  setTcurr(Tfirst);
}

void t_randomwalk::setq(double q)
{
  this->_dSig = q;
}


void t_randomwalk::setTprev(const t_gtime& Tprev)
{
  this->_Tprev = Tprev;
}

void t_randomwalk::setTcurr(const t_gtime& Tcurr)
{
  this->_Tcurr = Tcurr;
}

double t_randomwalk::getQ()
{
  double Q;
  double q = (_dSig*1e-3*_dSig*1e-3)/3600;
//   double q = (_dSig*1.e-3)/3600;

  Q = q*std::abs(_Tcurr - _Tprev);   
   
#ifdef DEBUG
   cout << "Current time: " << _Tcurr.str("%Y:%m:%d  %H:%M:%S") << endl;
   cout << "Previous time: " << _Tprev.str("%Y:%m:%d  %H:%M:%S") << endl;
   cout << "dT [s]: "<< _Tcurr - _Tprev << endl;
   cout << scientific << setprecision(3) << "_dSig: " << _dSig << endl;   
   cout << scientific << setprecision(3) << "q: " << q << endl;
   cout << scientific << setprecision(3) << "Q: " << Q << endl;
#endif   
   
   return Q;
}

double t_randomwalk::get_dt()
{
  return (_Tcurr - _Tprev);
}


void t_randomwalk::updateTime(const t_gtime& Tnew)
{
  setTprev(_Tcurr);
  setTcurr(Tnew);
}

// White noise stochastic model

t_whitenoise::t_whitenoise(double var) : t_stochastic()
{
  this->_var = var;
}

void t_whitenoise::setVar(double var)
{
  this->_var = var;
}

double t_whitenoise::getQ()
{
  return this->_var;
}

} // namespace
