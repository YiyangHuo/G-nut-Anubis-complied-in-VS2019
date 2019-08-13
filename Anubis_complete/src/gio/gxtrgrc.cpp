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

#include "gio/gxtrgrc.h"
#include "gset/gsetqc.h"


using namespace std;
using namespace pugi;

namespace gnut {

// Constructors
// ---------
t_gxtrgrc::t_gxtrgrc()
{
}
   
t_gxtrgrc::t_gxtrgrc(t_gsetbase* set,
             const string& pgm, const t_gtime& dt)
 : t_gxtrqc(set, pgm, dt)
{
  gtrace("t_gxtrgrc::t_gxtrgrc");

  _get_settings();      
  _qcdata = make_shared<t_gqcdata>();

}

// Destructor
// ---------
t_gxtrgrc::~t_gxtrgrc()
{
  gtrace("t_gxtrgrc::desctructor");
}


// ===================
// PROTECTED FUNCTIONS
// ===================
   
// Get settings
// -----------------------------
void t_gxtrgrc::_get_settings()
{
  gtrace("t_gxtrgrc::_get_settings");

  t_gxtrqc::_get_settings();
   
  _gkpi             = dynamic_cast<t_gsetqc*>(_set)->gkpi();

  _dV_lim           = dynamic_cast<t_gsetqc*>(_set)->dV_lim();
  _dH_lim           = dynamic_cast<t_gsetqc*>(_set)->dH_lim();
  _dG_lim           = dynamic_cast<t_gsetqc*>(_set)->dG_lim();

  if(  _gkpi && !_band ) _band = -1; // do always band testing
  if(  _gkpi && !_calc ) _calc = -1; // do always position calculation
  if(  _gkpi && !_summ ) _summ = -1; // maybe separate in the future (counts of epoch)

}

   
// Section: KPI pars
// -----------------
void t_gxtrgrc::_grckpi(ostringstream& os)
{
  gtrace("t_gxtrgrc::_grckpi");

  if( _gkpi > 0 ){ _section(os, "Key-parameter indicators", _gkpi);
  }else if( !_gkpi ) return;

  set<GSYS> set_sys = _gnav->systems();

  string tt(_ref.str_ymdhms() + " ");
  
  ostringstream os_v2;  

  // loop GSYS
  for(auto itSYS = set_sys.begin(); itSYS != set_sys.end(); ++itSYS ){
    string gs = t_gsys::gsys2str(*itSYS);
//  char   ch = t_gsys::gsys2char(*itSYS);

    // SETTINGS INFORMATION
    t_gtriple _pos_ref(_qcdata->qc_est.xyz[*itSYS]);
    shared_ptr<t_grec> setrec = _rec->grec(_site);
    if( setrec == 0 ){
       if( _log ) _log->comment(0,"gxtrgrc","Error: KPI no user receiver settings "+_site);
       else               cerr << "gxtrgrc - Error: KPI no user receiver settings "+_site << endl;
       setrec = make_shared<t_grec>();      
    }
    else{ _pos_ref = setrec->crd(_obs_beg); }
    
    t_gtriple _blh_ref;
    xyz2ell(_pos_ref, _blh_ref, false);

    string key = gs + "KPI";

    if( _gkpi >= 2 ){
      _setkey(os_v2, key, '#', tt);
      os_v2 << W(9)  << "dN[m]" <<  W(9) << "dE[m]"  << W(9) << "dU[m]"
            << W(9)  << "GDOP"  <<  W(9) << "FLAGS"  << endl;
    }

    for(auto itEPO  = _qcdata->qc_est.pos_line[*itSYS].begin();
             itEPO != _qcdata->qc_est.pos_line[*itSYS].end();
           ++itEPO )
    {
      t_gtriple dNEU;
      t_gtriple dXYZ( itEPO->second.xyz - _pos_ref );
      
      xyz2neu( itEPO->second.blh, dXYZ, dNEU);
     
      _qcdata->qc_kpi.nepo[*itSYS]++;

      string flag = "";

      if( itEPO->second.gdop > _dG_lim ){
             _qcdata->qc_kpi.GD_XX[*itSYS]++; flag += "G";
      }else{ _qcdata->qc_kpi.GD_OK[*itSYS]++; flag += "_"; }

      if( sqrt(dNEU[0]*dNEU[0] + dNEU[1]*dNEU[1]) > _dH_lim ){ 
             _qcdata->qc_kpi.dH_XX[*itSYS]++; flag += "H";
      }else{ _qcdata->qc_kpi.dH_OK[*itSYS]++; flag += "_"; }
                                                                      
      if( dNEU[2] > _dV_lim ){
             _qcdata->qc_kpi.dV_XX[*itSYS]++; flag += "V";
      }else{ _qcdata->qc_kpi.dV_OK[*itSYS]++; flag += "_"; }

      if( _gkpi >= 2 ){
        _setkey(os_v2, key, ' ', itEPO->first.str_ymdhms());
        os_v2 << " " << fixed
              << setprecision(3) << W(9) << dNEU[0]
                                 << W(9) << dNEU[1]
                                 << W(9) << dNEU[2]
              << setprecision(2) << W(9) << _qcdata->qc_est.pos_line[*itSYS][itEPO->first].gdop
                                 << W(9) << flag
              << endl;
      }
    }

    ostringstream ostr;
    _setkey(os, "GNSKPI", '#', tt);
    os << W(9) << "CountOK" << W(9) << "CountXX" << W(9) << "%_Expect" << W(9) << "%_Exist" << W(9) << "%_Posit" << endl;
 
    int posEpo = _qcdata->qc_kpi.GD_OK[*itSYS] + _qcdata->qc_kpi.GD_XX[*itSYS];     // tot # of estimable positions

//  cerr  << " Expect: " << _nEpoExpect[*itSYS] << "Exists: " << _nEpoExists[*itSYS] << endl;
    if( _nEpoExpect[*itSYS] == 0 || _nEpoExpect[*itSYS] )
    {

      _setkey(os, gs+"_EP", ' ', tt);
      os << fixed << setprecision(0)
         << W(9) << _qcdata->qc_est.epo_used[*itSYS]
         << W(9) << _qcdata->qc_est.epo_excl[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_est.epo_used[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_est.epo_used[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << "-" << endl;

      _setkey(os, gs+"_DH", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.dH_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.dH_XX[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_kpi.dH_OK[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.dH_OK[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.dH_OK[*itSYS]*100.0 /  posEpo) << endl;

      _setkey(os, gs+"_DV", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.dV_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.dV_XX[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_kpi.dV_OK[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.dV_OK[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.dV_OK[*itSYS]*100.0 /  posEpo) << endl;

      _setkey(os, gs+"_GD", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.GD_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.GD_XX[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_kpi.GD_OK[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.GD_OK[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.GD_OK[*itSYS]*100.0 /  posEpo) << endl;

      _setkey(os, gs+"_DF", ' ', tt);
      os << W(9) << _qcdata->qc_kpi.DF_OK[*itSYS]
         << W(9) << _qcdata->qc_kpi.DF_XX[*itSYS]
         << W(9) << dbl2str(_qcdata->qc_kpi.DF_OK[*itSYS]*100.0 / _nEpoExpect[*itSYS])
         << W(9) << dbl2str(_qcdata->qc_kpi.DF_OK[*itSYS]*100.0 / _nEpoExists[*itSYS])
         << W(9) << "-" << endl;

      os << endl;
    }
    else{
      if( _log ) _log->comment(0,"gxtrgrc","Warning: KPI ["+gs+"] skipped, no epoch available"); 
       os                              << "# Warning: KPI ["+gs+"] skipped, no epoch available\n";
    }
  }
  
  if(_gkpi >= 2) os << os_v2.str();

}

} // namespace
