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

#include <algorithm>

#include "gutils/gcommon.h"
#include "gutils/gconst.h"
#include "gmodels/gnwmbase.h"
#include "gset/gsetnwm.h"
#include "gset/gsetout.h"

using namespace std;

namespace gnut {   

// constructor
// ----------
t_gnwmbase::t_gnwmbase()
: t_gdata()
{
  gtrace("t_gnwmbase::constructor");
  id_type(t_gdata::NWMBASE);
  id_group(t_gdata::GRP_PRODUCT);
}

// constructor
// ----------
t_gnwmbase::t_gnwmbase(set<MET_DATA> data, set<MET_SURF> surf)
: t_gdata()
{
  gtrace("t_gnwmbase::constructor(set)");
  id_type(t_gdata::NWMSURF);
  id_group(t_gdata::GRP_PRODUCT);   
}


// destructor
// ----------
t_gnwmbase::~t_gnwmbase(){}


// get surf name id
// ----------
string t_gnwmbase::met_surf_str(const MET_SURF& type, bool full)
{
  switch( type )
  {    
    case LAT   : if (full) return "LAT";    else return "LAT";
    case LON   : if (full) return "LON";    else return "LON";
    case HEL   : if (full) return "HEL";    else return "H";
    case XPRJ  : if (full) return "XPRJ";   else return "X";
    case YPRJ  : if (full) return "YPRJ";   else return "Y";
    case GEOID : if (full) return "GEOID";  else return "N";
    case OROGR : if (full) return "OROGR";  else return "SF";
    case T2M   : if (full) return "T2M";    else return "T2";
    case TSK   : if (full) return "TSK";    else return "TSK";
    case TCWV  : if (full) return "TCWV";   else return "TCWV";
    case LSM   : if (full) return "LSM";    else return "LSM";
    case MSL   : if (full) return "MSL";    else return "MSL";

    // surface adjusted values
    case TSURF : if (full) return "TSURF";  else return "Ts";  // surface TEMP (from fitting)
    case ESURF : if (full) return "ESURF";  else return "Es";  // surface EPRE (from fitting)
    case QSURF : if (full) return "QSURF";  else return "Qs";  // surface SHUM (from fitting)

    // surface model values
    case MSURF : if (full) return "MSURF";  else return "Ms";  // surface TM   (model calculation)
    case WSURF : if (full) return "WSURF";  else return "Ws";  // surface ZWD  (model calculation)

    // adopted (default) values for lapse/decay rate (see settings: vert_fit + vert_adj)
    case dT    : if (full) return "DT";     else return "dT";  // temperature lapse rate  [K/km]
    case dM    : if (full) return "DM";     else return "dM";  // mean temper lapse rate  [K/km]
    case dE    : if (full) return "DE";     else return "dE";  // WV  exponential decay rate [-]
    case dW    : if (full) return "DW";     else return "dW";  // ZWD exponential decay rate [-]
    case dI    : if (full) return "DI";     else return "dI";  // IWV exponential decay rate [-]
    case dR    : if (full) return "DR";     else return "dR";  // ZWD/WV ratio decay rate [%/km]

    case eqH   : if (full) return "EQH";    else return "eqH"; // P/ZHD equiv height [km]
    case eqW   : if (full) return "EQW";    else return "eqW"; // ZWD   equiv height [km]
    case eqE   : if (full) return "EQE";    else return "eqE"; // WV    equiv height [km]

    case scH   : if (full) return "SCH";    else return "scH"; // P/ZHD scale height [km]
    case scW   : if (full) return "SCW";    else return "scW"; // ZWD   scale height [km]
    case scE   : if (full) return "SCE";    else return "scE"; // WV    scale height [km]

    // =================================
    // two parameter fit (incl. surface)
    // =================================
    case MS_Ep : if (full) return "MS_Ep";  else return "dE";
    case MS_Tp : if (full) return "MS_Tp";  else return "dT";
    // =================================

    case ME_Ph : if (full) return "ME_Ph";  else return "Ph"; // CONTROL ==> SHOULD BE 1
    case ME_Eh : if (full) return "ME_Eh";  else return "dE";
    case ME_Wh : if (full) return "ME_Wh";  else return "dW";

    case MP_Tp : if (full) return "MP_Tp";  else return "dT";
    case SP_Tp : if (full) return "SP_Tp";  else return "dT";
    case AP_Tp : if (full) return "AP_Tp";  else return "dT";
    case ML_Th : if (full) return "ML_Th";  else return "dT";
    case SL_Th : if (full) return "SL_Th";  else return "dT";
    case AL_Th : if (full) return "AL_Th";  else return "dT";

    case MP_Mp : if (full) return "MP_Mp";  else return "dTm";
    case SP_Mp : if (full) return "SP_Mp";  else return "dTm";
    case AP_Mp : if (full) return "AP_Mp";  else return "dTm";
    case ML_Mh : if (full) return "ML_Mh";  else return "dTm";
    case SL_Mh : if (full) return "SL_Mh";  else return "dTm";
    case AL_Mh : if (full) return "AL_Mh";  else return "dTm";

    case MP_Ep : if (full) return "MP_Ep";  else return "dE";
    case SP_Ep : if (full) return "SP_Ep";  else return "dE";
    case AP_Ep : if (full) return "AP_Ep";  else return "dE";
    case ML_Ep : if (full) return "ML_Ep";  else return "dE";
    case SL_Ep : if (full) return "SL_Ep";  else return "dE";
    case AL_Ep : if (full) return "AL_Ep";  else return "dE";

    case MP_Eh : if (full) return "MP_Eh";  else return "dE";
    case SP_Eh : if (full) return "SP_Eh";  else return "dE";
    case AP_Eh : if (full) return "AP_Eh";  else return "dE";
    case ML_Eh : if (full) return "ML_Eh";  else return "dE";
    case SL_Eh : if (full) return "SL_Eh";  else return "dE";
    case AL_Eh : if (full) return "AL_Eh";  else return "dE";

    case LR_Ep : if (full) return "LR_Ep";  else return "dE";
    case LM_Ep : if (full) return "LM_Ep";  else return "dE";
    case ER_Ep : if (full) return "ER_Ep";  else return "dE";
    case LL_Ep : if (full) return "LL_Ep";  else return "dE";
    case LR_Eh : if (full) return "LR_Eh";  else return "dE";
    case LM_Eh : if (full) return "LM_Eh";  else return "dE";
    case ER_Eh : if (full) return "ER_Eh";  else return "dE";
    case LL_Eh : if (full) return "LL_Eh";  else return "dE";

    case MP_Fp : if (full) return "MP_Fp";  else return "dE";
    case SP_Fp : if (full) return "SP_Fp";  else return "dE";
    case AP_Fp : if (full) return "AP_Fp";  else return "dE";
    case ML_Fp : if (full) return "ML_Fp";  else return "dE";
    case SL_Fp : if (full) return "SL_Fp";  else return "dE";
    case AL_Fp : if (full) return "AL_Fp";  else return "dE";

    case MP_Fh : if (full) return "MP_Fh";  else return "dE";
    case SP_Fh : if (full) return "SP_Fh";  else return "dE";
    case AP_Fh : if (full) return "AP_Fh";  else return "dE";
    case ML_Fh : if (full) return "ML_Fh";  else return "dE";
    case SL_Fh : if (full) return "SL_Fh";  else return "dE";
    case AL_Fh : if (full) return "AL_Fh";  else return "dE";

    case LR_Fp : if (full) return "LR_Fp";  else return "dE";
    case LM_Fp : if (full) return "LM_Fp";  else return "dE";
    case ER_Fp : if (full) return "ER_Fp";  else return "dE";
    case LL_Fp : if (full) return "LL_Fp";  else return "dE";
    case LR_Fh : if (full) return "LR_Fh";  else return "dE";
    case LM_Fh : if (full) return "LM_Fh";  else return "dE";
    case ER_Fh : if (full) return "ER_Fh";  else return "dE";
    case LL_Fh : if (full) return "LL_Fh";  else return "dE";

    case MP_Wp : if (full) return "MP_Wp";  else return "dW";
    case SP_Wp : if (full) return "SP_Wp";  else return "dW";
    case AP_Wp : if (full) return "AP_Wp";  else return "dW";
    case ML_Wp : if (full) return "ML_Wp";  else return "dW";
    case SL_Wp : if (full) return "SL_Wp";  else return "dW";
    case AL_Wp : if (full) return "AL_Wp";  else return "dW";

    case MP_Wh : if (full) return "MP_Wh";  else return "dW";
    case SP_Wh : if (full) return "SP_Wh";  else return "dW";
    case AP_Wh : if (full) return "AP_Wh";  else return "dW";
    case ML_Wh : if (full) return "ML_Wh";  else return "dW";
    case SL_Wh : if (full) return "SL_Wh";  else return "dW";
    case AL_Wh : if (full) return "AL_Wh";  else return "dW";

    case MP_Ip : if (full) return "MP_Ip";  else return "dI";
    case SP_Ip : if (full) return "SP_Ip";  else return "dI";
    case AP_Ip : if (full) return "AP_Ip";  else return "dI";
    case ML_Ip : if (full) return "ML_Ip";  else return "dI";
    case SL_Ip : if (full) return "SL_Ip";  else return "dI";
    case AL_Ip : if (full) return "AL_Ip";  else return "dI";

    case MP_Ih : if (full) return "MP_Ih";  else return "dI";
    case SP_Ih : if (full) return "SP_Ih";  else return "dI";
    case AP_Ih : if (full) return "AP_Ih";  else return "dI";
    case ML_Ih : if (full) return "ML_Ih";  else return "dI";
    case SL_Ih : if (full) return "SL_Ih";  else return "dI";
    case AL_Ih : if (full) return "AL_Ih";  else return "dI";

    case LR_Wp : if (full) return "LR_Wp";  else return "dW";
    case LM_Wp : if (full) return "LM_Wp";  else return "dW";
    case ER_Wp : if (full) return "ER_Wp";  else return "dW";
    case LL_Wp : if (full) return "LL_Wp";  else return "dW";
    case LR_Wh : if (full) return "LR_Wh";  else return "dW";
    case LM_Wh : if (full) return "LM_Wh";  else return "dW";
    case ER_Wh : if (full) return "ER_Wh";  else return "dW";
    case LL_Wh : if (full) return "LL_Wh";  else return "dW";

    case MP_WE : if (full) return "MP_WE";  else return "dE";
    case ML_WE : if (full) return "ML_WE";  else return "dE";
    case MP_Wi : if (full) return "MP_Wi";  else return "dE";
    case ML_Wi : if (full) return "ML_Wi";  else return "dE";
    case SP_Wi : if (full) return "SP_Wi";  else return "dE";
    case SL_Wi : if (full) return "SL_Wi";  else return "dE";
    case MP_Fi : if (full) return "MP_Fi";  else return "dE";
    case ML_Fi : if (full) return "ML_Fi";  else return "dE";
    case SP_Fi : if (full) return "SP_Fi";  else return "dE";
    case SL_Fi : if (full) return "SL_Fi";  else return "dE";

    case LR_WE : if (full) return "LR_WE";  else return "dE";
    case ER_Wi : if (full) return "ER_Wi";  else return "dE";
    case ER_Fi : if (full) return "ER_Fi";  else return "dE";
     
    case MP_Rt : if (full) return "MP_Rt";  else return "dR";
    case ML_Rt : if (full) return "ML_Rt";  else return "dR";
    case SP_Rt : if (full) return "SP_Rt";  else return "dR";
    case SL_Rt : if (full) return "SL_Rt";  else return "dR";
    case LR_Rt : if (full) return "LR_Rt";  else return "dR";
    case ER_Rt : if (full) return "ER_Rt";  else return "dR";
  }
  return "xxx";
}
// get surf name id
// ----------
string t_gnwmbase::met_surf_id(const MET_SURF& type)
{
  switch( type )
  {    
    case LAT   : return "lat[deg]";      //  north latitude [deg]
    case LON   : return "lon[deg]";      //  east longitude [deg]
    case HEL   : return "hel[m]";        //  surface elipsoidal height [m]
    case XPRJ  : return "X[m]";          //  projection coordinates [m]
    case YPRJ  : return "Y[m]";          //  projection coordinates [m]
    case GEOID : return "geo[m]";        //  geoid undulation (above WGS84) [m]
    case OROGR : return "oro[m]";        //  orography (above geoid) [m]
    case T2M   : return "t2m[K]";        //  temperature at 2 meters [K]
    case TSK   : return "tsk[K]";        //  skin temperature [K]
    case TCWV  : return "tcwv[mm]";      //  total column water vapour [mm]
    case LSM   : return "lsm[-]";        //  land-sea-mask [0/1]
    case MSL   : return "msl[hPa]";      //  mean-seal-level presssure [hPa]

    // adjusted surface values
    case TSURF : return "Ts[K]";         // surface TEMP (from fitting)  
    case ESURF : return "Es[hPa]";       // surface EPRE (from fitting)
    case QSURF : return "Qs[kg/kg]";     // surface SHUM (from fitting)

    // model calculation
    case MSURF : return "TMs[K]";        // surface TM   (model calculation)
    case WSURF : return "ZWDs[mm]";      // surface ZWD  (model calculation)

    // adopted (default) values for lapse/decay rate (see settings: vert_fit + vert_adj)
    case dT    : return "dT[K/km]";      // temperature lapse rate  [K/km]
    case dM    : return "dTm[K/km]";     // mean temper lapse rate  [K/km]
    case dE    : return "dE[-]";         // WV  exponential decay rate [-]
    case dW    : return "dW[-]";         // ZWD exponential decay rate [-]
    case dI    : return "dI[-]";         // IWV exponential decay rate [-]
    case dR    : return "dR[%/km]";      //     exponential decay rate [-]
     
    // scale heights
    case scH    : return "scH[km]";        // P/ZHD scale height [km]
    case scW    : return "scW[km]";        // ZWD   scale height [km]
    case scE    : return "scE[km]";        // WV    scale height [km]

    case eqH    : return "eqH[km]";        // P/ZHD equiv height [km]
    case eqW    : return "eqW[km]";        // ZWD   equiv height [km]
    case eqE    : return "eqE[km]";        // WV    equiv height [km]

    // two parameter fit (incl. surface)
    // =================================
    case MS_Tp : return "dT[K/km]";      //  temperature lapse rate [K/km]
    case MS_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
    // =================================
 
    case ME_Ph : return "dPh[-]";        //  pressure exponencial decay [-] --> 1 !
    case ME_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
    case ME_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]

    case MP_Tp : return "dT[K/km]";      //  temperature lapse rate [K/km]
    case SP_Tp : return "dT[K/km]";      //  temperature lapse rate [K/km]
    case AP_Tp : return "dT[K/km]";      //  temperature lapse rate [K/km]
    case ML_Th : return "dT[K/km]";      //  temperature lapse rate [K/km]
    case SL_Th : return "dT[K/km]";      //  temperature lapse rate [K/km]
    case AL_Th : return "dT[K/km]";      //  temperature lapse rate [K/km]

    case MP_Mp : return "dTm[K/km]";     //  mean temperature lapse rate [K/km]
    case SP_Mp : return "dTm[K/km]";     //  mean temperature lapse rate [K/km]
    case AP_Mp : return "dTm[K/km]";     //  mean temperature lapse rate [K/km]
    case ML_Mh : return "dTm[K/km]";     //  mean temperature lapse rate [K/km]
    case SL_Mh : return "dTm[K/km]";     //  mean temperature lapse rate [K/km]
    case AL_Mh : return "dTm[K/km]";     //  mean temperature lapse rate [K/km]

    case MP_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
    case SP_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
    case AP_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
    case ML_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
    case SL_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
    case AL_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]

    case MP_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
    case SP_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
    case AP_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
    case ML_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
    case SL_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
    case AL_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]

     case LR_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
     case LM_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
     case ER_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
     case LL_Ep : return "dEp[-]";        //  WV pressure exp. decay parameter [-]
     case LR_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
     case LM_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
     case ER_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]
     case LL_Eh : return "dEh[-]";        //  WV pressure exp. decay parameter [-]

    case MP_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
    case SP_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
    case AP_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
    case ML_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
    case SL_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
    case AL_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]

    case MP_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
    case SP_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
    case AP_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
    case ML_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
    case SL_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
    case AL_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
     
     case LR_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
     case LM_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
     case ER_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
     case LL_Fp : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
     case LR_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
     case LM_Fh : return "dFp[-]";        //  WV pressure exp. decay parameter [-]
     case ER_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
     case LL_Fh : return "dFh[-]";        //  WV pressure exp. decay parameter [-]
     
    case MP_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
    case SP_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
    case AP_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
    case ML_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
    case SL_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
    case AL_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]

    case MP_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
    case SP_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
    case AP_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
    case ML_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
    case SL_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
    case AL_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]

     case LR_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
     case LM_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
     case ER_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
     case LL_Wp : return "dWp[-]";        //  zenith wetD exp. decay parameter [-]
     case LR_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
     case LM_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
     case ER_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
     case LL_Wh : return "dWh[-]";        //  zenith wetD exp. decay parameter [-]
     
    case MP_Ip : return "dIp[-]";        //  integr.WV exp. decay parameter [-]
    case SP_Ip : return "dIp[-]";        //  integr.WV exp. decay parameter [-]
    case AP_Ip : return "dIp[-]";        //  integr.WV exp. decay parameter [-]
    case ML_Ip : return "dIp[-]";        //  integr.WV exp. decay parameter [-]
    case SL_Ip : return "dIp[-]";        //  integr.WV exp. decay parameter [-]
    case AL_Ip : return "dIp[-]";        //  integr.WV exp. decay parameter [-]

    case MP_Ih : return "dIh[-]";        //  integr.WV exp. decay parameter [-]
    case SP_Ih : return "dIh[-]";        //  integr.WV exp. decay parameter [-]
    case AP_Ih : return "dIh[-]";        //  integr.WV exp. decay parameter [-]
    case ML_Ih : return "dIh[-]";        //  integr.WV exp. decay parameter [-]
    case SL_Ih : return "dIh[-]";        //  integr.WV exp. decay parameter [-]
    case AL_Ih : return "dIh[-]";        //  integr.WV exp. decay parameter [-]

    case MP_WE : return "dWE[-]";        //  zenith wetD/e [-]
    case ML_WE : return "dWE[-]";        //  zenith wetD/e [-]
    case MP_Wi : return "dWi[-]";        //  zenith wetD exp. decay parameter [-]
    case ML_Wi : return "dWi[-]";        //  zenith wetD exp. decay parameter [-]
    case SP_Wi : return "dWi[-]";        //  zenith wetD exp. decay parameter [-]
    case SL_Wi : return "dWi[-]";        //  zenith wetD exp. decay parameter [-]
    case MP_Fi : return "dFi[km]";       //  FaF isotermic scale-height [km]
    case ML_Fi : return "dFi[km]";       //  FaF isotermic scale-height [km]
    case SP_Fi : return "dFi[km]";       //  FaF isotermic scale-height [km]
    case SL_Fi : return "dFi[km]";       //  FaF isotermic scale-height [km]
     
     case LR_WE : return "dWE[-]";        //  zenith wetD/e [-]
     case ER_Wi : return "dWi[-]";        //  zenith wetD exp. decay parameter [-]
     case ER_Fi : return "dFi[km]";       //  FaF isotermic scale-height [km]
     
    case MP_Rt : return "dRtX[%/m]";     //  Ex/Wx ratio exp.decay parameter [-/m]
    case ML_Rt : return "dRtL[%/m]";     //  El/Wl ratio exp.decay parameter [-/m]
    case SP_Rt : return "dRtX[%/m]";     //  Ex/Wx ratio exp.decay parameter [-/m]
    case SL_Rt : return "dRtL[%/m]";     //  El/Wl ratio exp.decay parameter [-/m]
      case LR_Rt : return "dRat[%/m]";     //  E/W ratio exp.decay parameter [-/m]
      case ER_Rt : return "dRat[%/m]";     //  Ex/Wx ratio exp.decay parameter [-/m]
  }
  return "xxx";
}


// get data name enum
// ----------
string t_gnwmbase::met_data_str(const MET_DATA& type, bool full)
{
  switch( type )
  {
    case GEOM  : if (full) return "GEOM"; else return "Ho";
    case GEOP  : if (full) return "GEOP"; else return "Hp";
    case GRAV  : if (full) return "GRAV"; else return "G";
    case PRES  : if (full) return "PRES"; else return "P";
    case EPRE  : if (full) return "EPRE"; else return "E";
    case SHUM  : if (full) return "SHUM"; else return "Q";
    case ZTD   : if (full) return  "ZTD"; else return "ZTD";
    case ZHD   : if (full) return  "ZHD"; else return "ZHD";
    case ZWD   : if (full) return  "ZWD"; else return "ZWD";
    case IWV   : if (full) return  "IWV"; else return "IWV";
    case TM    : if (full) return   "TM"; else return "Tm";
    case TEMP  : if (full) return "TEMP"; else return "T";
    case NGRD  : if (full) return "NGRD"; else return "NGRD";
    case EGRD  : if (full) return "EGRD"; else return "EGRD";
    case NGH   : if (full) return  "NGH"; else return "NGH";
    case NGW   : if (full) return  "NGW"; else return "NGW";
    case EGH   : if (full) return  "EGH"; else return "EGH";
    case EGW   : if (full) return  "EGW"; else return "EGW";
    case RATX  : if (full) return "RATX"; else return "RATIO";
    case RATL  : if (full) return "RATL"; else return "RATIO";
    case RATE  : if (full) return "RATE"; else return "RATIO";
    case RATW  : if (full) return "RATW"; else return "RATIO";
  }
  return "xxx";
}

// get data name id
// ----------
string t_gnwmbase::met_data_id(const MET_DATA& type)
{
  switch( type )
  {
    case GEOM  : return "Ho[m]";         //  geometric height [m]
    case GEOP  : return "G[gpm]";        //  geopotential height [gpm]
    case GRAV  : return "g[m/s2]";       //  gravity accelleration [m/s2]
    case ZWD   : return "ZWD[m]";        //  zenith wet delay (integrated) [m]
    case ZHD   : return "ZHD[m]";        //  zenith hydrostatic delay (integrated) [m]
    case ZTD   : return "ZTD[m]";        //  zenith total delay (integrated) [m]
    case IWV   : return "I[kg/m2]";      //  integrated water vapour [kg/m2]
    case TM    : return "Tm[K]";         //  mean temperature [K]
    case TEMP  : return "T[K]";          //  temperature [K]
    case PRES  : return "p[hPa]";        //  pressure [hPa]
    case EPRE  : return "e[hPa]";        //  water vapour pressure [hPa]
    case SHUM  : return "q[g/kg]";       //  specific humidity [g/kg]
    case NGRD  : return "NGR[mm]";       //  total tropospheric horizontal gradient [mm]
    case EGRD  : return "EGR[mm]";       //  total tropospheric horizontal gradient [mm]
    case NGH   : return "NGh[mm]";       //  hydrostatic tropospheric horizontal gradient [mm]
    case NGW   : return "NGw[mm]";       //  hydrostatic tropospheric horizontal gradient [mm]
    case EGH   : return "EGh[mm]";       //  wet tropospheric horizontal gradient [mm]
    case EGW   : return "EGw[mm]";       //  wet tropospheric horizontal gradient [mm]
    case RATL  : return "ratL[%]";       //  ratio of El/Wl exp.decay parameter
    case RATE  : return "ratE[%]";       //  ratio of Ex/Ex exp.decay parameter
    case RATW  : return "ratW[%]";       //  ratio of Wx/Wx exp.decay parameter
    case RATX  : return "ratX[%]";       //  ratio of Ex/Wx exp.decay parameter
  }
  return "xxx";
}


// get enum from data name
// ----------
int t_gnwmbase::met_data(string name)
{
  transform(name.begin(), name.end(), name.begin(), ::toupper);

  if(      name.compare("GEOM")  == 0 ) return GEOM;
  else if( name.compare("GEOP")  == 0 ) return GEOP;
  else if( name.compare("GRAV")  == 0 ) return GRAV;
  else if( name.compare("PRES")  == 0 || name.compare("P")  == 0 ) return PRES;
  else if( name.compare("EPRE")  == 0 || name.compare("E")  == 0 ) return EPRE;
  else if( name.compare("SHUM")  == 0 || name.compare("Q")  == 0 ) return SHUM;
  else if( name.compare("TEMP")  == 0 || name.compare("T")  == 0 ) return TEMP;
  else if( name.compare("TM")    == 0 || name.compare("M")  == 0 ) return TM;
  else if( name.compare("ZTD")   == 0 ) return ZTD;
  else if( name.compare("ZHD")   == 0 ) return ZHD;
  else if( name.compare("ZWD")   == 0 ) return ZWD;
  else if( name.compare("IWV")   == 0 ) return IWV;
  else if( name.compare("NGRD")  == 0 ) return NGRD;
  else if( name.compare("EGRD")  == 0 ) return EGRD;
  else if( name.compare("NGH")   == 0 ) return NGH;
  else if( name.compare("NGW")   == 0 ) return NGW;
  else if( name.compare("EGH")   == 0 ) return EGH;
  else if( name.compare("EGW")   == 0 ) return EGW;
  else if( name.compare("RATIO") == 0 || name.compare("RAT") == 0 || name.compare("RT") == 0 ) return RATX;
  return -1;
}


// get enum from surf name
// ----------
int t_gnwmbase::met_surf(string name)
{
  transform(name.begin(), name.end(), name.begin(), ::toupper);

  if(      name.compare("LAT")   == 0 ) return LAT;
  else if( name.compare("LON")   == 0 ) return LON;
  else if( name.compare("HEL")   == 0 ) return HEL;
  else if( name.compare("XPRJ")  == 0 ) return XPRJ;
  else if( name.compare("YPRJ")  == 0 ) return YPRJ;
  else if( name.compare("GEOID") == 0 ) return GEOID;
  else if( name.compare("OROGR") == 0 ) return OROGR;
  else if( name.compare("T2M")   == 0 ) return T2M;
  else if( name.compare("TSK")   == 0 ) return TSK;
  else if( name.compare("TCWV")  == 0 ) return TCWV;
  else if( name.compare("LSM")   == 0 ) return LSM;
  else if( name.compare("MSL")   == 0 ) return MSL;
  else if( name.compare("TSURF") == 0
	|| name.compare("TS")    == 0 ) return TSURF;
  else if( name.compare("MSURF") == 0
	|| name.compare("MS")    == 0 
        || name.compare("TMS")   == 0 ) return MSURF;
  else if( name.compare("ESURF") == 0
        || name.compare("ES")    == 0 ) return ESURF;
  else if( name.compare("WSURF") == 0 
	|| name.compare("WS")    == 0 ) return WSURF;
  else if( name.compare("QSURF") == 0 
	|| name.compare("QS")    == 0 ) return QSURF;
  else if( name.compare("DT")    == 0 ) return dT;
  else if( name.compare("DM")    == 0 
	|| name.compare("DTM")   == 0 ) return dM;
  else if( name.compare("DE")    == 0 ) return dE;
  else if( name.compare("DW")    == 0 ) return dW;
  else if( name.compare("DI")    == 0 ) return dI;
  else if( name.compare("DR")    == 0
	|| name.compare("DRAT")  == 0 ) return dR;
  else if (name.compare("EQH")   == 0 ) return eqH;
  else if (name.compare("EQW")   == 0 ) return eqW;
  else if (name.compare("EQE")   == 0 ) return eqE;
  else if (name.compare("SCH")   == 0 ) return scH;
  else if (name.compare("SCW")   == 0 ) return scW;
  else if (name.compare("SCE")   == 0 ) return scE;

  return -1;
}


// get data decimal resolution
// ----------
int t_gnwmbase::digit_dec(const MET_DATA& type)
{
  switch( type )
  {
    case GEOM : case GEOP : return 3;
    case GRAV :             return 5; 
    case ZTD  : case ZHD  : 
    case ZWD  : case IWV  : return 3;   // [m]
    case TM   : case TEMP : return 1;   // [K]
    case PRES : case EPRE : return 2;   // [hPa]
    case SHUM :             return 5; 
    case NGRD : case EGRD :
    case NGH  : case NGW  :
    case EGH  : case EGW  : return 3;   // [mm]
    case RATL : case RATE : return 3;   // [-]
    case RATW : case RATX : return 3;   // [-]
  }
  return 2;
}


// get data lenght
// ----------
int t_gnwmbase::digit_len(const MET_DATA& type)
{
  switch( type )
  {
    case GEOM : case GEOP : return 9; 
    case GRAV :             return 9; 
    case ZTD  : case ZHD  : 
    case ZWD  : case IWV  : return 6;   // [m]
    case TM   : case TEMP : return 6;   // [K]
    case PRES : case EPRE : return 8;   // [hPa]
    case SHUM :             return 6; 
    case NGRD : case EGRD :
    case NGH  : case NGW  :
    case EGH  : case EGW  : return 6;   // [mm]
    case RATL : case RATE : return 6;   // [-]
    case RATW : case RATX : return 6;   // [-]
  }
  return 9;
}


// get data decimal resolution
// ----------
int t_gnwmbase::digit_dec(const MET_SURF& type)
{
  switch( type )
  {	
   case LAT   : case LON  : return 3;
   case XPRJ  : case YPRJ : return 3;
   case OROGR : case HEL  : return 3;
   case GEOID : case MSL  : return 3;
   case T2M   : case TSK  : return 1;
   case TCWV  :             return 3;
   case LSM   :             return 0;
   case TSURF : case MSURF: return 1;
   case ESURF : case QSURF: return 2;
   case WSURF :             return 3;
   default    :             return 3;
  }
  return 2;
}


// get surf lenght
// ----------
int t_gnwmbase::digit_len(const MET_SURF& type)
{
  switch( type )
  {	
   case LAT   : case LON  : return 8;
   case XPRJ  : case YPRJ : return 8;
   case OROGR : case HEL  : return 9;
   case GEOID : case MSL  : return 7;
   case T2M   : case TSK  : return 6;
   case TCWV  :             return 6;
   case LSM   :             return 4;
   case TSURF : case MSURF: return 6;
   case ESURF : case QSURF: return 8;
   case WSURF :             return 6;
   default    :             return 6;
  }
  return 9;
}

// convert string set to MET_DATA & MET_SURF sets
// ----------
int t_gnwmbase::set_param(t_gsetbase* gset,
			  set<MET_SURF>& surf,
			  set<MET_DATA>& data )
{
  set<string> par = dynamic_cast<t_gsetnwm*>(gset)->param();        // selected parameters
  set<string> ass = dynamic_cast<t_gsetnwm*>(gset)->assess();       // selected assessment
  set<string> pro = dynamic_cast<t_gsetnwm*>(gset)->profile();      // selected profile
  
  NWMADJ      adj = dynamic_cast<t_gsetnwm*>(gset)->vert_adj();     // method of vertical adjustement
  NWMFIT      fit = dynamic_cast<t_gsetnwm*>(gset)->vert_fit();     // method of vertical approximation
  
  string     ofit = dynamic_cast<t_gsetout*>(gset)->outputs("fit"); // if fitting output requested only

#ifdef DEBUG
     cerr << "par_beg: " << par.size() << " ";
     for( set<string>::iterator it = par.begin(); it!= par.end(); ++it ) cerr << " " << *it;
     cerr << endl;
#endif

  set<string>::iterator it;
  for( it=par.begin(); it!=par.end(); ++it ){

    // all comparisons in CAPITAL!
    // ===========================
     
    // DATA - only if explicitely requested
    if( it->compare("ZWD")   == 0 ||   it->compare("ZTD")== 0    )  data.insert(ZWD);
    if( it->compare("ZHD")   == 0 ||   it->compare("ZTD")== 0    )  data.insert(ZHD);
    if( it->compare("ZTD")   == 0                                )  data.insert(ZTD);
     
    // SURF - only if explicitely requested
    if( it->compare("MS_TP") == 0 ||   fit == POW_SURF           )  surf.insert(MS_Tp);
    if( it->compare("MS_EP") == 0 ||   fit == POW_SURF           )  surf.insert(MS_Ep);

    if( it->compare("MP_TP") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Tp);
    if( it->compare("MP_MP") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Mp);
    if( it->compare("MP_EP") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Ep);
    if( it->compare("MP_EH") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Eh);
    if( it->compare("MP_WP") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Wp);
    if( it->compare("MP_WH") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Wh);
    if( it->compare("MP_IP") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Ip);
    if( it->compare("MP_IH") == 0 || ( fit == POW && adj == LMM ))  surf.insert(MP_Ih);

    if( it->compare("ML_TH") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Th);
    if( it->compare("ML_MH") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Mh);
    if( it->compare("ML_EP") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Ep);
    if( it->compare("ML_EH") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Eh);
    if( it->compare("ML_WP") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Wp);
    if( it->compare("ML_WH") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Wh);
    if( it->compare("ML_IP") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Ip);
    if( it->compare("ML_IH") == 0 || ( fit == LOG && adj == LMM ))  surf.insert(ML_Ih);

    if( it->compare("SP_TP") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Tp);
    if( it->compare("SP_MP") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Mp);
    if( it->compare("SP_EP") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Ep);
    if( it->compare("SP_EH") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Eh);
    if( it->compare("SP_WP") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Wp);
    if( it->compare("SP_WH") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Wh);
    if( it->compare("SP_IP") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Ip);
    if( it->compare("SP_IH") == 0 || ( fit == POW && adj == LSQ ))  surf.insert(SP_Ih);
     
    if( it->compare("SL_TH") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Th);
    if( it->compare("SL_MH") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Mh);
    if( it->compare("SL_EP") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Ep);
    if( it->compare("SL_EH") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Eh);
    if( it->compare("SL_WP") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Wp);
    if( it->compare("SL_WH") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Wh);
    if( it->compare("SL_IP") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Ip);
    if( it->compare("SL_IH") == 0 || ( fit == LOG && adj == LSQ ))  surf.insert(SL_Ih);

    if( it->compare("ME_PH") == 0 || ( fit == POW && adj == LMM ))  surf.insert(ME_Ph);
    if( it->compare("ME_EH") == 0 || ( fit == POW && adj == LMM ))  surf.insert(ME_Eh);
    if( it->compare("ME_WH") == 0 || ( fit == POW && adj == LMM ))  surf.insert(ME_Wh);

    if( it->compare("EQH")   == 0 || ( fit == POW && adj == LMM ))  surf.insert(eqH);
    if( it->compare("EQW")   == 0 || ( fit == POW && adj == LMM ))  surf.insert(eqW);
    if( it->compare("EQE")   == 0 || ( fit == POW && adj == LMM ))  surf.insert(eqE);

    if( it->compare("SCH")   == 0 || ( fit == POW && adj == LMM ))  surf.insert(scH);
    if( it->compare("SCW")   == 0 || ( fit == POW && adj == LMM ))  surf.insert(scW);
    if( it->compare("SCE")   == 0 || ( fit == POW && adj == LMM ))  surf.insert(scE);

    // RATIO
    if((it->compare("RATX")  == 0 ||
        it->compare("RATIO") == 0 ) && ( fit == POW || fit == POW_SURF ))  data.insert(RATX);
    if((it->compare("RATL")  == 0 ||
        it->compare("RATIO") == 0 ) && ( fit == LOG                    ))  data.insert(RATL);

    if((it->compare("MP_RT") == 0 ||
        it->compare("DR")    == 0 ) && ( fit == POW || fit == POW_SURF ))  surf.insert(MP_Rt);
    if((it->compare("ML_RT") == 0 ||
        it->compare("DR")    == 0 ) && ( fit == LOG                    ))  surf.insert(ML_Rt);

    // old style - on special request only!
    if( it->compare("ER_EP") == 0 )  surf.insert(ER_Ep);
    if( it->compare("LR_EP") == 0 )  surf.insert(LR_Ep);
    if( it->compare("ER_EH") == 0 )  surf.insert(ER_Eh);
    if( it->compare("LR_EH") == 0 )  surf.insert(LR_Eh);
     				                                             
    if( it->compare("ER_WP") == 0 )  surf.insert(ER_Wp);
    if( it->compare("LR_WP") == 0 )  surf.insert(LR_Wp);
    if( it->compare("ER_WH") == 0 )  surf.insert(ER_Wh);
    if( it->compare("LR_WH") == 0 )  surf.insert(LR_Wh);

//  if( it->compare("LR_WE") == 0 )  surf.insert(LR_WE);
 
    // only if FIT OUTPUT requested
    if( ! ofit.empty() ){
      surf.insert(ML_Th);  surf.insert(ML_Mh); 
      surf.insert(MP_Eh);  surf.insert(ML_Ep);  surf.insert(SP_Ep);
      surf.insert(MP_Wh);  surf.insert(ML_Wp);  surf.insert(SP_Wp);
    }

    // only if ASSESSMENT requested
    if( ass.size() == 1 && ass.find("FIT_RATIO") != ass.end() ){}
    else if( ass.size() > 0 ){
      surf.insert(ML_Th);  surf.insert(SL_Th);  surf.insert(SP_Tp);
      surf.insert(ML_Mh);  surf.insert(SL_Mh);  surf.insert(SP_Mp);
    
      surf.insert(MP_Ep);  surf.insert(ML_Eh);
      surf.insert(ML_Ep);  surf.insert(ML_Eh);
      surf.insert(MP_Wp);  surf.insert(ML_Wh);
      surf.insert(ML_Wp);  surf.insert(ML_Wh);
    
      surf.insert(MP_Eh);  surf.insert(MP_Wh);
      surf.insert(ME_Eh);  surf.insert(ME_Wh);  surf.insert(ME_Ph); // fit exponential via altitude
       
      surf.insert(scE);    surf.insert(scW);    surf.insert(scH);   // scale heights
      surf.insert(eqE);    surf.insert(eqW);    surf.insert(eqH);   // scale heights
  
      surf.insert(LR_Ep);  surf.insert(LR_Eh);
      surf.insert(LR_Wp);  surf.insert(LR_Wh);
      surf.insert(ER_Ep);  surf.insert(ER_Eh);
      surf.insert(ER_Wp);  surf.insert(ER_Wh);
	  
      surf.insert(SP_Ep);  surf.insert(SL_Eh);
      surf.insert(SL_Ep);  surf.insert(SL_Eh);
      surf.insert(SP_Wp);  surf.insert(SL_Wh);
      surf.insert(SL_Wp);  surf.insert(SL_Wh);
	    
      surf.insert(LL_Ep);
      surf.insert(LL_Wp);

//    surf.insert(ML_Rt); // NOT NECESSARY
    }
}
  return 0;
}

} // namespace
