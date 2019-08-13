
#ifndef GNWMBASE_H
#define GNWMBASE_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implementation of NWM base data class
  Version: $ Rev: $

  2013-06-25 /JD: created

-*/

#include <string>
#include <set>

#include "gdata/gdata.h"
#include "gset/gsetbase.h"

#define NWM_UNKNOWN     -999
#define NWM_DEF_T_LAPSE  6.5
#define NWM_DEF_I_DECAY  2.5
#define NWM_DEF_W_DECAY  2.5
#define NWM_DEF_E_DECAY  2.5

using namespace std;

namespace gnut {   

enum P2GEOP   { HYPS, PSAS, CUMUL };
enum LEVELS   { HIGH, FULL, HALF  };

enum NWMFIT   { POW, LOG, POW_SURF };
enum NWMADJ   { LMM, LSQ };
enum PROFTEMP { dT_ZERO, dT_LINEAR, dT_CONST };
enum PROFZWD  { dW_DE, dW_DE_MOPS, dW_MOPS };
enum SURFZWD  { W_INTEGR, W_MOPS, W_DE, W_DE_AUTO, W_DE_GAMMA, W_DE_LAMBDA, W_ESA, W_ESA_L };

// ZTD must be after ZHD, ZWD!!!
enum MET_DATA { GEOM, GEOP, GRAV, PRES, EPRE, TEMP, SHUM, ZHD, ZWD, ZTD, IWV, TM, NGRD, EGRD, NGH, NGW, EGH, EGW, RATX, RATL, RATE, RATW };
enum MET_SURF { LAT, LON, HEL, XPRJ, YPRJ, GEOID, OROGR, T2M, TSK, TCWV, LSM, MSL,
     
                // adjusted surface values
                TSURF, // as SURF(TEMP) = adjusted value    vs  DATA(TEMP) = extrapolated value !
		ESURF, // as SURF(EPRE) = adjusted value    vs  DATA(EPRE) = extrapolated value !
		QSURF, // as SURF(SHUM) = adjusted value    vs  DATA(SHUM) = extrapolated value !
 
                // model caculation
		MSURF, // as SURF(TM)   = model value       vs  DATA(TM)   = integrated value !
		WSURF, // as SURF(ZWD)  = model value       vs  DATA(ZWD)  = integrated value !

                // adopted (default) values for lapse/decay rate (see settings: vert_fit + vert_adj)
                dT, // temperature lapse rate  [K/km]
                dM, // mean temper lapse rate  [K/km]
                dE, // WV  exponential decay rate [-]
                dW, // ZWD exponential decay rate [-]
                dI, // IWV exponential decay rate [-]
                dR, // ZWD/WV ratio decay rate [%/km]

                eqH, // ZHD equiv height [km]  ~43km
                eqW, // ZWD equiv height [km]  ~12km
                eqE, // WV  equiv height [km]  ~12km
     
                scH, // ZHD scale height [km]  ~8km
                scW, // ZWD scale height [km]  ~2km
                scE, // WV  scale height [km]  ~2km
      
                // M.[LMM], S.[LSQ], Ax[AVE]  ====  .P[POW], .L[LOG], .E[EXP]
    	        // --------------------------------------------------
                // pow-fit             (+surf)    log-fit
	        // --------------------------------------------------
                MP_Tp, SP_Tp, AP_Tp, MS_Tp,    ML_Th, SL_Th, AL_Th, //      temperature (pow-fit, log-fit)
                MP_Mp, SP_Mp, AP_Mp,           ML_Mh, SL_Mh, AL_Mh, // mean temperature (pow-fit, log-fit)

                MP_Ep, SP_Ep, AP_Ep, MS_Ep,    ML_Ep, SL_Ep, AL_Ep, // water vapour pressure via E
                MP_Eh, SP_Eh, AP_Eh,           ML_Eh, SL_Eh, AL_Eh, // water vapour pressure via E
                MP_Fp, SP_Fp, AP_Fp,           ML_Fp, SL_Fp, AL_Fp, // water vapour pressure via ZWD
                MP_Fh, SP_Fh, AP_Fh,           ML_Fh, SL_Fh, AL_Fh, // water vapour pressure via ZWD

                MP_Wp, SP_Wp, AP_Wp,           ML_Wp, SL_Wp, AL_Wp, // zenith wet delay via ZWD
                MP_Wh, SP_Wh, AP_Wh,           ML_Wh, SL_Wh, AL_Wh, // zenith wet delay via ZWD

                MP_Ip, SP_Ip, AP_Ip,           ML_Ip, SL_Ip, AL_Ip, // integrated water vapour via IWV
                MP_Ih, SP_Ih, AP_Ih,           ML_Ih, SL_Ih, AL_Ih, // integrated water vapour via IWV

                MP_WE,                         ML_WE,               // GPTw
                
                ME_Ph, ME_Wh, ME_Eh,                                // P/ZHD, ZWD, WV decays

                MP_Wi, SP_Wi,                  ML_Wi, SL_Wi,
                MP_Fi, SP_Fi,                  ML_Fi, SL_Fi,
     
                MP_Rt, SP_Rt,                  ML_Rt, SL_Rt,        // ratio of decay rates
    
// TESTING/OBSOLETE/LEGACY
//
                LM_Ep, LM_Eh,   LM_Fp, LM_Fh,
                LR_Ep, LR_Eh,   LR_Fp, LR_Fh,
                ER_Ep, ER_Eh,   ER_Fp, ER_Fh,
                LL_Ep, LL_Eh,   LL_Fp, LL_Fh,

                LM_Wp, LM_Wh,
                LR_Wp, LR_Wh,
                ER_Wp, ER_Wh,
                LL_Wp, LL_Wh,

                LR_WE,
                ER_Wi,
                ER_Fi,

                LR_Rt, ER_Rt
           };


// ----------
class t_gnwmbase : public t_gdata {

 public:
   t_gnwmbase();
   t_gnwmbase(set<MET_DATA> data, set<MET_SURF> surf);
   virtual ~t_gnwmbase();

   static int set_param(t_gsetbase* gset, set<MET_SURF>& surf, set<MET_DATA>& data);

   static string met_surf_str(const MET_SURF& type, bool full = true );
   static string met_data_str(const MET_DATA& type, bool full = true );

   static string met_surf_id(const MET_SURF& type );
   static string met_data_id(const MET_DATA& type );
   
   static int met_data(string type );
   static int met_surf(string type );

   static int digit_dec(const MET_SURF& type );
   static int digit_dec(const MET_DATA& type );
   static int digit_len(const MET_SURF& type );
   static int digit_len(const MET_DATA& type );

};

} // namespace

#endif

