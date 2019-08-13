
#ifndef GSPP_H
#define GSPP_H

#include "gutils/gtime.h"
#include "gutils/gmutex.h"
#include "gprod/gprodcrd.h"
#include "gprod/gprodtrp.h"
#include "gset/gsetbase.h"
#include "gset/gsetout.h"
#include "gset/gsetgen.h"
#include "gset/gsetgnss.h"
#include "gset/gsetrec.h"
#include "gset/gsetproc.h"
#include "gio/glog.h"
#include "gio/giof.h"
#include "gall/gallobs.h"
#include "gall/gallnav.h"
#include "gall/gallprod.h"
#include "gall/gallobj.h"
#include "gall/gallbias.h"
#include "gmodels/gbancroft.h"
#include "gmodels/gmodel.h"
#include "gmodels/gpar.h"

namespace gnut {
   
class t_gspp
{
 public:
   t_gspp(string mark, t_gsetbase* set);
   virtual ~t_gspp();

   virtual int  processBatch(const t_gtime& beg, const t_gtime& end) = 0;
   
   virtual void setDAT(t_gallobs*   gobs, t_gallnav* gnav);   
           void setOUT(t_gallprod* products);
   virtual void setOBJ(t_gallobj*   gobj);
   virtual void setDCB(t_gallbias*  gbias);
   virtual void setFCB(t_gallbias*  gbias);
   /** @brief Set general log file. */
   virtual void glog(t_glog* glog) {_glog = glog;};
   
           void tropo(bool tropo);
           void tropo_slant(bool slant);
           void phase(bool phase);
           void setgnss(GSYS sys);   
   string  site(){ return _site; }   
     void  fixObsTypes();
   
 protected:
   shared_ptr<t_gobj> _grec;        /**< Transmitter/receiver object. */
   t_gallobj*    _gallobj;          /**< Objects transmitter/receiver. */
   t_gallbias*   _gallbias;         /**< Differential code biases. */
   t_gallbias*	 _gallfcb;			/**< Phase fractional cycle bias. */
   // constraining parameters
   CONSTRPAR     _crd_est;
   t_gtime       _crd_begStat;
   t_gtime       _crd_endStat;
   t_gtime       _ztd_begStat;
   t_gtime       _ztd_endStat;
   
   // a priory parameters
   double        _aprox_ztd_xml;    /**< Approximate ztd. */
   
   // observation weighting
   OBSWEIGHT _weight;
   
   // observation type model
   OBSCOMBIN _observ;
   
   /** @brief Get settings from XML file and set local variables. */
   virtual int   _get_settings();   
   
   /** @brief Set output file (empty). */
   virtual void  _setOut();
   
   /** @brief Set processing log file (<ppp>). */   
   virtual void  _setLog();
   
   
   bool          _valid_crd_xml;
   bool          _valid_ztd_xml;
   
   t_gmodel*     _gModel;           /**< models. */  
   string        _site;             /**< Site internal ID. */
   t_gsetbase*   _set;              /**< Base setting. */
   t_glog*       _log;              /**< Processing log output. */
   t_glog*       _glog;             /**< Genereal log output. */   
   t_giof*       _res;              /**< Residuals output. */
   t_gallobs*    _gobs;             /**< Observation data. */
   t_gallnav*    _gnav;             /**< Objects for ephemerides. */
   t_gallprod*   _gmet;	            /**< Troposphere products input*/
   t_gallprod*   _gion;	            /**< Ionosphere products input*/
   t_gallprod*   _allprod;          /**< Products output. */
   bool          _phase;            /**< Phase is used. */
   bool          _tropo_est;        /**< Tropo is estimated. */
   bool          _iono_est;         /**< Iono is estimated. */   
   bool          _tropo_grad;       /**< Tropo horizontal gradients. */
   bool          _tropo_slant;      /**< Tropo slants are produced. */
   GSYS          _gnss;             /**< GNSS system to be used. */
   ZTDMPFUNC     _ztd_mf;           /**< ZTD mapping function. */
   GRDMPFUNC     _grd_mf;           /**< GRD mapping function. */   
   double        _minElev;          /**< Elevation cut-off. */
   double        _sampling;         /**< Sampling interval. */
   double        _scale;            /**< Sampling scaling factor. */
   bool          _initialized;      /**< Initialized status. */
   //set<string>   _sats;             /**< Configured satellites to be used. */
   int           _nSat;             /**< Number of satellites comming into the processing */
   int           _nSat_excl;        /**< Number of excluded satellites due to various reason */   

   // init sigma
   double        _sig_init_crd;     /**< Initial coordinates sigma. */
   double        _sig_init_ztd;     /**< Initial ZTD sigma. */
   double        _sig_init_vion;    /**< Initial VION sigma. */
   double        _sig_init_grd;     /**< Initial GRD  sigma. */   
   
   // ISB sigma
   double        _sig_init_glo;     /**< Initial GLONASS ISB sigma. */
   double        _sig_init_gal;     /**< Initial Galileo ISB sigma. */
   double        _sig_init_bds;     /**< Initial BeiDou ISB sigma. */
   double        _sig_init_qzs;     /**< Initial QZSS ISB sigma. */
   
   // observations sigma

   double        _sigCodeGPS;       /**< Code sigma.  - GPS*/
   double        _sigPhaseGPS;      /**< Phase sigma. - GPS*/
   double        _sigCodeGLO;       /**< Code sigma.  - GLONASS*/
   double        _sigPhaseGLO;      /**< Phase sigma. - GLONASS*/
   double        _sigCodeGAL;       /**< Code sigma.  - Galileo*/
   double        _sigPhaseGAL;      /**< Phase sigma. - Galileo*/   
   double        _sigCodeBDS;       /**< Code sigma.  - BeiDou*/
   double        _sigPhaseBDS;      /**< Phase sigma. - BeiDou*/
   double        _sigCodeQZS;       /**< Code sigma.  - QZSS*/
   double        _sigPhaseQZS;      /**< Phase sigma. - QZSS*/   

   bool          _pos_kin;
   bool          _extern_log;
   
   bool          _use_ecl;
   bool          _success;

#ifdef BMUTEX
   boost::mutex  _mutex;
#endif
   t_gmutex      _gmutex;   

// signals used for processing
   map<string, map<GOBSBAND, GOBSATTR>>  _signals;

};

} // namespace

#endif

