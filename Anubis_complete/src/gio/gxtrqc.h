
#ifndef GXTRQC_H
#define GXTRQC_H

#define MISS     "-"     // print missing values
#define XML_NAN  "n/a"   // print missing values
#define STP   (15*60)    // [s] .. 15 min (STEP INTERVAL)
#define GAP   (10*60)    // [s] .. 10 min (GAP DEFINITION)
#define PCS   (30*60)    // [s] .. 30 min (SMALL PIECES)
#define MP_NEPOCHS  (15) // # epochs for multipath
#define MP_LIMIT   (5.0) // sigma-multiplicator for MP cycle-slip & outlier detection
#define MP_UNITS   (100) // multipath [m]->[cm]
#define DAY_EPOCHS (2880) // TEMP


/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements extration class
  Version: $ Rev: $

  2013-01-10 /JD: created

-*/

#include "gio/gxtr.h"
#include "gdata/grnxhdr.h"
#include "gdata/gqcdata.h"
#include "gall/galloqc.h"
#include "gproc/gpreproc.h"
#include "gutils/gsatview.h"

using namespace std;
using namespace pugi;

namespace gnut {

class t_gxtrqc : public t_gxtr
{
 public:
   t_gxtrqc(t_gsetbase* set, const string& pgm, const t_gtime& dt);
   t_gxtrqc();   
  virtual ~t_gxtrqc();

   typedef map<string, pair<double,int> > t_map_mpsat;
   typedef map<t_gtime, t_map_mpsat >     t_map_mpepo;
   typedef map<GOBS, t_map_mpepo >        t_map_mpobs;        // used only for re-ordering
   typedef map<GSYS, t_map_mpobs >        t_map_mpsys;        // used only for re-ordering

   virtual void setDAT(t_galloqc* gobs, t_gallnav* gnav);   
   
   virtual void summary(string site);
   virtual void multipath(const string& prn, GSYS gsys, const t_gtime& epo,
                          const unsigned int& nepo, t_map_mpobs& m_obs);

   virtual void extended_xml(bool b){ _xml_ext = b; };

 protected:
   virtual void _get_settings();
   
   virtual void _list_gap(ostringstream& os);
   virtual void _list_pcs(ostringstream& os);
   virtual void _list_smp(ostringstream& os);

   virtual void _list_band(ostringstream& os);
   virtual void _list_type(ostringstream& os);
   virtual void _list_typh(ostringstream& os);   
   virtual void _list_gsys(ostringstream& os);
// virtual void _list_gobs(ostringstream& os);
   virtual void _list_summ(ostringstream& os);
   virtual void _list_cslp(ostringstream& os);
   virtual void _list_jump(ostringstream& os);
   virtual void _list_prep(ostringstream& os);

   virtual void _header(ostringstream& os);
   virtual void _calcul(ostringstream& os);
   virtual void _observ(ostringstream& os);
   virtual void _nbands(ostringstream& os);
   virtual void _pieces(ostringstream& os);
   virtual void _skyplt(ostringstream& os);
   virtual void _mlpath(ostringstream& os);
   virtual void _snoise(ostringstream& os);      
   virtual void _svinfo(ostringstream& os);      
   virtual void _prepro(ostringstream& os);
   virtual void _grckpi(ostringstream& os){};
   virtual void _summar(ostringstream& os);   
   
   virtual void _xml_meta();
   virtual void _xml_head();
   virtual void _xml_navi();
   virtual void _xml_data();
   virtual void _xml_full();

 protected:
   virtual void _setOut();
   virtual void _clear();
   virtual int  _get_nsat(GSYS gsys);
   virtual int  _get_rnxhdr();
     
   xml_node    _XFMT , _XNAV , _XHDR , _XDAT , _XEXT,  _ROOT ;
// xml_node    _XFMT2, _XNAV2, _XHDR2, _XDAT2, _XEXT2, _ROOT2;

   t_galloqc*  _gobs;
   t_glog*     _xtr;

   t_gtime     _beg, _obs_beg, _obs_sync;
   t_gtime     _end, _obs_end;
   t_gtriple   _xyz_approx();
   t_gtriple   _xyz_est;
   t_gtriple   _xyz_rep;

   t_grec::t_maphdr  _m_rnxHdr;
   
   double      _smp_est;
   double      _min_ele;
   double      _cut_ele;
   double      _hours;
   double      _minInt;
   double      _maxInt;
   double      _pos_cut;
   int         _pos_int;
   bool        _xml;
   bool        _xml_ext;
// bool        _xml2;
   bool        _kinematic;
   bool        _ele_new;
   bool        _ele_app;
   int         _ngap; 
   int         _npcs;

   int         _summ;
   int         _head;
   int         _stat;
   int         _band;
   int         _gaps;
   int         _prep;
   int         _elev;
   int         _mult;
   int         _calc;
   int         _stnr;
   int         _sinf;
   int         _gkpi;

   double      _step;
   int         _tgap;
   int         _tpcs;
   int         _mp_nepochs;
   double      _mp_limit;
   bool        _mp_all;
   bool        _auto_band;
   
   double      _dV_lim;
   double      _dH_lim;
   double      _dG_lim;

   bool        _sat_rec;

   USE_HEALTH  _useHealth;

   t_map_sky_sat  _sat_view;               // satellite visibility (zero horizon)
   t_map_sky_sat  _sat_mask;               // satellite visibility (user cut-off)

   t_galloqc::t_map_stt_smp  _stat_smp;    // interval statistics
   t_galloqc::t_map_stt_sys  _stat_obs;    // observation type statistics GSYS/SAT/GOBS
   t_galloqc::t_map_sky_sys  _stat_ele;    // observation type statistics (elevation bins)
   t_map_mpsys               _m_mp;        // all multipath (sys/obs/epo/sat)

   map<GSYS,t_gtriple>  _m_xyz_est;        // estimated coordinates
   map<GSYS,t_gtriple>  _m_xyz_rep;        // estimated repeatabilities

   map<GSYS, map<GOBS,int> > _m_totAll;    // total interruptions (All below - related to ambiguity set up)
   map<GSYS, map<GOBS,int> > _m_totSlp;    // total interruptions (Cycle slips)
   map<GSYS, map<GOBS,int> > _m_totSat;    // total interruptions (Satellites)
   map<GSYS, map<GOBS,int> > _m_totSig;    // total interruptions (Signals)
   map<GSYS, map<GOBS,int> > _m_totEpo;    // total interruptions (Epochs)
   
   map<GSYS, map<GOBS,int> > _m_satHav;    // total satellite     (Have)
   map<GSYS, map<GOBS,int> > _m_obsHav;    // total observations  (Have)
   map<GSYS, map<GOBS,int> > _m_obsExp;    // total observations  (Expt)
   map<GSYS, map<GOBS,int> > _m_cutExp;    // total observations  (Have), mask cut-off
   map<GSYS, map<GOBS,int> > _m_cutHav;    // total observations  (Expt), mask cut-off

   map<GSYS, map<GOBS,int> > _m_okElev;    // total observations with    elevations
   map<GSYS, map<GOBS,int> > _m_woElev;    // total observations without elevations

   map<GSYS, pair<int,int> > _m_UnusEpo;   // unusable epochs pair<code,phase>
   map<GSYS, pair<int,int> > _m_UnusSat;   // unusable satellite epochs pair<code,phase>
   map<GSYS, int >           _m_FineEpo;   // fine epochs for the processing
   map<GSYS, int >           _nEpoExists;  // existing 
   map<GSYS, int >           _nEpoExpect;  // expected epochs for system

   map<t_gtime, int>         _m_Breaks;    // receiver clk jumps (pha/cod inconsistency)
   int                       _ClkSync;     // receiver clk jumps (clock synchronization)

   map<t_gtime, map<GSYS,int> >                   _m_nSatEpo;  // satellites in epoch per system
   map<t_gtime, map<string, map<GOBS, double> > > _m_Slips;    // cycle slips map
   map<t_gtime, map<string, map<GOBS, int> > >    _m_SlipsGap; // cycle slips map due to gap - 1(epoch gap), 2(sat gap), 3(GOBS gap)

   string                _vers;
   shared_ptr<t_gqcdata> _qcdata;
};

} // namespace

#endif
