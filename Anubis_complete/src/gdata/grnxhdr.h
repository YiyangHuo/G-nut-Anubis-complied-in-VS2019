#ifndef T_RNXHDR
#define T_RNXHDR

#include <string>
#include <vector>

#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gdata/gobsgnss.h"

namespace gnut {

class t_rnxhdr
{
 public:
   t_rnxhdr();
  ~t_rnxhdr();

   typedef vector< pair<GOBS,double> > t_vobstypes;
   typedef map<string, t_vobstypes>     t_obstypes;
   
   string                     path() const {return _path;}
   void                       path(string s) {_path = s;}
   
   char                       rnxsys() {return _rnxsys;}
   void                       rnxsys(char s) {_rnxsys = s;}
   
   string                     rnxver() const {return _rnxver;}
   void                       rnxver(string s) {_rnxver = s;}
   
   string                     program() const {return _program;}
   void                       program(string pgm) {_program = pgm;}
   
   string                     runby() const {return _runby;}
   void                       runby(string rnb) {_runby = rnb;}
   
   t_gtime                    gtime() const {return _gtime;}
   void                       gtime(t_gtime& t) {_gtime = t;}
   
   vector<string>             comment() {return _comment;}
   void                       comment(vector<string>& cmt) {_comment = cmt;}
   
   string                     markname() const {return _markname;}
   void                       markname(string mrk) {_markname = mrk;}
   
   string                     marknumb() const{return _marknumb;}
   void                       marknumb(string mnb) {_marknumb = mnb;}
   
   string                     marktype() const {return _marktype;}
   void                       marktype(string mtp) {_marktype = mtp;}
   
   string                     observer() const {return _observer;}
   void                       observer(string obsrv) {_observer = obsrv;}
   
   string                     agency() const {return _agency;}
   void                       agency(string agn) {_agency = agn;}
   
   string                     recnumb() const {return _recnumb;}
   void                       recnumb(string rnb) {_recnumb = rnb;}
   
   string                     rectype() const {return _rectype;}
   void                       rectype(string rtp) {_rectype = rtp;}
   
   string                     recvers() const {return _recvers;}
   void                       recvers(string rvs) {_recvers = rvs;}
   
   string                     antnumb() const {return _antnumb;}
   void                       antnumb(string ant) {_antnumb = ant;}
   
   string                     anttype() const {return _anttype;}
   void                       anttype(string ant) {_anttype = ant;}
   
   t_gtriple                  aprxyz() const {return _aprxyz;}
   void                       aprxyz(t_gtriple& apr) {_aprxyz = apr;}
   
   t_gtriple                  antxyz() const {return _antxyz;}
   void                       antxyz(t_gtriple& ecc) {_antxyz = ecc;}
   
   t_gtriple                  antneu() const {return _antneu;}
   void                       antneu(t_gtriple& ecc) {_antneu = ecc;}

   string                     strength() const {return _strength;}
   void                       strength(string s) {_strength = s;}
   
   double                     interval() const {return _interval;}
   void                       interval(double i) {_interval = i;}
   
   t_gtime                    first() const {return _first;}
   void                       first(t_gtime& frst) {_first = frst;}
   
   t_gtime                    last() const {return _last;}
   void                       last(t_gtime& lst) {_last = lst;}
   
   int                        leapsec() const {return _leapsec;}
   void                       leapsec(int ls) {_leapsec = ls;}
   
   int                        numsats() const {return _numsats;}
   void                       numsats(int nums) {_numsats = nums;}
   
   t_obstypes                 mapobs() const {return _mapobs;}
   void                       mapobs(t_obstypes types) {_mapobs = types;}
   
   t_obstypes                 mapcyc() const {return _mapcyc;}
   void                       mapcyc(t_obstypes types) {_mapcyc = types;}

   t_obstypes                 glofrq() const {return _glofrq;}
   void                       glofrq(t_obstypes types) {_glofrq = types;}
   
   t_vobstypes                globia() const {return _globia;}
   void                       globia(t_vobstypes types) {_globia = types;}

   void clear();

   friend ostream& operator<<(ostream& os, const t_rnxhdr& x);
   
 private:   
   char                       _rnxsys;     // G=GPS, R=GLO, E=GAL, S=SBAS, M=Mix
   string                     _path;       // rinex file path
   string                     _rnxver;     // rinex version
   string                     _program;    // name of program creating RINEX file
   string                     _runby;      // name of agency  creating RINEX file
   t_gtime                    _gtime;      // name of date and file of RINEX creation
   vector<string>             _comment;    // vector of comments
   string                     _markname;   // marker name
   string                     _marknumb;   // marker number
   string                     _marktype;   // marker type
   string                     _observer;   // name of observer
   string                     _agency;     // name of agency
   string                     _recnumb;    // receiver number
   string                     _rectype;    // receiver type
   string                     _recvers;    // receiver version
   string                     _antnumb;    // antenna number
   string                     _anttype;    // antenna type
   t_gtriple                  _aprxyz;     // approximate xyz position [m]
   t_gtriple                  _antxyz;     // antenna xyx eccentricities [m]
   t_gtriple                  _antneu;     // antenna north/east/up eccentricities [m]                
   string                     _strength;   // signal strength [DBHZ/...]
   double                     _interval;   // interval [sec]
   t_gtime                    _first;      // time of first observation
   t_gtime                    _last;       // time of last observation
   int                        _leapsec;    // leapseconds since 6-Jan-1980
   int                        _numsats;    // number of satellites
   t_obstypes                 _mapobs;     // map of GOBS and scaling factors
   t_obstypes                 _mapcyc;     // map of GOBS phase quater-cycle shifts
   t_obstypes                 _glofrq;     // map of GLONASS slot/frequency
   t_vobstypes                _globia;     // vec of GLONASS obs code-phase biases
   
};

} // namespace

#endif
