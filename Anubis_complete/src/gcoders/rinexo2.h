
#ifndef RINEXO2_H
#define RINEXO2_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: Clock RINEX encoder/decoder
  Version: $Rev:$

  2011-11-04 /JD: created
  2012-09-24 /JD: RINEX 2.11 implemented
  2012-09-27 /JD: support of multi-GNSS (phase not converted to [m])
  2013-03-08 /JD: filtering via general settings
  2014-05-03 /PV: gobj getting information from rinexo header
  2014-11-19 /JD: shared pointers
  2014-12-21 /JD: advanced request filtering
 
  Todo: header information completed 
        encoders implementation

-*/

#include <string> 
#include <vector> 

#include "gall/gallobs.h"
#include "gall/gallobj.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gutils/gsys.h"
#include "gutils/gobs.h"
#include "gutils/gtypeconv.h"
#include "gcoders/gcoder.h"
#include "gdata/grnxhdr.h"

using namespace std;

namespace gnut {

class t_rinexo2 : public t_gcoder {

 public:
   t_rinexo2( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
  ~t_rinexo2();

  virtual void clear();

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg) = 0;
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg) = 0;
  
  virtual int encode_head( char* buff, int sz,           vector<string>& errmsg ) = 0;
  virtual int encode_data( char* buff, int sz, int& cnt, vector<string>& errmsg ) = 0;

 protected:

  virtual int _decode_head();
  virtual int _decode_data();

  virtual int _read_epoch();                                               // read epoch & number of satellites, return flag   
  virtual int _read_syskey(string& sat, string& key);                      // read sat/sys-specific key to _mapobs (G,E,R,.. M, satellite)
  virtual int _read_satvec(vector<string>& satellites);                    // read satellite list (vector)
  virtual int _read_obstypes(const string& sat, const string& sys);        // read single satellite observation types
  virtual int _read_obs(const unsigned int& idx, 
		        const t_rnxhdr::t_vobstypes::const_iterator& it,
		              t_spt_gobs obs);                             // fill single observation element
  virtual int _fill_head();                                                // fill header information
  virtual int _fill_data();                                                // fill observation data structure
  virtual int _check_head();                                               // fill header information
  virtual int _fix_band(string sys, string& go);                           // check band (extended RINEX 2.11)
  virtual int _null_log(const string& sat, const string& obstype);         // write log for null observations (sat, obstype)
  virtual int _stop_read();                                                // common stop reading used in local subroutines

  string                     _csys;       // A1 (G=GPS, R=GLO, E=GAL, S=SBAS, M=Mix), but may be also G11, ..
  t_gtime                    _epoch;      // working epoch
  string                     _line;       // working line read from
  string                     _site;       // working site
  int                        _tmpsize;    // working amount of bytes processed
  int                        _consume;    // working total amount of bytes (return)
  bool                       _complete;   // working flag for completed epoch decoding
  char                       _flag;       // working special event flag
  vector<t_spt_gobs>         _vobs;       // working GNSS observation vector
  int                        _nsat;       // working # of satellites
  int                        _count;      // working # of epochs

  t_gtime                    _epo_beg;
  t_gtime                    _epo_end;

  int                        _xbeg;       // # epochs filtered before BEG
  int                        _xend;       // # epochs filtered after  END
  int                        _xsmp;       // # epochs filtered by sampling
  int                        _xsys;       // # obs data filtered by systems

  t_rnxhdr                   _rnxhdr;     // Rinex header instance
  vector<string>             _comment;    // Rinex comments
  vector<char>               _pcosat;     // GNSS (G/R/E/S/...) for antenna phase center vs. ARP
  vector<string>             _pcosys;     // XYZ or NEU system for pco
  vector<t_gtriple>          _pcoecc;     // PCO-ARP values in system above (XYZ, NEU) [m]
  t_rnxhdr::t_obstypes       _mapobs;     // GNSS/sat (G/R/E/S/...) observation group
                                          // vector due to saving order in header !
                                          // including scale factor (1000,100,10,1)
};

} // namespace

#endif
