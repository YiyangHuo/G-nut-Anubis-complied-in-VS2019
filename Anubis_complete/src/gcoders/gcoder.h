
#ifndef GCODER_H
#define GCODER_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2011-04-20 /JD: created
  2013-03-07 /JD: general settings supported

-*/

#include <set>
#include <map>
#include <vector>
#include <string>
#include <cstring>
#include <memory>
#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gio/gio.h"
#include "gio/gnote.h"
#include "gdata/gdata.h"
#include "gutils/gsys.h"
#include "gutils/gmutex.h"
#include "gutils/gcommon.h"
#include "gset/gsetbase.h"
#include "gset/gsetgnss.h"
#include "gset/gsetgen.h"
#include "gset/gsetvar.h"

#define BUFFER_INCREASE_FAC  1.5
#define DEFAULT_BUFFER_SIZE  4096
#define MAXIMUM_BUFFER_SIZE  240000 // 9000000
using	 namespace std;
namespace gnut {

// enum t_irc { failure = -1, success, fatal } ; // return code

using namespace std;

class t_gio;
class t_gcoder {

 public:
   t_gcoder(t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE, string id = "gcoder" );
   virtual ~t_gcoder();

   typedef map<int, pair<GOBS,int> > t_vobstypes;
   virtual void clear();

   virtual void   version(string s){ _version = s; }
   virtual string version(){ return _version; }

   virtual void    out_length(int len){ _out_len = len; }
   virtual int     out_length(){ return _out_len; }

   virtual void    out_sample(float smp){ _out_smp = smp; }
   virtual float   out_sample(){ return _out_smp; }

   virtual void    out_epoch(t_gtime epo){ _out_epo = epo; }
   virtual t_gtime out_epoch(){     return _out_epo; }

   virtual int decode_head( char* buff, int sz,           vector<string>& errmsg ){ return 0; } // = 0;
   virtual int decode_data( char* buff, int sz, int& cnt, vector<string>& errmsg ){ return 0; } // = 0;

   virtual int encode_head( char* buff, int sz,           vector<string>& errmsg ){ return 0; } // = 0;
   virtual int encode_data( char* buff, int sz, int& cnt, vector<string>& errmsg ){ return 0; } // = 0;

   int irc(){ return _irc; }
   
   void glog(t_glog* l){ _log = l; }             // set/get glog pointer
   t_glog* glog(){ return _log; }

   void path(string s);
   string path(){ return _fname; }

   void add_gio(weak_ptr<t_gio> p){ _gio_ptr = p; }

   int add_data( string data_id, t_gdata* data );
// int rem_data( string data_id );
   t_gdata* data( string data_id ){ return _data[data_id]; }

   void            mesg(t_note m, string s);     // set/get notes (messages/warning/errors)
   vector<t_gnote> mesg();

   void close_with_warning(bool b){ _close_with_warning = b; }

   void pgm(string pgm) {_pgm = pgm;}
   
 protected:
   int   endpos(){ return _endpos; }
   int     size(){ return _endpos; }             // number of char elements in buffer (excl \0)
   char* buffer(){ return _buffer; }

   virtual int  _gset(t_gsetbase* s);            // get settings (NO MORE PUBLIC, only via CONSTRUCT)
   virtual void _init();
   virtual void _add_data(string id,t_gdata* data){};
   virtual bool _filter_epoch(const t_gtime& epo);
   virtual bool _filter_gnss(const string& prn);

   int _getline(string& str, int from_pos = 0);
   int _add2buffer( char* buff, int sz );
   int _consume( int bytes_to_eat );

   weak_ptr<t_gio>         _gio_ptr;
// vector<string>          _err;                 // cummative error message default systems
   vector<t_gnote>         _notes;               // cummative notes message/warning/error
   string                  _fname;               // decoded file
   string                  _class;               // string for reporting
   bool                    _initialized;         // if initialized
   bool                    _gnss;                // if gnss definition is requested
   int                     _out_len;             // [min] encoder data batch length
   float                   _out_smp;             // [sec] encoder data batch sample
   t_gtime                 _out_epo;             // encoder refrence epoch
   string                  _version;             // format version
   string                  _initver;             // format initial version
   map<string,t_gdata*>    _data;                // data pointer
   t_glog*                 _log;                 // log pointer
   t_gsetbase*             _set;                 // set pointer
   char*                   _buffer;              // class buffer
   int                     _buffsz;              // size of buffer
   int                     _endpos;              // after last position
   int                     _irc;                 // IRC code   
   bool                    _close_with_warning;  // close with warnings (desctructor)
   
   // settings
   t_gtime     _beg;                             // default beg time
   t_gtime     _end;                             // default end time
   double      _int;                             // default interval
   int         _scl;                             // default scaling for decimation (if >= 1Hz)
   set<string> _sys;                             // default systems
   set<string> _rec;                             // default sites/receivers
   
   map<GSYS,set<string>> _sat;                   // default satellites
   map<GSYS,set<string>> _obs;                   // default observations (all mixed for single GNSS !!)
   map<GSYS,set<string>> _nav;                   // default navigation messages

  // ENCODING
   int _fill_buffer( char* buff, int sz );
   stringstream     _ss;
   long             _ss_position;
   string           _pgm;
   bool             _hdr;
   
#ifdef BMUTEX
   boost::mutex            _mutex;               // buffer mutex
#else
  t_gmutex                 _mutex;               // buffer mutex
#endif
 private:

};

} // namespace

#endif
