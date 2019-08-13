
#ifndef RINEXN_H
#define RINEXN_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2011-04-20 /JD: created

-*/

#include <vector> 
#include <sstream>

#include "gcoders/gcoder.h"
#include "gset/gsetinp.h"

using namespace std;

namespace gnut {

class t_rinexn : public t_gcoder {

 public:
   t_rinexn( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
  ~t_rinexn(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);
   
  virtual  int encode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int encode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

  void gnsssys(char s ){ _gnsssys = s; }
  char gnsssys(){ return _gnsssys; }

 protected:
//  virtual int  _gset(t_gsetbase* set);
          bool _filter_gnav(shared_ptr<t_gnav> geph, const string& prn);
   

 private:
// map<GSYS,set<GNAVTYPE>> _gnavtype;
   char     _gnsssys;
// bool     _check;
   t_gtime  _check_dt; // to validate the message
};

} // namespace

#endif
