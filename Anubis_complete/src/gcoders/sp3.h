
#ifndef SP3_H
#define SP3_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2011-04-20 /JD: created

-*/

#include <vector> 

#include "gcoders/gcoder.h"
#include "gutils/gtime.h"

using namespace std;

namespace gnut {

class t_sp3 : public t_gcoder {

 public:
   t_sp3( t_gsetbase* s, string version = "", int sz = DEFAULT_BUFFER_SIZE );
  ~t_sp3(){};

  virtual  int decode_head(char* buff, int sz,           vector<string>& errmsg);
  virtual  int decode_data(char* buff, int sz, int& cnt, vector<string>& errmsg);

 protected:

 private:
  t_gtime        _start;
  t_gtime        _lastepo;
  long           _orbintv; // [sec]
  int            _nepochs;  
  int            _nrecmax;
  int            _nrecord;
  string         _orbrefs;
  string         _orbtype;
  vector<string> _prn;    
  vector<int>    _acc;    
  vector<string> _timesys;
  vector<int>    _accbase;
};

} // namespace

#endif
