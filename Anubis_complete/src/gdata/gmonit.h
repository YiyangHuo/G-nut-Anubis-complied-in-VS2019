
#ifndef GMONIT_H
#define GMONIT_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $Rev:$

  2011-03-25 /JD: created

-*/

#include <string>
#include <sstream>

using namespace std;

namespace gnut {

class t_gmonit {

 public:
  t_gmonit( string id );
  virtual ~t_gmonit();

  virtual void  show( ostringstream& os, int verb );
//  virtual string show( int verb, int repeat, string& buff ) = 0;

 protected:
  string   _moni_id;

 private:

};

} // namespace

#endif
