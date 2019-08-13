
#ifndef GXTRGRC_H
#define GXTRGRC_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements extration class
  Version: $ Rev: $

  2013-01-10 /JD: created

-*/

#include "gio/gxtrqc.h"

using namespace std;
using namespace pugi;

namespace gnut {

class t_gxtrgrc : public t_gxtrqc
{
 public:
   t_gxtrgrc(t_gsetbase* set, const string& pgm, const t_gtime& dt);
   t_gxtrgrc();   
  
   virtual ~t_gxtrgrc();

 protected:
   virtual void _get_settings();
   
   virtual void _grckpi(ostringstream& os);      

};

} // namespace

#endif
