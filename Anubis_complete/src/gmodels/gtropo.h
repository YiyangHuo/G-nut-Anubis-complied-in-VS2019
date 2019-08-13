
#ifndef GTROPO_H
#define GTROPO_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements troposphere model class
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gall/gallsurf.h"
#include "gall/gallprod.h"
#include "gprod/gprodmet.h"
#include "gprod/gprodtrp.h"
#include "gprod/gprodcrd.h"
#include "gmodels/gnwmbase.h"
#include "gmodels/gnwmsurf.h"
#include "gmodels/ginterp.h"
#include "gmodels/ggpt.h"

using namespace std;

namespace gnut {   

class t_gtropo {
 public:
//t_gtropo(string site);
  t_gtropo();
  virtual ~t_gtropo();
  
  void nwm(t_gallsurf* nwm);                          // set pointer to gsurf data
  void met(t_gallprod* met);                          // set pointer to gproduct data

  virtual double getZHD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
  virtual double getZWD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
   
  virtual bool   init(string site = "");
  virtual string site()const{return _site;}

 protected:
  double _interp_temporal(string site, t_gtime epo, METEO_ID type);

  string       _site;
  t_gtriple    _ell;
  t_gpt        _gpt;
  t_gallsurf*  _nwm;
  t_gallprod*  _met;
  t_ginterp    _interp;

  map<string, t_gtriple>                              _m_sites;
  map<string, map< t_gtime, shared_ptr<t_gprodmet> >> _m_prods;
};

// ----------------------------- SITE MODELS ------------------------------------
 
class t_saast : public t_gtropo {
 public:
   t_saast(){}
  ~t_saast(){}

   virtual double getSTD(const double& ele,    const double& hel);   // ! Radians: elevation
   virtual double getZHD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& ell, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
};


class t_davis : public t_gtropo {
 public:
   t_davis(){}
  ~t_davis(){}
   
   virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]  
};


class t_hopf : public t_gtropo {
 public:
   t_hopf(){}
  ~t_hopf(){}
   
   virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
};


class t_baby : public t_gtropo {
 public:
   t_baby(){}
  ~t_baby(){}
   
   virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
};


class t_chao : public t_gtropo {
 public:
   t_chao(){}
  ~t_chao(){}
   
   virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
};


class t_ifad : public t_gtropo {
 public:
   t_ifad(){}
  ~t_ifad(){}

   virtual double getZHD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& ele, const t_gtime& epo); // ! Radians: Ell[0] and Ell[1]
};


class t_askn : public t_gtropo {
 public:
   t_askn(){}
  ~t_askn(){}
   
   virtual double getZHD(const t_gtriple& Ele, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& Ele, const t_gtime& epo);  // ! Radians: Ell[0] and Ell[1]
   virtual double getZWD(const t_gtriple& Ell,
		         double  T,          // temperature [K]           
		         double  E,	     // water vapour [Pa]         
     		         double  dT,	     // T lapse rate [K/km]     
		         double  dE,	     // E lapse rate [-]        
		         double  dW,	     // ZWD lapse rate [-]        
		         double  ref,        // ZWD reference valeu [m]
                         double& ratio,      // combi ratio of E/ZWD lps rares [%]
			 double  tm=NWM_UNKNOWN
		        );
};

} // namespace

#endif
