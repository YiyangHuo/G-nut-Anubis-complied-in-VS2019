
#ifndef STOCHASTIC_H
#define STOCHASTIC_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: stochastic models
  Version: $ Rev: $

  2012-05-02 /PV: created

-*/

#include "gutils/gtime.h"

using namespace std;

namespace gnut {   

class t_stochastic 
{
 public:
  t_stochastic();
  virtual ~t_stochastic() {};
  
  virtual double getQ()
  {
    return 0.0;
  }
  

 protected:
  
 private:

};

class t_randomwalk : public t_stochastic 
{
 public:
  t_randomwalk();
  virtual ~t_randomwalk() {};
   
  virtual double getQ();
  void setTprev(const t_gtime&);
  void setTcurr(const t_gtime&);
  void updateTime(const t_gtime&);
  void setq(double q);
  double get_dt();
  
 protected:
 private:
  t_gtime _Tprev;
  t_gtime _Tcurr;
  double _dSig;
};

class t_whitenoise : public t_stochastic
{
 public:
  t_whitenoise(double);
  virtual ~t_whitenoise() {};
   
   virtual double getQ();
   void setVar(double);
   
 private:
   double _var;
};

} // namespace

#endif // STOCHASTIC_H
