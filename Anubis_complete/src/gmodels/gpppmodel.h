
#ifndef GPPPMODEL_H
#define GPPPMODEL_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: various PPP models
  Version: $ Rev: $

  2012-09-20 /PV: created

-*/

#include <string>
#include <map>
#include <cmath>

#include "gmodels/gsppmodel.h"
#include "gutils/gtypeconv.h"
#include "gmodels/gephplan.h"
#include "gall/gallobj.h"

using namespace std;

namespace gnut {   

class t_gpppmodel : public t_gsppmodel
{
 public:
   t_gpppmodel(string site,t_glog* glog, t_gsetbase* settings);
   t_gpppmodel() {};
   virtual ~t_gpppmodel();
   
   virtual double windUp(t_gsatdata& satdata, const ColumnVector&);   
   virtual double cmpObs(t_gtime& epoch, t_gallpar& param, t_gsatdata&, t_gobs& gobs, bool com);
   virtual double tropoDelay(t_gtime& epoch, t_gallpar& param, t_gtriple ell, t_gsatdata& satdata);


   virtual void setOBJ(t_gallobj* obj);
   
   // attitude modeling - public interface
   int  attitude_old(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k); // from RTKlib (to remove)
   int  attitude(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);  
    
  protected:   
    
    // From RTKlib - needs to be removed
    int  _yaw(t_gsatdata& satdata,string antype, ColumnVector& xs, ColumnVector& ys, ColumnVector& zs);  
    
   // attitude niminal modeling
   void _ysm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);
   void _onm(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);

    void _noon_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    void _midnight_turn(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    
    // attitude for GPS Block IIA
    void _attitude_GPSIIA(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    void _midnight_turn_GPSIIA(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    
    // attitude for GPS Block IIR
    void _attitude_GPSIIR(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    
    // attitude for GPS Block IIR-M
    void _attitude_GPSIIRM(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    
    // attitude for GPS Block IIF
    void _attitude_GPSIIF(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    void _midnight_turn_GPSIIF(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
       
    // attitude for Galileo IOV
    void _attitude_GAL1(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    void _noon_turn_GAL1(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);

    // attitude for Galileo FOC
    void _attitude_GAL2(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    void _noon_turn_GAL2(t_gsatdata& satdata, ColumnVector& i, ColumnVector& j, ColumnVector& k);   
    
    // attitude for BeiDou
    void _attitude_BDS(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    
    // attitude for QZSS
    void _attitude_QZS(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);  

    // attitude for GLO
    void _attitude_GLO(t_gsatdata& satdata, string antype, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    void _midnight_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);
    void _noon_turn_GLOM(t_gsatdata& satdata, double R, ColumnVector& i, ColumnVector& j, ColumnVector& k);  
    
    void _yaw2ijk(t_gsatdata& satdata, double& yaw, ColumnVector& i, ColumnVector& j, ColumnVector& k);
       
   map<string, double> _windUpTime;
   map<string, double> _windUpSum;
   t_gephplan          _ephplan;
   t_gallobj*          _gallobj;
   GRDMPFUNC           _grad_mf;
   
    map<string, double>  _last_beta;
    map<string, double>  _last_yaw;
    map<string, t_gtime> _last_epo;
   
};

} // namespace

#endif //  GPPPMODEL_H
