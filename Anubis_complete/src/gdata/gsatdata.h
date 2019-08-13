
#ifndef  GSATDATA_H
#define  GSATDATA_H

#include <string>
#include <vector>

#include "../newmat/newmat.h"
#include "gdata/gdata.h"
#include "gdata/gobsgnss.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gall/gallnav.h"
#include "gset/gsetproc.h"

using namespace std;

namespace gnut {

class t_gsatdata : public t_gobsgnss {
   
 public:
   t_gsatdata(){};
   t_gsatdata(string site, string sat, const t_gtime& t);
   t_gsatdata(const t_gobsgnss& obs);
   virtual ~t_gsatdata();
   
   void addcrd( const t_gtriple& crd );        // add satellite position
   void addpco( const t_gtriple& pco );        // add satellite pco
   void addvel( const t_gtriple& vel );        // add satellite velocity   
   void addclk( double d );                    // add satellite clocks at the transmision time
   void addecl(map<string, t_gtime>& lastEcl); // determine wheather the satellite is eclipsed (postshadow is considered)
   void	setecl(bool ecl);						// reset the eclipsing time
   
   // add satellite pos, clk and ecl (corrTOT is correction of transmition time)
   // add satellite pos, clk and ecl (corrTOT is correction of transmition time) - low precision of sat pos
   int  addprd(    t_gallnav* gnav, bool corrTOT = true, bool msk_health = true  ); // true  to support PPP/SPP (set by INP:chk_health)
   int  addprd_nav(t_gallnav* gnav, bool corrTOT = true, bool msk_health = false ); // false to support QC (combines INP:chk_health+QC:use_health for Anubis) 
   int  cmpVal(t_gtriple& xyz);
   
   void addele( double d );                // add satellite elevation                           
   void addazi( double d );                // add satellite azimuth                             
   void addrho( double d );                // add satellite rho-vector
   void addres( RESIDTYPE restype, GOBSTYPE type, double res );
   void addmfH( const double& d);
   void addmfW( const double& d);
   void addmfG( const double& d);
   void addwind( const double& d);   
   
   t_gtriple satcrd();                     // get satellite position
   t_gtriple satpco();
   t_gtriple satvel();                     // get satellite velocity
   double  clk();			   // get satellite clocks at the transmision time 
   double  ele();			   // get satellite elevation
   double  ele_deg();			   // get satellite elevation [deg]
   double  azi();			   // get satellite azimuth                        
   double  rho();			   // get satellite rho-vector
   bool    ecl();                          // get eclipsing
   double  mfH();
   double  mfW();
   double  mfG();
   double  wind();
   void	   addslip(const bool& flag);
   bool    islip();
   double  beta();           // Sun elevation relative to orbital plane
   double  orb_angle();      // Orbit angle
   double  yaw(){return _yaw;}          // get yaw angle
   void    yaw(double yaw){_yaw = yaw;} // set yaw angle
   
   vector<double> residuals(RESIDTYPE restype, GOBSTYPE type);
   void           clear_res(RESIDTYPE restype);
   
   void clear();
   bool valid();   
   
 private:
   int  _addprd(t_gallnav* gnav, bool corrTOT = true, bool msk_health = true );          // add satellite pos, clk and ecl (corrTOT is correction of transmition time)
   
   double  _beta();       
   double  _orb_angle();
   
   virtual void _clear();
   virtual bool _valid() const;

   t_gtriple         _satcrd; // satellite position (X,Y,Z)
   t_gtriple         _satpco; // satellite pco
   t_gtriple         _satvel; // satellite velocity   
   double            _clk;    // satellite clocks (precise, time of transmision)
   double            _ele;    // satellite elevation
   double            _azi;    // satellite azimuth
   double            _rho;    // satellite-station geometrical distance
   bool              _eclipse; // eclipsed satellite
   double            _mfH;
   double            _mfW;
   double            _mfG;
   double            _wind; 
   bool              _low_prec; // low precision of sat pos
   bool              _slipf;	// cycle slip flag   
   
   // normalized residuals
   vector<double>    _code_res_norm;
   vector<double>    _phase_res_norm;
   
   // original residuals
   vector<double>    _code_res_orig;
   vector<double>    _phase_res_orig;
   
   double _beta_val;
   double _orb_angle_val;
   double _yaw;
   
};

} // namespace

#endif
