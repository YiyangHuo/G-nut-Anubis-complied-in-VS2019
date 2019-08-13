
#ifndef GPAR_H
#define GPAR_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: parametres class
  Version: $ Rev: $

  2011-04-18 /PV: created

-*/

#include <string>

#include "gdata/gsatdata.h"
#include "gutils/gtriple.h"
#include "gset/gsetproc.h"

using namespace std;

namespace gnut {   
   
class t_gpar
{
 public:
   
   enum t_type{	
     CRD_X, CRD_Y, CRD_Z,                                        // coordinates
     TRP, GRD_N, GRD_E, SION, VION,                              // atmospheric parameters
     CLK, CLK_SAT,                                               // clocks
     IFB_C3, IFB_C4, IFB_C5, IFCB_F3, IFCB_F4, IFCB_F5,          // inter-freq. code biases for FREQ_3, FREQ_4, FREQ_5, inter-freq. clock bias for GPS FREQ_3
     CLK_ICB, CLUSTERB,                                          // initial clock bias, cluster-dependent bias
     AMB_IF, AMB_L1, AMB_L2, AMB_L3, AMB_L4, AMB_L5,             // ambiguities for indiv. freq. (number indicates freq not band)
     FCBS_IF, FCBS_L1, FCBS_L2, FCBS_L3, FCBS_L4, FCBS_L5,       // satellite fractional cycle biases for indiv. freq.
     FCBR_IF, FCBR_L1, FCBR_L2, FCBR_L3, FCBR_L4, FCBR_L5,	     // receiver fractional cycle biases for indiv. freq.
     GLO_ISB, GLO_IFCB, GLO_IFPB, GAL_ISB, BDS_ISB, QZS_ISB,     // multi-GNSS
     P1P2G_REC, P1P2E_REC, P1P2R_REC, P1P2C_REC                  // GNSS-specific receiver code DCB P1-P2
   };
   
   t_gpar(string site, t_type t, int i, string p);
   t_gpar();
   ~t_gpar();

   bool  operator==(const t_gpar&) const;
   bool  operator <(const t_gpar&) const;
   t_gpar operator-(const t_gpar&) const;
   t_gpar operator+(const t_gpar&) const;   
   
   void setTime(const t_gtime&, const t_gtime&);
   
   double partial(t_gsatdata&, t_gtime&, t_gtriple, t_gobs& gobs);

   void   value( double val ){ _value = val; }
   double value()const{ return _value; }
   
   void   apriori( double apr ){ _apriori = apr; }
   double apriori()const{ return _apriori; }   

   t_type parType;
   int index;
   string prn;
   string site;
   t_gtime beg;
   t_gtime end;
   t_gtime stime;
   double  aprval;
   string str_type() const;
   
  void setMF(ZTDMPFUNC MF);
  void setMF(GRDMPFUNC MF);
   
 protected:
   double    _value;   // value
   ZTDMPFUNC _mf_ztd;      // mapping function for ZTD
   GRDMPFUNC _mf_grd;      // mapping function for GRD   
   double    _apriori;
   void      _getmf(t_gsatdata& satData, const t_gtriple& crd, const t_gtime& epoch, double& mfw, double& mfh, double& dmfw, double& dmfh);
};

} // namespace

#endif
