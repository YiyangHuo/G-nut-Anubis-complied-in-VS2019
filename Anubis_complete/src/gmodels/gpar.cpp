
/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 * 
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  This file is part of the G-Nut C++ library.
 
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 3 of
  the License, or (at your option) any later version.
 
  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses>.

-*/

#include <iostream>
#include <cmath>

#include "../../newmat/newmat.h"

#include "gmodels/gpar.h"
#include "gmodels/ggmf.h"

using namespace std;

namespace gnut {   

// constructor
// --------------------------------------------------------
t_gpar::t_gpar(string site, t_type t, int i, string p)
{
   beg = FIRST_TIME;
   end = LAST_TIME;

   this->site = site;
   parType = t;
   index = i;
   prn = p;
   value( 0.0 );
   apriori( 0.0 );
}

t_gpar::t_gpar()
{
   beg = FIRST_TIME;
   end = LAST_TIME;   
}

// t_gpar destructor
// ---------------------------------------------------
t_gpar::~t_gpar()
{
}

// setting mapping function for ZTD
// ---------------------------------------------------
void t_gpar::setMF(ZTDMPFUNC MF)
{
  _mf_ztd = MF;   
}

// setting mapping function for GRD
// ---------------------------------------------------
void t_gpar::setMF(GRDMPFUNC MF)
{
  _mf_grd = MF;
}
   
// Partial derivatives 
// ----------------------------------------------------
double t_gpar::partial(t_gsatdata& satData, t_gtime& epoch, t_gtriple ground, t_gobs& gobs)
{
  double mfw, dmfw, mfh, dmfh;
  mfw = dmfw = mfh = dmfh = 0.0;   
  
  switch( parType ){
    case CRD_X : return (value() - satData.satcrd().crd(0))/satData.rho();
    case CRD_Y : return (value() - satData.satcrd().crd(1))/satData.rho();
    case CRD_Z : return (value() - satData.satcrd().crd(2))/satData.rho();
    case CLK   : return 1.0;
    case IFCB_F3 : if (t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_3 
                       && prn == satData.sat()) return -1.0; else return 0.0;
    case IFCB_F4 : if (t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_4
                       && prn == satData.sat()) return -1.0; else return 0.0;
    case IFCB_F5 : if (t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_5 
                       && prn == satData.sat()) return -1.0; else return 0.0;    
    case IFB_C3 : if (gobs.is_code() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_3) return 1.0; else return 0.0;
    case IFB_C4 : if (gobs.is_code() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_4) return 1.0; else return 0.0;
    case IFB_C5 : if (gobs.is_code() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_5) return 1.0; else return 0.0;    
    case CLK_SAT: if (satData.sat() == this->prn) return -1.0; else return 0.0;
    case TRP   :
           _getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);
           return mfw;
    case SION  :
      {       
        double f1 = G01_F;
        double fk = satData.frequency(gobs.band());        
        double alfa = 0.0;
        if(gobs.is_phase() && prn == satData.sat()) {alfa =  - (f1*f1) / (fk*fk);}
        if(gobs.is_code()  && prn == satData.sat()) {alfa =    (f1*f1) / (fk*fk);}
        return alfa;
      }
    case VION  :
      {
        double f1 = G01_F;
        double fk = satData.frequency(gobs.band());
        double mf = 1.0 / sqrt( 1.0 - pow(R_SPHERE/(R_SPHERE+450000.0) * sin(G_PI/2.0 - satData.ele()),2) );
        double alfa = 0.0;
        if(gobs.is_phase() && prn == satData.sat()) {alfa =  - (f1*f1) / (fk*fk);}
        if(gobs.is_code()  && prn == satData.sat()) {alfa =    (f1*f1) / (fk*fk);}
//cout << "Partial " << satData.sat() << " " << prn << "  " << gobs2str(gobs.gobs()) << "  " << gobs.band() << "  " << alfa << "  " << f1 << "  " << fk << endl;
        return alfa * mf;
      }
    case P1P2G_REC :
      {
        double f1 = G01_F;
        double fk = satData.frequency(gobs.band());        
        double alfa = (f1*f1) / (fk*fk);        
        double beta = (G02_F*G02_F)/(G01_F*G01_F - G02_F*G02_F);        
        FREQ_SEQ freq = t_gsys::band2freq(satData.gsys(), gobs.band());        
        if( satData.gsys() == GPS && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2) ) {return - alfa * beta;}
        else return 0.0;
      }
    case P1P2E_REC :
      {
        double f1 = E01_F;
        double fk = satData.frequency(gobs.band());        
        double alfa = (f1*f1) / (fk*fk);        
        double beta = (E05_F*E05_F)/(E01_F*E01_F - E05_F*E05_F);        
        FREQ_SEQ freq = t_gsys::band2freq(satData.gsys(), gobs.band());        
        if( satData.gsys() == GAL && gobs.is_code() && (freq == FREQ_1 || freq == FREQ_2) ) {return - alfa * beta;}
        else return 0.0;
      }    
    case GRD_N :
        _getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);
     if (_mf_grd == CHEN_HERRING) {      
         double sinel = sin( satData.ele() );
         double tanel = tan( satData.ele() );
         double cosaz = cos( satData.azi() );
         return (1.0/(sinel*tanel+0.0032)) * cosaz;
      }else if (_mf_grd == TILTING)  {
         double cosaz = cos( satData.azi() );         
         return dmfw * cosaz;
      }else if (_mf_grd == BAR_SEVER)  {
         double tanel = tan( satData.ele() );
         double cosaz = cos( satData.azi() );
         return mfw * (1.0/tanel) * cosaz;
      }else cerr << "Grad N mapping function is not set up correctly!!!" << endl;
    case GRD_E :
         _getmf(satData, ground, epoch, mfw, mfh, dmfw, dmfh);     
         if (_mf_grd == CHEN_HERRING) {        
         double sinel = sin( satData.ele() );
         double tanel = tan( satData.ele() );         
         double sinaz = sin( satData.azi() );
         return (1.0/(sinel*tanel+0.0032)) * sinaz;
      }else if (_mf_grd == TILTING) {                     
         double sinaz = sin( satData.azi() );         
         return dmfw * sinaz;
      }else if (_mf_grd == BAR_SEVER) {
         double tanel = tan( satData.ele() );
         double sinaz = sin( satData.azi() );
         return mfw * (1/tanel) * sinaz;        
      } else cerr << "Grad E mapping function is not set up correctly!!!" << endl;     
     
    case AMB_IF  :
          if (gobs.is_phase() && prn == satData.sat()) return 1.0;
          else return 0.0;
    case AMB_L1  :
          if (gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_1 && prn == satData.sat()) return 1.0;
          else return 0.0;
    case AMB_L2  :
          if (gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_2 && prn == satData.sat()) return 1.0;
          else return 0.0;  
    case AMB_L3  :
          if (gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_3 && prn == satData.sat()) return 1.0;
          else return 0.0;     
    case AMB_L4  :
          if (gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_4 && prn == satData.sat()) return 1.0;
          else return 0.0;
    case AMB_L5  :
          if (gobs.is_phase() && t_gsys::band2freq(satData.gsys(), gobs.band()) == FREQ_5 && prn == satData.sat()) return 1.0;
          else return 0.0;        
    case GLO_ISB :
          if (satData.gsys() == GLO) return 1.0;
          else return 0.0;
   case GLO_IFCB :
          if (!gobs.is_phase() && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
          else return 0.0;
   case GLO_IFPB :
          if (gobs.is_phase() && satData.gsys() == GLO && prn == satData.sat()) return 1.0;
          else return 0.0;     
    case GAL_ISB :
          if (satData.gsys() == GAL) return 1.0;
          else return 0.0;
    case BDS_ISB :
          if (satData.gsys() == BDS) return 1.0;
          else return 0.0;
    case QZS_ISB :
          if (satData.gsys() == QZS) return 1.0;
          else return 0.0;
  }
  return 0.0;
}


// Operators for t_gpar
// -------------------------------------------------
bool t_gpar::operator==(const t_gpar& par) const
{
   if ( parType == par.parType      &&
      //  prn.compare(par.prn)   == 0 &&
      //  site.compare(par.site) == 0 &&
   beg == par.beg &&
   end == par.end ){
      return true;
   } else {
      return false;
   }
   
}

// Doplnil gabo aby par mohla byt pouzita ako klic v mape
bool t_gpar::operator <(const t_gpar& par) const
{
    if (( parType <  par.parType )                                        ||
        ( parType == par.parType && prn <  par.prn )                      ||
        ( parType == par.parType && prn == par.prn && site <  par.site )  ||
        ( parType == par.parType && prn == par.prn && site == par.site && beg < par.beg ) ||
        ( parType == par.parType && prn == par.prn && site == par.site && beg == par.beg && end < par.end )
        ){
       return true;
    }
    return false;
}

t_gpar t_gpar::operator-(const t_gpar& p) const
{
   t_gpar par = (*this);
   par.value( value() - p.value() );
   return par;
}

t_gpar t_gpar::operator+(const t_gpar& p) const
{
   t_gpar par = (*this);
   par.value( value() + p.value() );
   return par;
}

// Setting begin and end time for validity object
// -------------------------------------------------
void t_gpar::setTime(const t_gtime& t1, const t_gtime& t2)
{
   this->beg = t1;
   this->end = t2;
}

// get data type
// ----------
string t_gpar::str_type() const
{
  string type;
  switch( parType ){
   case  CRD_X    :  type = "CRD_X";     break;
   case  CRD_Y    :  type = "CRD_Y";     break;
   case  CRD_Z    :  type = "CRD_Z";     break;
   case  TRP      :  type = "TRP"  ;     break;
   case  SION     :  type = "SION_" + prn ;     break;
   case  VION     :  type = "VION_" + prn ;     break;  
   case  CLK      :  type = "CLK"  ;     break;
   case  IFCB_F3  :  type = "IFCB_F3_" + prn;     break;
   case  IFCB_F4  :  type = "IFCB_F4_" + prn;     break;
   case  IFCB_F5  :  type = "IFCB_F5_" + prn;     break;    
   case  IFB_C3   :  type = "IFB_F3_" + prn;     break;
   case  IFB_C4   :  type = "IFB_F4_" + prn;     break;
   case  IFB_C5   :  type = "IFB_F5_" + prn;     break;    
   case  CLK_SAT  :  type = "CLK_SAT_" + prn;   break;
   case  AMB_IF   :  type = "AMB_IF_" + prn ;   break;
   case  AMB_L1   :  type = "AMB_L1_" + prn;    break;
   case  AMB_L2   :  type = "AMB_L2_" + prn;    break;
   case  AMB_L3   :  type = "AMB_L3_" + prn;    break;
   case  AMB_L4   :  type = "AMB_L4_" + prn;    break;  
   case  AMB_L5   :  type = "AMB_L5_" + prn;    break;
   case  FCBS_IF  :  type = "FCBS_IF";   break;
   case  FCBS_L1  :  type = "FCBS_L1";   break;
   case  FCBS_L2  :  type = "FCBS_L2";   break;
   case  FCBS_L3  :  type = "FCBS_L3";   break;
   case  FCBS_L4  :  type = "FCBS_L4";   break;
   case  FCBS_L5  :  type = "FCBS_L5";   break;
   case  FCBR_IF  :  type = "FCBR_IF";   break;
   case  FCBR_L1  :  type = "FCBR_L1";   break;
   case  FCBR_L2  :  type = "FCBR_L2";   break;
   case  FCBR_L3  :  type = "FCBR_L3";   break;
   case  FCBR_L4  :  type = "FCBR_L4";   break;
   case  FCBR_L5  :  type = "FCBR_L5";   break;
   case  CLK_ICB  :  type = "CLK_ICB";   break;
   case  CLUSTERB :  type = "CLUSTERB";  break;
   case  GRD_N    :  type = "GRD_N";     break;
   case  GRD_E    :  type = "GRD_E";     break;
   case  GLO_ISB  :  type = "GLO_ISB";   break;
   case  GLO_IFCB :  type = "GLO_IFCB";  break;
   case  GLO_IFPB :  type = "GLO_IFPB";  break;
   case  GAL_ISB  :  type = "GAL_ISB";   break;
   case  BDS_ISB  :  type = "BDS_ISB";   break;
   case  QZS_ISB  :  type = "QZS_ISB";   break;
   case  P1P2G_REC :  type = "P1P2G_REC";  break;
   case  P1P2E_REC :  type = "P1P2E_REC";  break;    
   default        :  type = "UNDEF";
  }
   
return type;
}

// get ZTD mf according to settings   
void t_gpar::_getmf(t_gsatdata& satData, const t_gtriple& crd, const t_gtime& epoch, double& mfw, double& mfh, double& dmfw, double& dmfh)
{
   if(parType != TRP && parType != GRD_N && parType != GRD_E) return;

   double ele = satData.ele();
     
   if (_mf_ztd == COSZ) {
      mfw = mfh = 1.0/sin(ele);
   }else if (_mf_ztd == GMF)  {
      t_gmf mf;
      mf.gmf( epoch.mjd(), crd[0], crd[1], crd[2], G_PI/2.0-ele,
         mfh, mfw, dmfh, dmfw );
   } else cerr << "ZTD mapping function is not set up correctly!!!" << endl;
}   
   
} // namespace
