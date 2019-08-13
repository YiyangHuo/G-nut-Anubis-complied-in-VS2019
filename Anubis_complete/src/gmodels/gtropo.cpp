
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

#include <cmath>
#include <iomanip>

#include "gmodels/gtropo.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"
#include "gutils/gsysconv.h"

using namespace std;

namespace gnut {   

// ---------
// Base class
// ----------
/*
t_gtropo::t_gtropo(string site)
 : _site(site),
   _ell(0.0,0.0,0.0),
   _nwm(0),
   _met(0)
{
  if( !site.empty() ) cerr << "gtropo model site[" << site << "]\n"; 
}
*/

t_gtropo::t_gtropo() // : t_gtropo::t_gtropo("")
 : _site(""),
   _ell(0.0,0.0,0.0),
   _nwm(0),
   _met(0)
{}


t_gtropo::~t_gtropo()
{}

   
// settings   
// ----------
void t_gtropo::nwm( t_gallsurf* nwm ){ _nwm = nwm; }
void t_gtropo::met( t_gallprod* met ){ _met = met; }


// initialize site (from product)
// ----------
bool t_gtropo::init(string site)
{
  gtrace("t_gtropo::init");

  set<string> sites = _met->prod_sites();

  if( !_met         ){ cerr << "gtropo - Error: cannot be initiated, no product found!\n"; return false;}
  if( sites.size()<1){ cerr << "gtropo - Error: cannot be initiated, no sites found!\n";   return false;}

  if( !site.empty() ){
    _site = site;
    if( sites.find(site) == sites.end() ){
      cerr << "gtropo - Error: cannot be initiated, site ["+site+"] not found!\n"; return false;
    }
  }

  // collect sites & coordinates
  for( set<string>::iterator itSIT = sites.begin(); itSIT != sites.end(); ++itSIT ){
    string s = *itSIT;

    if( !site.empty() && s != site ) continue;
    set<t_gtime> epo_pos = _met->prod_epochs(s,t_gdata::POS); // get epocsh for coordinates
    
    // CRD
    if( epo_pos.size() == 0 ){ cerr << "gtropo_init - Error: CRD not in the product!\n"; return false; }
    for( set<t_gtime>::iterator itEPO = epo_pos.begin(); itEPO != epo_pos.end(); ++itEPO ){
      t_gtime epo = *itEPO;
      shared_ptr<t_gprodcrd> pt_crd = dynamic_pointer_cast<t_gprodcrd>(_met->get(s, t_gdata::POS, epo));
     
      if( !pt_crd ){ cerr << s << " reading CRD epochs -> not found\n"; continue; }
      else{ t_gtriple xyz = pt_crd->xyz(), ell; xyz2ell( xyz, ell, false);   // RADIANS
//        cerr << s  << " reading CRD epochs" << epo.str_ymdhms(" ");
        map<string,t_gtriple>::iterator itS = _m_sites.begin();
        if( (itS = _m_sites.find(s)) == _m_sites.end() ) _m_sites[s] = ell;  // first case
        else if( fabs(itS->second[2] - ell[2]) > 1.0 ) cerr << s + " - Warning: high difference between epochs\n";
      }
    }

    // MET (RINEX)
    set<t_gtime> epo_met = _met->prod_epochs(s,t_gdata::MET); // get epochs for coordinates
//  if( epo_met.size() < 0 ){ cerr << "gtropo_init - Error: MET not in the product!\n"; return false; }
    for( set<t_gtime>::iterator itEPO = epo_met.begin(); itEPO != epo_met.end(); ++itEPO ){
      t_gtime epo = *itEPO;
      shared_ptr<t_gprodmet> pt_trp = dynamic_pointer_cast<t_gprodmet>(_met->get(s, t_gdata::MET, epo));
      
      if( !pt_trp ){ cout << s << " reading MET epochs -> not found\n"; continue; }
#ifdef DEBUG
      else           cout << s << " reading MET epochs " + epo.str_ymdhms(" ")
                          << fixed << setprecision(3) << " " << pt_trp->pres() << " " << pt_trp->epre()
	                                              << " " << pt_trp->temp() << endl;
#endif       
      _m_prods[s][epo] = pt_trp;     // first or any case (just rewrite)       
    }

    // TRP (SINEX)
    set<t_gtime> epo_trp = _met->prod_epochs(s,t_gdata::TRP); // get epochs for coordinates
    if( epo_trp.size() == 0 && epo_met.size() == 0 ){ cerr << "gtropo_init - Error: TRP not in the product!\n"; return false; }
    for( set<t_gtime>::iterator itEPO = epo_trp.begin(); itEPO != epo_trp.end(); ++itEPO ){
      t_gtime epo = *itEPO;
      shared_ptr<t_gprodtrp> pt_trp = dynamic_pointer_cast<t_gprodtrp>(_met->get(s, t_gdata::TRP, epo));
      
      if( !pt_trp ){ cout << s << " reading TRP epochs -> not found\n"; continue; }
#ifdef DEBUG
      else           cout << s << " reading TRP epochs " + epo.str_ymdhms(" ")
	                  << fixed << setprecision(3) << " " << pt_trp->ztd() << " " << pt_trp->zhd()
	                                              << " " << pt_trp->zwd() << endl;
#endif
      _m_prods[s][epo] = pt_trp;     // first or any case (just rewrite)       
    }

  }
  return true;
}


// temporal interpolation
// ----------------------
double t_gtropo::_interp_temporal(string site, t_gtime epo, METEO_ID type)
{
  gtrace("t_gtropo::interp_temporal");
   
  double res = NWM_UNKNOWN;
  map<t_gtime,double> map_gtime;

  if( _m_prods.find(site) == _m_prods.end() ) return res;
  map<t_gtime, shared_ptr<t_gprodmet>>::iterator itEPO = _m_prods[site].upper_bound(epo);
  if( itEPO == _m_prods[site].end() || itEPO == _m_prods[site].begin() ) return res; // out of boundary
  set<METEO_ID> params = itEPO->second->get_params();
  if( params.find( type ) == params.end() ) return res;
  --itEPO; // 1st value before EPO

  if( distance(_m_prods[site].begin(),itEPO) < 1 || distance(itEPO,_m_prods[site].end()) < 3 ){ // 3 for spline: -1+3!
     switch(type){
      case METEO_ZHD : map_gtime[  (itEPO)->first] = itEPO->second->zhd();
                       map_gtime[(++itEPO)->first] = itEPO->second->zhd();
	               break;
      case METEO_ZWD : map_gtime[  (itEPO)->first] = itEPO->second->zwd();
                       map_gtime[(++itEPO)->first] = itEPO->second->zwd();
	               break;
      case METEO_PRES: map_gtime[  (itEPO)->first] = itEPO->second->pres();
                       map_gtime[(++itEPO)->first] = itEPO->second->pres();
	               break;
      case METEO_EPRE: map_gtime[  (itEPO)->first] = itEPO->second->epre();
                       map_gtime[(++itEPO)->first] = itEPO->second->epre();
	               break;
      case METEO_TEMP: map_gtime[  (itEPO)->first] = itEPO->second->temp();
                       map_gtime[(++itEPO)->first] = itEPO->second->temp();
	               break;
      default:         map_gtime[  (itEPO)->first] = itEPO->second->met(type);
                       map_gtime[(++itEPO)->first] = itEPO->second->met(type);
     }
     _interp.linear( map_gtime, epo, res );
#ifdef DEBUG     
     cout << type << epo.str_ymdhms( " interp[linear]: ") << " " << itEPO->first.str_ymdhms()  << " " << res
          << " " << map_gtime.size() << " " <<  distance(itEPO,_m_prods[site].end()) << endl;
#endif
  }else{
     switch(type){
      case METEO_ZHD : map_gtime[(--itEPO)->first] = itEPO->second->zhd();
                       map_gtime[(++itEPO)->first] = itEPO->second->zhd();
                       map_gtime[(++itEPO)->first] = itEPO->second->zhd();
                       map_gtime[(++itEPO)->first] = itEPO->second->zhd();
	               break;
      case METEO_ZWD : map_gtime[(--itEPO)->first] = itEPO->second->zwd();
                       map_gtime[(++itEPO)->first] = itEPO->second->zwd();
                       map_gtime[(++itEPO)->first] = itEPO->second->zwd();
	               map_gtime[(++itEPO)->first] = itEPO->second->zwd();
	               break;
      case METEO_PRES: map_gtime[(--itEPO)->first] = itEPO->second->pres();
                       map_gtime[(++itEPO)->first] = itEPO->second->pres();
                       map_gtime[(++itEPO)->first] = itEPO->second->pres();
                       map_gtime[(++itEPO)->first] = itEPO->second->pres();
	               break;
      case METEO_EPRE: map_gtime[(--itEPO)->first] = itEPO->second->epre();
                       map_gtime[(++itEPO)->first] = itEPO->second->epre();
                       map_gtime[(++itEPO)->first] = itEPO->second->epre();
                       map_gtime[(++itEPO)->first] = itEPO->second->epre();
	               break;
      case METEO_TEMP: map_gtime[(--itEPO)->first] = itEPO->second->temp();
                       map_gtime[(++itEPO)->first] = itEPO->second->temp();
                       map_gtime[(++itEPO)->first] = itEPO->second->temp();
                       map_gtime[(++itEPO)->first] = itEPO->second->temp();
	               break;
      default:         map_gtime[(--itEPO)->first] = itEPO->second->met(type);
                       map_gtime[(++itEPO)->first] = itEPO->second->met(type);
                       map_gtime[(++itEPO)->first] = itEPO->second->met(type);
                       map_gtime[(++itEPO)->first] = itEPO->second->met(type);
     }
     _interp.spline( map_gtime, epo, res );
#ifdef DEBUG
     cout << type << epo.str_ymdhms( " interp[spline]: ") << " " << itEPO->first.str_ymdhms()  << " " << res
          << " " << map_gtime.size() << " " <<  distance(itEPO,_m_prods[site].end()) << endl;
#endif     
  }
  return res;
}


// get ZHD
// ----------
double t_gtropo::getZHD(const t_gtriple& ell, const t_gtime& epo) // ell v RADIANECH !! TREBA SJEDNOTIT
{   
  gtrace("t_gtropo::getZHD");  

  // in future this may be significantly optimized!
  map<t_gtime, map<MET_DATA, double> > map_nwm;

  double res = NWM_UNKNOWN;

  // TROPO product (SINEX_TRO/RINEXM) ---> PRIORITY ZHD FROM IN SITU!!!
  if( _met ){
    // a nearby point 
    if( _site.empty() ){
      for(map<string, t_gtriple>::iterator itSIT = _m_sites.begin(); itSIT != _m_sites.end(); ++itSIT ){
        t_gtriple crd = itSIT->second;

	if( fabs(crd[2] - ell[2]) < 2.0           &&    // vertical:     2 m
	    fabs(crd[0] - ell[0]) < 10.0/R_SPHERE &&    // horizontal: ~10 m
	    fabs(crd[1] - ell[1]) < 10.0/R_SPHERE ){    // horizontal: ~10 m

	  res = _interp_temporal(itSIT->first,epo,METEO_ZHD); // ZHD
//          cout << epo.str_ymdhms("MET data: ") << " ZHD " << setw(12) << res << " kinematic" << endl;
	}
      }
    }else{
      res = _interp_temporal(_site,epo,METEO_ZHD);         // try ZHD
      if( res == NWM_UNKNOWN ){
	double P = _interp_temporal(_site,epo,METEO_PRES); // try PRES instead
	if( P != NWM_UNKNOWN ){
	   // Saastamoinen(1972)
           res = (0.002277 * P) / (1.0 - 0.00266 * cos(2.0*ell[0]) - 0.00000028 * ell[2]);
	}
      }
//      cout << epo.str_ymdhms("MET data: ") << " ZHD " << setw(12) << res << "  site: " << _site << endl;
    }
//    cout << epo.str_ymdhms("MET data: ")
//	 << fixed << setprecision(3) << " ZHD " << setw(12) << res << endl;
    if( res != NWM_UNKNOWN ) return res;
  }

  // NWM tropo/meteo grid
  if( _nwm ){
    t_gtriple tmp( ell[0]*R2D, ell[1]*R2D, ell[2]);             // TEMPORARILY TO DEGREE !!
//     cout << " ell = " << ell[0]*R2D << " " << ell[1]*R2D << " " << ell[2]  << " "
//          << _nwm->epochs().size() << " " << _nwm->pairs(epo).size() << endl;
    _nwm->interpolate_point(tmp, epo, epo, &map_nwm );
     
    if( map_nwm.find(epo)      != map_nwm.end() &&
        map_nwm[epo].find(ZHD) != map_nwm[epo].end() )
    {
//      cout << epo.str_ymdhms("NWM data: ")
//	   << fixed << setprecision(3) << " ZHD " << setw(12) << map_nwm[epo][ZHD] << endl;
      return map_nwm[epo][ZHD];
    }
  }

  cerr << "t_gtropo:getZHD not available data " << epo.str_ymdhms() << ", using default: 2.3m \n";

  return 2.3; // NWM_UNKNOWN;
}


// get ZWD
// ----------
double t_gtropo::getZWD(const t_gtriple& ell, const t_gtime& epo) // ell v RADIANECH !! TREBA SJEDNOTIT
{ 
  gtrace("t_gtropo::getZWD");

  // in future this may be significantly optimized!
  map<t_gtime, map<MET_DATA, double> > map_nwm;
   
  // NWM tropo/meteo grid -- PRIORITY ZWD FROM MODEL  !
  if( _nwm ){
    t_gtriple tmp( ell[0]*R2D, ell[1]*R2D, ell[2]);             // TEMPORARILY TO DEGREE !!
    _nwm->interpolate_point(tmp, epo, epo, &map_nwm );

    if( map_nwm.find(epo)      != map_nwm.end() &&
        map_nwm[epo].find(ZWD) != map_nwm[epo].end() )
    {
//      cout << epo.str_ymdhms("NWM data: ")
//	   << fixed << setprecision(3) << " ZWD " << setw(12) << map_nwm[epo][ZWD] << endl;
      return map_nwm[epo][ZWD];
    }
  }

  double res = NWM_UNKNOWN;

  // TROPO product (SINEX_TRO/RINEXM)
  if( _met ){

    // a nearby point 
    if( _site.empty() ){
      for(map<string, t_gtriple>::iterator itSIT = _m_sites.begin(); itSIT != _m_sites.end(); ++itSIT ){
        t_gtriple crd = itSIT->second;

	if( fabs(crd[2] - ell[2]) < 2.0           &&    // vertical:     2 m
	    fabs(crd[0] - ell[0]) < 10.0/R_SPHERE &&    // horizontal: ~10 m
	    fabs(crd[1] - ell[1]) < 10.0/R_SPHERE ){    // horizontal: ~10 m

	  res = _interp_temporal(itSIT->first,epo,METEO_ZWD); // ZWD
//          cout << epo.str_ymdhms("MET data: ") << " ZWD " << setw(12) << res << " kinematic" << endl;
	}
      }
    }else{
      res = _interp_temporal(_site,epo,METEO_ZWD);  // ZWD
//          cout << epo.str_ymdhms("MET data: ") << " ZWD " << setw(12) << res << endl;
      if( res == NWM_UNKNOWN ){
	double T = _interp_temporal(_site,epo,METEO_TEMP); // T
	double E = _interp_temporal(_site,epo,METEO_EPRE); // E
//        cout << epo.str_ymdhms("MET data: ") << fixed << setprecision(3) << " E: " << E << " T: " << T << endl;
	if( E != NWM_UNKNOWN && T != NWM_UNKNOWN ){
	   
	  // Askne-Nordius (1987)
          double dT = 6.5;                                    // [K/km] standard temperature lapse rate
	  double dE = 3.0;                                    // [-]    standard water vapor decay
          double gm = 9.784*( 1 - 0.00266 * cos(2*ell[0])
		                - 0.00028 * ell[2]/1000.0);   // [m/s] Askne-Nordius original
          double Tm = T * ( 1 - dT/1000.0 * Rd/gm/(dE + 1));  // [K]   Mean temperature
      
	  res = 1e-06 * (K2_THA - K1_THA*Eps + K3_THA/Tm) * Rd/gm/(dE + 1) * E; // [m]

//        cout << epo.str_ymdhms("MET data: ") << " ZWD " << setw(12) << res << "  site: " << _site
//  	       << " " << " Tm: " << Tm << "  hgt: " << ell[2] << " lat: " << ell[0] << endl;
	}
      }
    }
     
//    cout << epo.str_ymdhms("NWM data: ")
//	 << fixed << setprecision(3) << " ZWD " << setw(12) << res << endl;
    if( res != NWM_UNKNOWN ) return res;
  }
  
  cerr << "t_gtropo:getZWD not available data " << epo.str_ymdhms() << ", using default: 0.0m \n";

  return 0.0; // NWM_UNKNOWN;
}

// ----------------------------- SITE MODELS ------------------------------------

// ---------
// Saastamoinen 1972
// ----------
double t_saast::getSTD(const double& ele, const double& height)
{    
  gtrace("t_saast::getSTD"); 
   
  double pp =  1013.25 * pow(1.0 - 2.26e-5 * height, 5.225) ;
  double TT =  18.0 - height * 0.0065 + 273.15;
  double hh =  50.0 * exp(-6.396e-4 * height);
  double ee =  hh / 100.0 * exp(-37.2465 + 0.213166*TT - 0.000256908*TT*TT);

#ifdef DEBUG
  cout << "meteo data: " << fixed << setprecision(3)
       << " pp " << setw(12) << pp
       << " TT " << setw(12) << TT
       << " hh " << setw(12) << hh
       << " ee " << setw(12) << ee << endl;
#endif
      
  double h_km = height / 1000.0;

  if (h_km < 0.0) h_km = 0.0;
  if (h_km > 5.0) h_km = 5.0;
  int    ii   = int(h_km + 1);
  double href = ii - 1;
	
  double bCor[6];
  bCor[0] = 1.156;
  bCor[1] = 1.006;
  bCor[2] = 0.874;
  bCor[3] = 0.757;
  bCor[4] = 0.654;
  bCor[5] = 0.563;

  double BB = bCor[ii-1] + (bCor[ii]-bCor[ii-1]) * (h_km - href);
	
  double zen   = G_PI/2.0 - ele;
  double delay = (0.002277/cos(zen)) * (pp + ((1255.0/TT)+0.05)*ee - BB*(tan(zen)*tan(zen)));
  return delay;
}

// ----------
double t_saast::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{   
  gtrace("t_saast::getZHD");
   
  double P, T, N;
   
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P,T,N);
   
  double delay = (0.002277 * P) /
                 (1.0 - 0.00266 * cos(2.0*Ell[0]) - 0.00000028 * Ell[2]);

#ifdef DEBUG
  cout << "GPT data: " << fixed << setprecision(3)
       << " P " << setw(12) << P
       << " T " << setw(12) << T
       << " N " << setw(12) << N
       << " ZHD " << setw(12) << delay << endl;           
#endif   
   
  return delay;
}

// ----------
double t_saast::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{   
  gtrace("t_saast::getZWD");   
   
  double P, T, N;   
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P, T, N);

  T += 273.15;   
   
  double hh =  50.0 * exp(-6.396e-4 * Ell[2]);
  double e =  hh / 100.0 * exp(-37.2465 + 0.213166*T - 0.000256908*T*T);   

  double delay = (0.0022768*e*(1255.0/T + 0.05))
               / (1.0 - 0.00266 * cos(2.0*Ell[0]) - 0.00000028 * Ell[2]);
  return delay;
}


// ---------
// Davis 1985 model (Saast + Thayer 1974)
// ----------
double t_davis::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_davis::getZHD");   
   
  double P, T, N;
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P,T,N);
   
  double delay = (0.0022768 * P)
               / (1.0 - 0.0026 * cos(2.0*Ell[0]) - 0.00000028 * Ell[2]);

#ifdef DEBUG
  cout << "GPT data: " << fixed << setprecision(3)
       << " P " << setw(12) << P
       << " T " << setw(12) << T
       << " N " << setw(12) << N
       << " ZHD " << setw(12) << delay << endl;           
#endif

  return delay;
}

// ---------
double t_davis::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_davis::getZWD");   
   
  // not implemented yet
  return 0.0;
}


// ----------
// Hopfield model
// ----------
double t_hopf::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_hopf::getZHD");   
   
  double P, T, N;

  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P, T, N);
     
  double delay = (1e-06/5.0)*(77.64*(P/T))*(40136.0+148.72*T);
  return delay;
}

// ----------
double t_hopf::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_hopf::getZWD");   
   
  double P, T, N;
  double hh =  50.0 * exp(-6.396e-4 * Ell[2]);
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P, T, N);
      
  double e =  hh / 100.0 * exp(-37.2465 + 0.213166*T - 0.000256908*T*T);

  double delay = (1.e-06/5.0)*((-12.96)*(e/T)+(3.718*1.e05)*(e/(T*T)))*11000.0;
  return delay;
}


// ---------
// Baby model
// ---------
double t_baby::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_baby::getZHD");   
   
  double P, T, N;
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P, T, N);

  T += 273.15;    // [K]

  double gs = 9.81;         // [ms^2] surface gravity !!!
  double rs = A_WGS+Ell[2];
  double sigma = Eps/T;
  double mu = gs/(Rd*Eps)*(1.0-(2.0/(rs*sigma)));
      
  double delay = (0.022277*P/gs)*(1.0+(2.0/(rs*sigma*(mu+1.0))));     
  return delay;
}

// ---------
double t_baby::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_baby::getZWD");   
   
  // not implemeted yet
  return 0.0;
}


// Chao model
// --------
double t_chao::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_chao::getZHD");   
   
  // not implemeted yet
  return 0.0;
}


// Chao model (ZWD), Mendes pp. 85
// --------
double t_chao::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_chao::getZWD");   
   
  double P, T, N, alpha;
  double hh =  50.0 * exp(-6.396e-4 * Ell[2]);   
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P, T, N);

  T += 273.15;
  alpha = 0.0065; // model is not very sensitive to the temperature lapse rate   
     
  double e =  hh / 100.0 * exp(-37.2465 + 0.213166*T - 0.000256908*T*T);   
      
  double delay = 4.70*100.0*(pow(e,1.23)/(T*T))+1.71*1.e6*(pow(e,1.46)/(T*T*T))*alpha;
  return delay;
}


// ---------
// Ifadis model
// ---------
double t_ifad::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_ifad::getZHD");   
   
  // not implemeted yet
  return 0.0;
}

// ---------
double t_ifad::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_ifad::getZWD");   
   
  double P, T, N;   
  double hh =  50.0 * exp(-6.396e-4 * Ell[2]);   
  _gpt.gpt_v1( epoch.mjd(), Ell[0], Ell[1], Ell[2], P, T, N);
   
  T += 273.15;
      
  double e =  hh / 100.0 * exp(-37.2465 + 0.213166*T - 0.000256908*T*T);   
      
  double delay = 0.00554-0.0880*1.e-4*(P - 1000.0)+0.272*1.e-4*e+2.771*(e/T);
  return delay;
}


// ---------
// Askne Nordius model
// ---------
double t_askn::getZHD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_askn::getZHD");   
   
  return 0.0;   // not implemeted yet
}

// ---------
double t_askn::getZWD(const t_gtriple& Ell, const t_gtime& epoch)
{
  gtrace("t_askn::getZWD");   
   
  return 0.0;   // not implemeted yet (WV pressure not available)
}

// ---------  BACHA TADY JE TO V degrees !!!!!
double t_askn::getZWD(const t_gtriple& Ell,
		      double  T,             // temperature [K]
		      double  E,             // water vapour [hPa]
     		      double  dT,            // T lapse rate [K/km]
		      double  dE,            // E decay parameter [-]
		      double  dW,            // ZWD decay parameter [-]
		      double  ref,           // ZWD reference value [m]
                      double& ratio,         // combi ratio of E/ZWD decays [%]
		      double  tm
		     )
{
  gtrace("t_askn::getZWD");   
   
  double lat = Ell[0];                              // lat [deg]
  double hgt = Ell[2]/1000;                         // height [km]
  double gm = 9.784*( 1 - 0.00266 * cos(2*lat*D2R) 
		        - 0.00028 * hgt);           // [m/s] Askne-Nordius original

  double TmE = T * ( 1 - dT/1000 * Rd/gm/(dE + 1)); // [K] mean temperature
  double Tm  = TmE;

  double zwdO  = 0.0;
  double coef1 = 0.0, coef2 = 0.0, decay = 0.0, dO = 0.0;


  if( tm != NWM_UNKNOWN ){ Tm = tm;  }              // replace by argument
  if( dT <= 0 ) dT = 6.5;                           // [K/km]
   
  if( ratio == -1 && ref != 0.0 )
  {
    if( fabs(dE - dW) <= RATIO_LIM ){
       ratio = RATIO_EW;
    }else{
       
      if( tm != NWM_UNKNOWN ){      // use requested Tm (argument tm integrated)
         dO = 1e-06 * (K2_THA - K1_THA*Eps + K3_THA/tm) * Rd/gm/ref * E - 1;
         ratio = (dO-dW) / (dE-dW) * 100.0;
         if( ratio > 100.0 ) ratio = 100.0;
         if( ratio <   0.0 ) ratio =   0.0;
	 
      }else{                        // use requested Tm (argument tmE approximated)
         dO = 1e-06 * (K2_THA - K1_THA*Eps + K3_THA/TmE) * Rd/gm/ref * E - 1;
         ratio = (dO-dW) / (dE-dW) * 100.0;
         if( ratio > 100.0 ) ratio = 100.0;
         if( ratio <   0.0 ) ratio =   0.0;
      }
    }
  }

  // estimate ZWD (optimal decay & Tm)
  coef1 = ratio/100;
  coef2 = 1.0-coef1;
  decay = coef1*dE + coef2*dW; if( decay < 0.0 ) decay = 0.0;
  zwdO  = 1e-06 * (K2_THA - K1_THA*Eps + K3_THA/Tm) * Rd/gm/(decay + 1) * E; // [m]

#ifdef DEBUG
  double TmW   = T * ( 1 - dT/1000 * Rd/gm/(dW + 1)); // [K] mean temperature
  double zwdE  = 1e-06 * (K2_THA - K1_THA*Eps + K3_THA/TmE) * Rd/gm/(dE + 1) * E;  // [m]
  double zwdW  = 1e-06 * (K2_THA - K1_THA*Eps + K3_THA/TmW) * Rd/gm/(dW + 1) * E;  // [m]
  cerr << fixed << setprecision(3)
       << " ratio:" << setw(7) << ratio
       << " decay:" << setw(6) << decay
       << " zwdO:"  << setw(6) << zwdO
       << " zwdE:"  << setw(6) << zwdE
       << " zwdW:"  << setw(6) << zwdW
       << " zwdR:"  << setw(6) << ref
       << " zEzR:"  << setw(6) << (zwdE-zwdO)
       << " zWzR:"  << setw(6) << (zwdW-zwdO)
       << " dE:"    << setw(6) << dE
       << " dW:"    << setw(6) << dW
       << " dO:"    << setw(6) << dO
       << " Tm:"    << setw(6) << Tm
       << " e:"     << setw(6) << E
       << endl;
#endif

  return zwdO; // # [m]
}

} // namespace
