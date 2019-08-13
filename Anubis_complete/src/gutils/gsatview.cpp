
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

#include <math.h>

#include "gall/gallnav.h"
#include "gutils/gtime.h"
#include "gutils/gtriple.h"
#include "gutils/gsatview.h"
#include "gutils/gcommon.h"

using namespace std;

namespace gnut {

// calculate satellite visibility
// ----------
int sat_view( t_map_sky_sat&   sat_view,
	      t_gtriple        xyz_rec,
 	      t_gallnav*       all_nav,
	      const t_gtime&   beg,       // should be setup reasonably !
	      const t_gtime&   end,       // should be setup reasonably !
	      double mask
	    )
{
  gtrace("sat_view");

  if( beg > end || beg == FIRST_TIME || end == LAST_TIME ||
      ! all_nav ) return -1;

  t_gtriple xyz_sat;
  t_gtriple top_sat;
  double sec = 120.0; // add seconds as a function of elevation (120s ~ 1 deg)
  double tdif = 0.0;
  double r_radi, s_radi, s_topo;

  // loop over all satellites
  t_map_sats gnss_sats = GNSS_SATS();
  t_map_sats::const_iterator  itSYS;
  set<string>::const_iterator itSAT;

  for( itSYS = gnss_sats.begin(); itSYS != gnss_sats.end(); ++itSYS )
  { 
    GSYS gs = itSYS->first;

    for( itSAT = gnss_sats[gs].begin(); itSAT != gnss_sats[gs].end(); ++itSAT )
    { 
      string prn = *itSAT;
#ifdef DEBUG
      t_gtime beg_eph = all_nav->beg_time( prn );
      t_gtime end_eph = all_nav->end_time( prn );
      cerr << prn << beg_eph.str_ymdhms("  beg: ") << end_eph.str_ymdhms("  end: ") << endl;
       
      if( beg_eph == LAST_TIME  || end_eph == FIRST_TIME )
      {
	cerr << prn << " no eph - sat skipped "
	     << " per: " << beg_eph.str_ymdhms()
	     <<    " - " << end_eph.str_ymdhms() << endl;
        continue;
      }
#endif	 
       
      t_gtime beg_lim = beg;
      t_gtime end_lim = end;

#ifdef DEBUG
      double dt_lim = end_lim - beg_lim;
      if( dt_lim > 86400 || beg_lim == FIRST_TIME ) beg_lim = beg_eph;
      if( dt_lim > 86400 || end_lim == LAST_TIME  ) end_lim = end_eph;

      cout << " sat: " << prn
	   << " per: " << beg_lim.str_ymdhms() 
	   <<    " - " << end_lim.str_ymdhms()
	   << endl;
#endif

      t_gtime tt     = beg_lim;
      t_gtime tt_sav = beg_lim;
      double ele     = 0.0;
      double ele_sav = 0.0;
      int  direction = 0;
      int  iter = 0;
      bool init = true;
      
      while( tt < end_lim + 1 && ++iter < 75 )
      {
        if( all_nav->pos( prn, tt, xyz_sat.crd_array(), NULL, NULL ) < 0 ){
	  
	  // TEMPORARY SOLUTION FOR GLONASS (hope this all function will be re-implemented soon)
	  if( gs == GLO ){ tt.add_secs(3600); ele = 1.0; continue; } // risky to skip nav gap	 
	  break; // better STOP and skip this SAT, but a problem with GLONASS when missing NAV
	}

        top_sat = xyz_sat - xyz_rec;

        r_radi = xyz_rec.crd_cvect().norm_Frobenius()/1000.0; // km
        s_radi = xyz_sat.crd_cvect().norm_Frobenius()/1000.0; // km
        s_topo = top_sat.crd_cvect().norm_Frobenius()/1000.0; // km
     
        ele = -90 + acos( ( pow(s_topo, 2.0) + pow(r_radi, 2.0) - pow(s_radi, 2.0) )
	  		   / ( 2 * r_radi * s_topo ) ) / G_PI * 180;

        if( gs == QZS || gs == SBS || t_gsys::bds_geo(prn) ){
          // initiate new start with end time!
	  if( ele > mask ) sat_view[prn][beg_lim] = end_lim;
//	  cerr << prn << " - satview: skipping geostacionary satellites QZS, SBS or BDS, ele: " << ele << endl;
	  break;
	}
	
	// gradient
	tdif = tt - tt_sav;
	if( tdif != 0.0 ) sec = tdif / fabs(ele - ele_sav) / 2.0;
	if( sec > 240.0 ) sec = 240.0;

	// horizon control
//	bool horizon_found = ( fabs(ele) < 0.01 && !init ); // not below 0.01! for 1s increase at least
	bool horizon_found = ( fabs(tt-tt_sav) < 2 && !init ); // MORE ROBUST, GLONASS Particularly!
//	bool horizon_cross = ((ele > 0) - (ele_sav > 0) && !init ); // alternative to a sign check
	bool horizon_cross = ((ele > mask) - (ele_sav > mask) && !init ); // alternative to a sign check
	 
	if( !init && !horizon_found ) //  !horizon_cross &&
	{
  	  if( ele_sav < ele ) direction = +1;
	  if( ele_sav > ele ) direction = -1;
	}

 	init = false;

#ifdef DEBUG	 
	cerr << fixed << setprecision(3)
	     << setw(3) << iter
	     << "  " << prn  
	     << "  " << tt.str_ymdhms()  
	     << "  " << tt_sav.str_ymdhms()
	     << "  " << tdif
	     << " sec: " << setw(7) << sec
	     << " ele: " << setw(7) << ele
	     << " sav: " << setw(7) << ele_sav
	     << " " << init << ":" << horizon_cross << ":" << horizon_found
	     << " " << setw(2) << direction
             << "  " << endl;
#endif

	// horizon found
	if( horizon_found ){ // && direction != 0 ){
#ifdef DEBUG
	  cerr << "  " << prn << " add H 1800.000";
	  if( direction > 0 ){  cerr << "  ascend horizon\n\n"; }
	  else{                 cerr << " descend horizon\n\n"; }
#endif
	  // FILL 
          map<t_gtime, t_gtime>::reverse_iterator itEND = sat_view[prn].rbegin();
	  if( sat_view[prn].size() > 0 )                          // multiple visibility record
	  { if( direction < 0 )     itEND->second = tt;           // add visibility end time
	    else                sat_view[prn][tt] = end_lim;      // initiate new start with end time!
	  }
	  else                                                    // initial visibility record
	  { if( direction < 0 ) sat_view[prn][beg_lim] = tt;      // add visibility end time
	    else                sat_view[prn][tt] = end_lim;      // initiate new start with end time!
	  }

	  iter   = 0;
          tt_sav = tt;
	  tt.add_secs(300);
	  ele_sav = mask + direction*0.01;
	  continue;
 	}
	  
        // horizon roughly crossing
	else if( horizon_cross ){
	  
	  if( fabs(ele_sav-mask) > fabs(ele-mask) || fabs(ele-mask) < 1.0 ){
//            cerr << "  " << prn << "  add M " << setw(8) << -sec*fabs(ele-mask) << endl;
	    tt_sav = tt;
	    tt.add_secs( (int)(-sec*fabs(ele-mask)) );
	  }else{            
//            cerr << "  " << prn << "  add P " << setw(8) << sec*fabs(ele_sav-mask) << endl;
	    if( fabs(ele_sav-mask) < 0.2 ) tt = tt_sav + sec*fabs((ele_sav-mask)/2);
	    else                           tt = tt_sav + sec*fabs((ele_sav-mask));
	  }
	}

	else{
  	  ele_sav = ele;
          tt_sav = tt;

	  double add = sec*fabs(ele-mask);
	  if( fabs(ele-mask) > 10 ) add *= 0.6;

//	  cout << prn << " add A " << setw(8) << add << endl;
          tt.add_secs((int)add);

	}
      }       
    }
  }

#ifdef DEBUG 
    t_map_sky_sat::iterator itPRN;
    for( itPRN = sat_view.begin(); itPRN != sat_view.end(); ++itPRN )
    {	 
      string prn = itPRN->first;

      t_map_sky_epo::iterator itBEG;
      for( itBEG = sat_view[prn].begin(); itBEG != sat_view[prn].end(); ++itBEG )
      {
        cerr << " PRN: " << prn 
	     << "   "    << itBEG->first.str_ymdhms()
	     << " - "    << itBEG->second.str_ymdhms()
	     << " time:" << prn << " " << itBEG->second-itBEG->first
	     << endl;
      }  
      cerr << endl;
    }
#endif

  return 1;
}

   
   
// calculate satellite visibility
// ----------
int sat_elev( t_map_sky_epo&   sat_elev,
	      t_gtriple        xyz_rec,
 	      t_gallnav*       all_nav,
	      const t_gtime&   beg,       // should be setup reasonably !
	      const t_gtime&   end,       // should be setup reasonably !
	      const  string&   prn,
	      double mask
	    )
{
  gtrace("sat_elev");

  if( beg > end || beg == FIRST_TIME || end == LAST_TIME ||
      ! all_nav ) return -1;

  t_gtriple xyz_sat;
  t_gtriple top_sat;
  double sec = 120.0; // add seconds as a function of elevation (120s ~ 1 deg)
  double tdif = 0.0;
  double r_radi, s_radi, s_topo;

  // loop over all satellites
  t_map_sats gnss_sats = GNSS_SATS();
  t_map_sats::const_iterator  itSYS;
  set<string>::const_iterator itSAT;

  GSYS gs = t_gsys::str2gsys(prn);
   
#ifdef DEBUG
      t_gtime beg_eph = all_nav->beg_time( prn );
      t_gtime end_eph = all_nav->end_time( prn );
      cerr << prn << beg_eph.str_ymdhms("  beg: ") << end_eph.str_ymdhms("  end: ") << endl;
       
      if( beg_eph == LAST_TIME  || end_eph == FIRST_TIME )
      {
	cerr << prn << " no eph - sat skipped "
	     << " per: " << beg_eph.str_ymdhms()
	     <<    " - " << end_eph.str_ymdhms() << endl;
        return -1;
      }
#endif	 
       
      t_gtime beg_lim = beg;
      t_gtime end_lim = end;

#ifdef DEBUG
      double dt_lim = end_lim - beg_lim;
      if( dt_lim > 86400 || beg_lim == FIRST_TIME ) beg_lim = beg_eph;
      if( dt_lim > 86400 || end_lim == LAST_TIME  ) end_lim = end_eph;

      cout << " sat: " << prn
	   << " per: " << beg_lim.str_ymdhms() 
	   <<    " - " << end_lim.str_ymdhms()
	   << endl;
#endif

      t_gtime tt     = beg_lim;
      t_gtime tt_sav = beg_lim;
      double ele     = 0.0;
      double ele_sav = 0.0;
      int  direction = 0;
      int  iter = 0;
      bool init = true;
      
      while( tt < end_lim + 1 && ++iter < 75 )
      {
        if( all_nav->pos( prn, tt, xyz_sat.crd_array(), NULL, NULL ) < 0 ){
	  
	  // TEMPORARY SOLUTION FOR GLONASS (hope this all function will be re-implemented soon)
	  if( gs == GLO ){ tt.add_secs(3600); ele = 1.0; continue; } // risky to skip nav gap	 
	  break; // better STOP and skip this SAT, but a problem with GLONASS when missing NAV
	}

        top_sat = xyz_sat - xyz_rec;

        r_radi = xyz_rec.crd_cvect().norm_Frobenius()/1000.0; // km
        s_radi = xyz_sat.crd_cvect().norm_Frobenius()/1000.0; // km
        s_topo = top_sat.crd_cvect().norm_Frobenius()/1000.0; // km
     
        ele = -90 + acos( ( pow(s_topo, 2.0) + pow(r_radi, 2.0) - pow(s_radi, 2.0) )
	  		   / ( 2 * r_radi * s_topo ) ) / G_PI * 180;

        if( gs == QZS || gs == SBS || t_gsys::bds_geo(prn) ){
          // initiate new start with end time!
	  if( ele > mask ) sat_elev[beg_lim] = end_lim;
//	  cerr << prn << " - satview: skipping geostacionary satellites QZS, SBS or BDS, ele: " << ele << endl;
	  break;
	}
	
	// gradient
	tdif = tt - tt_sav;
	if( tdif != 0.0 ) sec = tdif / fabs(ele - ele_sav) / 2.0;
	if( sec > 240.0 ) sec = 240.0;

	// horizon control
//	bool horizon_found = ( fabs(ele) < 0.01 && !init ); // not below 0.01! for 1s increase at least
	bool horizon_found = ( fabs(tt-tt_sav) < 2 && !init ); // MORE ROBUST, GLONASS Particularly!
//	bool horizon_cross = ((ele > 0) - (ele_sav > 0) && !init ); // alternative to a sign check
	bool horizon_cross = ((ele > mask) - (ele_sav > mask) && !init ); // alternative to a sign check
	 
	if( !init && !horizon_found ) //  !horizon_cross &&
	{
  	  if( ele_sav < ele ) direction = +1;
	  if( ele_sav > ele ) direction = -1;
	}

 	init = false;

//#ifdef DEBUG	 
	cerr << fixed << setprecision(3)
	     << setw(3) << iter
	     << "  " << prn  
	     << "  " << tt.str_ymdhms()  
	     << "  " << tt_sav.str_ymdhms()
	     << "  " << tdif
	     << " sec: " << setw(7) << sec
	     << " ele: " << setw(7) << ele
	     << " sav: " << setw(7) << ele_sav
	     << " " << init << ":" << horizon_cross << ":" << horizon_found
	     << " " << setw(2) << direction
             << "  " << endl;
//#endif

	// horizon found
	if( horizon_found ){ // && direction != 0 ){
//#ifdef DEBUG
	  cerr << "  " << prn << " add H 1800.000";
	  if( direction > 0 ){  cerr << "  ascend horizon\n\n"; }
	  else{                 cerr << " descend horizon\n\n"; }
//#endif
	  // FILL 
          map<t_gtime, t_gtime>::reverse_iterator itEND = sat_elev.rbegin();
	  if( sat_elev.size() > 0 )                          // multiple visibility record
	  { if( direction < 0 )     itEND->second = tt;           // add visibility end time
	    else                sat_elev[tt] = end_lim;      // initiate new start with end time!
	  }
	  else                                                    // initial visibility record
	  { if( direction < 0 ) sat_elev[beg_lim] = tt;      // add visibility end time
	    else                sat_elev[tt] = end_lim;      // initiate new start with end time!
	  }

	  iter   = 0;
          tt_sav = tt;
	  tt.add_secs(300);
	  ele_sav = mask + direction*0.01;
	  continue;
 	}
	  
        // horizon roughly crossing
	else if( horizon_cross ){
	  
	  if( fabs(ele_sav-mask) > fabs(ele-mask) || fabs(ele-mask) < 1.0 ){
//            cerr << "  " << prn << "  add M " << setw(8) << -sec*fabs(ele-mask) << endl;
	    tt_sav = tt;
	    tt.add_secs( (int)(-sec*fabs(ele-mask)) );
	  }else{            
//            cerr << "  " << prn << "  add P " << setw(8) << sec*fabs(ele_sav-mask) << endl;
	    if( fabs(ele_sav-mask) < 0.2 ) tt = tt_sav + sec*fabs((ele_sav-mask)/2);
	    else                           tt = tt_sav + sec*fabs((ele_sav-mask));
	  }
	}

	else{
  	  ele_sav = ele;
          tt_sav = tt;

	  double add = sec*fabs(ele-mask);
	  if( fabs(ele-mask) > 10 ) add *= 0.6;

//	  cout << prn << " add A " << setw(8) << add << endl;
          tt.add_secs((int)add);

	}
      }

//#ifdef DEBUG 
  t_map_sky_epo::iterator itBEG;
  for( itBEG = sat_elev.begin(); itBEG != sat_elev.end(); ++itBEG )
  {
    cerr << " PRN: " << prn
	     << "   "    << itBEG->first.str_ymdhms()
	     << " - "    << itBEG->second.str_ymdhms()
	     << " time:" << prn << " " << itBEG->second-itBEG->first
	     << endl;
  }
  cerr << endl;
//#endif

  return 1;
}

} // namespace
