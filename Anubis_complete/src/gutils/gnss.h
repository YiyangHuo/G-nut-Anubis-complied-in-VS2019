
#ifndef GNSS_H
#define GNSS_H
 
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  
  (c) 2011-2017 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: definition of GNSS data
  Version: $ Rev: $

  2012-09-26 /JD: created

-*/

#include <map>
#include <set>
#include <vector>
#include <string>

#include "gutils/gtriple.h"
#include "gutils/gcommon.h"   // pragma

using namespace std;

// ------------------------------------------------------------------------------------------------------
// ENUMS
// ------------------------------------------------------------------------------------------------------

namespace gnut
{
   
 // GNSS systems and augmentations
 // ------------------------------
 enum GSYS { // GXX = -1,
             GPS, GAL, GLO, BDS, QZS, SBS, IRN, GNS };

 // GNSS freq Sequence ID
 // ---------------------
 enum FREQ_SEQ { FREQ_1 = 1, FREQ_2 = 2, FREQ_3 = 3, FREQ_4 = 4, FREQ_5 = 5, 
                 FREQ_X = 999 };
    
 // GNSS frequencies
 // ----------------
 enum GFRQ { // FXX = -1,
             G01 = 10, G02 = 11, G05 = 12,     // GPS
             R01 = 20, R02 = 21,               // GLONASS FDMA
             R01_CDMA = 30, R02_CDMA = 31,
             R03_CDMA = 32, R05_CDMA = 33,     // GLONASS CDMA
             E01 = 50, E05 = 51, E07 = 52,
                       E08 = 53, E06 = 54,     // Galileo
             C02 = 60, C07 = 61, C06 = 62,     // BeiDou
             J01 = 70, J02 = 71, J05 = 72, 
                                 J06 = 73,     // QZSS
             S01 = 80, S05 = 81,               // SBAS
             I05 = 90,  // I09,                // IRNSS
             LAST_GFRQ = 999
           };

 // GNSS receiver types
 // -------------------
 enum RECTYPE {  P1P2,    // receiver providing C1, P1, P2
                 C1X2,    // cross-correlation
                 C1P2     // modern receivers providing C1, P2
 };

 // Broadcast messages types
 // ------------------------
 // Broadcast messages types
 enum GNAVTYPE { FNAV, INAV, INAV_E01, INAV_E07, CNAV, NAV };

// enum GNAVTYPE {  NAV_G01, CNAV_G02, CNAV_G05,
//                  NAV_R01,
//                  NAV_C02,
//                  FNAV_E05, INAV_E01, INAV_E07, CNAV_E06, GNAV_E01, GNAV_E06,
//                  NAV_DEF
//              };

 // GNSS type/band/attr definitions
 // -------------------------------
 enum GOBSTYPE { TYPE_C = 1,   TYPE_L = 2,   TYPE_D = 3,   TYPE_S = 4, 
                 TYPE_P = 101, // only for P-code!
                 TYPE   = 999  // ""  UNKNOWN
               };
 enum GOBSBAND { BAND_1 = 1,   BAND_2 = 2,   BAND_3 = 3,   BAND_5 = 5,
                 BAND_6 = 6,   BAND_7 = 7,   BAND_8 = 8,
                 BAND_A = 101, BAND_B = 102, BAND_C = 103, BAND_D = 104,
                 BAND   = 999  // ""  UNKNOWN
               };
 enum GOBSATTR { ATTR_A, ATTR_B, ATTR_C, ATTR_D, ATTR_I, ATTR_L, ATTR_M, ATTR_N,
                 ATTR_P, ATTR_Q, ATTR_S, ATTR_W, ATTR_X, ATTR_Y, ATTR_Z,
                 ATTR_NULL,    // " " 2CHAR code
                 ATTR = 999    // ""  UNKNOWN
               };

 // GNSS observations
 // -----------------
 enum GOBS {

    // psedorange [in meters] (RINEX 3.x)
    C1A=0,   C1B, C1C,      C1I, C1L, C1M,      C1P, C1S, C1Q, C1W, C1X, C1Y, C1Z,
                  C2C, C2D, C2I, C2L, C2M,      C2P, C2S, C2Q, C2W, C2X, C2Y,
                            C3I,                          C3Q,      C3X,
    C5A,     C5B, C5C,      C5I,                          C5Q,      C5X,
    C6A,     C6B, C6C,      C6I, C6L,                C6S, C6Q,      C6X,      C6Z,
                            C7I,                          C7Q,      C7X,
                            C8I,                          C8Q,      C8X,

    // carrier phase [in whole cycles] (RINEX 3.x)
    L1A=100, L1B, L1C,      L1I, L1L, L1M, L1N, L1P, L1S, L1Q, L1W, L1X, L1Y, L1Z,
                  L2C, L2D, L2I, L2L, L2M, L2N, L2P, L2S, L2Q, L2W, L2X, L2Y,
                            L3I,                     L3Q,      L3X,
    L5A,     L5B, L5C,      L5I,                     L5Q,      L5X,
    L6A,     L6B, L6C,      L6I, L6L,           L6S, L6Q,      L6X,      L6Z,
                            L7I,                     L7Q,      L7X,
                            L8I,                     L8Q,      L8X,

    // doppler [cycles/sec] (RINEX 3.x)
    D1A=200, D1B, D1C,      D1I, D1L, D1M, D1N, D1P, D1S, D1Q, D1W, D1X, D1Y, D1Z,
                  D2C, D2D, D2I, D2L, D2M, D2N, D2P, D2S, D2Q, D2W, D2X, D2Y,
                            D3I,                     D3Q,      D3X,
    D5A,     D5B, D5C,      D5I,                     D5Q,      D5X,
    D6A,     D6B, D6C,      D6I, D6L,           D6S, D6Q,      D6X,      D6Z,
                            D7I,                     D7Q,      D7X,
                            D8I,                     D8Q,      D8X,

    // signal strength [DBHZ] (RINEX 3.x)
    S1A=300, S1B, S1C,      S1I, S1L, S1M, S1N, S1P, S1S, S1Q, S1W, S1X, S1Y, S1Z,
                  S2C, S2D, S2I, S2L, S2M, S2N, S2P, S2S, S2Q, S2W, S2X, S2Y,
                            S3I,                     S3Q,      S3X,
    S5A,     S5B, S5C,      S5I,                     S5Q,      S5X,
    S6A,     S6B, S6C,      S6I, S6L,           S6S, S6Q,      S6X,      S6Z,
                            S7I,                     S7Q,      S7X,
                            S8I,                     S8Q,      S8X,

    // special cases: v2.x or unknown tracking modes
    P1=1000, P2,         P5, C1, C2, C5, C6, C7, C8, CA, CB, CC, CD,
    L1=1100, L2,         L5, L6, L7, L8, LA, LB, LC, LD,
    D1=1200, D2,         D5, D6, D7, D8, DA, DB, DC, DD,
    S1=1300, S2,         S5, S6, S7, S8, SA, SB, SC, SD,

    X // LAST_GOBS
 };

 enum GOBS_LC {LC_UNDEF, LC_L1, LC_L2, LC_L3, LC_L4, LC_L5, LC_IF, LC_MW, LC_NL, LC_WL, LC_GF};
 
 // ------------------------------------------------------------------------------------------------------
 // TYPEDEF
 // ------------------------------------------------------------------------------------------------------
 typedef vector< GOBSATTR >          t_vec_attr;
 typedef vector< GOBSBAND >          t_vec_band;
 typedef vector< GFRQ     >          t_vec_freq;

 typedef map< GOBSTYPE, t_vec_attr > t_map_attr;
 typedef map< GOBSBAND, t_map_attr > t_map_type;

 typedef map< GSYS, set<string> >    t_map_sats;
 typedef map< GSYS, set<string> >    t_map_gnav;
 typedef map< GSYS, t_map_type  >    t_map_gnss;
 typedef map< GSYS, t_vec_band  >    t_map_band;
 typedef map< GSYS, t_vec_freq  >    t_map_freq;
   
 typedef map< GOBSBAND, t_gtriple >  t_map_pcos;  // triple: ATX  NORTH / EAST / UP
 typedef map< GSYS,    t_map_pcos >  t_map_offs;

 // ------------------------------------------------------------------------------------------------------
 // GLOBAL FUNCTIONS
 // ------------------------------------------------------------------------------------------------------
 GOBSATTR str2gobsattr( string s );           // get GOBSATTR enum from gobs string
 GOBSBAND str2gobsband( string s );           // get GOBSBAND enum from gobs string
 GOBSTYPE str2gobstype( string s );           // get GOBSTYPE enum from gobs string
 GNAVTYPE str2gnavtype( string s );           // get GNAVTYPE enum from gobs string

// GOBSATTR gobs2gobsattr( GOBS o );            // get GOBSATTR enum from gobs string
// GOBSBAND gobs2gobsband( GOBS o );            // get GOBSBAND enum from gobs string
// GOBSTYPE gobs2gobstype( GOBS o );            // get GOBSTYPE enum from gobs string
   
 GOBSATTR char2gobsattr( char c );            // get GOBSATTR enum from char
 GOBSBAND char2gobsband( char c );            // get GOBSBAND enum from char
 GOBSBAND  int2gobsband( int  c );            // get GOBSBAND enum from char   
 GOBSTYPE char2gobstype( char c );            // get GOBSTYPE enum from char
   
 string gobsattr2str( GOBSATTR e );           // get string enum from GOBSATTR 
 string gobsband2str( GOBSBAND e );           // get string enum from GOBSBAND
 string gobstype2str( GOBSTYPE e );           // get string enum from GOBSTYPE 
 string gnavtype2str( GNAVTYPE e );           // get string enum from GNAVTYPE

 string gobs2str( GOBS );                     // get string from GOBS enum
 GOBS str2gobs( string s );                   // get GOBS enum from string
 GOBS tba2gobs(GOBSTYPE t, GOBSBAND b, GOBSATTR a); // get GOBS from type, band, and attribute
 
 int gobs2band( GOBS o );                     // get band from GOBS enum

 GOBS pha2snr( GOBS o );                      // get GOBS enum (pha->snr)
 bool gobs_code( GOBS o );                    // get true for code obs
 bool gobs_phase( GOBS o );                   // get true for phase obs

 t_map_sats GNSS_SATS();                      // static map of default GNSS satellites                       
 t_map_gnav GNSS_GNAV();                      // static map of default GNSS navigation types
 t_map_freq GNSS_FREQ_PRIORITY();	            // static map of default GNSS freq priorities
 t_map_band GNSS_BAND_PRIORITY(); 	          // static map of default GNSS band priorities
 t_map_gnss GNSS_DATA_PRIORITY(); 	          // static map of default GNSS data types/bands/attrs priorities

 t_map_band GNSS_BAND_SORTED(); 	            // static map of sorted GNSS band w.r.t. wavelength
 vector<GOBSBAND> sort_band(GSYS gs, set<GOBSBAND>& bands);        // sort set of bands w.r.t. wavelength
	 
 t_map_offs GNSS_PCO_OFFSETS(); 	      // static map of default GNSS PCO offsets

 set<GSYS> GNSS_SUPPORTED();                  // supported GNSS
   
// ------------------------------------------------------------------------------------------------------
// STATIC MAPS
 // ------------------------------------------------------------------------------------------------------

/*
 // c++11 (problem with Debian 5.0)
 // -------------------------------
 static t_map_gnss GNSS_DATA_PRIORITY = {
  {GPS, { {BAND_1,{ {TYPE_C, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M,         ATTR_NULL }},
  	            {TYPE_L, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL }},
 	            {TYPE_D, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL }},
	            {TYPE_S, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_N, ATTR_NULL }},
                    {TYPE_P, { ATTR_NULL }}
          }},
          {BAND_2,{ {TYPE_C, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D,         ATTR_NULL }},
	            {TYPE_L, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D, ATTR_N, ATTR_NULL }},
	            {TYPE_D, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D, ATTR_N, ATTR_NULL }},
	            {TYPE_S, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_P, ATTR_W, ATTR_Y, ATTR_M, ATTR_D, ATTR_N, ATTR_NULL }},
                    {TYPE_P, { ATTR_NULL }}
          }},
          {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
	  }},
          {BAND_A,{ {TYPE_C, { ATTR_NULL }},
	            {TYPE_L, { ATTR_NULL }},
	            {TYPE_D, { ATTR_NULL }},
	            {TYPE_S, { ATTR_NULL }}
	  }},
          {BAND_B,{ {TYPE_C, { ATTR_NULL }},
	            {TYPE_L, { ATTR_NULL }},
	            {TYPE_D, { ATTR_NULL }},
	            {TYPE_S, { ATTR_NULL }}
	  }} 
  }},

  {GLO, { {BAND_1,{ {TYPE_C, { ATTR_C, ATTR_P, ATTR_NULL }},
	            {TYPE_L, { ATTR_C, ATTR_P, ATTR_NULL }},
	            {TYPE_D, { ATTR_C, ATTR_P, ATTR_NULL }},
	            {TYPE_S, { ATTR_C, ATTR_P, ATTR_NULL }}
                    {TYPE_P, { ATTR_NULL }}
          }},
          {BAND_2,{ {TYPE_C, { ATTR_C, ATTR_P, ATTR_NULL }},
	            {TYPE_L, { ATTR_C, ATTR_P, ATTR_NULL }},
	            {TYPE_D, { ATTR_C, ATTR_P, ATTR_NULL }},
	            {TYPE_S, { ATTR_C, ATTR_P, ATTR_NULL }}
                    {TYPE_P, { ATTR_NULL }}                   
          }},
          {BAND_3,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}  
	  }},
          {BAND_A,{ {TYPE_C, { ATTR_NULL }},
	            {TYPE_L, { ATTR_NULL }},
	            {TYPE_D, { ATTR_NULL }},
	            {TYPE_S, { ATTR_NULL }}
	  }},
          {BAND_B,{ {TYPE_C, { ATTR_NULL }},
	            {TYPE_L, { ATTR_NULL }},
	            {TYPE_D, { ATTR_NULL }},
	            {TYPE_S, { ATTR_NULL }}
	  }} 
  }},
   
  {GAL, { {BAND_1,{ {TYPE_C, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
	            {TYPE_L, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
	            {TYPE_D, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
	            {TYPE_S, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }}
          }},
          {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
	  }},
          {BAND_6,{ {TYPE_C, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
	            {TYPE_L, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
	            {TYPE_D, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }},
	            {TYPE_S, { ATTR_A, ATTR_B, ATTR_C, ATTR_X, ATTR_Z, ATTR_NULL }}
	  }},
          {BAND_7,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
	  }}, 
          {BAND_8,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
	  }},
          {BAND_A,{ {TYPE_C, { ATTR_NULL }},
	            {TYPE_L, { ATTR_NULL }},
	            {TYPE_D, { ATTR_NULL }},
	            {TYPE_S, { ATTR_NULL }}
	  }},
          {BAND_B,{ {TYPE_C, { ATTR_NULL }},
	            {TYPE_L, { ATTR_NULL }},
	            {TYPE_D, { ATTR_NULL }},
	            {TYPE_S, { ATTR_NULL }}
	  }}
  }},

  {BDS, { {BAND_1,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}, // WILL BE DELETED AFTER CONVERSION IN GCODERS
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}, // WILL BE DELETED AFTER CONVERSION IN GCODERS
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}, // WILL BE DELETED AFTER CONVERSION IN GCODERS
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}  // WILL BE DELETED AFTER CONVERSION IN GCODERS
          }},
          {BAND_2,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }} 
          }},
          {BAND_6,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
	  }},
          {BAND_7,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X, ATTR_NULL }}
	  }}
  }},
   
  {SBS, { {BAND_1,{ {TYPE_C, { ATTR_C }},
	            {TYPE_L, { ATTR_C }},
	            {TYPE_D, { ATTR_C }},
	            {TYPE_S, { ATTR_C }},
                    {TYPE_P, { ATTR_NULL }}
          }},
          {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X }}
	  }}
  }},

  {QZS, { {BAND_1,{ {TYPE_C, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }},
	            {TYPE_L, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }},
	            {TYPE_D, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }},
	            {TYPE_S, { ATTR_C, ATTR_S, ATTR_L, ATTR_X, ATTR_Z }}
          }},
          {BAND_2,{ {TYPE_C, { ATTR_S, ATTR_L, ATTR_X }},
	            {TYPE_L, { ATTR_S, ATTR_L, ATTR_X }},
	            {TYPE_D, { ATTR_S, ATTR_L, ATTR_X }},
	            {TYPE_S, { ATTR_S, ATTR_L, ATTR_X }}
	  }},
          {BAND_5,{ {TYPE_C, { ATTR_I, ATTR_Q, ATTR_X }},
	            {TYPE_L, { ATTR_I, ATTR_Q, ATTR_X }},
	            {TYPE_D, { ATTR_I, ATTR_Q, ATTR_X }},
	            {TYPE_S, { ATTR_I, ATTR_Q, ATTR_X }}
	  }},
          {BAND_6,{ {TYPE_C, { ATTR_S, ATTR_L, ATTR_X }},
	            {TYPE_L, { ATTR_S, ATTR_L, ATTR_X }},
	            {TYPE_D, { ATTR_S, ATTR_L, ATTR_X }},
	            {TYPE_S, { ATTR_S, ATTR_L, ATTR_X }}
	  }}
  }}
 };

// c++11 (problem with Debian 5.0)
// -------------------------------
static map<GSYS, set<string> > GNSS_SATS =
{
  { GPS, { "G01","G02","G03","G04","G05","G06","G07","G08","G09","G10",
           "G11","G12","G13","G14","G15","G16","G17","G18","G19","G20",
           "G21","G22","G23","G24","G25","G26","G27","G28","G29","G30",
           "G31","G32" } },
   
  { GLO, { "R01","R02","R03","R04","R05","R06","R07","R08","R09","R10",
           "R11","R12","R13","R14","R15","R16","R17","R18","R19","R20",
           "R21","R22","R23","R24" } },

  { GAL, { "E01","E02","E03","E04","E05","E06","E07","E08","E09","E10",
           "E11","E12","E13","E14","E15","E16","E17","E18","E19","E20",
           "E21","E22","E23","E24","E25","E26","E27","E28","E29","E30" } },

  { BDS, { "C01","C02","C03","C04","C05","C06","C07","C08","C09","C10",
           "C11","C12","C13","C14","C15","C16","C17","C18","C19","C20",
           "C21","C22","C23","C24","C25","C26","C27","C28","C29","C30" } },
   
  { SBS, { "S20","S24","S25","S26","S27","S28","S29",
           "S33","S35","S36","S37","S38",
           "S40",
           "S83" } },

  { QZS, { "J01" } }
};
*/

} // namespace

#endif // GOBS_H

