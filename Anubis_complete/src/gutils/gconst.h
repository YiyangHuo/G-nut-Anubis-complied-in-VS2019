
#ifndef CONST_H
#define CONST_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: defines constants
  Version: $ Rev: $

  2011-01-10 /JD: created

-*/

#include <map>
#include <string>
#include <vector>

using namespace std;

namespace gnut {
#define G_PI        3.1415926535897932    // pi
#define D2R         (G_PI/180.0)          // deg to rad
#define R2D         (180.0/G_PI)          // rad to deg
#define SEC2RAD     (D2R/3600.0)          // sec to rad
#define RAD2SEC     (R2D*3600.0)          // rad to sec
//#define RAD2SEC     (180.0*3600.0/G_PI)   // rad to sec
//#define SEC2RAD     (G_PI/(180.0*3600.0)) // sec to rad
#define RHO_SEC     3600.0*D2R
#define CLIGHT      2.99792458e+8         // speed of light [m/s]
#define OMEGA       7292115.1467e-11
#define Aell        6378137.000           // semi major axes [m]
#define Finv        298.2572236           // inverse flattening [-]
#define MJD_J2000   51544.5               // J2000 MJD date
#define R_SPHERE    6371000.000           // [m] approx. sphere radius
#define POST_SHADOW 1800.0                // satellite post shadow period [s]
   
// GPS   
#define A_WGS       6378137.000           // [m] WGS84 semi-major axis
#define B_WGS       6356752.300           // [m] WGS84 semi-minor axis
#define E_WGS       0.081819              // [-] WGS84 eccentricity
#define F_WGS       0.003352811           // [-] WGS84 flatenning
#define MUDOT_GPS   (0.00836*D2R)         // avarage GPS satellite angular velocity [rad]
#define EPS0_GPS    (13.5*D2R)            // maximal GPS satellites crossing angle [deg]
	 
// GLONASS
#define GM_PZ90     398.60044e12          // [] PZ90 earth's graviational constant
#define Aell_PZ90   6378136.000           // [m] PZ90 semi-major axes
#define C20_PZ90   -1082.62575e-6         // [] PZ90
#define MUDOT_GLO   (0.00888*D2R)         // avarage GLO satellite angular velocity [rad]   
#define EPS0_GLO    (14.2*D2R)            // maximal GLO satellites crossing angle [deg]	 
#define Finv_GLO    298.25784             // inverse flattening [-]
#define OMGE_DOT_GLO 7.292115e-5          // Mean angular velocity of the Earth [rad/sec]
   
// Galileo
#define GM_GAL       3.986004418e14      // Geocentric gravitational constant [m^3/s^2]
#define OMGE_DOT_GAL 7.2921151467e-5     // Mean angular velocity of the Earth [rad/sec]
   
// BeiDou
#define GM_CGCS      3.986004418e14      // []  CGCS2000 earth's graviational constant
#define OMGE_DOT_BDS 7.2921150e-5        // BDS value of the earth's rotation rate [rad/sec]// 
#define Aell_CGCS    6378137.000         // [m] CGCS2000 semi-major axes
   
// METEO
#define G_RATIO     0.003449787           // [-] gravity ratio
#define G_EQUA      9.7803253359          // [m/s2] equatorial gravity
#define G_POLE      9.8321849378          // [m/s2] polar gravity
#define G_WMO       9.80665               // [m/s2] gravity acceleration (WMO constant for lat=45)
#define KS_SOM      0.001931853           // [-] Somigliana's constant
#define Md          28.9644               // [g/mol] molar mass (mean molecular weight) of dry air
#define Mw          18.01528              // [g/mol] molar mass (mean molecular weight) of water
#define Ru          8.3144621             // [J/mol/K] universal gas constant
#define Rd          287.058               // [J/kg/K] spec. gas constant for Dry Air (Rd = 287.058)  Ru/(Md/1000)
#define Rv          461.495               // [J/kg/K] spec. gas constant for Wat Vap (Rv = 461.495)  Ru/(Mw/1000)
#define Eps         0.62198               // [-] molecular weight ratio (Mw/Md)     (Eps = 0.62198)
#define TC2TK       273.15                // conversion from deg of Celsius to Kelvins
#define TPOINT      273.16                // temperature of triple point [K], i.e. 0.01 Deg C
   
 typedef map<string, vector<double> > t_map_refr; // refractivity coefficients
 t_map_refr REFR_COEF();                          // static map of refractivity coefficients

#define K1_ESS      77.64                 // [K/hPa]   Essen and Froome (1951)
#define K2_ESS      64.68                 // [K/hPa]   Essen and Froome (1951)
#define K3_ESS      3.718e5               // [K^2/hPa] Essen and Froome (1951)
   
#define K1_SMW      77.607                // [K/hPa]   Smith and Weintraub (1953) // FULL RESOLUTION
#define K2_SMW      71.6                  // [K/hPa]   Smith and Weintraub (1953) // original paper 3-par formula!
#define K3_SMW      3.747e5               // [K^2/hPa] Smith and Weintraub (1953) // original paper 3-par formula!

#define K1_THA      77.604                // [K/hPa]   Thayer (1974)
#define K2_THA      64.79                 // [K/hPa]   Thayer (1974)
#define K3_THA      3.776e5               // [K^2/hPa] Thayer (1974)

#define K1_BEV      77.6                  // [K/hPa]   Bevis (1994)
#define K2_BEV      70.4                  // [K/hPa]   Bevis (1994)
#define K3_BEV      3.739e5               // [K^2/hPa] Bevis (1994)

#define K1_FOE      77.65                 // [K/hPa]   Foelsche (1999)
#define K2_FOE      65.99                 // [K/hPa]   Foelsche (1999)
#define K3_FOE      3.777e5               // [K^2/hPa] Foelsche (1999)

#define K1_RUE      77.6848               // [K/hPa]   Rueger (2002)
#define K2_RUE      71.2952               // [K/hPa]   Rueger (2002)
#define K3_RUE      3.75463e5             // [K^2/hPa] Rueger (2002)

#define K1          K1_THA
#define K2          K2_THA
#define K3          K3_THA
#define K22         K2 - K1*Eps

} // namespace

#endif

