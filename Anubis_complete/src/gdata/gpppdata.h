/* ----------------------------------------------------------------------
(c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
Research Institute of Geodesy, Topography and Cartography
Ondrejov 244, 251 65, Czech Republic

2018-01-01 /LEWEN: created
		  Purpose: Temporary file for DCB correction
-*/

#ifndef GPPPDATA_H
#define GPPPDATA_H

#include "string.h"

using namespace std;

namespace gnut {
	class t_gbiasDCB {
	public:
		t_gbiasDCB() {};
		virtual ~t_gbiasDCB() {};

		static double str2num(const char *s, int i, int n);
		static int	  readDCB(const char * file);
		static double dcbcorr(int sat, int type);

	private:
		static double cbias[32][4];	 /* code bias (0:p1-p2,1:p1-c1,2:p2-c2) (m):only for GPS currently*/
//		static double cbiasE[32][3];	 /* code bias (0:p1-p2,1:p1-c1,2:p2-c2) (m):only for GPS currently*/     
	};
}

#endif
