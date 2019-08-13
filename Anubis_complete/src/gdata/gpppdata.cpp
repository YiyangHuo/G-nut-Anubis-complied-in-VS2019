/* ----------------------------------------------------------------------
(c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
Research Institute of Geodesy, Topography and Cartography
Ondrejov 244, 251 65, Czech Republic

2018-01-01 /LEWEN: created
		  Purpose: Temporary file for DCB correction
-*/

#include "gdata/gpppdata.h"
#include "gutils/gconst.h"
#include "gutils/gtypeconv.h"

namespace gnut {

	double t_gbiasDCB::cbias[32][4] = { 0.0 };
	/* string to number ------------------------------------------------------------
	* convert substring in string to number
	* args   : char   *s        I   string ("... nnn.nnn ...")
	*          int    i,n       I   substring position and width
	* return : converted number (0.0:error)
	*-----------------------------------------------------------------------------*/
	double t_gbiasDCB::str2num(const char *s, int i, int n)
	{
		/* double value; */
		/* char str[256], *p = str; */

		/* if (i < 0 || (int)strlen(s) < i || (int)sizeof(str) - 1 < n) return 0.0; */
		/* for (s += i; *s&&--n >= 0; s++) {*p++ = *s == 'd' || *s == 'D' ? 'E' : *s; *p = '\0';} */
		/* return sscanf(str, "%lf", &value) == 1 ? value : 0.0; */

		string str(s, 47);
		string num_str = str.substr(26, 9);
		double num = str2dbl(num_str);
		//     cout << "str " << str << endl << "num_str " << num_str << endl << "num " << num << endl;
		return num;
	}

	double t_gbiasDCB::dcbcorr(int sat, int type)
	{
		return cbias[sat - 1][type];
	}

	int t_gbiasDCB::readDCB(const char * file)
	{
		FILE *fp;
		int type = 0;
		char buff[256], code;
		int prn;

		if (!(fp = fopen(file, "r"))) {
			printf("dcb parameters file open error: %s\n", file);

			return 0;
		}
		printf("read dcb parameters file: %s\n", file);
		while (fgets(buff, sizeof(buff), fp)) {

			if (strstr(buff, "DIFFERENTIAL (P1-P2) CODE BIASES")) type = 1;
			else if (strstr(buff, "DIFFERENTIAL (P1-C1) CODE BIASES")) type = 2;
			else if (strstr(buff, "DIFFERENTIAL (P2-C2) CODE BIASES")) type = 3;
			else if (strstr(buff, "DIFFERENTIAL (E1-E5) CODE BIASES")) type = 4;

			if (!type) continue;
			if (sscanf(buff, "%c%d", &code, &prn) == 2)
			{
				if (code == 'G' || code == 'E') {
					cbias[prn - 1][type - 1] = str2num(buff, 26, 9)*1E-9*CLIGHT; /* ns -> m */
																				 //              cout << "Ctu DCB " << prn << " " << type-1 << " " << cbias[prn - 1][type - 1] << endl;
				}
			}
		}
		fclose(fp);
		return 1;
	}
}