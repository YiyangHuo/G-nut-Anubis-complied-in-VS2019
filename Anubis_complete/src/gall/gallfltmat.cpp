
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

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "gall/gallfltmat.h"
#include "gutils/gmatrixconv.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut {  

// constructor
// ----------
t_gallfltmat::t_gallfltmat()
{
  clear_addAmb();
}


// destructor
// ----------
t_gallfltmat::~t_gallfltmat()
{
  data.clear();
  clear_addAmb();   
}

void t_gallfltmat::add(const t_gfltmat& mat)
{
   data.push_back(mat);
}


// sync par1 and par2 
// only common parameter remains
// -----------------------------
int t_gallfltmat::common(t_gallpar& par1, SymmetricMatrix& Q1, t_gallpar& par2, SymmetricMatrix& Q2, DiagonalMatrix* Noise)
{

   for (unsigned int i = 0; i < par1.parNumber(); i++){
      int idx = par2.getParam(par1[i].site, par1[i].parType, par1[i].prn);
      if (idx >= 0) continue;
      else {	 
//cout << "common: Budu mazat par1 " << par1[i].str_type() << " " << par1[i].prn << endl;
	 int index = par1[i].index;
       	 Matrix_remRC(Q1, index, index);
         par1.delParam(i);
	 par1.decIndex(index);
	 i--;
      }
   }

   for (unsigned int i = 0; i < par2.parNumber(); i++){
      int idx = par1.getParam(par2[i].site, par2[i].parType, par2[i].prn);
      if (idx >= 0) continue;
      else {
//cout << "common: Budu mazat par2 " << par2[i].str_type() << " " << par2[i].prn << endl;	 
	 int index = par2[i].index;
       	 Matrix_remRC(Q2, index, index);
	 if (Noise) Matrix_remRC(*Noise, index);
         par2.delParam(i);	 
	 par2.decIndex(index);
	 i--;
      }
   }   

   return 1;
}

// sync smooth vector according to par1
// -------------------------------
int t_gallfltmat::syncSMT(t_gallpar& Xu, SymmetricMatrix& Qu, t_gallpar& Xsm, SymmetricMatrix& Qsm, DiagonalMatrix* Noise)
{
   
   // removing from Xsm
   for (unsigned int i = 0; i < Xsm.parNumber(); i++){
      int idx = Xu.getParam(Xsm[i].site, Xsm[i].parType, Xsm[i].prn);
      if (idx >= 0) continue;
      else {
//cout << "syncSMT: Budu mazat par " << Xsm[i].str_type() << " " << Xsm[i].prn << endl;
	 int index = Xsm[i].index;
       	 Matrix_remRC(Qsm, index, index);
         Xsm.delParam(i);
	 Xsm.decIndex(index);
	 i--;
      }
   }
   
   // adding missing par to Xsm
   for (unsigned int i = 0; i < Xu.parNumber(); i++){ 
      int ii = Xu.getParam(i+1);
      t_gpar parTMP = Xu[ii];
      int idx = Xsm.getParam(parTMP.site, parTMP.parType, parTMP.prn);
      if (idx >= 0) continue;
      else {
//cout << "syncSMT: Budu pridavat par " << parTMP.str_type() << " " << parTMP.prn << endl;
	 // introduce new  parameter (Xu[i] -> Xsm)
	 Xsm.incIndex(parTMP.index);
	 Xsm.addParam(parTMP);
	 Matrix_addRC(Qsm, parTMP.index, parTMP.index);
/*	 
	 // copy data from Qu into expanded Qsm   
	 for (int i=1; i<=Qsm.Nrows(); i++){       
	    for (int j=1; j<=Qsm.Ncols(); j++){
	       if ( i == parTMP.index || j == parTMP.index ){
		  int iS = Xsm.getParam(i);
		  int jS = Xsm.getParam(j);
		  int iU = Xu.getParam(Xsm[iS].site, Xsm[iS].parType, Xsm[iS].prn);
		  int jU = Xu.getParam(Xsm[jS].site, Xsm[jS].parType, Xsm[jS].prn);
		  int r = Xu[iU].index;
		  int c = Xu[jU].index;	     
		  Qsm(i,j) = Qu(r,c);
	       }
	    }
	 }
*/	 
	 Qsm(parTMP.index, parTMP.index) = Qu(parTMP.index, parTMP.index);
	 if(parTMP.parType == t_gpar::AMB_IF) addAmb.insert(parTMP.index);
      }
   }

   return 1;
}

// checking reset ambiguities
// -----------------------
   int t_gallfltmat::checkSlp(t_gallpar& Xu, SymmetricMatrix& Qu, t_gallpar& Xp, SymmetricMatrix& Qp, t_gallpar& Xsm, SymmetricMatrix& Qsm, set<string>& slips)
{
   for (unsigned int i = 0; i < Xsm.parNumber(); i++){
      string prn = Xsm[i].prn;
      string site = Xsm[i].site;
      int idxu = Xu.getParam(site, t_gpar::AMB_IF, prn);
      int idxp = Xp.getParam(site, t_gpar::AMB_IF, prn);
      if(slips.find(prn) != slips.end()){
         Xp[idxp].value( Xu[idxu].value() );
         Xsm[i].value( Xu[idxu].value() );	 
	 Matrix_remRC(Qsm, Xp[idxp].index, Xp[idxp].index);
	 Matrix_addRC(Qsm, Xp[idxp].index, Xp[idxp].index);
	 Qsm(Xsm[i].index, Xsm[i].index) = Qu(Xu[idxu].index, Xu[idxu].index);
      }
   }         
   

   return 1;
}

// checking outliers
// -----------------------   
int t_gallfltmat::checkOutl(t_gallpar& Xp, SymmetricMatrix& Qp, t_gallpar& Xsm, SymmetricMatrix& Qu)
{
   t_gallpar diff = Xsm - Xp;
   for (unsigned int i = 0; i < diff.parNumber(); i++){
      if(diff[i].parType == t_gpar::AMB_IF && fabs(diff[i].value()) > 2){
	 int index = diff[i].index;
	 addAmb.insert(index);
      }
   }
   return 1;
}
   

// swap rows and cols for same order
// only common pars have to be included
// return: both X1 and X2 have the sama position in cov matrix
int t_gallfltmat::reorder(t_gallpar& X1, SymmetricMatrix& Q1, t_gallpar& X2, SymmetricMatrix& Q2)
{
   if ( X1.parNumber() != X2.parNumber() ) {
      cerr << "t_gallfltmat:reorder - not synchronized!";
      return -1;
   }   
   
   for (unsigned int i = 0; i < X1.parNumber(); i++){
      int idx = X2.getParam(X1[i].site, X1[i].parType, X1[i].prn);
      if (idx >= 0){
	 if (X1[i].index != X2[idx].index){
//cout << "reorder: Prehazuju " << X1[i].index << " za " << X2[idx].index << endl;
//cout << "Before: " << fixed << setprecision(6) << Q1 << endl;
	    Matrix_swap(Q1, X1[i].index, X2[idx].index);
//cout << "After: " << fixed << setprecision(6) << Q1 << endl;
	    int in = X1.getParam(X2[idx].index);
	    X1[in].index = X1[i].index;
	    X1[i].index = X2[idx].index;	    
	 }
      }else{
	 cerr << "t_gallfltmat::reorder - not common params!" << endl;
	 return -1;
      }
   }
   return 1;
}

int t_gallfltmat::reorder2(t_gallpar& X1, SymmetricMatrix& Q1, t_gallpar& X2, SymmetricMatrix& Q2)
{
   if ( X1.parNumber() != X2.parNumber() ) {
      cerr << "t_gallfltmat:reorder - not synchronized!";
      return -1;
   }   

   SymmetricMatrix Q_sav = Q1;
   
    for (size_t i=1; i<=Q1.Nrows(); i++){       
       for (size_t j=1; j<=Q1.Ncols(); j++){	  	  
	     int i1 = X1.getParam(i);
	     int j1 = X1.getParam(j);	  	  
	     int i2 = X2.getParam(X1[i1].site, X1[i1].parType, X1[i1].prn);
	     int j2 = X2.getParam(X1[j1].site, X1[j1].parType, X1[j1].prn);	  
	     int r = X2[i2].index;
	     int c = X2[j2].index;
	     Q1(i,j) = Q_sav(r,c);
       }
    }
   
   for (unsigned int i = 0; i < X1.parNumber(); i++){
      int idx = X2.getParam(X1[i].site, X1[i].parType, X1[i].prn);
      X1[i].index = X2[idx].index;
   }
   
   
   return 1;
}


// clear data
// -------------------------
void t_gallfltmat::clear()
{
   data.clear();
}

// clear data
// -------------------------
void t_gallfltmat::clear_addAmb()
{   
   addAmb.clear();
}

} // namespace
