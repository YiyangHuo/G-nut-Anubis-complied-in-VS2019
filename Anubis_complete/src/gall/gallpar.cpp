
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
#include <iomanip>
#include <cmath>

#include "gall/gallpar.h"

using namespace std;

namespace gnut {  

// Add t_gpar to t_gallpar
// Parameter is stored at the end of the vector
// -----------------------------------------------
void t_gallpar::addParam(t_gpar newPar)
{
   gtrace("t_gallpar::addParam");
   
   this->_vParam.push_back(newPar);
}

// Delete paremeter
// Parameter is deleted according to index value
// ----------------------------------------------------
void t_gallpar::delParam(int i)
{
   gtrace("t_gallpar::delParam");   
   
   _vParam.erase( _vParam.begin() + i );
      
   /* for(iter=_vParam.begin(); iter!=_vParam.end(); iter++) */
   /*   { */
   /*    if (iter->index == i+1){ */
   /*         _vParam.erase(iter); */
   /*         break; */
   /*    }   */
   /*    else continue; */
   /*   } */
}

// get position of item according: station name, par type, PRN, begin time, end time
// -----------------------------------------------------
int t_gallpar::getParam(string mark, t_gpar::t_type type, string prn, 
                t_gtime beg, t_gtime end) const
{
  gtrace("t_gallpar::getParam");
   
//   for ( iter = _vParam.begin(); iter != _vParam.end(); iter++ )
  for ( unsigned int i=0; i < _vParam.size(); i++)
    {
      if (_vParam[i].site.compare(mark) == 0 &&
          _vParam[i].parType == type         &&
          (_vParam[i].beg == beg || _vParam[i].beg > beg ) &&
          (_vParam[i].end == end || _vParam[i].end < end ) )
        {
          if(prn == ""){
            return i;
            break;
          }else if(_vParam[i].prn.compare(prn) == 0){
            return i;
            break;        
          } else continue;
        }   
    }
  return -1;
}

// get position of item according: index (position in covariance matrix)
// -----------------------------------------------------
int t_gallpar::getParam(int index)
{
   gtrace("t_gallpar::getParam(int)");   
   
   for ( unsigned int i=0; i<=_vParam.size()-1; i++)
     {
   if (_vParam[i].index == index)
     {
        return i;
        break;
     }   
     } 
   return -1;
}

// Reindexing parametres.
// New indexes are reordered form 1 to n
// -----------------------------------------------
void t_gallpar::reIndex()
{
   gtrace("t_gallpar::reIndex");   
   
   int index_new = 1;
   for (unsigned int iPar=0; iPar<=_vParam.size()-1; iPar++)
     {
   _vParam[iPar].index = index_new;
   index_new++;
     }
}

// Reindexing parametres.
// All indexes larger than "i" is decresed by 1
// -----------------------------------------------
void t_gallpar::decIndex(int i)
{
   gtrace("t_gallpar::decIndex(int)");   
   
   for (unsigned int iPar=0; iPar<=_vParam.size()-1; iPar++)
     {
   if (_vParam[iPar].index > i) _vParam[iPar].index -= 1;
     }   
}

// Reindexing parametres.
// All indexes larger than "i" is incresed by 1
// -----------------------------------------------
void t_gallpar::incIndex(int i)
{
   gtrace("t_gallpar::incIndex(int)");   
   
   for (unsigned int iPar=0; iPar<=_vParam.size()-1; iPar++)
     {
   if (_vParam[iPar].index >= i) _vParam[iPar].index += 1;
     }   
}

// Get number of parametres
// ---------------------------------------------
unsigned int t_gallpar::parNumber() const
{
   gtrace("t_gallpar::parNumber");
   
   return _vParam.size();
}
   
// Get number of ambiguities
// ---------------------------------------------
unsigned int t_gallpar::ambNumber() const
{
   gtrace("t_gallpar::ambNumber");
   
   int i = 0;
   for (unsigned int iPar=0; iPar<=_vParam.size()-1; iPar++){
      if(_vParam[iPar].parType == t_gpar::AMB_IF) i++;
   }
   
   return i;
}   

// Strore Coordinates parametres to t_triple crd
// ------------------------------------------------

int t_gallpar::getCrdParam(string station, t_gtriple& crd, t_gtime Tbeg, t_gtime Tend) const
{      
  gtrace("t_gallpar::getCrdParam");   
   
  int found = 0;
  vector<t_gpar>::const_iterator iter;
  for(iter=_vParam.begin(); iter!=_vParam.end(); iter++){
    if (iter->parType == t_gpar::CRD_X  &&  iter->site.compare(station) == 0 &&
        iter->beg     == Tbeg           &&  iter->end == Tend) {
      crd.set(0, iter->value());
      found++;
    }else if (iter->parType == t_gpar::CRD_Y  &&  iter->site.compare(station) == 0 && 
              iter->beg     == Tbeg           &&  iter->end == Tend) {
      crd.set(1, iter->value());
      found++;
    }else if (iter->parType == t_gpar::CRD_Z  &&  iter->site.compare(station) == 0 && 
              iter->beg     == Tbeg           &&  iter->end == Tend) {
      crd.set(2, iter->value());
      found++;
    }
  }
  if (found == 3)  return  1;    // all three crd were found
  if (found == 1)  return -1;    // just one crd was found
  if (found == 2)  return -2;    // just two crd was found

  if (found == 0) return  -3;    // crd not found
  
  return -1;
}


// Repair ambiguity parameters due to receiver clock jump
// Argumetn is clk jump value
// ---------------------------------------------------
/* void t_gallpar::repairClkJump(double jump) */
/* { */
/*    for( iter=_vParam.begin(); iter!=_vParam.end(); iter++ ) */
/*      { */
/*    if (iter->parType == t_gpar::AMB_IF) */
/*      { */
/*         iter->value += jump*CLIGHT; */
/*      }    */
/*      }    */
/* } */

// get single t_gpar element from container
// ----------------------------------------------
t_gpar& t_gallpar::operator[](const size_t idx)
{
  return _vParam[idx];
}

// Operator -
// ----------------------------------------
t_gallpar t_gallpar::operator-(t_gallpar& gallpar)
{   
   t_gallpar diff;
   if ( this->parNumber() != gallpar.parNumber() ){
      cerr << "t_gallpar::operator-: Incompatible dimension ("
      << this->parNumber() << ", " << gallpar.parNumber() << ")" 
      << endl;
      return diff;
   }
      
   diff = (*this);
   vector<t_gpar>::const_iterator iter;
   for (iter = _vParam.begin(); iter != _vParam.end(); iter++) {     
      int i = gallpar.getParam(iter->site, iter->parType, iter->prn, iter->beg ,iter->end);
      if (i >= 0) {
    diff[i] = (*iter) - gallpar[i];
      }       
//      else cerr << "NENASEL JSEM " << iter->str_type() << endl;
    } 
   return diff;
}

// Operator +
// ----------------------------------------
t_gallpar t_gallpar::operator+(t_gallpar& gallpar)
{   
   t_gallpar diff;
     
   if ( this->parNumber() != gallpar.parNumber() ){
      cerr << "t_gallpar::operator+: Incopatible dimension ("  
      << this->parNumber() << ", " << gallpar.parNumber() << ")" 
      << endl;
      return diff;
   }
   
     diff = (*this);
     vector<t_gpar>::const_iterator iter;
     for (iter = _vParam.begin(); iter != _vParam.end(); iter++) {     
       int i = gallpar.getParam(iter->site, iter->parType, iter->prn, iter->beg ,iter->end);
       if (i >= 0) {
     diff[i] = (*iter) + gallpar[i];
       }
   
     }
   return diff;
}

//Doplnil Gabo
void t_gallpar::delAllParam()
{
   gtrace("t_gallpar::delAllParam");
   
   _vParam.clear();
}

// Delete all ambiguity params
// ------------------------------------
vector<int> t_gallpar::delAmb()
{
   gtrace("t_gallpar::addAmb");   
   
   vector<int> ind;
   vector<t_gpar>::iterator iter;
   iter = _vParam.begin();
   while (iter != _vParam.end())
   {
      if (iter->parType == t_gpar::AMB_IF) {
    ind.push_back(iter->index);
    iter = _vParam.erase(iter);
      }else iter++;
   }
   return ind;
}

// Reset all parems value
// -------------------------------
void t_gallpar::resetAllParam()
{
   gtrace("t_gallpar::resetAllParam");   
   
   vector<t_gpar>::iterator iter;
   for (iter = _vParam.begin(); iter != _vParam.end(); iter++) iter->value(0);
}


// set site name for all pars in gallpar
// -----------------------------------------
void t_gallpar::setSite(string site)
{
   gtrace("t_gallpar::setSite");   
   
   vector<t_gpar>::iterator iter;
   for (iter = _vParam.begin(); iter != _vParam.end(); iter++)  iter->site = site;
   
}

// Multiple matrix and gallpar
// -----------------
int t_gallpar::mult(const Matrix& K, t_gallpar& mult, t_gallpar& res)
{   
  gtrace("t_gallpar::mult");   
   
  size_t parN = mult.parNumber();
  if ( parN != K.Ncols()) {
    cerr << "t_gallpar::mult - incorrect dimension " << parN << " " << K.Ncols() << endl;
    return -1;
  }

  res = mult;

  double c = 0;
  int pos = -1;
  for (size_t i = 1; i <= K.Nrows(); i++) {
    for (size_t j = 1; j <= K.Ncols(); j++) {
      pos = mult.getParam(j);
      if (pos < 0 ) { return -1; }
      c += K(i,j)*mult[j-1].value();
    }
    pos = mult.getParam(i);
    res[i-1].value(c);
    c = 0;
  }
  
  return 1;
}

// Get prn of all AMB params
// --------------
set<string> t_gallpar::ambs()
{
  gtrace("t_gallpar::ambs");   
   
  set<string> prns;
  vector<t_gpar>::const_iterator iter;
  for (iter = _vParam.begin(); iter != _vParam.end(); iter++) {
    if ( iter->parType == t_gpar::AMB_IF || iter->parType == t_gpar::AMB_L1 ) prns.insert(iter->prn);
  }
  return prns;
}

// overloading << operator
// -----------------------------
ostream& operator<<(ostream& os, t_gallpar& x)
{
  for (unsigned int i = 0; i < x.parNumber(); i++){
    if (x[i].parType == t_gpar::AMB_IF || 
        x[i].parType == t_gpar::AMB_L1 || 
        x[i].parType == t_gpar::AMB_L2 ||
        x[i].parType == t_gpar::AMB_L3 ||
        x[i].parType == t_gpar::AMB_L4 ||
        x[i].parType == t_gpar::AMB_L5 ||      
        x[i].parType == t_gpar::SION   ||
        x[i].parType == t_gpar::VION   ||
        x[i].parType == t_gpar::IFB_C3 ||
        x[i].parType == t_gpar::IFB_C4 ||
        x[i].parType == t_gpar::IFB_C5 ||
        x[i].parType == t_gpar::IFCB_F3 ||
        x[i].parType == t_gpar::IFCB_F4 ||        
        x[i].parType == t_gpar::IFCB_F5     ) os << x[i].str_type() << "_" << x[i].prn << " ";
    else os << x[i].str_type() << " ";
      
    if (x[i].parType == t_gpar::GRD_N || x[i].parType == t_gpar::GRD_E)
      os << "value: " << x[i].value()*1000 << " " << "index:" << x[i].index;   
    else os << "value: " << x[i].value() << " " << "index:" << x[i].index;
      os << endl;
   } 
  return os;        
}

// X1 + X2
//---------------------------------
int t_gallpar::sum(t_gallpar& X1, t_gallpar& X2)
{
  if (X1.parNumber() != X2.parNumber() ) return -1;
   
  for (unsigned int i = 0; i < _vParam.size(); i++){
    int id1 = X1.getParam(_vParam[i].site, _vParam[i].parType, _vParam[i].prn);
    int id2 = X2.getParam(_vParam[i].site, _vParam[i].parType, _vParam[i].prn);
    if (id1 >= 0 && id2 >= 0) _vParam[i].value(X1[id1].value() + X2[id2].value());
    else return -1;
  }
  return 1;
}

// get ColumnVector w.r.t par indexes
// ---------------------------
ColumnVector t_gallpar::get_cvect(t_gallpar& par)
{
  int n = par.parNumber();
  ColumnVector vec(n);
   
  for (int i = 1; i <= n; i++){
    int idx = this->getParam(i);
    vec(i) = _vParam[idx].value();
  }
  return vec;
}

} // namespace
