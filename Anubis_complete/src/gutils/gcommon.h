
#ifndef GCOMMON_H
#define GCOMMON_H 
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: 
  Version: $ Rev: $

  2011-02-14 /JD: created

-*/

#include <string>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <stdio.h>  // defines FILENAME_MAX
#include <sys/stat.h>

#ifdef _WIN32
#include <direct.h>
#include <Shlwapi.h>
#define GET_CURRENT_PATH _getcwd
#pragma warning(disable:4503)                   // suppress Visual Studio WARNINGS 4503 about truncated decorated names

#else
#include <unistd.h>
#define GET_CURRENT_PATH  getcwd
#endif

#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif

using namespace std;

namespace gnut {

static const char lf = 0x0a; // = "\n" .. line feed
static const char cr = 0x0d; // = "\r" .. carriage return

#if defined __linux__ || defined __APPLE__
static const string crlf = string("") + lf;   // "\n"
#endif

#if defined _WIN32 || defined _WIN64
//static const string crlf = string("") + cr + lf; // "\r\n"
static const string crlf = string("") + lf; // "\r\n"
#endif 

// ----------
inline string cut_crlf(string s)
{
  if( !s.empty() && s[s.length()-1] == lf ) s.erase(s.length()-1);
  return s;
}
   
   
// ----------
inline void gtrace(const string& str)
{
#ifdef TRACE
  cout << "#" << str << endl;
#endif
}


// ----------
inline string get_current_path()
{
  char cPath[FILENAME_MAX];
  if( ! GET_CURRENT_PATH(cPath, sizeof(cPath)) ) return "./";
   
  cPath[sizeof(cPath) - 1] = '\0'; // not really required
  return string(cPath);
}


// ----------
inline bool chk_directory(const string& dir)
{
struct stat st;

  if( stat(dir.c_str(), &st) == 0 && S_ISDIR(st.st_mode) ) return true;
// if( stat(dir.c_str(), &st) == 0 && (((st.st_mode) & _S_IFMT) == _S_IFDIR) ) return true;
//  if (stat(dir.substr(0,255).c_str(), &st) == 0 && st.st_mode == _S_IFDIR ) return true;

//  if( stat(dir.c_str(),&st) == 0 )
//    if( st.st_mode & S_IFDIR != 0 ) return true;  //  printf(" /tmp is present\n");

  return false;
}  

   
// ----------
inline void gclk(const string& str, time_t& clk)  // usage:  clock_t clk; gclk("TEXT",clk);
{
#ifdef CLOCK
  cout << "#"
       << fixed << setprecision(4)
       << setw(7) << difftime(time(0) - clk) << " sec: "
       << str << " " << endl;
  clk = time(0); // clock();
#endif
}

   
// this is designed for main program
// ----------
inline string clk_string( clock_t clk )
{
  clk = clock() - clk;
  ostringstream os;
  os << fixed << setprecision(3) << setw(7) << ((float)clk)/CLOCKS_PER_SEC << " sec ";
  return os.str();
}


/*
inline void gclk(const string& str, clock_t& clk)  // usage:  clock_t clk; gclk("TEXT",clk);
{
#ifdef CLOCK
  cout << "#"
       << fixed << setprecision(4)
       << setw(7) << (double)( clock()-clk )/CLOCKS_PER_SEC << " sec: "
       << str << " " << endl;
  clk = clock();
#endif
}
*/

#ifdef DEBUG
 #define ID_FUNC printf(" .. %s\n", __func__)
#else
 #define ID_FUNC 
#endif

#define SQR(x)   ((x)*(x))
#define SQRT(x)  ((x)<=0.0?0.0:sqrt(x))

// tail of a string
// ---------------
inline string tail(const string& str, const size_t length)
{
   if(length >= str.size()) return str;
   else return str.substr(str.size() - length);
}
   
   
} // namespace

#endif // # GCOMMON_H
