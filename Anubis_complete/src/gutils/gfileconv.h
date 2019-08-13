
#ifndef GFILECONV_H
#define GFILECONV_H
 
/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: file conversion utilities
  Version: $ Rev: $

  2016-06-26 /JD: created

-*/

#include <string>

#if defined _WIN32 || defined _WIN32
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

#define GFILE_PREFIX "file://"

using namespace std;

namespace gnut {


string base_name(const string& path);               // extract file base name
string  dir_name(const string& path);               // extract file dir  name
string file_name(const string& path);
bool  dir_exists(const string& path);               // check existance of path
int    make_path(const string& path);               // create path recursively
int    make_dir(const string& path);                // create single directory
} // namespace
   
#endif