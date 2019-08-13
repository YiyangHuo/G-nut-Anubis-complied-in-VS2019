
#ifndef GXML_H
#define GXML_H

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements base XML class
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <set>
#include <string>
#include <ctype.h>
#include <iostream>

#ifdef BMUTEX
#include <boost/thread/mutex.hpp>
#endif

#include "gio/glog.h"
#include "gutils/gmutex.h"
#include "gutils/gtypeconv.h"
#include "pugixml/src/pugixml.hpp"

using namespace std;
using namespace pugi;

namespace gnut {

class t_gxml : public pugi::xml_node
{
 public:
   t_gxml(string s, bool upper = true );
   virtual ~t_gxml();
  
   virtual int glog_set(t_glog* glog);           // read glog file and set it
   
   int read(const string& file);                 // read from file
   int read_istream(istream& is);                // read from istream
   int write(const string& file);                // write in file
   int write_ostream(ostream& os);               // write in ostream
   
   string root(){ return _root; }   
   
   virtual void check() {};                     // settings check
   virtual void help()  {};                     // settings help

   string         strval(const string& elem, const string& subelem); // public interface
   set<string>    setval(const string& elem, const string& subelem); // public interface
   vector<string> vecval(const string& elem, const string& subelem); // public interface

 protected:
   string         _strval(const string& elem, const string& subelem);
   set<string>    _setval(const string& elem, const string& subelem);
   vector<string> _vecval(const string& elem, const string& subelem);

   xml_node _default_node(xml_node& node, const char* n, const char* val="", bool reset=false);

   void _default_attr(xml_node& node, const char* n, const string& val, bool reset=false);
   void _default_attr(xml_node& node, const char* n, const bool& val,   bool reset=false);
   void _default_attr(xml_node& node, const char* n, const int& val,    bool reset=false);
   void _default_attr(xml_node& node, const char* n, const double& val, bool reset=false);

   virtual void help_header();                   // XML settings header
   virtual void help_footer();                   // XML settings footer

   xml_document     _doc;                        // root document
   xml_parse_result _irc;                        // result status
   
   t_glog*        _logxml;                          // pointer to log file
   string         _name;                         // file name
   string         _root;                         // root directory
   string         _delimiter;                    // delimiter for writting nodes/elements
   bool           _ucase;                        // upper/lower case for keywords
   t_gmutex       _gmutexxml;
#ifdef BMUTEX   
   boost::mutex   _mutex;                        // muttable access
#endif
   
 private:

};

} // namespace

#endif
