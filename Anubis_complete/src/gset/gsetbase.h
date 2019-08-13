
#ifndef GSETBASE_H
#define GSETBASE_H

#define XMLKEY_ROOT "config"

/* ----------------------------------------------------------------------
  (c) 2011 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: implements base setting class
  Version: $ Rev: $

  2012-10-23 /JD: created

-*/

#include <set>
#include <map>
#include <string>
#include <iostream>

#ifdef WIN32
#pragma warning(disable:4250)
#endif

#include "gio/glog.h"
#include "gutils/gtypeconv.h"
#include "pugixml/src/pugixml.hpp"

using namespace std;
using namespace pugi;

namespace gnut {

class t_gsetbase : public pugi::xml_node
{
 public:
   t_gsetbase();
   virtual ~t_gsetbase();
  
   int read(const string& file);                 // read configuration from file
   int read_istream(istream& is);                // read configuration from istream
   int write(const string& file);                // write configuration in file
   int write_ostream(ostream& os);               // write configuration in ostream

   virtual void glog(t_glog* l){ _log = l; }     // set/get glog pointer (VIRTUAL for gsetout!)
   t_glog* glog()const{ return _log; }
   
   virtual string pgm();                         // application name
   virtual string app();                         // application info/version
   virtual void app(const string& pgm, const string& ver,
	            const string& dat, const string& tim,
	            const string& rev, const string& own);

   virtual void usage();                         // application usage
   virtual void arg(int argc, char *argv[], 
               bool add=false, bool thin=false); // command line arguments

   virtual void help() = 0;                      // settings help
   virtual void check() = 0;                     // settings check

   set<string>    types(){ return _set; }        // get all types of actual settings

   string         strval(const string& elem, const string& subelem); // public interface
   set<string>    setval(const string& elem, const string& subelem); // public interface
   vector<string> vecval(const string& elem, const string& subelem); // public interface

   bool thin();                                  // thin execution

   void print_log();
 protected:
   double         _dblatt(const string& elem, const string& attrib);      
   double         _dblval(const string& elem, const string& subelem);
   string         _strval(const string& elem, const string& subelem);
   set<string>    _setval(const string& elem, const string& subelem);
   vector<string> _vecval(const string& elem, const string& subelem);
   
   xml_node _default_node(xml_node& node, const char* n, const char* val="", bool reset=false);

   void _default_attr(xml_node& node, const char* n, const string& v,     bool reset=false);
   void _default_attr(xml_node& node, const char* n, const bool& value,   bool reset=false);
   void _default_attr(xml_node& node, const char* n, const int& value,    bool reset=false);
   void _default_attr(xml_node& node, const char* n, const double& value, bool reset=false);

   void _add_log(string element, string msg);    // report element issues
     
   virtual void _upd_glog(){};                   // update glog (mask,verb)

   virtual void help_header();                   // XML settings header
   virtual void help_footer();                   // XML settings footer

   xml_document     _doc;                        // root document
   xml_parse_result _irc;                        // result status

   t_glog*        _log;                          // pointer to log file
   string         _name;                         // configuration name
   string         _delimiter;                    // delimiter for writting nodes/elements

   string         _pgm;
   string         _ver;
   string         _rev;
   string         _own;
   string         _dat;
   string         _tim;

   set<string>    _set;
   map<string, set<string> >    _chache_log;
   
   mutable t_gmutex        _gmutex;   
   
 private:

};

} // namespace

#endif
