
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

#include <iomanip>
#include <sstream>
#include <algorithm>

#include "gio/gxml.h"
#include "gutils/gfileconv.h"

using namespace std;
using namespace pugi;

namespace gnut {

// Constructor
// ----------
t_gxml::t_gxml(string s, bool upper)
{
  _root      = s;
  _logxml    = 0;
  _name      = "";
  _delimiter = "  "; // only for nodes/elements
  _ucase     = true;
}


// Destructor
// ----------
t_gxml::~t_gxml()
{}


// Read/set glog
int t_gxml::glog_set(t_glog* glog )
{
  _logxml = glog;

  return 0;
}


// Read from file
// ----------
int t_gxml::read(const string& file)
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutexxml.lock();

  _name = file;

  if( ! (_irc = _doc.load_file( _name.c_str() )) ){ 
    if( _logxml ) _logxml->comment(0, "gxml" , "XML-config not read file " + _name + " " + string(_irc.description()) );
    else                         cerr << "XML-config not read file " + _name + " " + string(_irc.description()) << endl;
    { _gmutexxml.unlock(); return -1; }
  }

  if( _logxml ) _logxml->comment(0, "gxml" , "XML-config read from file " + _name );

#ifdef BMUTEX
  lock.unlock();
#endif
   
  _gmutexxml.unlock();
  this->check();

  return 0;
}


// Read from istream
// ----------
int t_gxml::read_istream(istream& is)
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutexxml.lock();

  if( ! (_irc = _doc.load( is ) ) ){
    if( _logxml ) _logxml->comment(0, "gxml" , "XML-config not read istream: " + string(_irc.description()) );
    else                             cerr << "XML-config not read istream: " + string(_irc.description()) << endl;
    { _gmutexxml.unlock(); return -1; }
   
  }
   
  if( _logxml ) _logxml->comment(0, "gxml" , "XML-config read from istream");

#ifdef BMUTEX
  lock.unlock();
#endif
  _gmutexxml.unlock();
  this->check();
  return 0;
}


// Write to file
// ----------
int t_gxml::write(const string& file)
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutexxml.lock();
   
  ofstream of;
  string name(file);   // if( name.empty() ) name = _name;  // ???
  substitute(name,GFILE_PREFIX,"");     
   
  try{
      of.open(name.c_str());
  }catch( fstream::failure e ){
    cerr << "Warning: Exception opening file " << name << ": " << e.what() << endl;
    _gmutexxml.unlock(); return -1;
  }

  _doc.save_file(name.c_str(), _delimiter.c_str(), pugi::format_indent);
  of.close();

  if( _logxml ) _logxml->comment(0,"gxml","XML-file saved: " + name );
//else               cout << "gxml:  XML-file saved: " + name << endl;

  _gmutexxml.unlock(); return 0;
}


// Write to ostream
// ----------
int t_gxml::write_ostream(ostream& os)
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutexxml.lock();

  _doc.save(os, _delimiter.c_str(), pugi::format_indent);
  
  _gmutexxml.unlock(); return 0;
}
 
   
// Return first string of value based on key
// ----------
string t_gxml::strval( const string& elem, const string& subelem )
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutexxml.lock();
  string tmp = _strval( elem, subelem );
  _gmutexxml.unlock();
  return tmp;
}


// Return set of value based on key
// ----------
set<string> t_gxml::setval( const string& elem, const string& subelem )
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutexxml.lock();

  set<string> tmp = _setval( elem, subelem );
  _gmutexxml.unlock();
  return tmp;
}


// Return vector of value based on key
// ----------
vector<string> t_gxml::vecval( const string& elem, const string& subelem )
{
#ifdef BMUTEX
  boost::mutex::scoped_lock lock(_mutex);
#endif
  _gmutexxml.lock();

  vector<string> tmp = _vecval( elem, subelem );
  _gmutexxml.unlock(); 
  return tmp;
}


// Return first string of value based on key
// ----------
string t_gxml::_strval( const string& elem, const string& subelem )
{
  string word;

  istringstream is( _doc.child(_root.c_str()).child( elem.c_str() ).child_value( subelem.c_str() ));
   
  if( is >> word ) return word;
  return "";
}


// Return set of value based on key
// ----------
set<string> t_gxml::_setval( const string& elem, const string& subelem )
{
  set<string> vals;
  string word;

  istringstream is( _doc.child(_root.c_str()).child( elem.c_str() ).child_value( subelem.c_str() ));
  while( is >> word ) vals.insert( word );
               
  return vals;
}

// Return vector of value based on key
// ----------
vector<string> t_gxml::_vecval( const string& elem, const string& subelem )
{
  vector<string> vec;
  string str;

  istringstream is( _doc.child(_root.c_str()).child( elem.c_str() ).child_value( subelem.c_str() ));
  while( is >> str ) vec.push_back( str );
               
  return vec;
}


// check/set node
// ----------
xml_node t_gxml::_default_node(xml_node& node, const char* n, const char* val, bool reset)
{
#ifdef DEBUG
  cout << "creating node [" << n << "]  for parent [" << node.name() << "]  set with value [" << val << "]\n";
#endif

  string s(n);

//  if( _ucase )  transform(s.begin(), s.end(), s.begin(), ::toupper);
//  else          transform(s.begin(), s.end(), s.begin(), ::tolower);
   
  // remove if to be reset
//  if( strcmp(val,"") && reset ) node.remove_child(s.c_str());
  if( reset ) node.remove_child(s.c_str());

#ifdef DEBUG
  if( ! node ) cerr << " NODE[" << node.name() << "] " << node << " does'not exists !\n";
  else         cerr << " NODE[" << node.name() << "] " << node << " does exists .... ok \n";
#endif
   
  xml_node elem = node.child(s.c_str());

#ifdef DEBUG
  if( ! elem ) cerr << " ELEM[" << s << "] " << elem << " does'not exists !\n";
  else         cerr << " ELEM[" << s << "] " << elem << " does exists .... ok \n";
#endif

  if( ! elem ) elem = node.append_child(s.c_str());


  if( ! elem ){
    if( _logxml ) _logxml->comment(0,"gxml","warning - cannot create element " + s);
    else               cerr << "gxml:  warning - cannot create element " + s + "\n";
  }
   
//  if( elem && strcmp(val,"") && reset ){
  if( elem && (strcmp(val,"") || reset) ){
    elem.append_child(pugi::node_pcdata).set_value(val);
  }

  return elem;
}


// check attributes & set default
// ----------
void t_gxml::_default_attr(xml_node& node, const char* n, const string& val, bool reset)
{
  if( node.attribute(n).empty() ) node.append_attribute(n);
  if( strlen(node.attribute(n).value()) == 0 || reset ) node.attribute(n).set_value(val.c_str());
   
  string s = node.attribute(n).value();
//  if( _ucase ) transform(s.begin(), s.end(), s.begin(), ::toupper);
//  else         transform(s.begin(), s.end(), s.begin(), ::tolower);
  node.attribute(n).set_value(s.c_str());

//  cerr << " strcmp(val,) " << strcmp(val.c_str(),"") 
//       << " reset:" << reset << " val:" << val << "." << endl;

}


// check attributes & set default
// ----------
void t_gxml::_default_attr(xml_node& node, const char* n, const bool& val, bool reset)
{
  if( node.attribute(n).empty() ) node.append_attribute(n);
  if( strlen(node.attribute(n).value()) == 0 || reset ) node.attribute(n).set_value(val);
}


// check attributes & set default
// ----------
void t_gxml::_default_attr(xml_node& node, const char* n, const int& val, bool reset)
{ 
  if( node.attribute(n).empty() ) node.append_attribute(n);
  if( strlen(node.attribute(n).value()) == 0 || reset ) node.attribute(n).set_value(val);
}


// check attributes & set default
// ----------
void t_gxml::_default_attr(xml_node& node, const char* n, const double& val, bool reset)
{
  if( node.attribute(n).empty() ) node.append_attribute(n);
  if( strlen(node.attribute(n).value()) == 0 || reset ) node.attribute(n).set_value(val);
}


// XML help header (optional text comment)
// ----------
void t_gxml::help_header()
{
  cerr << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?> \n"
       << "<!DOCTYPE " << _root << ">\n\n"
       << "<" << _root << ">\n";
}


// XML help footer (optional text comment)
// ----------
void t_gxml::help_footer()
{
  cerr << "</" << _root << ">\n";
}

} // namespace
