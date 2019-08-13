
/* ----------------------------------------------------------------------
  (c) 2011-2016 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
      Research Institute of Geodesy, Topography and Cartography
      Ondrejov 244, 251 65, Czech Republic

  Purpose: QC RINEX format (in future modify, convert, concatenate obs, nav files)
  Version: $ Rev: $

  2011-10-10 /JD: created
  2013-08-16 /JD: released 1.0
  2014-06-26 /JD: released 1.1
  2014-08-12 /JD: released 1.2
  2015-03-13 /JD: released 1.3
  2016-01-27 /JD: released 1.4
  2016-10-05 /JD: released 2.0
  2017-12-10 /JD: released 2.1
  2018-08-01 /JD: released 2.2

-*/

#include <thread>
#include "gset/gcfg_anubis.h"

using namespace std;
using namespace gnut;

void catch_signal(int){ cout << "Program interrupted by Ctrl-C [SIGINT,2]\n"; exit(1); }

/* -----
 * MAIN
-*/
int main(int argc, char *argv[])
{   
  clock_t clktot = clock();
  clock_t clkxtr = clock();
  signal(SIGINT, catch_signal);
  vector<t_glog> vlogs;
  bool   thrd = false;
  bool   thin = false;
  bool   redu = false;
  bool   dat1 = false;
  bool   exml = false;
  bool   gkpi = false;
  double cfdi = -1;
  int    irc  = 0;
  string main = "main";
  t_glog glog; glog.mask("anubis.log"); glog.cache_size(30);

  // READ COMMAND-LINE
  for( int i=1; i<argc; i++ ){
    if(      !strcmp(argv[i],"-1") && i  <argc ){ dat1 = true; }
    else if( !strcmp(argv[i],"-m") && i  <argc ){ thrd = true; }
    else if( !strcmp(argv[i],"-r") && i  <argc ){ redu = true; }
    else if( !strcmp(argv[i],"-e") && i  <argc ){ exml = true; }
    else if( !strcmp(argv[i],"-g") && i  <argc ){ gkpi = true; }
    else if( !strcmp(argv[i],"-c") && i+1<argc ){ cfdi = atof(argv[++i]); }
    else if( !strcmp(argv[i],"-t") && i  <argc ){ thin = true; main = "thin-exec"; glog.verb(-1); }
  }

  t_gcfg_anubis gset;  gset.glog(&glog);
  gset.app("G-Nut/Anubis", "2.2.4", "$Rev: 2554 $","(gnss@pecny.cz)", __DATE__, __TIME__);
  gset.arg( argc, argv, true, thin);

  bool chk_health = dynamic_cast<t_gsetinp*>(&gset)->chkHealth();
  bool chk_navig  = dynamic_cast<t_gsetinp*>(&gset)->chkNavig();

  // SETTING
  // --------
  t_gdata*    gdata  = 0;
  t_gallobj*  gobj   = new t_gallobj;     gobj->glog( &glog );
  t_gallprec* gnav   = new t_gallprec;    gnav->glog( &glog ); gnav->multi(redu);
              gnav->use_clksp3(true);   // USE ALSO SP3 CLK!
              gnav->use_clknav(true);   // USE ALSO NAV CLK!
              gnav->use_posnav(true);   // USE ALSO NAV ORB!	 
              gnav->chk_health(chk_health);
              gnav->chk_navig(chk_navig);

  vector<t_galloqc*> gobs;
  vector<t_gcoder*> gcoder;
  vector<t_gfile*> gfile;
  vector<t_gxtrqc*> gxtr;

  vector<thread> gthread;

/*
  set<string>::const_iterator itSAT;
  map<GSYS,set<string>>::const_iterator itGNS;
  map<GSYS,set<string>> gsat;
 // TESTING
  gsat[GPS]  = dynamic_cast<t_gsetgnss*>(&gset)->sat(GPS);
  gsat[GLO]  = dynamic_cast<t_gsetgnss*>(&gset)->sat(GLO);
  gsat[GAL]  = dynamic_cast<t_gsetgnss*>(&gset)->sat(GAL);
  gsat[BDS]  = dynamic_cast<t_gsetgnss*>(&gset)->sat(BDS);
  gsat[SBAS] = dynamic_cast<t_gsetgnss*>(&gset)->sat(SBAS);
  gsat[QZSS] = dynamic_cast<t_gsetgnss*>(&gset)->sat(QZSS);

  for( itGNS = gsat.begin(); itGNS != gsat.end(); ++itGNS ){
    cout << "satellite[" << t_gsys::gsys2str(itGNS->first) << "]:";
    for( itSAT = itGNS->second.begin(); itSAT != itGNS->second.end(); ++itSAT ){
      cout << " " << *itSAT;
    }
    cout << endl;

//    cout << "type[" << t_gsys::gsys2str(itGNS->first) << "]:";
//    for( itSAT = itGNS->second.begin(); itSAT != itGNS->second.end(); ++itSAT ){
//      cout << " " << *itSAT;
//    }
//    cout << endl;
  }
*/
   
  set<string>::const_iterator itOBJ;
  set<string> sites = dynamic_cast<t_gsetgen*>(&gset)->rec();
   
  // CREATE USER REQUEST OBJECTs (TO BE COMPARED WITH RINEX HEADER)
  set<string> obj   = dynamic_cast<t_gsetrec*>(&gset)->objects();
  for( itOBJ = obj.begin(); itOBJ != obj.end(); ++itOBJ ){
    string name = *itOBJ;
    shared_ptr<t_grec> rec = dynamic_cast<t_gsetrec*>(&gset)->grec(name, &glog);
    gobj->add( rec );
  }

  // PREPARE INPUTS
  // ---------------
  string path;
  int idx = 0;
  multimap<IFMT,string> inp = gset.inputs_all();
  multimap<IFMT,string>::const_iterator itINP = inp.begin();

  if( thin && inp.size() > 1 ){ glog.comment(0,main,"No single file for ThinExec - stopped!") ; return -1; }

  for( size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP ){
     
     IFMT   ifmt(itINP->first);
     path = itINP->second;
     string id("ID" + int2str(i));
     gfile.push_back(new t_gfile);
     
     if(       ifmt == RINEXO_INP  ){ if(!dat1 || gobs.size()==0) gobs.push_back(new t_galloqc); // single/separate containers
	                              gobs[gobs.size()-1]->maxepoch(0);
	                              gobs[gobs.size()-1]->glog(&glog);
	                              gobs[gobs.size()-1]->gset(&gset);
//	                              gobs[gobs.size()-1]->gfile(gfile[idx]);
	                              gdata = gobs[gobs.size()-1];
	                              gcoder.push_back(new t_rinexo(&gset,"3.02",4096));
     }else if( ifmt == RINEXN_INP  ){ gdata = gnav;
	                              gcoder.push_back(new t_rinexn(&gset,"3.00",4096));
     }else if( ifmt == SP3_INP     ){ gdata = gnav;
	                              gcoder.push_back(new t_sp3(&gset,"",8172));
     }else{
       if( !thin ) glog.comment(0,main,"Error: decoder not supported : "+ int2str(ifmt) );
       gdata = 0;
       continue;
     }

     gcoder[idx]->clear();
     gcoder[idx]->glog(&glog);
     gcoder[idx]->add_data("OBJ", gobj);
     gcoder[idx]->add_data(id, gdata);
     gcoder[idx]->path(path);
     gcoder[idx]->close_with_warning(true);

     gfile[idx]->glog(&glog);
     gfile[idx]->coder(gcoder[i]);
     gfile[idx]->path(path);
     idx++;
  }

  // READ INPUT
  // -----------
  for( size_t i = 0; i < gfile.size(); ++i )
  {
    if( ! thrd ){
      clock_t tmpclk = clock();
      glog.comment(1,main,"READ: " + gfile[i]->path() + " started");
      gfile[i]->run_read();
      glog.comment(0,main,"READ: " + gfile[i]->path() + " " + clk_string(tmpclk));
      // THIN EXEC
      vector<t_gnote> mesg = gcoder[i]->mesg();
      irc += gcoder[i]->irc();
      irc += gfile[i]->irc();
      for( auto itMSG = mesg.begin(); itMSG != mesg.end(); ++itMSG ){
        if( thin ) cerr << "*** " << itMSG->str() << endl;
//      else       glog.comment(0,main,itMSG->str());
      }
    }else{
      glog.comment(1,main,"Multi-thread reading: " + gfile[i]->path() + " started");
      gthread.push_back( thread( &t_gfile::run_read, gfile[i]) );
    }
  }

  // CLEAN XTR
  //-----------
  for( size_t i = 0; i < gfile.size();  ++i ){ delete gfile[i];  }; gfile.clear();
  for( size_t i = 0; i < gcoder.size(); ++i ){ delete gcoder[i]; }; gcoder.clear();

  // check navigation messages
  if( cfdi >= 0.0 ){
    glog.comment(1,main,"Consolidating navigation data");
    gnav->consolidate(cfdi); gnav->clean_invalid();
  }

  // CLEAN READING
  // --------------
  clkxtr = clock();
  if( !thin && thrd ){
    for( size_t i = 0; i < gthread.size(); ++i ) {gthread[i].join();} 
    gthread.clear();
    glog.comment(1,main,"Multi-thread input finished: " + clk_string(clkxtr));
  }

  // WRITE RINEXN
  // ============
  set<string> outf = gset.oformats();
  string ofmt = "rinexn";
  if( !thin && outf.find( ofmt ) != outf.end() ){
    string path = gset.outputs( ofmt );
    t_rinexn* gencode = new t_rinexn(&gset,"3.02",4096);
    gencode->glog(&glog);
    gencode->add_data("OBJ", gnav);

    t_gfile* gout = new t_gfile;
    gout->glog(&glog);
    gout->path(path);
    gout->coder(gencode);
          
    clock_t tmpclk = clock();
    glog.comment(1,main,"SAVE: " + gout->path() + " started");
    gout->run_write();
    glog.comment(0,main,"SAVE: " + gout->path() + " " + clk_string(tmpclk));
     
    delete gout;
    delete gencode;
  }
  // RINEXN2
  // =======
  ofmt = "rinexn2";
  if( !thin && outf.find( ofmt ) != outf.end() ){
    string path = gset.outputs( ofmt );
    t_rinexn* gencode = new t_rinexn(&gset,"2.11",4096);
    gencode->glog(&glog);
    gencode->add_data("OBJ", gnav);

    t_gfile* gout = new t_gfile;
    gout->glog(&glog);
    gout->path(path);
    gout->coder(gencode);
          
    clock_t tmpclk = clock();
    glog.comment(1,main,"SAVE: " + gout->path() + " started");
    gout->run_write();
    glog.comment(0,main,"SAVE: " + gout->path() + " " + clk_string(tmpclk));
     
    delete gout;
    delete gencode;
  }
  // XTR SUMMARY
  // ------------
  string D(__DATE__); transform(D.begin(), D.end(), D.begin(), ::toupper);
  t_gtime tt(t_gtime::UTC); tt.from_str("%M %d %Y %H:%M:%S", string(D)+string(__TIME__));

  for( size_t i = 0; !thin && i < gobs.size(); ++i )
  {
     set<string> all_sites = gobs[i]->stations();
//  if( all_sites.size() == 0 ) continue;
     
    set<string>::const_iterator itSIT;
    for( itSIT = all_sites.begin(); itSIT != all_sites.end(); ++itSIT )
    {
      string site = *itSIT;

      // check correctly 4-ch vs 9-ch (file/settings)
      bool site_ok = false;
      if( sites.size() == 0 ) site_ok = true; // no specification
      for(set<string>::iterator itSIT = sites.begin(); itSIT != sites.end(); ++itSIT ){
        if(                                                    itSIT->compare(site)         == 0 // equal (4-char or 9-char)
	    || ( itSIT->length() == 4 && site.length() == 9 && itSIT->compare(0,4,site,0,4) == 0)
	    || ( itSIT->length() == 9 && site.length() == 4 && itSIT->compare(0,4,site,0,4) == 0) )
	{ site_ok = true; break; }
      }
      if( ! site_ok ) continue;
//      if( sites.size() != 0 && sites.find(site)             == sites.end() &&
// 	                       sites.find(site.substr(0,4)) == sites.end() ) continue; 

      size_t idx = gxtr.size();
      if( gkpi ){ gxtr.push_back(new t_gxtrgrc(&gset, gset.app(), tt)); }
      else{       gxtr.push_back(new t_gxtrqc( &gset, gset.app(), tt)); }
      gxtr[idx]->glog(&glog);
      gxtr[idx]->extended_xml(exml);
      gxtr[idx]->setOBJ(gobj);
      gxtr[idx]->setDAT(gobs[i], gnav);

      gxtr[idx]->setSIT(site);
      
      if( ! thrd ){
        clock_t tmpclk = clock();
        glog.comment(1,main,"Single-thread summary: " + site + " started");
        gxtr[idx]->summary( site );
        glog.comment(0,main,"Single-thread summary: " + site + " " + clk_string(tmpclk));
      }else{
        glog.comment(1,main,"Multi-thread summary: " + site + " started");
        gthread.push_back( thread(&t_gxtrqc::summary, gxtr[idx], site) );
      }
    }
  }

  // CLEAN XTR
  //-----------
  if( thrd ){ 
    for( size_t i=0; i<gthread.size(); ++i ) {gthread[i].join();}
    gthread.clear();
    glog.comment(1,main,"Multi-thread summary finished: " + clk_string(clkxtr));
  }
  for( size_t i = 0; i < gxtr.size();   ++i ){ delete gxtr[i];   }; gxtr.clear();
  for( size_t i = 0; i < gobs.size();   ++i ){ delete gobs[i];   }; gobs.clear();
  delete gobj;
  delete gnav;

  glog.comment(0,main, "total time: " + clk_string(clktot));
  return irc;
}
