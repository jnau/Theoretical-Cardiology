#ifndef ___RECORDER2CLASS_H_
#define ___RECORDER2CLASS_H_

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cfloat>
#include <climits>
#include <cstdlib>
using namespace std;
#include "tissue.h"


#define COLORMAPDEFAULT 0
#define COLORMAPJET 1
#define COLORMAPHOT 2
#define COLORMAPGRAY 3


class REC2PARAM{
  public:
    string filenamev;
    string filenamevmm;
    string filenamet;
    string filenamevmovie;
    #ifndef VOLTAGE_ONLY
    string filenamec;
    string filenamecmm;
    string filenamecmovie;
    string filenamecnsr;
    string filenamecnsrmm;
    #endif
    string  filenameekg;
    int step;
    int interval;
    bool gzip;

    string ffmpegoptions;
    int pictureformat;//0: PPM ; 1: PNG ; 2: webp

    REC2PARAM(void){interval=1;step=1;filenamev="";filenamevmm="";filenamet="";filenamevmovie="";filenameekg="";gzip=true;ffmpegoptions="";pictureformat=1;
    #ifndef VOLTAGE_ONLY
    filenamec="";filenamecmm="";filenamecnsr="";filenamecnsrmm="";filenamecmovie="";
    #endif
    }
};

class Recorder2{
public:
  void recdata(double t);
  Recorder2(class Tissue *tiss, class REC2PARAM *param);
  virtual ~Recorder2(void);
  void recekg(double t, double posx, double posy, double posz, bool first = true, bool last = true);

  void setvcolormap(int newcolormap){ vcolormap = newcolormap; }
  void setcicolormap(int newcolormap){ cicolormap = newcolormap; }


  int vfrm, cfrm;
  string picdir;
  bool savepics;
  double picvmax;
  double picvmin;
  double piccimax;
  double piccimin;

private:
  static const int filenamemaxchar = 255;

  REC2PARAM *rcprm;
  ofstream *tout; bool IsRectout;
  ofstream *vout; ofstream *vmmout; bool IsRecvout;
  bool IsRecvmovieout;
#ifndef VOLTAGE_ONLY
  ofstream *ciout; ofstream *cimmout; bool IsRecciout;
  ofstream *cnsrout; ofstream *cnsrmmout; bool IsReccnsrout;
  bool IsReccimovieout;
#endif
  ofstream *ekgout; bool IsRecekgout;
  Tissue *tis;
  int X, Y, step;

  bool gzip;

  int vcolormap;
  int cicolormap;
};
#endif /* ___RECORDER2CLASS_H_ */
