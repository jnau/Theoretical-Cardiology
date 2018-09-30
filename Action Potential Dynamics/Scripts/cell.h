#ifndef ___CELL_H_
#define ___CELL_H_

#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>

class CCell{

public:
  double jparam;//tauj*jparam
  double pace(double stim=0);
  double pacevclamp(double clampv);
  double setdt(double dtt){dt=dtt;return dt;}
  double getdt(void){return dt;}
  int getdim(void){return N;}
  double getvc(void){return vc;}
  double getstim(void){return stim;}
  double getstimduration(void){return stimduration;}
  void apclamp(double t, double bcl, double Vmin=-80.0, double Vmax=30.0, double APD=0);
  CCell(void);
  virtual ~CCell();
  CCell& operator=(const CCell& cell);
  void prepare(double BCL=300, int Iter=0);
  double *y;
  double &m,&h,&j,&xr,&xs1,&xs2,&xtos,&ytos,&v,&ci,&cs,&cnsr,&cjsr,&cp;
  double &xir,&c1,&c2,&xi1ca,&xi1ba,&xi2ca,&xi2ba,&nai,&xtof,&ytof,&tropi,&trops;

  double gca;//ica conductance
  double gtos;// ito slow conductance 
  double gtof;// ito fast conductance 
  double gnaca;// exchanger strength 
  double gks;
  double gkr;
  double vup;
  double gna;// sodium conductance (mS/micro F) 
  double gkix;// Ik1 conductance
  double gnak;

  double nao;//mM external Na
  double cao;//mM external Ca
  double ek;

  double taud;//  diffusional delay (ms)
  double taur;// spark lifetime (ms)
  double taua;// NSR-JSR diffusional delay (ms)
  double av;
  double cstar;
  double gleak;
  double knsr;
  double cup;// uptake threshold


  double inaca,ica,iks,ikr,itof,itos,ik1,ina,inak;

  double getvmax(void){return 10;}
  double getvmin(void){return -80;}
  double setki(double newki){ ki = newki; ek = (1.0 / frt)*log(ko / ki); return ki; }
  double setko(double newko){ ko = newko; ek = (1.0 / frt)*log(ko / ki); return ko; }
  double getki(void){return ki; }
  double getko(void){return ko; }

private:
  double pacex(double stim=0);
  static const int N=26;
  static const double vc;
  static const double stim;
  static const double stimduration;
  static const double temp;// temperature (K)
  static const double xxr;//
  static const double xf;// Faraday's constant
  static const double frt;

  double comp_ina (void);
  double comp_ikr(void);
  double comp_iks(void);
  double comp_ik1(void);
  double comp_ito(void);
  double comp_inak(void);
  double comp_jnaca(double csm);
  double comp_icalpo(void);
  double comp_jup(void);
  double comp_jleak(void);
  double comp_inst_buffer(double c);

  double comp_rxa(double csm);
  double comp_Q(void);
  double comp_dir(double po, double Qr, double rxa, double dcj);
  double comp_dcp(double po, double Qr, double rxa);
  double vold;
  double dt,dtx;

  double ki;//mM internal K
  double ko;//mM external K


};

#endif /* ___CELL_H_ */
