#ifndef ___TISSUE_H
#define ___TISSUE_H
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
using namespace std;

//1D and 2D tissue
#ifdef ___USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#include "cell.h"

#define STIM_2D_CORNER1 0
#define STIM_2D_PARALLEL1 1
#define STIM_2D_PARALLEL2 6
#define STIM_2D_CORNER2 2
#define STIM_2D_CORNER3 3
#define STIM_2D_CORNER4 4
#define STIM_2D_CENTER 5
#define STIM_1D_LEFT 0
#define STIM_1D_RIGHT 1
#define STIM_2D_POS 99
#define STIM_1D_POS 99

class Tissue{
 public:
  Tissue(int X,int Y=1, bool CaDiff=false);
  Tissue(const Tissue &oldtis);
  virtual ~Tissue();
  void pace(double stim=0.0,int pos=0,int posx=0 ,int posy=0);//0:cor1 1:para 2:cor2 for 2d, 0:left 1:right for 1d
  void paceall(double stim=0.0);
  CCell *cell;
  void copydata(CCell *cell,int x1=-1,int x2=-1, int y1=-1, int y2=-1);
  double setdt(double dtt);
  int setDifLoop(int DL);
  int getX(void){return X;}
  int getY(void){return Y;}
  int getGX(void){return GX;}
  int getGY(void){return GY;}
  double getdt(void){return dt;}
  double getdx(void){return dx;}
  double setdx(double newdx){dx=newdx; return dx;}
//    friend ostream &operator<<(ostream &stream,Tissue *tiss);
  Tissue& operator=(const Tissue& tis);
  void setR2(int newR2){R2=newR2;Stripe=R2;}
  void setDfuCai(double D){DfuCa=D;}
  void setdif(double D);
  double getdif(void){return Dfu;}



  // for phase field
  void solvephi(void);
  void setphi(int xx,int yy,double phii=1);
  void writephi(string filename="phi.txt",int step=1);
  void setDxx(double theta,double dxx, double dyy, int posx=-1, int posy=-1);

  double *Dxx;
  double *Dxy;
  double *Dyy;
  double *phi;
  void setxi(double newxi){xi=newxi;}
  double getxi(void){return xi;}
  void setcutoff(double newcutoff){cutoff=newcutoff;}
  double getcutoff(void){return cutoff;}
  // end for phase field

 private:
  double Dfu;//V Diffusion Cofficient
  double DfuCa;//Ca Diffusion Cofficient
  double dx;
  double dt;
  int R2;
  int Stripe;
  int X,Y;//tissue size in this process
  int GX,GY;//tissue size (total)
  bool IsCaDiff;


  double *tmpv;
  double ***diffwork;
  int DifLoop;
  #ifdef ___USE_MPI
  int num_procs;
  double *mpiwork_send1;
  double *mpiwork_send2;
  double *mpiwork_recv1;
  double *mpiwork_recv2;
  double *Dxxphi1;
  double *Dxxphi2;
  double *Dxyphi1;
  double *Dxyphi2;
  #endif
  int my_proc;


  // for phase field
  double cutoff;
  double xi;
  void diffusion(void);
  void diffusionf(void);
  // end for phase field


};



#endif /* ___TISSUE_H */

