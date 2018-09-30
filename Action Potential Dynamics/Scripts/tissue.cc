//1D and 2D tissue

#include "tissue.h"
#include "tissue_diff.cc"
Tissue::Tissue(int sizeX,int sizeY, bool CaDiff) {
  R2=100;
  Stripe=5;
  Dfu=0.001;// or 0.0005
  dx = 0.015;
  GX=X=sizeX;
  GY=Y=sizeY;
  DifLoop=4;
  DfuCa=3E-9;//300[um^2/s]=300*10^-8*10^-3= is upper limit
  my_proc=0;

  cutoff=-1;//1E-5;
  xi=0.025;
  
  if (Y==1) {
    IsCaDiff=CaDiff;
    R2=5;
    cell=new CCell  [X];
    tmpv=new double[X];
    dt=cell[0].getdt();
  }
  else {
  #ifdef ___USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (GX%num_procs)exit(EXIT_FAILURE);
    X=GX/num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
    mpiwork_send1=new double[Y];
    mpiwork_send2=new double[Y];
    mpiwork_recv1=new double[Y];
    mpiwork_recv2=new double[Y];
    Dxxphi1=new double[Y];
    Dxxphi2=new double[Y];
    Dxyphi1=new double[Y];
    Dxyphi2=new double[Y];
  #endif
    Dxx=new double[X*Y];
    Dxy=new double[X*Y];
    Dyy=new double[X*Y];
    phi=new double[X*Y];
    for(int i=0;i<X;i++) {
      for(int j=0;j<Y;j++) {
        Dxx[i*Y+j]=Dfu;
        Dxy[i*Y+j]=0;
        Dyy[i*Y+j]=Dfu;
        phi[i*Y+j]=1;
      }
    }

    diffwork=new double**[X];
    for (int i=0;i<X;i++)
      diffwork[i] = new double* [Y];
    for(int i=0;i<X;i++) {
      for(int j=0;j<Y;j++) {
        diffwork[i][j] = new double [DifLoop+1];
      }
    }
    cell=new CCell  [X*Y];
    tmpv=new double[X*Y];
    dt=cell[0].getdt();
  }
}
Tissue::~Tissue() {
  if (Y==1) {
    delete[] tmpv;
    delete[] cell;
  }
  else {
    for(int i=0;i<X;i++) {
      for(int j=0;j<Y;j++) {
        delete[] diffwork[i][j];
      }
    }
    for (int i=0;i<X;i++) {
      delete[] diffwork[i];
    }
    delete[] diffwork;
    delete[] cell;
    delete[] tmpv;
#ifdef ___USE_MPI
    delete[] mpiwork_send1;
    delete[] mpiwork_send2;
    delete[] mpiwork_recv1;
    delete[] mpiwork_recv2;
    delete[] Dxxphi1;
    delete[] Dxxphi2;
    delete[] Dxyphi1;
    delete[] Dxyphi2;
#endif
    delete[] Dxx;
    delete[] Dxy;
    delete[] Dyy;
    delete[] phi;
  }
}
Tissue::Tissue(const Tissue& tis) {
//constructor
  R2=tis.R2;
  Dfu=tis.Dfu;
  dx = tis.dx;
  dt=tis.dt;
  GX=tis.GX;
  X=tis.X;
  GY=tis.GY;
  Y=tis.Y;
  DifLoop=tis.DifLoop;
  DfuCa=tis.DfuCa;
  my_proc=tis.my_proc;

  IsCaDiff=tis.IsCaDiff;
  #ifdef ___USE_MPI
  X=GX/num_procs;
  mpiwork_send1=new double[Y];
  mpiwork_send2=new double[Y];
  mpiwork_recv1=new double[Y];
  mpiwork_recv2=new double[Y];
  #endif
  diffwork=new double**[X];
  for (int i=0;i<X;i++)
    diffwork[i] = new double* [Y];
  for(int i=0;i<X;i++) {
    for(int j=0;j<Y;j++) {
      diffwork[i][j] = new double [DifLoop+1];
    }
  }
  cell=new CCell  [X*Y];
  for (int i=0;i<X;i++) {
    for (int j=0;j<Y;j++) {
      cell[i*Y+j]=tis.cell[i*Y+j];
    }
  }
  tmpv=new double[X*Y];
}
Tissue& Tissue::operator=(const Tissue& tis) {
  if (&tis!=this) {
    if (Y==1) {
      delete[] tmpv;
      delete[] cell;
    }
    else {
      for(int i=0;i<X;i++) {
        for(int j=0;j<Y;j++) {
          delete[] diffwork[i][j];
        }
      }
      for (int i=0;i<X;i++) {
        delete[] diffwork[i];
      }
      delete[] diffwork;
      delete[] cell;
      delete[] tmpv;
  #ifdef ___USE_MPI
      delete[] mpiwork_send1;
      delete[] mpiwork_send2;
      delete[] mpiwork_recv1;
      delete[] mpiwork_recv2;
  #endif
    }

  //constructor
    R2=tis.R2;
    Dfu=tis.Dfu;
    dx = tis.dx;
    dt=tis.dt;
    GX=tis.GX;
    X=tis.X;
    GY=tis.GY;
    Y=tis.Y;
    DifLoop=tis.DifLoop;
    DfuCa=tis.DfuCa;
    my_proc=tis.my_proc;

    if (Y==1) {
      IsCaDiff=tis.IsCaDiff;
      R2=tis.R2;
      cell=new CCell  [X];
      for (int i=0;i<X;i++) {
        cell[i]=tis.cell[i];
      }
      tmpv=new double[X];
    }
    else {
    #ifdef ___USE_MPI
      X=GX/num_procs;
      mpiwork_send1=new double[Y];
      mpiwork_send2=new double[Y];
      mpiwork_recv1=new double[Y];
      mpiwork_recv2=new double[Y];
    #endif
      diffwork=new double**[X];
      for (int i=0;i<X;i++)
        diffwork[i] = new double* [Y];
      for(int i=0;i<X;i++) {
        for(int j=0;j<Y;j++) {
          diffwork[i][j] = new double [DifLoop+1];
        }
      }
      cell=new CCell  [X*Y];
      for (int i=0;i<X;i++) {
        for (int j=0;j<Y;j++) {
          cell[i*Y+j]=tis.cell[i*Y+j];
        }
      }
      tmpv=new double[X*Y];
    }
  }
  return(*this);
}


//Don't use with MPI
int Tissue::setDifLoop(int DL) {
  DifLoop=DL;
  return DifLoop;
}

void Tissue::copydata(CCell  *cellx,int x1,int x2, int y1, int y2) {
  if (Y==1) {
    if (x1>=0 && x1<=GX && x2>=0 && x2<=GX) {
      if (x1==GX)x1--;
      if (x2==GX)x2--;
      for (int c=x1;c<=x2;c++)
        cell[c]=*cellx;
    }
    else {
      for (int c=0;c<X;c++)
        cell[c]=*cellx;
    }
  }
  else {
    if (x1>=0 && x1<=GX && x2>=0 && x2<=GX &&
             y1>=0 && y1<=GY && y2>=0 && y2<=GY) {
      if (x1==GX)x1--;
      if (x2==GX)x2--;
      if (y1==GY)y1--;
      if (y2==GY)y2--;
      for (int c=x1;c<=x2;c++) {
        if (my_proc*X<=c && c< (my_proc+1)*X) {
          for (int d=y1;d<=y2;d++) {
            cell[(c-my_proc*X)*Y+d]=*cellx;
          }
        }
      }
    }
    else {
      for (int c=0;c<X;c++) {
        for (int d=0;d<Y;d++) {
          cell[c*Y+d]=*cellx;
        }
      }
    }
  }
}
double Tissue::setdt(double dtt) {
  dt=dtt;
  for (int id=0;id<X*Y;id++)
    cell[id].setdt(dtt);
  return dt;
}

void Tissue::pace(double stim, int pos, int posx ,int posy) {
  //1d code
  if (Y==1) {
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
      switch (pos){
        case STIM_1D_LEFT:
          #ifdef _OPENMP
          #pragma omp for
          #endif
          for(int i=0;i<X;i++) {
            if (i<R2)
              cell[i].pace(stim);
            else
              cell[i].pace();
          }
          break;
        case STIM_1D_RIGHT:
          #ifdef _OPENMP
          #pragma omp for
          #endif
          for(int i=0;i<X;i++) {
            if (i>=X-R2)
              cell[i].pace(stim);
            else
              cell[i].pace();
          }
          break;
        case STIM_1D_POS:
          #ifdef _OPENMP
          #pragma omp for
          #endif
          for(int i=0;i<X;i++) {
            if (abs(i-posx)<R2/2)
              cell[i].pace(stim);
            else
              cell[i].pace();
          }
          break;
      }

      const int DIFFLOOP=2;
      double Dfudtdx2=Dfu*dt/(dx*dx)/DIFFLOOP;
      tmpv[0]=cell[0].v+(cell[1].v+cell[1].v-2*cell[0].v)*Dfudtdx2;
      tmpv[X-1]=cell[X-1].v+(cell[X-2].v+cell[X-2].v-2*cell[X-1].v)*Dfudtdx2;
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++)
        tmpv[i]=cell[i].v+(cell[i-1].v+cell[i+1].v-2*cell[i].v)*Dfudtdx2;
      
      cell[0].v=tmpv[0]+(tmpv[1]+tmpv[1]-2*tmpv[0])*Dfudtdx2;
      cell[X-1].v=tmpv[X-1]+(tmpv[X-2]+tmpv[X-2]-2*tmpv[X-1])*Dfudtdx2;
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++)
        cell[i].v=tmpv[i]+(tmpv[i-1]+tmpv[i+1]-2*tmpv[i])*Dfudtdx2;

      #ifndef VOLTAGE_ONLY
      if (IsCaDiff) {
        double DfuCadtdx2=DfuCa*dt/(dx*dx)/DIFFLOOP;
        tmpv[0]=cell[0].ci+(cell[1].ci+cell[1].ci-2*cell[0].ci)*DfuCadtdx2;
        tmpv[X-1]=cell[X-1].ci+(cell[X-2].ci+cell[X-2].ci-2*cell[X-1].ci)*DfuCadtdx2;
        #ifdef _OPENMP
        #pragma omp for
        #endif
        for (int i=1;i<X-1;i++)
          tmpv[i]=cell[i].ci+(cell[i-1].ci+cell[i+1].ci-2*cell[i].ci)*DfuCadtdx2;
        
        cell[0].ci=tmpv[0]+(tmpv[1]+tmpv[1]-2*tmpv[0])*DfuCadtdx2;
        cell[X-1].ci=tmpv[X-1]+(tmpv[X-2]+tmpv[X-2]-2*tmpv[X-1])*DfuCadtdx2;
        #ifdef _OPENMP
        #pragma omp for
        #endif
        for (int i=1;i<X-1;i++)
          cell[i].ci=tmpv[i]+(tmpv[i-1]+tmpv[i+1]-2*tmpv[i])*DfuCadtdx2;
      }
      #endif
    }
  }
  else {//2D
#ifdef ___USE_MPI
    switch (pos){
      case STIM_2D_CORNER1:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i+my_proc*X)*(i+my_proc*X)+j*j<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_PARALLEL1:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i+my_proc*X)<Stripe)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_PARALLEL2:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if (j<Stripe)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CORNER2:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i+my_proc*X-GX)*(i+my_proc*X-GX)+j*j<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CORNER3:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i+my_proc*X)*(i+my_proc*X)+(j-GY)*(j-GY)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CORNER4:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i+my_proc*X-GX)*(i+my_proc*X-GX)+(j-GY)*(j-GY)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CENTER:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i+my_proc*X-GX/2)*(i+my_proc*X-GX/2)+(j-GY/2)*(j-GY/2)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_POS:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i+my_proc*X-posx)*(i+my_proc*X-posx)+(j-posy)*(j-posy)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
    }
#else //non-MPI code
    switch (pos){
      case STIM_2D_CORNER1:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if (i*i+j*j<R2) cell[i*Y+j].pace(stim);
            else cell[i*Y+j].pace(0);
        break;
      case STIM_2D_PARALLEL1:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if (i<Stripe)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_PARALLEL2:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if (j<Stripe)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CORNER2:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i-X)*(i-X)+j*j<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CORNER3:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if (i*i+(j-Y)*(j-Y)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CORNER4:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i-X)*(i-X)+(j-Y)*(j-Y)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_CENTER:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i-X/2)*(i-X/2)+(j-Y/2)*(j-Y/2)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
      case STIM_2D_POS:
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i=0;i<X;i++)
          for(int j=0;j<Y;j++)
            if ((i-posx)*(i-posx)+(j-posy)*(j-posy)<R2)
              cell[i*Y+j].pace(stim);
            else
              cell[i*Y+j].pace(0);
        break;
    }
#endif

    #ifdef ___FIBER
    diffusionf();
    #else
    diffusion();
    #endif
  }
}

void Tissue::paceall(double stim) {
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int i=0;i<X;i++) {
      cell[i].pace(stim);
    }
    double Dfudtdx2=Dfu*dt/(dx*dx)/2;
  #ifdef ___PERIODIC
    tmpv[0]=cell[0].v+(cell[1].v+cell[X-1].v-2*cell[0].v)*Dfudtdx2;
    tmpv[X-1]=cell[X-1].v+(cell[X-2].v+cell[0].v-2*cell[X-1].v)*Dfudtdx2;
  #else
    tmpv[0]=cell[0].v+(cell[1].v+cell[1].v-2*cell[0].v)*Dfudtdx2;
    tmpv[X-1]=cell[X-1].v+(cell[X-2].v+cell[X-2].v-2*cell[X-1].v)*Dfudtdx2;
  #endif
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (int i=1;i<X-1;i++)
      tmpv[i]=cell[i].v+(cell[i-1].v+cell[i+1].v-2*cell[i].v)*Dfudtdx2;
    
  #ifdef ___PERIODIC
    cell[0].v=tmpv[0]+(tmpv[1]+tmpv[X-1]-2*tmpv[0])*Dfudtdx2;
    cell[X-1].v=tmpv[X-1]+(tmpv[X-2]+tmpv[0]-2*tmpv[X-1])*Dfudtdx2;
  #else
    cell[0].v=tmpv[0]+(tmpv[1]+tmpv[1]-2*tmpv[0])*Dfudtdx2;
    cell[X-1].v=tmpv[X-1]+(tmpv[X-2]+tmpv[X-2]-2*tmpv[X-1])*Dfudtdx2;
  #endif
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (int i=1;i<X-1;i++)
      cell[i].v=tmpv[i]+(tmpv[i-1]+tmpv[i+1]-2*tmpv[i])*Dfudtdx2;

    #ifndef VOLTAGE_ONLY
    if (IsCaDiff) {
      double DfuCadtdx2=DfuCa*dt/(dx*dx)/2;
      #ifdef ___PERIODIC
      tmpv[0]=cell[0].ci+(cell[1].ci+cell[X-1].ci-2*cell[0].ci)*DfuCadtdx2;
      tmpv[X-1]=cell[X-1].ci+(cell[X-2].ci+cell[0].ci-2*cell[X-1].ci)*DfuCadtdx2;
      #else
      tmpv[0]=cell[0].ci+(cell[1].ci+cell[1].ci-2*cell[0].ci)*DfuCadtdx2;
      tmpv[X-1]=cell[X-1].ci+(cell[X-2].ci+cell[X-2].ci-2*cell[X-1].ci)*DfuCadtdx2;
      #endif
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++)
        tmpv[i]=cell[i].ci+(cell[i-1].ci+cell[i+1].ci-2*cell[i].ci)*DfuCadtdx2;
      
      #ifdef ___PERIODIC
      cell[0].ci=tmpv[0]+(tmpv[1]+tmpv[X-1]-2*tmpv[0])*DfuCadtdx2;
      cell[X-1].ci=tmpv[X-1]+(tmpv[X-2]+tmpv[0]-2*tmpv[X-1])*DfuCadtdx2;
      #else
      cell[0].ci=tmpv[0]+(tmpv[1]+tmpv[1]-2*tmpv[0])*DfuCadtdx2;
      cell[X-1].ci=tmpv[X-1]+(tmpv[X-2]+tmpv[X-2]-2*tmpv[X-1])*DfuCadtdx2;
      #endif
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++)
        cell[i].ci=tmpv[i]+(tmpv[i-1]+tmpv[i+1]-2*tmpv[i])*DfuCadtdx2;
    }
    #endif
  }
}
