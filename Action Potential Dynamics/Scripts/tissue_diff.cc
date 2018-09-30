void Tissue::diffusion(void) {
#ifdef ___USE_MPI
  double Dfudtdx2=Dfu*dt/(dx*dx*DifLoop);
  for (int lp=0;lp<DifLoop/2;lp++) {
    if (my_proc==0) { //1st node
      for (int j=0;j<Y;j++) {
        mpiwork_send2[j]=cell[(X-1)*Y+j].v;
      }
      MPI_Status status;
      MPI_Sendrecv(mpiwork_send2,Y,MPI_DOUBLE,1,0,mpiwork_recv2,Y,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&status);

      tmpv[0*Y+0]=cell[0*Y+0].v+(cell[1*Y+0].v+cell[1*Y+0].v+cell[0*Y+1].v+cell[0*Y+1].v-4*cell[0*Y+0].v)*Dfudtdx2;
      tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+(cell[1*Y+Y-1].v+cell[1*Y+Y-1].v+cell[0*Y+Y-2].v+cell[0*Y+Y-2].v-4*cell[0*Y+Y-1].v)*Dfudtdx2;

      tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+(cell[(X-2)*Y+0].v+mpiwork_recv2[0]+cell[(X-1)*Y+1].v+cell[(X-1)*Y+1].v-4*cell[(X-1)*Y+0].v)*Dfudtdx2;
      tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+(cell[(X-2)*Y+Y-1].v+mpiwork_recv2[Y-1]+cell[(X-1)*Y+Y-2].v+cell[(X-1)*Y+Y-2].v-4*cell[(X-1)*Y+Y-1].v)*Dfudtdx2;

      for (int i=1;i<Y-1;i++) {
        tmpv[0*Y+i]=cell[0*Y+i].v+(cell[0*Y+i-1].v+cell[0*Y+i+1].v+cell[1*Y+i].v+cell[1*Y+i].v-4*cell[0*Y+i].v)*Dfudtdx2;
        tmpv[(X-1)*Y+i]=cell[(X-1)*Y+i].v+(cell[(X-1)*Y+i-1].v+cell[(X-1)*Y+i+1].v+cell[(X-2)*Y+i].v+mpiwork_recv2[i]-4*cell[(X-1)*Y+i].v)*Dfudtdx2;
      }
    }
    else if (my_proc==num_procs-1) { //last node
      for (int j=0;j<Y;j++) {
        mpiwork_send1[j]=cell[0*Y+j].v;
      }
      MPI_Status status;
      MPI_Sendrecv(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);

      tmpv[0*Y+0]=cell[0*Y+0].v+(mpiwork_recv1[0]+cell[1*Y+0].v+cell[0*Y+1].v+cell[0*Y+1].v-4*cell[0*Y+0].v)*Dfudtdx2;
      tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+(mpiwork_recv1[Y-1]+cell[1*Y+Y-1].v+cell[0*Y+Y-2].v+cell[0*Y+Y-2].v-4*cell[0*Y+Y-1].v)*Dfudtdx2;

      tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+(cell[(X-2)*Y+0].v+cell[(X-2)*Y+0].v+cell[(X-1)*Y+1].v+cell[(X-1)*Y+1].v-4*cell[(X-1)*Y+0].v)*Dfudtdx2;
      tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+(cell[(X-2)*Y+Y-1].v+cell[(X-2)*Y+Y-1].v+cell[(X-1)*Y+Y-2].v+cell[(X-1)*Y+Y-2].v-4*cell[(X-1)*Y+Y-1].v)*Dfudtdx2;

      for (int i=1;i<Y-1;i++) {
        tmpv[0*Y+i]=cell[0*Y+i].v+(cell[0*Y+i-1].v+cell[0*Y+i+1].v+mpiwork_recv1[i]+cell[1*Y+i].v-4*cell[0*Y+i].v)*Dfudtdx2;
        tmpv[(X-1)*Y+i]=cell[(X-1)*Y+i].v+(cell[(X-1)*Y+i-1].v+cell[(X-1)*Y+i+1].v+cell[(X-2)*Y+i].v+cell[(X-2)*Y+i].v-4*cell[(X-1)*Y+i].v)*Dfudtdx2;
      }
    }
    else {
      for (int j=0;j<Y;j++) { //the other nodes
        mpiwork_send1[j]=cell[0*Y+j].v;
        mpiwork_send2[j]=cell[(X-1)*Y+j].v;
      }
      MPI_Status status;
      MPI_Sendrecv(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
      MPI_Sendrecv(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);

      tmpv[0*Y+0]=cell[0*Y+0].v+(mpiwork_recv1[0]+cell[1*Y+0].v+cell[0*Y+1].v+cell[0*Y+1].v-4*cell[0*Y+0].v)*Dfudtdx2;
      tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+(mpiwork_recv1[Y-1]+cell[1*Y+Y-1].v+cell[0*Y+Y-2].v+cell[0*Y+Y-2].v-4*cell[0*Y+Y-1].v)*Dfudtdx2;
      tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+(cell[(X-2)*Y+0].v+mpiwork_recv2[0]+cell[(X-1)*Y+1].v+cell[(X-1)*Y+1].v-4*cell[(X-1)*Y+0].v)*Dfudtdx2;
      tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+(cell[(X-2)*Y+Y-1].v+mpiwork_recv2[Y-1]+cell[(X-1)*Y+Y-2].v+cell[(X-1)*Y+Y-2].v-4*cell[(X-1)*Y+Y-1].v)*Dfudtdx2;
      for (int i=1;i<Y-1;i++) {
        tmpv[0*Y+i]=cell[0*Y+i].v+(cell[0*Y+i-1].v+cell[0*Y+i+1].v+mpiwork_recv1[i]+cell[1*Y+i].v-4*cell[0*Y+i].v)*Dfudtdx2;
        tmpv[(X-1)*Y+i]=cell[(X-1)*Y+i].v+(cell[(X-1)*Y+i-1].v+cell[(X-1)*Y+i+1].v+cell[(X-2)*Y+i].v+mpiwork_recv2[i]-4*cell[(X-1)*Y+i].v)*Dfudtdx2;
      }
    }
    // all nodes (inside tissue)
    for (int i=1;i<X-1;i++) {
      tmpv[i*Y+0]=cell[i*Y+0].v+(cell[(i-1)*Y+0].v+cell[(i+1)*Y+0].v+cell[i*Y+1].v+cell[i*Y+1].v-4*cell[i*Y+0].v)*Dfudtdx2;
      tmpv[i*Y+Y-1]=cell[i*Y+Y-1].v+(cell[(i-1)*Y+Y-1].v+cell[(i+1)*Y+Y-1].v+cell[i*Y+Y-2].v+cell[i*Y+Y-2].v-4*cell[i*Y+Y-1].v)*Dfudtdx2;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=1;i<X-1;i++)
      for (int j=1;j<Y-1;j++)
        tmpv[i*Y+j]=cell[i*Y+j].v+(cell[(i-1)*Y+j].v+cell[(i+1)*Y+j].v+cell[i*Y+j+1].v+cell[i*Y+j-1].v-4*cell[i*Y+j].v)*Dfudtdx2;

    if (my_proc==0) {//1st node
      for (int j=0;j<Y;j++) {
        mpiwork_send2[j]=tmpv[(X-1)*Y+j];
      }
      MPI_Status status;
      MPI_Sendrecv(mpiwork_send2,Y,MPI_DOUBLE,1,0,mpiwork_recv2,Y,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&status);

      cell[0*Y+0].v=tmpv[0*Y+0]+(tmpv[1*Y+0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtdx2;
      cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+(tmpv[1*Y+Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtdx2;

      cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+mpiwork_recv2[0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtdx2;
      cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtdx2;


      for (int i=1;i<Y-1;i++) {
        cell[0*Y+i].v=tmpv[0*Y+i]+(tmpv[0*Y+i-1]+tmpv[0*Y+i+1]+tmpv[1*Y+i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtdx2;
        cell[(X-1)*Y+i].v=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+i-1]+tmpv[(X-1)*Y+i+1]+tmpv[(X-2)*Y+i]+mpiwork_recv2[i]-4*tmpv[(X-1)*Y+i])*Dfudtdx2;
      }
    }
    else if (my_proc==num_procs-1) {//last node
      for (int j=0;j<Y;j++) {
        mpiwork_send1[j]=tmpv[0*Y+j];
      }
      MPI_Status status;
      MPI_Sendrecv(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);

      cell[0*Y+0].v=tmpv[0*Y+0]+(mpiwork_recv1[0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtdx2;
      cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+(mpiwork_recv1[Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtdx2;

      cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+tmpv[(X-2)*Y+0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtdx2;
      cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+tmpv[(X-2)*Y+Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtdx2;
      for (int i=1;i<Y-1;i++) {
        cell[0*Y+i].v=tmpv[0*Y+i]+(tmpv[0*Y+i-1]+tmpv[0*Y+i+1]+mpiwork_recv1[i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtdx2;
        cell[(X-1)*Y+i].v=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+i-1]+tmpv[(X-1)*Y+i+1]+tmpv[(X-2)*Y+i]+tmpv[(X-2)*Y+i]-4*tmpv[(X-1)*Y+i])*Dfudtdx2;
      }
    }
    else {//the other nodes
      for (int j=0;j<Y;j++) {
        mpiwork_send1[j]=tmpv[0*Y+j];
        mpiwork_send2[j]=tmpv[(X-1)*Y+j];
      }
      MPI_Status status;
      MPI_Sendrecv(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
      MPI_Sendrecv(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);

      cell[0*Y+0].v=tmpv[0*Y+0]+(mpiwork_recv1[0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtdx2;
      cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+(mpiwork_recv1[Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtdx2;
      cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+mpiwork_recv2[0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtdx2;
      cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtdx2;
      for (int i=1;i<Y-1;i++)
      {
        cell[0*Y+i].v=tmpv[0*Y+i]+(tmpv[0*Y+i-1]+tmpv[0*Y+i+1]+mpiwork_recv1[i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtdx2;
        cell[(X-1)*Y+i].v=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+i-1]+tmpv[(X-1)*Y+i+1]+tmpv[(X-2)*Y+i]+mpiwork_recv2[i]-4*tmpv[(X-1)*Y+i])*Dfudtdx2;
      }
    }
    //all nodes
    for (int i=1;i<X-1;i++) {
      cell[i*Y+0].v=tmpv[i*Y+0]+(tmpv[(i-1)*Y+0]+tmpv[(i+1)*Y+0]+tmpv[i*Y+1]+tmpv[i*Y+1]-4*tmpv[i*Y+0])*Dfudtdx2;
      cell[i*Y+Y-1].v=tmpv[i*Y+Y-1]+(tmpv[(i-1)*Y+Y-1]+tmpv[(i+1)*Y+Y-1]+tmpv[i*Y+Y-2]+tmpv[i*Y+Y-2]-4*tmpv[i*Y+Y-1])*Dfudtdx2;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=1;i<X-1;i++)
      for (int j=1;j<Y-1;j++)
        cell[i*Y+j].v=tmpv[i*Y+j]+(tmpv[(i-1)*Y+j]+tmpv[(i+1)*Y+j]+tmpv[i*Y+j+1]+tmpv[i*Y+j-1]-4*tmpv[i*Y+j])*Dfudtdx2;
  }
#else //non-MPI code
  double Dfudtdx2=Dfu*dt/(dx*dx*DifLoop);
  for (int lp=0;lp<DifLoop/2;lp++) {
    //corners
    tmpv[0*Y+0]=cell[0*Y+0].v+(cell[1*Y+0].v+cell[1*Y+0].v+cell[0*Y+1].v+cell[0*Y+1].v-4*cell[0*Y+0].v)*Dfudtdx2;
    tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+(cell[(X-2)*Y+0].v+cell[(X-2)*Y+0].v+cell[(X-1)*Y+1].v+cell[(X-1)*Y+1].v-4*cell[(X-1)*Y+0].v)*Dfudtdx2;
    tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+(cell[1*Y+Y-1].v+cell[1*Y+Y-1].v+cell[0*Y+Y-2].v+cell[0*Y+Y-2].v-4*cell[0*Y+Y-1].v)*Dfudtdx2;
    tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+(cell[(X-2)*Y+Y-1].v+cell[(X-2)*Y+Y-1].v+cell[(X-1)*Y+Y-2].v+cell[(X-1)*Y+Y-2].v-4*cell[(X-1)*Y+Y-1].v)*Dfudtdx2;
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++) {
        tmpv[i*Y+0]=cell[i*Y+0].v+(cell[(i-1)*Y+0].v+cell[(i+1)*Y+0].v+cell[i*Y+1].v+cell[i*Y+1].v-4*cell[i*Y+0].v)*Dfudtdx2;
        tmpv[i*Y+Y-1]=cell[i*Y+Y-1].v+(cell[(i-1)*Y+Y-1].v+cell[(i+1)*Y+Y-1].v+cell[i*Y+Y-2].v+cell[i*Y+Y-2].v-4*cell[i*Y+Y-1].v)*Dfudtdx2;
      }
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<Y-1;i++) {
        tmpv[0*Y+i]=cell[0*Y+i].v+(cell[0*Y+i-1].v+cell[0*Y+i+1].v+cell[1*Y+i].v+cell[1*Y+i].v-4*cell[0*Y+i].v)*Dfudtdx2;
        tmpv[(X-1)*Y+i]=cell[(X-1)*Y+i].v+(cell[(X-1)*Y+i-1].v+cell[(X-1)*Y+i+1].v+cell[(X-2)*Y+i].v+cell[(X-2)*Y+i].v-4*cell[(X-1)*Y+i].v)*Dfudtdx2;
      }
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++)
        for (int j=1;j<Y-1;j++)
          tmpv[i*Y+j]=cell[i*Y+j].v+(cell[(i-1)*Y+j].v+cell[(i+1)*Y+j].v+cell[i*Y+j+1].v+cell[i*Y+j-1].v-4*cell[i*Y+j].v)*Dfudtdx2;

      //corners
      cell[0*Y+0].v=tmpv[0*Y+0]+(tmpv[1*Y+0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtdx2;
      cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+tmpv[(X-2)*Y+0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtdx2;
      cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+(tmpv[1*Y+Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtdx2;
      cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+tmpv[(X-2)*Y+Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtdx2;
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++) {
        cell[i*Y+0].v=tmpv[i*Y+0]+(tmpv[(i-1)*Y+0]+tmpv[(i+1)*Y+0]+tmpv[i*Y+1]+tmpv[i*Y+1]-4*tmpv[i*Y+0])*Dfudtdx2;
        cell[i*Y+Y-1].v=tmpv[i*Y+Y-1]+(tmpv[(i-1)*Y+Y-1]+tmpv[(i+1)*Y+Y-1]+tmpv[i*Y+Y-2]+tmpv[i*Y+Y-2]-4*tmpv[i*Y+Y-1])*Dfudtdx2;
      }
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<Y-1;i++) {
        cell[0*Y+i].v=tmpv[0*Y+i]+(tmpv[0*Y+i-1]+tmpv[0*Y+i+1]+tmpv[1*Y+i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtdx2;
        cell[(X-1)*Y+i].v=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+i-1]+tmpv[(X-1)*Y+i+1]+tmpv[(X-2)*Y+i]+tmpv[(X-2)*Y+i]-4*tmpv[(X-1)*Y+i])*Dfudtdx2;
      }
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for (int i=1;i<X-1;i++)
        for (int j=1;j<Y-1;j++)
          cell[i*Y+j].v=tmpv[i*Y+j]+(tmpv[(i-1)*Y+j]+tmpv[(i+1)*Y+j]+tmpv[i*Y+j+1]+tmpv[i*Y+j-1]-4*tmpv[i*Y+j])*Dfudtdx2;
    }
  }

#endif
}
//diffusion with fiber orientation
void Tissue::diffusionf(void) {
#ifdef ___USE_MPI
  if (my_proc==0) {//1st node
    #pragma ivdep
    #pragma vector always
    for (int j=0;j<Y;j++) {
      mpiwork_send2[j]=Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j];
    }
    MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
    MPI_Status status;
    MPI_Recv(Dxxphi2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
    #pragma ivdep
    #pragma vector always
    for (int j=0;j<Y;j++) {
      mpiwork_send2[j]=Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j];
    }
    MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
    MPI_Recv(Dxyphi2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
  }
  else if (my_proc==num_procs-1) { //last node
    #pragma ivdep
    #pragma vector always
    for (int j=0;j<Y;j++) {
      mpiwork_send1[j]=Dxx[0*Y+j]*phi[0*Y+j];
    }
    MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
    MPI_Status status;
    MPI_Recv(Dxxphi1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
    #pragma ivdep
    #pragma vector always
    for (int j=0;j<Y;j++)
    {
      mpiwork_send1[j]=Dxy[0*Y+j]*phi[0*Y+j];
    }
    MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
    MPI_Recv(Dxyphi1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
  }
  else { //the other nodes
    #pragma ivdep
    #pragma vector always
    for (int j=0;j<Y;j++) {
      mpiwork_send1[j]=Dxx[0*Y+j]*phi[0*Y+j];
      mpiwork_send2[j]=Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j];
    }
    MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
    MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
    MPI_Status status;
    MPI_Recv(Dxxphi1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
    MPI_Recv(Dxxphi2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
    #pragma ivdep
    #pragma vector always
    for (int j=0;j<Y;j++) {
      mpiwork_send1[j]=Dxy[0*Y+j]*phi[0*Y+j];
      mpiwork_send2[j]=Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j];
    }
    MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
    MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
    MPI_Recv(Dxyphi1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
    MPI_Recv(Dxyphi2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
  }
  double Dfudtdx2=Dfu*dt/(dx*dx*DifLoop);
  for (int i=0;i<X;i++) {
    #pragma ivdep
    #pragma vector always
    for (int j=0;j<Y;j++) {
      tmpv[i*Y+j]=cell[i*Y+j].v;
    }
  }
  for (int lp=0;lp<DifLoop/2;lp++) {
    if (my_proc==0) { //1st node
      #pragma ivdep
      #pragma vector always
      for (int j=0;j<Y;j++) {
        mpiwork_send2[j]=cell[(X-1)*Y+j].v;
      }
      MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
      MPI_Status status;
      MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);

      tmpv[0*Y+0]=cell[0*Y+0].v+(cell[1*Y+0].v+cell[1*Y+0].v+cell[0*Y+1].v+cell[0*Y+1].v-4*cell[0*Y+0].v)*Dfudtdx2;
      tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+(cell[1*Y+Y-1].v+cell[1*Y+Y-1].v+cell[0*Y+Y-2].v+cell[0*Y+Y-2].v-4*cell[0*Y+Y-1].v)*Dfudtdx2;

      tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+(cell[(X-2)*Y+0].v+mpiwork_recv2[0]+cell[(X-1)*Y+1].v+cell[(X-1)*Y+1].v-4*cell[(X-1)*Y+0].v)*Dfudtdx2;
      tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+(cell[(X-2)*Y+Y-1].v+mpiwork_recv2[Y-1]+cell[(X-1)*Y+Y-2].v+cell[(X-1)*Y+Y-2].v-4*cell[(X-1)*Y+Y-1].v)*Dfudtdx2;

      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        tmpv[0*Y+j]=cell[0*Y+j].v+(cell[0*Y+j-1].v+cell[0*Y+j+1].v+cell[1*Y+j].v+cell[1*Y+j].v-4*cell[0*Y+j].v)*Dfudtdx2;
        if (phi[(X-1)*Y+j]>cutoff)
        {
          double ldv=(cell[(X-1-1)*Y+j].v+cell[(X-1-1)*Y+(j-1)].v+cell[(X-1)*Y+(j-1)].v+cell[(X-1)*Y+j].v)/4;
          double luv=(cell[(X-1-1)*Y+j].v+cell[(X-1-1)*Y+(j+1)].v+cell[(X-1)*Y+(j+1)].v+cell[(X-1)*Y+j].v)/4;
          double rdv=(mpiwork_recv2[j]+mpiwork_recv2[j-1]+cell[(X-1)*Y+(j-1)].v+cell[(X-1)*Y+j].v)/4;
          double ruv=(mpiwork_recv2[j]+mpiwork_recv2[j+1]+cell[(X-1)*Y+(j+1)].v+cell[(X-1)*Y+j].v)/4;

          double J11p=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxxphi2[j])/2*(mpiwork_recv2[j]-cell[(X-1)*Y+j].v)/dx;
          double J11m=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxx[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(cell[(X-1)*Y+j].v-cell[(X-1-1)*Y+j].v)/dx;

          double J12p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxyphi2[j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(luv-ldv)/dx;

          double J21p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j+1]*phi[(X-1)*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j-1]*phi[(X-1)*Y+j-1])/2*(rdv-ldv)/dx;

          double J22p=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j+1)]*phi[(X-1)*Y+(j+1)])/2*(cell[(X-1)*Y+(j+1)].v-cell[(X-1)*Y+j].v)/dx;
          double J22m=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j-1)]*phi[(X-1)*Y+(j-1)])/2*(cell[(X-1)*Y+j].v-cell[(X-1)*Y+(j-1)].v)/dx;

          tmpv[(X-1)*Y+j]=cell[(X-1)*Y+j].v+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[(X-1)*Y+j]);
        }
      }
    }
    else if (my_proc==num_procs-1) { //last node
      #pragma ivdep
      #pragma vector always
      for (int j=0;j<Y;j++) {
        mpiwork_send1[j]=cell[0*Y+j].v;
      }
      MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
      MPI_Status status;
      MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);

      tmpv[0*Y+0]=cell[0*Y+0].v+(mpiwork_recv1[0]+cell[1*Y+0].v+cell[0*Y+1].v+cell[0*Y+1].v-4*cell[0*Y+0].v)*Dfudtdx2;
      tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+(mpiwork_recv1[Y-1]+cell[1*Y+Y-1].v+cell[0*Y+Y-2].v+cell[0*Y+Y-2].v-4*cell[0*Y+Y-1].v)*Dfudtdx2;

      tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+(cell[(X-2)*Y+0].v+cell[(X-2)*Y+0].v+cell[(X-1)*Y+1].v+cell[(X-1)*Y+1].v-4*cell[(X-1)*Y+0].v)*Dfudtdx2;
      tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+(cell[(X-2)*Y+Y-1].v+cell[(X-2)*Y+Y-1].v+cell[(X-1)*Y+Y-2].v+cell[(X-1)*Y+Y-2].v-4*cell[(X-1)*Y+Y-1].v)*Dfudtdx2;

      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        tmpv[(X-1)*Y+j]=cell[(X-1)*Y+j].v+(cell[(X-1)*Y+j-1].v+cell[(X-1)*Y+j+1].v+cell[(X-2)*Y+j].v+cell[(X-2)*Y+j].v-4*cell[(X-1)*Y+j].v)*Dfudtdx2;
        if (phi[0*Y+j]>cutoff) {
          double ldv=(mpiwork_recv1[j]+mpiwork_recv1[j-1]+cell[0*Y+(j-1)].v+cell[0*Y+j].v)/4;
          double luv=(mpiwork_recv1[j]+mpiwork_recv1[j+1]+cell[0*Y+(j+1)].v+cell[0*Y+j].v)/4;
          double rdv=(cell[(0+1)*Y+j].v+cell[(0+1)*Y+(j-1)].v+cell[0*Y+(j-1)].v+cell[0*Y+j].v)/4;
          double ruv=(cell[(0+1)*Y+j].v+cell[(0+1)*Y+(j+1)].v+cell[0*Y+(j+1)].v+cell[0*Y+j].v)/4;

          double J11p=(Dxx[0*Y+j]*phi[0*Y+j]+Dxx[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(cell[(0+1)*Y+j].v-cell[0*Y+j].v)/dx;
          double J11m=(Dxx[0*Y+j]*phi[0*Y+j]+Dxxphi1[j])/2*(cell[0*Y+j].v-mpiwork_recv1[j])/dx;
          double J12p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxyphi1[j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j+1]*phi[0*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j-1]*phi[0*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j+1)]*phi[0*Y+(j+1)])/2*(cell[0*Y+(j+1)].v-cell[0*Y+j].v)/dx;
          double J22m=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j-1)]*phi[0*Y+(j-1)])/2*(cell[0*Y+j].v-cell[0*Y+(j-1)].v)/dx;
          tmpv[0*Y+j]=cell[0*Y+j].v+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[0*Y+j]);
        }
      }
    }
    else {
      #pragma ivdep
      #pragma vector always
      for (int j=0;j<Y;j++) {
        mpiwork_send1[j]=cell[0*Y+j].v;
        mpiwork_send2[j]=cell[(X-1)*Y+j].v;
      }
      MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
      MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
      MPI_Status status;
      MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
      MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
      tmpv[0*Y+0]=cell[0*Y+0].v+(mpiwork_recv1[0]+cell[1*Y+0].v+cell[0*Y+1].v+cell[0*Y+1].v-4*cell[0*Y+0].v)*Dfudtdx2;
      tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+(mpiwork_recv1[Y-1]+cell[1*Y+Y-1].v+cell[0*Y+Y-2].v+cell[0*Y+Y-2].v-4*cell[0*Y+Y-1].v)*Dfudtdx2;
      tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+(cell[(X-2)*Y+0].v+mpiwork_recv2[0]+cell[(X-1)*Y+1].v+cell[(X-1)*Y+1].v-4*cell[(X-1)*Y+0].v)*Dfudtdx2;
      tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+(cell[(X-2)*Y+Y-1].v+mpiwork_recv2[Y-1]+cell[(X-1)*Y+Y-2].v+cell[(X-1)*Y+Y-2].v-4*cell[(X-1)*Y+Y-1].v)*Dfudtdx2;
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        if (phi[0*Y+j]>cutoff) {
          double ldv=(mpiwork_recv1[j]+mpiwork_recv1[j-1]+cell[0*Y+(j-1)].v+cell[0*Y+j].v)/4;
          double luv=(mpiwork_recv1[j]+mpiwork_recv1[j+1]+cell[0*Y+(j+1)].v+cell[0*Y+j].v)/4;
          double rdv=(cell[(0+1)*Y+j].v+cell[(0+1)*Y+(j-1)].v+cell[0*Y+(j-1)].v+cell[0*Y+j].v)/4;
          double ruv=(cell[(0+1)*Y+j].v+cell[(0+1)*Y+(j+1)].v+cell[0*Y+(j+1)].v+cell[0*Y+j].v)/4;

          double J11p=(Dxx[0*Y+j]*phi[0*Y+j]+Dxx[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(cell[(0+1)*Y+j].v-cell[0*Y+j].v)/dx;
          double J11m=(Dxx[0*Y+j]*phi[0*Y+j]+Dxxphi1[j])/2*(cell[0*Y+j].v-mpiwork_recv1[j])/dx;
          double J12p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxyphi1[j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j+1]*phi[0*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j-1]*phi[0*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j+1)]*phi[0*Y+(j+1)])/2*(cell[0*Y+(j+1)].v-cell[0*Y+j].v)/dx;
          double J22m=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j-1)]*phi[0*Y+(j-1)])/2*(cell[0*Y+j].v-cell[0*Y+(j-1)].v)/dx;

          tmpv[0*Y+j]=cell[0*Y+j].v+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[0*Y+j]);
        }
        if (phi[(X-1)*Y+j]>cutoff) {
          double ldv=(cell[(X-1-1)*Y+j].v+cell[(X-1-1)*Y+(j-1)].v+cell[(X-1)*Y+(j-1)].v+cell[(X-1)*Y+j].v)/4;
          double luv=(cell[(X-1-1)*Y+j].v+cell[(X-1-1)*Y+(j+1)].v+cell[(X-1)*Y+(j+1)].v+cell[(X-1)*Y+j].v)/4;
          double rdv=(mpiwork_recv2[j]+mpiwork_recv2[j-1]+cell[(X-1)*Y+(j-1)].v+cell[(X-1)*Y+j].v)/4;
          double ruv=(mpiwork_recv2[j]+mpiwork_recv2[j+1]+cell[(X-1)*Y+(j+1)].v+cell[(X-1)*Y+j].v)/4;

          double J11p=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxxphi2[j])/2*(mpiwork_recv2[j]-cell[(X-1)*Y+j].v)/dx;
          double J11m=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxx[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(cell[(X-1)*Y+j].v-cell[(X-1-1)*Y+j].v)/dx;
          double J12p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxyphi2[j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j+1]*phi[(X-1)*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j-1]*phi[(X-1)*Y+j-1])/2*(rdv-ldv)/dx;

          double J22p=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j+1)]*phi[(X-1)*Y+(j+1)])/2*(cell[(X-1)*Y+(j+1)].v-cell[(X-1)*Y+j].v)/dx;
          double J22m=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j-1)]*phi[(X-1)*Y+(j-1)])/2*(cell[(X-1)*Y+j].v-cell[(X-1)*Y+(j-1)].v)/dx;

          tmpv[(X-1)*Y+j]=cell[(X-1)*Y+j].v+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[(X-1)*Y+j]);
        }
      }
    }
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<X-1;i++) {
      tmpv[i*Y+0]=cell[i*Y+0].v+(cell[(i-1)*Y+0].v+cell[(i+1)*Y+0].v+cell[i*Y+1].v+cell[i*Y+1].v-4*cell[i*Y+0].v)*Dfudtdx2;
      tmpv[i*Y+Y-1]=cell[i*Y+Y-1].v+(cell[(i-1)*Y+Y-1].v+cell[(i+1)*Y+Y-1].v+cell[i*Y+Y-2].v+cell[i*Y+Y-2].v-4*cell[i*Y+Y-1].v)*Dfudtdx2;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=1;i<X-1;i++) {
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        if (phi[i*Y+j]>cutoff) {
          double ldv=(cell[(i-1)*Y+j].v+cell[(i-1)*Y+(j-1)].v+cell[i*Y+(j-1)].v+cell[i*Y+j].v)/4;
          double luv=(cell[(i-1)*Y+j].v+cell[(i-1)*Y+(j+1)].v+cell[i*Y+(j+1)].v+cell[i*Y+j].v)/4;
          double rdv=(cell[(i+1)*Y+j].v+cell[(i+1)*Y+(j-1)].v+cell[i*Y+(j-1)].v+cell[i*Y+j].v)/4;
          double ruv=(cell[(i+1)*Y+j].v+cell[(i+1)*Y+(j+1)].v+cell[i*Y+(j+1)].v+cell[i*Y+j].v)/4;

          double J11p=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(cell[(i+1)*Y+j].v-cell[i*Y+j].v)/dx;
          double J11m=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(cell[i*Y+j].v-cell[(i-1)*Y+j].v)/dx;
          double J12p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j+1]*phi[i*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j-1]*phi[i*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j+1)]*phi[i*Y+(j+1)])/2*(cell[i*Y+(j+1)].v-cell[i*Y+j].v)/dx;
          double J22m=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j-1)]*phi[i*Y+(j-1)])/2*(cell[i*Y+j].v-cell[i*Y+(j-1)].v)/dx;

          tmpv[i*Y+j]=cell[i*Y+j].v+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[i*Y+j]);
        }
      }
    }
//
    if (my_proc==0) {
      #pragma ivdep
      #pragma vector always
      for (int j=0;j<Y;j++) {
        mpiwork_send2[j]=tmpv[(X-1)*Y+j];
      }
      MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
      MPI_Status status;
      MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);

      cell[0*Y+0].v=tmpv[0*Y+0]+(tmpv[1*Y+0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtdx2;
      cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+(tmpv[1*Y+Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtdx2;

      cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+mpiwork_recv2[0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtdx2;
      cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtdx2;


      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        cell[0*Y+j].v=tmpv[0*Y+j]+(tmpv[0*Y+j-1]+tmpv[0*Y+j+1]+tmpv[1*Y+j]+tmpv[1*Y+j]-4*tmpv[0*Y+j])*Dfudtdx2;
        if (phi[(X-1)*Y+j]>cutoff) {
          double ldv=(tmpv[(X-1-1)*Y+j]+tmpv[(X-1-1)*Y+(j-1)]+tmpv[(X-1)*Y+(j-1)]+tmpv[(X-1)*Y+j])/4;
          double luv=(tmpv[(X-1-1)*Y+j]+tmpv[(X-1-1)*Y+(j+1)]+tmpv[(X-1)*Y+(j+1)]+tmpv[(X-1)*Y+j])/4;
          double rdv=(mpiwork_recv2[j]+mpiwork_recv2[j-1]+tmpv[(X-1)*Y+(j-1)]+tmpv[(X-1)*Y+j])/4;
          double ruv=(mpiwork_recv2[j]+mpiwork_recv2[j+1]+tmpv[(X-1)*Y+(j+1)]+tmpv[(X-1)*Y+j])/4;

          double J11p=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxxphi2[j])/2*(mpiwork_recv2[j]-tmpv[(X-1)*Y+j])/dx;
          double J11m=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxx[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(tmpv[(X-1)*Y+j]-tmpv[(X-1-1)*Y+j])/dx;
          double J12p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxyphi2[j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j+1]*phi[(X-1)*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j-1]*phi[(X-1)*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j+1)]*phi[(X-1)*Y+(j+1)])/2*(tmpv[(X-1)*Y+(j+1)]-tmpv[(X-1)*Y+j])/dx;
          double J22m=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j-1)]*phi[(X-1)*Y+(j-1)])/2*(tmpv[(X-1)*Y+j]-tmpv[(X-1)*Y+(j-1)])/dx;

          cell[(X-1)*Y+j].v=tmpv[(X-1)*Y+j]+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[(X-1)*Y+j]);
        }
      }
    }
    else if (my_proc==num_procs-1) {
      #pragma ivdep
      #pragma vector always
      for (int j=0;j<Y;j++) {
        mpiwork_send1[j]=tmpv[0*Y+j];
      }
      MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
      MPI_Status status;
      MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);

      cell[0*Y+0].v=tmpv[0*Y+0]+(mpiwork_recv1[0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtdx2;
      cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+(mpiwork_recv1[Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtdx2;
      cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+(tmpv[X-2*Y+0]+tmpv[X-2*Y+0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtdx2;
      cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+(tmpv[X-2*Y+Y-1]+tmpv[X-2*Y+Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtdx2;
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        cell[(X-1)*Y+j].v=tmpv[(X-1)*Y+j]+(tmpv[(X-1)*Y+j-1]+tmpv[(X-1)*Y+j+1]+tmpv[(X-2)*Y+j]+tmpv[(X-2)*Y+j]-4*tmpv[(X-1)*Y+j])*Dfudtdx2;
        if (phi[0*Y+j]>cutoff) {
          double ldv=(mpiwork_recv1[j]+mpiwork_recv1[j-1]+tmpv[0*Y+(j-1)]+tmpv[0*Y+j])/4;
          double luv=(mpiwork_recv1[j]+mpiwork_recv1[j+1]+tmpv[0*Y+(j+1)]+tmpv[0*Y+j])/4;
          double rdv=(tmpv[(0+1)*Y+j]+tmpv[(0+1)*Y+(j-1)]+tmpv[0*Y+(j-1)]+tmpv[0*Y+j])/4;
          double ruv=(tmpv[(0+1)*Y+j]+tmpv[(0+1)*Y+(j+1)]+tmpv[0*Y+(j+1)]+tmpv[0*Y+j])/4;

          double J11p=(Dxx[0*Y+j]*phi[0*Y+j]+Dxx[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(tmpv[(0+1)*Y+j]-tmpv[0*Y+j])/dx;
          double J11m=(Dxx[0*Y+j]*phi[0*Y+j]+Dxxphi1[j])/2*(tmpv[0*Y+j]-mpiwork_recv1[j])/dx;
          double J12p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxyphi1[j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j+1]*phi[0*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j-1]*phi[0*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j+1)]*phi[0*Y+(j+1)])/2*(tmpv[0*Y+(j+1)]-tmpv[0*Y+j])/dx;
          double J22m=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j-1)]*phi[0*Y+(j-1)])/2*(tmpv[0*Y+j]-tmpv[0*Y+(j-1)])/dx;

          cell[0*Y+j].v=tmpv[0*Y+j]+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[0*Y+j]);
        }
      }
    }
    else {
      #pragma ivdep
      #pragma vector always
      for (int j=0;j<Y;j++) {
        mpiwork_send1[j]=tmpv[0*Y+j];
        mpiwork_send2[j]=tmpv[(X-1)*Y+j];
      }
      MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
      MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
      MPI_Status status;
      MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
      MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
      cell[0*Y+0].v=tmpv[0*Y+0]+(mpiwork_recv1[0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtdx2;
      cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+(mpiwork_recv1[Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtdx2;
      cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+mpiwork_recv2[0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtdx2;
      cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtdx2;
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        if (phi[0*Y+j]>cutoff) {
          double ldv=(mpiwork_recv1[j]+mpiwork_recv1[j-1]+tmpv[0*Y+(j-1)]+tmpv[0*Y+j])/4;
          double luv=(mpiwork_recv1[j]+mpiwork_recv1[j+1]+tmpv[0*Y+(j+1)]+tmpv[0*Y+j])/4;
          double rdv=(tmpv[(0+1)*Y+j]+tmpv[(0+1)*Y+(j-1)]+tmpv[0*Y+(j-1)]+tmpv[0*Y+j])/4;
          double ruv=(tmpv[(0+1)*Y+j]+tmpv[(0+1)*Y+(j+1)]+tmpv[0*Y+(j+1)]+tmpv[0*Y+j])/4;

          double J11p=(Dxx[0*Y+j]*phi[0*Y+j]+Dxx[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(tmpv[(0+1)*Y+j]-tmpv[0*Y+j])/dx;
          double J11m=(Dxx[0*Y+j]*phi[0*Y+j]+Dxxphi1[j])/2*(tmpv[0*Y+j]-mpiwork_recv1[j])/dx;
          double J12p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[(0+1)*Y+j]*phi[(0+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxyphi1[j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j+1]*phi[0*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[0*Y+j]*phi[0*Y+j]+Dxy[0*Y+j-1]*phi[0*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j+1)]*phi[0*Y+(j+1)])/2*(tmpv[0*Y+(j+1)]-tmpv[0*Y+j])/dx;
          double J22m=(Dyy[0*Y+j]*phi[0*Y+j]+Dyy[0*Y+(j-1)]*phi[0*Y+(j-1)])/2*(tmpv[0*Y+j]-tmpv[0*Y+(j-1)])/dx;

          cell[0*Y+j].v=tmpv[0*Y+j]+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[0*Y+j]);
        }
        if (phi[(X-1)*Y+j]>cutoff) {
          double ldv=(tmpv[(X-1-1)*Y+j]+tmpv[(X-1-1)*Y+(j-1)]+tmpv[(X-1)*Y+(j-1)]+tmpv[(X-1)*Y+j])/4;
          double luv=(tmpv[(X-1-1)*Y+j]+tmpv[(X-1-1)*Y+(j+1)]+tmpv[(X-1)*Y+(j+1)]+tmpv[(X-1)*Y+j])/4;
          double rdv=(mpiwork_recv2[j]+mpiwork_recv2[j-1]+tmpv[(X-1)*Y+(j-1)]+tmpv[(X-1)*Y+j])/4;
          double ruv=(mpiwork_recv2[j]+mpiwork_recv2[j+1]+tmpv[(X-1)*Y+(j+1)]+tmpv[(X-1)*Y+j])/4;

          double J11p=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxxphi2[j])/2*(mpiwork_recv2[j]-tmpv[(X-1)*Y+j])/dx;
          double J11m=(Dxx[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxx[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(tmpv[(X-1)*Y+j]-tmpv[(X-1-1)*Y+j])/dx;
          double J12p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxyphi2[j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1-1)*Y+j]*phi[(X-1-1)*Y+j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j+1]*phi[(X-1)*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dxy[(X-1)*Y+j-1]*phi[(X-1)*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j+1)]*phi[(X-1)*Y+(j+1)])/2*(tmpv[(X-1)*Y+(j+1)]-tmpv[(X-1)*Y+j])/dx;
          double J22m=(Dyy[(X-1)*Y+j]*phi[(X-1)*Y+j]+Dyy[(X-1)*Y+(j-1)]*phi[(X-1)*Y+(j-1)])/2*(tmpv[(X-1)*Y+j]-tmpv[(X-1)*Y+(j-1)])/dx;

          cell[(X-1)*Y+j].v=tmpv[(X-1)*Y+j]+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[(X-1)*Y+j]);
        }
      }
    }
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<X-1;i++) {
      cell[i*Y+0].v=tmpv[i*Y+0]+(tmpv[(i-1)*Y+0]+tmpv[(i+1)*Y+0]+tmpv[i*Y+1]+tmpv[i*Y+1]-4*tmpv[i*Y+0])*Dfudtdx2;
      cell[i*Y+Y-1].v=tmpv[i*Y+Y-1]+(tmpv[(i-1)*Y+Y-1]+tmpv[(i+1)*Y+Y-1]+tmpv[i*Y+Y-2]+tmpv[i*Y+Y-2]-4*tmpv[i*Y+Y-1])*Dfudtdx2;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=1;i<X-1;i++) {
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        if (phi[i*Y+j]>cutoff) {
          double ldv=(tmpv[(i-1)*Y+j]+tmpv[(i-1)*Y+(j-1)]+tmpv[i*Y+(j-1)]+tmpv[i*Y+j])/4;
          double luv=(tmpv[(i-1)*Y+j]+tmpv[(i-1)*Y+(j+1)]+tmpv[i*Y+(j+1)]+tmpv[i*Y+j])/4;
          double rdv=(tmpv[(i+1)*Y+j]+tmpv[(i+1)*Y+(j-1)]+tmpv[i*Y+(j-1)]+tmpv[i*Y+j])/4;
          double ruv=(tmpv[(i+1)*Y+j]+tmpv[(i+1)*Y+(j+1)]+tmpv[i*Y+(j+1)]+tmpv[i*Y+j])/4;

          double J11p=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(tmpv[(i+1)*Y+j]-tmpv[i*Y+j])/dx;
          double J11m=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(tmpv[i*Y+j]-tmpv[(i-1)*Y+j])/dx;
          double J12p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j+1]*phi[i*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j-1]*phi[i*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j+1)]*phi[i*Y+(j+1)])/2*(tmpv[i*Y+(j+1)]-tmpv[i*Y+j])/dx;
          double J22m=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j-1)]*phi[i*Y+(j-1)])/2*(tmpv[i*Y+j]-tmpv[i*Y+(j-1)])/dx;

          cell[i*Y+j].v=tmpv[i*Y+j]+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[i*Y+j]);
        }
      }
    }
  }
#else // non-MPI code (diffusion with fiber)
  double Dfudtdx2=Dfu*dt/(dx*dx*DifLoop);
  double dtdx2=dt/(dx*dx*DifLoop);
  for (int lp=0;lp<DifLoop/2;lp++) {
    for (int i=1;i<(X-1);i++) {
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        tmpv[i*Y+j]=cell[i*Y+j].v;
      }
    }
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
    tmpv[0*Y+0]=cell[0*Y+0].v+
                ((Dxx[0*Y+0]*phi[0*Y+0]+Dxx[(0+1)*Y+0]*phi[(0+1)*Y+0])*
                (cell[(0+1)*Y+0].v-cell[0*Y+0].v)+
                (Dyy[0*Y+0]*phi[0*Y+0]+Dyy[0*Y+1]*phi[0*Y+1])*
                (cell[0*Y+1].v-cell[0*Y+0].v))*dtdx2/phi[0*Y+0];
                
    tmpv[(X-1)*Y+0]=cell[(X-1)*Y+0].v+
                ((Dxx[(X-1)*Y+0]*phi[(X-1)*Y+0]+Dxx[(X-2)*Y+0]*phi[(X-2)*Y+0])*
                (cell[(X-2)*Y+0].v-cell[(X-1)*Y+0].v)+
                (Dyy[(X-1)*Y+0]*phi[(X-1)*Y+0]+Dyy[(X-1)*Y+1]*phi[(X-1)*Y+1])*
                (cell[(X-1)*Y+1].v-cell[(X-1)*Y+0].v))*dtdx2/phi[(X-1)*Y+0];
                
    tmpv[0*Y+Y-1]=cell[0*Y+Y-1].v+
                ((Dxx[0*Y+Y-1]*phi[0*Y+Y-1]+Dxx[1*Y+Y-1]*phi[1*Y+Y-1])*
                (cell[1*Y+Y-1].v-cell[0*Y+Y-1].v)+
                (Dyy[0*Y+Y-1]*phi[0*Y+Y-1]+Dyy[0*Y+Y-2]*phi[0*Y+Y-2])*
                (cell[0*Y+Y-2].v-cell[0*Y+Y-1].v))*dtdx2/phi[0*Y+Y-1];
                
    tmpv[(X-1)*Y+Y-1]=cell[(X-1)*Y+Y-1].v+
                ((Dxx[(X-1)*Y+Y-1]*phi[(X-1)*Y+Y-1]+Dxx[(X-2)*Y+Y-1]*phi[(X-2)*Y+Y-1])*
                (cell[(X-2)*Y+Y-1].v-cell[(X-1)*Y+Y-1].v)+
                (Dyy[(X-1)*Y+Y-1]*phi[(X-1)*Y+Y-1]+Dyy[(X-1)*Y+Y-2]*phi[(X-1)*Y+Y-2])*
                (cell[(X-1)*Y+Y-2].v-cell[(X-1)*Y+Y-1].v))*dtdx2/phi[(X-1)*Y+Y-1];
                
                
                
    #ifdef _OPENMP
    #pragma omp for
    #endif
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<(X-1);i++) {
      tmpv[i*Y+0]=cell[i*Y+0].v+
                  ((Dxx[i*Y+0]*phi[i*Y+0]+Dxx[(i-1)*Y+0]*phi[(i-1)*Y+0])/2*
                  (cell[(i-1)*Y+0].v-cell[i*Y+0].v)+
                  (Dxx[i*Y+0]*phi[i*Y+0]+Dxx[(i+1)*Y+0]*phi[(i+1)*Y+0])/2*
                  (cell[(i+1)*Y+0].v-cell[i*Y+0].v)+
                  (Dyy[i*Y+0]*phi[i*Y+0]+Dyy[i*Y+1]*phi[i*Y+1])*
                  (cell[i*Y+1].v-cell[i*Y+0].v))*dtdx2/phi[i*Y+0];

      tmpv[i*Y+Y-1]=cell[i*Y+Y-1].v+
                    ((Dxx[i*Y+Y-1]*phi[i*Y+Y-1]+Dxx[(i-1)*Y+Y-1]*phi[(i-1)*Y+Y-1])/2*
                    (cell[(i-1)*Y+Y-1].v-cell[i*Y+Y-1].v)+
                    (Dxx[i*Y+Y-1]*phi[i*Y+Y-1]+Dxx[(i+1)*Y+Y-1]*phi[(i+1)*Y+Y-1])/2*
                    (cell[(i+1)*Y+Y-1].v-cell[i*Y+Y-1].v)+
                    (Dyy[i*Y+Y-1]*phi[i*Y+Y-1]+Dyy[i*Y+Y-2]*phi[i*Y+Y-2])*
                    (cell[i*Y+Y-2].v-cell[i*Y+Y-1].v))*dtdx2/phi[i*Y+Y-1];
                    

    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<Y-1;i++) {
      tmpv[0*Y+i]=cell[0*Y+i].v+
                  ((Dyy[0*Y+i]*phi[0*Y+i]+Dyy[0*Y+(i-1)]*phi[0*Y+(i-1)])/2*
                  (cell[0*Y+(i-1)].v-cell[0*Y+i].v)+
                  (Dyy[0*Y+i]*phi[0*Y+i]+Dyy[0*Y+(i+1)]*phi[0*Y+(i+1)])/2*
                  (cell[0*Y+(i+1)].v-cell[0*Y+i].v)+
                  (Dxx[0*Y+i]*phi[0*Y+i]+Dxx[1*Y+i]*phi[1*Y+i])*
                  (cell[1*Y+i].v-cell[0*Y+i].v))*dtdx2/phi[0*Y+i];

      tmpv[(X-1)*Y+i]=cell[(X-1)*Y+i].v+
                  ((Dyy[(X-1)*Y+i]*phi[(X-1)*Y+i]+Dyy[(X-1)*Y+(i-1)]*phi[(X-1)*Y+(i-1)])/2*
                  (cell[(X-1)*Y+(i-1)].v-cell[(X-1)*Y+i].v)+
                  (Dyy[(X-1)*Y+i]*phi[(X-1)*Y+i]+Dyy[(X-1)*Y+(i+1)]*phi[(X-1)*Y+(i+1)])/2*
                  (cell[(X-1)*Y+(i+1)].v-cell[(X-1)*Y+i].v)+
                  (Dxx[(X-1)*Y+i]*phi[(X-1)*Y+i]+Dxx[(X-2)*Y+i]*phi[(X-2)*Y+i])*
                  (cell[(X-2)*Y+i].v-cell[(X-1)*Y+i].v))*dtdx2/phi[(X-1)*Y+i];

    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (int i=1;i<(X-1);i++) {
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        if (phi[i*Y+j]>cutoff) {
          double ldv=(cell[(i-1)*Y+j].v+cell[(i-1)*Y+(j-1)].v+cell[i*Y+(j-1)].v+cell[i*Y+j].v)/4;
          double luv=(cell[(i-1)*Y+j].v+cell[(i-1)*Y+(j+1)].v+cell[i*Y+(j+1)].v+cell[i*Y+j].v)/4;
          double rdv=(cell[(i+1)*Y+j].v+cell[(i+1)*Y+(j-1)].v+cell[i*Y+(j-1)].v+cell[i*Y+j].v)/4;
          double ruv=(cell[(i+1)*Y+j].v+cell[(i+1)*Y+(j+1)].v+cell[i*Y+(j+1)].v+cell[i*Y+j].v)/4;

          double J11p=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(cell[(i+1)*Y+j].v-cell[i*Y+j].v)/dx;
          double J11m=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(cell[i*Y+j].v-cell[(i-1)*Y+j].v)/dx;
          double J12p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j+1]*phi[i*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j-1]*phi[i*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j+1)]*phi[i*Y+(j+1)])/2*(cell[i*Y+(j+1)].v-cell[i*Y+j].v)/dx;
          double J22m=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j-1)]*phi[i*Y+(j-1)])/2*(cell[i*Y+j].v-cell[i*Y+(j-1)].v)/dx;

          tmpv[i*Y+j]=cell[i*Y+j].v+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[i*Y+j]);
        }
      }
    }

    cell[0*Y+0].v=tmpv[0*Y+0]+
                ((Dxx[0*Y+0]*phi[0*Y+0]+Dxx[(0+1)*Y+0]*phi[(0+1)*Y+0])*
                (tmpv[(0+1)*Y+0]-tmpv[0*Y+0])+
                (Dyy[0*Y+0]*phi[0*Y+0]+Dyy[0*Y+1]*phi[0*Y+1])*
                (tmpv[0*Y+1]-tmpv[0*Y+0]))*dtdx2/phi[0*Y+0];
                
    cell[(X-1)*Y+0].v=tmpv[(X-1)*Y+0]+
                ((Dxx[(X-1)*Y+0]*phi[(X-1)*Y+0]+Dxx[(X-2)*Y+0]*phi[(X-2)*Y+0])*
                (tmpv[(X-2)*Y+0]-tmpv[(X-1)*Y+0])+
                (Dyy[(X-1)*Y+0]*phi[(X-1)*Y+0]+Dyy[(X-1)*Y+1]*phi[(X-1)*Y+1])*
                (tmpv[(X-1)*Y+1]-tmpv[(X-1)*Y+0]))*dtdx2/phi[(X-1)*Y+0];
                
    cell[0*Y+Y-1].v=tmpv[0*Y+Y-1]+
                ((Dxx[0*Y+Y-1]*phi[0*Y+Y-1]+Dxx[1*Y+Y-1]*phi[1*Y+Y-1])*
                (tmpv[1*Y+Y-1]-tmpv[0*Y+Y-1])+
                (Dyy[0*Y+Y-1]*phi[0*Y+Y-1]+Dyy[0*Y+Y-2]*phi[0*Y+Y-2])*
                (tmpv[0*Y+Y-2]-tmpv[0*Y+Y-1]))*dtdx2/phi[0*Y+Y-1];
                
    cell[(X-1)*Y+Y-1].v=tmpv[(X-1)*Y+Y-1]+
                ((Dxx[(X-1)*Y+Y-1]*phi[(X-1)*Y+Y-1]+Dxx[(X-2)*Y+Y-1]*phi[(X-2)*Y+Y-1])*
                (tmpv[(X-2)*Y+Y-1]-tmpv[(X-1)*Y+Y-1])+
                (Dyy[(X-1)*Y+Y-1]*phi[(X-1)*Y+Y-1]+Dyy[(X-1)*Y+Y-2]*phi[(X-1)*Y+Y-2])*
                (tmpv[(X-1)*Y+Y-2]-tmpv[(X-1)*Y+Y-1]))*dtdx2/phi[(X-1)*Y+Y-1];
                
                
                
    #ifdef _OPENMP
    #pragma omp for
    #endif
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<(X-1);i++) {
      
      cell[i*Y+0].v=tmpv[i*Y+0]+
                  ((Dxx[i*Y+0]*phi[i*Y+0]+Dxx[(i-1)*Y+0]*phi[(i-1)*Y+0])/2*
                  (tmpv[(i-1)*Y+0]-tmpv[i*Y+0])+
                  (Dxx[i*Y+0]*phi[i*Y+0]+Dxx[(i+1)*Y+0]*phi[(i+1)*Y+0])/2*
                  (tmpv[(i+1)*Y+0]-tmpv[i*Y+0])+
                  (Dyy[i*Y+0]*phi[i*Y+0]+Dyy[i*Y+1]*phi[i*Y+1])*
                  (tmpv[i*Y+1]-tmpv[i*Y+0]))*dtdx2/phi[i*Y+0];

      cell[i*Y+Y-1].v=tmpv[i*Y+Y-1]+
                    ((Dxx[i*Y+Y-1]*phi[i*Y+Y-1]+Dxx[(i-1)*Y+Y-1]*phi[(i-1)*Y+Y-1])/2*
                    (tmpv[(i-1)*Y+Y-1]-tmpv[i*Y+Y-1])+
                    (Dxx[i*Y+Y-1]*phi[i*Y+Y-1]+Dxx[(i+1)*Y+Y-1]*phi[(i+1)*Y+Y-1])/2*
                    (tmpv[(i+1)*Y+Y-1]-tmpv[i*Y+Y-1])+
                    (Dyy[i*Y+Y-1]*phi[i*Y+Y-1]+Dyy[i*Y+Y-2]*phi[i*Y+Y-2])*
                    (tmpv[i*Y+Y-2]-tmpv[i*Y+Y-1]))*dtdx2/phi[i*Y+Y-1];

    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<Y-1;i++) {
      
      cell[0*Y+i].v=tmpv[0*Y+i]+
                  ((Dyy[0*Y+i]*phi[0*Y+i]+Dyy[0*Y+(i-1)]*phi[0*Y+(i-1)])/2*
                  (tmpv[0*Y+(i-1)]-tmpv[0*Y+i])+
                  (Dyy[0*Y+i]*phi[0*Y+i]+Dyy[0*Y+(i+1)]*phi[0*Y+(i+1)])/2*
                  (tmpv[0*Y+(i+1)]-tmpv[0*Y+i])+
                  (Dxx[0*Y+i]*phi[0*Y+i]+Dxx[1*Y+i]*phi[1*Y+i])*
                  (tmpv[1*Y+i]-tmpv[0*Y+i]))*dtdx2/phi[0*Y+i];

      cell[(X-1)*Y+i].v=tmpv[(X-1)*Y+i]+
                  ((Dyy[(X-1)*Y+i]*phi[(X-1)*Y+i]+Dyy[(X-1)*Y+(i-1)]*phi[(X-1)*Y+(i-1)])/2*
                  (tmpv[(X-1)*Y+(i-1)]-tmpv[(X-1)*Y+i])+
                  (Dyy[(X-1)*Y+i]*phi[(X-1)*Y+i]+Dyy[(X-1)*Y+(i+1)]*phi[(X-1)*Y+(i+1)])/2*
                  (tmpv[(X-1)*Y+(i+1)]-tmpv[(X-1)*Y+i])+
                  (Dxx[(X-1)*Y+i]*phi[(X-1)*Y+i]+Dxx[(X-2)*Y+i]*phi[(X-2)*Y+i])*
                  (tmpv[(X-2)*Y+i]-tmpv[(X-1)*Y+i]))*dtdx2/phi[(X-1)*Y+i];

    }
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (int i=1;i<X-1;i++) {
      #pragma ivdep
      #pragma vector always
      for (int j=1;j<Y-1;j++) {
        if (phi[i*Y+j]>cutoff) {
          double ldv=(tmpv[(i-1)*Y+j]+tmpv[(i-1)*Y+(j-1)]+tmpv[i*Y+(j-1)]+tmpv[i*Y+j])/4;
          double luv=(tmpv[(i-1)*Y+j]+tmpv[(i-1)*Y+(j+1)]+tmpv[i*Y+(j+1)]+tmpv[i*Y+j])/4;
          double rdv=(tmpv[(i+1)*Y+j]+tmpv[(i+1)*Y+(j-1)]+tmpv[i*Y+(j-1)]+tmpv[i*Y+j])/4;
          double ruv=(tmpv[(i+1)*Y+j]+tmpv[(i+1)*Y+(j+1)]+tmpv[i*Y+(j+1)]+tmpv[i*Y+j])/4;

          double J11p=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(tmpv[(i+1)*Y+j]-tmpv[i*Y+j])/dx;
          double J11m=(Dxx[i*Y+j]*phi[i*Y+j]+Dxx[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(tmpv[i*Y+j]-tmpv[(i-1)*Y+j])/dx;
          double J12p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i+1)*Y+j]*phi[(i+1)*Y+j])/2*(ruv-rdv)/dx;
          double J12m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[(i-1)*Y+j]*phi[(i-1)*Y+j])/2*(luv-ldv)/dx;
          double J21p=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j+1]*phi[i*Y+j+1])/2*(ruv-luv)/dx;
          double J21m=(Dxy[i*Y+j]*phi[i*Y+j]+Dxy[i*Y+j-1]*phi[i*Y+j-1])/2*(rdv-ldv)/dx;
          double J22p=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j+1)]*phi[i*Y+(j+1)])/2*(tmpv[i*Y+(j+1)]-tmpv[i*Y+j])/dx;
          double J22m=(Dyy[i*Y+j]*phi[i*Y+j]+Dyy[i*Y+(j-1)]*phi[i*Y+(j-1)])/2*(tmpv[i*Y+j]-tmpv[i*Y+(j-1)])/dx;

          cell[i*Y+j].v=tmpv[i*Y+j]+(J11p-J11m+J12p-J12m+J21p-J21m+J22p-J22m)*dt/(dx*DifLoop*phi[i*Y+j]);
        }
      }
    }
    }
  }
#endif
}

void Tissue::solvephi(void)
{
  double dtx=0.1*2*dx*dx/(xi*xi);
#ifdef ___USE_MPI
  for (double t=0;t<10;t+=dtx) {
    double Dfudtxdx2=xi*xi*dtx/(dx*dx)/2;
    for (int lp=0;lp<DifLoop/2;lp++) {
      if (my_proc==0) {
        for (int j=0;j<Y;j++) {
          mpiwork_send2[j]=phi[(X-1)*Y+j];
        }
        MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);

        tmpv[0*Y+0]=phi[0*Y+0]+(phi[1*Y+0]+phi[1*Y+0]+phi[0*Y+1]+phi[0*Y+1]-4*phi[0*Y+0])*Dfudtxdx2-(8*tmpv[0*Y+0]*tmpv[0*Y+0]*tmpv[0*Y+0]-12*tmpv[0*Y+0]*tmpv[0*Y+0]+4*tmpv[0*Y+0])*dtx/2;
        tmpv[0*Y+Y-1]=phi[0*Y+Y-1]+(phi[1*Y+Y-1]+phi[1*Y+Y-1]+phi[0*Y+Y-2]+phi[0*Y+Y-2]-4*phi[0*Y+Y-1])*Dfudtxdx2-(8*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]-12*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]+4*tmpv[0*Y+Y-1])*dtx/2;
        tmpv[(X-1)*Y+0]=phi[(X-1)*Y+0]+(phi[(X-2)*Y+0]+mpiwork_recv2[0]+phi[(X-1)*Y+1]+phi[(X-1)*Y+1]-4*phi[(X-1)*Y+0])*Dfudtxdx2-(8*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]-12*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]+4*tmpv[(X-1)*Y+0])*dtx/2;
        tmpv[(X-1)*Y+Y-1]=phi[(X-1)*Y+Y-1]+(phi[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+phi[(X-1)*Y+Y-2]+phi[(X-1)*Y+Y-2]-4*phi[(X-1)*Y+Y-1])*Dfudtxdx2-(8*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]-12*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]+4*tmpv[(X-1)*Y+Y-1])*dtx/2;

        for (int i=1;i<Y-1;i++) {
          tmpv[0*Y+i]=phi[0*Y+i]+(phi[0*Y+i-1]+phi[0*Y+i+1]+phi[1*Y+i]+phi[1*Y+i]-4*phi[0*Y+i])*Dfudtxdx2-(8*tmpv[0*Y+i]*tmpv[0*Y+i]*tmpv[0*Y+i]-12*tmpv[0*Y+i]*tmpv[0*Y+i]+4*tmpv[0*Y+i])*dtx/2;
          tmpv[(X-1)*Y+i]=phi[(X-1)*Y+i]+(phi[(X-1)*Y+i-1]+phi[(X-1)*Y+i+1]+phi[(X-2)*Y+i]+mpiwork_recv2[i]-4*phi[(X-1)*Y+i])*Dfudtxdx2;-(8*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]-12*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]+4*tmpv[(X-1)*Y+i])*dtx/2;
        }
      }
      else if (my_proc==num_procs-1) {
        for (int j=0;j<Y;j++) {
          mpiwork_send1[j]=phi[0*Y+j];
        }
        MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);

        tmpv[0*Y+0]=phi[0*Y+0]+(mpiwork_recv1[0]+phi[1*Y+0]+phi[0*Y+1]+phi[0*Y+1]-4*phi[0*Y+0])*Dfudtxdx2-(8*tmpv[0*Y+0]*tmpv[0*Y+0]*tmpv[0*Y+0]-12*tmpv[0*Y+0]*tmpv[0*Y+0]+4*tmpv[0*Y+0])*dtx/2;
        tmpv[0*Y+Y-1]=phi[0*Y+Y-1]+(mpiwork_recv1[Y-1]+phi[1*Y+Y-1]+phi[0*Y+Y-2]+phi[0*Y+Y-2]-4*phi[0*Y+Y-1])*Dfudtxdx2-(8*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]-12*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]+4*tmpv[0*Y+Y-1])*dtx/2;
        tmpv[(X-1)*Y+0]=phi[(X-1)*Y+0]+(phi[(X-2)*Y+0]+phi[(X-2)*Y+0]+phi[(X-1)*Y+1]+phi[(X-1)*Y+1]-4*phi[(X-1)*Y+0])*Dfudtxdx2-(8*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]-12*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]+4*tmpv[(X-1)*Y+0])*dtx/2;
        tmpv[(X-1)*Y+Y-1]=phi[(X-1)*Y+Y-1]+(phi[(X-2)*Y+Y-1]+phi[(X-2)*Y+Y-1]+phi[(X-1)*Y+Y-2]+phi[(X-1)*Y+Y-2]-4*phi[(X-1)*Y+Y-1])*Dfudtxdx2-(8*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]-12*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]+4*tmpv[(X-1)*Y+Y-1])*dtx/2;

        for (int i=1;i<Y-1;i++)
        {
          tmpv[0*Y+i]=phi[0*Y+i]+(phi[0*Y+i-1]+phi[0*Y+i+1]+mpiwork_recv1[i]+phi[1*Y+i]-4*phi[0*Y+i])*Dfudtxdx2-(8*tmpv[0*Y+i]*tmpv[0*Y+i]*tmpv[0*Y+i]-12*tmpv[0*Y+i]*tmpv[0*Y+i]+4*tmpv[0*Y+i])*dtx/2;
          tmpv[(X-1)*Y+i]=phi[(X-1)*Y+i]+(phi[(X-1)*Y+i-1]+phi[(X-1)*Y+i+1]+phi[(X-2)*Y+i]+phi[(X-2)*Y+i]-4*phi[(X-1)*Y+i])*Dfudtxdx2-(8*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]-12*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]+4*tmpv[(X-1)*Y+i])*dtx/2;
        }
      }
      else {
        for (int j=0;j<Y;j++) {
          mpiwork_send1[j]=phi[0*Y+j];
          mpiwork_send2[j]=phi[(X-1)*Y+j];
        }
        MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
        MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
        MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
        tmpv[0*Y+0]=phi[0*Y+0]+(mpiwork_recv1[0]+phi[1*Y+0]+phi[0*Y+1]+phi[0*Y+1]-4*phi[0*Y+0])*Dfudtxdx2-(8*tmpv[0*Y+0]*tmpv[0*Y+0]*tmpv[0*Y+0]-12*tmpv[0*Y+0]*tmpv[0*Y+0]+4*tmpv[0*Y+0])*dtx/2;
        tmpv[0*Y+Y-1]=phi[0*Y+Y-1]+(mpiwork_recv1[Y-1]+phi[1*Y+Y-1]+phi[0*Y+Y-2]+phi[0*Y+Y-2]-4*phi[0*Y+Y-1])*Dfudtxdx2-(8*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]-12*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]+4*tmpv[0*Y+Y-1])*dtx/2;
        tmpv[(X-1)*Y+0]=phi[(X-1)*Y+0]+(phi[(X-2)*Y+0]+mpiwork_recv2[0]+phi[(X-1)*Y+1]+phi[(X-1)*Y+1]-4*phi[(X-1)*Y+0])*Dfudtxdx2-(8*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]-12*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]+4*tmpv[(X-1)*Y+0])*dtx/2;
        tmpv[(X-1)*Y+Y-1]=phi[(X-1)*Y+Y-1]+(phi[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+phi[(X-1)*Y+Y-2]+phi[(X-1)*Y+Y-2]-4*phi[(X-1)*Y+Y-1])*Dfudtxdx2-(8*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]-12*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]+4*tmpv[(X-1)*Y+Y-1])*dtx/2;
        for (int i=1;i<Y-1;i++)
        {
          tmpv[0*Y+i]=phi[0*Y+i]+(phi[0*Y+i-1]+phi[0*Y+i+1]+mpiwork_recv1[i]+phi[1*Y+i]-4*phi[0*Y+i])*Dfudtxdx2-(8*tmpv[0*Y+i]*tmpv[0*Y+i]*tmpv[0*Y+i]-12*tmpv[0*Y+i]*tmpv[0*Y+i]+4*tmpv[0*Y+i])*dtx/2;
          tmpv[(X-1)*Y+i]=phi[(X-1)*Y+i]+(phi[(X-1)*Y+i-1]+phi[(X-1)*Y+i+1]+phi[(X-2)*Y+i]+mpiwork_recv2[i]-4*phi[(X-1)*Y+i])*Dfudtxdx2-(8*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]-12*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]+4*tmpv[(X-1)*Y+i])*dtx/2;
        }
      }
      for (int i=1;i<X-1;i++) {
        tmpv[i*Y+0]=phi[i*Y+0]+(phi[(i-1)*Y+0]+phi[(i+1)*Y+0]+phi[i*Y+1]+phi[i*Y+1]-4*phi[i*Y+0])*Dfudtxdx2-(8*tmpv[i*Y+0]*tmpv[i*Y+0]*tmpv[i*Y+0]-12*tmpv[i*Y+0]*tmpv[i*Y+0]+4*tmpv[i*Y+0])*dtx/2;
        tmpv[i*Y+Y-1]=phi[i*Y+Y-1]+(phi[(i-1)*Y+Y-1]+phi[(i+1)*Y+Y-1]+phi[i*Y+Y-2]+phi[i*Y+Y-2]-4*phi[i*Y+Y-1])*Dfudtxdx2-(8*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]-12*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]+4*tmpv[i*Y+Y-1])*dtx/2;
      }
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int i=1;i<X-1;i++)
        for (int j=1;j<Y-1;j++)
          tmpv[i*Y+j]=phi[i*Y+j]+(phi[(i-1)*Y+j]+phi[(i+1)*Y+j]+phi[i*Y+j+1]+phi[i*Y+j-1]-4*phi[i*Y+j])*Dfudtxdx2-(8*tmpv[i*Y+j]*tmpv[i*Y+j]*tmpv[i*Y+j]-12*tmpv[i*Y+j]*tmpv[i*Y+j]+4*tmpv[i*Y+j])*dtx/2;

      if (my_proc==0) {
        for (int j=0;j<Y;j++) {
          mpiwork_send2[j]=tmpv[(X-1)*Y+j];
        }
        MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);

        phi[0*Y+0]=tmpv[0*Y+0]+(tmpv[1*Y+0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtxdx2-(8*tmpv[0*Y+0]*tmpv[0*Y+0]*tmpv[0*Y+0]-12*tmpv[0*Y+0]*tmpv[0*Y+0]+4*tmpv[0*Y+0])*dtx/2;
        phi[0*Y+Y-1]=tmpv[0*Y+Y-1]+(tmpv[1*Y+Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtxdx2-(8*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]-12*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]+4*tmpv[0*Y+Y-1])*dtx/2;
        phi[(X-1)*Y+0]=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+mpiwork_recv2[0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtxdx2-(8*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]-12*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]+4*tmpv[(X-1)*Y+0])*dtx/2;
        phi[(X-1)*Y+Y-1]=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtxdx2-(8*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]-12*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]+4*tmpv[(X-1)*Y+Y-1])*dtx/2;
        for (int i=1;i<Y-1;i++) {
          phi[0*Y+i]=tmpv[0*Y+i]+(tmpv[0*Y+i-1]+tmpv[0*Y+i+1]+tmpv[1*Y+i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtxdx2-(8*tmpv[0*Y+i]*tmpv[0*Y+i]*tmpv[0*Y+i]-12*tmpv[0*Y+i]*tmpv[0*Y+i]+4*tmpv[0*Y+i])*dtx/2;
          phi[(X-1)*Y+i]=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+i-1]+tmpv[(X-1)*Y+i+1]+tmpv[(X-2)*Y+i]+mpiwork_recv2[i]-4*tmpv[(X-1)*Y+i])*Dfudtxdx2-(8*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]-12*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]+4*tmpv[(X-1)*Y+i])*dtx/2;
        }
      }
      else if (my_proc==num_procs-1) {
        for (int j=0;j<Y;j++) {
          mpiwork_send1[j]=tmpv[0*Y+j];
        }
        MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
        phi[0*Y+0]=tmpv[0*Y+0]+(mpiwork_recv1[0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtxdx2-(8*tmpv[0*Y+0]*tmpv[0*Y+0]*tmpv[0*Y+0]-12*tmpv[0*Y+0]*tmpv[0*Y+0]+4*tmpv[0*Y+0])*dtx/2;
        phi[0*Y+Y-1]=tmpv[0*Y+Y-1]+(mpiwork_recv1[Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtxdx2-(8*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]-12*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]+4*tmpv[0*Y+Y-1])*dtx/2;
        phi[(X-1)*Y+0]=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+tmpv[(X-2)*Y+0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtxdx2-(8*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]-12*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]+4*tmpv[(X-1)*Y+0])*dtx/2;
        phi[(X-1)*Y+Y-1]=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+tmpv[(X-2)*Y+Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtxdx2-(8*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]-12*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]+4*tmpv[(X-1)*Y+Y-1])*dtx/2;
        for (int i=1;i<Y-1;i++) {
          phi[0*Y+i]=tmpv[0*Y+i]+(tmpv[0*Y+i-1]+tmpv[0*Y+i+1]+mpiwork_recv1[i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtxdx2-(8*tmpv[0*Y+i]*tmpv[0*Y+i]*tmpv[0*Y+i]-12*tmpv[0*Y+i]*tmpv[0*Y+i]+4*tmpv[0*Y+i])*dtx/2;
          phi[(X-1)*Y+i]=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+i-1]+tmpv[(X-1)*Y+i+1]+tmpv[(X-2)*Y+i]+tmpv[(X-2)*Y+i]-4*tmpv[(X-1)*Y+i])*Dfudtxdx2-(8*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]-12*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]+4*tmpv[(X-1)*Y+i])*dtx/2;
        }
      }
      else {
        for (int j=0;j<Y;j++) {
          mpiwork_send1[j]=tmpv[0*Y+j];
          mpiwork_send2[j]=tmpv[(X-1)*Y+j];
        }
        MPI_Send(mpiwork_send1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD);
        MPI_Send(mpiwork_send2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(mpiwork_recv1,Y,MPI_DOUBLE,my_proc-1,0,MPI_COMM_WORLD,&status);
        MPI_Recv(mpiwork_recv2,Y,MPI_DOUBLE,my_proc+1,0,MPI_COMM_WORLD,&status);
        phi[0*Y+0]=tmpv[0*Y+0]+(mpiwork_recv1[0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtxdx2-(8*tmpv[0*Y+0]*tmpv[0*Y+0]*tmpv[0*Y+0]-12*tmpv[0*Y+0]*tmpv[0*Y+0]+4*tmpv[0*Y+0])*dtx/2;
        phi[0*Y+Y-1]=tmpv[0*Y+Y-1]+(mpiwork_recv1[Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtxdx2-(8*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]-12*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]+4*tmpv[0*Y+Y-1])*dtx/2;
        phi[(X-1)*Y+0]=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+mpiwork_recv2[0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtxdx2-(8*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]-12*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]+4*tmpv[(X-1)*Y+0])*dtx/2;
        phi[(X-1)*Y+Y-1]=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+mpiwork_recv2[Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtxdx2-(8*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]-12*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]+4*tmpv[(X-1)*Y+Y-1])*dtx/2;
        for (int i=1;i<Y-1;i++) {
          phi[0*Y+i]=tmpv[0*Y+i]+(tmpv[0*Y+i-1]+tmpv[0*Y+i+1]+mpiwork_recv1[i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtxdx2-(8*tmpv[0*Y+i]*tmpv[0*Y+i]*tmpv[0*Y+i]-12*tmpv[0*Y+i]*tmpv[0*Y+i]+4*tmpv[0*Y+i])*dtx/2;
          phi[(X-1)*Y+i]=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+i-1]+tmpv[(X-1)*Y+i+1]+tmpv[(X-2)*Y+i]+mpiwork_recv2[i]-4*tmpv[(X-1)*Y+i])*Dfudtxdx2-(8*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]-12*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]+4*tmpv[(X-1)*Y+i])*dtx/2;
        }
      }
      for (int i=1;i<X-1;i++) {
        phi[i*Y+0]=tmpv[i*Y+0]+(tmpv[(i-1)*Y+0]+tmpv[(i+1)*Y+0]+tmpv[i*Y+1]+tmpv[i*Y+1]-4*tmpv[i*Y+0])*Dfudtxdx2-(8*tmpv[i*Y+0]*tmpv[i*Y+0]*tmpv[i*Y+0]-12*tmpv[i*Y+0]*tmpv[i*Y+0]+4*tmpv[i*Y+0])*dtx/2;
        phi[i*Y+Y-1]=tmpv[i*Y+Y-1]+(tmpv[(i-1)*Y+Y-1]+tmpv[(i+1)*Y+Y-1]+tmpv[i*Y+Y-2]+tmpv[i*Y+Y-2]-4*tmpv[i*Y+Y-1])*Dfudtxdx2-(8*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]-12*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]+4*tmpv[i*Y+Y-1])*dtx/2;
      }
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      #pragma ivdep
      #pragma vector always
      for (int i=1;i<X-1;i++)
        for (int j=1;j<Y-1;j++)
          phi[i*Y+j]=tmpv[i*Y+j]+(tmpv[(i-1)*Y+j]+tmpv[(i+1)*Y+j]+tmpv[i*Y+j+1]+tmpv[i*Y+j-1]-4*tmpv[i*Y+j])*Dfudtxdx2-(8*tmpv[i*Y+j]*tmpv[i*Y+j]*tmpv[i*Y+j]-12*tmpv[i*Y+j]*tmpv[i*Y+j]+4*tmpv[i*Y+j])*dtx/2;
    }
  }

#else

  for (double t=0;t<10;t+=dtx) {
    double Dfudtxdx2=xi*xi*dtx/(dx*dx)/2;
    tmpv[0*Y+0]=phi[0*Y+0]+(phi[1*Y+0]+phi[1*Y+0]+phi[0*Y+1]+phi[0*Y+1]-4*phi[0*Y+0])*Dfudtxdx2-(8*phi[0*Y+0]*phi[0*Y+0]*phi[0*Y+0]-12*phi[0*Y+0]*phi[0*Y+0]+4*phi[0*Y+0])*dtx/2;
    tmpv[(X-1)*Y+0]=phi[(X-1)*Y+0]+(phi[(X-2)*Y+0]+phi[(X-2)*Y+0]+phi[(X-1)*Y+1]+phi[(X-1)*Y+1]-4*phi[(X-1)*Y+0])*Dfudtxdx2-(8*phi[(X-1)*Y+0]*phi[(X-1)*Y+0]*phi[(X-1)*Y+0]-12*phi[(X-1)*Y+0]*phi[(X-1)*Y+0]+4*phi[(X-1)*Y+0])*dtx/2;
    tmpv[0*Y+Y-1]=phi[0*Y+Y-1]+(phi[1*Y+Y-1]+phi[1*Y+Y-1]+phi[0*Y+Y-2]+phi[0*Y+Y-2]-4*phi[0*Y+Y-1])*Dfudtxdx2-(8*phi[0*Y+Y-1]*phi[0*Y+Y-1]*phi[0*Y+Y-1]-12*phi[0*Y+Y-1]*phi[0*Y+Y-1]+4*phi[0*Y+Y-1])*dtx/2;
    tmpv[(X-1)*Y+Y-1]=phi[(X-1)*Y+Y-1]+(phi[(X-2)*Y+Y-1]+phi[(X-2)*Y+Y-1]+phi[(X-1)*Y+Y-2]+phi[(X-1)*Y+Y-2]-4*phi[(X-1)*Y+Y-1])*Dfudtxdx2-(8*phi[(X-1)*Y+Y-1]*phi[(X-1)*Y+Y-1]*phi[(X-1)*Y+Y-1]-12*phi[(X-1)*Y+Y-1]*phi[(X-1)*Y+Y-1]+4*phi[(X-1)*Y+Y-1])*dtx/2;
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<(X-1);i++) {
      tmpv[i*Y+0]=phi[i*Y+0]+(phi[(i-1)*Y+0]+phi[(i+1)*Y+0]+phi[i*Y+1]+phi[i*Y+1]-4*phi[i*Y+0])*Dfudtxdx2-(8*phi[i*Y+0]*phi[i*Y+0]*phi[i*Y+0]-12*phi[i*Y+0]*phi[i*Y+0]+4*phi[i*Y+0])*dtx/2;
      tmpv[i*Y+Y-1]=phi[i*Y+Y-1]+(phi[(i-1)*Y+Y-1]+phi[(i+1)*Y+Y-1]+phi[i*Y+Y-2]+phi[i*Y+Y-2]-4*phi[i*Y+Y-1])*Dfudtxdx2-(8*phi[i*Y+Y-1]*phi[i*Y+Y-1]*phi[i*Y+Y-1]-12*phi[i*Y+Y-1]*phi[i*Y+Y-1]+4*phi[i*Y+Y-1])*dtx/2;
    }
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<Y-1;i++) {
      tmpv[0*Y+i]=phi[0*Y+i]+(phi[0*Y+(i-1)]+phi[0*Y+(i+1)]+phi[1*Y+i]+phi[1*Y+i]-4*phi[0*Y+i])*Dfudtxdx2-(8*phi[0*Y+i]*phi[0*Y+i]*phi[0*Y+i]-12*phi[0*Y+i]*phi[0*Y+i]+4*phi[0*Y+i])*dtx/2;
      tmpv[(X-1)*Y+i]=phi[(X-1)*Y+i]+(phi[(X-1)*Y+(i-1)]+phi[(X-1)*Y+(i+1)]+phi[(X-2)*Y+i]+phi[(X-2)*Y+i]-4*phi[(X-1)*Y+i])*Dfudtxdx2-(8*phi[(X-1)*Y+i]*phi[(X-1)*Y+i]*phi[(X-1)*Y+i]-12*phi[(X-1)*Y+i]*phi[(X-1)*Y+i]+4*phi[(X-1)*Y+i])*dtx/2;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<(X-1);i++)
      for (int j=1;j<Y-1;j++)
        tmpv[i*Y+j]=phi[i*Y+j]+(phi[(i-1)*Y+j]+phi[(i+1)*Y+j]+phi[i*Y+j+1]+phi[i*Y+j-1]-4*phi[i*Y+j])*Dfudtxdx2-(8*phi[i*Y+j]*phi[i*Y+j]*phi[i*Y+j]-12*phi[i*Y+j]*phi[i*Y+j]+4*phi[i*Y+j])*dtx/2;
    phi[0*Y+0]=tmpv[0*Y+0]+(tmpv[1*Y+0]+tmpv[1*Y+0]+tmpv[0*Y+1]+tmpv[0*Y+1]-4*tmpv[0*Y+0])*Dfudtxdx2-(8*tmpv[0*Y+0]*tmpv[0*Y+0]*tmpv[0*Y+0]-12*tmpv[0*Y+0]*tmpv[0*Y+0]+4*tmpv[0*Y+0])*dtx/2;
    phi[(X-1)*Y+0]=tmpv[(X-1)*Y+0]+(tmpv[(X-2)*Y+0]+tmpv[(X-2)*Y+0]+tmpv[(X-1)*Y+1]+tmpv[(X-1)*Y+1]-4*tmpv[(X-1)*Y+0])*Dfudtxdx2-(8*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]-12*tmpv[(X-1)*Y+0]*tmpv[(X-1)*Y+0]+4*tmpv[(X-1)*Y+0])*dtx/2;
    phi[0*Y+Y-1]=tmpv[0*Y+Y-1]+(tmpv[1*Y+Y-1]+tmpv[1*Y+Y-1]+tmpv[0*Y+Y-2]+tmpv[0*Y+Y-2]-4*tmpv[0*Y+Y-1])*Dfudtxdx2-(8*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]-12*tmpv[0*Y+Y-1]*tmpv[0*Y+Y-1]+4*tmpv[0*Y+Y-1])*dtx/2;
    phi[(X-1)*Y+Y-1]=tmpv[(X-1)*Y+Y-1]+(tmpv[(X-2)*Y+Y-1]+tmpv[(X-2)*Y+Y-1]+tmpv[(X-1)*Y+Y-2]+tmpv[(X-1)*Y+Y-2]-4*tmpv[(X-1)*Y+Y-1])*Dfudtxdx2-(8*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]-12*tmpv[(X-1)*Y+Y-1]*tmpv[(X-1)*Y+Y-1]+4*tmpv[(X-1)*Y+Y-1])*dtx/2;
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<(X-1);i++) {
      phi[i*Y+0]=tmpv[i*Y+0]+(tmpv[(i-1)*Y+0]+tmpv[(i+1)*Y+0]+tmpv[i*Y+1]+tmpv[i*Y+1]-4*tmpv[i*Y+0])*Dfudtxdx2-(8*tmpv[i*Y+0]*tmpv[i*Y+0]*tmpv[i*Y+0]-12*tmpv[i*Y+0]*tmpv[i*Y+0]+4*tmpv[i*Y+0])*dtx/2;
      phi[i*Y+Y-1]=tmpv[i*Y+Y-1]+(tmpv[(i-1)*Y+Y-1]+tmpv[(i+1)*Y+Y-1]+tmpv[i*Y+Y-2]+tmpv[i*Y+Y-2]-4*tmpv[i*Y+Y-1])*Dfudtxdx2-(8*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]-12*tmpv[i*Y+Y-1]*tmpv[i*Y+Y-1]+4*tmpv[i*Y+Y-1])*dtx/2;
    }
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<Y-1;i++) {
      phi[0*Y+i]=tmpv[0*Y+i]+(tmpv[0*Y+(i-1)]+tmpv[0*Y+(i+1)]+tmpv[1*Y+i]+tmpv[1*Y+i]-4*tmpv[0*Y+i])*Dfudtxdx2-(8*tmpv[0*Y+i]*tmpv[0*Y+i]*tmpv[0*Y+i]-12*tmpv[0*Y+i]*tmpv[0*Y+i]+4*tmpv[0*Y+i])*dtx/2;
      phi[(X-1)*Y+i]=tmpv[(X-1)*Y+i]+(tmpv[(X-1)*Y+(i-1)]+tmpv[(X-1)*Y+(i+1)]+tmpv[(X-2)*Y+i]+tmpv[(X-2)*Y+i]-4*tmpv[(X-1)*Y+i])*Dfudtxdx2-(8*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]-12*tmpv[(X-1)*Y+i]*tmpv[(X-1)*Y+i]+4*tmpv[(X-1)*Y+i])*dtx/2;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    #pragma ivdep
    #pragma vector always
    for (int i=1;i<X-1;i++)
      for (int j=1;j<Y-1;j++)
        phi[i*Y+j]=tmpv[i*Y+j]+(tmpv[(i-1)*Y+j]+tmpv[(i+1)*Y+j]+tmpv[i*Y+j+1]+tmpv[i*Y+j-1]-4*tmpv[i*Y+j])*Dfudtxdx2-(8*tmpv[i*Y+j]*tmpv[i*Y+j]*tmpv[i*Y+j]-12*tmpv[i*Y+j]*tmpv[i*Y+j]+4*tmpv[i*Y+j])*dtx/2;
  }
#endif

}
void Tissue::writephi(string filename,int step) {
#ifdef ___USE_MPI
  double *tmpphi=new double[X*Y];
  double *tmpphi2=new double[X*Y*num_procs];
  for(int i=0;i<X;i++) {
    for(int j=0;j<Y;j++) {
      tmpphi[i*Y+j]=phi[i*Y+j];
    }
  }
  MPI_Gather(tmpphi,X*Y,MPI_DOUBLE,tmpphi2,X*Y,MPI_DOUBLE,0,MPI_COMM_WORLD);
  if (my_proc==0) {
    ofstream os(filename.c_str());
    for (int i=0;i<GX;i+=step) {
      for (int j=0;j<GY;j+=step) {
        if(fabs(tmpphi2[i*GY+j])<0.000001)tmpphi2[i*GY+j]=0;
        os<<tmpphi2[i*GY+j]<<"\t";
      }
      os<<endl;
    }
  }
  delete[] tmpphi,tmpphi2;
#else
  ofstream os(filename.c_str());
  for (int i=0;i<X;i+=step) {
    for (int j=0;j<Y;j+=step) {
      if(fabs(phi[i*Y+j])<0.000001)
        os<<0<<"\t";
      else
        os<<phi[i*Y+j]<<"\t";
    }
    os<<endl;
  }
#endif
}
void Tissue::setphi(int xx,int yy,double phii) {
#ifdef ___USE_MPI
  int x2=xx-my_proc*X;
  if (x2>=0 && x2<X) {
    phi[x2*Y+yy]=phii;
  }
#else
  phi[xx*Y+yy]=phii;
#endif
}
void Tissue::setdif(double D) {
  Dfu=D;
  if (Y!=1) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i=0;i<X;i++) {
      for(int j=0;j<Y;j++) {
        Dxx[i*Y+j]=Dfu;
        Dxy[i*Y+j]=0;
        Dyy[i*Y+j]=Dfu;
        phi[i*Y+j]=1;
      }
    }
  }
}
void Tissue::setDxx(double theta,double dxx, double dyy, int posx, int posy) {
  if (posx==-1) {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i=0;i<X;i++) {
      #pragma ivdep
      #pragma vector always
      for (int j=0;j<Y;j++) {
        Dxx[i*Y+j]=(dxx-dyy)/2*cos(2*theta)+(dxx+dyy)/2;
        Dxy[i*Y+j]=-(dxx-dyy)/2*sin(2*theta);
        Dyy[i*Y+j]=-(dxx-dyy)/2*cos(2*theta)+(dxx+dyy)/2;
      }
    }
  }
  else {
    Dxx[posx*Y+posy]=(dxx-dyy)/2*cos(2*theta)+(dxx+dyy)/2;
    Dxy[posx*Y+posy]=-(dxx-dyy)/2*sin(2*theta);
    Dyy[posx*Y+posy]=-(dxx-dyy)/2*cos(2*theta)+(dxx+dyy)/2;
  }
}
