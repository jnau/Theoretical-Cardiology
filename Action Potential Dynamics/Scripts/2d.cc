#include <iostream>
#include <fstream>
using namespace std;

#include "cell.h"
#include "tissue.h"
#include "recorder2.h"
#include "log.h"

int main(void) {
  Logging lg;
  //CCell cell;

  //change this to 10
  const int X = 1000;
  const int Y = 1000;
  Tissue tis(X, Y);
  ofstream os1("voltks130.txt"); // voltage exported onto text file.
  ofstream os2("apdks130.txt"); // checking the timing of AP.


  int t1 = 0;
  int t2 = 0;
  int t3 = 0;
  double volt1 = 0;
  //double volt2 = 0;
  //int tm = 0;


  for (int i = 0; i < X; i++) {
    for (int j = 0; j < Y; j++) {
      tis.cell[i*Y + j].gna *= 1.0;//Na conductance
      tis.cell[i*Y + j].gca *= 1.0;//Ca conductance
      tis.cell[i*Y + j].gkr *= 1.0;//K conductance
      tis.cell[i*Y + j].gks *= 1.3;//K conductance
      tis.cell[i*Y + j].gtos *= 1.0;//K conductance
      tis.cell[i*Y + j].gtof *= 1.0;//K conductance
      tis.cell[i*Y + j].gkix *= 1.0;//K conductance
      tis.cell[i*Y + j].gnak *= 1.0;//K conductance K1
      tis.cell[i*Y + j].gnaca *= 1.0;//K conductance

      }

  }

  
  const double dt = 0.1;// time step (0.1 ms)
  REC2PARAM param;
  param.filenamev="vks130.dat";
  param.filenamevmm="vmmks130.txt";
  param.filenamevmovie = "vks130.mp4";
  param.step=3;
//  param.gzip = false;
  Recorder2 rec(&tis, &param);

  int durn = 2.0 / dt;//duration of stimulation (1.0 ms)



  //main loop
  int ttn = 0;

// timing base it on voltage. you want when it is -80.

  
  int count = 0;
  double tmax = 300;
  int tnmax = tmax / dt;//convert to integer
  int apd = 0; // initializing variables
  int apdn = 0;
  int adj = tnmax;


  for (int tn = 0; tn < tnmax; tn++,ttn++) {

    double stim = 0;//stimulation current
    if (tn < durn) {
      stim = 80;
    }

    tis.pace(stim,STIM_2D_PARALLEL2);


    // checking voltage here, after the pace function
    volt1 = tis.cell[tn].v;
    os1 << tn << ' ' << volt1 << endl;



    // if-else statements to find when the voltage is > 10 and less than -80
    // made to find the timing
    if (volt1 > 10 && count == 0){
      count = 1;
      t1 = tn;
    }

    if (volt1 < -80 && count == 1){
      count = 2;
      t2 = tn;
      apd = t2 - t1; //find apd
      

      // the APD is then multiplied by a ratio to keep consistency. (~200-300)
      if (apd < 25){
        apd = apd*5;
        apdn = apd/dt;
      }

      else if (apd >= 25 && apd < 56){
        apd = apd*3.5;
        apdn = apd/dt;
      }

      else if (apd >= 57 && apd < 60){
        apd = apd*3;
        apdn = apd/dt;
      }

      else if (apd >= 60 && apd < 75){
        apd = apd*2.5;
        apdn = apd/dt;
      }

      else if (apd >= 75 && apd < 90){
        apd = apd*2;
        apdn = apd/dt;
      }

      else if (apd >= 90 && apd < 120){
        apd = apd*1.5;
        apdn = apd/dt;
      }

      else if (apd >= 120 && apd < 141){
        apd = apd*1.25;
        apdn = apd/dt;
      } 
      
      else if (apd >= 142 && apd < 195){
        apd = apd*1.1;
        apdn = apd/dt;
      } 

      else if (apd >= 195 && apd < 250){
        apd = apd/1.25;
        apdn = apd/dt;
      } //keep an eye on

      else if (apd >= 250 && apd < 310){

        apd = apd/1.5;
        apdn = apd/dt;
      }
/*
      else if (apd >= 275 && apd < 310){
        apd = apd/1.75;
        apdn = apd/dt;
      }
*/
      else if (apd >= 310){

        apd = apd/2;
        apdn = apd/dt;
      }

      else{

        apdn = apd/dt;
      }

      adj = apdn + t1; // to take into account for the timing in the beginning 
      os2 << apdn << ' ' << t1 << ' ' << t2 << ' '<< adj << endl; 
    }

    if (adj == tn){
      break;
    }

    //if (tn % 30 == 0 && count != 4){//write results every 3 ms
      if (tn % 30 == 0){
      double t = ttn*dt;
      rec.recdata(t);
    }

    
  }

  for (int i = 0; i < X/2; i++) {
    for (int j = 0; j < Y; j++) {
        tis.cell[i*Y + j].v = 0.0;
    }
  }

  // tmax here is adjusted to be the tmax from above, apdn.
  //tmax = adj;
  tmax = 3000;
  tnmax = tmax / dt;

  for (int tn = 0; tn < tnmax; tn++, ttn++) {
    tis.pace();

    if (tn % 30 == 0) {//write results every 3 ms
      double t = ttn*dt;
      rec.recdata(t);
    }
  }

  return 0;
}

