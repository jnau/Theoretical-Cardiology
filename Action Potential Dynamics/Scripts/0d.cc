#include <iostream>
#include <fstream>
using namespace std;

#include "cell.h"
#include "log.h"

double APDt1 = 0;

int main(void) {

  Logging lg;
  CCell cell;

  
  //double pcl = 400;//pacing cycle length (400 ms)
  double tmax = 0;
  int tnmax = 0;//convert to integer
  //int pcln = pcl / dt;//convert to integer
  int durn = 0;//duration of stimulation (2.0 ms)


  ofstream os("ks130.txt");
  //ofstream res("volts.txt");
  //ofstream apds("apd.txt");
  

  for (double i = 400; i > 169; i = i-1){
  	
    const int itr = 200;//the number of beats
    const double dt = 0.1;// time step (0.1 ms)

    tmax = i*itr;
    tnmax = tmax/dt;
    durn = 2.0 / dt;
  		
    int in = i/dt;
    double oldv = 0;
    double newv = 0;
  		
  		//main loop
  	for (int tn = 0; tn < tnmax; tn++) {

    	double stim = 0;//stimulation current
    	if (tn%in < durn) {
      		stim = 80.0;
    		}
    
    		oldv = cell.v;
    		cell.pace(stim);
        newv = cell.v;

        double aV = 0;
        double bV = 0;
        double APDt2 = 0;
        double di = 0;
        double apd = 0;
        
		
			if (newv > -81.5 && oldv < -81.5) {
				aV = (oldv-newv)/dt;
				bV = (oldv-aV*tn*dt);
				APDt1 = (-81.5-bV)/aV;	 
					
				}

    	else if (newv < -81.5 && oldv > -81.5) {
				aV = (oldv-newv)/dt;
				bV = (oldv-aV*tn*dt);
				APDt2 = (-81.5-bV)/aV;
					
				apd = APDt2 - APDt1; 
        di = i - apd; 
					

      				if (tn*dt > i*(itr-2)) {	 
               			 					
      					os << di << "," << apd << "," << i << endl;
    //            res << "Old v: " << oldv << "New v: " << newv << "Time: " << tn << "Pcl: " << i << endl;

      //          apds << "APD: " << apd << "APDt2: " << APDt2 << "APDt1 " << APDt1 << "Pcl: " << i << "Time :" << tn << endl;
      				}
     						

			} 

					

  
        

  			}
		}

  return 0;
}

