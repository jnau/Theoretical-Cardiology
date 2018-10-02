/*
CaDiffusion.cpp

This C++ script is meant to implement a calcium wave in one dimension using a stochastic process.
We will be using random walks to simulate a wave with similarities to the diffusion PDE. 
Since random walks are the fundementals of pdes, we want to study this behavior. 

*/

// Preprocessor Directives
#include <iostream>
#include <sstream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <random>
#include <new>
using namespace std;

// Probability Function
double fprob (double n)

{
	double pmin, pmax, pval, phalf, hc;

	pmin = 0.0001;
	pmax = 0.3;
	phalf = 100;
	hc = 10;
	pval = (double) pmin + (double)(pmax*(pow(n, hc)))/ ((double)(pow(n, hc)) + (pow(phalf, hc)));
	return pval;

}

// Random Walk simulation
int walk (double p)
{
	int w;

	if (p <= 0.5000){

		w = 1;
	}

	else {

		w = -1;
	}

	return w;

}


int main()
{

	// Initialization
	int len, ions, boundary, T, n_op, div;
	double incre, tau_op, tau_in;

	ions = 100000;
	boundary = 100;
	T = 10000;
	incre = (double) (2*boundary)/ (double) ions;
	n_op = 20;
	div = (2*boundary)/n_op;
	len = n_op+1;

	double * sr_ions, * cyto_ions;
	sr_ions = new double [ions];
	cyto_ions = new double [ions];

	// outfiles
	ofstream positions("s_positions_new.txt");
	ofstream positionc("c_positions_new.txt");
	ofstream State("States_100.txt"); 
	ofstream Location("Location.txt");

	for (int i = 0; i <= ions; i++) {

		sr_ions[i] = round((-boundary) + incre*i);
		cyto_ions[i] = -1000;
	
	}
	
	double * cru, * cru_state, * t_open, * t_inactive;
	cru = new double [len];
	cru_state = new double [len];
	t_open = new double [len];
	t_inactive = new double [len];
	
	tau_op = 50;
	tau_in = 200;

	// Initializing CRU locations and state
	for (int i = 0; i < len; i++) {

		cru[i] = -boundary + (div*i);
		cru_state[i] = 0;
		t_open[i] = 0;
		t_inactive[i] = 0;

		Location << cru[i] << endl;

	}

	// force open one cru to begin simulation
	cru_state[0] = 1;
	t_open[0] = 25;

	default_random_engine gen;
	uniform_real_distribution<double> dis(0,1);
	
	for (int time = 0; time < T; time++) {

		for (int j = 0; j < len; j++) {
			
			double prand = dis(gen);
			double num_cyions = 0;

			// Determine number of ions within range of influence
			for (int k = 0; k <= ions; k++) {

				if ((cyto_ions[k] > (cru[j]-div)) && (cyto_ions[k] < cru[j]+div)) {

					num_cyions++;

				}

			} 

			// Check state of calcium release units (CRUs), 0 close, 1 open, -1 inactive
			if (cru_state[j] == 0) {

				double actualProb = fprob(num_cyions);

				if (prand < actualProb) {
					cru_state[j] = 1;
					t_open[j] = time;
				}
			}

			else if (cru_state[j] ==1) {

				//cout<< cru[j]<< endl;
				if (time >= (t_open[j] + tau_op)) {
					cru_state[j] = -1;
					t_inactive[j] = time;
				}
			}

			else if (cru_state[j] == -1) {

				if (time >= (t_inactive[j] + tau_in)) {
					cru_state[j] = 0;
				}
			}		
		}

		// Random Walk Simulation
		// make for loop, generate random -1, or 1, and add to Cytosol (c), or Sarcoplasmic Reticulum (s). 

		for (int j = 0; j<= ions; j++) {

			if (cyto_ions[j] > -1000) {

				double rand_walk = dis(gen);
				int w = walk(rand_walk);
				cyto_ions[j] = cyto_ions[j] + w;


				if (cyto_ions[j] < -boundary) {

					cyto_ions[j] = -boundary + abs(-boundary-(cyto_ions[j]));
				}

				else if (cyto_ions[j] > boundary) {

					cyto_ions[j] = boundary - abs(cyto_ions[j] - boundary);

				}
			}

			if (sr_ions[j] > -1000) {

				double rand_walk = dis(gen);
				int w = walk(rand_walk);

				sr_ions[j] = sr_ions[j] + w;

				if (sr_ions[j] < -boundary) {
					sr_ions[j] = -boundary + abs(-boundary-(sr_ions[j]));
				}

				else if (sr_ions[j] > boundary) {
					sr_ions[j] = boundary - abs(sr_ions[j] - boundary);
				}
			}
		}

		for (int j = 0; j<len; j++) {

			for (int k = 0; k <= ions; k++) {

				if ((cru[j] == sr_ions[k]) && (cru_state[j]==1)) {

					cyto_ions[k] = sr_ions[k];
					sr_ions[k] = -1000;

				}
			}
		}

		// uptake
		for (int u = 0; u <= ions; u++) {

			double r = dis(gen);

			if ((r<0.01) && (cyto_ions[u] > -1000)) {

				sr_ions[u] = cyto_ions[u];
				cyto_ions[u] = -1000;
			}
		}

		// save data every 5 us
		for (int j = 0; j <= ions; j++) {

			if (time%5==0) {

				positions << sr_ions[j] << '\t';
				positionc << cyto_ions[j] << '\t';

			}

		}

		if (time%5==0) {
			positions << endl;
			positionc << endl;
		}

		for (int j = 0; j <len; j++) {

			if (time%5==0) {

				State << cru_state[j] << '\t';
			}

		}

		if (time%5==0) {
			State << endl;
			
		}
	}
	
	return 0;

}