// UpdatedMain.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <new>
#include <cstring>
#include <climits>
#include <omp.h>

using namespace std;

// Pseudo-Random Number Generator Xorwow Algorithm
unsigned int xorwow(void) {

	static unsigned int x = 123456789;
	static unsigned int y = 362436069;
	static unsigned int z = 521288629;
	static unsigned int w = 88675123;
	static unsigned int v = 5783321;
	static unsigned int d = 6615241;
	unsigned int t;

	t = (x ^ (x >> 2));
	x = y;
	y = z;
	z = w;
	w = v;
	v = (v ^ (v << 4)) ^ (t ^ (t << 1));

	return (d += 362437) + v;
}

// Uniform Distribution Mapping
double uni_dist(double num) {
	double fixednumber = double(num) / UINT_MAX;
	return fixednumber; 
}

// Random Walk
void RandomWalk(double *positionx, double *positiony, double r_move) {
	double prand = uni_dist(xorwow());
	double radius = prand * 2 * M_PI; // random direction in polar coordinates

	// Given the direction and position of the ion in question, we calculate its random walk. 
	*positionx += r_move * cos(radius);
	*positiony += r_move * sin(radius);
}

// Check Boundary
void CheckBoundary(double *positionx, double *positiony, int *location, double sr_domainx, double sr_domainy, double c_domainx, double c_domainy) {

	if (*location == 0) {
		if (*positionx < -sr_domainx) {
			*positionx = -sr_domainx * 2 - *positionx;
		}
		else if (*positionx > sr_domainx) {
			*positionx = sr_domainx * 2 - *positionx;
		}

		if (*positiony < -sr_domainy) {
			*positiony = -sr_domainy * 2 - *positiony;
		}
		else if (*positiony > sr_domainy) {
			*positiony = sr_domainy * 2 - *positiony;
		}
	}
	else if (*location != 0) {

		if (*positionx < -c_domainx) {
			*positionx = -c_domainx * 2 - *positionx;
		}
		else if (*positionx > c_domainx) {
			*positionx = c_domainx * 2 - *positionx;
		}

		if (*positiony < -c_domainy) {
			*positiony = -c_domainy * 2 - *positiony;
		}
		else if (*positiony > c_domainy) {
			*positiony = c_domainy * 2 - *positiony;
		}

	}
}

double Ions2Concentration(int ions, double volume) {

	double AvoConstant = 6.022;
	double ConCoeff = volume * AvoConstant * 100;
	double concentration = ions / ConCoeff;

	return concentration;
}

double Concentration2Ions(double concentration, double volume) {

	double AvoConstant = 6.022;
	double ConCoeff = volume * AvoConstant * 100;
	double ions = concentration * ConCoeff;

	return ions;
}

void circshift(int *grid, int *NewGrid, int x, int y, int size) {

	for (int i = 0; i < size; i++) {

		for (int j = 0; j < size; j++) {

			if (*((grid + i * size) + j) != 0) {
				*((NewGrid + (i + y)*size) + (j + x)) = *((grid + i * size) + j);
			}
		}
	}
}

void addSixGrids(int *Newgrid, int *G1, int *G2, int *G3, int *G4, int *G5, int *G6, int size) {

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			*((Newgrid + i * size) + j) = *((G1 + i * size) + j) + *((G2 + i * size) + j) + *((G3 + i * size) + j) + *((G4 + i * size) + j) + *((G5 + i * size) + j) + *((G6 + i * size) + j);
		}
	}

}

void SingleCluster(const int numRyRs, double radius_influence, double *Xpos, double *Ypos) {

	const int gridsize = 100;
	int grid[gridsize][gridsize];

	for (int i = 0; i < gridsize; i++) {
		for (int j = 0; j < gridsize; j++) {
			grid[i][j] = 0;
		}
	}

	grid[gridsize / 2][gridsize / 2] = 1;

	double prob_g = 0.025;
	double prob_d = 0.001; // probability to delete isolated ryr
	double prob_r = 0.999;

	int sum = 0;

	while (sum != numRyRs) {

		int grid1[gridsize][gridsize];
		int grid2[gridsize][gridsize];
		int grid3[gridsize][gridsize];
		int grid4[gridsize][gridsize];
		int grid5[gridsize][gridsize];
		int grid6[gridsize][gridsize];

		for (int i = 0; i < gridsize; i++) {
			for (int j = 0; j < gridsize; j++) {
				grid1[i][j] = 0;
				grid2[i][j] = 0;
				grid3[i][j] = 0;
				grid4[i][j] = 0;
				grid5[i][j] = 0;
				grid6[i][j] = 0;
			}
		}

		circshift((int *)grid, (int *)grid1, 0, 1, gridsize);
		circshift((int *)grid, (int *)grid2, 0, -1, gridsize);
		circshift((int *)grid, (int *)grid3, 1, 0, gridsize);
		circshift((int *)grid, (int *)grid4, -1, 0, gridsize);
		circshift((int *)grid, (int *)grid5, -1, -1, gridsize);
		circshift((int *)grid, (int *)grid6, -1, 1, gridsize);

		int prob_growth[gridsize][gridsize];

		addSixGrids((int *)prob_growth, (int *)grid1, (int *)grid2, (int *)grid3, (int *)grid4, (int *)grid5, (int *)grid6, gridsize);

		// while going through prob_g, generate a random number and compare. if r is less, then add 1 to original grid.

		for (int i = 0; i < gridsize; i++) {
			for (int j = 0; j < gridsize; j++) {
				double rand_num = xorwow();
				double prand = uni_dist(rand_num);
				if (prand < prob_growth[i][j] * prob_g) {
					grid[i][j] = 1;
				}

				grid1[i][j] = 0;
				grid2[i][j] = 0;
				grid3[i][j] = 0;
				grid4[i][j] = 0;
				grid5[i][j] = 0;
				grid6[i][j] = 0;

			}
		}

		circshift((int *)grid, (int *)grid1, 0, 1, gridsize);
		circshift((int *)grid, (int *)grid2, 0, -1, gridsize);
		circshift((int *)grid, (int *)grid3, 1, 0, gridsize);
		circshift((int *)grid, (int *)grid4, -1, 0, gridsize);
		circshift((int *)grid, (int *)grid5, -1, -1, gridsize);
		circshift((int *)grid, (int *)grid6, -1, 1, gridsize);

		int prob_delete[gridsize][gridsize];
		addSixGrids((int *)prob_delete, (int *)grid1, (int *)grid2, (int *)grid3, (int *)grid4, (int *)grid5, (int *)grid6, gridsize);

		for (int i = 0; i < gridsize; i++) {
			for (int j = 0; j < gridsize; j++) {
				double rand_num = xorwow();
				double prand = uni_dist(rand_num);
				if (prand < (6 - prob_delete[i][j])*prob_d) {
					grid[i][j] = 0;
				}
			}
		}

		for (int i = 0; i < gridsize; i++) {
			for (int j = 0; j < gridsize; j++) {
				double rand_num = xorwow();
				double prand = uni_dist(rand_num);
				if (prand > prob_r) {
					grid[i][j] = 0;
				}
			}
		}

		grid[gridsize / 2][gridsize / 2] = 1;

		int count = 0;

		for (int i = 0; i < gridsize; i++) {
			for (int j = 0; j < gridsize; j++) {
				if (grid[i][j] != 0) {

					count++;
				}
			}
		}

		sum = count;
		while (sum > numRyRs) {

			for (int i = 0; i < gridsize; i++) {
				for (int j = 0; j < gridsize; j++) {
					double rand_num = xorwow();
					double prand = uni_dist(rand_num);
					if (prand > prob_r) {
						grid[i][j] = 0;
						grid[gridsize / 2][gridsize / 2] = 1;
					}
				}
			}

			int count = 0;

			for (int i = 0; i < gridsize; i++) {
				for (int j = 0; j < gridsize; j++) {
					if (grid[i][j] != 0) {

						count++;
					}
				}
			}

			sum = count;
		}

		count = 0;

		for (int i = 0; i < gridsize; i++) {
			for (int j = 0; j < gridsize; j++) {
				if (grid[i][j] != 0) {

					count++;
				}
			}
		}

		sum = count;

	}

	int count = 0;
	for (int px = 0; px < gridsize; px++) {
		for (int py = 0; py < gridsize; py++) {

			if (grid[px][py] == 1) {

				if ((py % 2) == 1) {
					*(Xpos + count) = 2 * radius_influence*(px - (gridsize / 2)) + radius_influence;
				}
				else {
					*(Xpos + count) = 2 * radius_influence*(px - (gridsize / 2));
				}

				*(Ypos + count) = sqrt(3)*radius_influence*(py - (gridsize / 2));
				count++;
			}
		}
	}
}

int CaRelease(int open, const int numRyRs, double sr, double cyto, double sr_volume, double cleft_volume, int opentime) {
	const double ks = 0.005; // us
	// (open/(double)NumRyRs)
	double Stochastic_NumMove = Concentration2Ions((ks *(open/((double)numRyRs))*(sr - cyto)*(cleft_volume / sr_volume)), sr_volume);
	//cout << sr << '\t' << cyto << '\t' << (cleft_volume / sr_volume) << '\t' << ks * (open / (double)numRyRs) << endl;
	int numMove = floor(Stochastic_NumMove);
	//cout << Stochastic_NumMove << '\t' << numMove << endl;
	Stochastic_NumMove = Stochastic_NumMove - numMove;
	
	double prand = uni_dist(xorwow());

	if (prand < Stochastic_NumMove) {
		numMove++;
	}

	if (Stochastic_NumMove < 0) {
		numMove = 0;
	}
	//cout << sr << '\t' << cyto << '\t' << (ks * (open / (double)numRyRs)*(sr- cyto)*(cleft_volume / sr_volume)) << '\t' << numMove << endl;
	
	return numMove;
}


void LocateMoveIon(double *RyRx, double *RyRy, int *Location, double *XIons, double *YIons, const int numRyRs, double radius_influence, int *state, int numMove, int *count) {
	for (int i = 0; i < numRyRs; i++) {

		double dist = (*XIons - (*(RyRx + i))) * (*XIons - (*(RyRx + i))) + (*YIons - (*(RyRy + i))) * (*YIons - (*(RyRy + i)));

		if (dist < radius_influence*radius_influence) {
			*Location = 1;
#pragma omp atomic
			(*count)++;
		}
	}
}

void Uptake(int *Location, double uptake, const int numIons, double c_rangex, double c_rangey, double sr_rangex, double sr_rangey, int sumCyto, double* XIon, double* YIon) {

	double prand = uni_dist(xorwow());
	vector<int> v;
	std::vector<int>::iterator it;

	for (int i = 0; i < numIons; i++) {

		if (*(Location + i) == 1) {
			v.push_back(i);
		}
	}

	random_shuffle(v.begin(), v.end());

	double cyto_volume = 2 * c_rangex * 2 * c_rangey*0.012;
	double up_ions = Concentration2Ions(uptake, cyto_volume);
	int total = up_ions;
	double diff = up_ions - total;

	if (prand < diff) {
		total++;
	}

	if (total < sumCyto) {
		for (int i = 0; i < total; i++) {

			double prandx = uni_dist(xorwow());
			double prandy = uni_dist(xorwow());

			*(XIon + v[i]) = prandx * 2 * sr_rangex - sr_rangex;
			*(YIon + v[i]) = prandy * 2 * sr_rangey - sr_rangey;
			*(Location + v[i]) = 0;
		}
	}
	else {
		for (int i = 0; i < v.size(); i++) {

			double prandx = uni_dist(xorwow());
			double prandy = uni_dist(xorwow());

			*(XIon + v[i]) = prandx * 2 * sr_rangex - sr_rangex;
			*(YIon + v[i]) = prandy * 2 * sr_rangey - sr_rangey;
			*(Location + v[i]) = 0;
		}

	}

	v.erase(v.begin(), v.end());
}

void Buffer(int sumCyto, int sumBound, double volume, double *TroponinBuff, const int numIons, int *loc) {

	double k2_f = 0.0001;
	double k2_b = 0.0001;
	int total_b = 123;
	double prand = uni_dist(xorwow());

	double Cyto_c = Ions2Concentration(sumCyto, volume);
	double d_TrpnBuff = (k2_f*Cyto_c*(total_b - (*TroponinBuff))) - (k2_b*(*TroponinBuff));
	int d_CaIons = floor(-(Concentration2Ions(d_TrpnBuff, volume)));
	double NumSIons = abs(-(Concentration2Ions(d_TrpnBuff, volume)) - d_CaIons);
	if (prand < NumSIons) {
		d_CaIons++;
	}

	*TroponinBuff = *TroponinBuff - Ions2Concentration(d_CaIons, volume);
	if (d_CaIons < 0) {

		vector<int> v;
		std::vector<int>::iterator it;
		int count = 0;

		for (int i = 0; i < numIons; i++) {
			if (*(loc + i) == 1) {
				v.push_back(i);
				count++;
			}
		}
		random_shuffle(v.begin(), v.end());

		if (abs(d_CaIons) <= sumCyto) {

			for (int i = 0; i < abs(d_CaIons); i++) {
				*(loc + v[i]) = 2;
			}
		}

		else {
			for (int j = 0; j < count; j++) {
				*(loc + v[j]) = 2;
			}
		}

	}
	else if (d_CaIons > 0) {
		vector<int> v;
		std::vector<int>::iterator it;
		int count = 0;

		for (int i = 0; i < numIons; i++) {
			if (*(loc + i) == 2) {
				v.push_back(i);
				count++;
			}
		}
		random_shuffle(v.begin(), v.end());

		if (d_CaIons < count) {
			for (int j = 0; j < d_CaIons; j++) {
				*(loc + v[j]) = 1;
			}
		}
		else if (d_CaIons > count) {
			for (int j = 0; j < count; j++) {
				*(loc + v[j]) = 1;
			}
		}
	}
}


inline double k_CaSR(double concentration) {
	int SR_Max = 15;
	int SR_Min = 1;
	int ECSR = 450; // uM
	double H = 2.5;

	return (SR_Max - ((SR_Max - SR_Min) / (1 + pow(ECSR / concentration, H))));
}

double k0SR(double concentration) {
	//double k0Ca = 0.0000002; // uM-2 us-1
	double k0Ca = 0.0000001;
	return (k0Ca / k_CaSR(concentration));
}

double kiSR(double concentration) {
	//double kiCa = 0.00000005; // uM-1 us-1
	double kiCa = 0.000000005;
	return (kiCa*(k_CaSR(concentration)));
}

// Four State Shannon Model
void CalculateState(int *state, double k0, double ki, double SumCyto) {


	double kim = 0.000005;
	double kom = 0.00006;
	double dt = 1;

	if (*state == 0) {
		double prand = uni_dist(xorwow());
		if (prand < (k0*SumCyto*SumCyto)*dt) {
			*state = 1;
			//cout << *prob << '\t' << "State 0" << endl;
			//*prob = *prob + 1;
		}
		else if (prand < ((ki*SumCyto) + (k0*SumCyto*SumCyto))*dt) {
			*state = 3;
		}
	}
	else if (*state == 1) {
		double prand = uni_dist(xorwow());
		//(*prob)++;
		if (prand < (ki*SumCyto)*dt) {
			*state = 2;
		}
		else if (prand < (kom + (ki*SumCyto))*dt) {
			*state = 0;
		}
	}
	else if (*state == 2) {
		double prand = uni_dist(xorwow());
		if (prand < kim*dt) {
			*state = 1;
		}
		else if (prand < (kom + kim)*dt) {
			*state = 3;
		}
	}
	else if (*state == 3) {
		double prand = uni_dist(xorwow());
		if (prand < kim*dt) {
			*state = 0;
		}
		else if (prand < ((k0*SumCyto*SumCyto) + kim)*dt) {
			*state = 2;
		}
	}
}

void RecordConcentration(double *XIons, double *YIons, int *Location, double sr_volume, double c_volume, double d_volume, double c_xrange, double c_yrange, double dyadradius, const int numberIons, int time, ostream& File_SR, ostream& F_Cyto, ostream& F_local, ostream& F_buff, ostream& F_dyad) {

	const int sizex = 50;
	const int sizey = 50;

	int Cytosol_count[sizex][sizey] = {};
	double cyto_concentration[sizex*sizey] = {};

	for (int i = 0; i < sizex; i++) {
		for (int j = 0; j < sizey; j++) {
			Cytosol_count[i][j] = 0;
		}
	}

	double c_dx = (2 * c_xrange) / sizex;
	double c_dy = (2 * c_yrange) / sizey;

	int numbercytosol = 0;
	int numSR = 0;
	int numberbound = 0;
	int sumDyad = 0;

	for (int k = 0; k < numberIons; k++) {

		int xid = 0;
		int yid = 0;

		int lid = Location[k];

		if (lid == 1) {
			xid = int((XIons[k] - (-c_xrange)) / c_dx);
			yid = int((YIons[k] - (-c_yrange)) / c_dy);
			if (xid < 0 || xid >= sizex || yid < 0 || yid >= sizey) { break; }
		}

		switch (lid) {
		case 0:
			numSR++;
			break;
		case 1:
			Cytosol_count[xid][yid]++;
			numbercytosol++;
			break;
		case 2:
			numberbound++;
			break;
		}

		if (lid == 1) {
			double distx = (XIons[k] * XIons[k]);
			double disty = (YIons[k] * YIons[k]);
			if (distx + disty < dyadradius*dyadradius) {
#pragma omp atomic
				sumDyad++;
			}
		}
	}
	/*
	if (time % 1000 == 0) {

		for (int i = 0; i < sizex; i++) {
			for (int j = 0; j < sizey; j++) {

				double cgrid_volume = 2 * c_dx * 2 * c_dy*0.012;
				cyto_concentration[sizey*(i)+j] = Ions2Concentration(Cytosol_count[i][j], cgrid_volume);
			}
		}

		F_Cyto.write((char *)cyto_concentration, sizeof(double)*sizex*sizey);
	}
	
	*/
	

	double localbuff = Ions2Concentration(numberbound, c_volume);
	double localconcentration = Ions2Concentration(numbercytosol, c_volume);
	double dyadconcentration = Ions2Concentration(sumDyad, d_volume);
	double localSR = Ions2Concentration(numSR, sr_volume);

	File_SR.write((char *)&localSR, sizeof(double));
	F_local.write((char *)&localconcentration, sizeof(double));
	F_buff.write((char *)&localbuff, sizeof(double));
	F_dyad.write((char *)&dyadconcentration, sizeof(double));

}

void RecordPosition(double XIons[], double YIons[], int Location[], const int numIons, ostream& xpos, ostream& ypos, float X[], float Y[]) {

	for (int i = 0; i < numIons; i++) {

		if (Location[i] == 1) {
			*(X + i) = floor(XIons[i] * 100 + 0.5) / 100;
			*(Y + i) = floor(YIons[i] * 100 + 0.5) / 100;
			//cout << *(X + i) << '\t' << *(Y + i) << endl;

		}
		else {
			*(X + i) = 4;
			*(Y + i) = 4; 
		}
	}

	xpos.write((char *)X, sizeof(float)*numIons);
	ypos.write((char *)Y, sizeof(float)*numIons);
}

int main()
{

	// Initialize Constants
	double sr_xrange = 1; // um
	double sr_yrange = 1; // um
	double sr_volume = 4 * sr_xrange*sr_yrange*0.048;
	double cyto_xrange = 3; //um
	double cyto_yrange = 3; //um
	double cyto_volume = 4 * cyto_xrange*cyto_yrange*0.012;
	const int gridsizex = 50;
	const int gridsizey = 50;
	const int numIons = 80000; // 14500 ~500uM,  700uM =  20k
	const int TotalTime = 1000000; // us
	const int numRyRs = 50;
	//const double dyad_Xsize = 0.16;
	//const double dyad_Ysize = 0.16;
	double dyadradius = 0.183;
	double dyadheight = 0.012;
	int tau_uptake = 1000; // us
	int dt = 1;
	int OpenTime = 1;
	int Release = 0;
	int Prev_SumSR = 0;
	int numMove = 0;
	int OpenCount = 0;
	double radius_influence = 0.015;
	double dyadvolume = M_PI * dyadradius*dyadradius*dyadheight;
	//double volume_ryr = (4 / 3)*(M_PI*radius_influence*radius_influence*radius_influence);
	double uptake = 0.3;
	double TroponinBuff = 0;
	double OpenProb = 0;
	double r_move = sqrt(223 * 4 * dt*pow(10, -6)); // um, length of a random walk in dyadic space in 2D.

	ofstream ofs_OpenRyR("OpenRyR.txt");
	//ofstream ofs_Xpos("posx.dat", ios::binary);
	//ofstream ofs_Ypos("posy.dat", ios::binary);
	ofstream ofs_RyRLocation("RyRLocation.dat", ios::binary);
	ofstream ofs_SR("LocalSR.dat", ios::binary);
	ofstream ofs_Cytosol("Cyto_concentration.dat", ios::binary);
	ofstream ofs_ion("IonInfo.txt");
	ofstream ofs_LocalCyto("LocalCytosol.dat", ios::binary);
	ofstream ofs_State("State.dat", ios::binary);
	ofstream ofs_Buff("Buffer.dat", ios::binary);
	ofstream ofs_NumState("StateFreq.txt");
	ofstream ofs_dyad("dyad.dat", ios::binary);
	ofstream ofs_flux("Flux.txt");

	double XIons[numIons] = {};
	double YIons[numIons] = {};
	//float RecordX[numIons] = {};
	//float RecordY[numIons] = {};
	double RyR_posx[numRyRs] = {};
	double RyR_posy[numRyRs] = {};
	int state[numRyRs] = {};
	int LocationIons[numIons] = {};
	//int RyRCount[numRyRs] = {};

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < numIons; i++) {

		double prandx = uni_dist(xorwow());
		double prandy = uni_dist(xorwow());

		XIons[i] = prandx * 2 * sr_xrange - sr_xrange;
		YIons[i] = prandy * 2 * sr_yrange - sr_yrange;
		LocationIons[i] = 0;
		//RecordX[i] = 0;
		//RecordY[i] = 0;
	}

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < numRyRs; i++) {	
		state[i] = 0;
	}

	state[0] = 1;

	SingleCluster(numRyRs, radius_influence, RyR_posx, RyR_posy);
	ofs_RyRLocation.write((char *)RyR_posx, sizeof(double)*numRyRs);
	ofs_RyRLocation.write((char *)RyR_posy, sizeof(double)*numRyRs);
	ofs_RyRLocation.close();

	for (int time = 0; time < TotalTime; time++) {

			RecordConcentration(XIons, YIons, LocationIons, sr_volume, cyto_volume, dyadvolume, cyto_xrange, cyto_yrange, dyadradius, numIons, time, ofs_SR, ofs_Cytosol, ofs_LocalCyto, ofs_Buff, ofs_dyad);
			
			//if (time % 100 == 0) {
				//RecordPosition(XIons, YIons, LocationIons, numIons, ofs_Xpos, ofs_Ypos, RecordX, RecordY);	
			//}

			int open = 0; // 1
			int close_cyto = 0; // 0
			int inactive = 0; // 3
			int close_SR = 0; //2

#pragma omp parallel for
			for (int k = 0; k < numRyRs; k++) {
				switch (state[k]) {
				case 0:
#pragma omp atomic
					close_cyto++;
					break;
				case 1:
#pragma omp atomic
					open++;
					break;

				case 2:
#pragma omp atomic
					close_SR++;
					break;
				case 3:
#pragma omp atomic
					inactive++;
					break;
				}
			}

			ofs_OpenRyR << open << endl;
			ofs_NumState << close_cyto << '\t' << open << '\t' << close_SR << '\t' << inactive << endl;
			ofs_State.write((char *)state, sizeof(int)*numRyRs);

			int SumDyad = 0;
			int SumSR = 0;

#pragma omp parallel for
			for (int i = 0; i < numIons; i++) {

				if (LocationIons[i] == 0) {
#pragma omp atomic 
					SumSR++;
				}

				else if (LocationIons[i] == 1) {

					double distx = (XIons[i] * XIons[i]);
					double disty = (YIons[i] * YIons[i]);
					if (distx+disty < dyadradius*dyadradius) {
#pragma omp atomic
						SumDyad++;
					}
				}
			}

			double DyadConcentration = Ions2Concentration(SumDyad, dyadvolume);
			double SRConcentration = Ions2Concentration(SumSR, sr_volume);
			
			if (open > 0){
				OpenTime++;
			} else {
				OpenTime = 1;
			}

			numMove = CaRelease(open, numRyRs, SRConcentration, DyadConcentration, sr_volume, dyadvolume, OpenTime);

			SumDyad = 0;
			SumSR = 0;
			int count = 0;
			int SumCyto = 0;
			int SumBound = 0;

#pragma omp parallel 
			{

#pragma omp for schedule(dynamic)
				for (int i = 0; i < numIons; i++) {
					RandomWalk(&XIons[i], &YIons[i], r_move);
					CheckBoundary(&XIons[i], &YIons[i], &LocationIons[i], sr_xrange, sr_yrange, cyto_xrange, cyto_yrange);
					
					if (LocationIons[i] == 0 && count < numMove) {
						LocateMoveIon(RyR_posx, RyR_posy, &LocationIons[i], &XIons[i], &YIons[i], numRyRs, radius_influence, state, numMove, &count);
					}
					
					switch (LocationIons[i]) {
					case 0:
#pragma omp atomic
						SumSR++;
						break;
					case 1:
#pragma omp atomic 
						SumCyto++;
						break;

					case 2:
#pragma omp atomic
						SumBound++;
						break;
					}

					if (LocationIons[i] == 1) {

						double distx = (XIons[i] * XIons[i]);
						double disty = (YIons[i] * YIons[i]);
						if (distx + disty < dyadradius*dyadradius) {
#pragma omp atomic
							SumDyad++;
						}
					}
				}
			}
			
			double SR_con = Ions2Concentration(SumSR, sr_volume);
			double Dyad_con = Ions2Concentration(SumDyad, dyadvolume);
			double k0 = k0SR(SR_con);
			double ki = kiSR(SR_con);

#pragma omp parallel for
			for (int j = 0; j < numRyRs; j++) {
				//double Cyto_con = Ions2Concentration(RyRCount[j], volume_ryr);
				CalculateState(&state[j], k0, ki, Dyad_con);
				//RyRCount[j] = 0;
			}

			Release = Release + count;

			if (time % 1000 == 0) {
				Uptake(LocationIons, uptake, numIons, cyto_xrange, cyto_yrange, sr_xrange, sr_yrange, SumCyto, XIons, YIons);
				ofs_flux << Ions2Concentration(Release, sr_volume) << '\t' << Release << endl;
				Release = 0;
			}

			SumCyto = 0;
			#pragma omp parallel for
			for (int ion = 0; ion < numIons; ion++) {
				if (LocationIons[ion] == 1) {
				#pragma omp atomic
					SumCyto++;
				}
			}

			Buffer(SumCyto, SumBound, cyto_volume, &TroponinBuff, numIons, LocationIons);
	}

	ofs_SR.close();
	ofs_LocalCyto.close();
	ofs_NumState.close();
	ofs_State.close();
	ofs_OpenRyR.close();
	ofs_flux.close();
	ofs_Cytosol.close();

	ofs_ion << numIons << '\t' << TotalTime << '\t' << numRyRs << endl;
	ofs_ion << sr_xrange << '\t' << sr_yrange << '\t' << cyto_xrange << endl;
	ofs_ion << cyto_yrange << '\t' << gridsizex << '\t' << gridsizey << endl;

	ofs_ion.close();
	return 0;
}