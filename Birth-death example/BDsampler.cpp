
// Add important stuff overall

#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <thread>
using namespace std;

// Structure to track MJP states in epidemic

struct BDstate {
	double time;
	int X;
};

// Include functions

#include "functions.h"


/***************/
/* Main method */
/***************/

int main(int argc, char** argv) 
{

	// BASICS REGARDING DATA AND ALGORITHM

	const int N = 50; //3000 ; //1000; //500; //200; //100; //50; //stoi(argv[1]);
	const double birthRate = 0.5; //30.0; // 10.0; //5.0; //2.0; // 1.0; //0.5; // trivial
	const double scaling = 1.5; // Only for uniformization procedure; //stold(argv[2]);
	const double K = 1.0; // 0.5;; //1.0; // Percentage of candidate jumps to add to any current trajectory; // stold(argv[3]);
	const double avDev =  pow(N,0.5); // pow(N,0.5); //  N/10.0; //  Average deviation from ODE approximation
	const double sd = 1.5 * (1.0 + 0.05); // 0.65 * (1.0 + 0.1); // For deviation in Aux vars, enough, since jump are of a single unit
	const double obsSd = N/50.0; // Measurement error for observations
	const int mcmcIter = 50000; // per chain //stold(argv[4]);
	const int chains = 8; //stold(argv[5]);
	const int seed = 152153; // 156416156; //stold(argv[6]);

	// Populate look-up table for exponentials at key-points

	for (int i = 0; i < 1000001; ++i)
	{
		expoFloors[i] = exp(-i/1000.0);
	}


	// READ DATA AND STORE INPUT VARIABLES

	vector<vector<double>> observations = fileToVector("data_50.txt");
	vector<vector<double>> odeApprox = fileToVector("ode_50.txt");


	/*******************/
	/* MCMC ITERATIONS */
	/*******************/
    
	// Number of cores

	int cores = thread::hardware_concurrency();
	cout << "\nNumber of threads found in your system: "<< cores << endl;
	cout << "Number of chains requested: "<< chains << endl;

	// Run the chains

	cout << "\n##########################################\n" << endl;
	cout << "Max population: "<< N << endl;
	cout << "Number of observations: "<< observations.size() << endl;
	cout << "Time interval: [0,"<< observations.back()[0] << "]"<< endl;

	cout << "\n##########################################\n" << endl;
	cout << "Sampling in progress:" << endl;

	if (cores >= chains){

		//Launch a group of threads

		thread t[chains];

		for (int i = 1; i <= chains; ++i) {
			mt19937 generator(seed*i);
			// t[i-1] = thread(do_mcmc_DepThinning, observations, N, obsSd, birthRate, mcmcIter, K, generator, "mcmc_chain_" + to_string(i) + ".csv",i);
			//// t[i-1] = thread(do_mcmc_DepThinningCounts, observations, N, obsSd, birthRate, mcmcIter, K, generator, "mcmc_chain_" + to_string(i) + ".csv",i);
			t[i-1] = thread(do_mcmc_AuxVar, observations,odeApprox, N, obsSd, birthRate, mcmcIter, K, avDev, sd, generator, "mcmc_chain_" + to_string(i) + ".csv",i);
			// t[i-1] = thread(do_mcmc_Unif, observations, N, obsSd, birthRate, mcmcIter, scaling, generator, "mcmc_chain_" + to_string(i) + ".csv",i);
		}

        //Join the threads with the main thread

		for (int i = 0; i < chains; ++i) {
			t[i].join();
		}

	} else {

		cout << "ERROR: You do not have enough threads for the chains requested; specify fewer chains" << endl;

	}

	cout << "\n\n";

    return 0;
}

