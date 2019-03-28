
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

struct mjpState {
	double time;
	int S, I, R, O;
};

// Include functions

#include "functions.h"


/***************/
/* Main method */
/***************/

int main(int argc, char** argv) 
{

	// BASICS REGARDING DATA AND ALGORITHM

	const int population = 50; //stoi(argv[1]);
	const int initials = 1; 
	const double scaling = 1.5; //unif 1.0 DP 1.5; //stold(argv[2]);
	const int mcmcIter = 50000; // per chain //stold(argv[4]);
	const int chains = 8; //stold(argv[5]);
	const int seed = 152153; // 156416156; //stold(argv[6]);

	// Populate look-up table for exponentials at key-points

	for (int i = 0; i < 1000001; ++i)
	{
		expoFloors[i] = exp(-i/1000.0);
	}

	// READ DATA AND STORE INPUT VARIABLES

	vector<double> removalData = fileToVector("removals50.txt");
	vector<vector<double>> odeApprox = fileToVector2("ode_50.txt");

	/*******************/
	/* MCMC ITERATIONS */
	/*******************/
    
	// Number of cores

	int cores = thread::hardware_concurrency();
	cout << "\nNumber of threads found in your system: "<< cores << endl;
	cout << "Number of chains requested: "<< chains << endl;

	// Run the chains

	cout << "\n##########################################\n" << endl;
	cout << "Population: "<< population << endl;
	cout << "Final full amount of removals: "<< removalData.size() << endl;

	cout << "\n##########################################\n" << endl;
	cout << "Sampling in progress:" << endl;

	if (cores >= chains){

		//Launch a group of threads

		thread t[chains];

		for (int i = 1; i <= chains; ++i) {
			mt19937 generator(seed*i);
			t[i-1] = thread(do_mcmc_AuxVar, removalData, odeApprox, population, initials, mcmcIter, scaling, generator, "mcmc_chain_" + to_string(i) + ".csv",i);
			// t[i-1] = thread(do_mcmc_DP, removalData, population, initials, mcmcIter, scaling, generator, "mcmc_chain_" + to_string(i) + ".csv",i);
			// t[i-1] = thread(do_mcmc_Unif, removalData, population, initials, mcmcIter, scaling, generator, "mcmc_chain_" + to_string(i) + ".csv",i);
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

