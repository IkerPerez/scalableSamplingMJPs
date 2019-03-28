
/*****************************/
/* Stuff to read a text file */
/*****************************/

vector<double> fileToVector(string filename)
{
	std::vector<double> aVector;

	string line;
	ifstream myfile (filename);
	if (myfile.is_open())
	{
		int i = 0;
		while ( getline (myfile,line) )
		{
			aVector.push_back(stold(line));
		}
		myfile.close();
	}

	else cout << "Unable to open file\n"; 

	return aVector;  
}

vector<vector<double>> fileToVector2(string filename)
{
	vector<vector<double>> aVector;

	string line;
	ifstream myfile (filename);
	if (myfile.is_open())
	{
		while ( getline (myfile,line,' ') )
		{
			vector<double> dummy;
			dummy.push_back(stold(line));
			getline (myfile,line);
			dummy.push_back(stold(line));

			aVector.push_back(dummy);
		}
		myfile.close();
	}

	else cout << "Unable to open file\n"; 

	return aVector;  
}


/********************************/
/* Some functions/macros of use */
/********************************/

double ipow(double base, int exp) 
{
    double result = 1.0;
    while (exp > 0) {
        if (exp & 1) result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

#define i0max(i) ((i > 0) ? i : 0)

// Look-up tables for exponentials over range -0:1000; note rates always negative

double expoFloors[1000001];
#define fastExp(x) ((x > -1000) ? expoFloors[ (int) (-(x)*1000.0) ] : 0) // enough quality approximtion here


/******************************/
/* To produce a starting Path */
/******************************/

void mcmc_startingPathVector(vector<mjpState>& startingPath, vector<double>& removalData, int pop, int init, mt19937& gen)
{
	// Create and fill a frame to operate within
	startingPath.push_back({0.0,pop-init, init,0,0}); 
	for (int i = 0; i < removalData.size(); ++i){
		startingPath.push_back({removalData[i],0,0,0,3});
	} 

	// Associated uniform infection times
	for (int i = 0; i < removalData.size()-init; ++i){
		uniform_real_distribution<double> unif(0.0, removalData[i]);
		startingPath.push_back({unif(gen),0,0,0,2});
	} 

	// Sort according to times
	sort(startingPath.begin(), startingPath.end(), [](mjpState & i1, mjpState & i2) {return i1.time < i2.time;});

	// Populate starting path and organise data to feed into later procedures
	for (int i = 1; i < startingPath.size(); ++i)
	{
		if (startingPath[i].O == 2){
			startingPath[i].S = startingPath[i-1].S - 1;
			startingPath[i].I = startingPath[i-1].I + 1;
			startingPath[i].R = startingPath[i-1].R;
			startingPath[i].O = 0;
		} else if (startingPath[i].O == 3){
			startingPath[i].S = startingPath[i-1].S;
			startingPath[i].I = startingPath[i-1].I - 1;
			startingPath[i].R = startingPath[i-1].R + 1;
			startingPath[i].O = 1;
		}
	}

}

void mcmc_startingPath(mjpState * startingPath, vector<double>& removalData, int pop, int init, mt19937& gen)
{
	// Create and fill a frame to operate within
	startingPath[0] = {0.0,pop-init, init,0,0}; 
	for (int i = 0; i < removalData.size(); ++i){
		startingPath[i+1] = {removalData[i],0,0,0,3};
	} 

	// Associated uniform infection times
	for (int i = 0; i < removalData.size()-init; ++i){
		uniform_real_distribution<double> unif(0.0, removalData[i]);
		startingPath[removalData.size()+1+i] = {unif(gen),0,0,0,2};
	} 

	// Sort according to times
	sort(startingPath, startingPath + 2*removalData.size(), [](mjpState & i1, mjpState & i2) {return i1.time < i2.time;});

	// Populate starting path and organise data to feed into later procedures
	for (int i = 1; i < 2*removalData.size(); ++i)
	{
		if (startingPath[i].O == 2){
			startingPath[i].S = startingPath[i-1].S - 1;
			startingPath[i].I = startingPath[i-1].I + 1;
			startingPath[i].R = startingPath[i-1].R;
			startingPath[i].O = 0;
		} else if (startingPath[i].O == 3){
			startingPath[i].S = startingPath[i-1].S;
			startingPath[i].I = startingPath[i-1].I - 1;
			startingPath[i].R = startingPath[i-1].R + 1;
			startingPath[i].O = 1;
		}
	}

}


/**************************/
/* Vanilla Uniformization */
/**************************/

void mcmc_newPathUnif(vector<mjpState>& startingPath, int sizeMJPpath, int population, double scaling, double infectionRate, double removalRate, mt19937& generator, vector<double>& ffArray)
{

	// Dominating rate

	double omega = scaling * ( infectionRate * ipow(population/2.0,2) + removalRate * population );
	double invOmega = 1.0/omega;

    // Fill the new states, only times for now

	for (int i = 0; i < sizeMJPpath-1; ++i)
	{

		// Compute virtual rate

		double virtualRate = omega - infectionRate * startingPath[i].S * startingPath[i].I - removalRate * startingPath[i].I;

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval

		for (int j = 0; j < dummy0; ++j)
		{
			startingPath.push_back({unif(generator) , 0, 0, 0, 2}) ; // the two is to keep track what we need to fill-up later
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](mjpState & a, mjpState & b) {return a.time < b.time;});

	// populate, add to stack some important information from the vector; speeds things up, apparently

	for (int i = 0; i < startingPath.size()-1; ++i)
	{
		if (startingPath[i+1].O ==2)
		{
			startingPath[i+1].S = startingPath[i].S;
			startingPath[i+1].I = startingPath[i].I;
			startingPath[i+1].R = startingPath[i].R;
			startingPath[i+1].O = 0;
		}
	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Resize array for FF

	ffArray.resize(startingPath.size() * (population+1));

	// Forward filtering steps

	for (int i = 1; i < startingPath.size(); ++i){

		if (startingPath[i].O != 1) {  // No removal is possible, only a infection or virtual

			for (int j = population; j > population - startingPath[i].R; --j) ffArray[(population+1)*i + j] = 0.0;
			for (int j = population - startingPath[i].R; j > 0; --j)
			{

				double virtualJumpProb = 1 - (infectionRate * (population - j - startingPath[i-1].R) + removalRate)* j * invOmega;

				// Transition with weighting

				ffArray[(population+1)*i + j] = ffArray[(population+1)*(i-1) + j-1] * infectionRate * (j-1.0) * (population - j + 1 - startingPath[i-1].R ) * invOmega + ffArray[(population+1)*(i-1) + j] * virtualJumpProb;

			}

			// Special case at j=0
			ffArray[(population+1)*i] = ffArray[(population+1)*(i-1)];

		} else { // If removal observation exists at this time step

			for (int j = 0; j < population - startingPath[i].R + 1; ++j)
			{

				// Transition with weighting

				ffArray[(population+1)*i + j] = ffArray[(population+1)*(i-1) + j + 1] * removalRate * (j+1.0) * invOmega;

			}
			for (int j = population - startingPath[i].R + 1; j < population + 1; ++j) ffArray[(population+1)*i + j] = 0.0;

		}

		// Normalise
		double invSum = 1.0 / accumulate (ffArray.begin() + (population+1)* i, ffArray.begin() + (population+1)*(i+1), 0.0) ;
		for (int j = 0; j <population+1; ++j) ffArray[(population+1)*i + j] *= invSum;

	}


	/*********************/
	/* Backward Sampling */
	/*********************/
        
    // Assume end point is known - the epidemic dies off
    
    startingPath[startingPath.size() - 1].I = 0;
	startingPath[startingPath.size() - 1].S = population - 0 - startingPath[startingPath.size() - 1].R;

	// Keep track of indexes to keep
	vector<int> toKeep; toKeep.reserve(sizeMJPpath);

    // Loop
    
    for (int i = startingPath.size() - 2; i >= 0; --i){
      
      if (startingPath[i+1].O == 1){ // If previous in backward step was a removal observation time...
        
        startingPath[i].I = startingPath[i+1].I +1;
        startingPath[i].S = population - startingPath[i].I - startingPath[i].R;
        toKeep.push_back(i+1);
 
      } else{ // either infection or virtual time
        
        int candidates[] = { startingPath[i+1].I - 1, startingPath[i+1].I };

        double beta[2] = {};
        double sum = 0.0;

        // Start with infection

        beta[0] = infectionRate * candidates[0] * i0max(population - candidates[0] - startingPath[i].R) * invOmega * ffArray[(population+1)*i + candidates[0]];
        sum = sum + beta[0];

        // Continue with virtual

 		double virtualJumpProb = 1 - (infectionRate * (population - candidates[1] - startingPath[i].R) + removalRate)* candidates[1] * invOmega;
        beta[1] = virtualJumpProb * ffArray[(population+1)*i + candidates[1]];
        sum = sum + beta[1];

        // Normalise and sample

        beta[0] = beta[0] / sum;
		uniform_real_distribution<double> unif(0.0, 1.0);

        if (unif(generator) <= beta[0]){
        	startingPath[i].I = candidates[0];
        	toKeep.push_back(i+1);
        } else {
        	startingPath[i].I = candidates[1];
        }
        startingPath[i].S = population - startingPath[i].I - startingPath[i].R;

      }
      
    }
    toKeep.push_back(0);

    // Remove virtual jumps and return; done from the back as these are cheaper operations
    for (int i = 0; i < sizeMJPpath; ++i) startingPath[i] = startingPath[toKeep[sizeMJPpath - i - 1]];
    startingPath.resize(sizeMJPpath);
}


/******************/
/* Double Poisson */
/******************/

void mcmc_newPathAuxVar(mjpState * startingPath,  vector<vector<double>>& odeApprox, int sizeMJPpath, int population, double scaling, double infectionRate, double removalRate, mt19937& generator, double * ffArray)
{

	/*****************/
	/* Tune Aux Vars */
	/*****************/
	
	int auxVarLag = 25; //25; // all results so far are on 25
	uniform_int_distribution<int> uniDist(1,auxVarLag);
 	int counterLag = uniDist(generator);

	double meanNormals = log(population/10); // log(population/10); // log(sqrt(population));
	//double autoregrVal = 1.0; //0.5;
	//double sdNormals = 0.25 ; //* ( 1-ipow(1 - autoregrVal, 2 * auxVarLag) ) / ( 1-ipow(1 - autoregrVal, 2) );
	int gammaAlpha = 2;// 2; all results so far are on 2

	// Dominating rate

	double omega = scaling * ( infectionRate * ipow(population/2.0,2) + removalRate * population );
	double invOmega = 1.0/omega;

    // Fill the new states

	int finalSize = sizeMJPpath;
	for (int i = 0; i < sizeMJPpath-1; ++i){

		// Get dom rate for virtuals based on "scaling" proportionality

		double virtualRate = (scaling - 1.0) * (infectionRate * startingPath[i].S * startingPath[i].I + removalRate * startingPath[i].I); 

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval

		// Times need sampling 
		for (int j = 0; j < dummy0; ++j) {

			startingPath[finalSize] = {unif(generator) , startingPath[i].S, startingPath[i].I, startingPath[i].R, 0};
			finalSize++;

		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath, startingPath + finalSize, [](mjpState & a, mjpState & b) {return a.time < b.time;});

	// Populate mapping to corresponding ODE value

	double ODEvalues[finalSize];
	double endTime = startingPath[finalSize-1].time;

	for (int i = 0; i < finalSize; ++i) ODEvalues[i] = odeApprox[i0max(startingPath[i].time / endTime * 1000.0)][1];

	// Add weights and auxiliary variables

	int weightsPoisson[finalSize];

	double invgammaBeta[finalSize];
	double auxVars[finalSize];
	int range[finalSize*2]; //normal_distribution<double> normal;

	// First instance

	double dummyMeanNormals = meanNormals;
	auxVars[0] = 0.0;
	range[0] = 1;
	range[finalSize] = 1;

	double virtualRate = omega - scaling * (infectionRate * startingPath[0].S * startingPath[0].I + removalRate * startingPath[0].I);
	poisson_distribution<int> pois((startingPath[1].time - startingPath[0].time) * virtualRate);
	weightsPoisson[0] = pois(generator);

	for (int i = 1; i < finalSize-1; ++i)
	{
		// Weights

		double virtualRate = omega - scaling * (infectionRate * startingPath[i].S * startingPath[i].I + removalRate * startingPath[i].I);
		pois = poisson_distribution<int> ((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		weightsPoisson[i] = pois(generator);

		// Auxiliary Variable

 		if (counterLag < auxVarLag){

			auxVars[i] = 0.0;
			invgammaBeta[i] = 0.0;
 			range[i] = i0max( range[i-1] - 1 );
 			range[finalSize + i] = min(population - startingPath[i].R, range[finalSize + i-1] + 1);

 			counterLag++;

 		} else {

 			//dummyMeanNormals = meanNormals + ( dummyMeanNormals - meanNormals ) * ipow(1 - autoregrVal, auxVarLag);
 			//normal = normal_distribution<double> (dummyMeanNormals,sdNormals); 
 			//dummyMeanNormals = normal(generator);

			invgammaBeta[i] = exp(dummyMeanNormals)/gammaAlpha;

			gamma_distribution<double> gamma(gammaAlpha,invgammaBeta[i]);
 			auxVars[i] =  abs(startingPath[i].I-ODEvalues[i]) + gamma(generator) ;

 			range[i] = i0max((int) ceil(ODEvalues[i]-auxVars[i]));
 			range[finalSize + i] = min(range[finalSize + i-1] + 1 , min(population - startingPath[i].R,(int) floor(ODEvalues[i]+auxVars[i])) );

 			counterLag = 1;

 		}

	}

	// Add finals

	weightsPoisson[finalSize-1] = 0;
	auxVars[finalSize-1] = 0.0;
	range[finalSize-1] = i0max( range[finalSize-2] - 1 );
	range[finalSize + finalSize-1] = min(population - startingPath[finalSize-1].R, range[finalSize + finalSize-2] + 1);


	/*********************/
	/* Forward Filtering */
	/*********************/

	for (int i = 1; i < finalSize; ++i){

		if (startingPath[i].O != 1) {  // No removal is possible, only a infection or virtual

			for (int j = range[finalSize+i]+2; j >range[finalSize+i]; --j) ffArray[(population+1)*i + j] = 0.0;

			for (int j = range[finalSize+i]; j >= max(1,range[i]); --j)
			{

				double virtualProb = (scaling - 1.0) * (infectionRate * (population - j - startingPath[i-1].R ) * j + removalRate * j); 
				double weighting = 1 - scaling * (infectionRate * (population - j - startingPath[i].R) + removalRate)* j * invOmega; 

				// Transition with weighting

				ffArray[(population+1)*i + j] = ipow(weighting,weightsPoisson[i]) * 
					(ffArray[(population+1)*(i-1) + j-1] * infectionRate * (j-1.0) * (population - j + 1 - startingPath[i-1].R ) + 
						ffArray[(population+1)*(i-1) + j] * virtualProb) ;

			}

			for (int j = i0max( range[i]-2 ); j < max(1,range[i]); ++j) ffArray[(population+1)*i + j] = 0.0;

		} else { // If removal observation exists at this time step

			for (int j = i0max( range[i]-2 ); j < range[i]; ++j) ffArray[(population+1)*i + j] = 0.0;

			for (int j = range[i]; j <= range[finalSize+i]; ++j)
			{

				double weighting = 1 - scaling * (infectionRate * (population - j - startingPath[i].R) + removalRate)* j * invOmega; 

				// Transition with weighting

				ffArray[(population+1)*i + j] = ipow(weighting,weightsPoisson[i]) * 
					ffArray[(population+1)*(i-1) + j + 1] * removalRate * (j+1.0);

			}

			for (int j = range[finalSize+i]+2; j >range[finalSize+i]; --j) ffArray[(population+1)*i + j] = 0.0;

		}

		// Auxiliary variable effect

		if (auxVars[i] != 0) {

			for (int j = range[i]; j <= range[finalSize+i]; ++j)
			{

				double gammaVar = auxVars[i] - abs(j-ODEvalues[i]);
				ffArray[(population+1)*i + j] *=  ipow(gammaVar,gammaAlpha - 1 )  *  fastExp(  - gammaVar / invgammaBeta[i] ) ;

			}

		}

		// Normalise

		double invSum = 1.0 / accumulate (ffArray + (population+1)* i + range[i], ffArray+(population+1)*i + range[finalSize+i]+1,0.0) ;
		for (int j = range[i]; j < range[finalSize+i]+1 ; ++j) ffArray[(population+1)*i + j] *= invSum;

	}


	/*********************/
	/* Backward Sampling */
	/*********************/
        
    // Assume end point is known - the epidemic dies off
    
    startingPath[finalSize - 1].I = 0;
	startingPath[finalSize - 1].S = population - 0 - startingPath[finalSize - 1].R;

	// Keep track of indexes to keep
	vector<int> toKeep; toKeep.reserve(sizeMJPpath);

    // Loop
    
    for (int i = finalSize - 2; i >= 0; --i){
      
      if (startingPath[i+1].O == 1){ // If previous in backward step was a removal observation time...
        
        startingPath[i].I = startingPath[i+1].I +1;
        startingPath[i].S = population - startingPath[i].I - startingPath[i].R;
        toKeep.push_back(i+1);
 
      } else{ // either infection or virtual time
        
        int candidates[] = { startingPath[i+1].I - 1, startingPath[i+1].I }; // check range issues!!!!!

        double beta[2] = {};
        double sum = 0.0;

        // Start with infection

        beta[0] = infectionRate * candidates[0] * i0max(population - candidates[0] - startingPath[i].R) * ffArray[(population+1)*i + candidates[0]];
	    sum = sum + beta[0];

        // Continue with virtual

       	beta[1] = (scaling - 1.0) * (infectionRate * i0max(population - candidates[1] - startingPath[i].R ) * candidates[1] + removalRate * candidates[1]) * ffArray[(population+1)*i + candidates[1]];
        sum = sum + beta[1];

        // Normalise and sample

        beta[0] = beta[0] / sum;
		uniform_real_distribution<double> unif(0.0, 1.0);

        if (unif(generator) <= beta[0]){
        	startingPath[i].I = candidates[0];
        	toKeep.push_back(i+1);
        } else {
        	startingPath[i].I = candidates[1];
        }
        startingPath[i].S = population - startingPath[i].I - startingPath[i].R;

      }
      
    }
    toKeep.push_back(0);

    // Remove virtual jumps and return; done from the back as these are cheaper operations
    for (int i = 0; i < sizeMJPpath; ++i) startingPath[i] = startingPath[toKeep[sizeMJPpath - i - 1]];

}


void mcmc_newPathDP(mjpState * startingPath, int sizeMJPpath, int population, double scaling, double infectionRate, double removalRate, mt19937& generator, double * ffArray)
{

	// Dominating rate

	double omega = scaling * ( infectionRate * ipow(population/2.0,2) + removalRate * population );
	double invOmega = 1.0/omega;

    // Fill the new states

	int finalSize = sizeMJPpath;
	for (int i = 0; i < sizeMJPpath-1; ++i){

		// Get dom rate for virtuals based on "scaling" proportionality

		double virtualRate = (scaling - 1.0) * (infectionRate * startingPath[i].S * startingPath[i].I + removalRate * startingPath[i].I); 

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval

		// Times need sampling 
		for (int j = 0; j < dummy0; ++j) {

			startingPath[finalSize] = {unif(generator) , startingPath[i].S, startingPath[i].I, startingPath[i].R, 0};
			finalSize++;

		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath, startingPath + finalSize, [](mjpState & a, mjpState & b) {return a.time < b.time;});

	// Add weights; add to stack some important information from the vector; speeds things up, apparently

	int weightsPoisson[finalSize];

	for (int i = 0; i < finalSize-1; ++i)
	{
		// Weights

		double virtualRate = omega - scaling * (infectionRate * startingPath[i].S * startingPath[i].I + removalRate * startingPath[i].I);
		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		weightsPoisson[i] = pois(generator);

	}

	// Add final weight

	weightsPoisson[finalSize-1] = 0;

	/*********************/
	/* Forward Filtering */
	/*********************/

	for (int i = 1; i < finalSize; ++i){

		if (startingPath[i].O != 1) {  // No removal is possible, only a infection or virtual

			for (int j = population; j > population - startingPath[i].R; --j) ffArray[(population+1)*i + j] = 0.0;
			for (int j = population - startingPath[i].R; j > 0; --j)
			{

				double virtualProb = (scaling - 1.0) * (infectionRate * (population - j - startingPath[i-1].R ) * j + removalRate * j); 
				double weighting = 1 - scaling * (infectionRate * (population - j - startingPath[i].R) + removalRate)* j * invOmega; 

				// Transition with weighting

				ffArray[(population+1)*i + j] = ipow(weighting,weightsPoisson[i]) * 
					(ffArray[(population+1)*(i-1) + j-1] * infectionRate * (j-1.0) * (population - j + 1 - startingPath[i-1].R ) + 
						ffArray[(population+1)*(i-1) + j] * virtualProb) ;

			}

			// Special case at j=0
			ffArray[(population+1)*i] = 0;

		} else { // If removal observation exists at this time step

			for (int j = 0; j < population - startingPath[i].R + 1; ++j)
			{

				double weighting = 1 - scaling * (infectionRate * (population - j - startingPath[i].R) + removalRate)* j * invOmega; 

				// Transition with weighting

				ffArray[(population+1)*i + j] = ipow(weighting,weightsPoisson[i]) * 
					ffArray[(population+1)*(i-1) + j + 1] * removalRate * (j+1.0);

			}
			for (int j = population - startingPath[i].R + 1; j < population + 1; ++j) ffArray[(population+1)*i + j] = 0.0;

		}

		// Normalise
		double invSum = 1.0 / accumulate (ffArray + (population+1)* i, ffArray+(population+1)*(i+1),0.0) ;
		for (int j = 0; j <population+1; ++j) ffArray[(population+1)*i + j] *= invSum;

	}


	// /*********************/
	// /* Backward Sampling */
	// *******************
        
    // Assume end point is known - the epidemic dies off
    
    startingPath[finalSize - 1].I = 0;
	startingPath[finalSize - 1].S = population - 0 - startingPath[finalSize - 1].R;

	// Keep track of indexes to keep
	vector<int> toKeep; toKeep.reserve(sizeMJPpath);

    // Loop
    
    for (int i = finalSize - 2; i >= 0; --i){
      
      if (startingPath[i+1].O == 1){ // If previous in backward step was a removal observation time...
        
        startingPath[i].I = startingPath[i+1].I +1;
        startingPath[i].S = population - startingPath[i].I - startingPath[i].R;
        toKeep.push_back(i+1);
 
      } else{ // either infection or virtual time
        
        int candidates[] = { startingPath[i+1].I - 1, startingPath[i+1].I };

        double beta[2] = {};
        double sum = 0.0;

        // Start with infection

        beta[0] = infectionRate * candidates[0] * i0max(population - candidates[0] - startingPath[i].R) * ffArray[(population+1)*i + candidates[0]];
        sum = sum + beta[0];

        // Continue with virtual

        beta[1] = (scaling - 1.0) * (infectionRate * i0max(population - candidates[1] - startingPath[i].R ) * candidates[1] + removalRate * candidates[1]) * ffArray[(population+1)*i + candidates[1]];
        sum = sum + beta[1];

        // Normalise and sample

        beta[0] = beta[0] / sum;
		uniform_real_distribution<double> unif(0.0, 1.0);

        if (unif(generator) <= beta[0]){
        	startingPath[i].I = candidates[0];
        	toKeep.push_back(i+1);
        } else {
        	startingPath[i].I = candidates[1];
        }
        startingPath[i].S = population - startingPath[i].I - startingPath[i].R;

      }
      
    }
    toKeep.push_back(0);

    // Remove virtual jumps and return; done from the back as these are cheaper operations
    for (int i = 0; i < sizeMJPpath; ++i) startingPath[i] = startingPath[toKeep[sizeMJPpath - i - 1]];
}


/*******************/
/* MCMC Procedures */
/*******************/

void do_mcmc_Unif(vector<double> removalData, int population, int initials, int mcmcIter, double scaling, mt19937 generator, string filename, int chainNumber)
{
	// Create starting trajectory

	const int sizeMJPpath = (int) 2*removalData.size();
	vector<mjpState> startingPath; 
	mcmc_startingPathVector(startingPath, removalData, population, initials, generator);

	// Priors and sampler for rates

	double priorShapeRemoval = 1.0, priorRateRemoval = pow(10.0,-3.0), priorShapeInfec = 1.0, priorRateInfec = pow(10.0,-3.0);

	// Build a forward probabilities matrix to reuse, reserve enough space and then pass the pointer to functions

	vector<double> ffArray(population + 1); for (int i = 0; i < population + 1; ++i) ffArray[i] = 0.0; ffArray[1] = 1.0;

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "beta, gamma, t0 \n";

	// ITERATE CHOICE OF ALGORITHM - Assign below

	for (int i = 0; i < mcmcIter; ++i)
	{

		// Rates update

		double shapeRemove = startingPath[sizeMJPpath-1].R + priorShapeRemoval;
		double shapeInfec = startingPath[0].S - startingPath[sizeMJPpath-1].S + priorShapeInfec;
		double rateRemove = priorRateRemoval;
		double rateInfec = priorRateInfec;
		for (int j = 0; j < sizeMJPpath-1; ++j)
		{
			rateRemove = rateRemove + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].I;
			rateInfec = rateInfec + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].I * startingPath[j].S;
		}
		gamma_distribution<double> gamma(shapeRemove,1.0/rateRemove);
		double removalRate = gamma(generator);
		gamma = gamma_distribution<double>(shapeInfec,1.0/rateInfec);
		double infectionRate = gamma(generator);

		// First infection time update

		exponential_distribution<double> expon( startingPath[0].S * startingPath[0].I * infectionRate + startingPath[0].I * removalRate);
		startingPath[0].time = startingPath[1].time - expon(generator);

      	// Save parameters and update trajectory

		myfileTrace << infectionRate << ", " << removalRate << ", " <<  startingPath[0].time <<"\n";

        // Update trajectory

		mcmc_newPathUnif(startingPath, sizeMJPpath, population, scaling, infectionRate, removalRate, generator, ffArray);

		// Progress bar

		if (mcmcIter > 49 && (i+1) % (mcmcIter/50) == 0 && chainNumber==1){

			float progress = (float) (i+1) / (float) mcmcIter;
			int barWidth = 50;
			cout << "\r[";
			int pos = (int) barWidth * progress;
			for (int j = 0; j < barWidth; ++j) {
				if (j <= pos) cout << "*";
				else cout << " ";
			}
			cout << "] " << int(progress * 100.0) << " %";
			cout.flush();
		}

	}
	myfileTrace.close();

	// // IF WANT TO OUTPUT THAT PATH TO HAVE A LOOK IN R
	// ofstream myfile;
	// myfile.open ("startingPath.csv");
	// myfile << "Time, Susceptible, Infected, Removed, Observation \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].S << ", " << startingPath[i].I << ", " << startingPath[i].R << ", " << startingPath[i].O << "\n";
	// }
	// myfile.close();

}

void do_mcmc_DP(vector<double> removalData, int population, int initials, int mcmcIter, double scaling, mt19937 generator, string filename, int chainNumber)
{
	// Create starting trajectory

	const int sizeMJPpath = (int) 2*removalData.size();
	mjpState * startingPath = new mjpState [(int) 3*2*removalData.size()]; // put more space to make sure the virtual jumps fit
	mcmc_startingPath(startingPath, removalData, population, initials, generator);

	// Priors and sampler for rates

	double priorShapeRemoval = 1.0, priorRateRemoval = pow(10.0,-3.0), priorShapeInfec = 1.0, priorRateInfec = pow(10.0,-3.0);

	// Build a forward probabilities matrix to reuse, reserve enough space and then pass the pointer to functions

	double * ffArray = new double [ (int) (3*2*removalData.size() * (population+1)) ];
	for (int i = 0; i < population+1; ++i) ffArray[i] = 0.0; // Initialize all to 0, as we will later avoid accessing every index always
	ffArray[1] = 1.0;

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "beta, gamma, t0 \n";

	// ITERATE CHOICE OF ALGORITHM - Assign below

	for (int i = 0; i < mcmcIter; ++i)
	{

		// Rates update

		double shapeRemove = startingPath[sizeMJPpath-1].R + priorShapeRemoval;
		double shapeInfec = startingPath[0].S - startingPath[sizeMJPpath-1].S + priorShapeInfec;
		double rateRemove = priorRateRemoval;
		double rateInfec = priorRateInfec;
		for (int j = 0; j < sizeMJPpath-1; ++j)
		{
			rateRemove = rateRemove + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].I;
			rateInfec = rateInfec + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].I * startingPath[j].S;
		}
		gamma_distribution<double> gamma(shapeRemove,1.0/rateRemove);
		double removalRate = gamma(generator);
		gamma = gamma_distribution<double>(shapeInfec,1.0/rateInfec);
		double infectionRate = gamma(generator);

		// First infection time update

		exponential_distribution<double> expon( startingPath[0].S * startingPath[0].I * infectionRate + startingPath[0].I * removalRate);
		startingPath[0].time = startingPath[1].time - expon(generator);

      	// Save parameters and update trajectory

		myfileTrace << infectionRate << ", " << removalRate << ", " <<  startingPath[0].time <<"\n";

		
        // Update trajectory

		mcmc_newPathDP(startingPath, sizeMJPpath, population, scaling, infectionRate, removalRate, generator, ffArray);

		// Progress bar

		if (mcmcIter > 49 && (i+1) % (mcmcIter/50) == 0 && chainNumber==1){

			float progress = (float) (i+1) / (float) mcmcIter;
			int barWidth = 50;
			cout << "\r[";
			int pos = (int) barWidth * progress;
			for (int j = 0; j < barWidth; ++j) {
				if (j <= pos) cout << "*";
				else cout << " ";
			}
			cout << "] " << int(progress * 100.0) << " %";
			cout.flush();
		}

	}
	myfileTrace.close();

	// // IF WANT TO OUTPUT THAT PATH TO HAVE A LOOK IN R
	// ofstream myfile;
	// myfile.open ("startingPath.csv");
	// myfile << "Time, Susceptible, Infected, Removed, Observation \n";
	// for (int i = 0; i < sizeMJPpath; ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].S << ", " << startingPath[i].I << ", " << startingPath[i].R << ", " << startingPath[i].O << "\n";
	// }
	// myfile.close();

}

void do_mcmc_AuxVar(vector<double> removalData, vector<vector<double>> odeApprox, int population, int initials, int mcmcIter, double scaling, mt19937 generator, string filename, int chainNumber)
{
	// Create starting trajectory

	const int sizeMJPpath = (int) 2*removalData.size();
	mjpState * startingPath = new mjpState [(int) 3*2*removalData.size()]; // put more space to make sure the virtual jumps fit
	mcmc_startingPath(startingPath, removalData, population, initials, generator);

	// Priors and sampler for rates

	double priorShapeRemoval = 1.0, priorRateRemoval = pow(10.0,-3.0), priorShapeInfec = 1.0, priorRateInfec = pow(10.0,-3.0);

	// Build a forward probabilities matrix to reuse, reserve enough space and then pass the pointer to functions

	double * ffArray = new double [ (int) (3*2*removalData.size() * (population+1)) ];
	for (int i = 0; i < population+1; ++i) ffArray[i] = 0.0; // Initialize all to 0, as we will later avoid accessing every index always
	ffArray[1] = 1.0;

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "beta, gamma, t0 \n";

	// ITERATE CHOICE OF ALGORITHM - Assign below

	for (int i = 0; i < mcmcIter; ++i)
	{

		// Rates update

		double shapeRemove = startingPath[sizeMJPpath-1].R + priorShapeRemoval;
		double shapeInfec = startingPath[0].S - startingPath[sizeMJPpath-1].S + priorShapeInfec;
		double rateRemove = priorRateRemoval;
		double rateInfec = priorRateInfec;
		for (int j = 0; j < sizeMJPpath-1; ++j)
		{
			rateRemove = rateRemove + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].I;
			rateInfec = rateInfec + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].I * startingPath[j].S;
		}
		gamma_distribution<double> gamma(shapeRemove,1.0/rateRemove);
		double removalRate = gamma(generator);
		gamma = gamma_distribution<double>(shapeInfec,1.0/rateInfec);
		double infectionRate = gamma(generator);

		// First infection time update

		exponential_distribution<double> expon( startingPath[0].S * startingPath[0].I * infectionRate + startingPath[0].I * removalRate);
		startingPath[0].time = startingPath[1].time - expon(generator);

      	// Save parameters and update trajectory

		myfileTrace << infectionRate << ", " << removalRate << ", " <<  startingPath[0].time <<"\n";

		
        // Update trajectory

		mcmc_newPathAuxVar(startingPath, odeApprox, sizeMJPpath, population, scaling, infectionRate, removalRate, generator, ffArray);

		// Progress bar

		if (mcmcIter > 49 && (i+1) % (mcmcIter/50) == 0 && chainNumber==1){

			float progress = (float) (i+1) / (float) mcmcIter;
			int barWidth = 50;
			cout << "\r[";
			int pos = (int) barWidth * progress;
			for (int j = 0; j < barWidth; ++j) {
				if (j <= pos) cout << "*";
				else cout << " ";
			}
			cout << "] " << int(progress * 100.0) << " %";
			cout.flush();
		}

	}
	myfileTrace.close();

	// IF WANT TO OUTPUT THAT PATH TO HAVE A LOOK IN R
	// ofstream myfile;
	// myfile.open ("startingPath.csv");
	// myfile << "Time, Susceptible, Infected, Removed, Observation \n";
	// for (int i = 0; i < sizeMJPpath; ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].S << ", " << startingPath[i].I << ", " << startingPath[i].R << ", " << startingPath[i].O << "\n";
	// }
	// myfile.close();

}