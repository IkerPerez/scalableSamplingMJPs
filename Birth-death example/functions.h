
/*****************************/
/* Stuff to read a text file */
/*****************************/

vector<vector<double>> fileToVector(string filename)
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

// Seasonality functions

double seasonality(double t) 
{
    return 1.5 + cos(2.0*M_PI*t/100.0)/2.0;
}

double seasonalityIntegrated(double t) 
{
    return 1.5*t + sin(2.0*M_PI*t/100.0)/2.0 * 100.0/2.0/M_PI;
}

double normalCDF(double value)
{
   return 0.5 * erfc(-value * M_SQRT1_2);
}

// Look-up tables for exponentials over range -0:1000; note rates always negative

double expoFloors[1000001];
#define fastExp(x) ((x > -1000) ? expoFloors[ (int) (-(x)*1000.0) ] : 0) // enough quality approximtion here


/******************************/
/* To produce a starting Path */
/******************************/

void mcmc_startingPathVector(vector<BDstate>& startingPath, vector<vector<double>>& observations, int N, mt19937& gen)
{
	// Create and fill a frame to operate within

	double sum = 0;
	for (int i = 0; i < observations.size(); ++i) sum += observations[i][1];
	int avgState = (int) sum / observations.size();
	startingPath.push_back({0.0,avgState});
	double jumpswithin =  max(1.0,N/100.0);
	for (int i = 1; i < observations.size(); ++i) {
		for (int j = 0; j < jumpswithin; ++j)
		{
			startingPath.push_back({observations[i][0] + ( - observations[i][0] + observations[i-1][0]) * (j+1)/(jumpswithin+1),avgState});
		}
	}

}

void mcmc_startingPathVector2(vector<BDstate>& startingPath, vector<vector<double>>& observations, vector<vector<double>> odeApprox, int N, mt19937& gen)
{
	// Create and fill a frame to operate within

	double sum = 0;
	for (int i = 0; i < observations.size(); ++i) sum += observations[i][1];
	int avgState = (int) sum / observations.size();
	startingPath.push_back({0.0, (int) odeApprox[0][1]});
	double jumpswithin =  max(1.0,N/100.0);
	for (int i = 1; i < observations.size(); ++i) {
		for (int j = 0; j < jumpswithin; ++j)
		{
			startingPath.push_back({observations[i][0] + ( - observations[i][0] + observations[i-1][0]) * (j+1)/(jumpswithin+1),
				(int) odeApprox[1 + 10 * (observations[i][0] + ( - observations[i][0] + observations[i-1][0]) * (j+1)/(jumpswithin+1))][1]});
		}
	}

}

/**************************/
/* Vanilla Uniformization */
/**************************/

void mcmc_newPathUnif(vector<BDstate>& startingPath, vector<vector<double>>& observations, int N, double obsSd, double birthRate, double scaling, double deathRate, mt19937& generator)
{

	// Dominating rate

	double omega = scaling * (birthRate + deathRate * 2 * N);
	double invOmega = 1.0/omega;

    // Fill the new states

	startingPath.push_back({observations.back()[0] , -1}); // to make sure this goes on until last observation time
	int sizePathOriginal = startingPath.size();
	for (int i = 0; i < sizePathOriginal-1; ++i)
	{

		// Compute dominating virtual rate

		double virtualRate = omega - birthRate * (startingPath[i].X < N)  - deathRate * startingPath[i].X; // then thin seasonality

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);
		uniform_real_distribution<double> unifProb(0,1);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval, before thinning

		// Times need sampling and thinning, since non-homogeneous process!
		for (int j = 0; j < dummy0; ++j)
		{

			double candidateTime = unif(generator);
			double candidateProb = (omega - birthRate * (startingPath[i].X < N)  - deathRate * startingPath[i].X * seasonality(candidateTime) ) / virtualRate;
			if (unifProb(generator)<candidateProb) startingPath.push_back({candidateTime , startingPath[i].X}) ; 
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](BDstate & a, BDstate & b) {return a.time < b.time;});

	// Remove last entry, which I added to ensure we reach the end of observation time

	startingPath.pop_back() ;

	// Account for observations, at the latest jump right before each observation time, there could be more than one!

	vector<vector<double>> vectorObs(startingPath.size());
	int counter = 0;
	for (int i = 0; i < observations.size(); ++i)
	{
		while(startingPath[counter].time < observations[i][0] && counter < startingPath.size()) counter++;
		vectorObs[counter-1].push_back(observations[i][1]);
	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Build array for FF

	vector<double> ffArray(startingPath.size() * (N+1), 0);

	// First iteration

	for (int i = 0; i < N + 1; ++i) ffArray[i] = 1.0/( (double) N+1);

	// Weight by observation likelihood; if any

	for (int k = 0; k < vectorObs[0].size(); ++k) {

		for (int j = 0; j < N+1; ++j) ffArray[j] *=  exp( - ipow((vectorObs[0][k]-j)/obsSd,2) / 2.0 );

	}

	// Normalise
	double invSum = 1.0 / accumulate (ffArray.begin(), ffArray.begin() + (N+1), 0.0) ;
	for (int j = 0; j <N+1; ++j) ffArray[j] *= invSum;

	// Forward filtering steps

	for (int i = 1; i < startingPath.size(); ++i){

		for (int j = 1; j < N; ++j)
		{

			double virtualJumpProb = 1 - (birthRate + deathRate * j * seasonality(startingPath[i].time)) * invOmega;

			// Transition 

			ffArray[(N+1)*i + j] = ffArray[(N+1)*(i-1) + j-1] * birthRate * invOmega + ffArray[(N+1)*(i-1) + j + 1] * deathRate * (j+1) * seasonality(startingPath[i].time) * invOmega + ffArray[(N+1)*(i-1) + j] * virtualJumpProb;

		}

		// Special case at j=0

		double virtualJumpProb = 1 - birthRate * invOmega;
		ffArray[(N+1)*i] = ffArray[(N+1)*(i-1) + 1] * deathRate * seasonality(startingPath[i].time) * invOmega + ffArray[(N+1)*(i-1)] * virtualJumpProb;

		// Special case at j=N

		virtualJumpProb = 1 - deathRate * N * seasonality(startingPath[i].time) * invOmega;
		ffArray[(N+1)*i + N] = ffArray[(N+1)*(i-1) + N - 1] * birthRate * invOmega + ffArray[(N+1)*(i-1) + N] * virtualJumpProb;

		
		// Weight entire row by observation likelihood; if any
		for (int k = 0; k < vectorObs[i].size(); ++k) {

			for (int j = 0; j < N+1; ++j) ffArray[(N+1)*i + j] *=  exp( - ipow((vectorObs[i][k]-j)/obsSd,2) / 2.0 );

		}

		// Normalise
		double invSum = 1.0 / accumulate (ffArray.begin() + (N+1)* i, ffArray.begin() + (N+1)*(i+1), 0.0) ;
		for (int j = 0; j <N+1; ++j) ffArray[(N+1)*i + j] *= invSum;

	}


	/*********************/
	/* Backward Sampling */
	/*********************/
        
	// Sample end value from distribution at end point
	    
	uniform_real_distribution<double> unif(0.0, 1.0);
	double dummyUnif = unif(generator);
	int dummyIndex = -1;
	double dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += ffArray[(N+1)*(startingPath.size()-1) + dummyIndex];
   	}
    startingPath[startingPath.size()-1].X = dummyIndex;
	
    // Loop
    
    for (int i = startingPath.size() - 2; i >= 0; --i){
              
        int candidates[] = { startingPath[i+1].X - 1, startingPath[i+1].X , startingPath[i+1].X + 1};

        double beta[3] = {0.0,0.0,0.0};
        double sum = 0.0;

        // Start with arrival
        if( startingPath[i+1].X > 0){
        	
        	beta[0] = birthRate * invOmega * ffArray[(N+1)*i + candidates[0]];
        	sum = sum + beta[0];

    	}

        // Continue with virtual

 		double virtualJumpProb = 1 - (birthRate * (candidates[1] < N ? 1:0) + deathRate * candidates[1] * seasonality(startingPath[i+1].time)) * invOmega;
        beta[1] = virtualJumpProb * ffArray[(N+1)*i + candidates[1]];
        sum = sum + beta[1];

        // Finally death
        if( startingPath[i+1].X < N){
        	
        	beta[2] =  deathRate * candidates[2] * seasonality(startingPath[i+1].time) * invOmega * ffArray[(N+1)*i + candidates[2]];
        	sum = sum + beta[2];

    	}

        // Normalise and sample

        for (int j = 0; j < 2; ++j) beta[j] = beta[j] / sum;

		dummyUnif = unif(generator);

        if (dummyUnif <= beta[0]){
        	startingPath[i].X = candidates[0];
        } else if (dummyUnif <= beta[0]+beta[1]){
        	startingPath[i].X = candidates[1];
        } else {
        	startingPath[i].X = candidates[2];
        }
      
    }

    // Remove virtual jumps and return
    counter = 1;
    while (counter < startingPath.size()) {
    
    	if(startingPath[counter].X == startingPath[counter-1].X) {
    		startingPath.erase(startingPath.begin() + counter);
    	} else{
    		counter++;
    	}

    }

}


/******************/
/* Depen Thinning */
/******************/

void mcmc_newPathDepThinning(vector<BDstate>& startingPath, vector<vector<double>>& observations, int N, double obsSd, double birthRate, double K, double deathRate, mt19937& generator)
{

    // Fill the new states

	startingPath.push_back({observations.back()[0] , -1}); // to make sure this goes on until last observation time
	int sizePathOriginal = startingPath.size();
	for (int i = 0; i < sizePathOriginal-1; ++i)
	{

		// Get dom rate for virtuals based on K proportionality; then we will thin to get non-homogeneous behwaviour

		double virtualRate = K * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X * 2); // then thin seasonality

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);
		uniform_real_distribution<double> unifProb(0,1);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval, before thinning

		// Times need sampling and thinning, since non-homogeneous process!
		for (int j = 0; j < dummy0; ++j)
		{
			double candidateTime = unif(generator);
			double candidateProb = K * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X * seasonality(candidateTime)) /virtualRate;
			if (unifProb(generator)<candidateProb) startingPath.push_back({candidateTime , startingPath[i].X}) ; 
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](BDstate & a, BDstate & b) {return a.time < b.time;});

	// Remove last entry, which I added to ensure we reach the end of observation time

	startingPath.pop_back() ;
	
	// Store final size of augmented frame
	int finalSize = startingPath.size();

	// populate time differences

	double timeDiffs[finalSize];
	double timeDiffsSeas[finalSize];

	for (int i = 0; i < finalSize-1; ++i)
	{

		timeDiffs[i] = startingPath[i+1].time - startingPath[i].time;
		timeDiffsSeas[i] = seasonalityIntegrated(startingPath[i+1].time) - seasonalityIntegrated(startingPath[i].time);

	}

	// Add final

	timeDiffs[finalSize-1] = observations.back()[0] - startingPath[finalSize-1].time;
	timeDiffsSeas[finalSize-1] = seasonalityIntegrated(observations.back()[0]) - seasonalityIntegrated(startingPath[finalSize-1].time);

	// Account for observations, at the latest jump right before each observation time

	vector<vector<double>> vectorObs(startingPath.size());
	int counter = 0;
	for (int i = 0; i < observations.size(); ++i)
	{
		while(startingPath[counter].time < observations[i][0] && counter < startingPath.size()) counter++;
		vectorObs[counter-1].push_back(observations[i][1]);
	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Stuff for weights and transitions; to speed things up

	double weighting1[N+1];
	double weighting2[N+1];
	for (int i = 0; i < N; ++i)	{
		weighting1[i] = - birthRate; 
		weighting2[i] = - deathRate * i; 
	}
	weighting1[N] = 0; 
	weighting2[N] = - deathRate * N; 

	// Build array for FF

	vector<double> ffArray(finalSize * (N+1), 0);

	// First iteration

	for (int i = 0; i < N + 1; ++i) ffArray[i] = fastExp( (1+K) * (weighting1[i]*timeDiffs[0] + weighting2[i]*timeDiffsSeas[0]) );

	// Weight by observation likelihood; if any
	for (int k = 0; k < vectorObs[0].size(); ++k) {

		// till upper bound only
		for (int j = 0; j < N+1; ++j) ffArray[j] *=  fastExp( - ipow((vectorObs[0][k]-j)/obsSd,2) / 2.0 );

	}

	// Normalise
	double invSum = 1.0 / accumulate (ffArray.begin(), ffArray.begin() + (N+1), 0.0) ;
	for (int j = 0; j <N+1; ++j) ffArray[j] *= invSum;

	// Forward filtering iterations

	for (int i = 1; i < finalSize; ++i){

		for (int j = 1; j < N; ++j)
		{

			double virtualProb = K * ( birthRate + deathRate * j * seasonality(startingPath[i].time) );

			// Transition with weighting

			ffArray[(N+1)*i + j] = fastExp( (1+K) * (weighting1[j]*timeDiffs[i] + weighting2[j]*timeDiffsSeas[i]) ) *  (ffArray[(N+1)*(i-1) + j-1] * birthRate + ffArray[(N+1)*(i-1) + j + 1] * deathRate * (j+1) * seasonality(startingPath[i].time) + ffArray[(N+1)*(i-1) + j] * virtualProb);

		}

		// Special case at j=0

		ffArray[(N+1)*i] = fastExp( (1+K) * (weighting1[0]*timeDiffs[i] + weighting2[0]*timeDiffsSeas[i]) ) * (ffArray[(N+1)*(i-1) + 1] * deathRate * seasonality(startingPath[i].time) + ffArray[(N+1)*(i-1)] * K * birthRate);

		// Special case at j=N

		ffArray[(N+1)*i + N] = fastExp( (1+K) * (weighting1[N]*timeDiffs[i] + weighting2[N]*timeDiffsSeas[i]) ) * (ffArray[(N+1)*(i-1) + N - 1] * birthRate + ffArray[(N+1)*(i-1) + N] * K * deathRate * N * seasonality(startingPath[i].time) );

		// Weight entire row by observation likelihood; if any
		for (int k = 0; k < vectorObs[i].size(); ++k) {

			for (int j = 0; j < N+1; ++j) ffArray[(N+1)*i + j] *=  fastExp( - ipow((vectorObs[i][k]-j)/obsSd,2) / 2.0 );

		}

		// Normalise
		double invSum = 1.0 / accumulate (ffArray.begin() + (N+1)* i, ffArray.begin() + (N+1)*(i+1), 0.0) ;
		for (int j = 0; j <N+1; ++j) ffArray[(N+1)*i + j] *= invSum;

	}

	/*********************/
	/* Backward Sampling */
	/*********************/
        
	// Sample end value from distribution at end point
	    
	uniform_real_distribution<double> unif(0.0, 1.0);
	double dummyUnif = unif(generator);
	int dummyIndex = -1;
	double dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += ffArray[(N+1)*(finalSize-1) + dummyIndex];
   	}
    startingPath[finalSize-1].X = dummyIndex;

    // Loop
    
    for (int i = finalSize - 2; i >= 0; --i){

    	int candidates[] = { startingPath[i+1].X - 1, startingPath[i+1].X , startingPath[i+1].X + 1};

    	double beta[3] = {0.0,0.0,0.0};
    	double sum = 0.0;

        // Start with arrival
    	if( startingPath[i+1].X > 0){

    		beta[0] = birthRate * ffArray[(N+1)*i + candidates[0]];
    		sum = sum + beta[0];

    	}

        // Continue with virtual

    	beta[1] = K * ( birthRate * (startingPath[i+1].X < N) + deathRate * startingPath[i+1].X * seasonality(startingPath[i+1].time) ) * ffArray[(N+1)*i + candidates[1]];
    	sum = sum + beta[1];

        // Finally death
        if( startingPath[i+1].X < N){
        	
        	beta[2] =  deathRate * candidates[2] * seasonality(startingPath[i+1].time) * ffArray[(N+1)*i + candidates[2]];
        	sum = sum + beta[2];

    	}

        // Normalise and sample

    	for (int j = 0; j < 2; ++j) beta[j] = beta[j] / sum;

    	dummyUnif = unif(generator);

    	if (dummyUnif <= beta[0]){
    		startingPath[i].X = candidates[0];
    	} else if (dummyUnif <= beta[0]+beta[1]){
    		startingPath[i].X = candidates[1];
    	} else {
    		startingPath[i].X = candidates[2];
    	}

    }

    // Remove virtual jumps and return
    counter = 1;
    while (counter < startingPath.size()) {
    
    	if(startingPath[counter].X == startingPath[counter-1].X) {
    		startingPath.erase(startingPath.begin() + counter);
    	} else{
    		counter++;
    	}

    }

}

void mcmc_newPathDepThinningCounts(vector<BDstate>& startingPath, vector<vector<double>>& observations, int N, double obsSd, double birthRate, double K, double deathRate, mt19937& generator)
{

	// Global dominating rate

	double omega = (1+K) * (birthRate + deathRate * 2 * N);

    // Fill the new states

	startingPath.push_back({observations.back()[0] , -1}); // to make sure this goes on until last observation time
	int sizePathOriginal = startingPath.size();
	for (int i = 0; i < sizePathOriginal-1; ++i){

		// Get dom rate for virtuals based on K proportionality; then we will thin to get non-homogeneous behwaviour

		double virtualRate = K * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X * 2); // then thin seasonality

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);
		uniform_real_distribution<double> unifProb(0,1);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval, before thinning

		// Times need sampling and thinning, since non-homogeneous process!
		for (int j = 0; j < dummy0; ++j)
		{
			double candidateTime = unif(generator);
			double candidateProb = K * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X * seasonality(candidateTime)) /virtualRate;
			if (unifProb(generator)<candidateProb) startingPath.push_back({candidateTime , startingPath[i].X}) ; 
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](BDstate & a, BDstate & b) {return a.time < b.time;});

	// Now compute vectors of weights of length (N+1) up until the previous to the last epoch
	
	vector<double> weights((startingPath.size()-1) * (N+1), 1.0);

	for (int i = 0; i < (startingPath.size()-1); ++i){

		// Get dom rate

		double virtualRate = omega - (1+K) * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X); // then thin seasonality

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);
		uniform_real_distribution<double> unifProb(0,1);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval, before thinning

		// now accept times and if so, recalculate weights
		for (int j = 0; j < dummy0; ++j)
		{
			double candidateTime = unif(generator);
			double candidateProb = (omega - (1+K) * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X * seasonality(candidateTime))) /virtualRate;
			if (unifProb(generator)<candidateProb) {

				for (int k = 0; k < N+1; ++k)
				{
					weights[(N+1)*i + k] *= (1 - (1+K) * (birthRate * (k < N)  + deathRate * k * seasonality(candidateTime)) /omega);
				}

			}
		}
	
	}

	// Remove last entry, which I added to ensure we reach the end of observation time

	startingPath.pop_back() ;
	
	// Store final size of augmented frame
	int finalSize = startingPath.size();

	// Account for observations, at the latest jump right before each observation time

	vector<vector<double>> vectorObs(startingPath.size());
	int counter = 0;
	for (int i = 0; i < observations.size(); ++i)
	{
		while(startingPath[counter].time < observations[i][0] && counter < startingPath.size()) counter++;
		vectorObs[counter-1].push_back(observations[i][1]);
	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Build array for FF

	vector<double> ffArray(finalSize * (N+1), 0);

	// First iteration

	for (int i = 0; i < N + 1; ++i) ffArray[i] = weights[i];

	// Weight by observation likelihood; if any
	for (int k = 0; k < vectorObs[0].size(); ++k) {

		// till upper bound only
		for (int j = 0; j < N+1; ++j) ffArray[j] *=  fastExp( - ipow((vectorObs[0][k]-j)/obsSd,2) / 2.0 );

	}

	// Normalise
	double invSum = 1.0 / accumulate (ffArray.begin(), ffArray.begin() + (N+1), 0.0) ;
	for (int j = 0; j <N+1; ++j) ffArray[j] *= invSum;

	// Forward filtering iterations

	for (int i = 1; i < finalSize; ++i){

		for (int j = 1; j < N; ++j)
		{

			double virtualProb = K * ( birthRate + deathRate * j * seasonality(startingPath[i].time) );

			// Transition with weighting

			ffArray[(N+1)*i + j] = weights[(N+1)*i + j] *  (ffArray[(N+1)*(i-1) + j-1] * birthRate + ffArray[(N+1)*(i-1) + j + 1] * deathRate * (j+1) * seasonality(startingPath[i].time) + ffArray[(N+1)*(i-1) + j] * virtualProb);

		}

		// Special case at j=0

		ffArray[(N+1)*i] = weights[(N+1)*i] * (ffArray[(N+1)*(i-1) + 1] * deathRate * seasonality(startingPath[i].time) + ffArray[(N+1)*(i-1)] * K * birthRate);

		// Special case at j=N

		ffArray[(N+1)*i + N] = weights[(N+1)*i + N] * (ffArray[(N+1)*(i-1) + N - 1] * birthRate + ffArray[(N+1)*(i-1) + N] * K * deathRate * N * seasonality(startingPath[i].time) );

		// Weight entire row by observation likelihood; if any
		for (int k = 0; k < vectorObs[i].size(); ++k) {

			for (int j = 0; j < N+1; ++j) ffArray[(N+1)*i + j] *=  fastExp( - ipow((vectorObs[i][k]-j)/obsSd,2) / 2.0 );

		}

		// Normalise
		double invSum = 1.0 / accumulate (ffArray.begin() + (N+1)* i, ffArray.begin() + (N+1)*(i+1), 0.0) ;
		for (int j = 0; j <N+1; ++j) ffArray[(N+1)*i + j] *= invSum;

	}

	/*********************/
	/* Backward Sampling */
	/*********************/
        
	// Sample end value from distribution at end point
	    
	uniform_real_distribution<double> unif(0.0, 1.0);
	double dummyUnif = unif(generator);
	int dummyIndex = -1;
	double dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += ffArray[(N+1)*(finalSize-1) + dummyIndex];
   	}
    startingPath[finalSize-1].X = dummyIndex;

    // Loop
    
    for (int i = finalSize - 2; i >= 0; --i){

    	int candidates[] = { startingPath[i+1].X - 1, startingPath[i+1].X , startingPath[i+1].X + 1};

    	double beta[3] = {0.0,0.0,0.0};
    	double sum = 0.0;

        // Start with arrival
    	if( startingPath[i+1].X > 0){

    		beta[0] = birthRate * ffArray[(N+1)*i + candidates[0]];
    		sum = sum + beta[0];

    	}

        // Continue with virtual

    	beta[1] = K * ( birthRate * (startingPath[i+1].X < N) + deathRate * startingPath[i+1].X * seasonality(startingPath[i+1].time) ) * ffArray[(N+1)*i + candidates[1]];
    	sum = sum + beta[1];

        // Finally death
        if( startingPath[i+1].X < N){
        	
        	beta[2] =  deathRate * candidates[2] * seasonality(startingPath[i+1].time) * ffArray[(N+1)*i + candidates[2]];
        	sum = sum + beta[2];

    	}

        // Normalise and sample

    	for (int j = 0; j < 2; ++j) beta[j] = beta[j] / sum;

    	dummyUnif = unif(generator);

    	if (dummyUnif <= beta[0]){
    		startingPath[i].X = candidates[0];
    	} else if (dummyUnif <= beta[0]+beta[1]){
    		startingPath[i].X = candidates[1];
    	} else {
    		startingPath[i].X = candidates[2];
    	}

    }

    // Remove virtual jumps and return
    counter = 1;
    while (counter < startingPath.size()) {
    
    	if(startingPath[counter].X == startingPath[counter-1].X) {
    		startingPath.erase(startingPath.begin() + counter);
    	} else{
    		counter++;
    	}

    }

}


/*************************/
/* AuxVar for inf states */
/*************************/

void mcmc_newPathAuxVar(vector<BDstate>& startingPath, vector<vector<double>>& observations, vector<vector<double>>& odeApprox, int N, double obsSd, double birthRate, double K, double avDev, double sd, double deathRate, mt19937& generator)
{

	// Fill the new states

	startingPath.push_back({observations.back()[0] , -1}); // to make sure this goes on until last observation time
	int sizePathOriginal = startingPath.size();

	for (int i = 0; i < sizePathOriginal-1; ++i)
	{

		// Get dom rate for virtuals based on K proportionality; then we will thin to get non-homogeneous behwaviour

		double virtualRate = K * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X * 2); // then thin seasonality

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);
		uniform_real_distribution<double> unifProb(0,1);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval, before thinning

		// Times need sampling and thinning, since non-homogeneous process!
		for (int j = 0; j < dummy0; ++j)
		{
			double candidateTime = unif(generator);
			double candidateProb = K * (birthRate * (startingPath[i].X < N)  + deathRate * startingPath[i].X * seasonality(candidateTime)) /virtualRate;
			if (unifProb(generator)<candidateProb) startingPath.push_back({candidateTime , startingPath[i].X}) ; 
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](BDstate & a, BDstate & b) {return a.time < b.time;});

	// Remove last entry, which I added to ensure we reach the end of observation time

	startingPath.pop_back() ;
	
	// Store final size of augmented frame
	int finalSize = startingPath.size();

	// Populate mapping to corresponding ODE value, timeDiffs 

	double timeDiffs[finalSize];
	double timeDiffsSeas[finalSize];
	double ODEvalues[finalSize];

	for (int i = 0; i < finalSize-1; ++i)
	{

		// Time differences recorded

		timeDiffs[i] = startingPath[i+1].time - startingPath[i].time;
		timeDiffsSeas[i] = seasonalityIntegrated(startingPath[i+1].time) - seasonalityIntegrated(startingPath[i].time);

		// Associated ODE

		ODEvalues[i] = odeApprox[1 + startingPath[i].time * 10][1];

	}

	// Add final

	timeDiffs[finalSize-1] = observations.back()[0] - startingPath[finalSize-1].time;
	timeDiffsSeas[finalSize-1] = seasonalityIntegrated(observations.back()[0]) - seasonalityIntegrated(startingPath[finalSize-1].time);
	ODEvalues[finalSize-1] = odeApprox[1 + startingPath[finalSize-1].time * 10][1];

	// Account for observations, at the latest jump right before each observation time

	vector<vector<double>> vectorObs(startingPath.size());
	int counter = 0;
	for (int i = 0; i < observations.size(); ++i)
	{
		while(startingPath[counter].time < observations[i][0] && counter < startingPath.size()) counter++;
		vectorObs[counter-1].push_back(observations[i][1]);
	}

	// Auxiliary variables, plus their means for later use

	double meansAuxvar[finalSize];
	double auxVars[finalSize];
	int range[finalSize*2];

	double invSd = 1/sd;
	normal_distribution<double> normal; // may need to tinker with this
	//double autoregrVal = 0.95; //0.995; // 0.98; //0.95;
	double kappa = 0.05; //0.1;

	// First instance, altered so as to ensure algo is not stuck

	meansAuxvar[0] = max(avDev,abs(startingPath[0].X-ODEvalues[0]));
	normal = normal_distribution<double> (meansAuxvar[0],sd); 
	auxVars[0] = normal(generator) ;
	while(auxVars[0] < abs(startingPath[0].X-ODEvalues[0])) auxVars[0] = normal(generator) ;
	range[0] = max(0,(int) ceil(ODEvalues[0]-auxVars[0]));
	range[finalSize] = min(N,(int) floor(ODEvalues[0]+auxVars[0]));

	for (int i = 1; i < finalSize; ++i)
	{
		// Auxiliary Variable

		// meansAuxvar[i] = avDev * (1-autoregrVal) + auxVars[i-1] * autoregrVal;
		meansAuxvar[i] = max(avDev,  auxVars[i-1] - kappa);

		normal = normal_distribution<double> (meansAuxvar[i],sd); 
		auxVars[i] = normal(generator) ;
		while(auxVars[i] < abs(startingPath[i].X-ODEvalues[i])) auxVars[i] = normal(generator) ;
		range[i] = max(0,(int) ceil(ODEvalues[i]-auxVars[i]));
		range[finalSize + i] = min(N,(int) floor(ODEvalues[i]+auxVars[i]));

	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Stuff for weights

	double weighting1[N+1];
	double weighting2[N+1];
	for (int i = 0; i < N; ++i)	{
		weighting1[i] = - birthRate; 
		weighting2[i] = - deathRate * i; 
	}
	weighting1[N] = 0; 
	weighting2[N] = - deathRate * N; 

	// Build array for FF

	vector<double> ffArray(finalSize * (N+1), 0.0);


	// First iteration
	for (int i = range[0]; i < range[finalSize]+1; ++i) ffArray[i] =  fastExp((1+K) * (weighting1[i]*timeDiffs[0] + weighting2[i]*timeDiffsSeas[0])) / normalCDF( (meansAuxvar[0] - abs(i-ODEvalues[0])) * invSd );

	// Weight by observation likelihood; if any
	for (int k = 0; k < vectorObs[0].size(); ++k) {

		// till upper bound only
		for (int j = range[0]; j < range[finalSize]+1; ++j) ffArray[j] *=  fastExp( - ipow((vectorObs[0][k]-j)/obsSd,2) / 2.0 );

	}

	// Normalise
	double invSum = 1.0 / accumulate (ffArray.begin() + range[0], ffArray.begin() + range[finalSize]+1, 0.0) ;
	for (int j = range[0]; j < range[finalSize]+1; ++j) ffArray[j] *= invSum;

	// Forward filtering steps

	for (int i = 1; i < finalSize; ++i){

		for (int j = range[i]; j < range[finalSize+i]+1 ; ++j)
		{

			double virtualProb = K * ( birthRate + deathRate * j * seasonality(startingPath[i].time) );

			// Transition with weighting

			ffArray[(N+1)*i + j] = fastExp( (1+K) * (weighting1[j]*timeDiffs[i] + weighting2[j]*timeDiffsSeas[i]) ) *  (ffArray[(N+1)*(i-1) + j-1] * birthRate + ffArray[(N+1)*(i-1) + j + 1] * deathRate * (j+1) * seasonality(startingPath[i].time) + ffArray[(N+1)*(i-1) + j] * virtualProb) /
				normalCDF( (meansAuxvar[i] - abs(j-ODEvalues[i])) * invSd );

		}

		// Special case at j=0

		if (range[i] == 0){

			ffArray[(N+1)*i] = fastExp( (1+K) * (weighting1[0]*timeDiffs[i] + weighting2[0]*timeDiffsSeas[i]) ) *
			 	(ffArray[(N+1)*(i-1) + 1] * deathRate * seasonality(startingPath[i].time) + ffArray[(N+1)*(i-1)] * K * birthRate)/
				normalCDF( (meansAuxvar[i] - abs(-ODEvalues[i])) * invSd );

		}

		// Special case at j=N

		if(range[finalSize+i] == N){

			ffArray[(N+1)*i + N] = fastExp( (1+K) * (weighting1[N]*timeDiffs[i] + weighting2[N]*timeDiffsSeas[i]) ) *
			 	(ffArray[(N+1)*(i-1) + N - 1] * birthRate + ffArray[(N+1)*(i-1) + N] * K * deathRate * N * seasonality(startingPath[i].time) )/
				normalCDF( (meansAuxvar[i] - abs(N-ODEvalues[i])) * invSd );

		}	

		// Weight entire row by observation likelihood; if any
		for (int k = 0; k < vectorObs[i].size(); ++k) {

			for (int j = range[i]; j < range[finalSize+i]+1 ; ++j) ffArray[(N+1)*i + j] *=  fastExp( - ipow((vectorObs[i][k]-j)/obsSd,2) / 2.0 );

		}

		// Normalise
		double invSum = 1.0 / accumulate (ffArray.begin() + (N+1)*i + range[i], ffArray.begin() + (N+1)*i + range[finalSize+i]+1, 0.0) ;
		for (int j = range[i]; j < range[finalSize+i]+1 ; ++j) ffArray[(N+1)*i + j] *= invSum;

	}

	/*********************/
	/* Backward Sampling */
	/*********************/
        
	// Sample end value from distribution at end point
	    
	uniform_real_distribution<double> unif(0.0, 1.0);
	double dummyUnif = unif(generator);
	int dummyIndex = -1;
	double dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += ffArray[(N+1)*(finalSize-1) + dummyIndex];
   	}
    startingPath[finalSize-1].X = dummyIndex;

    // Loop
    
    for (int i = finalSize - 2; i >= 0; --i){
        
        int candidates[] = { startingPath[i+1].X - 1, startingPath[i+1].X , startingPath[i+1].X + 1};

        double beta[3] = {0.0,0.0,0.0};
        double sum = 0.0;

        // Start with arrival
        if( startingPath[i+1].X > 0){
        	
        	beta[0] = birthRate * ffArray[(N+1)*i + candidates[0]];
        	sum = sum + beta[0];

    	}

        // Continue with virtual

        beta[1] = K * ( birthRate * (startingPath[i+1].X < N) + deathRate * startingPath[i+1].X * seasonality(startingPath[i+1].time) ) * ffArray[(N+1)*i + candidates[1]];
        sum = sum + beta[1];

        // Finally death
        if( startingPath[i+1].X < N){
        	
        	beta[2] =  deathRate * candidates[2] * seasonality(startingPath[i+1].time) * ffArray[(N+1)*i + candidates[2]];
        	sum = sum + beta[2];

    	}

        // Normalise and sample

        for (int j = 0; j < 2; ++j) beta[j] = beta[j] / sum;

		dummyUnif = unif(generator);

        if (dummyUnif <= beta[0]){
        	startingPath[i].X = candidates[0];
        } else if (dummyUnif <= beta[0]+beta[1]){
        	startingPath[i].X = candidates[1];
        } else {
        	startingPath[i].X = candidates[2];
        }
      
    }

    // Remove virtual jumps and return
    counter = 1;
    while (counter < startingPath.size()) {
    
    	if(startingPath[counter].X == startingPath[counter-1].X) {
    		startingPath.erase(startingPath.begin() + counter);
    	} else{
    		counter++;
    	}

    }

}


/*******************/
/* MCMC Procedures */
/*******************/

void do_mcmc_Unif(vector<vector<double>> observations, int N, double obsSd, double birthRate, int mcmcIter, double scaling, mt19937 generator, string filename, int chainNumber)
{
	// Create starting trajectory

	vector<BDstate> startingPath; 
	mcmc_startingPathVector(startingPath, observations, N, generator);

	// Prior for death rate, birth assumed known

	double priorShape = 1.0, priorRate = pow(10.0,-3.0);

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "mu \n";

	// ITERATE ALGORITHM

	for (int i = 0; i < mcmcIter; ++i)
	{

		// Rates update

		double shape = priorShape; // amount of deaths in trajectory + prior
		double rate = priorRate; // scaled time allowing a "single" death to happen + prior
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X == startingPath[j].X - 1) shape++;
			rate = rate + (seasonalityIntegrated(startingPath[j+1].time) - seasonalityIntegrated(startingPath[j].time)) * startingPath[j].X;
		}

		gamma_distribution<double> gamma(shape,1.0/rate);
		double deathRate = gamma(generator);

      	// Save parameters and update trajectory

		myfileTrace << deathRate <<"\n";

        // Update trajectory

		mcmc_newPathUnif(startingPath, observations, N, obsSd, birthRate, scaling, deathRate, generator);

		// // Progress bar

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
	// myfile << "Time, X \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].X << "\n";
	// }
	// myfile.close();

}

void do_mcmc_DepThinning(vector<vector<double>> observations, int N, double obsSd, double birthRate, int mcmcIter, double K, mt19937 generator, string filename, int chainNumber)
{

	// Create starting trajectory

	vector<BDstate> startingPath; 
	mcmc_startingPathVector(startingPath, observations, N, generator);

	// Prior for death rate, birth assumed known

	double priorShape = 1.0, priorRate = pow(10.0,-3.0);

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "mu \n";

	// ITERATE ALGORITHM

	for (int i = 0; i < mcmcIter; ++i)
	{


		// Rates update

		double shape = priorShape; // amount of deaths in trajectory + prior
		double rate = priorRate; // scaled time allowing a "single" death to happen + prior
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X == startingPath[j].X - 1) shape++;
			rate = rate + (seasonalityIntegrated(startingPath[j+1].time) - seasonalityIntegrated(startingPath[j].time)) * startingPath[j].X;
		}

		gamma_distribution<double> gamma(shape,1.0/rate);
		double deathRate = gamma(generator);

      	// Save parameters and update trajectory

		myfileTrace << deathRate <<"\n";

        // Update trajectory

		mcmc_newPathDepThinning(startingPath, observations, N, obsSd, birthRate, K, deathRate, generator);

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
	// myfile << "Time, X \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].X << "\n";
	// }
	// myfile.close();

}

void do_mcmc_DepThinningCounts(vector<vector<double>> observations, int N, double obsSd, double birthRate, int mcmcIter, double K, mt19937 generator, string filename, int chainNumber)
{

	// Create starting trajectory

	vector<BDstate> startingPath; 
	mcmc_startingPathVector(startingPath, observations, N, generator);

	// Prior for death rate, birth assumed known

	double priorShape = 1.0, priorRate = pow(10.0,-3.0);

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "mu \n";

	// ITERATE ALGORITHM

	for (int i = 0; i < mcmcIter; ++i)
	{


		// Rates update

		double shape = priorShape; // amount of deaths in trajectory + prior
		double rate = priorRate; // scaled time allowing a "single" death to happen + prior
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X == startingPath[j].X - 1) shape++;
			rate = rate + (seasonalityIntegrated(startingPath[j+1].time) - seasonalityIntegrated(startingPath[j].time)) * startingPath[j].X;
		}

		gamma_distribution<double> gamma(shape,1.0/rate);
		double deathRate = gamma(generator);

      	// Save parameters and update trajectory

		myfileTrace << deathRate <<"\n";

        // Update trajectory

		mcmc_newPathDepThinningCounts(startingPath, observations, N, obsSd, birthRate, K, deathRate, generator);

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
	// myfile << "Time, X \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].X << "\n";
	// }
	// myfile.close();

}

void do_mcmc_AuxVar(vector<vector<double>> observations, vector<vector<double>> odeApprox, int N, double obsSd, double birthRate, int mcmcIter, double K, double avDev, double sd, mt19937 generator, string filename, int chainNumber)
{

	// Create starting trajectory

	vector<BDstate> startingPath; 
	mcmc_startingPathVector2(startingPath, observations, odeApprox, N, generator);

	// Prior for death rate, birth assumed known

	double priorShape = 1.0, priorRate = pow(10.0,-3.0);

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "mu \n";

	// First iteration is just to ensure we get a starting trajectory, since starting path gives an inappropriate thing

	double shape = priorShape; // amount of deaths in trajectory + prior
	double rate = priorRate; // scaled time allowing a "single" death to happen + prior
		
	for (int j = 0; j < startingPath.size()-1; ++j)
	{
		if (startingPath[j+1].X == startingPath[j].X - 1) shape++;
		rate = rate + (seasonalityIntegrated(startingPath[j+1].time) - seasonalityIntegrated(startingPath[j].time)) * startingPath[j].X;
	}

	gamma_distribution<double> gamma(shape,1.0/rate);
	double deathRate = gamma(generator);

	mcmc_newPathDepThinning(startingPath, observations, N, obsSd, birthRate, K, deathRate, generator);

	// ITERATE ALGORITHM

	for (int i = 0; i < mcmcIter; ++i)
	{

		// Rates update

		shape = priorShape; // amount of deaths in trajectory + prior
		rate = priorRate; // scaled time allowing a "single" death to happen + prior
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X == startingPath[j].X - 1) shape++;
			rate = rate + (seasonalityIntegrated(startingPath[j+1].time) - seasonalityIntegrated(startingPath[j].time)) * startingPath[j].X;
		}

		gamma_distribution<double> gamma(shape,1.0/rate);
		double deathRate = gamma(generator);

      	// Save parameters and update trajectory

		myfileTrace << deathRate <<"\n";

        // Update trajectory

		mcmc_newPathAuxVar(startingPath, observations, odeApprox, N, obsSd, birthRate, K, avDev, sd, deathRate, generator);

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
	// myfile << "Time, X \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].X << "\n";
	// }
	// myfile.close();

}