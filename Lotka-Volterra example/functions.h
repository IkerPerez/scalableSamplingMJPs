
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
			getline (myfile,line, ' ');
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

void mcmc_startingPathVector(vector<LVstate>& startingPath, vector<vector<double>>& observations, int N, mt19937& gen)
{
	// Create and fill a frame to operate within

	double sum1 = 0, sum2 =0;
	for (int i = 0; i < observations.size(); ++i) {
		sum1 += observations[i][1];
		sum2 += observations[i][2];
	}
	int avgState1 = (int) sum1 / observations.size(), avgState2 = (int) sum2 / observations.size();
	startingPath.push_back({0.0,avgState1,avgState2});
	double jumpswithin =  max(1.0,N/100.0);
	for (int i = 1; i < observations.size(); ++i) {
		for (int j = 0; j < jumpswithin; ++j)
		{
			startingPath.push_back({observations[i][0] + ( - observations[i][0] + observations[i-1][0]) * (j+1)/(jumpswithin+1),avgState1,avgState2});
		}
	}

}


/**************************/
/* Vanilla Uniformization */
/**************************/

void mcmc_newPathUnif(vector<LVstate>& startingPath, vector<vector<double>>& observations, int N, double obsSd, double alpha, double beta, double delta, double gamma, double scaling, mt19937& generator, 	vector<vector<vector<double>>> obsWeightsMat)
{

	// Dominating rate

	double omega = scaling * (N * (alpha+gamma) + (beta + delta) * ipow(N,2));
	double invOmega = 1.0/omega;

    // Fill the new states

	startingPath.push_back({observations.back()[0] , -1, -1}); // to make sure this goes on until last observation time
	int sizePathOriginal = startingPath.size();
	for (int i = 0; i < sizePathOriginal-1; ++i)
	{

		// Compute dominating virtual rate

		double virtualRate = omega - alpha * startingPath[i].X1 * (startingPath[i].X1 < N) - beta * startingPath[i].X1 * startingPath[i].X2	- delta * startingPath[i].X1 * startingPath[i].X2 * (startingPath[i].X2 < N) - gamma* startingPath[i].X2;

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval, before thinning

		for (int j = 0; j < dummy0; ++j)
		{
			double candidateTime = unif(generator);
			startingPath.push_back({candidateTime , startingPath[i].X1, startingPath[i].X2}) ; 
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](LVstate & a, LVstate & b) {return a.time < b.time;});

	// Remove last entry, which I added to ensure we reach the end of observation time

	startingPath.pop_back() ;

	// Account for observations, at the latest jump right before each observation time, there could be more than one!

	vector<vector<double>> indexesObs(startingPath.size());
	int counter = 0;
	for (int i = 0; i < observations.size(); ++i)
	{
		while(startingPath[counter].time < observations[i][0] && counter < startingPath.size()) counter++;
		indexesObs[counter-1].push_back(i);
	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Build array of matrices for FF

	vector<vector<vector<double>>> ffArray(startingPath.size());
	for ( int i = 0 ; i < startingPath.size() ; i++ ){
		ffArray[i].resize(N+1);
		for (int j = 0; j < N+1; ++j)
		{
			ffArray[i][j].resize(N+1);
		}
	}

	// First iteration

	for (int i = 0; i < N + 1; ++i) {
		for (int j = 0; j < N+1; ++j)
		{
			ffArray[0][i][j] = 1.0/( (double) N+1)/( (double) N+1);
		}
		
	}

	// Weight by observation likelihood; if any

	for (int k = 0; k < indexesObs[0].size(); ++k) {

		for (int i = 0; i < N+1; ++i){

			for (int j = 0; j < N+1; ++j){

				ffArray[0][i][j] *=  obsWeightsMat[indexesObs[0][k]][i][j];

			}

		}

	}

	// Normalise
	double sum = 0;
	for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) sum += ffArray[0][i][j];
	for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) ffArray[0][i][j] *= 1.0/sum;

	// Matrices of transition intensities for each jump type

	vector< vector<double> > virtualFire(N+1, vector<double>(N+1));
	for (int i = 0; i < N+1; ++i)	{
		for (int j = 0; j < N+1; ++j)	{
			virtualFire[i][j] = 1 - (alpha * i * (i < N) + beta * i * j	+ delta * i * j * (j < N) + gamma * j) * invOmega ; 
		}
	}

	// Forward filtering steps

	for (int r = 1; r < startingPath.size(); ++r){

		for (int i = 1; i < N; ++i){

			for (int j = 1; j < N; ++j){

				ffArray[r][i][j] =  virtualFire[i][j] * ffArray[r-1][i][j] + invOmega* (

					alpha * (i-1) * ffArray[r-1][i-1][j] + beta * (i+1) * j * ffArray[r-1][i+1][j] + delta * i * (j-1) * ffArray[r-1][i][j-1] +	gamma * (j+1) * ffArray[r-1][i][j+1] );

			}

		}

		// Special cases i = 0

		for (int j = 1; j < N; ++j){

			ffArray[r][0][j] =  virtualFire[0][j] * ffArray[r-1][0][j] + invOmega* ( beta * j * ffArray[r-1][1][j] + gamma * (j+1) * ffArray[r-1][0][j+1] );

		}
		ffArray[r][0][0] = ffArray[r-1][0][0] + invOmega* gamma * ffArray[r-1][0][1];
		ffArray[r][0][N] = (1 - gamma * N * invOmega) * ffArray[r-1][0][N] + invOmega * beta * N * ffArray[r-1][1][N];

		// Special cases i = N

		for (int j = 1; j < N; ++j){

			ffArray[r][N][j] =  virtualFire[N][j] * ffArray[r-1][N][j] + invOmega* (

				alpha * (N-1) * ffArray[r-1][N-1][j] + delta * N * (j-1) * ffArray[r-1][N][j-1] + gamma * (j+1) * ffArray[r-1][N][j+1] );

		}
		ffArray[r][N][0] =  ffArray[r-1][N][0] + invOmega* (alpha * (N-1) * ffArray[r-1][N-1][0] + gamma * ffArray[r-1][N][1] );
		ffArray[r][N][N] =  (1 - ( beta * N * N + gamma * N ) * invOmega) * ffArray[r-1][N][N] + invOmega* (alpha * (N-1) * ffArray[r-1][N-1][N] + delta * N * (N-1) * ffArray[r-1][N][N-1] );

		// Weight by observation likelihood; if any

		for (int k = 0; k < indexesObs[r].size(); ++k) {

			for (int i = 0; i < N+1; ++i){

				for (int j = 0; j < N+1; ++j){

					ffArray[r][i][j] *=  obsWeightsMat[indexesObs[r][k]][i][j];

				}

			}

		}

		// Normalise
		double sum = 0;
		for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) sum += ffArray[r][i][j];
		for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) ffArray[r][i][j] *= 1.0/sum;

	}


	/*********************/
	/* Backward Sampling */
	/*********************/
        
	// Sample end value from distribution at end point
	    
	uniform_real_distribution<double> unif(0.0, 1.0);
	double dummyUnif = unif(generator);
	int dummyIndex = -1;

	vector<double> marginalPrey(N+1); for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) marginalPrey[i] += ffArray[startingPath.size()-1][i][j];

	double dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += marginalPrey[dummyIndex];
   	}
    startingPath[startingPath.size()-1].X1 = dummyIndex;
	
    double normConst = 0;
	vector<double> marginalPred(N+1); 
	for (int j = 0; j < N+1; ++j) {
		marginalPred[j] = ffArray[startingPath.size()-1][dummyIndex][j];
		normConst += marginalPred[j];
	}
	for (int j = 0; j < N+1; ++j) marginalPred[j] *= 1.0/normConst; 

	dummyUnif = unif(generator);
	dummyIndex = -1;
	dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += marginalPred[dummyIndex];
   	}
    startingPath[startingPath.size()-1].X2 = dummyIndex;

    // Loop
    
    for (int i = startingPath.size() - 2; i >= 0; --i){
   
      	vector<LVstate> candidates;
      	vector<double> probsPrev;
      	double sum = 0.0;

      	// Start with virtual

      	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2});

		double virtualJumpProb = 1 - ( alpha * startingPath[i+1].X1 * (startingPath[i+1].X1 < N) +
			beta * startingPath[i+1].X1 * startingPath[i+1].X2 + delta * startingPath[i+1].X1 * startingPath[i+1].X2 * (startingPath[i+1].X2 < N) + 
 			gamma* startingPath[i+1].X2 ) * invOmega;

      	probsPrev.push_back( virtualJumpProb * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2] );

        sum = sum + probsPrev.back();

        // Prey born

        if(startingPath[i+1].X1>0){

        	candidates.push_back({0,startingPath[i+1].X1-1,startingPath[i+1].X2});

        	probsPrev.push_back( alpha * invOmega * (startingPath[i+1].X1-1) * ffArray[i][startingPath[i+1].X1-1][startingPath[i+1].X2] );

        	sum = sum + probsPrev.back();

        }

        // Prey death

        if(startingPath[i+1].X1<N){

        	candidates.push_back({0,startingPath[i+1].X1+1,startingPath[i+1].X2});

        	probsPrev.push_back( beta * invOmega * (startingPath[i+1].X1+1) * startingPath[i+1].X2 * ffArray[i][startingPath[i+1].X1+1][startingPath[i+1].X2] );

        	sum = sum + probsPrev.back();

        }
      
       	// Pred birth

        if(startingPath[i+1].X2>0){

        	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2-1});

        	probsPrev.push_back( delta * invOmega * startingPath[i+1].X1 * (startingPath[i+1].X2-1) * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2-1] );

        	sum = sum + probsPrev.back();

        }

  		// Pred death

        if(startingPath[i+1].X2<N){

        	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2+1});

        	probsPrev.push_back( gamma * invOmega * (startingPath[i+1].X2+1) * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2+1] );

        	sum = sum + probsPrev.back();

        }

        // Normalise and sample

        for (int j = 0; j < probsPrev.size(); ++j) probsPrev[j] = probsPrev[j] / sum;

		dummyUnif = unif(generator);
		dummyIndex = -1;
		dummySum = 0;
		while(dummySum < dummyUnif) {
			dummyIndex++;
      		dummySum += probsPrev[dummyIndex];
   		}
    	startingPath[i].X1 = candidates[dummyIndex].X1;
    	startingPath[i].X2 = candidates[dummyIndex].X2;
      
    }

    // Remove virtual jumps and return
    counter = 1;
    while (counter < startingPath.size()) {
    
    	if(startingPath[counter].X1 == startingPath[counter-1].X1 && startingPath[counter].X2 == startingPath[counter-1].X2) {
    		startingPath.erase(startingPath.begin() + counter);
    	} else{
    		counter++;
    	}

    }

}


/******************/
/* Depen Thinning */
/******************/

void mcmc_newPathDepThinning(vector<LVstate>& startingPath, vector<vector<double>>& observations, int N, double obsSd, double alpha, double beta, double delta, double gamma, double K, mt19937& generator, vector<vector<vector<double>>> obsWeightsMat)
{

    // Fill the new states

	startingPath.push_back({observations.back()[0] , -1, -1}); // to make sure this goes on until last observation time
	int sizePathOriginal = startingPath.size();
	for (int i = 0; i < sizePathOriginal-1; ++i)
	{

		// Get dom rate for virtuals based on K proportionality

		double virtualRate = K * ( alpha * startingPath[i].X1 * (startingPath[i].X1 < N) + beta * startingPath[i].X1 * startingPath[i].X2	+ 
			delta * startingPath[i].X1 * startingPath[i].X2 * (startingPath[i].X2 < N) + gamma* startingPath[i].X2); 

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval

		// add times

		for (int j = 0; j < dummy0; ++j)
		{
			double candidateTime = unif(generator);
			startingPath.push_back({candidateTime , startingPath[i].X1, startingPath[i].X2}) ; 
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](LVstate & a, LVstate & b) {return a.time < b.time;});

	// Remove last entry, which I added to ensure we reach the end of observation time

	startingPath.pop_back();
	
	// Store final size of augmented frame
	int finalSize = startingPath.size();

	// populate time differences

	double timeDiffs[finalSize];

	for (int i = 0; i < finalSize-1; ++i) timeDiffs[i] = startingPath[i+1].time - startingPath[i].time;

	// Add final

	timeDiffs[finalSize-1] = observations.back()[0] - startingPath[finalSize-1].time;

	// Account for observations, at the latest jump right before each observation time, there could be more than one!

	vector<vector<double>> indexesObs(startingPath.size());
	int counter = 0;
	for (int i = 0; i < observations.size(); ++i)
	{
		while(startingPath[counter].time < observations[i][0] && counter < startingPath.size()) counter++;
		indexesObs[counter-1].push_back(i);
	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Stuff for weights and transitions; to speed things up

	vector< vector<double> > weighting(N+1, vector<double>(N+1));
	for (int i = 0; i < N+1; ++i)	{
		for (int j = 0; j < N+1; ++j)	{
			weighting[i][j] = - (alpha * i * (i < N) + beta * i * j	+ delta * i * j * (j < N) + gamma * j); 
		}
	}

	// Build array of matrices for FF

	vector<vector<vector<double>>> ffArray(startingPath.size());
	for ( int i = 0 ; i < startingPath.size() ; i++ ){
		ffArray[i].resize(N+1);
		for (int j = 0; j < N+1; ++j)
		{
			ffArray[i][j].resize(N+1);
		}
	}

	// First iteration

	for (int i = 0; i < N + 1; ++i) {
		for (int j = 0; j < N+1; ++j)
		{
			ffArray[0][i][j] = fastExp((1+K) * weighting[i][j] * timeDiffs[0]);
		}
		
	}

	// Weight by observation likelihood; if any

	for (int k = 0; k < indexesObs[0].size(); ++k) {

		for (int i = 0; i < N+1; ++i){

			for (int j = 0; j < N+1; ++j){

				ffArray[0][i][j] *=  obsWeightsMat[indexesObs[0][k]][i][j];

			}

		}

	}

	// Normalise
	double sum = 0;
	for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) sum += ffArray[0][i][j];
	for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) ffArray[0][i][j] *= 1.0/sum;

	// Forward filtering iterations

	for (int r = 1; r < finalSize; ++r){

		for (int i = 1; i < N; ++i){

			for (int j = 1; j < N; ++j){

				ffArray[r][i][j] =  fastExp( (1+K) * weighting[i][j] * timeDiffs[r]) * 
					(- K * weighting[i][j] * ffArray[r-1][i][j] + alpha * (i-1) * ffArray[r-1][i-1][j] + 
					beta * (i+1) * j * ffArray[r-1][i+1][j] + delta * i * (j-1) * ffArray[r-1][i][j-1] +	gamma * (j+1) * ffArray[r-1][i][j+1] );

			}

		}

		// Special cases i = 0

		for (int j = 1; j < N; ++j){

			ffArray[r][0][j] =  fastExp( (1+K) * weighting[0][j] * timeDiffs[r] ) *
				( - K * weighting[0][j] * ffArray[r-1][0][j] + beta * j * ffArray[r-1][1][j] + gamma * (j+1) * ffArray[r-1][0][j+1] );

		}
		ffArray[r][0][0] = fastExp( (1+K) * weighting[0][0] * timeDiffs[r] ) * (- K * weighting[0][0] * ffArray[r-1][0][0] + gamma * ffArray[r-1][0][1]);
		ffArray[r][0][N] = fastExp( (1+K) * weighting[0][N] * timeDiffs[r] ) * (- K * weighting[0][N] * ffArray[r-1][0][N] + beta * N * ffArray[r-1][1][N] );

		// Special cases i = N

		for (int j = 1; j < N; ++j){

			ffArray[r][N][j] =  fastExp( (1+K) * weighting[N][j] * timeDiffs[r] ) *
				(- K * weighting[N][j] * ffArray[r-1][N][j] + alpha * (N-1) * ffArray[r-1][N-1][j] + delta * N * (j-1) * ffArray[r-1][N][j-1] + gamma * (j+1) * ffArray[r-1][N][j+1] );

		}
		ffArray[r][N][0] =  fastExp( (1+K) * weighting[N][0] * timeDiffs[r] ) * (- K * weighting[N][0] * ffArray[r-1][N][0] + alpha * (N-1) * ffArray[r-1][N-1][0] + gamma * ffArray[r-1][N][1] );
		ffArray[r][N][N] =  fastExp( (1+K) * weighting[N][N] * timeDiffs[r] ) * (- K * weighting[N][N] * ffArray[r-1][N][N] + alpha * (N-1) * ffArray[r-1][N-1][N] + delta * N * (N-1) * ffArray[r-1][N][N-1] );

		// Weight by observation likelihood; if any

		for (int k = 0; k < indexesObs[r].size(); ++k) {

			for (int i = 0; i < N+1; ++i){

				for (int j = 0; j < N+1; ++j){

					ffArray[r][i][j] *=  obsWeightsMat[indexesObs[r][k]][i][j];

				}

			}

		}

		// Normalise
		double sum = 0;
		for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) sum += ffArray[r][i][j];
		for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) ffArray[r][i][j] *= 1.0/sum;

	}


	/*********************/
	/* Backward Sampling */
	/*********************/
        
	// Sample end value from distribution at end point
	    
	uniform_real_distribution<double> unif(0.0, 1.0);
	double dummyUnif = unif(generator);
	int dummyIndex = -1;

	vector<double> marginalPrey(N+1); for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) marginalPrey[i] += ffArray[startingPath.size()-1][i][j];

	double dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += marginalPrey[dummyIndex];
   	}
    startingPath[startingPath.size()-1].X1 = dummyIndex;
	
    double normConst = 0;
	vector<double> marginalPred(N+1); 
	for (int j = 0; j < N+1; ++j) {
		marginalPred[j] = ffArray[startingPath.size()-1][dummyIndex][j];
		normConst += marginalPred[j];
	}
	for (int j = 0; j < N+1; ++j) marginalPred[j] *= 1.0/normConst; 

	dummyUnif = unif(generator);
	dummyIndex = -1;
	dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += marginalPred[dummyIndex];
   	}
    startingPath[startingPath.size()-1].X2 = dummyIndex;

    // Loop
    
    for (int i = finalSize - 2; i >= 0; --i){

		vector<LVstate> candidates;
      	vector<double> probsPrev;
      	double sum = 0.0;

      	// Start with virtual

      	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2});

      	probsPrev.push_back(- K * weighting[startingPath[i+1].X1][startingPath[i+1].X2] * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2] );

        sum = sum + probsPrev.back();

        // Prey born

        if(startingPath[i+1].X1>0){

        	candidates.push_back({0,startingPath[i+1].X1-1,startingPath[i+1].X2});

        	probsPrev.push_back( alpha * (startingPath[i+1].X1-1) * ffArray[i][startingPath[i+1].X1-1][startingPath[i+1].X2] );

        	sum = sum + probsPrev.back();

        }

        // Prey death

        if(startingPath[i+1].X1<N){

        	candidates.push_back({0,startingPath[i+1].X1+1,startingPath[i+1].X2});

        	probsPrev.push_back( beta * (startingPath[i+1].X1+1) * startingPath[i+1].X2 * ffArray[i][startingPath[i+1].X1+1][startingPath[i+1].X2] );

        	sum = sum + probsPrev.back();

        }
      
       	// Pred birth

        if(startingPath[i+1].X2>0){

        	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2-1});

        	probsPrev.push_back( delta * startingPath[i+1].X1 * (startingPath[i+1].X2-1) * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2-1] );

        	sum = sum + probsPrev.back();

        }

  		// Pred death

        if(startingPath[i+1].X2<N){

        	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2+1});

        	probsPrev.push_back( gamma * (startingPath[i+1].X2+1) * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2+1] );

        	sum = sum + probsPrev.back();

        }

        // Normalise and sample

        for (int j = 0; j < probsPrev.size(); ++j) probsPrev[j] = probsPrev[j] / sum;

		dummyUnif = unif(generator);
		dummyIndex = -1;
		dummySum = 0;
		while(dummySum < dummyUnif) {
			dummyIndex++;
      		dummySum += probsPrev[dummyIndex];
   		}
    	startingPath[i].X1 = candidates[dummyIndex].X1;
    	startingPath[i].X2 = candidates[dummyIndex].X2;

    }

    // Remove virtual jumps and return
    counter = 1;
    while (counter < startingPath.size()) {
    
    	if(startingPath[counter].X1 == startingPath[counter-1].X1 && startingPath[counter].X2 == startingPath[counter-1].X2) {
    		startingPath.erase(startingPath.begin() + counter);
    	} else{
    		counter++;
    	}

    }

}


/*************************/
/* AuxVar for inf states */
/*************************/

void mcmc_newPathAuxVar(vector<LVstate>& startingPath, vector<vector<double>>& observations, int N, double obsSd, double alpha, double beta, double delta, double gamma, double K, int lag, mt19937& generator, vector<vector<vector<double>>> obsWeightsMat)
{

	// Fill the new states

	startingPath.push_back({observations.back()[0] , -1, -1}); // to make sure this goes on until last observation time
	int sizePathOriginal = startingPath.size();

	for (int i = 0; i < sizePathOriginal-1; ++i)
	{

		// Get dom rate for virtuals based on K proportionality

		double virtualRate = K * ( alpha * startingPath[i].X1 * (startingPath[i].X1 < N) + beta * startingPath[i].X1 * startingPath[i].X2	+ 
			delta * startingPath[i].X1 * startingPath[i].X2 * (startingPath[i].X2 < N) + gamma* startingPath[i].X2); 

		// define distributions

		poisson_distribution<int> pois((startingPath[i+1].time - startingPath[i].time) * virtualRate);
		uniform_real_distribution<double> unif(startingPath[i].time, startingPath[i+1].time);

		int dummy0 = pois(generator); // this is amount of nodes/times we will have in this interval

		// add times

		for (int j = 0; j < dummy0; ++j)
		{
			double candidateTime = unif(generator);
			startingPath.push_back({candidateTime , startingPath[i].X1, startingPath[i].X2}) ; 
		}
	
	}

	// sort the data out according to jump times; use lambda expressions
	
	sort(startingPath.begin(), startingPath.end(), [](LVstate & a, LVstate & b) {return a.time < b.time;});

	// Remove last entry, which I added to ensure we reach the end of observation time

	startingPath.pop_back() ;
	
	// Store final size of augmented frame

	int finalSize = startingPath.size();

	// Populate mapping to corresponding ODE value, timeDiffs 

	double timeDiffs[finalSize];
	for (int i = 0; i < finalSize-1; ++i) timeDiffs[i] = startingPath[i+1].time - startingPath[i].time;

	// Add final

	timeDiffs[finalSize-1] = observations.back()[0] - startingPath[finalSize-1].time;

	// Account for observations, at the latest jump right before each observation time, there could be more than one!

	vector<vector<double>> indexesObs(startingPath.size());
	int counter = 0;
	for (int i = 0; i < observations.size(); ++i)
	{
		while(startingPath[counter].time < observations[i][0] && counter < startingPath.size()) counter++;
		indexesObs[counter-1].push_back(i);
	}


	/**********************/
	/* Auxiliary Variables*/
	/**********************/

	int rangePrey[finalSize*2];
	int rangePred[finalSize*2];

	uniform_int_distribution<int> startEpoch(0,lag);
	int dummyEpoch = startEpoch(generator);

	if(dummyEpoch < lag){
		rangePrey[0] = 0; rangePrey[finalSize] = N;
		rangePred[0] = 0; rangePred[finalSize] = N;
	} else{
		rangePrey[0] = startingPath[0].X1; rangePrey[finalSize] = startingPath[0].X1;
		rangePred[0] = startingPath[0].X2; rangePred[finalSize] = startingPath[0].X2;
	}

	for (int i = 1; i < finalSize; ++i)
	{
		
		dummyEpoch++;

		if(dummyEpoch < lag){
			rangePrey[i] = max(0,rangePrey[i-1]-1); rangePrey[finalSize + i] = min(N,rangePrey[finalSize + i-1]+1);
			rangePred[i] = max(0,rangePred[i-1]-1); rangePred[finalSize + i] = min(N,rangePred[finalSize + i-1]+1);
		} else{
			dummyEpoch = 0;
			rangePrey[i] = startingPath[i].X1; rangePrey[finalSize + i] = startingPath[i].X1;
			rangePred[i] = startingPath[i].X2; rangePred[finalSize + i] = startingPath[i].X2;
		}

	}


	/*********************/
	/* Forward Filtering */
	/*********************/

	// Stuff for weights and transitions; to speed things up

	vector< vector<double> > weighting(N+1, vector<double>(N+1));
	for (int i = 0; i < N+1; ++i)	{
		for (int j = 0; j < N+1; ++j)	{
			weighting[i][j] = - (alpha * i * (i < N) + beta * i * j	+ delta * i * j * (j < N) + gamma * j); 
		}
	}

	// Build array of matrices for FF

	vector<vector<vector<double>>> ffArray(startingPath.size());
	for ( int i = 0 ; i < startingPath.size() ; i++ ){
		ffArray[i].resize(N+1);
		for (int j = 0; j < N+1; ++j)
		{
			ffArray[i][j].resize(N+1);
		}
	}

	// First iteration

	for (int i = rangePrey[0]; i < rangePrey[finalSize] + 1; ++i) {
		for (int j = rangePred[0]; j < rangePred[finalSize]+1; ++j)
		{
			ffArray[0][i][j] = fastExp((1+K) * weighting[i][j] * timeDiffs[0]);
		}
		
	}

	// Weight by observation likelihood; if any

	for (int k = 0; k < indexesObs[0].size(); ++k) {

		for (int i = rangePrey[0]; i < rangePrey[finalSize] + 1; ++i){

			for (int j = rangePred[0]; j < rangePred[finalSize]+1; ++j){

				ffArray[0][i][j] *=  obsWeightsMat[indexesObs[0][k]][i][j];

			}

		}

	}

	// Normalise
	double sum = 0;
	for (int i = rangePrey[0]; i < rangePrey[finalSize] + 1; ++i) for (int j = rangePred[0]; j < rangePred[finalSize]+1; ++j) sum += ffArray[0][i][j];
	for (int i = rangePrey[0]; i < rangePrey[finalSize] + 1; ++i) for (int j = rangePred[0]; j < rangePred[finalSize]+1; ++j) ffArray[0][i][j] *= 1.0/sum;

	// Forward filtering steps

	for (int r = 1; r < finalSize; ++r){

		for (int i = max(rangePrey[r],1); i < min(rangePrey[finalSize + r] + 1,N); ++i){

			for (int j = max(rangePred[r],1); j < min(rangePred[finalSize + r] + 1,N); ++j){

				ffArray[r][i][j] =  fastExp( (1+K) * weighting[i][j] * timeDiffs[r]) * 
					(- K * weighting[i][j] * ffArray[r-1][i][j] + alpha * (i-1) * ffArray[r-1][i-1][j] + 
					beta * (i+1) * j * ffArray[r-1][i+1][j] + delta * i * (j-1) * ffArray[r-1][i][j-1] +	gamma * (j+1) * ffArray[r-1][i][j+1] );

			}

		}

		// Special cases i = 0

		if (rangePrey[r] == 0){

			for (int j = max(rangePred[r],1); j < min(rangePred[finalSize + r] + 1,N); ++j){

				ffArray[r][0][j] =  fastExp( (1+K) * weighting[0][j] * timeDiffs[r] ) *
				( - K * weighting[0][j] * ffArray[r-1][0][j] + beta * j * ffArray[r-1][1][j] + gamma * (j+1) * ffArray[r-1][0][j+1] );

			}
			if (rangePred[r] == 0) ffArray[r][0][0] = fastExp( (1+K) * weighting[0][0] * timeDiffs[r] ) * (- K * weighting[0][0] * ffArray[r-1][0][0] + gamma * ffArray[r-1][0][1]);
			if (rangePred[finalSize + r] == N) ffArray[r][0][N] = fastExp( (1+K) * weighting[0][N] * timeDiffs[r] ) * (- K * weighting[0][N] * ffArray[r-1][0][N] + beta * N * ffArray[r-1][1][N] );

		}


		// Special cases i = N

		if (rangePrey[finalSize + r] == N){

			for (int j = max(rangePred[r],1); j < min(rangePred[finalSize + r] + 1,N); ++j){

				ffArray[r][N][j] =  fastExp( (1+K) * weighting[N][j] * timeDiffs[r] ) *
				(- K * weighting[N][j] * ffArray[r-1][N][j] + alpha * (N-1) * ffArray[r-1][N-1][j] + delta * N * (j-1) * ffArray[r-1][N][j-1] + gamma * (j+1) * ffArray[r-1][N][j+1] );

			}
			if (rangePred[r] == 0) ffArray[r][N][0] =  fastExp( (1+K) * weighting[N][0] * timeDiffs[r] ) * (- K * weighting[N][0] * ffArray[r-1][N][0] + alpha * (N-1) * ffArray[r-1][N-1][0] + gamma * ffArray[r-1][N][1] );
			if (rangePred[finalSize + r] == N) ffArray[r][N][N] =  fastExp( (1+K) * weighting[N][N] * timeDiffs[r] ) * (- K * weighting[N][N] * ffArray[r-1][N][N] + alpha * (N-1) * ffArray[r-1][N-1][N] + delta * N * (N-1) * ffArray[r-1][N][N-1] );
		}

		// Weight by observation likelihood; if any

		for (int k = 0; k < indexesObs[r].size(); ++k) {

			for (int i = rangePrey[r]; i < rangePrey[finalSize + r] + 1; ++i){

				for (int j = rangePred[r]; j < rangePred[finalSize + r]+1; ++j){

					ffArray[r][i][j] *=  obsWeightsMat[indexesObs[r][k]][i][j];

				}

			}

		}

		// Normalise
		double sum = 0;
		for (int i = rangePrey[r]; i < rangePrey[finalSize+r] + 1; ++i) for (int j = rangePred[r]; j < rangePred[finalSize+r]+1; ++j) sum += ffArray[r][i][j];
		for (int i = rangePrey[r]; i < rangePrey[finalSize+r] + 1; ++i) for (int j = rangePred[r]; j < rangePred[finalSize+r]+1; ++j) ffArray[r][i][j] *= 1.0/sum;

	}

	/*********************/
	/* Backward Sampling */
	/*********************/
        
	// Sample end value from distribution at end point
	    
	uniform_real_distribution<double> unif(0.0, 1.0);
	double dummyUnif = unif(generator);
	int dummyIndex = -1;

	vector<double> marginalPrey(N+1); for (int i = 0; i < N+1; ++i) for (int j = 0; j < N+1; ++j) marginalPrey[i] += ffArray[startingPath.size()-1][i][j];

	double dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += marginalPrey[dummyIndex];
   	}
    startingPath[startingPath.size()-1].X1 = dummyIndex;
	
    double normConst = 0;
	vector<double> marginalPred(N+1); 
	for (int j = 0; j < N+1; ++j) {
		marginalPred[j] = ffArray[startingPath.size()-1][dummyIndex][j];
		normConst += marginalPred[j];
	}
	for (int j = 0; j < N+1; ++j) marginalPred[j] *= 1.0/normConst; 

	dummyUnif = unif(generator);
	dummyIndex = -1;
	dummySum = 0;
	while(dummySum < dummyUnif) {
		dummyIndex++;
      	dummySum += marginalPred[dummyIndex];
   	}
    startingPath[startingPath.size()-1].X2 = dummyIndex;

    // Loop
    
    for (int i = finalSize - 2; i >= 0; --i){

		vector<LVstate> candidates;
      	vector<double> probsPrev;
      	double sum = 0.0;

      	// Start with virtual

      	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2});

      	probsPrev.push_back(- K * weighting[startingPath[i+1].X1][startingPath[i+1].X2] * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2] );

        sum = sum + probsPrev.back();

        // Prey born

        if(startingPath[i+1].X1>0){

        	candidates.push_back({0,startingPath[i+1].X1-1,startingPath[i+1].X2});

        	probsPrev.push_back( alpha * (startingPath[i+1].X1-1) * ffArray[i][startingPath[i+1].X1-1][startingPath[i+1].X2] );

        	sum = sum + probsPrev.back();

        }

        // Prey death

        if(startingPath[i+1].X1<N){

        	candidates.push_back({0,startingPath[i+1].X1+1,startingPath[i+1].X2});

        	probsPrev.push_back( beta * (startingPath[i+1].X1+1) * startingPath[i+1].X2 * ffArray[i][startingPath[i+1].X1+1][startingPath[i+1].X2] );

        	sum = sum + probsPrev.back();

        }
      
       	// Pred birth

        if(startingPath[i+1].X2>0){

        	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2-1});

        	probsPrev.push_back( delta * startingPath[i+1].X1 * (startingPath[i+1].X2-1) * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2-1] );

        	sum = sum + probsPrev.back();

        }

  		// Pred death

        if(startingPath[i+1].X2<N){

        	candidates.push_back({0,startingPath[i+1].X1,startingPath[i+1].X2+1});

        	probsPrev.push_back( gamma * (startingPath[i+1].X2+1) * ffArray[i][startingPath[i+1].X1][startingPath[i+1].X2+1] );

        	sum = sum + probsPrev.back();

        }

        // Normalise and sample

        for (int j = 0; j < probsPrev.size(); ++j) probsPrev[j] = probsPrev[j] / sum;

		dummyUnif = unif(generator);
		dummyIndex = -1;
		dummySum = 0;
		while(dummySum < dummyUnif) {
			dummyIndex++;
      		dummySum += probsPrev[dummyIndex];
   		}
    	startingPath[i].X1 = candidates[dummyIndex].X1;
    	startingPath[i].X2 = candidates[dummyIndex].X2;

    }

    // Remove virtual jumps and return
    counter = 1;
    while (counter < startingPath.size()) {
    
    	if(startingPath[counter].X1 == startingPath[counter-1].X1 && startingPath[counter].X2 == startingPath[counter-1].X2) {
    		startingPath.erase(startingPath.begin() + counter);
    	} else{
    		counter++;
    	}

    }
		
}


/*******************/
/* MCMC Procedures */
/*******************/

void do_mcmc_Unif(vector<vector<double>> observations, int N, double obsSd, int mcmcIter, double scaling, mt19937 generator, string filename, int chainNumber)
{
	// Create starting trajectory

	vector<LVstate> startingPath; 
	mcmc_startingPathVector(startingPath, observations, N, generator);

	// Prior for death rate, birth assumed known

	double priorShape = 1.0, priorRate = pow(10.0,-3.0);

	// Build observation weight matrices

	vector<vector<vector<double>>> obsWeightsMat(observations.size());
	for ( int i = 0 ; i < observations.size() ; i++ ){
		obsWeightsMat[i].resize(N+1);
		for (int j = 0; j < N+1; ++j)
		{
			obsWeightsMat[i][j].resize(N+1);
		}
	}

	for (int k = 0; k < observations.size(); ++k) {

		for (int i = 0; i < N+1; ++i){

			for (int j = 0; j < N+1; ++j){

				obsWeightsMat[k][i][j] =  exp( - ipow((observations[k][1]-i)/obsSd,2) / 2.0 ) * exp( - ipow((observations[k][2]-j)/obsSd,2) / 2.0 );

			}

		}

	}

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "alpha, beta, delta, gamma \n";

	// ITERATE ALGORITHM

	for (int i = 0; i < mcmcIter; ++i)
	{

		// Rates update

		double shape = priorShape; 
		double rate = priorRate;
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 + 1) shape++;
			if (startingPath[j+1].X1 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1;
		}

		gamma_distribution<double> gammaD(shape,1.0/rate);
		double alpha = gammaD(generator);

      	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double beta = gammaD(generator);
   
		shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 + 1) shape++;
			if (startingPath[j+1].X2 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double delta = gammaD(generator);
 
    	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double gamma = gammaD(generator);
  
       	// Save parameters and update trajectory

		myfileTrace << alpha << ", " << beta << ", " << delta << ", " << gamma <<"\n";

        // Update trajectory

		mcmc_newPathUnif(startingPath, observations, N, obsSd, alpha, beta, delta, gamma, scaling, generator, obsWeightsMat);

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
	// myfile << "Time, X1 , X2 \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].X1 << ", " << startingPath[i].X2 << "\n";
	// }
	// myfile.close();

}

void do_mcmc_DepThinning(vector<vector<double>> observations, int N, double obsSd, int mcmcIter, double K, mt19937 generator, string filename, int chainNumber)
{

	// Create starting trajectory

	vector<LVstate> startingPath; 
	mcmc_startingPathVector(startingPath, observations, N, generator);

	// Prior for death rate, birth assumed known

	double priorShape = 1.0, priorRate = pow(10.0,-3.0);

	// Build observation weight matrices

	vector<vector<vector<double>>> obsWeightsMat(observations.size());
	for ( int i = 0 ; i < observations.size() ; i++ ){
		obsWeightsMat[i].resize(N+1);
		for (int j = 0; j < N+1; ++j)
		{
			obsWeightsMat[i][j].resize(N+1);
		}
	}

	for (int k = 0; k < observations.size(); ++k) {

		for (int i = 0; i < N+1; ++i){

			for (int j = 0; j < N+1; ++j){

				obsWeightsMat[k][i][j] =  exp( - ipow((observations[k][1]-i)/obsSd,2) / 2.0 ) * exp( - ipow((observations[k][2]-j)/obsSd,2) / 2.0 );

			}

		}

	}

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "alpha, beta, delta, gamma \n";

	// ITERATE ALGORITHM

	for (int i = 0; i < mcmcIter; ++i)
	{

		// Rates update

		double shape = priorShape; 
		double rate = priorRate;
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 + 1) shape++;
			if (startingPath[j+1].X1 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1;
		}

		gamma_distribution<double> gammaD(shape,1.0/rate);
		double alpha = gammaD(generator);

      	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double beta = gammaD(generator);
   
		shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 + 1) shape++;
			if (startingPath[j+1].X2 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double delta = gammaD(generator);
 
    	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double gamma = gammaD(generator);
  
       	// Save parameters and update trajectory

		myfileTrace << alpha << ", " << beta << ", " << delta << ", " << gamma <<"\n";

        // Update trajectory

		mcmc_newPathDepThinning(startingPath, observations, N, obsSd, alpha, beta, delta, gamma, K, generator, obsWeightsMat);

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
	// myfile << "Time, X1 , X2 \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].X1 << ", " << startingPath[i].X2 << "\n";
	// }
	// myfile.close();

}

void do_mcmc_AuxVar(vector<vector<double>> observations, int N, double obsSd, int mcmcIter, double K, int lag,mt19937 generator, string filename, int chainNumber)
{

	
	// Create starting trajectory

	vector<LVstate> startingPath; 
	mcmc_startingPathVector(startingPath, observations, N, generator);

	// Prior for death rate, birth assumed known

	double priorShape = 1.0, priorRate = pow(10.0,-3.0);

	// Build observation weight matrices

	vector<vector<vector<double>>> obsWeightsMat(observations.size());
	for ( int i = 0 ; i < observations.size() ; i++ ){
		obsWeightsMat[i].resize(N+1);
		for (int j = 0; j < N+1; ++j)
		{
			obsWeightsMat[i][j].resize(N+1);
		}
	}

	for (int k = 0; k < observations.size(); ++k) {

		for (int i = 0; i < N+1; ++i){

			for (int j = 0; j < N+1; ++j){

				obsWeightsMat[k][i][j] =  exp( - ipow((observations[k][1]-i)/obsSd,2) / 2.0 ) * exp( - ipow((observations[k][2]-j)/obsSd,2) / 2.0 );

			}

		}

	}

	// Open file to store output chain

	ofstream myfileTrace;
	myfileTrace.open (filename);
	myfileTrace << "alpha, beta, delta, gamma \n";

	// ITERATE ALGORITHM

	for (int i = 0; i < min(100,mcmcIter); ++i)
	{

		// Rates update

		double shape = priorShape; 
		double rate = priorRate;
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 + 1) shape++;
			if (startingPath[j+1].X1 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1;
		}

		gamma_distribution<double> gammaD(shape,1.0/rate);
		double alpha = gammaD(generator);

      	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double beta = gammaD(generator);
   
		shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 + 1) shape++;
			if (startingPath[j+1].X2 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double delta = gammaD(generator);
 
    	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double gamma = gammaD(generator);
  
       	// Save parameters and update trajectory

		myfileTrace << alpha << ", " << beta << ", " << delta << ", " << gamma <<"\n";

        // Update trajectory
		
		mcmc_newPathDepThinning(startingPath, observations, N, obsSd, alpha, beta, delta, gamma, K, generator, obsWeightsMat);

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

	for (int i = min(100,mcmcIter); i < mcmcIter; ++i)
	{

		// Rates update

		double shape = priorShape; 
		double rate = priorRate;
		
		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 + 1) shape++;
			if (startingPath[j+1].X1 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1;
		}

		gamma_distribution<double> gammaD(shape,1.0/rate);
		double alpha = gammaD(generator);

      	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X1 == startingPath[j].X1 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double beta = gammaD(generator);
   
		shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 + 1) shape++;
			if (startingPath[j+1].X2 < N) rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X1 * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double delta = gammaD(generator);
 
    	shape = priorShape; 
		rate = priorRate;

		for (int j = 0; j < startingPath.size()-1; ++j)
		{
			if (startingPath[j+1].X2 == startingPath[j].X2 - 1) shape++;
			rate = rate + (startingPath[j+1].time - startingPath[j].time) * startingPath[j].X2;
		}

		gammaD = gamma_distribution<double>(shape,1.0/rate);
		double gamma = gammaD(generator);
  
       	// Save parameters and update trajectory

		myfileTrace << alpha << ", " << beta << ", " << delta << ", " << gamma <<"\n";

        // Update trajectory

		mcmc_newPathAuxVar(startingPath, observations, N, obsSd, alpha, beta, delta, gamma, K, lag, generator, obsWeightsMat);

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
	// myfile << "Time, X1 , X2 \n";
	// for (int i = 0; i < startingPath.size(); ++i)
	// {
	// 		myfile << startingPath[i].time << ", " << startingPath[i].X1 << ", " << startingPath[i].X2 << "\n";
	// }
	// myfile.close();

}