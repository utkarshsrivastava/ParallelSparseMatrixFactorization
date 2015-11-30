/*
 * mainDriver.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: afarasat
 */
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <cstdlib>
#include "solution.h"
#include "dataPreparation.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;




int main(int argc, char* args[]) {
	// *****************************
	// Parallel Setup
	int NUM_PROCS, NUM_THREADS;
	NUM_PROCS = omp_get_num_procs();
	//NUM_PROCS = (int)(args[0]);
	NUM_THREADS = NUM_PROCS;
	omp_set_num_threads(NUM_THREADS); 
	
	// Parameters Setting
	double mu (10.0); // This represent parameter of l2 regularization
	double init_mu = mu;
	double lambda (1.0); // This represent parameter of l1 regularization
	double gamma (10.0); // This represent parameter of Frobenius norm
	double init_gamma = gamma;
	static int numOfIteration (200); // Number of Iteration of main optimization loop
	double sparcityParameter = 0.5;
	double stepSize (10); // The stepsize in the optimization process 
	double init_stepSize = stepSize;
	int numOfNodes;     // number of processes requested
	std::vector<double> I_index;
	std::vector<double> K_index;
 	std::vector<double> J_index;
	
	// Generating Random Indexes for whole otimization
	
	//*****************************
	int ** data; // full data
	int ** observation; // observation = data-missing values
	int Q(5),N(4), K(2);
	// Generating Random Indexes for whole otimization
	std::vector<double> W;
	std::vector<double> C;

	dataPreparation dataObj(Q,N,K);
	dataObj.initialization();
	data = dataObj.getMainData();
	observation = dataObj.getObservation();
	cout<<"****************"<<" OBSERVATION "<<endl;
//for (int i=0; i < Q; i++){
//for (int j =0; j < N; j++){
		//	cout << observation[i][j] << " ";
//		}
		//cout << '\n';
//	}
	int index;
	double timeStart,timeEnd;
	//pragma omp for private(index) 	
        for (index = 0; index < numOfIteration*numOfNodes; index++){
                I_index[index] =(rand() % (Q-1));
                K_index[index] = (rand() % (K-1));
                J_index[index] = (rand() % (N-1));
        }

	timeStart = omp_get_wtime();
	solution solutionObj(observation, sparcityParameter, Q, N, K);
	W = solutionObj.getW();
	C = solutionObj.getC();

	cout << "Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;

	cout<<"****************"<<" Optimization "<<endl;
	int currentIteration = 0;
	int counter = 1;
	while (currentIteration < numOfIteration){
//		pragma omp parallel
		{
		int tid; 
		int iW, kW, kC, jC;
//		pragma omp for schedule(dynamic) private(iW, kW, kC, jC)
		for (int  procIterator = 0; procIterator < NUM_PROCS; procIterator++){
			tid = omp_get_thread_num();
			iW =(rand() % (Q-1));
			kW = (rand() % (K-1));
			//printf("process id: %d Iterator: %d iW: %d kW: %d\n",tid,procIterator,iW,kW); 
			W[iW*K+kW] = solutionObj.updateWij(iW,kW,mu,stepSize);
			//W[iW]  = solutionObj.updateRowWij(iW,mu,stepSize);
			kC = (rand() % (K-1)); jC = (rand() % (N-1));
			C[kC*N+jC] = solutionObj.updateCij(kC,jC,gamma,stepSize);
			//C[jC] = solutionObj.updateColumnCij(jC,gamma,stepSize);
			//cout << "C[i][j]: " << C[kC][jC] << endl;
		}
		}
		cout << "Iteration: " << currentIteration << ", Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;
		stepSize /= 1.5;
		mu = mu/1.5;
		gamma = gamma/1.5;
		if(currentIteration % 100 == 0 && currentIteration != 0 ){
			stepSize = (double)(1/++counter) * init_stepSize;;
			mu = (double)(1/++counter) * init_mu;;
			gamma = (double)(1/++counter) * init_gamma;
		} 
		currentIteration++;
	}

	timeEnd =omp_get_wtime();
//	for (int i=0; i < Q; i++){
//		for (int j =0; j < K; j++){
		//	cout << W[i*K+j] << " ";
//		}
		//cout << '\n';
//	}
	std::cout << "Num of Procs: " <<  NUM_PROCS  <<" Total time: " << timeEnd-timeStart << std::endl;
return -1;
}



