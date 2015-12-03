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
#include <time.h>

using namespace std;




int main(int argc, char* args[]) {
	// *****************************
	// Parallel Setup
	int NUM_PROCS, NUM_THREADS;
	//NUM_PROCS = omp_get_num_procs(); // number of processes requested
	NUM_PROCS = 4;
	//NUM_PROCS = (int)(args[0]);
	NUM_THREADS = NUM_PROCS;
	omp_set_num_threads(NUM_THREADS); 
	
	// Parameters Setting
	double mu (10.0); // This represent parameter of l2 regularization
	double init_mu = mu;
	double lambda (1.0); // This represent parameter of l1 regularization
	double gamma (10.0); // This represent parameter of Frobenius norm
	double init_gamma = gamma;
	static int numOfIteration (2); // Number of Iteration of main optimization loop
	double sparcityParameter = 0.05;
	double stepSize (1000); // The stepsize in the optimization process 
	double init_stepSize = stepSize;
	// number of processes requested
	
	// Generating Random Indexes for whole otimization
	
	//*****************************
	int ** data; // full data
	int ** observation; // observation = data-missing values
	int Q(2000),N(1000), K(10);
	
	std::vector<double> W;
	std::vector<double> C;

	dataPreparation dataObj(Q,N,K);
	dataObj.initialization();
	data = dataObj.getMainData();
	observation = dataObj.getObservation();
	cout<<"****************"<<" OBSERVATION "<<endl;
	int index;
	double timeStart,timeEnd;
	int num_rand = numOfIteration*NUM_PROCS;
	 // Generating Random Indexes for whole otimization
	std::vector<double> I_index(num_rand);
        std::vector<double> K_index(num_rand);
        std::vector<double> J_index(num_rand);
	
	timeStart = omp_get_wtime();
	
	#pragma omp parallel shared(I_index,K_index,J_index,numOfIteration,NUM_PROCS)
	{
	srand(time(NULL)+omp_get_thread_num()); 
	#pragma omp for private(index) 	
        for (index = 0; index < num_rand; index++){
                I_index[index] =(rand() % (Q));
                K_index[index] = (rand() % (K));
                J_index[index] = (rand() % (N));
        }
	}

	solution solutionObj(observation,  NUM_PROCS,sparcityParameter, Q, N, K);
	W = solutionObj.getW();
	C = solutionObj.getC();	
	double timeS_Obj, timeE_Obj;
	timeS_Obj = omp_get_wtime();
	cout << "Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;
	timeE_Obj = omp_get_wtime();
	std::cout << "Time Objective Function: " << timeE_Obj - timeS_Obj <<  std::endl;
	cout<<"****************"<<" Optimization Started..."<<endl;
	int currentIteration = 0;
	int counter = 1;
	while (currentIteration < numOfIteration){
		//NUM_PROCS = 1;
		omp_set_num_threads(NUM_PROCS);
		int tid;
		#pragma omp parallel shared(W,C)
		{ 
		int  procIterator;
		int iW, kW, kC, jC;
		//#pragma omp for schedule(dynamic) private(iW, kW, kC, jC,procIterator)
		//for (procIterator = 0; procIterator < NUM_PROCS; procIterator++){
			tid = omp_get_thread_num();
			//std::cout<< "ThreadID: " << tid << " rand: " << I_index[currentIteration*NUM_PROCS+tid]<< std::endl;
			iW = I_index[currentIteration*NUM_PROCS+tid];
			kW = K_index[currentIteration*NUM_PROCS+tid];
			jC = J_index[currentIteration*NUM_PROCS+tid];
			//iW =(rand() % (Q));
			//kW = (rand() % (K));
			//printf("process id: %d Iterator: %d iW: %d kW: %d\n",tid,procIterator,iW,kW); 
			
			W[iW*K+kW] = solutionObj.updateWij(iW,kW,mu,stepSize);
			//W[iW]  = solutionObj.updateRowWij(iW,mu,stepSize);
			kC = kW; jC = (rand() % (N));
			C[kC*N+jC] = solutionObj.updateCij(kC,jC,gamma,stepSize);
			//C[jC] = solutionObj.updateColumnCij(jC,gamma,stepSize);
			//cout << "C[i][j]: " << C[kC][jC] << endl;
		//}
		}
		cout << "Iteration: " << currentIteration << ", Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;
		stepSize /= 1.2;
		mu = mu/1.2;
		gamma = gamma/1.2;
		if(currentIteration % 10 == 0 && currentIteration != 0 ){
			stepSize = (double)(1/++counter) * init_stepSize;;
			mu = (double)(1/++counter) * init_mu;;
			gamma = (double)(1/++counter) * init_gamma;
		} 
		currentIteration++;
	}

	timeEnd =omp_get_wtime();
	for (int i=0; i < Q; i++){
		for (int j =0; j < K; j++){
		//	cout << W[i*K+j] << " ";
		}
		//cout << '\n';
	}
	std::cout << "Num of Procs: " <<  NUM_PROCS  <<" Total time: " << timeEnd-timeStart << std::endl;
return -1;
}



