/*
 * solution.cpp
 *
 *  Created on: Oct 10, 2015
 *      Author: afarasat
 */

#include "solution.h"
# include <iostream>
#include <cstdlib>
#include <math.h>
//#include <omp.h>
#include <time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

using namespace std;
std::vector<double> _W;
std::vector<double> _C;
int _Q, _N, _K;
solution::solution() {


}
solution::solution(int ** observation, double & sparsityThreshold, int & Q, int  & N, int & K) {
	_observation = observation;
	cout<<"Initial W: " << '\n';
	_Q = Q; _N = N; _K = K;
	double seed, randt;
        int i, j, m, n;
	_W.resize(_Q*_K);
	_C.resize(_K*_N);
int numofproc=100;// WHY ARE YOU AGAINS SETTING  NUMBER OF PROCS

//	omp_set_num_threads(omp_get_num_procs());
long KQ=_K*_Q;

//	{
	srand(time(NULL)+numofproc);	
	cilk_for (i = 0; i <KQ; i++) {
	//	{
//		cilk_for (j = 0; j < _K; j++) {
			_W[i]=0;
			randt = ((double) rand() / (RAND_MAX));
			if (randt < sparsityThreshold) {
				_W[i]= 5 * ((double) rand() / (RAND_MAX));
			}
			//cout << _W[i*_K+j]<< " ";
//		}
	//	}
		//cout << '\n';
	}
	
//	std::cout<<"Initial C: "<<endl;	
        srand(time(NULL)+numofproc);
	cilk_for (n = 0; n < _K*_N; n++) {
		//for (int m = 0; m < _N; m++) {
			_C[n]=2 * ((double) rand() / (RAND_MAX)) - 1;
			//cout << _C[i*_N+j]<< " ";
		}
		//cout << '\n';
	}
	}
	std::cout<<"Initial C: "<<endl;
}
solution::~solution() {
	// TODO Auto-generated destructor stub
}
std::vector<double> solution::getW(){
	return _W;
}
std::vector<double> solution::getC(){
	return _C;
}
double solution::objectiveFunction (const double & lambda, const double & mu, const double & gamma){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumWl1(0.0);  double sumWl2(0.0); double sumC(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (_observation[i][j] != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i*_K+k]*_C[j+k*_N];
				}
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),_observation[i][j])*pow(1-(1/(1+exp(-WiCj))),(1-_observation[i][j])));
				sumLog += _observation[i][j]* log(0.01+1/(1+exp(-WiCj))) + (1-_observation[i][j])*log(1.01-1/(1+exp(-WiCj)));
			}
		}
		double sumTemp(0.0);
		for (int k = 0; k < _K; k++){

			sumWl1 += (_W[i*_K+k] < 0) ? -_W[i*_K+k] : _W[i*_K+k];
			sumTemp += _W[i*_K+k]*_W[i*_K+k];
		}
		sumWl2 += pow(sumTemp,0.5);
	}
	for (int i = 0; i < _K; i++){
		for (int j = 0; j < _N; j++){
			sumC += _C[i*_N+j]*_C[i*_N+j];
		}
	}
	objectiveFunction = -sumLog + lambda*sumWl1 + mu*sumWl2 + gamma*sqrt(sumC);
	return objectiveFunction;
}
double solution::objectiveFunctionW (double  & lambda, double & mu){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumWl1(0.0);  double sumWl2(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (_observation[i][j] != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i*_K+k]*_C[j+k*_N];
				}
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),_observation[i][j])*pow(1-(1/(1+exp(-WiCj))),(1-_observation[i][j])));
				sumLog += _observation[i][j]* log(1/(1+exp(-WiCj))) + (1-_observation[i][j])*log(1-1/(1+exp(-WiCj)));
			}
		}
		double sumTemp(0.0);
		for (int k = 0; k < _K; k++){

			sumWl1 += (_W[i*_K+k] < 0) ? -_W[i*_K+k] : _W[i*_K+k];
			sumTemp += _W[i*_K+k]*_W[i*_K+k];
		}
		sumWl2 += pow(sumTemp,0.5);
	}

	objectiveFunction = -sumLog + lambda*sumWl1 + mu*sumWl2;
	return objectiveFunction;
}
double solution::objectiveFunctionC (double & gamma){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumC(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (_observation[i][j] != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i*_K+k]*_C[j+k*_N];
				}
				sumLog+= (_observation[i][j]* log(1/(1+exp(-WiCj)))+(1-_observation[i][j])*log(1-(1/(1+exp(-WiCj)))));
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),observation[i][j])   *    pow(1-(1/(1+exp(-WiCj))),(1-observation[i][j])));
			}
		}
	}
	for (int i = 0; i < _K; i++){
		for (int j = 0; j < _N; j++){
			sumC += _C[i*_N+j]*_C[i*_N+j];
		}
	}
	objectiveFunction = -sumLog + gamma*sqrt(sumC);
	return objectiveFunction;
}
double solution::updateWij(const int & index_i, const int  & index_k, const double & mu, const double & stepSize){
	double delF(0.0);
	double WikCkj (0.0);
	double yipi[_N];
	double sum (0.0);
	for (int j = 0; j < _N; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[index_i*_K+k]*_C[j+k*_N];
			}
		if (_observation[index_i][j] != -1){
			yipi[j] =_observation[index_i][j] - (1/(1+exp(-WikCkj)));
		}else{
			yipi[j] =(rand()%2) - (1/(1+exp(-WikCkj)));
		}

		sum += _C[index_k*_N+j]*yipi[j];
	}
	//cout << "Sum->W: " << sum << endl;

	delF = -sum + mu*_W[index_i*_K+index_k];
	//cout << "DeltaF_w: " << delF << endl;
	_W[index_i*_K+index_k] =(_W[index_i*_K+index_k]-stepSize*delF) < 0 ? 0 : (_W[index_i*_K+index_k]-stepSize*delF);

	return _W[index_i*_K+index_k];
}
double solution::updateCij(const int & index_k, const int & index_j, const double & gamma,const double & stepSize){
	double delF(0.0);
	double WikCkj (0.0);
	double yipi[_Q];
	double sum (0.0);
	for (int j = 0; j < _Q; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[j*_K+k]*_C[index_j+k*_N];
			}
		if (_observation[j][index_j] != -1){
			yipi[j] = _observation[j][index_j]-(1/(1+exp(-WikCkj)));
		}else{
			yipi[j] = (rand()%2)-(1/(1+exp(-WikCkj)));
		}
		sum += _W[j*_K+index_k]*yipi[j];
	}
	//cout << "Sum->C: " << sum << endl;
	delF = -sum + gamma*_C[index_k*_N+index_j];
	//cout << "DeltaF_w: " << delF << endl;
	
	_C[index_k*_N+index_j] = (1/(1+gamma*stepSize))*(_C[index_k*_N+index_j]-stepSize*delF);
	//_C[index_k][index_j] = _C[index_k][index_j]-stepSize*delF;

	return _C[index_k*_N+index_j];
}
/*
 *

double * solution::updateRowWij(int & index_i, double & mu,double & stepSize){
	double delF [_K];
	double WikCkj (0.0);
	double yipi[_N];
	for (int j = 0; j < _N; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[index_i][k]*_C[k][j];
			}
		if (_observation[index_i][j] != -1){
			yipi[j] =_observation[index_i][j] - (1/(1+exp(-WikCkj)));
		}else{
			yipi[j] =(rand()%2) - (1/(1+exp(-WikCkj)));
		}
	}
	double sum [_K];
	for (int k = 0; k < _K; k++){
		for (int n = 0; n < _N; n++){
			sum[k] += _C[k][n]*yipi[n];
		}
		delF[k] = -sum[k] + mu*_W[index_i][k];
		_W[index_i][k] =(_W[index_i][k]-stepSize*delF[k]) < 0 ? 0 : (_W[index_i][k]-stepSize*delF[k]);
	}

	return _W[index_i];
}

double * solution::updateColumnCij(int & index_j, double & gamma,double & stepSize){
	double delF[_K];
	double WikCkj (0.0);
	double yipi[_Q];
	for (int j = 0; j < _Q; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[j][k]*_C[k][index_j];
			}
		if (_observation[j][index_j] != -1){
			yipi[j] = _observation[j][index_j]-(1/(1+exp(-WikCkj)));
		}else{
			yipi[j] = (rand()%2)-(1/(1+exp(-WikCkj)));
		}
	}
	double sum [_K];
	for (int k = 0; k < _K; k++){
		for (int n = 0; n < _Q; n++){
			sum[k] += _W[n][k]*yipi[n];
		}
		delF[k] = -sum[k] + gamma*_C[k][index_j];
		_C[k][index_j] = (1/(1+gamma*stepSize))*(_C[k][index_j]-stepSize*delF[k]);
	}

	//_C[index_k][index_j] = _C[index_k][index_j]-stepSize*delF;

	return _C[index_j];
}
 */
