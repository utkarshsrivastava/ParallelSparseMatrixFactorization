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
#include "DataRetriever.h"

using namespace std;

solution::solution() {
}

solution::solution(DataRetriever *DRobj, double & sparsityThreshold, int & Q, int  & N, int & K) {
	_DRobj = DRobj;
	cout<<"Initial W: " << '\n';
	_Q = Q; _N = N; _K = K;
	 _W.resize(_Q*_K); // Set the size of _W
        _C.resize(_K*_N); // Set the size of _C
	int interator_1, interator_2;
	int i,j,m,n; 
	for (i = 0; i < _Q; i++) {
		interator = i * _K;
		for (j = 0; j < _K; j++) {
			double randt = ((double) rand() / (RAND_MAX));
			if (randt < sparsityThreshold) {
				_W[interator++]= 5 * ((double) rand() / (RAND_MAX));
			}else{
			 _W[interator++] = 0.0;
			}
		//	cout << _W[i*_K+j]<< " ";
		}
		//cout << '\n';
	}
	cout<<"Initial C: "<<endl;
	for (n = 0; n < _K; n++) {
		interator_2 = n*_N;
		for (m = 0; m < _N; m++) {
			_C[interator_2++] = 2 * ((double) rand() / (RAND_MAX)) - 1;
			//cout << _C[i*_N+j]<< " ";
		}
		//cout << '\n';
	}

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
unsigned long solution::objectiveFunction (const double & lambda, const double & mu, const double & gamma){
	unsigned long objectiveFunction(0.0);
	unsigned long sumLog (0.0); unsigned long sumWl1(0.0);  unsigned long sumWl2(0.0); unsigned long sumC(0.0);
	double WiCj;
	int obs_ij;		//observation data
	double logesticFunc;	// = 1/1+exp(...)
	int iter_1;		//i*_K+k
	int iter_2;		//j+k*_N
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			 obs_ij=_DRobj->observation(i,j);
			if (obs_ij != -1){
				WiCj=0.0;
				iter_1=i*_K;
				iter_2=j;
				for (int k=0; k < _K; k++){
					WiCj += _W[iter_1++]*_C[iter_2];
					iter_2+=_N;
				}
				logesticFunc=1/(1+exp(-WiCj));
				sumLog += obs_ij*log(0.01+logesticFunc) + (1-obs_ij)*log(1.01-logesticFunc);
			}
		}
		unsigned long sumTemp(0.0);
		iter_1=i*_K;
		for (int k = 0; k < _K; k++){

			sumWl1 += (_W[iter_1] < 0) ? -_W[iter_1] : _W[iter_1];
			sumTemp += _W[iter_1]*_W[iter_1];
			iter_1++;
		}
		sumWl2 += sqrt(sumTemp);
	}
	for (int k = 0; k < _K; k++){
		iter_1=k*_N;
		for (int j = 0; j < _N; j++){
			sumC += _C[iter_1]*_C[iter_1];
			iter_1++;
		}
	}
	std::cout << "SumLog: "<< sumLog <<"SumWl1: "<<sumWl1 << "SumWl2: "<<sumWl2 <<"Sumc: " << sumC << std::endl; 
	objectiveFunction = -sumLog + lambda*sumWl1 + mu*sumWl2 + gamma*sqrt(sumC);
	return objectiveFunction;
}
double solution::objectiveFunctionW (double  & lambda, double & mu){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumWl1(0.0);  double sumWl2(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (_DRobj->observation(i,j) != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i*_K+k]*_C[j+k*_N];
				}
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),_DRobj.observation(i,j))*pow(1-(1/(1+exp(-WiCj))),(1-_DRobj.observation(i,j))));
				sumLog += _DRobj->observation(i,j)* log(1/(1+exp(-WiCj))) + (1-_DRobj->observation(i,j))*log(1-1/(1+exp(-WiCj)));
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
			if (_DRobj->observation(i,j) != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i*_K+k]*_C[j+k*_N];
				}
				sumLog+= (_DRobj->observation(i,j)* log(1/(1+exp(-WiCj)))+(1-_DRobj->observation(i,j))*log(1-(1/(1+exp(-WiCj)))));
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
	double logisticFun;
	int obs_index_ij;
	for (int j = 0; j < _N; j++){
		obs_index_ij = _DRobj->observation(index_i,j);
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[index_i*_K+k]*_C[j+k*_N];
		}
		logisticFun = (1/(1+exp(-WikCkj)));
		if (obs_index_ij != -1){
			yipi[j] =obs_index_ij - logisticFun;
		}else{
			yipi[j] =(rand()%2) - logisticFun;
		}

		sum += _C[index_k*_N+j]*yipi[j];
	}
	//cout << "Sum->W: " << sum << endl;

	delF = -sum + mu*_W[index_i*_K+index_k];
	cout << "DeltaF_w: " << delF << endl;
	_W[index_i*_K+index_k] =(_W[index_i*_K+index_k]-stepSize*delF) < 0 ? 0 : (_W[index_i*_K+index_k]-stepSize*delF);

	return _W[index_i*_K+index_k];
}
double solution::updateCij(const int & index_k, const int & index_j, const double & gamma,const double & stepSize){
	double delF(0.0);
	double WikCkj (0.0);
	double yipi[_Q];
	double sum (0.0);
	int obs_index_j
	for (int j = 0; j < _Q; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[j*_K+k]*_C[index_j+k*_N];
			}
		if (_DRobj->observation(j,index_j) != -1){
			yipi[j] = _DRobj->observation(j,index_j)-(1/(1+exp(-WikCkj)));
		}else{
			yipi[j] = (rand()%2)-(1/(1+exp(-WikCkj)));
		}
		sum += _W[j*_K+index_k]*yipi[j];
	}
	//cout << "Sum->C: " << sum << endl;
	delF = -sum + gamma*_C[index_k*_N+index_j];
	cout << "DeltaF_w: " << delF << endl;

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
		if (_DRobj.observation(index_i,j) != -1){
			yipi[j] =_DRobj.observation(index_i,j) - (1/(1+exp(-WikCkj)));
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
