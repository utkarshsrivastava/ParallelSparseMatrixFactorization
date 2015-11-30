/*
 * dataPreparation.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: afarasat
 */
#include <iostream>
#include <math.h>
#include <cstdlib>
#include "dataPreparation.h"
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <ctime>

using namespace std;






dataPreparation::dataPreparation(int a, int b, int c) {
	_Q = a, _N = b, _K =c;
	_dataThresold = 0.5;
	_obsThresold = 0.5;
	_Y = new int *[_Q];
	for (int k = 0; k < _Q; k++){
		_Y[k] = new int[_N];
	}
	_observation =  new int *[_Q];
	for (int k = 0; k < _Q; k++){
		_observation[k] = new int[_N];
	}





}

dataPreparation::~dataPreparation() {
	// TODO Auto-generated destructor stub
}
 void dataPreparation::initialization() {

    std::vector<double> W(_Q*_K);
    std::vector<double> C(_K*_N);

	/*vector<vector<double> >  C1( _K, std::vector<double> ( _N, 0 ) );
	vector<vector<double> > W1( _Q, std::vector<double> ( _K, 0 ) );
	vector<vector<int> > Y ( _Q, std::vector<int> ( _N, 0 ) );*/



	srand(time(NULL));
	cout<< "W: "<< '\n';
	for (int i = 0; i < _Q; i++) {
		for (int j = 0; j < _K; j++) {
			W[i*_K+j]= 0.0;
			double randt = ((double) rand() / (RAND_MAX));
			if (randt < _dataThresold) {
				W[i*_K+j] = 5 * ((double) rand() / (RAND_MAX));
			}
		//	cout << W[i*_K+j] << " ";
		}
		//cout << '\n';
	}
	cout<<"C: "<<endl;
	for (int i = 0; i < _K; i++) {
		for (int j = 0; j < _N; j++) {
			C[i*_N+j] = 2 * ((double) rand() / (RAND_MAX)) - 1;
			//cout << C[i*_N+j]<< " ";
		}
		//cout << '\n';
	}
	//cout<< "Main Data: "<< '\n';
	for (int i = 0; i < _Q; i++) {
		for (int j = 0; j < _N; j++) {
			double sum = 0;
			for (int k = 0; k < _K; k++) {
				sum += W[i*_K+k]*C[j+k*+_N];
			}
			if (sum > 0){
				_Y[i][j] = 1;
				_observation[i][j] = 1;
				if ((double) rand() / (RAND_MAX)<_obsThresold/5 ){
					_observation[i][j] = -1;
				}
			}else{
				_Y[i][j] = 0;
				_observation[i][j] = 0;
				if ((double) rand() / (RAND_MAX)<_obsThresold ){
					_observation[i][j] = -1;
				}

			}
			//cout<<_Y[i][j]<<" ";
		}
		//cout << endl;

	}

	//_Y =  & Y;
}

int ** dataPreparation::getMainData(){
	return  _Y;
}
int ** dataPreparation::getObservation(){
	return  _observation;
}

