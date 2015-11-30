/*
 * solution.h
 *
 *  Created on: Oct 10, 2015
 *      Author: afarasat
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_
#include <math.h>
#include<iostream>
#include<vector>
#include "DataRetriever.h"

class solution {
public:
	solution();
	solution(DataRetriever * ,double &,int &,int &,int &);
	std::vector<double> getW();
	std::vector<double> getC();
	unsigned long objectiveFunction(const double &, const double &,const double &);
	double objectiveFunctionW(double &, double &);
	double objectiveFunctionC(double &);
	double updateWij(const int &,const int &,const double &,const double &);
	double updateCij(const int &,const int &, const double &,const double &);
	double * updateRowWij(int &,double &,double &);
	double * updateColumnCij(int &,double &,double &);

	virtual ~solution();
private:
	int _Q, _N, _K;
	std::vector<double> _W;
	std::vector<double> _C;
	DataRetriever *_DRobj;
};

#endif /* SOLUTION_H_ */
