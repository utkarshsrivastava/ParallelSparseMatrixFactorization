/*
 * dataPreparation.h
 *
 *  Created on: Oct 8, 2015
 *      Author: afarasat
 */
#include <vector>
#ifndef DATAPREPARATION_H_
#define DATAPREPARATION_H_

class dataPreparation {
public:
	dataPreparation(int, int, int);
	virtual ~dataPreparation();
	void initialization();
	//std::vector<std::vector<int> > * getMainData();
	int ** getMainData();
	int ** getObservation();
private:
	int _Q,_N,_K;
	double _dataThresold;
	double _obsThresold;
	int ** _Y;
	int ** _observation;
	//std::vector<std::vector<int> > * _Y;


};

#endif /* DATAPREPARATION_H_ */
