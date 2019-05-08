#include "NeuronLayer.h"
#include"ConvelutionLayer.h"
#include"PoolingLayer.h"
#pragma once
class network
{
public:
	int size;
	int convSize;
	int x[40005][28][28];
	int y[40005];
	double studyrate;
	double lamda;
	NeuronLayer layers[10];
	ConvelutionLayer clayers[2];
	PoolingLayer players[2];
	network(int size,int *l,double x[],double y[]);
	network();
	void BP(int testsize);
	void fullBP(int testsize);
	void init();
	void test();
	~network();
	void fullBPTest(int testsize);
};

