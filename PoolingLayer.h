#include "Matrix.h"
#pragma once
class PoolingLayer
{
public:
	int sizex,sizey,presizex,presizey,num;
	Matrix *w[10];
	double bias[10];
	Matrix *val[10];
	Matrix *preval[10];
	Matrix *q[10]; //the inaccuracy
	PoolingLayer(int num,int sizex,int sizey,int presizex,int presizey,Matrix preval[]);
	PoolingLayer();
	~PoolingLayer();
	void output();
};

