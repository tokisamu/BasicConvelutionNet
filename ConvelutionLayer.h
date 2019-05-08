#include "Matrix.h"
#pragma once
class ConvelutionLayer
{
public:
	int sizex,sizey,presizex,presizey,num,distance;
	Matrix *w[10];
	double bias[10];
	Matrix *val[10];
	Matrix preval;
	Matrix *q[10]; //the inaccuracy
	ConvelutionLayer(int num,int distance,int sizex,int sizey,int presizex,int presizey,Matrix preval);
	ConvelutionLayer();
	~ConvelutionLayer();
	void output();
};

