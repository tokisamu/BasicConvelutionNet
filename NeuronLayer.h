#include "Matrix.h"
#pragma once
class NeuronLayer
{
public:
	int size,presize;
	Matrix w;
	Matrix bias;
	double val[2005];
	Matrix origin; //没有施加冲激函数的原输出值
	double preval[2005];
	Matrix q; //the inaccuracy
	NeuronLayer();
	NeuronLayer(int size,int lastsize,double preval[]);
	~NeuronLayer();
	void output();
};

