#include "Matrix.h"
#pragma once
class NeuronLayer
{
public:
	int size,presize;
	Matrix w;
	Matrix bias;
	double val[2005];
	Matrix origin; //û��ʩ�ӳ弤������ԭ���ֵ
	double preval[2005];
	Matrix q; //the inaccuracy
	NeuronLayer();
	NeuronLayer(int size,int lastsize,double preval[]);
	~NeuronLayer();
	void output();
};

