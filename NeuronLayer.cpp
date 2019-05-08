#include "NeuronLayer.h"
#include<math.h>
#include <stdlib.h>
#include <iostream>



NeuronLayer::NeuronLayer()
{
	this->size = 1;
}

NeuronLayer::NeuronLayer(int size,int lastsize,double preval[])
{
	this->size = size;
	this->presize = lastsize;
	this->w = Matrix(presize, size);
	this->w.init(size);
	this->bias = Matrix(1,size);
	this->origin =Matrix(1, size);
	for (int i = 0; i < sizeof(preval)/sizeof(double); i++)
	{
		this->preval[i] = preval[i];
	}
	this->q = Matrix(1, size);
}

NeuronLayer::~NeuronLayer()
{
}

//use sigmoid to calculate output
void NeuronLayer::output()
{
	Matrix temp(1, presize);
	Matrix ans(1, size);
	for (int i = 0; i < presize; i++) //将输入转换为行矩阵
		temp.M[0][i] = this->preval[i];
	ans = this->w.multi(temp, this->w) + this->bias; //z=w*x+b
	for (int i = 0; i < this->size; i++) //sigmoid
	{
		this->origin.M[0][i] = ans.M[0][i];
		//std::cout << ans.M[0][i] << ' ';
		ans.M[0][i] = 1.0/(1.0+exp(-ans.M[0][i]));
		this->val[i] = ans.M[0][i];
	}
}
