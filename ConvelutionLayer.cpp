#include "ConvelutionLayer.h"
#include"Matrix.h"
#include<math.h>
#include <iostream>

double sigmoid(double x)
{
	return 1.0 / (1.0 + exp(-x));
}

ConvelutionLayer::ConvelutionLayer(int num, int distance, int sizex, int sizey, int presizex, int presizey, Matrix preval)
{
	this->distance = distance; //卷积间距
	this->num = num; //卷积层数
	this->sizex = sizex;
	this->sizey = sizey;
	this->presizex = presizex;
	this->presizey = presizey;
	this->preval = preval;
	for (int i = 0; i < num; i++)
	{
		this->w[i] =new Matrix(sizey, sizex);
		this->w[i]->init(sizex * sizey);
		this->val[i] =new Matrix(presizey - sizey + 1, presizex - sizex + 1);
		this->q[i] =new Matrix(presizey - sizey + 1, presizex - sizex + 1);
		bias[i] = 0.0;
	}
	
}

ConvelutionLayer::ConvelutionLayer()
{
}

ConvelutionLayer::~ConvelutionLayer()
{
}

//use siamoid to calculate output
void ConvelutionLayer::output()
{
	double cnt = 0.0;
	for (int i = 0; i < this->presizey - this->sizey + 1; i+=distance)
	{
		for (int j = 0; j < this->presizex - this->sizex + 1; j+=distance)
		{
			for (int nn = 0; nn < this->num; nn++)
			{
				double cnt = 0.0;
				for (int k = 0; k < this->sizey; k++)
				{
					for (int l = 0; l < this->sizex; l++)
					{
						cnt += preval.M[i * this->distance + k][j * this->distance + l]* this->w[nn]->M[k][l];
					}
				}
				cnt += this->bias[nn];
				//std::cout << i << ' ' << cnt<< ' ' << sigmoid( cnt )<< '\n';
				this->val[nn]->M[i][j] = sigmoid(cnt);
			}
		}
	}
	//this->val[0]->print();
}
