#include "PoolingLayer.h"
#include <iostream>


PoolingLayer::PoolingLayer(int num, int sizex, int sizey, int presizex, int presizey, Matrix preval[])
{
	this->num = num; //¾í»ý²ãÊý
	this->sizex = sizex;
	this->sizey = sizey;
	this->presizex = presizex;
	this->presizey = presizey;
	for (int i = 0; i < num; i++)
	{
		this->bias[i] = 0.0;
		this->preval[i] = &preval[i];
		this->val[i] =new Matrix(presizey / sizey, presizex / sizex);
		this->q[i] =new Matrix(presizey / sizey, presizex / sizex);
	}
}

PoolingLayer::PoolingLayer()
{
	this->num = 0;
	this->sizex = 0;
	this->sizey = 0;
	this->presizex = 0;
	this->presizey = 0;
	for (int i = 0; i < num; i++)
	{

	}
}


PoolingLayer::~PoolingLayer()
{
}

//use mean-pooling
void PoolingLayer::output()
{
	for (int i = 0; i < this->val[0]->row; i++)
	{
		for (int j = 0; j < this->val[0]->column; j++)
		{
			for (int k = 0; k < this->num; k++)
			{
				double cnt = 0.0;
				int flag = 0;
				for (int l = 0; l < this->sizey; l++)
				{
					for (int m = 0; m < this->sizex; m++)
					{
						//std::cout << i << ' ' << j;
						/*
						if (flag == 0)
						{
							cnt = this->preval[k]->M[i * this->sizey + l][j * this->sizex + m];
							flag = 1;
						}
						else
						{
							if (this->preval[k]->M[i*this->sizey+l][j*this->sizex+m] > cnt)
								cnt = this->preval[k]->M[i * this->sizey + l][j * this->sizex + m];
						}
						*/
						cnt += this->preval[k]->M[i * this->sizey + l][j * this->sizex + m];
					}
				}
				this->val[k]->M[i][j] = cnt/(this->sizex*this->sizey);
			}
		}
	}
}
