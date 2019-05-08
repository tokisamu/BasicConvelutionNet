#include "network.h"
#include"NeuronLayer.h"
#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include <iostream>

network::network(int size,int *l, double x[], double y[])
{
	this->size = size;
	/*
	for (int i = 0; i < size; i++)
	{
		this->layers[i] = NeuronLayer(l[i], num, pre);
		layers[i].output();
		num = l[i];
		pre = layers[i].val;
	}*/
}

network::network()
{
}

double sigmoid2(double x)
{
	return 1.0 / (1.0 + exp(-x));
}

void network::BP(int minitest)
{
	int correct = 0;
	for (int tt = 0; tt < minitest; tt++)
	{
		//前向传播
		int t = rand() % 10000;
		double inp[28 * 28+10];
		int pos = 0;
		for (int i = 0; i < 28; i++) //输入矩阵
		{
			for (int j = 0; j < 28; j++)
				//for(int k=0;k<clayers->num;k++)
				inp[pos++] = x[t][i][j];
		}
		for(int i=0;i<pos;i++)
			layers[0].preval[i] = inp[i];
		for (int i = 0; i < this->size; i++) //前一层输入作为后一层输出
		{
			layers[i].output();
			if (i + 1 < this->size)
				for(int j=0;j<layers[i].size;j++)
					layers[i + 1].preval[j] = layers[i].val[j];
		}
		//for (int i = 0; i < pos; i++)
		//	std::cout << layers[0].preval[i] << ' ';
		//std::cout << '\n';
		//计算最后一层梯度
		//for (int i = 0; i < 10; i++)
		//	std::cout << layers[size - 1].val[i] << ' ';
		Matrix delta(1, 10);
		//使用交叉熵进行误差计算
		double entropy = 0.0;
		for (int i = 0; i < 10; i++)
		{
			if(fabs(y[t]-i)<0.000001)
			//entropy += (0.5 * (y[t][i] - layers[size - 1].val[i]) * (y[t][i] - layers[size - 1].val[i]));
				entropy += ( - (1) * log(1 - layers[size - 1].val[i]));
			else
				entropy += (-(1 * log(layers[size - 1].val[i])));
			//std::cout << -( log(1 + layers[size - 1].val[i])) <<' '<<  log(1 - layers[size - 1].val[i]) << '\n';
		}
		//std::cout << entropy << '\n';
		//测试输出 
		double ma = -10.0;
		int ans=-99.99;
		for (int xx = 0; xx < 10; xx++)
		{
			//std::cout << layers[1].val[xx] << ' ';
			if (layers[1].val[xx] >= ma)
			{
				ma = layers[1].val[xx];
				ans = xx;
			}
		}
		//std::cout << '\n';
		int ans2 = y[t];
		//if(tt%10==7)
			//std::cout << entropy/10 << ' '<<ans<<' '<<ans2<<'\n';
		//测试完毕
		if (fabs(ans - ans2) <= 0.001)
			correct++;
		//std::cout << ans << " " << ans2 << ' ' << entropy << '\n';
		for (int i = 0; i < 10; i++)
		{
			if (fabs(y[t] - i) < 0.000001)
			{
				//delta.M[0][i]= layers[size - 1].val[i] - y[t][i];
				delta.M[0][i] = (layers[size - 1].val[i]-1) / 10;
			}
			else
			{
				delta.M[0][i] = (layers[size - 1].val[i]) / 10;
			}
		}
		//delta.print();
		layers[size - 1].bias = layers[size - 1].bias - delta*studyrate;
		Matrix deltaw(layers[size - 2].size, layers[size - 1].size);
		Matrix tempw(layers[size - 2].size, layers[size - 1].size);
		for (int i = 0; i < layers[size - 2].size; i++)
		{
			for (int j = 0; j < layers[size - 1].size; j++)
			{
				deltaw.M[i][j] = layers[size - 2].origin.M[0][i] *delta.M[0][j];
				tempw.M[i][j] = layers[size - 1].w.M[i][j];
				//cout << n1.val.M[0][i] << " " << delta.M[0][j] << endl;
			}
		}
		//delta.print();
		//deltaw.print();
		layers[size - 1].w = layers[size - 1].w - deltaw*studyrate;
		//layers[size - 1].w.print();
		//反向传导
		for (int i = this->size - 2; i >= 0; i--)
		{
			Matrix tranDelta=delta.tran();
			//tranDelta.print();
			delta = delta.multi(tempw, tranDelta);
			//delta.print();
			//std::cout << delta.row << ' ';
			//std::cout << sigmoid2(layers[0].origin.M[0][1]) * (1 - sigmoid2(layers[0].origin.M[0][1]))<<' ';
			tranDelta.column = layers[i].size;
			tranDelta.row = 1;
			//delta.print();
			for (int j = 0; j < layers[i].size; j++)
				tranDelta.M[0][j] = delta.M[j][0]* sigmoid2(layers[i].origin.M[0][j]) * (1 - sigmoid2(layers[i].origin.M[0][j]));
			//tranDelta.print();
			delta = tranDelta;
			//delta.print();
			layers[i].bias = layers[i].bias - delta*studyrate;
			deltaw.reset();
			deltaw.column = layers[i].size;
			deltaw.row = layers[i].presize;
			tempw.row = layers[i].presize;
			tempw.column = layers[i].size;
			for (int ii = 0; ii < layers[i].presize; ii++)
			{
				for (int j = 0; j < layers[i].size; j++)
				{
					tempw.M[i][j] = layers[i].w.M[ii][j];
					deltaw.M[ii][j] = layers[i].preval[ii] * delta.M[0][j];
					//cout << n1.val.M[0][i] << " " << delta.M[0][j] << endl;
				}
			}
			//deltaw.print();
			layers[i].w = layers[i].w - deltaw*studyrate;
			//if (tt == 9) std::cout << ' ' << layers[i+1].w.M[45][7];
		}
		//layers[0].bias.print();
		//layers[1].bias.print();
	}
	std::cout << correct << '/' << minitest;
}
void network::fullBP(int testsize)
{
	int correct = 0;
	for (int tt = 0; tt < testsize; tt++) //对一个小批样本进行训练,一轮100次，样本随机
	{
		int t = rand() % 40000;
		for (int i = 0; i < 28; i++) //输入矩阵
		{
			for (int j = 0; j < 28; j++)
				//for(int k=0;k<clayers->num;k++)
					clayers[0].preval.M[i][j] = this->x[t][i][j];
		}
		//前向传播
		//卷积层
		clayers[0].output();
		//clayers[0].val[0]->print();
		//clayers[0].w[0]->print();
		//std::cout << clayers[0].bias[0] << ' ';
		for (int i = 0; i < clayers[0].num; i++)
		{
			players[0].preval[i] = clayers[0].val[i];
		}
		players[0].output();
		//players[0].val[0]->print();
		//players[0].val[1]->print();
		//players[0].val[2]->print();
		//池化层->输出层
		double poolOut[144*10]; //将输出矩阵映射到样本空间
		int pos = 0;
		for (int k = 0; k < players->num; k++)
		{
			for (int i = 0; i < 12; i++)
			{
				for (int j = 0; j < 12; j++)
				{
					poolOut[pos++] = players[0].val[k]->M[i][j];
					//std::cout << poolOut[pos - 1] << ' ';
				}
			}
		}
		for(int i=0;i<pos;i++)
			layers[0].preval[i] = poolOut[i];
		//for (int i = 0; i < pos; i++)
		//	std::cout << layers[0].preval[i] << ' ';
		//std::cout << '\n';
		/*/ 神经层单元测试
		double inpp[10005];
		pos = 0;
		for (int i = 0; i < 28; i++)
		{
			for (int j = 0; j < 28; j++)
			{
				layers[0].preval[pos++] = x[t][i][j];
			}
		}
		*/
		//隐层和输出层
		for (int i = 0; i < size; i++) //前一层输入作为后一层输出
		{
			layers[i].output();
			if (i + 1 < size)
				for (int j = 0; j < layers[i].size; j++)
					layers[i + 1].preval[j] = layers[i].val[j];
		}
		//for (int i = 0; i < pos; i++)
		//	std::cout << layers[0].preval[i] << ' ';
		//std::cout << '\n';
		//计算最后一层梯度
		//for (int i = 0; i < 10; i++)
		//	std::cout << layers[size - 1].val[i] << ' ';
		Matrix delta(1, 10);
		//使用交叉熵进行误差计算
		double entropy = 0.0;
		for (int i = 0; i < 10; i++)
		{
			if (fabs(y[t] - i) < 0.000001)
				//entropy += (0.5 * (y[t][i] - layers[size - 1].val[i]) * (y[t][i] - layers[size - 1].val[i]));
				entropy += (-(1) * log(1 - layers[size - 1].val[i]));
			else
				entropy += (-(1 * log(layers[size - 1].val[i])));
			//std::cout << -( log(1 + layers[size - 1].val[i])) <<' '<<  log(1 - layers[size - 1].val[i]) << '\n';
		}
		//std::cout << entropy << '\n';
		//测试输出 
		double ma = -10.0;
		int ans = -99.99;
		for (int xx = 0; xx < 10; xx++)
		{
			//std::cout << layers[1].val[xx] << ' ';
			if (layers[size-1].val[xx] >= ma)
			{
				ma = layers[size-1].val[xx];
				ans = xx;
			}
		}
		//std::cout << '\n';
		int ans2 = y[t];
		//if(tt%10==7)
			//std::cout << entropy/10 << ' '<<ans<<' '<<ans2<<'\n';
		//测试完毕
		if (fabs(ans - ans2) <= 0.001)
			correct++;
		//std::cout << ans << " " << ans2 << ' ' << entropy << '\n';
		for (int i = 0; i < 10; i++)
		{
			if (fabs(y[t] - i) < 0.000001)
			{
				//delta.M[0][i]= layers[size - 1].val[i] - y[t][i];
				delta.M[0][i] = (layers[size - 1].val[i] - 1) / 10;
			}
			else
			{
				delta.M[0][i] = (layers[size - 1].val[i]) / 10;
			}
		}
		//delta.print();
		layers[size - 1].bias = layers[size - 1].bias - delta * studyrate;
		Matrix deltaw(12 * 12 *players->num, layers[size - 1].size);
		//std::cout <<12* 12* players->num << '\n';
		//Matrix deltaw(layers[size - 2].size, layers[size - 1].size);
		//Matrix tempw(layers[size - 2].size, layers[size - 1].size);
		Matrix tempw(12 * 12 * players->num, layers[size - 1].size);
		for (int i = 0; i < layers[size - 1].presize; i++)
		{
			for (int j = 0; j < layers[size - 1].size; j++)
			{
				deltaw.M[i][j] = layers[size - 1].preval[i] * delta.M[0][j];
				tempw.M[i][j] = layers[size - 1].w.M[i][j];
				//cout << n1.val.M[0][i] << " " << delta.M[0][j] << endl;
			}
		}
		//delta.print();
		//deltaw.print();
		layers[size - 1].w = layers[size - 1].w - deltaw * studyrate;
		//layers[size - 1].w.print();
		//反向传导
		for (int i = size - 2; i >= 0; i--)
		{
			Matrix tranDelta = delta.tran();
			//tranDelta.print();
			delta = delta.multi(tempw, tranDelta);
			//delta.print();
			//std::cout << delta.row << ' ';
			//std::cout << sigmoid2(layers[0].origin.M[0][1]) * (1 - sigmoid2(layers[0].origin.M[0][1]))<<' ';
			tranDelta.column = layers[i].size;
			tranDelta.row = 1;
			//delta.print();
			for (int j = 0; j < layers[i].size; j++)
				tranDelta.M[0][j] = delta.M[j][0] * sigmoid2(layers[i].origin.M[0][j]) * (1 - sigmoid2(layers[i].origin.M[0][j]));
			//tranDelta.print();
			delta = tranDelta;
			//delta.print();
			layers[i].bias = layers[i].bias - delta * studyrate;
			deltaw.reset();
			deltaw.column = layers[i].size;
			deltaw.row = layers[i].presize;
			tempw.reset();
			tempw.row = layers[i].presize;
			tempw.column = layers[i].size;
			for (int ii = 0; ii < layers[i].presize; ii++)
			{
				for (int j = 0; j < layers[i].size; j++)
				{
					tempw.M[ii][j] = layers[i].w.M[ii][j];
					deltaw.M[ii][j] = layers[i].preval[ii] * delta.M[0][j];
					//cout << n1.val.M[0][i] << " " << delta.M[0][j] << endl;
				}
			}
			//deltaw.print();
			layers[i].w = layers[i].w - deltaw * studyrate;
			//if (tt == 9) std::cout << ' ' << layers[i+1].w.M[45][7];
		}
		Matrix tranDelta(5, 5);
		//接着是池化层
		tranDelta.row = delta.column;
		tranDelta.column = 1;
		for (int i = 0; i < delta.column; i++)
			tranDelta.M[i][0] = delta.M[0][i];
		//tranDelta.print();
		//layers[0].w.print();
		delta = delta.multi(tempw, tranDelta);
		//delta.print();
		for (int j = 0; j < 144*players->num; j++)
			tranDelta.M[0][j] = delta.M[j][0];
		tranDelta.column = 144* players->num;
		tranDelta.row = 1;
		delta = tranDelta;
		//delta.print();
		pos=0;
		for (int k = 0; k < players->num; k++)
		{
			for (int i = 0; i < 12; i++)
			{
				for (int j = 0; j < 12; j++)
				{
					players[0].q[k]->M[i][j] = delta.M[0][pos];
					pos++;
				}
			}
		}
		//players[0].q[0]->print();
		//players[0].q[1]->print();
		//最后是卷积层
		for (int nn = 0; nn < clayers->num; nn++)
		{
			for (int i = 0; i < 12; i++)
			{
				for (int j = 0; j < 12; j++)
				{
					for (int k = 0; k < 2; k++)
					{
						for (int l = 0; l < 2; l++)
						{
							clayers[0].q[nn]->M[i * 2 + k][j * 2 + l] = players[0].q[nn]->M[i][j] / 4 * sigmoid2(clayers[0].val[nn]->M[i * 2 + k][j * 2 + l]) * (1 - sigmoid2(clayers[0].val[nn]->M[i * 2 + k][j * 2 + l]));
						}
					}
				}
			}
		}
		//clayers[0].q[0]->print();
		//clayers[0].q[1]->print();
		//卷积层参数更新,首先是对池化层误差超采样，再求和得到卷积层误差
		double del[30];
		for (int i = 0; i < 10; i++) del[i] = 0.0;
		for (int k = 0; k < clayers->num; k++)
		{
			for (int i = 0; i < 24; i++)
			{
				for (int j = 0; j < 24; j++)
				{
					del[k] += clayers[0].q[k]->M[i][j];
				}
			}
			//std::cout << ' ' << del[k];
		}
		//std::cout << '\n';
		//std::cout << ' '<<del[0];
		for(int k=0;k<clayers->num;k++)
			clayers[0].bias[k] = clayers[0].bias[k]-del[k]*studyrate;
		//std::cout << clayers[0].bias[0]<< ' ' << clayers[0].bias[1]<<' '<< clayers[0].bias[2]<<'\n';
		//模仿在一般神经元的BP推导，delta w = delta l * x(l-1) , 将x (l -1 )置换为原来的输入矩阵，并求和
		Matrix delW(5, 5);
		for (int nn = 0; nn < clayers->num; nn++)
		{
			delW.reset();
			//delW.print();
			for (int i = 0; i < 24; i++)
			{
				for (int j = 0; j < 24; j++)
				{
					for (int k = 0; k < 5; k++)
					{
						for (int l = 0; l < 5; l++)
						{
							delW.M[k][l] += clayers[0].q[nn]->M[i][j] * clayers[0].preval.M[i+k][j+l];
						}
					}
				}
			}
			for (int ii = 0; ii < 5; ii++)
			{
				for (int jj = 0; jj < 5; jj++)
				{
					clayers->w[nn]->M[ii][jj] = clayers->w[nn]->M[ii][jj] - studyrate * delW.M[ii][jj];
				}
			}
			//delW.print();
			//std::cout << '\n';
		}
		//delW.print();
		//clayers->w[0] = clayers->w[0]-delW;
		//clayers->w[0]->print();
		//std::cout << clayers->bias[0] << '\n';
	}
	std::cout <<correct << ' ' << testsize;
}
void network::init()
{
	this->studyrate = 0.25;
	this->size = 2;
	this->convSize = 1;
	// 卷积网络
	this->size = 2;
	//layers[1] = NeuronLayer(10, 100, *x[0]);
	//layers[1].w.init(layers[1].size);
	//layers[1] = NeuronLayer(10, 120, layers[0].val);
	//(int num, int distance, int sizex, int sizey, int presizex, int presizey, Matrix * preval)
	clayers[0] = ConvelutionLayer(8, 1, 5, 5, 28, 28, layers[0].bias);
	for (int i = 0; i < 8; i++)
		clayers->w[i]->init(24);
	//PoolingLayer(int num, int sizex, int sizey, int presizex, int presizey, Matrix** preval)
	players[0] =PoolingLayer(8, 2, 2, 24, 24, *clayers[0].val);
	double ss[144*8];
	for (int i = 0; i < 144 * 8; i++)
		ss[i] = 0.0;
	layers[0] = NeuronLayer(100, 144*8, ss);
	layers[0].w.init(100);
	layers[1] = NeuronLayer(10, 100, ss);
	layers[1].w.init(10);
	//layers[2]= NeuronLayer(10, 20, layers[1].val);
	/*
	this->size = 2;
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 28; j++)
		{
			for (int k = 0; k < 28; k++)
			{
				//std::cout << x[i][j][k] << ' ';
			}
			//std::cout << '\n';
		}
		//std::cout << y[i] << '\n';
	}
	layers[0] = NeuronLayer(100, 28*28, *x[0]);
	layers[0].w.init(layers[0].size);
	layers[1]= NeuronLayer(10, 100, *x[0]);
	layers[1].w.init(layers[1].size);
	//layers[0].w.print();
	*/
}
network::~network()
{
}
void network::fullBPTest(int testsize)
{
	int correct = 0;
	for (int t = 0; t < testsize; t++) //对一个小批样本进行训练,一轮100次，样本随机
	{
		//int t = rand() % 5000;
		for (int i = 0; i < 28; i++) //输入矩阵
		{
			for (int j = 0; j < 28; j++)
				//for(int k=0;k<clayers->num;k++)
				clayers[0].preval.M[i][j] = this->x[t][i][j];
		}
		//前向传播
		//卷积层
		clayers[0].output();
		//clayers[0].val[0]->print();
		//clayers[0].w[0]->print();
		//std::cout << clayers[0].bias[0] << ' ';
		for (int i = 0; i < clayers[0].num; i++)
		{
			players[0].preval[i] = clayers[0].val[i];
		}
		players[0].output();
		//players[0].val[0]->print();
		//players[0].val[1]->print();
		//players[0].val[2]->print();
		//池化层->输出层
		double poolOut[144 * 10]; //将输出矩阵映射到样本空间
		int pos = 0;
		for (int k = 0; k < players->num; k++)
		{
			for (int i = 0; i < 12; i++)
			{
				for (int j = 0; j < 12; j++)
				{
					poolOut[pos++] = players[0].val[k]->M[i][j];
					//std::cout << poolOut[pos - 1] << ' ';
				}
			}
		}
		for (int i = 0; i < pos; i++)
			layers[0].preval[i] = poolOut[i];
		//for (int i = 0; i < pos; i++)
		//	std::cout << layers[0].preval[i] << ' ';
		//std::cout << '\n';
		/* 神经层单元测试
		double inpp[10005];
		int pos = 0;
		for (int i = 0; i < 28; i++)
		{
			for (int j = 0; j < 28; j++)
			{
				layers[0].preval[pos++] = x[t][i][j];
			}
		}
		*/
		//隐层和输出层
		for (int i = 0; i < this->size; i++) //前一层输入作为后一层输出
		{
			layers[i].output();
			if (i + 1 < this->size)
				for (int j = 0; j < layers[i].size; j++)
					layers[i + 1].preval[j] = layers[i].val[j];
		}
		Matrix delta(1, 10);
		double entropy = 0.0;
		for (int i = 0; i < 10; i++)
		{
			if (fabs(y[t] - i) < 0.000001)
				//entropy += (0.5 * (y[t][i] - layers[size - 1].val[i]) * (y[t][i] - layers[size - 1].val[i]));
				entropy += (-(1) * log(1 - layers[size - 1].val[i]));
			else
				entropy += (-(1 * log(layers[size - 1].val[i])));
			//std::cout << -( log(1 + layers[size - 1].val[i])) <<' '<<  log(1 - layers[size - 1].val[i]) << '\n';
		}
		double ma = -10.0;
		int ans = -99.99;
		for (int xx = 0; xx < 10; xx++)
		{
			//std::cout << layers[1].val[xx] << ' ';
			if (layers[size-1].val[xx] >= ma)
			{
				ma = layers[size-1].val[xx];
				ans = xx;
			}
		}
		int ans2 = y[t];
		if (fabs(ans - ans2) <= 0.001)
			correct++;
	}
	std::cout << correct << ' ' << testsize;
}