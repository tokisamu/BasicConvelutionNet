#include "Matrix.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
using namespace std;

//gaussrand v is variance
double gaussrand(int V)
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
	double aaa = 1.0001 / (V + 1);
	double sq = sqrt(aaa);
	if (phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
	double ans = X * sq;
	return ans;
}

Matrix::Matrix(int row,int column)
{
	this->row = row;
	this->column = column;
	for (int i = 0; i < row; i++)
		for (int j = 0; j < column; j++)
			this->M[i][j] = 0.0;
}
Matrix::Matrix()
{
	this->row = 10;
	this->column = 10;
	for (int i = 0; i < 10; i++)
		for (int j = 0; j < 10; j++)
			this->M[i][j] = 0;
}

Matrix::~Matrix()
{
}

Matrix Matrix::multi(Matrix a, Matrix b)
{
	Matrix c(a.row, b.column);
	double cnt = 0.0;
	for (int i = 0; i < a.row; i++)
	{
		for (int j = 0; j < b.column; j++)
		{
			for (int k = 0; k < a.column; k++)
			{
				cnt += a.M[i][k] * b.M[k][j];
			}
			c.M[i][j] = cnt;
			cnt = 0.0;
		}
	}
	return c;
}

Matrix Matrix::multi(double a[], Matrix b)
{
	Matrix c(1, b.column);
	double cnt = 0.0;
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < b.column; j++)
		{
			for (int k = 0; k < b.row; k++)
			{
				cnt += a[k] * b.M[k][j];
			}
			c.M[i][j] = cnt;
			cnt = 0.0;
		}
	}
	return c;
}

void Matrix::print()
{
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->column; j++)
			cout << this->M[i][j] << ' ' ;
		cout << endl;
	}
	cout << endl;
}

void Matrix::reset()
{
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->column; j++)
		{
			this->M[i][j] = 0.0;
		}
	}
}

void Matrix::init(int x)
{
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->column; j++)
		{
			this->M[i][j] = gaussrand(this->row);
		}
	}
}

Matrix Matrix::tran()
{
	Matrix ans(this->column, this->row);
	for (int i = 0; i < this->column; i++)
	{
		for (int j = 0; j < this->row; j++)
		{
			ans.M[i][j] = this->M[j][i];
		}
	}
	return ans;
}
