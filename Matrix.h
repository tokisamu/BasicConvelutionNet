#pragma once
class Matrix
{
public:
	int row, column;
	double M[1153][1153];
	Matrix(int row,int column);
	Matrix();
	~Matrix();
	Matrix multi(Matrix a, Matrix b);
	Matrix multi(double a[], Matrix b);
	void print();
	void reset();
	void init(int x);
	Matrix tran();
	Matrix operator +(Matrix b)
	{
		Matrix c(this->row, this->column);
		int ro, col;
		ro = this->row;
		col = this->column;
		if (ro > b.row) ro = b.row;
		if (col > b.column) col = b.column;
		for (int i = 0; i <ro; i++)
		{
			for (int j = 0; j < col; j++)
			{
				c.M[i][j] = this->M[i][j] + b.M[i][j];
			}
		}
		return c;
	}
	Matrix operator *(double b)
	{
		Matrix c(this->row, this->column);
		int ro, col;
		ro = this->row;
		col = this->column;
		for (int i = 0; i < ro; i++)
		{
			for (int j = 0; j < col; j++)
			{
				c.M[i][j] = this->M[i][j] * b;
			}
		}
		return c;
	}
	Matrix operator -(Matrix b)
	{
		int ro, col;
		ro = this->row;
		col = this->column;
		if (ro < b.row) ro = b.row;
		if (col < b.column) col = b.column;
		Matrix c(this->row, this->column);
		for (int i = 0; i < ro; i++)
		{
			for (int j = 0; j < col; j++)
			{
				c.M[i][j] = this->M[i][j] - b.M[i][j];
			}
		}
		return c;
	}
};

