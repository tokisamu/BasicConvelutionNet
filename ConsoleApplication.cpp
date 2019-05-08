// ConsoleApplication.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "Matrix.h"
#include"NeuronLayer.h"
#include"network.h"
#include"ConvelutionLayer.h"
#include"PoolingLayer.h"
#include <iostream>
#include<vector>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;
int ReverseInt(int i)
{
	unsigned char ch1, ch2, ch3, ch4;
	ch1 = i & 255;
	ch2 = (i >> 8) & 255;
	ch3 = (i >> 16) & 255;
	ch4 = (i >> 24) & 255;
	return((int)ch1 << 24) + ((int)ch2 << 16) + ((int)ch3 << 8) + ch4;
}

void read_Mnist_Label(string filename, vector<double> & labels)
{
	ifstream file(filename, ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		file.read((char*)& magic_number, sizeof(magic_number));
		file.read((char*)& number_of_images, sizeof(number_of_images));
		magic_number = ReverseInt(magic_number);
		number_of_images = ReverseInt(number_of_images);
		cout << "magic number = " << magic_number << endl;
		cout << "number of images = " << number_of_images << endl;


		for (int i = 0; i < number_of_images; i++)
		{
			unsigned char label = 0;
			file.read((char*)& label, sizeof(label));
			labels.push_back((double)label);
		}

	}
}

void read_Mnist_Images(string filename, vector<vector<double>> & images)
{
	ifstream file(filename, ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int number_of_images = 0;
		int n_rows = 0;
		int n_cols = 0;
		unsigned char label;
		file.read((char*)& magic_number, sizeof(magic_number));
		file.read((char*)& number_of_images, sizeof(number_of_images));
		file.read((char*)& n_rows, sizeof(n_rows));
		file.read((char*)& n_cols, sizeof(n_cols));
		magic_number = ReverseInt(magic_number);
		number_of_images = ReverseInt(number_of_images);
		n_rows = ReverseInt(n_rows);
		n_cols = ReverseInt(n_cols);

		cout << "magic number = " << magic_number << endl;
		cout << "number of images = " << number_of_images << endl;
		cout << "rows = " << n_rows << endl;
		cout << "cols = " << n_cols << endl;

		for (int i = 0; i < number_of_images; i++)
		{
			vector<double>tp;
			for (int r = 0; r < n_rows; r++)
			{
				for (int c = 0; c < n_cols; c++)
				{
					unsigned char image = 0;
					file.read((char*)& image, sizeof(image));
					tp.push_back(image);
				}
			}
			images.push_back(tp);
		}
	}
}

int main()
{
	//建立神经元网络
	network n1;
	//读取mnist
	vector<vector<double>>images;
	//read_Mnist_Images("t10k-images.idx3-ubyte", images);
	read_Mnist_Images("train-images.idx3-ubyte", images);
	for (int i = 0; i < 40000; i++)
	{
		for (int j = 0; j < images[0].size(); j++)
		{
			n1.x[i][j / 28][j % 28] = images[i][j];
			//cout << images[i][j] << " ";
			//if (j % 28 == 0) cout << endl;
		}
	}
	vector<double>labels;
	//read_Mnist_Label("t10k-labels.idx1-ubyte", labels);
	read_Mnist_Label("train-labels.idx1-ubyte", labels);
	for (int i = 0; i < 40000; i++)
	{
		n1.y[i] = labels[i];
		//cout << labels[i] << endl;
	}
	n1.init();
	for (int i = 0; i <0; i++) //朴素神经网络
	{
		cout << "Round" << ' ' << i << ' ';
		n1.BP(100);
		cout << endl;
	}
	for (int i = 0; i < 1000; i++) //卷积神经网络
	{
		cout << "Round" << ' ' << i << ": ";
		n1.fullBP(100);
		cout << endl;
	}
	images.clear();
	read_Mnist_Images("t10k-images.idx3-ubyte", images);
	for (int i = 0; i < 10000; i++)
	{
		for (int j = 0; j < images[0].size(); j++)
		{
			n1.x[i][j / 28][j % 28] = images[i][j];
			//cout << images[i][j] << " ";
			//if (j % 28 == 0) cout << endl;
		}
	}
	labels.clear();
	read_Mnist_Label("t10k-labels.idx1-ubyte", labels);
	for (int i = 0; i < 10000; i++)
	{
		n1.y[i] = labels[i];
		//cout << labels[i] << endl;
	}
	//n1.BP(2000);
	n1.fullBPTest(10000);
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
