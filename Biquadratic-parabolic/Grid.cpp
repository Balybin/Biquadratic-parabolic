#include "Grid.h"
#include <fstream>
#include <iostream>
#include <math.h>;

double log(double a, double b)
{
	return log(b) / log(a);
}

int Grid::calculatePosistion(int i, int j)
{
	if (i >= n || j >= m || i < 0 || j < 0 )
		return -1;
	int k = j * n + i;
	if (k >= 0)
		return k;
	else
		return -1;
}


void Grid::inputGrid()
{
	ifstream file("grid.txt");
	vector<double> toBorder, toMinStep, toDecompactionCoef;
	int count, ibuf, k;
	double buf;
	file >> buf;//for X
	toBorder.push_back(buf);
	file >> count;
	for (int i = 1; i <= count; i++)
	{
		file >> buf;
		toBorder.push_back(buf);
	}
	for (int i = 0; i < count; i++)
	{
		file >> buf;
		toMinStep.push_back(buf);
	}
	for (int i = 0; i < count; i++)
	{
		file >> buf;
		toDecompactionCoef.push_back(buf);
	}
	X.push_back(toBorder[0]);
	for (int i = 0; i < count; ++i)
	{
		double Xk = toBorder[i + 1] - toBorder[i];
		if (toDecompactionCoef[i] == 1)
			k = toDecompactionCoef[i] / toMinStep[i];
		else
			k = log(toDecompactionCoef[i], Xk*(toDecompactionCoef[i] - 1) / toMinStep[i]);
		int j = 1;
		if (toDecompactionCoef[i] == 1)
			buf = toMinStep[i] * j;
		else
			buf = toBorder[i] + toMinStep[i]*(pow(toDecompactionCoef[i],j) - 1)/(toDecompactionCoef[i] - 1);
		while (buf < toBorder[i + 1])
		{
			X.push_back(buf);
			++j;
			if (toDecompactionCoef[i] == 1)
				buf = toMinStep[i] * j;
			else
				buf = toBorder[i] + toMinStep[i] * (pow(toDecompactionCoef[i], j) - 1) / (toDecompactionCoef[i] - 1);
		}
		if (toBorder[i + 1] - X[X.size() - 1] < (X[X.size() - 1] - X[X.size() - 2]) / 2)
			X.pop_back();
		X.push_back(toBorder[i + 1]);
	}
	//now for Y
	toBorder.clear();
	toMinStep.clear();
	toDecompactionCoef.clear();
	file >> buf;
	toBorder.push_back(buf);
	file >> count;
	for (int i = 1; i <= count; i++)
	{
		file >> buf;
		toBorder.push_back(buf);
	}
	for (int i = 0; i < count; i++)
	{
		file >> buf;
		toMinStep.push_back(buf);
	}
	for (int i = 0; i < count; i++)
	{
		file >> buf;
		toDecompactionCoef.push_back(buf);
	}
	Y.push_back(toBorder[0]);
	for (int i = 0; i < count; ++i)
	{
		double Xk = toBorder[i + 1] - toBorder[i];
		if (toDecompactionCoef[i] == 1)
			k = toDecompactionCoef[i] / toMinStep[i];
		else
			k = log(toDecompactionCoef[i], Xk*(toDecompactionCoef[i] - 1) / toMinStep[i]);
		int j = 1;
		if (toDecompactionCoef[i] == 1)
			buf = toMinStep[i] * j;
		else
			buf = toBorder[i] + toMinStep[i] * (pow(toDecompactionCoef[i], j) - 1) / (toDecompactionCoef[i] - 1);
		while (buf < toBorder[i + 1])
		{
			Y.push_back(buf);
			++j;
			if (toDecompactionCoef[i] == 1)
				buf = toMinStep[i] * j;
			else
				buf = toBorder[i] + toMinStep[i] * (pow(toDecompactionCoef[i], j) - 1) / (toDecompactionCoef[i] - 1);
		}
		if (toBorder[i + 1] - Y[Y.size() - 1] < (Y[Y.size() - 1] - Y[Y.size() - 2]) / 2)
			Y.pop_back();
		Y.push_back(toBorder[i + 1]);
	}
	toBorder.clear();
	toMinStep.clear();
	toDecompactionCoef.clear();
	//now for T
	file >> buf;
	toBorder.push_back(buf);
	file >> count;
	for (int i = 1; i <= count; i++)
	{
		file >> buf;
		toBorder.push_back(buf);
	}
	for (int i = 0; i < count; i++)
	{
		file >> buf;
		toMinStep.push_back(buf);
	}
	for (int i = 0; i < count; i++)
	{
		file >> buf;
		toDecompactionCoef.push_back(buf);
	}
	T.push_back(toBorder[0]);
	for (int i = 0; i < count; ++i)
	{
		double Xk = toBorder[i + 1] - toBorder[i];
		if (toDecompactionCoef[i] == 1)
			k = toDecompactionCoef[i] / toMinStep[i];
		else
			k = log(toDecompactionCoef[i], Xk*(toDecompactionCoef[i] - 1) / toMinStep[i]);
		int j = 1;
		if (toDecompactionCoef[i] == 1)
			buf = toMinStep[i] * j;
		else
			buf = toBorder[i] + toMinStep[i] * (pow(toDecompactionCoef[i], j) - 1) / (toDecompactionCoef[i] - 1);
		while (buf < toBorder[i + 1])
		{
			T.push_back(buf);
			++j;
			if (toDecompactionCoef[i] == 1)
				buf = toMinStep[i] * j;
			else
				buf = toBorder[i] + toMinStep[i] * (pow(toDecompactionCoef[i], j) - 1) / (toDecompactionCoef[i] - 1);
		}
		if (toBorder[i + 1] - T[T.size() - 1] < (T[T.size() - 1] - T[T.size() - 2]) / 2)
			T.pop_back();
		T.push_back(toBorder[i + 1]);
	}
	n = X.size();
	m = Y.size();
}

void Grid::doubleGrid()
{
	int buf = X.size() - 1;
	int i = 0;
	while (i < buf)
	{
		X.insert(X.begin() + 1 + i, X[i] + (X[i + 1] - X[i]) / 2);
		buf++;
		i += 2;
	}
	n = buf + 1;
	buf = Y.size() - 1;
	i = 0;
	while (i < buf)
	{
		Y.insert(Y.begin() + 1 + i, Y[i] + (Y[i + 1] - Y[i]) / 2);
		buf++;
		i += 2;
	}
	m = buf + 1;
}
