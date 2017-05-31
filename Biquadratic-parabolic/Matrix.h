#pragma once
#include "ListOfAdjacency.h"
#include "Grid.h"
#include "GlobalMatrix.h"
#include "LOS.h"

using namespace std;

extern double F(double x, double y, double t);
extern double U(double x, double y, double t);
extern double Fvector(double Hx, double Hy, int i, double *f);
extern double Mmatrix(double Hx, double Hy, int i, int j);
extern double Gmatrix(double Hx, double Hy, int i, int j);

class Matrix
{
public:
	int n;	
	vector<double> di, al, au;
	vector<int> ia, ja;
	
	const Matrix operator+(const Matrix &B)
	{
		Matrix result;
		result.resize(*this);
		for (int i = 0; i < al.size(); ++i)
		{
			result.al[i] = al[i] + B.al[i];
			result.au[i] = au[i] + B.au[i];
		}
		for (int i = 0; i < di.size(); ++i)
		{
			result.di[i] = di[i] + B.di[i];
		}
		return result;
	};
	Matrix& operator=(const Matrix& B)
	{
			al = B.al;
			au = B.au;
			di = B.di;
			ia = B.ia;
			ja = B.ja;
			n = B.n;
			return *this;
	}
	const Matrix operator*(double num)
	{
		Matrix result;
		result.resize(*this);
		for (int i = 0; i < al.size(); ++i)
		{
			result.al[i] = al[i] * num;
			result.au[i] = al[i] * num;
		}
		for (int i = 0; i < di.size(); ++i)
		{
			result.di[i] = al[i] * num;
		}
		return result;
	}

	const vector <double> operator*(vector <double> B)
	{
		return this->mul(B);
	}

	void resize(Matrix _A)
	{
		di.resize(_A.di.size(), 0);
		al.resize(_A.al.size(), 0);
		au.resize(_A.au.size(), 0);
		ia = _A.ia;
		ja = _A.ja;
		n = _A.n;
	};

	int setEl(int i, int j, double El);
	vector <double> mul(vector <double> B);
	double GetEl(int i, int j);
	void profileDefining(ListOfAdjacency listOfAdjacency, Grid grid);
	void nullMatrix();
	int setEl_BC(int i, int j, double El);
	void outMatrix(int t);
};