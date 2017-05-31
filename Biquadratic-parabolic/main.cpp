#include <iostream>
#include <fstream>
#include "Matrix.h"
#include "LOS.h"
#include "Slau.h"
#include <cmath>
#include <string>
//i	| j  
//n	| m  <--вот их вообще логично было бы поменять местами, но
//x	| y  я глупенький и мне лень переписывать, поэтому будет так
using namespace std;

double lyambda=1, gamma=1;

double U(double x, double y, double t)
{
	return x + y;					//F=U
//	return x * x + y * y;			//F=x^2+y^2-4
//	return x * x * x + y * y * y;	//F=x^3-6x+y^3-6y
//	return x*x*x*x + y*y*y*y;		//F=x^4-12x^2+y^4-12y^2
//	return exp(x) + exp(y);			//F = 0
//	return sin(x + y);				//F = 3 * sin(x + y);
}

double F(double x, double y, double t)
{
	return (x + y);										//U=F
//	return x * x + y * y - 4;							//U=x^2+y^2
//	return x*x*x - 6*x + y*y*y - 6*y;			//U=x^3+y^3
//	return x*x*x*x + y*y*y*y - 12 * x*x - 12 * y*y;	//U=x^4+y^4
//	return 0;											//U = e^x+e^y
//	return 3 * sin(x + y);								//U = sin(x + y);	
}
double Gmatrix(double Hx, double Hy, int i, int j)
{
	int G1[4][4] = {
	2, -2, 1, -1,
	-2, 2, -1, 1,
	1, -1, 2, -2,
	-1, 1, -2, 2 };
	int G2[4][4] = {
	2, 1, -2, -1,
	1, 2, -1, -2,
	-2, -1, 2, 1,
	-1, -2, 1, 2 };
	return lyambda / 6 * Hy / Hx * G1[i][j] + lyambda / 6 * Hx / Hy * G2[i][j];
}
double Mmatrix(double Hx, double Hy, int i, int j)
{
	int M[4][4] = {
		4, 2, 2, 1,
		2, 4, 1, 2,
		2, 1, 4, 2,
		1, 2, 2, 4 };
	return gamma * Hx * Hy / 36. * (double)M[i][j];
}
double Fvector(double Hx, double Hy, int i, double *f)
{
	int M[4][4] = {
		4, 2, 2, 1,
		2, 4, 1, 2,
		2, 1, 4, 2,
		1, 2, 2, 4 };

	double sum = 0;
	for (int j = 0; j < 4; j++)
	{
		sum += M[i][j] * f[j];
	}
	return Hx*Hy / 36 * sum;
}

void outvector(int i, vector<double> a)
{
	ofstream out("t=" + to_string(i) + "Result.txt");
	for (int i = 0; i < a.size(); ++i)
		out << a[i] << endl;
	out.close();
}

void main()
{
	LOS los;
	Slau slau;
	slau.input();
	outvector(0, slau.qprev2);
	outvector(1, slau.qprev1);

	for (slau.t = 2; slau.t < slau.grid.T.size(); ++slau.t)
	{
		slau.make();
		outvector(slau.t, slau.los.x0.V);
		slau.nullMatrix();
	}
}