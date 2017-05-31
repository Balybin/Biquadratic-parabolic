#pragma once
#include "grid.h"
#include "Matrix.h"
#include "LOS.h"

extern double F(double x, double y);
extern double U(double x, double y, double t);

class Slau
{
public:
	Matrix A, M, G;
	ListOfAdjacency listOfAdjacency;
	Grid grid;
	LOS los;
	double deltaT, deltaT0, deltaT1;
	vector<double> f, qprev1, qprev2, mq2, mq1;
	double LocM[9][9], LocG[9][9], LocF[9], f1[9], Hx, Hy;
	int index[9], t;

	void input();
	void matrixFilling();
	void make();
	void firstBoundaryConditions();
	void nullMatrix();
};