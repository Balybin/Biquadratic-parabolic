#pragma once
#include "Slau.h"

void Slau::input()
{
	grid.inputGrid();
	grid.doubleGrid();
	A.profileDefining(listOfAdjacency, grid);
	M.resize(A); G.resize(A);
	qprev2.reserve(grid.X.size()*grid.Y.size());
	qprev1.reserve(grid.X.size()*grid.Y.size());
	f.resize(A.n); mq2.resize(A.n); mq1.resize(A.n);
	for (int j = 0; j < grid.Y.size(); ++j)
	{
		for (int i = 0; i < grid.X.size(); ++i)
		{
			int k = grid.calculatePosistion(i, j);
			qprev2.push_back(U(grid.X[i], grid.Y[j], grid.T[0]));
			qprev1.push_back(U(grid.X[i], grid.Y[j], grid.T[1]));
		}
	}
};


void Slau::matrixFilling()
{
	int n = grid.n - 2, m = grid.m - 2;
	for (int j = 0; j < m; j += 2) //по конечным элементам
	{
		for (int i = 0; i < n; i += 2)
		{
			index[0] = grid.calculatePosistion(i, j);
			index[1] = grid.calculatePosistion(i + 1, j);
			index[2] = grid.calculatePosistion(i + 2, j);
			index[3] = grid.calculatePosistion(i, j + 1);
			index[4] = grid.calculatePosistion(i + 1, j + 1);
			index[5] = grid.calculatePosistion(i + 2, j + 1);
			index[6] = grid.calculatePosistion(i, j + 2);
			index[7] = grid.calculatePosistion(i + 1, j + 2);
			index[8] = grid.calculatePosistion(i + 2, j + 2);
			f1[0] = F(grid.X[i], grid.Y[j]);
			f1[1] = F(grid.X[i + 1], grid.Y[j]);
			f1[2] = F(grid.X[i + 2], grid.Y[j]);
			f1[3] = F(grid.X[i], grid.Y[j + 1]);
			f1[4] = F(grid.X[i + 1], grid.Y[j + 1]);
			f1[5] = F(grid.X[i + 2], grid.Y[j + 1]);
			f1[6] = F(grid.X[i], grid.Y[j + 2]);
			f1[7] = F(grid.X[i + 1], grid.Y[j + 2]);
			f1[8] = F(grid.X[i + 2], grid.Y[j + 2]);
			Hx = grid.X[i + 1] - grid.X[i];
			Hy = grid.Y[j + 1] - grid.Y[j];
			for (int i1 = 0; i1 < 9; i1++)
			{
				for (int j1 = 0; j1 < 9; j1++)
				{
					LocG[i1][j1] = Gmatrix(Hx, Hy, i1, j1);// +(deltaT + deltaT0) / (deltaT * deltaT0)*
					LocM[i1][j1] = Mmatrix(Hx, Hy, i1, j1);
				}
				LocF[i1] = Fvector(Hx, Hy, i1, f1);
			}

			for (int i1 = 0; i1 < 9; i1++)
			{
				for (int j1 = 0; j1 < 9; j1++)
				{
					G.setEl(index[i1], index[j1], LocG[i1][j1]);
					M.setEl(index[i1], index[j1], LocM[i1][j1]);
				}
				f[index[i1]] += LocF[i1];
			}

		}
	}
	M.au = M.al;
	G.au = G.al;
	A = G + M*((deltaT + deltaT0) / (deltaT * deltaT0));

	mq2 = M*(qprev2);
	for (int i = 0; i < f.size(); ++i)
	{
		f[i] -= deltaT0*mq2[i] / (deltaT*deltaT1);
	}
	mq1 = M*(qprev1);
	for (int i = 0; i < f.size(); ++i)
	{
		f[i] += deltaT*mq1[i] / (deltaT1*deltaT0);
	}
};

void Slau::make()
{
	deltaT = grid.T[t] - grid.T[t - 2];
	deltaT1 = grid.T[t - 1] - grid.T[t - 2];
	deltaT0 = grid.T[t] - grid.T[t - 1];
	matrixFilling();
	firstBoundaryConditions();
	los.A.di = A.di;
	los.A.ggl = A.al;
	los.A.ggu = A.au;
	los.A.ig = A.ia;
	los.A.jg = A.ja;
	los.A.n = A.n;
	los.A.F.V = f;
	los.A.normF = sqrt(los.A.F*los.A.F);
	los.make();
	los.LUdec();
	los.LOS_LU();
	qprev2 = qprev1;
	qprev1 = los.x0.V;
};

void Slau::firstBoundaryConditions()
{
	int gridn = grid.n, gridm = grid.m;
	for (int i = 0; i < gridn; ++i)
	{//для нижней границы
		int k = grid.calculatePosistion(i, 0);
		for (int j = 0; j < A.n; ++j)
		{
			A.setEl_BC(k, j, 0);
			A.setEl_BC(j, k, 0);
		}
		A.setEl_BC(k, k, 1);
		f[k] = U(grid.X[k], grid.Y[0], grid.T[t]);
	}
	for (int j = 0; j < gridm; ++j)
	{//для правой границы
		int k = grid.calculatePosistion(gridn - 1, j);
		for (int i = 0; i < A.n; ++i)
		{
			A.setEl_BC(k, i, 0);
			A.setEl_BC(i, k, 0);
		}
		A.setEl_BC(k, k, 1);
		f[k] = U(grid.X[gridn - 1], grid.Y[j], grid.T[t]);
	}
	for (int i = 0; i < gridn; ++i)
	{//для верхней границы
		int k = grid.calculatePosistion(i, gridm - 1);
		for (int j = 0; j < A.n; ++j)
		{
			A.setEl_BC(k, j, 0);
			A.setEl_BC(j, k, 0);
		}
		A.setEl_BC(k, k, 1);
		f[k] = U(grid.X[i], grid.Y[gridm - 1], grid.T[t]);
	}
	for (int j = 0; j < gridm; ++j)
	{//для левой границы, надеюсь, я ничего не перепутал
		int k = grid.calculatePosistion(0, j);
		for (int i = 0; i < A.n; ++i)
		{
			A.setEl_BC(k, i, 0);
			A.setEl_BC(i, k, 0);
		}
		A.setEl_BC(k, k, 1);
		f[k] = U(grid.X[0], grid.Y[j], grid.T[t]);
	}
}

void Slau::nullMatrix()
{
	for (int i = 0; i < f.size(); ++i)
		f[i] = 0;
	A.nullMatrix();
	G.nullMatrix();
	M.nullMatrix();
};