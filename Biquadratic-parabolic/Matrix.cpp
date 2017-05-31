#include <string>
#include "Matrix.h"


void Matrix::outMatrix(int t)
{
	ofstream out("Mmatrix"+std::to_string(t)+".txt");
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			out << GetEl(i, j) << "\t";
		}
		out << endl;
	}
}

//void Matrix::input()
//{
	
//	grid.inputGrid();
//	profileDefining();
//	qprev2.reserve(grid.X.size()*grid.Y.size());
//	qprev1.reserve(grid.X.size()*grid.Y.size());
//	for (int j = 0; j < grid.Y.size(); ++j)
//	{
	//	for (int i = 0; i < grid.X.size(); ++i)
	//	{
	//		int k = grid.calculatePosistion(i, j);
//			qprev2.push_back(U(grid.X[i], grid.Y[j], grid.T[0]));
//			qprev1.push_back(U(grid.X[i], grid.Y[j], grid.T[1]));
//		}
//	}
//}
void Matrix::nullMatrix()
{
	for (int i = 0; i < al.size(); ++i)
	{
		al[i] = 0;
		au[i] = 0;
	}
	for (int i = 0; i < di.size(); ++i)
	{
		di[i] = 0;
	}
}


int Matrix::setEl_BC(int i, int j, double El)
{
	if (i == j)
	{
		di[i] = El;
		return 0;
	}
	if (i > j)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k<i1; k++)
			if (ja[k] == j)
			{
				al[k] = El;
				return 0;
			}
		return -1;
	}
	else
	{
		swap(i, j);
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k<i1; k++)
			if (ja[k] == j)
			{
				au[k] = El;
				return 0;
			}
		return -1;
	}
}

vector <double> Matrix::mul(vector <double> B)
{
	vector <double> result;
	int n = B.size();
	result.resize(n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			result[i] += GetEl(i, j)*B[j];
		}
	}
	return result;
}

double Matrix::GetEl(int i, int j)
{
	if (i == j)
	{
		return di[i];
	}
	if (i > j)
	{
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k<i1; k++)
			if (ja[k] == j)
			{
				return al[k];
			}
		return 0;
	}
	else
	{
		swap(i, j);
		int i0 = ia[i];
		int i1 = ia[i + 1];
		for (int k = i0; k<i1; k++)
			if (ja[k] == j)
			{
				return au[k];
			}
		return 0;
	}
}
int Matrix::setEl(int i, int j, double El)
{
	if (i == j)
	{
		di[i] += El;
		return 0;
	}
	if(i > j)
	{
		int i0 = ia[i];
		int i1 = ia[i+1];
		for(int k=i0; k<i1; k++)
			if (ja[k] == j)
			{
				al[k] += El;
				return 0;
			}
		return -1;
	}
}


void Matrix::profileDefining(ListOfAdjacency listOfAdjacency, Grid grid)
{
	n = grid.n * grid.m;
	int k = 0;
	listOfAdjacency.fillingList(grid);
	ia.reserve(n + 1);
	ja.reserve(n);
	ia.push_back(0);
	di.reserve(n);
	al.reserve(4 * n);
	for (int i = 0; i < n; i++)
	{
		di.push_back(0);
		for (int j = 0; j < listOfAdjacency.list[i].size(); j++)
		{
			ja.push_back(listOfAdjacency.list[i][j]);
			al.push_back(0);
			au.push_back(0);
			k++;
		}
		ia.push_back(k);
	}
}