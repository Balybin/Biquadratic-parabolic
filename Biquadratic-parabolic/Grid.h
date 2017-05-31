#pragma once
#include <vector>
using namespace std;
class Grid
{
public:
	vector<double> X, Y, T;
	int n, m;
	int calculatePosistion(int i, int j);
	void inputGrid();
	void doubleGrid();
};