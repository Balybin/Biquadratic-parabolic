#pragma once
#include "Grid.h"
#include <vector>

using namespace std;

class ListOfAdjacency
{
public:
	vector<vector<int> > list;
	void fillingList(Grid grid);
};