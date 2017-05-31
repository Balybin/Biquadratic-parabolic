#include "ListOfAdjacency.h"

void ListOfAdjacency::fillingList(Grid grid)
{
	int n = grid.n, m = grid.m;
	list.resize(m * n);
	for (int j = 0; j < m; ++j)
	{
		for (int i = 0; i < n; ++i)
		{
			int k = grid.calculatePosistion(i, j);
			for (int j1 = -1; j1 < 2; ++j1)
			{
				for (int i1 = -1; i1 < 2; ++i1)
				{
					int ii = i + i1, jj = j + j1;
					int k1 = grid.calculatePosistion(ii, jj);
					if (k1 != -1 && k1 < k)
					{
						list[k].push_back(k1);
					}
				}
			}
		}
	}
}