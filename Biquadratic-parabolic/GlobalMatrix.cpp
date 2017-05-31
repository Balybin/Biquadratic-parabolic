#include "GlobalMatrix.h"

// построение портрета матрицы
void GlobalMatrix::make(std::vector<pt> &point, std::vector<std::vector<int>> &finElem)
{
	n = point.size();
	list.reserve(n);
	di.resize(n);
	F.make(n);
	ig.resize(n + 1);
	for (int i = 0; i < finElem.size(); i++)
	{
		addToList(finElem[i]);
	}
	for (int i = 0; i < list.size(); i++)
	{
		ig[list[i].x + 1]++;
	}
	for (int i = 1; i <= n; i++)
	{
		ig[i] += ig[i - 1];
	}
	jg.resize(ig[n]);
	sortList();
	for (int i = 0; i < ig[n]; i++)
	{
		jg[i] = list[i].y;
	}
	ggl.resize(ig[n]);
	ggu.resize(ig[n]);
}

void GlobalMatrix::addElem(std::vector<int> el, std::vector<std::vector<double>> &G, std::vector<std::vector<double>> &M, std::vector<double> &F)
{
	list.resize(0);
	list.reserve(n);
	addToList(el);
	sortList();
	for (int i = 0, k = 0; i < 4; i++)
	{
		di[el[i]] += G[i][i] + M[i][i];
		this->F.V[el[i]] += F[i];
		for (int j = 0; j < i; k++,j++)
		{
			int x = ig[list[k].x];
			for (; x < ig[list[k].x + 1]; x++)
			{
				if (list[k].y == jg[x])
				{
					ggl[x] += M[i][j] + G[i][j];
					ggu[x] += M[j][i] + G[j][i];
				}
			}
		}
	}
}

// обработка одного конечного элемента, составляем координат матрицы
void GlobalMatrix::addToList(std::vector<int> &el)
{
	dot buf;
	int min;
	int pos;
    for(int i = 0; i < 4; ++i) 
    { 
        pos = i; 
        min = el[i];
        for(int j = i + 1; j < 4; j++)
        {
           if (el[j] < min) 
           {
               pos = j; 
               min = el[j]; 
           }
        }
        std::swap(el[i],el[pos]);
    }
	for (int i = 3; i >= 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			buf.x = el[i];
			buf.y = el[j];
			if (findInList(buf) == 0)
			{
				list.push_back(buf);
			}
		}
	}
}

// поиск одинаковых элементов в списке
int GlobalMatrix::findInList(dot x)
{
	for (int i = 0; i < list.size(); i++)
	{
		if (list[i].x == x.x && list[i].y == x.y)
		{
			return 1;
		}
	}
	return 0;
}

// сортировка списка по первому и второму элементу
void GlobalMatrix::sortList()
{
	int min, pos, minx;
	for(int i = 0; i < list.size(); i++) 
    { 
        pos = i; 
        min = list[i].x;
        for(int j = i + 1; j < list.size(); j++)
        {
           if (list[j].x < min) 
           {
               pos = j; 
               min = list[j].x; 
           }
        }
        std::swap(list[i],list[pos]);
    }
	for (int i = 0; i < list.size(); i++)
	{
		pos = i; 
		min = list[i].y;
		minx = list[i].x;
		for(int j = i + 1; j < list.size(); j++)
		{
			if (list[j].y < min && minx == list[j].x) 
			{
				pos = j; 
				min = list[j].y; 
			}
		}
		std::swap(list[i],list[pos]);
	}
}

void GlobalMatrix::addFirstBound(std::vector<int> &el,std::vector<pt> &point)
{
	di[el[0]] = 1;
	di[el[1]] = 1;
	F.V[el[0]] = first(point[el[0]], el[3]);
	F.V[el[1]] = first(point[el[1]], el[3]);
	// зануляем строки
	for (int i = ig[el[0]]; i < ig[el[0] + 1]; i++)
	{
		ggl[i] = 0;
	}
	for (int i = ig[el[1]]; i < ig[el[1] + 1]; i++)
	{
		ggl[i] = 0;
	}
	for (int i = 0; i < jg.size(); i++)
	{
		if (jg[i] == el[0])
		{
			ggu[i] = 0;
		}
		if (jg[i] == el[1])
		{
			ggu[i] = 0;
		}
	}
}

void GlobalMatrix::addSecondBound(std::vector<int> &el,std::vector<pt> &point)
{

}
void GlobalMatrix::addThirdBound(std::vector<int> &el,std::vector<pt> &point)
{

}

GlobalMatrix::GlobalMatrix(void)
{
}

GlobalMatrix::~GlobalMatrix(void)
{
}

// умножение матрицы на вектор. 1: вектор 2: результат
void GlobalMatrix::multMV(Vector &X, Vector &Y)
{
        int k1,k2;
        for (int i = 0; i < n; i++)
        {
                Y.V[i] = di[i] * X.V[i];
                k1 = ig[i];
                k2 = ig[i+1];
                for(int k = k1; k < k2; k++)
                {
                        int j = jg[k];
                        Y.V[i] += ggl[k] * X.V[j];
                        Y.V[j] += ggu[k] * X.V[i];
                }
        }
}