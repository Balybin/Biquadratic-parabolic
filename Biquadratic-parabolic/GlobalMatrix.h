#pragma once
#include "Vector.h"


class GlobalMatrix
{
	public:
        std::vector<double> di; //диагональные элементы
        std::vector<double> ggl; // элементы нижней диагонали
        std::vector<double> ggu; // элементы верхней диагонали
        std::vector<int> jg; // массив с индексами начала строк(столбцов) в массиве ggl(ggu) 
        std::vector<int> ig; // номера столбцов(строк) элементов в массиве ggl(ggu)
		Vector F; // вектор правой части
		std::vector<dot> list; //список смежных элементов, храню, потому что могу
        int n; // размерность матрицы
        double normF; // норма вектора F
		GlobalMatrix(void);
		~GlobalMatrix(void);
		void addElem(std::vector<int> el, std::vector<std::vector<double>> &G, std::vector<std::vector<double>> &M, std::vector<double> &F);
        void make(std::vector<pt> &point, std::vector<std::vector<int>> &finElem);
        void multMV(Vector &X, Vector &Y);
		void addToList(std::vector<int> &el);
		int findInList(dot x);
		void sortList();
		void addFirstBound(std::vector<int> &el,std::vector<pt> &point);
		void addSecondBound(std::vector<int> &el,std::vector<pt> &point);
		void addThirdBound(std::vector<int> &el,std::vector<pt> &point);
};

