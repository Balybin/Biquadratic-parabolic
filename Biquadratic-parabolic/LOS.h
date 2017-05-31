#pragma once
#include "GlobalMatrix.h"
class LOS
{
public:
		GlobalMatrix A;
        Vector p;
        Vector z;
        Vector r;
        Vector x0;
        Vector Ar;
        Vector y;
		std::vector<double> L;
		std::vector<double> U;
		std::vector<double> D;
        int maxiter;
        double eps;
        int iter;
        double normR;
		
		LOS(void);
        ~LOS(void);
		void add(std::vector<int> el, std::vector<std::vector<double>> &G, std::vector<std::vector<double>> &M, std::vector<double> &F);
		void firstBound(std::vector<int> &el,std::vector<pt> &point);
		void secondBound(std::vector<int> &el,std::vector<pt> &point);
		void thirdBound(std::vector<int> &el,std::vector<pt> &point);
		void make();
        void LUdec();
        void Direct(Vector &y, Vector &F);
        void Reverse(Vector &x, Vector &y);
        void LOS_LU();
		void outputMatrix(std::ofstream &X);
        void output(std::ofstream &X);
};

