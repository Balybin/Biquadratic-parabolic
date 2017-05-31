#include "LOS.h"
#include <stdlib.h>
LOS::LOS(void)
{
}

LOS::~LOS(void)
{
}

void LOS::make()
{
        x0.make(A.n);
        r.make(A.n);
        z.make(A.n);
        p.make(A.n);
        Ar.make(A.n);
        y.make(A.n);
		y.V = A.F.V;
		D.resize(A.n);
		L.resize(A.ggl.size());
		U.resize(A.ggu.size());
		maxiter = 10000;
		eps = pow(10.0,-13);
}

void LOS::add(std::vector<int> el, std::vector<std::vector<double>> &G, std::vector<std::vector<double>> &M, std::vector<double> &F)
{
	A.addElem(el, G, M, F);
}

void LOS::firstBound(std::vector<int> &el,std::vector<pt> &point)
{
	A.addFirstBound(el, point);
}
void LOS::secondBound(std::vector<int> &el,std::vector<pt> &point)
{
	A.addSecondBound(el, point);
}
void LOS::thirdBound(std::vector<int> &el,std::vector<pt> &point)
{
	A.addThirdBound(el, point);
}

void LOS::LOS_LU()
{
        double scalPP = 0,scalPR = 0, scalRR = 0, scalPAr = 0;
        double a = 0, b = 0;
        bool exit = 1;
        A.multMV(x0,y); // y = A * x
        A.F - y;                // y = F - y
        Direct(r,y);// r = L-1 * y
        Reverse(z,r); // z = U-1 * r
        A.multMV(z,y);  // y = A * z
        Direct(p,y);// p = L-1 * y 
        scalRR = r * r;
        normR = sqrt(scalRR)/A.normF;
        for (iter = 1; iter < maxiter && exit != 0; iter++)
        {
                if( normR < eps)
                {
                        exit = 0;
                        break;
                }
                scalPP = p * p;
                scalPR = p * r;
                        a = scalPR/scalPP;
                for (int i = 0; i < A.n; i++)
                {
                        x0.V[i] = x0.V[i] + z.V[i] * a;
                        r.V[i] = r.V[i] - p.V[i] * a;
                }
                Reverse(y, r);
                A.multMV(y, Ar);
                Direct(Ar, Ar);
                scalPAr = p * Ar;
                b = -(scalPAr/scalPP);
                for (int i = 0; i < A.n; i++)
                {
                        z.V[i] = y.V[i] + z.V[i] * b;
                        p.V[i] = Ar.V[i] + p.V[i] * b;
                }       
                scalRR = r * r;
                normR = sqrt(scalRR)/A.normF;
        }
}


void LOS::LUdec()
{
        L = A.ggl;
        U = A.ggu;
        D = A.di;
        double l, u, d;
        for(int k = 0; k < A.n; k++)
        {
                d = 0;
                int i0 = A.ig[k], i1 = A.ig[k + 1];
                int i = i0;
                for(; i0 < i1; i0++)
                {
                        l = 0;
                        u = 0;
                        int j0 = i, j1 = i0;
                        for(; j0 < j1; j0++)
                        {
                                int t0 = A.ig[A.jg[i0]], t1 = A.ig[A.jg[i0]+1];
                                for (;   t0 < t1; t0++)
                                {
                                        if (A.jg[j0] == A.jg[t0])
                                        {
                                                l += L[j0] * U[t0];
                                                u += L[t0] * U[j0];
                                        }
                                }
                        }
                        L[i0] -= l;
                        U[i0] -= u;
                        U[i0] /= D[A.jg[i0]];
                        d += L[i0] * U[i0];
                }
                D[k] -= d;
        }
}

// прямой ход Ly=F
void LOS::Direct(Vector &y, Vector &F)
{
        y = F;
        for (int i = 0; i < A.n; i++)
        {
                double sum = 0;
                int k0 = A.ig[i], k1 = A.ig[i + 1];
                int j;
                for (; k0 < k1; k0++)
                {
                        j = A.jg[k0];
                        sum += y.V[j] * L[k0];
                }
                double buf = y.V[i] - sum;
                y.V[i] = buf / D[i];
        }
}

// обратный ход Ux=y
void LOS::Reverse(Vector &x, Vector &y) 
{
        x = y;
        for(int i = A.n - 1; i >= 0; i--)
        {
                int k0 = A.ig[i], k1 = A.ig[i + 1];
                int j;
                for (; k0 < k1; k0++)
                {
                        j = A.jg[k0];
                        x.V[j] -= x.V[i]*U[k0];
                }
        }
}

void LOS::outputMatrix(std::ofstream &X)
{
	for (int i = 0; i < A.n; ++i)
	{
		for (int j = 0; j < A.n; ++j)
		{
			if (i == j)
			{
				X << A.di[i] << "\t";
			}
			else
			{
				if (i > j)
				{
					int i0 = A.ig[i];
					int i1 = A.ig[i + 1];
					for (int k = i0; k < i1; k++)
						if (A.jg[k] == j)
						{
							X << A.ggl[k] << "\t";
							goto xxx;
						}
					X << 0 << "\t";
				}
				else
				{
					int i2 = j, j2 = i;
					int i0 = A.ig[i2];
					int i1 = A.ig[i2 + 1];
					for (int k = i0; k < i1; k++)
						if (A.jg[k] == j2)
						{
							X << A.ggu[k] << "\t";
							goto xxx;
						}
					X << 0 << "\t";
				}
			}
		xxx:;
		}
		X << "\n";
	}
}

void LOS::output(std::ofstream &X)
{
        X << iter << std::endl;
        X.precision(17);
        X << normR << std::endl;
        for (int i = 0; i < A.n; i++)
        {
                X << x0.V[i] << std::endl;
        }
}

