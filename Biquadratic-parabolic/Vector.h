#pragma once
#include "func.h"

class Vector
{
	public:
        std::vector<double> V;
        Vector(void);
        Vector(int n);
        void read(std::ifstream &vect);
        void make(int n);
        Vector& operator=(const Vector& newVect);
        double operator* (const Vector& newVect);
        Vector& operator-(Vector& newVect);
        Vector& operator+(Vector& newVect); 
        double norm();
	~Vector(void);
};

