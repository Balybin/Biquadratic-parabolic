#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>

struct pt
{
	double x;
	double y;
};

struct dot
{
	int x;
	int y;
};

//double lambda(pt x);
//double gamma(pt x);
double func(pt x);
double first(pt x, int choise);
double second(pt x, int choise);
double third(pt x);