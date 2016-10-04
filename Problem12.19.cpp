#include <iostream>
#include <cmath>
#include "GaussQuadrature.h"

using namespace std;
using uint = size_t;
double f1(double x)
{
	return sqrt(0.4*x*x + 1.5) / (2.5 + sqrt(2 * x + 0.8));
}
double f2(double x)
{
	return tan(x*x + 0.5) / (1 + 2*x*x);
}double f3(double x)
{
	return (1+0.3*x*x) / (0.8+sqrt(0.6*x*x+1.3));
}double f4(double x)
{
	return x*x/sqrt(x+1);
}
double SimpsonIntegration(double (*func)(double), double minValue, double maxValue, uint steps)
{
	double res{};
	if (steps % 2)
	{
		steps++;
	}
	double h = (maxValue - minValue) / steps;
	res = func(minValue) + func(maxValue);
	double x{ minValue + h };
	for (uint i{1}; i < steps; i++, x += h)
	{
		if (i % 2) 
		{
			res += 4 * func(x);
		}
		else
		{
			res += 2 * func(x);
		}
	}
	res *= h / 3;
	return res;
}

double ChebyshevIntegration(double(*func)(double), double minValue, double maxValue)
{
	// n = 9
	static const double roots[] {0.167906184214804, 0.528761783057880, 0.601018655380238, 0.911589307728434};
	double middle = (minValue + maxValue) / 2;
	double interval = (maxValue - minValue) / 2;
	double res = func(middle);
	for (uint i = 0; i < 4; i++)
	{
		res += func(middle + interval * roots[i]);
		res += func(middle - interval * roots[i]);
	}
	res *= 2 * interval / 9;
	return res;
}


int main()
{
	double eps;
	cout << "Enter presicious: \n";
	cin >> eps;
	cout.precision(16);
	GaussQuadrature::Initialize();
	cout << SimpsonIntegration(f1, 0.8, 2.4, 10000) << endl;
	cout << ChebyshevIntegration(f1, 0.8, 2.4) << endl;
	cout << GaussQuadrature::Integration(f1, 0.8, 2.4) << endl;
	cout << endl;
	cout << SimpsonIntegration(f2, 0.4, 0.8, 10000) << endl; 
	cout << ChebyshevIntegration(f2, 0.4, 0.8) << endl;
	cout << GaussQuadrature::Integration(f2, 0.4, 0.8) << endl;
	cout << endl;
	cout << SimpsonIntegration(f3, 0.4, 2.6, 10000) << endl;
	cout << ChebyshevIntegration(f3, 0.4, 2.6) << endl;
	cout << GaussQuadrature::Integration(f3, 0.4, 2.6) << endl;
	cout << endl;
	cout << SimpsonIntegration(f4, 2.4, 3.2, 10000) << endl;
	cout << ChebyshevIntegration(f4, 2.4, 3.2) << endl;
	cout << GaussQuadrature::Integration(f4, 2.4, 3.2) << endl;

	return 0;
}

