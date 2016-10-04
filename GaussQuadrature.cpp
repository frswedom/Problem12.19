#include "GaussQuadrature.h"

const int GaussQuadrature::n = 128;
std::vector<double> GaussQuadrature::weights = std::vector<double>(n);
std::vector<double> GaussQuadrature::abscissa = std::vector<double>(n);

void GaussQuadrature::Initialize()
{
	std::ifstream in_abscissa("D:\\Work\\NumericalMethods\\Problem12.19\\_abscissa.txt");
	std::ifstream in_weights("D:\\Work\\NumericalMethods\\Problem12.19\\_weights.txt");
	weights = std::vector<double>(n);
	abscissa = std::vector<double>(n);

	for (int i = 0; i < n; i++)
	{
		in_abscissa >> abscissa[i] ;
		in_weights >> weights[i];
	}
	in_abscissa.close();
	in_weights.close();
}

 double GaussQuadrature::Integration(double(*func)(double), double minValue, double maxValue)
{
	double res = 0.0;
	double middle = (minValue + maxValue) / 2;
	double interval = (maxValue - minValue) / 2;
	for (int i = 0; i < n; i++)
	{
		res += weights[i]*func(middle + interval * abscissa[i]);
	}
	res *= interval;
	return res;
}

GaussQuadrature::GaussQuadrature()
{
}


GaussQuadrature::~GaussQuadrature()
{
}
