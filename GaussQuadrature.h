#include<vector>
#include <fstream>
#pragma once
class GaussQuadrature
{
public:
	static const int n;
	static std::vector<double> weights;
	static std::vector<double> abscissa;
	static void Initialize();
	static double Integration(double(*func)(double), double minValue, double maxValue);
private:
	GaussQuadrature();
	~GaussQuadrature();
};

