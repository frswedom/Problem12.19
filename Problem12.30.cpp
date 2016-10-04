#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;
const int n = 128;
vector<double> weights(n);
vector<double> abscissa(n);

double f1(double x)
{
	return sin(0.6*x*x+0.3) / (2.4 + cos(x + 0.5));
}
double f2(double x)
{
	return 2.0 / sqrt(1.5 * x * x + 0.7);
}
double f3(double x)
{
	return (1 + 0.6*x*x) / (2.5 + sqrt(0.3*x*x + 1.6));
}
double f4(double x)
{
	return (x+2.2) / sqrt(x * x + 1);
}

void Initialize()//файлы, из которых считываются координаты точек и веса
{
	ifstream in_abscissa("_abscissa.txt");
	ifstream in_weights("_weights.txt");

	for (int i = 0; i < n; i++)
	{
		in_abscissa >> abscissa[i];
		in_weights >> weights[i];
	}

	in_abscissa.close();
	in_weights.close();
}

double SimpsonIntegration(double (*func)(double), double minValue, double maxValue, int steps)//интегрирование с помощью формулы Симпсона
/*согласно формуле Симпсона интеграл на одном промежутке от х(і) до х(і+1) равен h/3*(y(0)+4y(1)+y(2),
 тогда весь интергал равен h/3*(func(minValue))+sum(4*y(2n+1)+2*y(2n))+func(maxValue)), 
 где y(2n+1) - сумма всех у с нечетными индексами, а y(2n) - сумма всех у с четными индексами.*/
	double res = 0.0;
	if (steps % 2)
	{
		steps++;
	}
	double h = (maxValue - minValue) / steps;
	res = func(minValue) + func(maxValue);
	double x = minValue + h;
	for (int i = 1; i < steps; i++, x += h)
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
double SimpsonIntegrationWithPrecision(double (*func)(double), double minValue, double maxValue, double eps)//интегрирование с помощью формулы Симпсона с определенной точностью
{
	//интегрируем с различным кол-вом шагов пока разность между полученными интегралами не будет меньше точности
	int k = 16;
	double res_1 = SimpsonIntegration(func, minValue, maxValue, k);
	k *= 2;
	double res_2 = = SimpsonIntegration(func, minValue, maxValue, k);

	while (abs(res_1 - res_2) > eps)
	{
		k *= 2;
		res_1 = res_2;
		res_2 = SimpsonIntegration(func, minValue, maxValue, k);
	}
	return res_2;
}

double ChebyshevIntegration(double (*func)(double), double minValue, double maxValue)//интегрирование с помощью формул Чебышева
/*интеграл с пределами [-1;1] равняется сумме значений функции в точках(узлах) с определенными коэффициентами(весами)
узлы сетки находим, а веса у всех точек одинаковы
меняем пределы интегрировани с данных нам на [-1;1] заменой переменных
*/
	// n = 9
	static const double roots[]{ 0.167906184214804, 0.528761783057880, 0.601018655380238, 0.911589307728434 };
	double middle = (minValue + maxValue) / 2;
	double interval = (maxValue - minValue) / 2;
	double res = func(middle);
	for (int i = 0; i < 4; i++)
	{
		res += func(middle + interval * roots[i]);
		res += func(middle - interval * roots[i]);
	}
	res *= 2 * interval / 9;
	return res;
}

double GaussIntegration(double(*func)(double), double minValue, double maxValue)// интегрирование с помощью формул Гаусса
/*в этом методе выбираются и веса и точки оптимальным путем
координаты точек - нули полинома Лежандра
*/
{
	double res = 0.0;
	double middle = (minValue + maxValue) / 2;
	double interval = (maxValue - minValue) / 2;
	for (int i = 0; i < n; i++)
	{
		res += weights[i] * func(middle + interval * abscissa[i]);
	}
	res *= interval;
	return res;
}

int main()//подсчет и вывод результатов на экран, проверка правильности нахождения интеграла производится с помощью взятия каждого интеграла тремя способами
{
	double eps;
	cout << "Enter presicious: \n";
	cin >> eps;
	cout.precision(16);
	Initialize();
	cout << SimpsonIntegrationWithPrecision(f1, 0.3, 1.1, eps) << endl;
	cout << ChebyshevIntegration(f1, 0.3, 1.1) << endl;
	cout << GaussIntegration(f1, 0.3, 1.1) << endl;
	cout << endl;
	cout << SimpsonIntegrationWithPrecision(f2, 1.4, 2.6, eps) << endl;
	cout << ChebyshevIntegration(f2, 1.4, 2.6) << endl;
	cout << GaussIntegration(f2, 1.4, 2.6) << endl;
	cout << endl;
	cout << SimpsonIntegrationWithPrecision(f3, 0.5, 2.3, eps) << endl;
	cout << ChebyshevIntegration(f3, 0.5, 2.3) << endl;
	cout << GaussIntegration(f3, 0.5, 2.3) << endl;
	cout << endl;
	cout << SimpsonIntegrationWithPrecision(f4, 0.4, 1.7, eps) << endl;
	cout << ChebyshevIntegration(f4, 0.4, 1.7) << endl;
	cout << GaussIntegration(f4, 0.4, 1.7) << endl;

	return 0;
}

