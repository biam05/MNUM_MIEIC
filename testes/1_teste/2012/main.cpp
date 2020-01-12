#include <cmath>
#include <iostream>

using namespace std;

double f(double x)
{
	return pow(x, 3) + 2 * pow(x, 2) + 10 * x - 17;
}
double df(double x)
{
	return 3 * pow(x, 2) + 4 * x + 10;
}
void mnewton(double x)
{
	double xseguinte = x;
	for (int i = 0; i <= 3; i++)
	{
		x = xseguinte;
		xseguinte = x - f(x) / df(x);
		cout << "x = " << x << endl;
	}
}

//#############################################//

double f1(double x, double y) { return y - log(x - 1); }
double df1x(double x, double y) { return -1 / (x - 1); }
double df1y(double x, double y) { return 1; }
double f2(double x, double y) { return pow(y,2) + pow(x - 3, 2) - pow(2,2); }
double df2x(double x, double y) { return 2 * (x - 3); }
double df2y(double x, double y) { return 2 * y; }

void mnewtonsistemas(double x, double y)
{
	
	for (int i = 0; i <= 3; i++)
	{
		double jacobian = df1x(x, y) * df2y(x, y) - df2x(x, y) * df1y(x, y);
		double hn = f1(x, y) * df2y(x, y) - df1y(x, y) * f2(x, y);
		double kn = df1x(x, y) * f2(x, y) - df2x(x, y) * f1(x, y);
		x -= hn / jacobian;
		y -= kn / jacobian;

		cout << "x = " << x << " | y = " << y << endl;
	}
}

int main()
{
	//mnewton(0);
	mnewtonsistemas(1.50000, 1.30000);
	return 0;
}