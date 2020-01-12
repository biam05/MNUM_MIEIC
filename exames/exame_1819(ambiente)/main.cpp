#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double f1(double x, double y) { return exp(x - y) - sin(x+y); }
double f1x(double x, double y) { return exp(x - y) - cos(x + y); }
double f1y(double x, double y) { return -exp(x - y) - cos(x + y); }
double f2(double x, double y) { return pow(x, 2) * pow(y, 2) - cos(x + y); }
double f2x(double x, double y) { return sin(y + x) + 2 * x * pow(y, 2); }
double f2y(double x, double y) { return sin(y + x) + 2 * y * pow(x, 2); }

double f(double x, double y)
{
	if (x == 0.0)
	{
		if (y == 0.0) return 1.1;
		else if (y == 1.0) return 2.1;
		else if (y == 2.0) return 7.8;
	}
	else if (x == 1.0)
	{
		if (y == 0.0) return 1.4;
		else if (y == 1.0) return 4.0;
		else if (y == 2.0) return 1.5;
	}
	else if (x == 2.0)
	{
		if (y == 0.0) return 9.8;
		else if (y == 1.0) return 2.2;
		else if (y == 2.0) return 1.2;
	}
	return 0;
}

double dydx(double x, double y, double z) { return z; }
double dzdx(double x, double y, double z) { return -7 * z - 4 * y; }

void mnewton_sistemas(double x, double y)
{
	double jacobian, hn, kn;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tx: " << x << "\ty: " << y << endl;
		jacobian = f1x(x, y) * f2y(x, y) - f1y(x, y) * f2x(x, y);
		hn = -(f1(x, y) * f2y(x, y) - f1y(x, y) * f2(x, y)) / jacobian;
		kn = -(f1x(x, y) * f2(x, y) - f1(x, y) * f2x(x, y)) / jacobian;

		x += hn;
		y += kn;
	}
}

void msimpson_cubatura(double a0, double b0, double a1, double b1)
{
	double hx = (a1 - a0) / 2.0;
	double hy = (b1 - b0) / 2.0;

	double s = f(a0, b0) + f(a0, b1) + f(a1, b0) + f(a1, b1);
	s += 4 * (f(a0 + hx, b0) + f(a0 + hx, b1) + f(a0, b0 + hy) + f(a1, b0 + hy));
	s += 16 * f(a0 + hx, b0 + hy);

	s *= hx * hy / 9.0;
	cout << "S: " << s << endl;
}

void meuler(double x, double y, double z, double h)
{
	double yprev;
	for (int i = 0; i < 4; i++)
	{
		yprev = y;
		cout << "iteracao: " << i << "\tx: " << x << "\ty: " << y << "\ty': " << z << endl;
		y += h * dydx(x, y, z);
		z += h * dzdx(x, yprev, z);
		x += h;
	}
}

int main()
{
	cout << "\t\tMETODO DE NEWTON SISTEMAS - EX1" << endl << endl;
	mnewton_sistemas(0.5, 0.25);

	cout << "\n\t\tCUBATURA PELO METODO DE SIMPSON - EX3" << endl << endl;
	msimpson_cubatura(0, 0, 2, 2);

	cout << "\n\t\tMETODO DE EULER - EX4" << endl << endl;
	meuler(0.4, 2.0, 1.0, 0.2);
	return 0;
}