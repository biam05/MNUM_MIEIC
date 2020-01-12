#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double f1(double x, double y, double z, double w)
{
	return -(-y - z + w - 1) / 4.5;
}

double f2(double x, double y, double z, double w)
{
	return -(-x + z - w + 1) / 4.5;
}

double f3(double x, double y, double z, double w)
{
	return -(-x + 2 * y - w + 1) / 4.5;
}

double f4(double x, double y, double z, double w)
{
	return -(2 * x - y - z) / 4.5;
}

double dxdt(double x, double t)
{
	return sin(2 * x) + sin(2 * t);
}



void gauss_jacobi(double x0, double y0, double z0, double w0)
{
	double x, y, z, w;
	for (int i = 0; i < 3; i++)
	{
		cout << "It: " << i << "\t|x: " << x0 << "\t|y: " << y0 << "\t|z: "
			<< z0 << "\t|w: " << w0 << endl;
		x = f1(x0, y0, z0, w0);
		y = f2(x0, y0, z0, w0);
		z = f3(x0, y0, z0, w0);
		w = f4(x0, y0, z0, w0);
		x0 = x; y0 = y; z0 = z; w0 = w;
	}
}

double rk4edo(double x, double t, double h)
{
	double x1, x2, x3, x4;
	while(x < 1.5)
	{
		cout << "x: " << t << "\t|t: " << x << endl;
		x1 = h * dxdt(x, t);
		x2 = h * dxdt(x + h / 2, t + x1 / 2);
		x3 = h * dxdt(x + h / 2, t + x2 / 2);
		x4 = h * dxdt(x + h, t + x3);

		t += x1 / 6 + x2 / 3 + x3 / 3 + x4 / 6;
		x += h;

	}
	return t;
}

double simpson(double a, double b, double h)
{

	return 0;
}

int main()
{
	cout << ":::::GAUSS-JACOBI (EX1A):::::" << endl;
	gauss_jacobi(0.25, 0.25, 0.25, 0.25);
	cout << "\n:::::RK4 (EX2):::::" << endl;
	double rk4a = rk4edo(1, 1, 0.125);
	double rk4 = rk4edo(1, 1, 0.5/8);
	double rk41 = rk4edo(1, 1, 0.5/16);
	double rk42 = rk4edo(1, 1, 0.5/32);
	double qc = (rk41 - rk4) / (rk42 - rk41);
	cout << "Qc: " << qc << endl;
	cout << "\n:::::SIMPSON VS TRAPEZIOS (EX2):::::" << endl;



	return 0;
}