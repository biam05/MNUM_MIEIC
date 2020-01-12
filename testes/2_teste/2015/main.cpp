#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double f1(double x, double y, double z, double w)
{
	return -(- y - z + w - 1) / 4.5;
}

double f2(double x, double y, double z, double w)
{
	return -(-x + z -w + 1) / 4.5;
}

double f3(double x, double y, double z, double w)
{
	return -(-x + 2 * y -w +1) / 4.5;
}

double f4(double x, double y, double z, double w)
{
	return -(2*x -y -z) / 4.5;
}

double dtdt(double T, double Ta)
{
	return -0.25*(T - Ta);
}

double dydx(double z, double y)
{
	return z;
}

double dzdx(double z, double y)
{
	return -7 * z + 2 * y;
}

double simpson(double a, double b, double h)
{
	vector<double> f = { 0.36, 1.19, 1.32, 0.21, 1.15, 1.39, 0.12, 1.22, 0.60 };
	double n = (b - a) / h;
	double s = f[a];
	for (int i = 1; i < n; i++)
	{
		if (i % 2 == 0)
			s += 2 * f[(a + i * h) / 0.25];
		else
			s += 4 * f[(a + i * h) / 0.25];
	}
	s += f[b / 0.25];

	s = (h / 3) * s;
	return s;
}

void gauss_jacobi(double x0, double y0, double z0, double w0)
{
	double x, y, z, w;
	for (int i = 0; i < 3; i++)
	{
		cout << "Iteracao: " << i << "\t|x: " << x0 << "\t|y: " << y0 << "\t|z: " << z0
			<< "\t|w: " << w0 << endl;
		x = f1(x0, y0, z0, w0);
		y = f2(x0, y0, z0, w0);
		z = f3(x0, y0, z0, w0);
		w = f4(x0, y0, z0, w0);
		x0 = x; y0 = y; z0 = z; w0 = w;
	}
}

void euler(double t, double T, double Ta, double h)
{
	double deltaT, deltaTa;
	for (int i = 0; i < 3; i++)
	{
		cout << "Iteracao: " << i << "\t|T: " << T << "\t|t: " << t << endl;
		deltaT = dtdt(T, Ta);
		t += h;
		T += deltaT * h;
	}
}

void euleredo(double z, double y, double h)
{
	double deltaZ, deltaY;
	for (int i = 0; i < 3; i++)
	{
		cout << "Iteracao: " << i << "\t|y:" << y << "\t|ylinha: " << z << endl;
		deltaZ = h * dzdx(z, y);
		deltaY = h * dydx(z, y);
		z += deltaZ;
		y += deltaY;
	}
}

int main()
{
	cout << ":::::METODO DE SIMPSON (EX1):::::" << endl;
	double S = simpson(0, 2, 1);
	double S1 = simpson(0, 2, 0.5);
	double S2 = simpson(0, 2, 0.25);
	double E = abs(S2 - S1) / 15;
	cout << "S: " << S << "\t|E: " << E << endl;
	cout << "\n:::::GAUSS-JACOBI (EX2A):::::" << endl;
	gauss_jacobi(0.25, 0.25, 0.25, 0.25);
	cout << "\n:::::METODO DE EULER (EX4):::::" << endl;
	euler(1, 23, 45, 0.4);
	cout << "\n:::::METODO DE EULER EDOS (EX5):::::" << endl;
	euleredo(1, 0, 0.25);
	return 0;
}