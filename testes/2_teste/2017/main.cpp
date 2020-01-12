#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double f1(double x, double y, double z, double w)
{
	return -(0.5*y + 3 * z + 0.25*w - 25) / 6;
}

double f2(double x, double y, double z, double w)
{
	return -(1.2*x + 0.25*z + 0.20*w - 10) / 3;
}

double f3(double x, double y, double z, double w)
{
	return -(-1 * x + 0.25*y + 2 * w - 7) / 4;
}

double f4(double x, double y, double z, double w)
{
	return -(2 * x + 4 * y + 1 * z + 12) / 8;
}

double dzdt(double t, double z)
{
	return 2 + pow(t, 2) + t * z;
}

double dydt(double t, double z)
{
	return z;
}

void eliminacao_ext(double x1, double y1, double z1)
{
	vector<vector<double>> matrix = {
   {18,-1,1, 0.1 * (1 - 0.552949 + 0.15347 + 0.10655)},
   {3,-5,4, 0.1 * (1 - 0.552949 + 0.15347 + 0.10655)},
   {6,8,29, 0.1 * (1 - 0.552949 + 0.15347 + 0.10655)}
	};

	double pivot;

	for (int i = 0; i < matrix.size(); i++) {
		pivot = matrix[i][i];
		for (int j = i; j < matrix[0].size(); j++) {
			matrix[i][j] /= pivot;
		}
		for (int j = i + 1; j < matrix.size(); j++) {
			pivot = -matrix[j][i];
			for (int k = i; k < matrix[0].size(); k++)
				matrix[j][k] += matrix[i][k] * pivot;
		}
	}

	double z = matrix[2][3] / matrix[2][2];
	double y = (matrix[1][3] - z * matrix[1][2]) / matrix[1][1];
	double x = (matrix[0][3] - z * matrix[0][2] - y * matrix[0][1]) / matrix[0][0];
	cout << "x: " << x << "\t| y: " << y << "\t| z: " << z << endl;
}

void gauss_seidel(double x, double y, double z, double w)
{
	x = f1(x, y, z, w);
	y = f2(x, y, z, w);
	z = f3(x, y, z, w);
	w = f4(x, y, z, w);

	cout << "x: " << x << "\t|y: " << y << "\t|z: " << z << "\t|w: " << w << endl;
}

double simpson(double a, double b, double h)
{
	vector<double> f = { 1.04,0.37,0.38,1.49,1.08,0.13,0.64,0.84,0.12 };
	double s = f[a];
	double n = (b - a) / h;
	for (int i = 1; i < n; i++)
	{
		if (i % 2 == 0)
			s += 2 * f[(a + i * h) / 0.25];
		else
			s += 4 * f[(a + i * h) / 0.25];
	}
	s += f[b / 0.25];

	s = (h / 3)*s;
	return s;
}

double cubatura_trapezios()
{
	vector<vector<double>> xy = { {1.1,1.4,7.7},{2.1,3.1,2.2},{7.3,1.5,1.2} };
	double s = 0;
	for (int i = 0; i < xy.size(); i++)
	{
		for (int j = 0; j < xy.size(); j++)
		{
			if ((i == 0 || i == xy.size() - 1) && (j == 0 || j == xy.size() - 1))
				s += xy[i][j];
			
			else if (i == 0 || j == 0 || i == xy.size() - 1 || j == xy.size() - 1)
				s += 2 * xy[i][j];
			
			else
				s += 4 * xy[i][j];
			
		}
	}
	s /= 4;
	cout << "S: " << s << endl;
	return 0;
}

void euler(double t, double y, double z, double h)
{
	cout << "\tEULER" << endl;
	double deltaz;
	for (int i = 0; i < 3; i++)
	{
		cout << "It: " << i << "\t|t: " << t << "\t|y: " << y << endl;
		deltaz = dzdt(t, z);
		t += h;
		y += h * z;
		z += deltaz * h;
	}
}

void rk4(double t, double y, double z, double h)
{
	cout << "\tRK4" << endl;
	double y1, y2, y3, y4, z1, z2, z3, z4;
	for (int i = 0; i < 3; i++)
	{
		cout << "It: " << i << "\t|t: " << t << "\t|y: " << y << endl;
		y1 = h * dydt(t, z);
		z1 = h * dzdt(t, z);
		y2 = h * dydt(t + h / 2, z + z1 / 2);
		z2 = h * dzdt(t + h / 2, z + z1 / 2);
		y3 = h * dydt(t + h / 2, z + z2 / 2);
		z3 = h * dzdt(t + h / 2, z + z2 / 2);
		y4 = h * dydt(t + h, z + z3);
		z4 = h * dzdt(t + h, z + z3);
		y += y1 / 6 + y2 / 3 + y3 / 3 + y4 / 6;
		z += z1 / 6 + z2 / 3 + z3 / 3 + z4 / 6;
		t += h;
	}

}

int main()
{
	cout << ":::::ESTABILIDADE EXTERNA (EX1):::::" << endl;
	eliminacao_ext(0,0,0);
	cout << "\n:::::GAUSS SEIDEL (EX2):::::" << endl;
	gauss_seidel(2.12687,2.39858,3.99517,-3.73040);
	cout << "\n:::::SIMPSON (EX3):::::" << endl;
	double S = simpson(0, 2.0, 1);
	double S1 = simpson(0, 2.0, 0.5);
	double S2 = simpson(0, 2.0, 0.25);
	cout << "S: " << S;
	cout << endl << "E: " << abs(S2 - S1) / 15.0 << endl;
	cout << "\n:::::CUBATURA_TRAPEZIOS (EX4):::::" << endl;
	cubatura_trapezios();
	cout << "\n:::::EULER VS RK4 (EX6):::::" << endl;
	euler(1, 1, 0, 0.25);
	rk4(1, 1, 0, 0.25);
	return 0;
}