#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double g(double x) { return pow(4 * pow(x, 3) - x + 1, 0.25); }

double y(double x) { return x + pow(x - 2, 2) / (sin(x) + 4); }

double f(double x) { return -x + 60 * cos(sqrt(x)) + 2; }
double df(double x) { return -(30 * sin(sqrt(x))) / (sqrt(x)) - 1; }

void mpicardpeano(double x)
{
	double x1;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tx: " << x << endl;
		x1 = g(x);
		x = x1;
	}
}

void meliminacaodegauss()
{
	vector<vector<double>> matrix = { {0.1,0.5,3.0,0.25},{1.2,0.2,0.25,0.2},{-1.0,0.25,0.3,2.0},{2.0,0.00001,1.0,0.4} }, err_matrix = matrix;
	vector<double> b = { 0,1.0,2.0,3.0 };
	vector<double> sol = { 0,0,0,0 };

	for (int i = 0; i < 4; i++)
	{
		double div = matrix[i][i];
		for (int j = i; j < 4; j++)
		{
			matrix[i][j] /= div;
		}
		b[i] /= div;

		for (int j = i + 1; j < 4; j++)
		{
			double mul = matrix[j][i];
			for (int k = 0; k < 4; k++)
			{
				matrix[j][k] -= matrix[i][k] * mul;
			}
			b[j] -= b[i] * mul;
		}	
	}

	for (int i = 0; i < 4; i++)
	{
		cout << "[";
		for (int j = 0; j < 4; j++)
		{
			cout << matrix[i][j] << "| ";
		}
		cout << b[i] << "]" << endl;
	}

	double result;

	for (int i = 3; i > -1; i--)
	{
		result = b[i];
		for (int j = i + 1; j < 4; j++)
		{
			result -= sol[j] * matrix[i][j];
		}
		sol[i] = result;
	}
	cout << "\nx1: " << sol[0] << "\tx2: " << sol[1] << "\tx3: " << sol[2] << "\tx4: " << sol[3] << endl << endl;

	vector<double> err =
	{ 0.3 - 0.3*sol[0] - 0.3 * sol[1] - 0.3 * sol[2] - 0.3 * sol[3],
	0.3 - 0.3*sol[0] - 0.3 * sol[1] - 0.3 * sol[2] - 0.3 * sol[3],
	0.3 - 0.3*sol[0] - 0.3 * sol[1] - 0.3 * sol[2] - 0.3 * sol[3],
	0.3 - 0.3*sol[0] - 0.3 * sol[1] - 0.3 * sol[2] - 0.3 * sol[3] };

	for (int i = 0; i < 4; i++)
	{
		double div = err_matrix[i][i];
		for (int j = i; j < 4; j++)
		{
			err_matrix[i][j] /= div;
		}
		err[i] /= div;

		for (int j = i + 1; j < 4; j++)
		{
			double mul = err_matrix[j][i];
			for (int k = 0; k < 4; k++)
			{
				err_matrix[j][k] -= err_matrix[i][k] * mul;
			}
			err[j] -= err[i] * mul;
		}
	}

	for (int i = 0; i < 4; i++)
	{
		cout << "[";
		for (int j = 0; j < 4; j++)
		{
			cout << err_matrix[i][j] << "| ";
		}
		cout << err[i] << "]" << endl;
	}
	sol = { 0,0,0,0 };
	result = 0;
	for (int i = 3; i > -1; i--)
	{
		result = err[i];
		for (int j = i + 1; j < 4; j++)
		{
			result -= sol[j] * err_matrix[i][j];
		}
		sol[i] = result;
	}
	cout << "\nox1: " << sol[0] << "\tox2: " << sol[1] << "\tox3: " << sol[2] << "\tox4: " << sol[3] << endl << endl;

}

void mseccaoaurea(double x1, double x2)
{
	double B = (sqrt(5) - 1) / 2.0;
	double A = pow(B, 2);
	double x3 = A * (x2 - x1) + x1;
	double x4 = B * (x2 - x1) + x1;
	for (int i = 0; i < 3; i++)
	{
		cout << endl << "x1: " << x1 << "\tx2: " << x2 << "\tx3: " << x3 << "\tx4: " << x4 << endl <<
			"\tf(x1): " << y(x1) << "\tf(x2): " << y(x2) << "\tf(x3): " << y(x3) << "\tf(x4): " << y(x4) << endl;

		if (y(x3) < y(x4))
		{
			x2 = x4;
			x4 = x3;
			x3 = B * (x4 - x1) + x1;
		}
		else
		{
			x1 = x3;
			x3 = x4;
			x4 = B * (x2 - x3) + x3;
		}
	}
}

void mnewton(double x)
{
	double x1;
	for (int i = 0; i < 3; i++)
	{
		cout << "xn: " << x << "\tg(xn): " << f(x) << endl;
		x1 = x - f(x) / df(x);
		x = x1;
	}
}

int main()
{
	cout << "\t\tMETODO DE PICARD PEANO - EX1" << endl << endl;
	mpicardpeano(4.0);

	
	cout << "\n\t\tMETODO DE ELIMINACAO DE GAUSS - EX2" << endl << endl;
	meliminacaodegauss();


	cout << "\n\t\tMETODO DA SECCAO AUERA - EX6" << endl << endl;
	mseccaoaurea(-1, 1.5);

	
	cout << "\n\t\tMETODO DE NEWTON - EX7" << endl << endl;
	mnewton(1.8);

	return 0;
}