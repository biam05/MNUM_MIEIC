#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double f(double x) { return (x - 2.6) + pow(cos(x + 1.1), 3); }
double dfx(double x) { return 1 - 3 * pow(cos(x + 1.1), 2)*sin(x + 1.1); }

double y(double x) { return 5 * cos(x) - sin(x); }

double dxt(double x, double t) { return sin(1 * x) + sin(2 * t); }

void mnewton(double x)
{
	for (int i = 0; i < 2; i++)
	{
		cout << "iteracao: " << i << "\tx: " << x << endl;
		x -= f(x) / dfx(x);
	}
}

void meliminacaodegauss()
{
	vector<vector<double>> matrix = { {0.1,0.5,3.0,0.25},{1.2,0.2,0.25,0.2},{-1,0.25,0.3,2},{2,0.00001,2,0.4} }, err_matrix = matrix;
	vector<double> b = { 0,1,2,3 };
	vector<double> sol = { 0,0,0,0 };
	double result;

	for (int i = 0; i < 4; i++) //filas
	{
		double div = matrix[i][i];
		for (int j = i; j < 4; j++) //colunas
		{
			matrix[i][j] /= div;
		}
		b[i] /= div;

		for (int j = i + 1; j < 4; j++)
		{
			double mul = matrix[j][i];
			for (int k = 0; k < 4; k++)
			{
				matrix[j][k] -= mul * matrix[i][k];
			}
			b[j] -= mul * b[i];
		}

	}
	for (int i = 0; i < 4; i++)
	{
		cout << "[";
		for (int j = 0; j < 4; j++)
		{
			cout << matrix[i][j] << "|";
		}
		cout << "|" << b[i];
		cout << "]" << endl;
	}

	for (int i = 3; i >= 0; i--)
	{
		result = b[i];
		for (int j = i; j < 4; j++)
		{
			result -= sol[j] * matrix[i][j];
		}
		sol[i] = result;
	}
	cout << "\nx1: " << sol[0] << "\nx2: " << sol[1] << "\nx3: " << sol[2] << "\nx4: " << sol[3] << endl;

	vector<double> err = { 0.05 - 0.05*sol[0] - 0.05*sol[1] - 0.05*sol[2] - 0.05*sol[3] , 0.05 - 0.05*sol[0] - 0.05*sol[1] - 0.05*sol[2] - 0.05*sol[3] , 0.05 - 0.05*sol[0] - 0.05*sol[1] - 0.05*sol[2] - 0.05*sol[3], 0.05 - 0.05*sol[0] - 0.05*sol[1] - 0.05*sol[2] - 0.05*sol[3] };
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
	cout << endl;
	//prints
	for (int i = 0; i < 4; i++)
	{
		cout << "[";
		for (int j = 0; j < 4; j++)
		{
			cout << err_matrix[i][j] << "\t|";
		}
		cout << err[i] << "]" << endl;
	}

	//substituir tudo
	for (int i = 3; i > -1; i--)
	{
		result = err[i];
		for (int j = i + 1; j < 4; j++)
		{
			result -= sol[j] * err_matrix[i][j];
		}
		sol[i] = result;
	}
	cout << "\ndeltax1 = " << sol[0] << endl << "deltax2 = " << sol[1] << endl << "deltax3 = " << sol[2] << endl << "deltax4 = " << sol[3] << endl;
}

void mseccaoaurea(double x1, double x2)
{
	double B = (sqrt(5) - 1) / 2.0;
	double A = pow(B, 2);
	double x3 = x1 + A * (x2 - x1);
	double x4 = x1 + B * (x2 - x1);
	for (int i = 0; i < 4; i++)
	{
		cout << "iteracao: " << i << "\tx1: " << x1 << "\tx2: " << x2 << "\tx3: " << x3 << "\tx4: " << x4 << endl
			<< "\tf(x1): " << y(x1) << "\tf(x2): " << y(x2) << "\tf(x3): " << y(x3) << "\tf(x4): " << y(x4) << endl;
		cout << "amplitude: " << x2 - x1 << endl << endl;

		if (y(x3) < y(x4))
		{
			x2 = x4;
			x4 = x3;
			x3 = x1 + A * (x2 - x1);
		}
		else
		{
			x1 = x3;
			x3 = x4;
			x4 = x1 + B * (x2 - x1);
		}
	}
}

double mrk4(double x, double t, double h)
{
	double xprev, o1, o2, o3, o4, deltax;
	while (t <= 1.5)
	{
		xprev = x;
		cout << "t: " << t << "\tx: " << x << endl;
		o1 = h * dxt(x, t);
		o2 = h * dxt(x + o1 / 2, t + h / 2);
		o3 = h * dxt(x + o2 / 2, t + h / 2);
		o4 = h * dxt(x + o3, t + h);
		deltax = o1 / 6 + o2 / 3 + o3 / 3 + o4 / 6;
		x += deltax;
		t += h;
	}
	return xprev;
}

double mtrapezios(double a, double b, double h)
{
	vector<double> f = { 5, 5.1, 5.6, 5.9, 6.2, 7, 7.8, 8, 8.5 };
	double s = f[0];
	double n = (b - a) / h;
	for (int i = 1; i < n; i++)
	{
		s += 2 * f[(i * h) / 0.1];
	}
	s+= f[(b-a) / 0.1];

	s *= h / 2;
	
	return s;
}

int main()
{
	cout << "\t\tMETODO DE NEWTON - EX1" << endl << endl;
	mnewton(1.8);
	cout << "\n\t\tMETODO DE ELIMINACAO DE GAUSS - EX3" << endl << endl;
	meliminacaodegauss();
	cout << "\n\t\tMETODO DA SECCAO AUREA - EX4" << endl << endl;
	mseccaoaurea(2, 4);
	cout << "\n\t\tMETODO DE RK4 - EX5" << endl << endl;
	double s2 = mrk4(1, 1, 0.125);
	double qc = (1.768150 - 1.767816) / (s2 - 1.768150);
	cout << "\nQC: " << qc << endl;
	cout << "\n\t\tMETODO DOS TRAPEZIOS - EX6" << endl << endl;
	double t = mtrapezios(1, 1.8, 0.4);
	double t1 = mtrapezios(1, 1.8, 0.2);
	double t2 = mtrapezios(1, 1.8, 0.1);
	cout << "I'': " << t2 << "\tI': " << t1 << "\tI: " << t << endl;
	double qct = (t1 - t) / (t2 - t1);
	cout << "Qc: " << qct << endl;
	double e = (t2 - t1) / 3.0;
	cout << "E: " << e << endl;
	return 0;
}