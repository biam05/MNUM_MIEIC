#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double dT(double T) { return -0.25 * (T - 37); }

double f(double x) { return 2 * log(2 * x); }

double y(double x) { return sqrt(1 + pow(2.5 * exp(2.5*x),2)); }

double w(double x) { return pow(x, 3) - 10 * sin(x) + 2.8; }

void meuler(double h, double T, double t)
{
	double T1;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\tT: " << T << endl;
		t += h;
		T1 = dT(T);
		T += h * T1;
	}
}

void meliminacaodegaus()
{
	vector<vector<double>> matrix = { {1,0.5,1 / 3.0},{0.5,1 / 3.0,0.25},{1 / 3.0,0.25,0.2} }, err_matrix = matrix;
	vector<double> b = { -1,1,1 }, sol = { 0,0,0 };

	// OBTER A MATRIZ TRIANGULAR SUPERIOR
	for (int i = 0; i < 3; i++)
	{
		//div = elemento da diagonal principal de cada linha
		double div = matrix[i][i];

		//dividir a linha toda pelo elemento da diagonal princpal para ficar a 1, inclusive a resposta (b)
		for (int j = i; j < 3; j++)
		{

			matrix[i][j] /= div;
		}
		b[i] /= div;

		//subtrair as linhas seguintes a linha que coloquei a diagonal a 1
		for (int j = i + 1; j < 3; j++)
		{
			double mul = matrix[j][i];
			for (int k = 0; k < 3; k++)
			{
				matrix[j][k] -= matrix[i][k] * mul;
			}
			b[j] -= b[i] * mul;
		}
	}

	//prints
	for (int i = 0; i < 3; i++)
	{
		cout << "[";
		for (int j = 0; j < 3; j++)
		{
			cout << matrix[i][j] << "\t|";
		}
		cout << b[i] << "]" << endl;
	}


	double result;

	//começar a substituir os valores dados em b pelos valores na matriz ( z = b[2]; etc)
	for (int i = 2; i > -1; i--)
	{
		result = b[i];
		for (int j = i + 1; j < 3; j++)
		{
			result -= sol[j] * matrix[i][j];
		}
		sol[i] = result;
	}
	cout << "x1 = " << sol[0] << endl << "x2 = " << sol[1] << endl << "x3 = " << sol[2] << endl;

	//resolver tudo igual so q agora o b é o erro (estabiliidade externa)
	vector<double> err = { 0.05 - 0.05*sol[0] - 0.05*sol[1] - 0.05*sol[2] , 0.05 - 0.05*sol[0] - 0.05*sol[1] - 0.05*sol[2] , 0.05 - 0.05*sol[0] - 0.05*sol[1] - 0.05*sol[2] };
	for (int i = 0; i < 3; i++)
	{
		double div = err_matrix[i][i];
		for (int j = i; j < 3; j++)
		{
			err_matrix[i][j] /= div;
		}
		err[i] /= div;
		for (int j = i + 1; j < 3; j++)
		{
			double mul = err_matrix[j][i];
			for (int k = 0; k < 3; k++)
			{
				err_matrix[j][k] -= err_matrix[i][k] * mul;
			}
			err[j] -= err[i] * mul;
		}
	}

	//prints
	for (int i = 0; i < 3; i++)
	{
		cout << "[";
		for (int j = 0; j < 3; j++)
		{
			cout << err_matrix[i][j] << "\t|";
		}
		cout << err[i] << "]" << endl;
	}

	//substituir tudo
	for (int i = 2; i > -1; i--)
	{
		result = err[i];
		for (int j = i + 1; j < 3; j++)
		{
			result -= sol[j] * err_matrix[i][j];
		}
		sol[i] = result;
	}
	cout << "deltax1 = " << sol[0] << endl << "deltax2 = " << sol[1] << endl << "deltax3 = " << sol[2] << endl;
}

void mpicardpeano(double x)
{
	double x1;
	for (int i = 0; i < 2; i++)
	{
		cout << "iteracao: " << i << "\tx: " << x << endl;
		x1 = f(x);
		x = x1;
	}
}

double mtrapezios(double a, double b, double h)
{
	double s = y(a) + y(b);
	double n = (b - a) / h;
	for (int i = 1; i < n; i++)
	{
		s += 2 * y(a + i * h);
	}
	s *= h / 2;
	return s;
}

double msimpson(double a, double b, double h)
{
	double s = y(a) + y(b);
	double n = (b - a) / h;
	for (int i = 1; i < n; i++)
	{
		if( i % 2 == 0)
			s += 2 * y(a + i * h);
		else
			s += 4 * y(a + i * h);
	}
	s *= h / 3;
	return s;
}

void mbisseccao(double a, double b)
{
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\ta: " << a << "\tb: " << b << endl;
		double m = (a + b) / 2;
		if (w(a) * w(m) < 0)
			b = m;
		else
			a = m;
	}
}

int main()
{
	cout << "\t\tMETODO DE EULER - EX1" << endl << endl;
	meuler(0.4, 3, 5);


	cout << "\n\t\tMETODO DE ELIMINACAO DE GAUSS - EX3" << endl << endl;
	meliminacaodegaus();


	cout << "\n\t\tMETODO DE PICARD PEANO - EX4" << endl << endl;
	mpicardpeano(1.1);


	cout << "\n\t\tMETODO DE SIMPSON VS METODO DOS TRAPEZIOS - EX5" << endl << endl;

	cout << "\t\t\tMETODO DOS TRAPEZIOS" << endl << endl;

	double ht0 = mtrapezios(0, 1, 0.125);
	double ht1 = mtrapezios(0, 1, 0.125 / 2.0);
	double ht2 = mtrapezios(0, 1, 0.125 / 4.0);

	cout << "h: " << ht0 << "\th': " << ht1 << "\th'': " << ht2 << endl;

	double qct = (ht1 - ht0) / (ht2 - ht1);

	cout << "Qc: " << qct << endl;

	double et = (ht2 - ht1) / 3.0;

	cout << "E: " << et << endl;

	cout << "\n\t\t\tMETODO DE SIMPSON" << endl << endl;

	double hs0 = msimpson(0, 1, 0.125);
	double hs1 = msimpson(0, 1, 0.125 / 2.0);
	double hs2 = msimpson(0, 1, 0.125 / 4.0);

	cout << "h: " << hs0 << "\th': " << hs1 << "\th'': " << hs2 << endl;

	double qcs = (hs1 - hs0) / (hs2 - hs1);

	cout << "Qc: " << qcs << endl;

	double es = (hs2 - hs1) / 15.0;

	cout << "E: " << es << endl;


	cout << "\n\t\tMETODO DA BISSECAO - EX7" << endl << endl;
	mbisseccao(1.5, 4.2);

	return 0;
}
