#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double f(double x) { return 2 * log(2 * x); }

double dxt(double x, double t) { return sin(x) + sin(2 * t); }

double z(double x, double y) { return 6 * pow(x, 2) - x * y + 12 * y + pow(y, 2) - 8 * x; }
double zdx(double x, double y) { return -y + 12 * x - 8; }
double zdy(double x, double y) { return 2 * y - x + 12; }

void mpicardpeano(double x)
{
	for (int i = 0; i < 2; i++)
	{
		cout << "iteracao: " << i << "\tx: " << x << endl;
		x = f(x);
	}
}

double mrk4(double t, double h, double x)
{
	double o1, o2, o3, o4, deltax, xprev;
	while (t <= 1.5)
	{
		xprev = x;

		cout << "t: " << t << "\tx:" << x << endl;
		o1 = h * dxt(x, t);
		o2 = h * dxt(x + o1 / 2, t + h / 2);
		o3 = h * dxt(x + o2 / 2, t + h / 2);
		o4 = h * dxt(x + o3, t + h);

		deltax = o1 / 6.0 + o2 / 3.0 + o3 / 3.0 + o4 / 6.0;

		x += deltax;
		t += h;
	}
	return xprev;
}

vector<double> gradiente(double x, double y)
{
	vector<double> gradiente = { 0, 0 };
	gradiente[0] = zdx(x, y);
	gradiente[1] = zdy(x, y);
	return gradiente;
}

void mgradiente(double h, double x, double y)
{
	vector<double> Xn = { 1000, 1000 };
	vector<double> Xn1 = { x,y };
	for (int i = 0; i < 2; i++)
	{
		cout << "iteracao: " << i << "\tx: " << Xn1[0] << "\ty: " << Xn1[1] << "\tf(x,y): " << z(Xn1[0], Xn1[1]) <<
			"\tgradientex: " << zdx(Xn1[0], Xn1[1]) << "\tgradientey: " << zdy(Xn1[0], Xn1[1]) << endl;

		Xn[0] = Xn1[0];
		Xn[1] = Xn1[1];
		Xn1[0] = Xn[0] - h * gradiente(Xn[0], Xn[1])[0];
		Xn1[1] = Xn[1] - h * gradiente(Xn[0], Xn[1])[1];

		if (z(Xn1[0], Xn1[1]) < z(Xn[0], Xn[1]))
		{
			continue;
		}
		else if (z(Xn1[0], Xn1[1]) > z(Xn[0], Xn[1]))
		{
			h = h / 2;
			Xn1[0] = Xn[0] - h * gradiente(x, y)[0];
			Xn1[1] = Xn[0] - h * gradiente(x, y)[1];
		}
	}


}

void meliminacaodegauss()
{
	vector<vector<double>> matrix_err = { {18,-1,1},{3,-5,4},{6,8,29} };
	vector<double> sol = { 0,0,0 };
	vector<double> res= { 0.552949,-0.15347,-0.10655 };
	vector<double> err = { 0.1 - 0.1*res[0] - 0.1*res[1] - 0.1*res[2],0.1 - 0.1*res[0] - 0.1*res[1] - 0.1*res[2] ,0.1 - 0.1*res[0] - 0.1*res[1] - 0.1*res[2] };

	for (int i = 0; i < 3; i++)
	{
		double div = matrix_err[i][i];
		for (int j = i; j < 3; j++)
		{
			matrix_err[i][j] /= div;
		}
		err[i] /= div;

		for (int j = i + 1; j < 3; j++)
		{
			double mul = matrix_err[j][i];
			for (int k = 0; k < 3; k++)
			{
				matrix_err[j][k] -= mul * matrix_err[i][k];
			}
			err[j] -= err[i] * mul;
		}
	}

	double result;

	for (int i = 2; i > -1; i--)
	{
		result = err[i];
		for (int j = i + 1; j < 3; j++)
		{
			result -= sol[j] * matrix_err[i][j];
		}
		sol[i] = result;
	}

	cout << "deltax1 = " << sol[0] << endl << "deltax2 = " << sol[1] << endl << "deltax3 = " << sol[2] << endl;
}

int main()
{
	cout << "\t\tMETODO PICARD PEANO - EX1" << endl << endl;
	mpicardpeano(0.9);


	cout << "\n\t\tMETODO RK4 - EX3" << endl << endl;

	double s1 = mrk4(1, 0.5, 0);
	cout << endl;
	double s2 = mrk4(1, 0.5 / 2.0, 0);
	cout << endl;
	double s3 = mrk4(1, 0.5 / 4.0, 0);
	cout << "\nS: " << s1 << "\tS': " << s2 << "\tS'': " << s3;
	cout << "\nQc: " << (s2 - s1) / (s3 - s2);


	cout << "\n\t\tMETODO DO GRADIENTE - EX5" << endl << endl;
	mgradiente(0.25, 0, 0);


	cout << "\n\t\tESTUDO DA ESTABILIDADE EXTERNA - EX6" << endl << endl;
	meliminacaodegauss();

	return 0;
}