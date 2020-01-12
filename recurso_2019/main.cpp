#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double dyt(double y, double t) { return y / (t - 1); }

double w(double x, double y) { return -1.7 * x * y + 12 * y + 7 * pow(x, 2) - 8 * x; }
double wdx(double x, double y) { return -1.7*y + 14 * x - 8; }
double wdy(double x, double y) { return -1.7*x + 12; }

double f1(double x, double y) { return pow(x, 2) - y - 2; }
double df1x(double x, double y) { return 2 * x; }
double df1y(double x, double y) { return -1; }
double f2(double x, double y) { return -x * pow(y, 2) - 2; }
double df2x(double x, double y) { return -1; }
double df2y(double x, double y) { return 2 * y; }

double f(double x) { return exp(1.5*x); }

void meuler(double t, double y, double h)
{
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\ty: " << y << endl;
		y += h * dyt(y, t);
		t += h;
	}
}

void mrk4(double t, double y, double h)
{
	double o1, o2, o3, o4, deltay;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\ty: " << y;
		o1 = h * dyt(y, t);
		o2 = h * dyt(y + o1 / 2.0, t + h / 2.0);
		o3 = h * dyt(y + o2 / 2.0, t + h / 2.0);
		o4 = h * dyt(y + o3, t + h);
		cout << "\tdy1: " << o1 << "\tdy2: " << o2 << "\tdy3: " << o3 << "\tdy4: " << o4 << endl;
		deltay = o1 / 6.0 + o2 / 3.0 + o3 / 3.0 + o4 / 6.0;
		y += deltay;
		t += h;
		
	}
}

vector<double> gradiente(double x, double y)
{
	vector<double> gradiente;
	gradiente.push_back(wdx(x, y));
	gradiente.push_back(wdy(x, y));
	return gradiente;
}

void mgradiente(double x, double y, double lambda)
{
	vector<double> Xn = { 1000,1000 };
	vector<double> Xn1 = { x,y };
	for (int i = 0; i < 2; i++)
	{
		cout << "iteracao: " << i << "\tw(x,y): " << w(Xn1[0], Xn1[1]) << endl;
		Xn[0] = Xn1[0];
		Xn[1] = Xn1[1];
		Xn1[0] = Xn[0] - lambda * gradiente(Xn[0], Xn[1])[0];
		Xn1[1] = Xn[1] - lambda * gradiente(Xn[0], Xn[1])[1];

		if (w(Xn1[0], Xn1[1]) < w(Xn[0], Xn[0]))
			lambda = 2 * lambda;
		else
		{
			Xn1[0] = Xn[0];
			Xn1[1] = Xn[1];
			lambda = lambda / 2.0;
		}
	}
}

void mnewton_sistemas(double x, double y)
{
	double hn, kn, jacobian;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tx: " << x << "\ty: " << y << endl;
		jacobian = df1x(x, y) * df2y(x, y) - df1y(x, y) * df2x(x, y);
		hn = - (f1(x, y) * df2y(x, y) - df1y(x, y) * f2(x, y)) / jacobian;
		kn = - (df1x(x, y) * f2(x, y) - f1(x, y) * df2x(x, y)) / jacobian;

		x += hn;
		y += kn;
	}
}

double msimpson(double a, double b, double h)
{
	double s = f(a) + f(b);
	double n = (b - a) / h;
	for (int i = 1; i < n; i++)
	{
		if (i % 2 == 0)
			s += 2 * f(a + i * h);
		else
			s += 4 * f(a + i * h);
	}
	s *= h / 3;
	return s;
}

int main()
{
	cout << "\t\tMETODO DE EULER VS METODO RK4 - EX1" << endl << endl;
	cout << "\t\t\tMETODO DE EULER" << endl << endl;
	meuler(2, 2, 0.25);
	cout << "\n\t\t\tMETODO DE RK4" << endl << endl;
	mrk4(2, 2, 0.25);
	cout << "\n\t\tMETODO DO GRANDIENTE - EX3" << endl << endl;
	mgradiente(2.4, 4.3, 0.1);
	cout << "\n\t\tMETODO DE NEWTON SISTEMAS - EX5" << endl << endl;
	mnewton_sistemas(1.5, 0.8);
	cout << "\n\t\tMETODO DE SIMPSON - EX7" << endl << endl;
	double s = msimpson(2.5, 3, 0.125);
	double s1 = msimpson(2.5, 3, 0.125/2.0);
	double s2 = msimpson(2.5, 3, 0.125/4.0);
	cout << "S: " << s << "\tS': " << s1 << "\tS'': " << s2 << endl;
	double qc = (s1 - s) / (s2 - s1);
	cout << "Qc: " << qc << endl;
	double e = (s2 - s1) / 15.0;
	cout << "E: " << e << endl;
	return 0;
}