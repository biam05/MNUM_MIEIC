#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double dy(double t, double y, double z) { return z; }
double dz(double t, double y, double z) { return 0.5 + pow(t, 2) + t * dy(t, y, z); }

double f(double x, double y)
{
	return 3 * pow(x, 2) - x * y + 11 * y + pow(y, 2) - 8 * x;
}

double dfx(double x, double y) { return -y + 6*x - 8; }
double dfy(double x, double y) { return 2*y - x + 11; }

double g(double x) { return 1.5*exp(1.5*x); }

double w(double x) { return (x - 3.7) + pow(cos(x + 1.2), 3); }
double dwx(double x) {return 1 - 3 * pow(cos(x + 1.2), 2) * sin(x + 1.2);}

void meuler(double t, double y, double z, double h)
{
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\ty: " << y << endl;
		y += dy(t, y, z) * h;
		z += dz(t, y, z)* h;
		t += h;

	}
}

void mrk4(double t, double y, double z, double h)
{
	double deltay1, deltay2, deltay3, deltay4, deltaz1, deltaz2, deltaz3, deltaz4, deltay, deltaz;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\ty: " << y << endl;
		deltay1 = h * dy(t + h, y, z);
		deltaz1 = h * dz(t + h, y, z);
		deltay2 = h * dy(t + h/2, y + deltay1 / 2, z + deltaz1 / 2);
		deltaz2 = h * dz(t + h/2, y + deltay1 / 2, z + deltaz1 / 2);
		deltay3 = h * dy(t + h/2, y + deltay2 / 2, z + deltaz2 / 2);
		deltaz3 = h * dz(t + h/2, y + deltay2 / 2, z + deltaz2 / 2);
		deltay4 = h * dy(t + h, y + deltay3, z + deltaz3);
		deltaz4 = h * dz(t + h, y + deltay3, z + deltaz3);

		deltay = deltay1 / 6.0 + deltay2 / 3.0 + deltay3 / 3.0 + deltay4 / 6.0;
		deltaz = deltaz1 / 6.0 + deltaz2 / 3.0 + deltaz3 / 3.0 + deltaz4 / 6.0;

		y += deltay;
		z += deltaz;
		t += h;
	}
}

vector<double> gradiente(double x, double y)
{
	vector<double> gradiente;
	gradiente.push_back(dfx(x, y));
	gradiente.push_back(dfy(x, y));
	return gradiente;
}

void mgradiente(double h, double x, double y)
{
	vector<double>Xn;
	Xn.push_back(1000); Xn.push_back(1000);
	vector<double> Xn1;
	Xn1.push_back(x); Xn1.push_back(y);
	for(int i = 0; i < 2; i++)
	{

		cout << "iteracao: " << i << "\tx: " << Xn1[0] << "\ty: " << Xn1[1] << "\tf(x,y): " << f(Xn1[0], Xn1[1]) << 
			"\tgradientex: " << dfx(Xn1[0],Xn1[1]) << "\tgradientey: " << dfy(Xn1[0],Xn1[1]) << endl;

		Xn[0] = Xn1[0];
		Xn[1] = Xn1[1];
		Xn1[0] = Xn[0] - h * gradiente(Xn[0], Xn[1])[0];
		Xn1[1] = Xn[1] - h * gradiente(Xn[0], Xn[1])[1];

		if (f(Xn1[0], Xn1[1]) < f(Xn[0], Xn[1]))
		{
			continue;
		}
		else if (f(Xn1[0], Xn1[1]) > f(Xn[0], Xn[1]))
		{
			h = h / 2;
			Xn1[0] = Xn[0] - h * gradiente(x, y)[0];
			Xn1[1] = Xn[0] - h * gradiente(x, y)[1];
		}
	}
}

double msimpson(double a, double b, double h)
{
	double s = g(a) + g(b);
	double n = (b - a) / h;
	for (int i = 1; i < n; i++)
	{
		if (i % 2 == 0)
			s += 2 * g(a + i * h);
		else
			s += 4 * g(a + i * h);
	}
	s *= h / 3;

	return s;
}

void mnewton(double x)
{
	x -= w(x) / dwx(x);
	cout << "x: " << x << endl;
}

int main()
{
	cout << "\t\tMETODOS DE INTEGRACAO - EX1" << endl << endl;

	cout << "\t\t\tMETODO DE EULER" << endl << endl;
	meuler(0, 0, 1, 0.25);

	cout << "\n\t\t\tMETODO DE RK4" << endl << endl;
	mrk4(0, 0, 1, 0.25);


	cout << "\n\t\tMETODO DO GRADIENTE - EX3" << endl << endl;
	mgradiente(0.5, 2, 2);


	cout << "\n\t\tMETODO DE SIMPSON - EX4" << endl << endl;

	double s1 = msimpson(1, 1.5, 0.125);
	double s2 = msimpson(1, 1.5, 0.125 / 2.0);
	double s3 = msimpson(1, 1.5, 0.125 / 4.0);
	cout << "S: " << s1 << "\tS': " << s2 << "\tS'': " << s3 << endl;
	double qc = (s2 - s1) / (s3 - s2);
	cout << "QC: " << qc << endl;
	double e = (s3 - s2) / 15.0;
	cout << "E: " << e << endl;


	cout << "\n\t\tMETODO DE NEWTON - EX5" << endl << endl;
	mnewton(3.8);

	return 0;
}