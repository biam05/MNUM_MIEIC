#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

double f(double x){	return pow(x - 1, 2) + pow(x, 4);}

double y(double x) { return sqrt(1 + pow(2.5*exp(2.5*x), 2)); }

double f1(double x) { return exp(x) - x - 5; }
double df1x(double x) { return exp(x) - 1; }
double f2(double x) { return exp(x) - 5; }
double f3(double x) { return log(5 + x); }

double dCdt(double C, double T) { return -exp(-0.5 / (T + 273))*C; }
double dTdt(double C, double T) { return 30 * exp(-0.5 / (T + 273)) * C - 0.5 * (T - 20); }

double w(double x, double y) { return -1.1*x*y + 12 * y + 7 * pow(x, 2) - 8 * x; }
double wdx(double x, double y) { return -1.1 * y + 14 * x - 8; }
double wdy(double x, double y) { return -1.1*x + 12; }

void aurea(double x1, double x2, double precision)
{
	const double B = (sqrt(5) - 1) / 2;
	const double A = B * B;
	double x3 = x1 + A * (x2 - x1);
	double x4 = x1 + B * (x2 - x1);

	while (abs(x1 - x4) > precision && abs(x3 - x2) > precision)
	{
		if (f(x3) < f(x4))
		{
			x2 = x4;
			x4 = x3;
			x3 = x1 + B * (x4 - x1);
		}
		else
		{
			x1 = x3;
			x3 = x4;
			x4 = x3 + B * (x2 - x3);
		}
	}

	cout << "x1: " << x1 << "\tx2: " << x2 << "\tx3: " << x3 << "\tx4: " << x4 << endl;

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
		if(i % 2 == 0)
			s += 2 * y(a + i * h);
		else
			s += 4 * y(a + i * h);
	}
	s *= h / 3;
	return s;
}

void mnewton(double x)
{
	int i = 0;
	double x1 = 100;
	while(abs(x1 - x) > 0.00001)
	{
		x1 = x;
		x = x1 - f1(x1) / df1x(x1);
		i++;
	}
	cout <<	"i: " << i << "\tx: " << x << endl;
}

void mpicardpeano(double x)
{
	int i = 0;
	double x1 = 100;
	while (abs(x1 - x) > 0.00001)
	{
		x1 = x;
		x = f3(x1);
		i++;
	}
	cout << "i: " << i << "\tx: " << x << endl;
}

double meuler(double h, double C, double T)
{
	double t = 0, T1, C1, i = 0, tprev;
	while( t <= 0.5)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\tC: " << C << "\tT: " << T << endl;
		tprev = T;
		t += h;
		T1 = T + dTdt(C, T) * h;
		C1 = C + dCdt(C, T) * h;
		T = T1;
		C = C1;
		i++;
	}
	return tprev;
}

double mrk4(double h, double C, double T)
{
	double o1c, o2c, o3c, o4c, o1t, o2t, o3t, o4t, t = 0;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\tC: " << C << "\tT: " << T << endl;
		t += h;
		o1c = h * dCdt(C, T);
		o1t = h * dTdt(C, T);
		o2c = h * dCdt(C + o1c / 2, T + o1t / 2);
		o2t = h * dTdt(C + o1c / 2, T + o1t / 2);
		o3c = h * dCdt(C + o2c / 2, T + o2t / 2);
		o3t = h * dTdt(C + o2c / 2, T + o2t / 2);
		o4c = h * dCdt(C + o3c, T + o3t);
		o4t = h * dTdt(C + o3c, T + o3t);
		C += (o1c / 6.0 + o2c / 3.0 + o3c / 3.0 + o4c / 6.0);
		T += (o1t / 6.0 + o2t / 3.0 + o3t / 3.0 + o4t / 6.0);
	}
	return T;
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

		if (Xn1[0] < Xn[0] && Xn1[1] < Xn[1])
			lambda = 2 * lambda;
		else
		{
			Xn1[0] = Xn[0];
			lambda = lambda / 2.0;
		}		
	}
}

int main()
{
	cout << "\t\tMETODO DA SECCAO AUREA - EX1" << endl << endl;
	aurea(1.0, 3.0, 0.001);


	cout << "\n\t\tQUADRATURA - EX2" << endl << endl;

	cout << "\t\t\tMETODO DOS TRAPEZIOS" << endl << endl;

	double it1 = mtrapezios(0, 1.0, 0.125);
	double it2 = mtrapezios(0, 1.0, 0.125 / 2.0);
	double it3 = mtrapezios(0, 1.0, 0.125 / 4.0);
	cout << "I: " << it1 << "\tI': " << it2 << "\tI'': " << it3 << endl;
	double qct = (it2 - it1) / (it3 - it2);
	cout << "Qc: " << qct << endl;
	double et = (it3 - it2) / 3;
	cout << "E: " << et << endl;

	cout << "\n\t\t\tMETODO DE SIMPSON" << endl << endl;
	double is1 = msimpson(0, 1.0, 0.125);
	double is2 = msimpson(0, 1.0, 0.125 / 2.0);
	double is3 = msimpson(0, 1.0, 0.125 / 4.0);
	cout << "I: " << is1 << "\tI': " << is2 << "\tI'': " << is3 << endl;
	double qcs = (is2 - is1) / (is3 - is2);
	cout << "Qc: " << qcs << endl;
	double es = (is3 - is2) / 15;
	cout << "E: " << es << endl;


	cout << "\n\t\tPICARD-PEANO VS NEWTON - EX3" << endl << endl;

	cout << "\t\t\tMETODO DE NEWTON" << endl << endl;
	mnewton(2);

	cout << "\n\t\t\tMETODO DE PICARD PEANO" << endl << endl;
	mpicardpeano(2);


	cout << "\n\t\tMETODO DE EULER vs METODO RK4 - EX4" << endl << endl;

	cout << "\t\t\tMETODO DE EULER" << endl << endl;
	double Teuler = meuler(0.25, 2.5, 25);

	cout << "\n\t\t\tMETODO DE RK4" << endl << endl;
	double Trk41 = mrk4(0.25, 2.5, 25);

	cout << "\n\t\t\tQC E ERRO NO METODO DE EULER (C)" << endl << endl;
	double Teuler1 = meuler(0.25, 2.5, 25.0);
	cout << endl;
	double Teuler2 = meuler(0.25 / 2.0, 2.5, 25.0);
	cout << endl;
	double Teuler3 = meuler(0.25 / 4.0, 2.5, 25.0);
	cout << endl << "T': " << Teuler2 << "\tT'': " << Teuler3 << endl;
	double QCTeuler = (Teuler2 - Teuler1) / (Teuler3 - Teuler2);
	cout << "Qc: " << QCTeuler << endl;
	double ETeuler = (Teuler3 - Teuler2);
	cout << "E: " << ETeuler << endl;

	cout << "\n\t\tMETODO DO GRADIENTE - EX5" << endl << endl;
	mgradiente(3, 1, 0.1);

	return 0;
}