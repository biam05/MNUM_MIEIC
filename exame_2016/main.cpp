#include <cmath>
#include <iostream>

using namespace std;

double f1(double x, double y){ return pow(x, 2) - y - 1.2;}
double f2(double x, double y){ return -x + pow(y, 2) - 1.0;}
double f1x(double x, double y) { return 2 * x; }
double f1y(double x, double y) { return -1; }
double f2x(double x, double y) { return -1; }
double f2y(double x, double y) { return 2*y; }

double f(double x) { return 1 * pow(x, 7) + 0.5*x - 0.5; }

double fy(double t, double y, double z) {return z;}
double fz(double t, double y, double z) {return 0.5 + pow(t,2) + t*fy(t,y,z);}

double f6(double x) { return sqrt(1 + pow(1.5*exp(1.5*x), 2));}

void mnewton(double xn, double yn)
{
	for (int i = 0; i < 3; i++)
	{
		cout << "xn: " << xn << "\t yn: " << yn << endl;
		double jacobian = f1x(xn, yn)* f2y(xn, yn) - f1y(xn, yn) * f2x(xn, yn);
		double hn = (f1(xn, yn)*f2y(xn, yn) - f1y(xn, yn)*f2(xn, yn)) / jacobian;
		double kn = (f1x(xn, yn)*f2(xn, yn) - f1(xn, yn)*f2x(xn, yn)) / jacobian;
		xn -= hn;
		yn -= kn;		
	}
}

void mcorda(double a, double b)
{
	for (int i = 0; i < 4; i++)
	{
		double rr = (a*f(b) - b * f(a)) / (f(b) - f(a));

		cout << "iteracao: " << i + 1 << "\ta (xe): " << a << "\tb (xd) " << b << "\tx (xn): " << rr << "\tf(x): " << f(rr) << endl;
		if (f(b) - f(a) == 0)
			break;
		
		if (f(a) * f(rr) < 0)
			b = rr;
		else
			a = rr;
		
	}
}

void meuler_2a(double t, double y, double z, double h)
{
	double t1, y1, z1;
	for (int i = 0; i < 3 ; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\ty: " << y << endl;
		t1 = t + h;
		y1 = y + z * h;
		z1 = fz(t, y, z)*h + z;
		t = t1;
		y = y1;
		z = z1;
	}
}

void mrk4_2a(double t, double y, double z, double h)
{
	double y1, y2, y3, y4, z1, z2, z3, z4, zn, yn, zn1, yn1;
	for (int i = 0; i < 3; i++)
	{
		cout << "iteracao: " << i << "\tt: " << t << "\ty: " << y << endl;

		z1 = h * fz(t, y, z);
		y1 = h * fy(t, y, z);
		z2 = h * fz(t + h / 2, y + y1 / 2, z + z1 / 2);
		y2 = h * fy(t + h / 2, y + y1 / 2, z + z1 / 2);
		z3 = h * fz(t + h / 2, y + y2 / 2, z + z2 / 2);
		y3 = h * fy(t + h / 2, y + y2 / 2, z + z2 / 2);
		z4 = h * fz(t + h, y + y3, z + z3);
		y4 = h * fy(t + h, y + y3, z + z3);

		zn = z1 / 6.0 + z2 / 3.0 + z3 / 3.0 + z4 / 6.0;
		yn = y1 / 6.0 + y2 / 3.0 + y3 / 3.0 + y4 / 6.0;

		yn1 = y + yn;
		zn1 = z + zn;
		t += h;
		y = yn1;
		z = zn1;
	}
}

double mtrapezios(double a, double b, double h)
{
	double n = (b - a) / h; 
	double s = f6(a) + f6(b);
	for (int i = 1; i < n; i++)
	{
		s += 2 * f6(a + i * h);
	}
	s *= h / 2;
	return s;
}

double msimpson(double a, double b, double h)
{
	double n = (b - a) / h;
	double s = f6(a) + f6(b);
	for (int i = 1; i < n; i++)
	{
		if(i % 2 == 0)
			s += 2 * f6(a + i * h);
		else
			s += 4 * f6(a + i * h);
	}
	s *= h / 3;
	return s;
}

int main()
{
	cout << "\t\tMETODO DE NEWTON (SISTEMAS) - EX 2" << endl << endl;
	mnewton(1, 1);
	cout << "\n\t\tMETODO DA CORDA - EX4" << endl << endl;
	mcorda(0, 0.8);
	cout << "\n\t\tMETODOS DE INTEGRACAO - EX5" << endl << endl;
	cout << "\t\t\tMETODO DE EULER" << endl << endl;
	meuler_2a(0, 0, 1, 0.25);
	cout << "\n\t\t\tMETODO DE RK4" << endl << endl;
	mrk4_2a(0, 0, 1, 0.25);
	cout << "\n\t\tQUADRATURA - EX6" << endl << endl;
	cout << "\t\t\tMETODO DOS TRAPEZIOS" << endl << endl;
	double it1 = mtrapezios(0, 2.0, 0.5);
	double it2 = mtrapezios(0, 2.0, 0.5 / 2.0);
	double it3 = mtrapezios(0, 2.0, 0.5 / 4.0);
	cout << "I: " << it1 << "\tI': " << it2 << "\tI'': " << it3 << endl;
	double qct = (it2 - it1) / (it3 - it2);
	cout << "Qc: " << qct << endl;
	double et = (it3 - it2) / 3;
	cout << "E: " << et << endl;
	cout << "\n\t\t\tMETODO DE SIMPSON" << endl << endl;
	double is1 = msimpson(0, 2.0, 0.5);
	double is2 = msimpson(0, 2.0, 0.5 / 2.0);
	double is3 = msimpson(0, 2.0, 0.5 / 4.0);
	cout << "I: " << is1 << "\tI': " << is2 << "\tI'': " << is3 << endl;
	double qcs = (is2 - is1) / (is3 - is2);
	cout << "Qc: " << qcs << endl;
	double es = (is3 - is2) / 15;
	cout << "E: " << es << endl;
	return 0;
}
