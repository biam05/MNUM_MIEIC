#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

/*funções do enunciado*/

double dCdt(double C, double T) {
	return -exp(-0.5 / (T + 273))*C;
}

double dTdt(double C, double T) {
	return 20.00000 * exp(-0.50000 / (T + 273))*C - 0.50000 * (T - 20);
}

double dTdt4(double T, double Ta) { 
	return -0.25 *(T - Ta); 
}

double f1(double x, double y, double z, double w)
{
	return -( 0.5 * y + 3 * z + 0.25 * w - 2.5)/6;
}

double f2(double x, double y, double z, double w)
{
	return -(1.2 * x + 0.25 * z + 0.2 * w - 3.8)/3;
}

double f3(double x, double y, double z, double w)
{
	return -(-1 * x + 0.25 * y  + 2* w - 10)/4;
}

double f4(double x, double y, double z, double w)
{
	return -(2 * x + 4 * y + 1 * z - 7)/8;
}

/*algoritmos*/

double eulersistemasedo(double t, double C, double T, double h, int n)
{
	double deltaC, deltaT, cres;
	for (int i = 0; i <= n; i++)
	{
		cres = C;
		cout << "it: " << i << "\t| t: " << t << "\t| C: " << C << 
			"\t| T: " << T << endl;
		deltaC = h * dCdt(C, T);
		deltaT = h * dTdt(C, T);
		C += deltaC;
		T += deltaT;
		t += h;
	}
	return cres;
}

void rk4sistemasedo(double t, double C, double T, double h)
{
	double c1, t1;
	double dc1, dc2, dc3, dc4, dt1, dt2, dt3, dt4;

	for (int i = 0; i < 3; i++)
	{
		cout << "it: " << i << "\t| t: " << t << "\t| C: " << C <<
			"\t| T: " << T << endl;
		dc1 = h * dCdt(C, T);
		dt1 = h * dTdt(C, T);
		dc2 = h * dCdt(C + dc1/2, T + dt1/2);
		dt2 = h * dTdt(C + dc1/2, T + dc1/2);
		dc3 = h * dCdt(C + dc2/2, T + dt2/2);
		dt3 = h * dTdt(C + dc2/2, T + dc2/2);
		dc4 = h * dCdt(C + dc3, T + dt3);
		dt4 = h * dTdt(C + dc3, T + dt3);
		C += dc1 / 6 + dc2 / 3 + dc3 / 3 + dc4 / 6;
		T += dt1 / 6 + dt2 / 3 + dt3 / 3 + dt4 / 6;
		t += h;
	}
}

double simpson(double a, double b, double h)
{
	//vector<double> f = { 0.18,0.91,0.83,1.23,0.88,1.37,0.80,1.34,0.43 };
	vector<double> f = { 1.02, 1.21, 1.45, 0.89, 0.62, 1.46, 0.74, 0.36, 0.87 };
	double s = f[a];
	double n = (b - a) / h;
	for (int i = 1; i < n; i++)
	{
		if (i % 2 == 0)
			s += 2 * f[(a + i * h)/0.2];
		else
			s += 4 * f[(a + i * h)/0.2];
	}
	s += f[b/0.2];
	
	s = (h / 3) * s;
	return s;
	
}

double euler(double t, double T, double Ta, double h)
{
	double deltaC, deltaT, cres;
	for (int i = 0; i <= 2; i++)
	{
		cres = T;
		cout << "it: " << i << "\t| t: " << t << "\t| T: " << T << endl;
		deltaT = h * dTdt4(T, Ta);
		T += deltaT;
		t += h;
	}
	return cres;
}

void gauss_seidel(double x, double y, double z, double w) 
{
	
	x = f1(x, y, z, w);
	y = f2(x, y, z, w);
	z = f3(x, y, z, w);
	w = f4(x, y, z, w);

	cout << "x: " << x << "\t|y: " << y << "\t|z: " << z << "\t|w: " << w << endl;
}

int main()
{
	cout << "::::: METODO DE EULER (EX1A) :::::" << endl;
	eulersistemasedo(0, 1.00000, 15.00000, 0.25, 2);
	cout << "\n::::: METODO DE RK4 (EX1B) :::::" << endl;
	rk4sistemasedo(0, 1.00000, 15.00000, 0.25);
	cout << "\n::::: METODO DE EULER - QC & E (EX1C) :::::" << endl;
	double C = eulersistemasedo(0, 1.00000, 15.00000, 0.25, 2);
	double C1 = eulersistemasedo(0, 1.00000, 15.00000, 0.25/2, 4);
	double C2 = eulersistemasedo(0, 1.00000, 15.00000, 0.25/4, 8);
	double Qc = (C1 - C) / (C2 - C1);
	double E = abs(C2 - C1);
	cout << "C:" << C << "\t| C1: " << C1 << "\t| C2: " << C2 << endl;
	cout << "Qc: " << Qc << "\t| E:" << E << endl;
	cout << "\n::::: SIMPSON (EX2) :::::" << endl;
	double S = simpson(0, 1.6, 0.8);
	double S1 = simpson(0, 1.6, 0.4);
	double S2 = simpson(0, 1.6, 0.2);
	cout << "S2: " << S2 << endl;
	double Esimpson = (S2 - S1)/15;
	cout << "E: " << Esimpson << endl;
	cout << "\n::::: METODO DE EULER (EX4) :::::" << endl;
	euler(5, 10, 42, 0.4);
	cout << "\n::::: METODO DE GAUSS-SEIDEL (EX5) :::::" << endl;
	gauss_seidel(0, 0, 0, 0);

	return 0;
}