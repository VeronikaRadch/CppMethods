#include <iostream> 
#include <math.h> 
#include <iomanip> 
using namespace std;
const double eps1 = 0.0001;
const double eps2 = 0.00001;

double function(double x)
{
	return sqrt(x + x * x * x);
}

double function(double x, double y)
{
	return (x * x + 2 * y);
}

double MethodOfTrap(double a, double b)
{
	double sum = 0; double sum_prev = 0;
	double h = (b - a); 
	int kol_vo = 0; 

	do {
		sum_prev = sum;                               
		sum = 0;                                     
		for (int i = 1; i <= ((b - a) / h) - 1; i++)      
		{
			sum += 2 * function(a + h * i);
		}           
		sum = (function(a) + sum + function(b)) * h / 2;
		h /= 2;
		kol_vo++;
	} while (abs(sum - sum_prev) > 3 * eps1);
	cout << "Trampezija: " << endl << kol_vo << " iter" << endl;
	return sum;
}

double MethodOfSimpson(double a, double b)
{
	double sum = 0; double sum_prev = 0;
	double h = (b - a) / 2;                            
	int kol_vo = 0;

	do {
		sum_prev = sum;
		sum = 0;
		for (int i = 1; i <= ((b - a) / h); i += 2)
		{
			sum += 4 * function(a + h * i);
		}
		for (int i = 2; i <= ((b - a) / h) - 1; i += 2)
		{
			sum += 2 * function(a + h * i);
		}
		sum = (function(a) + sum + function(b)) * h / 3;
		h /= 2;                                    
		kol_vo++;
	} while (abs(sum - sum_prev) > 15 * eps2);
	cout << "Simpsson: " << endl << kol_vo << " iter" << endl;
	return sum;
}
double MethodOfCubeSimpson(double a, double b, double c, double d)
{
	int m = 10; int n = 2 * m;
	double hx = (b - a) / (2 * n);
	double hy = (d - c) / (2 * m);
	double sum = 0;

	double xi = a;
	double yi = c;

	double* Xi = new double[2 * n + 1];
	Xi[0] = xi;

	for (int i = 1; i <= 2 * n; i++)
		Xi[i] = Xi[i - 1] + hx;;

	double* Yi = new double[2 * m + 1];
	Yi[0] = yi;

	for (int j = 1; j <= 2 * m; j++)
		Yi[j] = Yi[j - 1] + hy;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			sum += function(Xi[2 * i], Yi[2 * j]);
			sum += 4 * function(Xi[2 * i + 1], Yi[2 * j]);
			sum += function(Xi[2 * i + 2], Yi[2 * j]);
			sum += 4 * function(Xi[2 * i], Yi[2 * j + 1]);
			sum += 16 * function(Xi[2 * i + 1], Yi[2 * j + 1]);
			sum += 4 * function(Xi[2 * i + 2], Yi[2 * j + 1]);
			sum += function(Xi[2 * i], Yi[2 * j + 2]);
			sum += 4 * function(Xi[2 * i + 1], Yi[2 * j + 2]);
			sum += function(Xi[2 * i + 2], Yi[2 * j + 2]);
		}
	}
	sum *= (hx * hy / 9);
	return sum;
}

void main()
{
	double a = 0.6;
	double b = 1.724;

	double trapp;
	double simp;
	double cube_simp;
	trapp = MethodOfTrap(a, b);
	cout << trapp << endl << endl;
	simp = MethodOfSimpson(a, b);
	cout << simp << endl << endl;

	double a1 = 0.0;
	double b1 = 2.0;
	double c1 = 0.5;
	double d1 = 1.5;
	cube_simp = MethodOfCubeSimpson(a1, b1, c1, d1);
	cout << "Cube Simpsson: " << endl << cube_simp << endl;
}