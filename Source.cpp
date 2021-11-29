#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

using namespace std;

const int m = 2; const int N = 11;
double** Initialization(int n);
void Input_X(double* A);
void Input_Y(double* A);
void Output(double** A, int n);
void Copy(double** A, double** copyA, int n);
void Zero(double** A, int n);
void Answer(double* A, int n);
double* Gauss(double** Array, int n, int m);
double Dispers(double* Dx, double* Dy, double* An, double n);

void main()
{
	double* Dx = new double[N];
	double* Dy = new double[N];
	Input_X(Dx);
	Input_Y(Dy);

	for (int j = 0; j < N; j++)
		cout << Dx[j] << " " << Dy[j] << endl;
	cout << endl;

	double* Powx = new double[2 * m];
	for (int i = 0; i < 2 * m; i++)
		Powx[i] = 0;

	for (int k = 1; k <= 2 * m; k++)
		for (int i = 0; i < N; i++)
			Powx[k - 1] += pow(Dx[i], k);

	for (int i = 0; i < 2 * m; i++)
		cout << Powx[i] << " ";
	cout << endl << endl;

	double** sumx = Initialization(m + 1);
	Zero(sumx, m + 1);

	for (int i = 1; i <= m + 1; i++) 
	{
		for (int j = 1; j <= m + 1; j++)
			if ((i + j) >= 3)
				sumx[i - 1][j - 1] = Powx[i + j - 3];
	}
	sumx[0][0] = N;

	Output(sumx, m + 1);
	for (int k = 0; k <= m; k++)
	{
		for (int i = 0; i < N; i++)
			sumx[k][m + 1] += (Dy[i] * pow(Dx[i], k));
	}

	Output(sumx, m + 1);
	double* An = new double[m + 1];
	double** CopyA = Initialization(m + 1);
	Copy(sumx, CopyA, m + 1);

	An = Gauss(CopyA, m + 1, m + 2);
	Answer(An, m + 1);

	cout << "SKO=" << sqrt(Dispers(Dx, Dy, An, m + 1)) << endl;
	cout << "y=" << An[2] << "x^2+" << An[1] << "x" << An[0] << endl;
}

double** Initialization(int n)
{
	double** A = new double* [n];
	for (int i = 0; i < n; i++)
		A[i] = new double[n + 1];
	return A;
}
void Input_X(double* A)
{
	A[0] = 0; A[1] = 1; A[2] = 2; A[3] = 3; A[4] = 4; A[5] = 5;
	A[6] = 6; A[7] = 7; A[8] = 8; A[9] = 9; A[10] = 10;
}
void Input_Y(double* A)
{
	A[0] = 3; A[1] = 87; A[2] = 156; A[3] = 210; A[4] = 238; A[5] = 252;
	A[6] = 239; A[7] = 211; A[8] = 158; A[9] = 90; A[10] = -5;
}
void Output(double** A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < (n + 1); j++)
			cout << setw(15) << A[i][j];
		cout << endl;
	}
	cout << endl;
}
void Copy(double** A, double** copyA, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n + 1; j++)
			copyA[i][j] = A[i][j];
	}
}
void Zero(double** A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < (n + 1); j++)
		{
			A[i][j] = 0;
		}
	}
}
void Answer(double* A, int n)
{
	for (int i = 0; i < n; i++)
		cout << "a[" << i << "]=" << A[i] << endl;
}
double* Gauss(double** Array, int n, int m)
{
	double elem;
	for (int j = 0; j < n; j++)
	{
		double max = 0;
		int str = 0;
		for (int t = j; t < n; t++)
		{
			if (abs(Array[t][j]) > max)
			{
				max = abs(Array[t][j]); str = t;
			}
		}
		if (max > abs(Array[j][j]))
		{
			double* ptr = Array[j];
			Array[j] = Array[str];
			Array[str] = ptr;
		}
		elem = Array[j][j];
		for (int c = j; c < m; c++)
		{
			Array[j][c] /= elem;
		}

		for (int i2 = j + 1; i2 < n; i2++)
		{
			elem = Array[i2][j];
			for (int k = j; k < m; k++)
				Array[i2][k] -= elem * Array[j][k];
		}
	}
	double* xx = new double[m];
	xx[n - 1] = Array[n - 1][n];
	for (int i = n - 2; i >= 0; i--)
	{
		xx[i] = Array[i][n];
		for (int j = i + 1; j < n; j++)
			xx[i] -= Array[i][j] * xx[j];
	}

	for (int i = 0; i < n; i++)
		cout << xx[i] << " ";
	cout << endl;

	return xx;
}
double Dispers(double* Dx, double* Dy, double* An, double n)
{
	double S = 0, temp;
	for (int i = 0; i < N; i++)
	{
		temp = Dy[i];
		for (int j = 0; j < n; j++)
			temp -= An[j] * pow(Dx[i], j);
		S += temp * temp;
	}
	S /= (N - m - 1);
	return S;
}