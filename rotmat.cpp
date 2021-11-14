#include "headMat.h"

using namespace std;



int uptriangleMat(double(*mat), double(*tam), double(*m), double(*t), int n, int I, int J, double eps)
{
	double x = mat[I * n + I];
	double y = mat[J * n + I];
	if ((abs(x) > 1e-16*eps) || (abs(y) > 1e-16*eps))
	{
		double Cos = x / sqrt(x * x + y * y);
		double Sin = -y / sqrt(x * x + y * y);
		for (int j = I; j < n; j++)
		{

			m[j] = Cos * mat[I * n + j] - Sin * mat[J * n + j];
			m[n + j] = Sin * mat[I * n + j] + Cos * mat[J * n + j];
		}
		for (int j = I; j < n; j++)
		{
			mat[I * n + j] = m[j];
			mat[J * n + j] = m[n + j];
		}
		for (int j = 0; j < n; j++)
		{
			t[j] = Cos * tam[I * n + j] - Sin * tam[J * n + j];
			t[n + j] = Sin * tam[I * n + j] + Cos * tam[J * n + j];
		}
		for (int j = 0; j < n; j++)
		{
			tam[I * n + j] = t[j];
			tam[J * n + j] = t[n + j];
		}

		return 0;
	}

	return 1;
}

//����� �������
double normMat(double(*mat), double(*tam), int n)
{
	double a = 0, A = 0;
	for (int j = 0; j < n; j++)
	{
		a +=  mat[j] - tam[j];
	}
	A = abs(a);
	for (int i = 1; i < n; i++)
	{
		a = 0;
		for (int j = 0; j < n; j++)
		{
			a +=  mat[i * n + j] - tam[i * n + j];
		}
		a = abs(a);
		if (a > A)
		{
			A = a;
		}
	}
	return A;
}
/*
int rotMat(double(*mat), double(*T), int n, int I, int J) //����� ������ ������ ����� 
{
	double x = mat[J * n + J];
	double y = mat[I * n + J];

	idMat(T, n);
	if ((abs(x)  < 1e-16) || (abs(y)  < 1e-16))
	{
		double Cos = x / sqrt(x * x + y * y);
		double Sin = -y / sqrt(x * x + y * y);
		T[J * n + J] = Cos;
		T[J * n + I] = -Sin;
		T[I * n + J] = Sin;
		T[I * n + I] = Cos;
		return 0;
	}
	return 1;
}
*/
int MatInverse(double mat[], double(*tam), double(*m), int n)
{
	double eps = normByMaxMat(mat, n);

	
	for (int i = 0; i < n; i++)
	{
		int bad = 1;
		for (int j = i + 1; j < n; j++)
		{
			//bad += uptriangleMat(mat, tam, m, t, n, i, j, eps);
			int I = i;
			int J = j;
			double x = mat[I * n + I];
			double y = mat[J * n + I];
			if ((abs(x) > 1e-25*eps) || (abs(y) > 1e-25*eps))
			{
				double Cos = x / sqrt(x * x + y * y);
				double Sin = -y / sqrt(x * x + y * y);
				for (int j1 = I; j1 < n; j1++)
				{

					m[j1] = (Cos) * mat[I * n + j1] - (Sin) * mat[J * n + j1];
					m[n + j1] = (Sin) * mat[I * n + j1] + (Cos) * mat[J * n + j1];
				}
				for (int j1 = I; j1 < n; j1++)
				{
					mat[I * n + j1] = m[j1];
					mat[J * n + j1] = m[n + j1];
				}
				for (int j1 = 0; j1 < n; j1++)
				{
					m[j1] = (Cos) * tam[I * n + j1] - (Sin) * tam[J * n + j1];
					m[n + j1] = (Sin) * tam[I * n + j1] + (Cos) * tam[J * n + j1];
				}
				for (int j1 = 0; j1 < n; j1++)
				{
					tam[I * n + j1] = m[j1];
					tam[J * n + j1] = m[n + j1];
				}

			}else
			{
				cout << i << ' ' << j << endl;
				bad++;
			}
		}
		if (bad == n)
		{
			cout << "Your matrix is uninvertable" << endl;
			return 1;
		}
	}

	for (int i = n - 1; i > -1; i--)
	{
		if (abs(mat[i * n + i]) < eps*1e-25)
		{
			cout << "UNINVERTABLE " << endl;
			return -2;
		}
		for (int k = 0; k < n; k++)
		{
			tam[i * n + k] /= mat[i * n + i];
		}
		for (int k = i + 1; k < n; k++)
		{
			mat[i * n + k] /= mat[i * n + i];
		}
		mat[i * n + i] = 1;
	}
	for (int J = n - 1; J > 0; J--)
	{
		for (int I = J - 1; I > -1; I--)
		{
			for (int j = 0; j < n; j++)
			{
				tam[I * n + j] = tam[I * n + j] - mat[I * n + J] * tam[J * n + j];
			}
		}
	}
	return 0;
}

double normByMaxMat(double(*mat), int n)
{
	double max = 0;
	double a;
	for (int i=0; i < n; i++)
	{
		a = 0;
		for (int j = 0; j < n; j++)
		{
			a = a + abs(mat[i * n + j]);

		}
		if(a>max)
		{
			max = a; 
		}
	}
	return max;
}

