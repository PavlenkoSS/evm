#include "headMat.h"
using namespace std;

// Набор матричного джентльмена//////////////////////////
 
//
//// Заполняет единичную матрицу
//void idMat(double (*mat), int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			mat[i*n+j] = 0;
//		}
//	}
//	for (int i = 0; i < n; i++)
//	{
//		mat[i*(n+1)] = 1;
//	}
//} 
//// Заполняет матрицу по желаемому правилу
//int fulMat(double (*mat), int n, int k, string filename)
//{
//	ifstream fin;
//	double I = 0, J = 0;
//	double N = n;
//	switch (k)
//	{
//	case 0:
//		fin.open(filename);
//		if (!fin.is_open())
//		{
//			cout << "Open was not opened";
//			return -1; // не открылся файл
//		}
//		for (int i = 0; i < n; i++)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				double a;
//				fin >> a;
//				mat[i * n + j] = a;
//			}
//			cout << endl;
//		}
//		cout << endl << endl;
//		fin.close();
//		break;
//	case 1: // corrected
//		for (int i = 0; i < n; i++, I++)
//		{
//			J = 0;
//			for (int j = 0; j < n; j++, J++)
//			{
//				mat[i*n+j] = N - max(I + 1, J + 1) + 1; // шо делать блинб
//			}
//		}
//		break;
//	case 2: // есть разница (функция как в задании)
//		for (int i = 0; i < n;I++, i++)
//		{
//			J = 0;
//			for (int j = 0; j < n; J++, j++)
//			{
//				mat[i*n+j]=max(I + 1, J + 1);
//			}
//		}
//		break;
//	case 3: // без разницы (функция как в задании)
//		for (int i = 0; i < n; I++, i++)
//		{
//			J = 0;
//			for (int j = 0; j < n; J++, j++)
//			{
//				mat[i*n+j]=abs(I - J);
//			}
//
//		}
//		break;
//	case 4: // есть разница (функция изменена)
//		for (int  i = 0; i < n;I++, i++)
//		{
//			J = 0;
//			for (int j = 0; j < n;J++,  j++)
//			{
//				mat[i*n+j]= (1 / (I + J + 1));
//			}
//		}
//		break;
//	case 5:
//		for (int i = 0; i < n; I++, i++)
//		{
//			J = 0;
//			for (int j = 0; j < n; J++, j++)
//			{
//				mat[i*n+j] = J + I;
//			}
//		}
//		break;
//	case 6:
//		for (int i = 0; i < n; I++ ,i++)
//		{
//			for (int j = 0; j < n; J++, j++)
//			{
//				mat[i * n + j] = 0;
//			}
//			mat[i * n + i] = I+1;
//		}
//
//		break;
//	case 7:
//		for (int i = 0; i < n; I++, i++)
//		{
//			J = 0;
//			for (int j = 0; j < i+2; J++, j++)
//			{
//				mat[i * n + j] = abs(I-2*J)+1;
//			}
//			for (int j = i+1; j < n; J++, j++)
//			{
//				mat[i * n + j] = 0;
//			}
//		}
//		break;
//	case 8:
//		for (int i = 0; i < n; I++, i++)
//		{
//			J = 0;
//			for (int j = 0; j < i + 2; J++, j++)
//			{
//		
//				mat[i * n + j] = abs(I * I + J) + 1;
//			}
//		}
//		for (int i = 0; i < n; I++, i++)
//		{
//			for (int j = 0; j < i + 2; J++, j++)
//			{
//				if (i < j)
//				{
//					mat[i * n + j] = 0;
//				}
//			}
//		}
//
//		I=0;
//		break;
//	}
//
//	return 0;
//}
//// Выводит матрицу
//int outMat(double mat[], int n, int m)
//{
//	if (n < m)
//	{
//		return -1;
//	}
//	n = m;
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			cout << mat[i*n+j] << ' ';
//		}
//		cout << endl;
//	}
//	cout << endl << endl;
//}
//// Выводит все возможные заполнения (кроме текста пока что)

//void testFillAndOut(double (*M), int n)
//{
//	for (int K = 0; K < 6; K++)
//	{
//		cout << K << endl;
//		fulMat(M, n, K, "file.txt");
//		cout << endl;
//		outMat(M, n);
//	}
//}

////домножаем mat2 := mat1 * mat2
//void multMat(double(*mat1), double(*mat2), double(*m), int n)
//{
//
//	double a = 0;
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			a = 0;
//			for (int k = 0; k < n; k++)
//			{
//				a = a + mat1[i * n + k] * mat2[k * n + j];
//			}
//			m[i * n + j] = a;
//		}
//	}
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			mat2[i * n + j] = m[i * n + j];
//		}
//	}
//}
//// умножение на матрицу поворота
//
//// //////////////////////////////////////////////////////
//int uptriangleMat(double(*mat), double(*tam), double (*m), double (*t),  int n, int I, int J)
//{
//	double x = mat[I * n + I];
//	double y = mat[J * n + I];
//	//cout << "mat[" << I << "][" << I << "] = " << x << endl;
//	//cout << "mat[" << J << "][" << I << "] = " << y << endl;
//	
//	if ((abs(x) > 1e-16) || (abs(y) > 1e-16))
//	{
//		double Cos = x / sqrt(x * x + y * y);
//		double Sin = -y / sqrt(x * x + y * y);
//		//cout << "Cos = " << Cos << " Sin = " << Sin << endl;
//		for (int j = I; j < n; j++)
//		{
//			double a,b;
//			/*a = Cos * mat[I * n + j] - Sin * mat[J * n + j]; 
//			b = Sin * mat[I * n + j] + Cos * mat[J * n + j];
//			m[I * n + j] = a;
//			m[J * n + j] = b;
//			a = Cos * tam[I * n + j] - Sin * tam[J * n + j];
//			b = Sin * tam[I * n + j] + Cos * tam[J * n + j];
//			t[I * n + j] = a;
//			t[J * n + j] = b;*/
//			m[j] = Cos * mat[I * n + j] - Sin * mat[J * n + j];
//			m[n + j] = Sin * mat[I * n + j] + Cos * mat[J * n + j];
//		}
//		for (int j = I; j < n; j++)
//		{
//			double a, b;
//			/*a = Cos * mat[I * n + j] - Sin * mat[J * n + j];
//			b = Sin * mat[I * n + j] + Cos * mat[J * n + j];
//			m[I * n + j] = a;
//			m[J * n + j] = b;
//			a = Cos * tam[I * n + j] - Sin * tam[J * n + j];
//			b = Sin * tam[I * n + j] + Cos * tam[J * n + j];
//			t[I * n + j] = a;
//			t[J * n + j] = b;*/
//			mat[I * n + j] = m[j];
//			mat[J * n + j] = m[n + j];
//		}
//		for (int j = 0; j < n; j++)
//		{
//			t[j] = Cos * tam[I * n + j] - Sin * tam[J * n + j];
//			t[n + j] = Sin * tam[I * n + j] + Cos * tam[J * n + j];
//		}
//		for (int j = 0; j < n; j++)
//		{
//			tam[I * n + j] = t[j];
//			tam[J * n + j] = t[n + j];
//		}
//
//		return 0;
//	}
//
//	return 1;
//}
//
//
//// не работает
//
////int downtriangleMat(double(*mat), double(*tam), int n, int I, int J)
////{
////	double x = mat[I * n + I];
////	double y = mat[J * n + I];
////	//cout << "mat[" << I << "][" << I << "] = " << x << endl;
////	//cout << "mat[" << J << "][" << I << "] = " << y << endl;
////	if ((x != 0) || (y != 0))
////	{
////		double Cos = x / sqrt(x * x + y * y);
////		double Sin = -y / sqrt(x * x + y * y);
////		//cout << "Cos = " << Cos << " Sin = " << Sin << endl;
////		for (int j = I; j < n; j++)
////		{
////			double a, b;
////			a = Cos * mat[j * n + I] + Sin * mat[j * n + J];
////			b = -Sin * mat[j * n + I] + Cos * mat[j * n + J];
////			mat[j * n + I] = a;
////			mat[j * n + J] = b;
////			a = Cos * tam[j * n + I] + Sin * tam[j * n + J];
////			b = -Sin * tam[j * n + I] + Cos * tam[j * n + J];
////			tam[j * n + I] = a;
////			tam[j * n + J] = b;
////			//a = Cos * mat[j * n + I] + Sin * mat[j * n + J];
////			//b = -Sin * mat[j * n + I] + Cos * mat[j * n + J];
////			//mat[j * n + I] = a;
////			//mat[j * n + J] = b;
////			//a = Cos * tam[j * n + I] + Sin * tam[j * n + J];
////			//b = -Sin * tam[j * n + I] + Cos * tam[j * n + J];
////			//tam[j * n + I] = a;
////			//tam[j * n + J] = b;
////		}
////		return 0;
////	}
////	cout << " err";
////	return 1;
////}
//
//
//int rotMat(double (*mat),double(*T), int n, int I, int J) //Можно убрать шелуху здесь 
//{
//	double x = mat[J * n + J];
//	double y = mat[I * n + J];
//
//	idMat(T, n);
//	if ((x != 0) || (y != 0))
//	{
//		double Cos = x / sqrt(x * x + y * y);
//		double Sin = -y / sqrt(x * x + y * y);
//		T[J * n + J] = Cos;
//		T[J * n + I] = -Sin;
//		T[I * n + J] = Sin;
//		T[I * n + I] = Cos;
//		return 0;
//	}
//	return 1;
//}
//int MatInverse(double mat[], double (*tam), double (*m), double (*t), int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		int bad = 1;
//		for (int j = i+1; j < n; j++)
//		{
//			bad = bad + uptriangleMat(mat, tam, m, t, n, i, j);
//		}
//		if (bad == n)
//		{
//			cout << "Your matrix is uninvertable" << endl;
//			return 1;
//		}
//	}
//
//
//	for (int i = n - 1; i > -1; i--)
//	{
//		if (abs(mat[i * n + i]) < 1e-15)
//		{
//			cout << "UNINVERTABLE " << endl;
//			return 1;
//			break;
//		}
//		for (int k = 0; k < n; k++)
//		{
//			tam[i * n + k] = tam[i * n + k] / mat[i * n + i];
//		}
//		
//		//for (int k = 0; k < i; k++)
//		//{
//		//	mat[i * n + k] = 0;
//		//}
//		for (int k = i+1; k < n; k++)
//		{
//			mat[i * n + k] = mat[i * n + k] / mat[i * n + i];
//		}
//		mat[i * n + i] = 1;
//	}
//	for (int J = n - 1; J > 0; J--)
//	{
//		for (int I = J-1; I >-1; I--)
//		{
//			for (int j = 0; j < n; j++)
//			{
//				tam[I * n + j] = tam[I * n + j] - mat[I * n + J] * tam[J * n + j];
//			}
//		}
//	}
//	return 0;
//}
////норма невязки
//double normMat(double(*mat), double(*tam), int n)
//{
//	double a = 0, A=0;
//	for (int j = 0; j < n; j++)
//	{
//		a = a + mat[j] - tam[j];
//	}
//	A = abs(a);
//	for (int i = 1; i < n; i++)
//	{
//		a = 0;
//		for (int j = 0; j < n;j++)
//		{
//			a = a + mat[i * n + j] - tam[i * n + j];
//		}
//		a = abs(a);
//		if (a > A)
//		{
//			A = a;
//		}
//	}
//	return A;
//}


// переделать массив если ноль, надо делать физически 
// в методе гаусса выбирать главный элемент -самый большой по модулю 

int main(int nargs, char** args)
{

	if (nargs < 2)
	{
		cout << "Not enought args: launch as following..." << endl;
		return -1;
	}
	int n = 0;
	int m = 0;
	int k = 1;
	char* filename = new char[128];

	if (sscanf(args[1], "%d", &n) != 1 || sscanf(args[2], "%d", &m) != 1 || sscanf(args[3], "%d", &k) != 1)
	{
		cout << "bad params" << endl;
		return -1;
	}
	if(n==0)
	{
		cout << "bad param n" << endl;
		return -1;
	}
	 if(args[4] == nullptr || sscanf(args[4], "%s", filename) != 1)
            {
                cout << "No filename" << endl;
                return -1;
            }
	double* M;
	double* M_copy;
	double* E;
	double* mem;
	double* mem2;
	mem = new double[n * n];
	M = new double[n * n];
	mem2 = new double[n + n];
	M_copy = new double[n * n];
	E = new double[n * n];	

	idMat(E, n);

	
	if (fulMat(M, n, k, string(filename)) == -1)
	{
		return -1;
	}

	eqMat(M_copy, M, n, n);
	cout << "M = " << endl;
	outMat1(M, n,m);

	int start_time = clock();
	if (MatInverse(M, E, mem2, n) == -2)
	{
		delete[]M_copy;
		delete[]E;
		delete[]M;
		delete[]mem;
		delete[]mem2;
		return -2;
	}
	int end_time = clock();
	cout <<"Time of inversing = " << (-start_time + end_time)/1e6 << endl;

	cout << "E = " << endl;
	if (outMat1(E, n, m) == -2)
	{
		delete[]M_copy;
		delete[]E;
		delete[]M;
		delete[]mem;
		delete[]mem2;
		return -2;
	}
	multMat(E, M_copy, mem, n);
	cout << "ANS = ";
	idMat(E, n);
	cout << normMat(M_copy, E , n)<< endl;

	delete[]M;
	delete[]M_copy;
	delete[]E;
	delete[]mem;
	delete[]mem2;
	return 0;
}
