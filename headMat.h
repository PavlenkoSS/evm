#ifndef headMat
#define headMat
#include <iostream>
#include <vector>
#include <fstream>
#include <typeinfo>
#include <ctime>
#include<math.h>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>

using namespace std;

double normByMaxMat(double(*mat), int n);
int uptriangleMat(double(*mat), double(*tam), double(*m), double(*t), int n, int I, int J, double eps);
double normMat(double(*mat), double(*tam), int n);
int rotMat(double(*mat), double(*T), int n, int I, int J);
int MatInverse(double mat[], double(*tam), double(*m), int n);
double maxMat(double(*mat), int n);
void normalizeMat(double(*mat), int n, double nor);
void idMat(double(*mat), int n);
int fulMat(double(*mat), int n, int k, string filename);
void multMat(double(*mat1), double(*mat2), double(*m), int n);
int outMat1(double mat[], int n, int m);
int eqMat(double(*mat1), double(*mat2), int n, int m);
double smartNormMat(double (*mat1), double (*mat2), int n);

#endif
