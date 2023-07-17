#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
using std::vector;

double cal_R2(vector<double> vx, vector<double> vy, double coefficient[]);

double sum(vector<double> Vnum, int n);

double MutilSum(vector<double> Vx, vector<double> Vy, int n);

double RelatePow(vector<double> Vx, int n, int ex);

double RelateMutiXY(vector<double> Vx, vector<double> Vy, int n, int ex);

void EMatrix(vector<double> Vx, vector<double> Vy, int n, int ex, double coefficient[]);

void CalEquation(int exp, double coefficient[]);

double F(double c[], int l, int m);