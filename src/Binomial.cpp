#include"../include/Binomial.h"
using namespace std;

double Em[6][4];

double cal_R2(vector<double> vx, vector<double> vy, double coefficient[])
{
    double a = coefficient[3];
    double b = coefficient[2];
    double c = coefficient[1];
    vector<double> vz(vx.size());
    double SST = 0, SSE = 0;
    double average = accumulate(vy.begin(), vy.end(), 0.0) / vy.size();
    for (int i = 0; i < vx.size(); i++)
    {
        SST += (vy[i] - average) * (vy[i] - average);
        vz[i] = a * vx[i] * vx[i] + b * vx[i] + c;
        SSE += (vz[i] - vy[i]) * (vz[i] - vy[i]);
    }
    return (1 - SSE / SST);
}

double sum(vector<double> Vnum, int n)
{
    double dsum = 0;
    for (int i = 0; i < n; i++)
    {
        dsum += Vnum[i];
    }
    return dsum;
}

double MutilSum(vector<double> Vx, vector<double> Vy, int n)
{
    double dMultiSum = 0;
    for (int i = 0; i < n; i++)
    {
        dMultiSum += Vx[i] * Vy[i];
    }
    return dMultiSum;
}

double RelatePow(vector<double> Vx, int n, int ex)
{
    double ReSum = 0;
    for (int i = 0; i < n; i++)
    {
        ReSum += pow(Vx[i], ex);
    }
    return ReSum;
}

double RelateMutiXY(vector<double> Vx, vector<double> Vy, int n, int ex)
{
    double dReMultiSum = 0;
    for (int i = 0; i < n; i++)
    {
        dReMultiSum += pow(Vx[i], ex) * Vy[i];
    }
    return dReMultiSum;
}

void EMatrix(vector<double> Vx, vector<double> Vy, int n, int ex, double coefficient[])
{
    for (int i = 1; i <= ex; i++)
    {
        for (int j = 1; j <= ex; j++)
        {
            Em[i][j] = RelatePow(Vx, n, i + j - 2);
        }
        Em[i][ex + 1] = RelateMutiXY(Vx, Vy, n, i - 1);
    }
    Em[1][1] = n;
    CalEquation(ex, coefficient);
}

void CalEquation(int exp, double coefficient[])
{
    for (int k = 1; k < exp; k++) 
    {
        for (int i = k + 1; i < exp + 1; i++)
        {
            double p1 = 0;

            if (Em[k][k] != 0)
                p1 = Em[i][k] / Em[k][k];

            for (int j = k; j < exp + 2; j++)
                Em[i][j] = Em[i][j] - Em[k][j] * p1;
        }
    }
    coefficient[exp] = Em[exp][exp + 1] / Em[exp][exp];
    for (int l = exp - 1; l >= 1; l--)   
        coefficient[l] = (Em[l][exp + 1] - F(coefficient, l + 1, exp)) / Em[l][l];
}

double F(double c[], int l, int m)
{
    double sum = 0;
    for (int i = l; i <= m; i++)
        sum += Em[l - 1][i] * c[i];
    return sum;
}