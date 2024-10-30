#pragma once
#ifndef _TOOLS_
#define _TOOLS_
#include<stdlib.h>
#include<stdio.h>
#include<string>
#include<algorithm>
#include<string>
#include<vector>
#include<string.h>
#include<map>
#include<math.h>
#include<read_write.h>
#include<stack>
#include<set>
using namespace std;
string GetInfoINCAR(char key[]);

int IsNum(string symbol);

vector<double> AddVec(vector<double> v1, vector<double> v2);

vector<double> operator-(vector<double> v1, double value);

vector<double> operator+(vector<double> v1, vector<double> v2);

vector<double> operator*(vector<double> v1, vector<double> v2);

vector<double> Intergral(vector<double> v1, vector<double> v2);

double SumVec(vector<double> v);

double volume(double a[3][3]);

void RecMat(double POSI[3][3], double ret[3][3]);

void copy_vec3d(double vec1[3][3], double vec2[3][3]);

void multi_mat(double mat1[3][3], double mat2[3][3], double multi_mat[3][3]);

bool Isequal(double a, double b, double SYMPREC);

vector<int> ExtractNumbersFromString(char a[], int len);

double maxvalue(vector<double> v);

double minvalue(vector<double> v);

double maxvalue(vector<vector<double> > v);

double minvalue(vector<vector<double> > v);

string RemovePostfixNum(string elem);

bool mysort(pair<int, int> a, pair<int, int> b);

class closest {
	int globaldelta;
	vector<int> result;
	vector<int> path;
public:
	vector<int> getclosest(vector<pair<int, int> > v, int target);
	void backtrack(vector<pair<int, int> > v, int start, int sum, int delta, int target);
};

#endif