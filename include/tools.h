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
using namespace std;
string GetInfoINCAR(const char key[]);

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

int zdgys(int a, int b);

string TranDecimalTofraction(string x);

string RemovePostfixNum(string elem);

class Node {
public:
	double num ;
	char op ;
	bool isop;

	Node() {
		num = 0;
		op = '+';
	};
	Node(int num, bool isop) {
		this->num = num;
		this->isop = false;
	};
	Node(char op, bool isop) {
		this->op = op;
		this->isop = true;
	}
};

vector<Node> midToBack(string mid);

double getValue(vector<Node> back);

template<typename T>
bool IsSameVector(vector<T> v1, vector<T> v2)
{
	if (v1.size() != v2.size())
		return false;
	for (int i = 0; i < v1.size(); i++)
	{
		if (fabs(v1[i] - v2[i]) >= 1e-4)
			return false;
	}
	return true;
}

void copy_file(const  char* target, const  char* source);

int len(char buf[]);

void getlength(double vec[3][3], double length[3]);

void brinv(double vec[3][3], double inv_vec[3][3]);

void transpose_matrix(vector<vector<double> > v, vector<vector<double> >& tran_v);

vector<vector<double>> mutil(vector<vector<double>> m1, vector<vector<double>> m2);

void cross(double a[3], double b[3], double c[3]);

void ITP_french(int init_size[3], vector<double> init_pot, int out_size[3], vector<double> &out_pot, int scale);
#endif