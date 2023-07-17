#pragma once
#ifndef _SPAPLOT_
#define _SPAPLOT_
#include"sv_Chart.h"
#include"outcar.h"
#include<string>
#include<stdio.h>
#include<stdlib.h>
#include<algorithm>
using std::string;
using std::min_element;
using std::max_element;
string convert_name(const char file[]);
string convert_TDOS(const char file[], int ISPIN);
string convert_ADOS(const char file[]);
string convert_ADOS(const char file_up[], const char file_dw[]); 
string convert_basic_band(const char file[], int ISPIN);
string convert_fat_band(const char file[]);
double get_Enthalpy_formation(vector<double> en);
class hull_point
{
public:
	hull_point() {}
	static bool mysort(hull_point p1, hull_point p2) { return p1.x < p2.x; };
	void spa_convexhull(const char data[] = "data");
private:
	double x, y;
};
#endif