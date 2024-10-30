#pragma once
#ifndef _ELAS_
#define _ELAS_
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<string>
#include<unistd.h>
#include<functional>
#include <fstream>
#include"../include/VASPMATE_include/tools.h"
#include"../include/VASPMATE_include/outcar.h"
#include"../include/VASPMATE_include/structure_operator.h"
#include"../include/VASPMATE_include/read_write.h"
#include"../spglib/spglib.h"
#include"../include/VASPMATE_include/structure_operator.h"
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::function;
namespace elas
{
	enum LAUE
	{
		C,
		H,
		R1,
		R2,
		T1,
		T2,
		O,
		M,
		N
	};
	class elastic
	{
	public:
		elastic(vector<double> strain = { -0.005, -0.003, -0.001, 0.000, 0.001, 0.003, 0.005 }, const char file[] = "INPOS");
		void generatefile(const char filename[], vector<vector<double> > matrix);
		void getstress(const char file[] = "Stress_Strain");
		void generate();
		void calculate();
		void PrintElasticProperty(const char file[] = "ELAS_INFO.dat");
		void cal_cub();
		void cal_hexa();
		void cal_trig1();
		void cal_trig2();
		void cal_tetra1();
		void cal_tetra2();
		void cal_ortho();
		void cal_mono();
		void cal_tric();
		vector<function<void(elastic*)> > cal_method;
	private:
		POSCAR pos;
		LAUE Laue;
		vector<double> strain;
		vector<vector<vector<double> > > stress;
		vector<vector<double> > coeff;
		vector<vector<double> > Celas;
		vector<int> lists;
	};
	class elastic_en
	{
	public:
		elastic_en(vector<double> strain_energy = {-0.018, -0.015, -0.012, -0.009, -0.006, -0.003, 0.000, 0.003, 0.006, 0.009, 0.012, 0.015, 0.018}, const char file[] = "INPOS");
		void generatefile(const char filename[], vector<vector<double> > matrix);
		void getenergy(const char file[] = "Energy_Strain");
		void generate();
		void calculate();
		void PrintElasticProperty(const char file[] = "ELAS_INFO.dat");
		void cal_cub();
		void cal_hexa();
		void cal_trig1();
		void cal_trig2();
		void cal_tetra1();
		void cal_tetra2();
		void cal_ortho();
		void cal_mono();
		void cal_tric();
		vector<function<void(elastic_en*)>> cal_method_energy;
		void matrix_cub(int j, double val);
		void matrix_hexa(int j, double val);
		void matrix_trig1(int j, double val);
		void matrix_trig2(int j, double val);
		void matrix_tetra1(int j, double val);
		void matrix_tetra2(int j, double val);
		void matrix_ortho(int j, double val);
		void matrix_mono(int j, double val);
		void matrix_tric(int j, double val);
		double ela_get_energy_oszi();
		//vector<function<void(elastic_energy*,int,double)>> cal_method_energy;
	private:
		POSCAR pos;
		LAUE Laue;
		vector<double> strain_energy;
		vector<vector<double> > defMat;
		vector<vector<double> > engxy;
		vector<vector<double> > dataxy;
		vector<vector<double> > Celas;
		vector<double> defVect;
		vector<double> conf_e2s;
		vector<int> lists;
		int nelastic;
		double vol;
	};
};
#endif
