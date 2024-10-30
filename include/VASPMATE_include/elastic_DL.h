#pragma once
#ifndef _ELASDL_
#define _ELASDL_
#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<string>
#include<unistd.h>
#include<functional>
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
namespace elas_DL
{
	enum LAUE
	{
		C,
		H,
		T,
		O,
		M,
		R
	};
    class elastic_DL
	{
	public:
		elastic_DL(int numDL = 0 , vector<double> strain_energyDL = {-0.018, -0.015, -0.012, -0.009, -0.006, -0.003, 0.000, 0.003, 0.006, 0.009, 0.012, 0.015, 0.018}, const char file[] = "INPOS");
		void generatefile(const char filename[], vector<vector<double> > matrix);
		void Recell(int numDL , const char file1[], const char file2[]);
		void Show_cell(double lattice[3][3], double position[][3],	const int types[],	const int num_atom);
		void getenergy(const char file[] = "Energy_Strain");
		void generate();
		void calculate();
		void PrintElasticProperty(const char file[] = "ELAS_INFO.dat");
		void cal_cub();
		void cal_hexa();
		void cal_trig();
		void cal_tetra();
		void cal_ortho();
		void cal_mono();
		vector<function<void(elastic_DL*)>> cal_method_energy;
		void matrix_cub(int j, double val);
		void matrix_hexa(int j, double val);
		void matrix_trig(int j, double val);
		void matrix_tetra(int j, double val);
		void matrix_ortho(int j, double val);
		void matrix_mono(int j, double val);
		//vector<function<void(elastic_energy*,int,double)>> cal_method_energy;
	private:
		POSCAR pos;
		LAUE Laue;
		vector<double> strain_energyDL;
		vector<vector<double> > defMat;
		vector<vector<double> > engxy;
		vector<vector<double> > dataxy;
		vector<vector<double> > Celas;
		vector<double> defVect;
		vector<double> conf_e2s;
		vector<int> lists;
		int nelastic;
		double vol;
		int numDL;
	};
};
#endif
