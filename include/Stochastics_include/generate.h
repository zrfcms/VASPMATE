#pragma once
#ifndef _GEN_
#define _GEN_
#include<stdio.h>
#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<sys/wait.h>
#include<dirent.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string>
#include<vector>
#include<spglib/spglib.h>
#include<algorithm>
#include"read_write.h"
#include"kpta.h"
#include"variant.h"
#include"tools.h"
using namespace std;
typedef struct
{
	char name[MAX_FNAME];//name of the stucture
	double energy;
	double volume;
	double formation_energy;
	double addition_data;
	double weight; // weight = w1*energy + w2*addition_data
	int method;//how we get the structure
}Account;

typedef struct
{
	int num;
	double(*pos)[3];
	int* wyckoff;
	int* equivalent_id;
	int* equivalent_num;
}SK_site_list;

string GetValueAfterEqual(char buf[]);

void read(const char poscar[], int atom[], int& num);

void write_POTCAT(int atom[], int num);

void gene_first_population();

void getbox(int mode, int R_flag, double lattice[3][3], double ia, double ib, double ic, double ialpha, double ibeta, double igamma);

void set_lattice(int type, double& a, double& b, double& c, double& alpha, double& beta, double& gamma);

void SK_malloc_sitelist(int n, SK_site_list* list);

void SK_get_sitelist_3D(SK_site_list* list,
	const int sym_size,
	int(*rot)[3][3],
	const double(*trans)[3],
	SPGCONST double bravais_lattice[3][3],
	const int hall_number,
	double symprec);

//vector<vector<double> > get_equiv_atoms(POSCAR &pos, int hall_number, int(*rot)[3][3], double(*tra)[3], double x, double y, double z, double symprec);

void gene_Random_structrue(vector<int> atom, vector<string> element, int formula_num, const char* file, double volume);

vector<string> gene_first_population_Random(vector<int> atom, vector<string> element, int POPULATION_NUM, int formula_num, double volume);

#endif

