#pragma once
#ifndef _BSKPT_
#define _BSKPT_

#include"structure_operator.h"
#include"kpta.h"
#include"tools.h"
#include<spglib/spacegroup.h>
#include<spglib/spg_database.h>
#include<spglib/spglib.h>
#include<spglib/debug.h>
#include<vector>
#include<algorithm>
#include<string>
#include<map>
#include<vector>
#include<iostream>
using namespace std;
typedef enum
{
	TYPE_WRONG,
	CUB,//CUB cp
	FCC,//Face-centered cubic cF
	BCC,//Body-centered cubic cI
	TET,//Tetraonal tP
	BCT,//Body-centered tetragonal tI
	ORC,//Orthorhombic oP
	ORCF,//Face-centered orthorhombic oF
	ORCI,//Body-centered orthorhombic oI
	ORCC,//C-centered orthorhombic os
	HEX,//Hexagonal hP
	RHL,//Rhombohedral hR
	MCL,//Monoclinic mP
	MCLC,//C-centered monoclinic ms
	TRL,//Triclinic ap
}Bravais;

class Symmpoint
{
public:
	Symmpoint() {}
	Symmpoint(double b1, double b2, double b3, string symbol)
	{
		this->b1 = b1;
		this->b2 = b2;
		this->b3 = b3;
		this->symbol = symbol;
	}
	double b1;
	double b2;
	double b3;
	string symbol;
};

void get_cell_params(double vec[3][3], double& a, double& b, double& c, double& alpha, double& beta, double& gamma);

Cell* Transform_to_Primitive(const Cell* cell,
	SPGCONST double trans_mat[3][3],
	const Centering centering,
	int spacegroup_number,
	const double symprec);

static int Find_Primitive_Cell(double lattice[3][3],
	double position[][3],
	int types[],
	const int num_atom,
	const double symprec,
	double angle_tolerance);

void bandskpt_3d(const char* file1, const char* file2, int kppra);

Bravais find_bravais(Centering centering, int spacegroup_number);

vector<Symmpoint> get_path(string path, map<char, Symmpoint> path_point);

void write_kpt(vector<Symmpoint> highsymmpoint, int kppra, bool time_reversal);

void write_highsymmpoint(vector<Symmpoint> point, string path, map<char, Symmpoint> path_point,bool time_reversal);

int* ir_reciprocal_mesh(const char* file, int ikppra, double kspac, char ikmesh, int kpt[3], char label[], vector<int> &weight, vector<vector<double> > &irr_coor);

void HSE_mesh(const char file[], int ikppra, double kspac, int kpt[], char label[], double kpra, char ikmesh);

vector<Symmpoint> find_highsymmpoint(Bravais bravais, double vec[3][3], int spacegroup_num, bool time_reversal);

vector<Symmpoint> get_2Dpath(double a1, double a2, double gamma);

void bandskpt_2d(const char* file1, const char* file2, int kppra);

void GetFermiMesh(const char file[], int kppra, double kspac, int mesh[3], char label[], char ikmesh);

void Get3DbandMesh(const char file[], int kppra, double kspac, int kpt[3], char label[], char ikmesh);

void getunfoldkpoints(const char file[], int ikppra, double kspac, int kpt[], char label[], double kpra, char ikmesh,const char tranfile[]);

void get_tran_matrix(const char tranfile[], double tran[3][3]);
#endif