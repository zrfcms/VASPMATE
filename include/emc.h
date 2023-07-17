#pragma once
#ifndef _EMC_
#define _EMC_
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string>
#include<vector>
#include"Eigenval.h"
#include"read_write.h"
#include"bskpt.h"
#include"tools.h"
#include"Binomial.h"
#include"dos.h"

static const double Bohr = 0.5291772108;
static const double Hartree = 0.036749;

using std::vector;
using std::string;
class EMC
{
public:
	EMC(const char file[]);
	void get_emc_kpath(const char file[]);
	void get_emc(int band_index = -1);
	void get_emc_option(int band_index = -1);
	void read_level_emc(double efermi, double vec[3][3], int &cbm, int &vbm, int band_index = -1);
	vector<vector<double> > fitting(double efermi);
	void write_emc(int cbm, int vbm);
private:
	typedef struct path
	{
		double start[3];
		double end[3];
		char direct[100];
	}path;
	char ikmesh;
	double kspacing;
	int mode;
	int fit_point;
	double cutoff;
	int task;
	int ispin;
	vector<path> pa;
	vector<vector<vector<double> > > nkpts_emc;
	vector<vector<vector<vector<double> > > > level_emc; //2*ispin*task*fit_points
	vector<vector<double> > kstep_emc;
	vector<vector<double> >mass, r2;
};

#endif