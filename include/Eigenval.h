#pragma once
#ifndef _EIGENVAL_
#define _EIGENVAL_

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<set>
#include"plotWave.h"
#include"tools.h"
#include"structure_operator.h"
#include"read_write.h"
#include"dos.h"
#include"bskpt.h"
#include"band.h"
using std::vector;

class EIGENVAL
{
public:
	EIGENVAL(const char file[], double efermi);
	void TranEigenToXcrysden(set<int> select_band_index, double efermi);
	void TranEigenToFermiSurface(double efermi, string LORBIT, vector<vector<vector<vector<double> > > > ion_dos, 
	vector<vector<vector<vector<double> > >> ion_dos_up, vector<vector<vector<vector<double> > > > ion_dos_dw, int argc, char* argv[]);
	void get3dband(vector<int> band_index, double efermi);
	void read_kpts_mapping_table(vector<vector<double> > &kmesh_primcell, vector<vector<double> > &kmesh_supercell, vector<vector<int> > &gvector_supercell);
	void getunfold();
	int get_ispin() { return ISPIN; };
	int get_bands() { return num_bands; };
	int get_cbm() { return cbm; };
	int get_vbm() { return vbm; };
	vector<vector<vector<double> > > get_eigenvalue() { return eigenvalue; };
	vector<double> getkpointsline(double vec[3][3]);
private:
	int vbm, cbm;
	int ISPIN;
	int value_electron, num_kpts, num_bands;
	vector<vector<double>> position;
	vector<double> weight;
	vector<vector<vector<double> > > eigenvalue;  //ispin*num_kpts*num_bands
	set<int> band_index_up, band_index_dw;
};

void get_overlap_G_single_kpt(int npws_single_kpt, ionizing::MatrixX3i G_single_kpt, vector<int>& index_G_primcell, int& ng, double tran[3][3]);
#endif