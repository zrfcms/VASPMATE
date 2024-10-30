#pragma once
#define _POSIX_C_SOURCE 2
#ifndef _MAGN_
#define _MAGN_

#define MAX_ATOM_MAG    100   // max number of atoms in amag
#define MIN_EPSILON     1e-6

#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<cmath>
#include<regex>
#include<unordered_map>
#include"../spglib/spglib.h"
#include"../include/VASPMATE_include/read_write.h"
#include"../include/VASPMATE_include/structure_operator.h"
#include"../include/VASPMATE_include/incar.h"

typedef struct MAGTYPE
{
	vector< vector<string > > order_mode;  // different ordering strategies to use, like:
	//ferromagnetic, antiferromagnetic, antiferromagnetic_by_motif,
    //ferrimagnetic_by_motif and ferrimagnetic_by_species (here, "motif",
    //means to use a different ordering parameter for symmetry inequivalent sites)
	//Cited from: High-throughput prediction of the ground-state collinear magnetic order of inorganic materials using Density Functional Theory
	vector< vector<double> > order_mag;  // MAGMON in INCAR
}MAGTYPE;


class magorder
{
public:
	magorder(const char file[] = "INPOS",const char file2[] = "INCAR", const char table[] = nullptr);
	void generate(char mode[] = nullptr, int cell_supernumb = 1);
	void incar_generate(char mode[] = nullptr);
	void derive();
	void derive_POS(string folderName);
	void read_pos_mag(const char file[] = "INPOS", char mode[] = nullptr);
	//void write_file_mag();
	void generatefile_MAGMON_list(FILE* fp_wri, string j_number, string mode, vector<double> magmon, int nant_0);
	void generatefile_INCAR(string mode , vector<double> magmon, int nant_0);
	void generatefile_INCAR_imag(string mode, string j_number, vector<double> magmon, int nant_0);
	void cal_order_mag(POSCAR pos);
	void cal_order_mode(POSCAR pos);
	//void ferromagnetic();
	//void ferrimagnetic_wyckoff();
	//void ferrimagnetic_element();
	//void antiferromagnetic();
	//void antiferromagnetic_wyckoff();
	vector<int> Similarity_check(vector<vector<double>> check_vector);

private:
	MAGTYPE magtype;
	int diffposition = 0; //Whenever the Wyckoff position and the element are not the same, one will be added.
    static const int MY_num_orderings = 64;
    static const int MY_max_unique_sites = 8;
	vector<string> amag_wyckoffs;
	vector<string> amag_ele;
	vector<int> amag_equivalent_atoms;
	vector<vector<int>> combination;
};

#endif