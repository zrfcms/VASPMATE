#pragma once
#ifndef _STRUCTURE_
#define _STRUCTURE_
#include<spglib/spglib.h>
#include<read_write.h>
#include<string>
#include<vector>
#define SYMPREAC 1e-4
using namespace std;
void EV_cross(const double a[], const  double b[], double c[]);

void get_unitvector(double vec[3][3], double dx[3], double dy[3], double dz[3]);

static void show_cell(double lattice[3][3],
	double position[][3],
	const int types[],
	const int num_atom);

static void show_spacegroup_type(const SpglibSpacegroupType spgtype);

static int show_spg_dataset(double lattice[3][3],
	const double origin_shift[3],
	double position[][3],
	const int num_atom,
	const int types[]);

void transpose_matrix(double vec[3][3], double n_vec[3][3]);

void translate_typenum_type(int* type, int atom_type, int typenum[]);//typenum->type

void translate_type_typenum(int* type, int nant[], int typenum[]);//type->typenum

void translate_typenum_type(int* type, int atom_type, int typenum[]);

int EV_primitive(const char file1[], const char file2[]);

int EV_unitcell(const char file1[], const char file2[]);

void EV_get_symmetry(const char file[]);

void EV_supercell(const char file1[], const char file2[], int super[3]);

void get_defvect(char mode[], double current_strain, double defvect[3][3]);

void EV_affine(const char file1[], char model[], char mode[], double init_value, double step_length, int step_num);

void EV_alias_tensile(const char file[], char mode[], double istart, double iend, double ispacing, double value);

void EV_alias_shear(const char file[], char mode[], double istart1, double iend1, double ispacing1,
	double istart2, double iend2, double ispacing2, double value);

void EV_rotproj(const char file1[], const char file2[], double rot[3]);

void EV_indproj(const char file1[], const char file2[], int pvh, int pvk, int pvl, int uvu, int uvv, int uvw);

void EV_matproj(const char file1[], const char file2[], double projmat[3][3]);

void EV_redefine(const char file1[], const char file2[], int rot[3][3], double symprec, int a1, int x1, int a2, int x2, bool fixdirect);

void EV_recell(const char file1[], const char file2[]);

void EV_fixatomcoor(const char file1[], const char file2[], char axis[], double start, double end, char addfix[3]);

void EV_fixatomindex(const char file1[], const char file2[], int argc, char* argv[], char addfix[3]);

void EV_fixatomele(const char file1[], const char file2[], vector<string> elem, char addfix[3]);

void EV_cleanfix(const char file[]);

void EV_carts_to_direct(const char file1[], const char file2[]);

void EV_direct_to_carts(const char file1[], const char file2[]);

void EV_atomsort_coord(const char file1[], const char file2[], char mode[]);

void EV_atomsort_element(const char file1[], const char file2[], char ele1[], char ele2[]);

void EV_move(const char file1[], const char file2[], char label[], char axis[], double move[3], double mmin, double mmax);

void adsorbent(char basic[], char clus[], char destin[], double dis);
#endif
