#pragma once
#ifndef _KPTA_
#define _KPTA_
#include"read_write.h"
typedef struct matrix
{
	double mat[3][3];
}_matrix;

double K_volume(double a[3][3]);

_matrix GetVec(const char* poscar);

_matrix RecMat(_matrix POSI);

void cross(const double a[], const  double b[], double c[]);

int* kpta(_matrix rec, int kppra, int atom_num, char ikmesh);

void w_kpt(char ikmesh, int imesh[3]);

void E_ikppra(const char* file, char ikmesh, int ikppra);

int* kpvf(_matrix rec, double kspac, int atom_num, char ikmesh);

void E_ikspac(const char* file, char ikmesh, double kspac);

#endif