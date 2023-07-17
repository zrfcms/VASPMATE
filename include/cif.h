#pragma once
#ifndef _CIF_
#define _CIF_

#include<stdlib.h>
#include<stdio.h>
#include<vector>
#include<string>
#include<queue>
#include<stack>
#include<iostream>
#include<algorithm>
#include<sstream>
#include<map>
#include"read_write.h"
using namespace std;
typedef struct CIF
{
	double length_a, length_b, length_c;
	double alpha, beta, gamma;
	vector<string> elemsymbol;
	vector < vector<vector<string> >>frac_xyz;
	vector<vector < vector<double> > > pos_xyz; //elem*typenum*3
	vector<vector<string> >  equiv_pos;
	vector<int> cif_typenum;
	vector<int> pos_typenum;
	int pos_atom_num;
	double vec[3][3];
}CIF;

void CIF_init(CIF& cif);

void TranCellToVec(CIF& cif);

template<typename T>
T TranEquivToCal(string input, string equiv);

void getpos_xyz(CIF& cif);

void readcif(FILE* fp, CIF& cif);

void printinfo(CIF cif);

void TranCifToPOS(CIF cif, POSCAR& pos);

void TranCIFToPOSCAR(const char cif[], const char poscar[]);

void TranPOSCARToCIF(const char poscar[], const char cif[]);

#endif