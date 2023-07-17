#pragma once
#ifndef _DOS_
#define _DOS_

#include<stdlib.h>
#include<stdio.h>
#include<string>
#include<algorithm>
#include<string>
#include<vector>
#include<string.h>
#include<map>
using namespace std;
typedef struct COMB
{
	vector<int> index;
	vector<string> orb;
}COMB;
class DOS
{
public:
	void readDoscar(int ISPIN, int LORBIT);
	double fermi_energy();
	void GetDos(int argc, char* argv[]);
	void TDos();
	void ADos(int num, const char elemlabel[], int headind, int ISPIN, int LORBIT, int max_element);
	void EDos(vector<string> elem, int typenum[], char elemsym[][3], int nant[], int max_element, int ISPIN, int LORBIT);
	void SDos(vector<int> num, string postfix, int max_element, int ISPIN, int LORBIT);
	void ODos(vector<vector<string> > ind, vector<vector<string> > orb, vector<COMB> orbit, int typenum[], char elemsym[][3], int nant[],
		int max_element, int ISPIN, int LORBIT);
	void BandCenter(vector<int> num, string atomid, int max_element, int ISPIN, int LORBIT);
private:
	int IonWithSphere;
	int Ion;
	int Pdos;	// no 0 incl 1
	int NCDIJ;	// currently not used
	double volume;
	double vec[3];
	double POTIM;
	double init_T;
	string car;
	string system;
	double Emax;
	double Emin;
	int NEDOS;	// default 301
	double efermi;
	double dosnum;	// 1.0000
	vector<double> energy, dos_nospin, integ_dos_nospin, dos_up, dos_dw, integ_dos_up, integ_dos_dw;
	vector<vector<double> >
		//ISPIN = 1 LORBIT = 10 
		energy_ion, s_dos, p_dos, d_dos, f_dos,
		//ISPIN = 1 LORBIT = 11
		s0_dos, py_dos, pz_dos, px_dos, dxy_dos, dyz_dos, dz2_dos, dxz_dos, dx2_dos,
		f1_dos, f2_dos, f3_dos, f4_dos, f5_dos, f6_dos, f7_dos,
		//ISPIN = 2 LORBIT = 10
		s_dos_up, p_dos_up, d_dos_up, f_dos_up, s_dos_dw, p_dos_dw, d_dos_dw, f_dos_dw,
		//ISPIN = 2 LORBIT = 11
		s0_dos_up, py_dos_up, pz_dos_up, px_dos_up, dxy_dos_up, dyz_dos_up, dz2_dos_up, dxz_dos_up, dx2_dos_up,
		f1_dos_up, f2_dos_up, f3_dos_up, f4_dos_up, f5_dos_up, f6_dos_up, f7_dos_up,
		s0_dos_dw, py_dos_dw, pz_dos_dw, px_dos_dw, dxy_dos_dw, dyz_dos_dw, dz2_dos_dw, dxz_dos_dw, dx2_dos_dw,
		f1_dos_dw, f2_dos_dw, f3_dos_dw, f4_dos_dw, f5_dos_dw, f6_dos_dw, f7_dos_dw;
	// intergral of all orbits
	vector<vector<vector<double > > > orbit_dos;
	vector<vector<vector<double > > > Idos;
};

int IsRightOrbit(string orbit, int LORBIT, int max_element);

vector<int> TranArgvToAtomIndex(int argc, char* argv[], int nant[], int typenum[], char elemsym[][3], double xyz[][3], char label[]);

#endif
