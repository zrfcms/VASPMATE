#pragma once
#ifndef _CHGCAR_
#define _CHGCAR_

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<string.h>
#include"tools.h"
#include"read_write.h"
#include"thermo_constants.h"
using namespace std;
class CHGCAR
{
public:
	void readchgcar(const char file[]);
	void operator_bader(int argc, char* argv[]);
	CHGCAR operator_cdd(CHGCAR chg1, CHGCAR chg2, bool oper);  // 0:- 1:+
	void series_oper(int argc, char* argv[], bool oper);
	void writechgcar(const char file[], bool label);
	void writespincar(const char file[]);
	void writespinUp_Dwcar(const char up_file[], const char dw_file[]);
	void ChgcarToCube(const char file1[], const char file2[]);
	POSCAR getchgcarpos() { return this->pos; }
	void getchg_size(int _chg_size[3]) {
		for (int i = 0; i < 3; i++)
			_chg_size[i] = this->chg_size[i];
	}
	vector<double> getchg() { return this->chg; }
	vector<vector<vector<double> > > getChg() { return this->Chg; }
	vector<vector<double> > getchg_aug() { return this->chg_aug; }
	void get_Separator(char _Separator[100]) { strcpy(_Separator, this->Separator); }
	void getmag_size(int _mag_size[3]) {
		for (int i = 0; i < 3; i++)
			_mag_size[i] = this->mag_size[i];
	}
	vector<double> getmag() { return this->mag; }
	vector<vector<double> > getmag_aug() { return this->mag_aug; }
private:
	POSCAR pos;
	//in order of x y z
	int chg_size[3];
	//in order of z y x 
	vector<double> chg;
	vector<vector<vector<double> > > Chg;
	vector<vector<double> >chg_aug;
	char Separator[100];
	int mag_size[3];
	vector<double> mag;
	vector<vector<double> >mag_aug;
};

void writechgcar(const char file[], POSCAR pos, int chg_size[3], vector<double> chg, const vector<vector<double> >& chg_aug = vector<vector<double> >(),
	char Separator[100] = NULL, int mag_size[3] = NULL, const vector<double>& mag = vector<double>(),
	const vector<vector<double> >& mag_aug = vector<vector<double> >(), bool label = 0);

void Chgcarinterpolation(const char file1[], const char file2[], int scale = 2);

class Cube
{
public:
	Cube(const char file[]);
	void CubeToChgcar(const char file[], int orbit_num);
private:
	int atom_num;
	double org[3]; // origin point
	double vec[3][3];
	int cub_size[3];
	vector<vector<double> > xyz;
	vector<int> elemnum;
	vector < vector<int> >typenum;
	vector<double> charge;
	int nmo; // number of orbit
	vector<int> orbit_index;
	vector<vector<vector<vector< double> > > >  cub; //orbit_num*i*j*k
};
#endif