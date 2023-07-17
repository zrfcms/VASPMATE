#pragma once
#ifndef _BAND_
#define _BAND_

#include<stdlib.h>
#include<stdio.h>
#include<string>
#include<algorithm>
#include<string>
#include<vector>
#include<string.h>
#include<map>
#include"read_write.h"
#include"emc.h"
#include"tools.h"
#include"Binomial.h"
#include"dos.h"
#include"spa_plot.h"
using namespace std;
class BAND
{
public:
	void readPROCAR(int ISPIN, int LORBIT);
	void getband(int argc, char* argv[]);
	vector<double> getkpointsline(double vec[3][3]);
	vector<double> hse_getkpointsline(int nkpts, vector<vector<double> > kcoor, double vec[3][3]);
	void basicband(int ISPIN, vector<double> kline, double efermi);
	void Aband(int num, const char elemlabel[], int headind, int ISPIN, int LORBIT, int max_element, vector<double> kline, double efermi);
	void Eband(vector<string> element, int typenum[], char elemsym[][3], int nant[], int max_element,
		int ISPIN, int LORBIT, vector<double> kline, double efermi);
	void Sband(vector<int> num, string postfix, int max_element, int ISPIN, int LORBIT, vector<double> kline, double efermi);
	void Oband(vector<int> num, vector<int> orbit_ind, vector<double> kline, int ISPIN, int LORBIT, double efermi);
	vector<vector<int> > BandGap(int ISPIN, vector<double> kline, double efermi, bool write = 1);
	double Effective_Mass(int fitting_point, int position, int band_index, bool direct, vector<double> kline, double efermi);
	vector<vector<vector<vector<double> > >	> getion_dos() { return this->ion_dos; };
	vector<vector<vector<vector<double> > >	> getion_dos_up() { return this->ion_dos_up; };
	vector<vector<vector<vector<double> > >	> getion_dos_dw() { return this->ion_dos_dw; };
	int getibz_kpt() { return this->start_point; };
private:
	string title;
	int num_kpts;
	int num_bands;
	int num_ions;
	vector<vector<double> > k_coor; // kpoints coordinates num_kpts*3
	vector<double> weight; // num_kpt
	vector<vector<double> > energy; // num_kpts*num_bands
	vector<vector<double> > occ; // num_kpts*num_bands
	vector<vector<vector<vector<double> > >	> ion_dos; // num_kpts*num_bands*(num_ion+1)*orbit
	vector<vector<double> > energy_up;
	vector<vector<double> > energy_dw;
	vector<vector<double> > occ_up;
	vector<vector<double> > occ_dw;
	vector<vector<vector<vector<double> > >	> ion_dos_up;
	vector<vector<vector<vector<double> > >	> ion_dos_dw;
	int start_point;
};

class point
{
public:
	point(double pos[3], char label[10])
	{
		this->pos.resize(3);
		for (int i = 0; i < 3; i++)
			this->pos[i] = pos[i];
		this->label = label;
	}
	vector<double> pos;
	string label;
};
#endif