#ifndef _SHERMO_
#define _SHERMO_
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<Eigen/Eigen>
#include"thermo_constants.h"
#include"read_write.h"
#include "msym.h"
using namespace Eigen;
using namespace std;

int TranposcarTomsym_element_t(POSCAR &pos, msym_element_t** element);

class init
{
public:
	init(const char poscar[], const char outcar[]);
	vector<double> getmass() { return mass; };
	char* getpointgroup() { return pointgroup; };
	inline double getE() { return E; }
	vector<double> gettot_wavenum() { return tot_wavenum; };
	vector < vector<double> > getxyz() { return xyz; };
private:
	vector<double> mass;
	vector<vector<double> > xyz;
	vector<double> tot_wavenum;
	char pointgroup[6];
	double E; //T = 0K energy
};

class InertiaValue
{
public:
	InertiaValue(vector<double> mass, vector<vector<double>> xyz);
	InertiaValue() { totmass = 0; }
	double gettotmass() { return totmass; };
	VectorXcd getinert() { return inert; };
private:
	Matrix3d inertmat;
	VectorXcd inert;
	double totmass;
};

class thermol :public InertiaValue
{
public:
	thermol(init _in, int imode, int p, double T, int ilowfreq,double ravib, int nSpinmulti)
		:InertiaValue(_in.getmass(), _in.getxyz()) {
		this->imode = imode; this->p = p; this->T = T; ilinear = getilinear();
		this->ilowfreq = ilowfreq; this->ravib = ravib; this->nSpinmulti = nSpinmulti;
		U = S = Q = H = G = Cv = Cp = 0;
		sclZPE = sclheat = sclS = sclCV = 1;
		rotsym = getrotsym(_in.getpointgroup());
	}
	thermol() {
		this->imode = 0; this->p = 0; this->T = 0; ilinear = 0;
		U = S = Q = H = G= Cv = Cp = rotsym = ilowfreq = ravib = 0;
		sclZPE = sclheat = sclS = sclCV = 1;
	}
	inline int getimode() { return imode; }
	inline int getp() { return p; }
	inline double getT() { return T; }
	inline int getrotsym() { return rotsym; };
	inline int getilowfreq() { return ilowfreq; };
	inline double getravib() { return ravib; };
	inline int getilinear()
	{
		VectorXcd inert = getinert();
		for (int i = 0; i < 3; i++)
			if (inert[0].real() < 0.001)
				return 1;
		return 0;
	}
	int getrotsym(char pointgroup[6]);
	int getnSpinmulti() { return nSpinmulti; };
	inline int getsclZPE() { return sclZPE; }
	inline int getsclheat() { return sclheat; }
	inline int getsclS() { return sclS; }
	inline int getsclCV() { return sclCV; }
public:
	double U, S, Q, H, G, Cv, Cp;
private:
	double sclZPE , sclheat, sclS, sclCV ;
	int imode; // Mode of evaluating thermodynamic quantities. 0: Consider all terms. 1: Ignore translation and rotation
	int p; // P = p * Pa
	double T; // temperature
	int ilinear; //1:linear 0:nolinear
	int rotsym; //rotation symmetry number
	int ilowfreq; //Treatment of low frequencies. 0: Harmonic. 1: Raising low frequencies. 2: Grimme's entropy interpolation
	double ravib; //Raising lower frequencies to this value (cm^-1) when ilowfreq=1
	int nSpinmulti;//spin multiplicity
};

void getthermol_data(thermol& th, vector<double> tot_wavenum,double E);

void shermo(int argc, char* argv[]);
#endif
#pragma once
