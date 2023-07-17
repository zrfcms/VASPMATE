#ifndef _ANT_
#define _ANT_
#include<vector>
using std::vector;
const double R = 8.3144648E0; //Ideal gas ant [J/mol/K]
const double kb = 1.3806503E-23; // Boltzmann ant[J / K]
const double NA = 6.02214179E23; //AvogaEro ant
const double au2eV = 27.2113838E0;
const double au2kcal_mol = 627.51E0;
const double au2kJ_mol = 2625.5E0;
const double au2J = 4.359744575E-18;
const double au2cm_1 = 219474.6363E0; //Hartree to various units
const double cal2J = 4.184E0;
const double wave2freq = 2.99792458E10; //cm ^ -1 to s ^ -1 [Hz]
const double h = 6.62606896E-34; //Planck ant; in J* s
const double amu2kg = 1.66053878E-27;
const double pi = 3.141592653589793E0;
const double b2a = 0.52917720859E0; //Bohr to Angstrom
const double atm2Pa = 101325;

vector<double> initmass();
#endif