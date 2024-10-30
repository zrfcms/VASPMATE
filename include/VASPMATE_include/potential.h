#include"read_write.h"

double Tersoff_describe(const char file[], double cutoff);

double bond_length(double atomi[3], double atomj[3]);

double bond_angle(double atomi[3], double atomj[3], double atomk[3]); // i->j,i->k,angle = i