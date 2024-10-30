#include"potential.h"

double bond_length(double atomi[3], double atomj[3])
{
	double length = 0;
	for (int i = 0; i < 3; i++)
		length += (atomi[i] - atomj[i]) * (atomi[i] - atomj[i]);
	return sqrt(length);
}

double bond_angle(double atomi[3], double atomj[3], double atomk[3]) // i->j,i->k,angle = i
{
	double Rij[3], Rik[3];
	double rij = 0, rik = 0, rijk = 0;
	for (int i = 0; i < 3; i++)
	{
		Rij[i] = atomi[i] - atomj[i];
		Rik[i] = atomi[i] - atomk[i];
		rij += Rij[i] * Rij[i];
		rik += Rik[i] * Rik[i];
		rijk += Rij[i] * Rik[i];
	}
	return rijk / sqrt(rij) / sqrt(rik);
}

double Tersoff_describe(const char file[], double cutoff)
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return -1;
	}
	POSCAR pos;
	readposcar(fp, pos);
	direct_to_carts(pos.iflg, pos.vec, pos.xyz, pos.nant);
	fclose(fp);
	double describe = 0;
	double m = 3;
	double gamma = 1;
	double lambda3 = 1;
	double c = 19981;
	double d = 7.034;
	double costheta0 = -0.33953;
	for (int i = 0; i < pos.nant[0]; i++)
	{
		for (int j = 0; j < pos.nant[0]; j++)
		{
			for (int k = j + 1; k < pos.nant[0]; k++)
			{
				if (i == j || i == k || j == k)
					continue;
				double length_ij = bond_length(pos.xyz[i], pos.xyz[j]);
				double length_ik = bond_length(pos.xyz[i], pos.xyz[k]);
				if (length_ik > cutoff || length_ij > cutoff)
					continue;
				double theta_ijk = bond_angle(pos.xyz[i], pos.xyz[j], pos.xyz[k]);
				double g = gamma * (1 + c * c / d / d - c * c / (d * d + (cos(theta_ijk) - costheta0) * (cos(theta_ijk) - costheta0)));
				describe += g * exp(pow(lambda3 * (length_ij - length_ik), 3));
			}
		}
	}
	return describe;
}