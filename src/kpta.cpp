#include"../include/kpta.h"
int atom_num = 0;

_matrix GetVec(const char* poscar)
{
	char   title[MAX_FNAME];                   // title of the system
	double latt;                     // Scaling Lattice parameter 
	int    ifix;                     // Select 0 or normal 1
	int    iflg;                     // Direct 0 or Cartes 1
	int    nant[2];                      // nant[0];nant[1] number of atoms and types
	int    typenum[MAX_NELEM];                 // number of atoms for each ion type
	int    elemnum[MAX_NELEM];                 // atomic number of each element
	char   elemsym[MAX_NELEM][3];              // atomic symbol of each element
	double vec[3][3];                 // lattice vectors
	double xyz[MAX_NATOM][3];                  // atomic coordinates
	char   fix[MAX_NATOM][3];                   // Selective fix on each atom
	FILE* fp = fopen(poscar, "r");
	if (fp == NULL) {
		printf("%s cannot open \n", poscar);
		exit(1);
	}
	readposcar(fp,
		title,                   // title of the system
		latt,                     // Scaling Lattice parameter 
		ifix,                     // Select 0 or normal 1
		iflg,                     // Direct 0 or Cartes 1
		nant,                      // nant[0],nant[1] number of atoms and types
		typenum,                 // number of atoms for each ion type
		elemnum,                 // atomic number of each element
		elemsym,              // atomic symbol of each element
		vec,                 // lattice vectors
		xyz,                  // atomic coordinates
		fix);                 // Selective fix on each atom
	fclose(fp);
	atom_num = nant[0];
	_matrix ret;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ret.mat[i][j] = vec[i][j] * latt;
	return ret;
}

double K_volume(double a[3][3])
{
	double odd = 0;
	double even = 0;
	for (int i = 0; i < 3; i++)
		odd += a[0][i] * a[1][(i + 1) % 3] * a[2][(i + 2) % 3];
	for (int i = 0; i < 3; i++)
		even += a[0][i] * a[1][(i + 2) % 3] * a[2][(i + 1) % 3];
	return fabs(odd - even);
}

_matrix RecMat(_matrix POSI)
{
	double v = K_volume(POSI.mat);
	if (v == 0)
		v = 1;
	_matrix ret;
	ret.mat[0][0] = (POSI.mat[1][1] * POSI.mat[2][2] - POSI.mat[1][2] * POSI.mat[2][1]) * (2 * Pi) / v;
	ret.mat[0][1] = -(POSI.mat[1][0] * POSI.mat[2][2] - POSI.mat[1][2] * POSI.mat[2][0]) * (2 * Pi) / v;
	ret.mat[0][2] = (POSI.mat[1][0] * POSI.mat[2][1] - POSI.mat[1][1] * POSI.mat[2][0]) * (2 * Pi) / v;
	ret.mat[1][0] = -(POSI.mat[0][1] * POSI.mat[2][2] - POSI.mat[0][2] * POSI.mat[2][1]) * (2 * Pi) / v;
	ret.mat[1][1] = (POSI.mat[0][0] * POSI.mat[2][2] - POSI.mat[0][2] * POSI.mat[2][0]) * (2 * Pi) / v;
	ret.mat[1][2] = -(POSI.mat[0][0] * POSI.mat[2][1] - POSI.mat[0][1] * POSI.mat[2][0]) * (2 * Pi) / v;
	ret.mat[2][0] = (POSI.mat[0][1] * POSI.mat[1][2] - POSI.mat[1][1] * POSI.mat[0][2]) * (2 * Pi) / v;
	ret.mat[2][1] = -(POSI.mat[0][0] * POSI.mat[1][2] - POSI.mat[1][0] * POSI.mat[0][2]) * (2 * Pi) / v;
	ret.mat[2][2] = (POSI.mat[1][1] * POSI.mat[0][0] - POSI.mat[1][0] * POSI.mat[0][1]) * (2 * Pi) / v;
	return ret;
}

void cross(const double a[], const  double b[], double c[])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

int* kpta(_matrix rec, int ikppra, int atom_num, char ikmesh, bool has_vacuum_slab[3])
{
	double N[3] = { 0 };
	double dotnorm[3] = { 0,0,0 };
	double crossvec[3][3];
	//crossvec:the cross product of rec
	for (int i = 0; i < 3; i++)
		cross(rec.mat[(i + 1) % 3], rec.mat[(i + 2) % 3], crossvec[i]);
	//dotnorm:the volume of rec
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			dotnorm[i] += rec.mat[i][j] * crossvec[i][j];
			N[i] += crossvec[i][j] * crossvec[i][j];
		}
	}
	double proj[3];
	for (int i = 0; i < 3; i++)
	{
		proj[i] = dotnorm[i] / sqrt(N[i]);
	}
	static int imesh[3] = { 0 };
	double fmesh[3] = { 0 };
	ikppra /= atom_num;
	double normalizer = pow(proj[0] * proj[1] * proj[2] / ikppra, 1.00 / 3.00);
	for (int i = 0; i < 3; i++)
	{
		imesh[i] = (int)(proj[i] / normalizer + 0.5);
		fmesh[i] = proj[i] / normalizer - (double)imesh[i];
	}
	if (ikmesh == 'M' || ikmesh == 'm')
	{
		for (int i = 0; i < 3; i++)
			if (imesh[i] % 2 == 1)
			{
				imesh[i]++;
				fmesh[i] = 0.00;
			}
	}
	double zero_tolerance = 1E-4;
	while ((double)(imesh[1] * imesh[2] * imesh[0]) < ikppra)
	{
		double temp = 0;
		for (int i = 0; i < 3; i++)
			if (fmesh[i] > temp)
				temp = fmesh[i];
		for (int i = 0; i < 3; i++)
			if (fabs(fmesh[i] - temp) < zero_tolerance)
			{
				imesh[i]++;
				if (ikmesh == 'M' || ikmesh == 'm')
					imesh[i]++;
				fmesh[i] = 0.000;
			}
	}
	for (int i = 0; i < 3; i++)
	{
		if (has_vacuum_slab[i] == 1)
			imesh[i] = 1;
	}
	return imesh;
}

int* kpvf(_matrix rec, double kspac, int atom_num, char ikmesh, bool has_vacuum_slab[3])
{
	static int imesh[3];
	double reciplgth[3] = { 0 };
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			reciplgth[i] += rec.mat[i][j] * rec.mat[i][j];
		}
		reciplgth[i] = sqrt(reciplgth[i]);
		imesh[i] = (int)(reciplgth[i] / kspac + 0.5);
		if (imesh[i] < 1)
			imesh[i] = 1;
		if (imesh[i] % 2 == 1)
			imesh[i]++;
	}
	if (ikmesh == 'M' || ikmesh == 'm')
	{
		for (int i = 0; i < 3; i++)
			if (imesh[i] % 2 == 1)
			{
				imesh[i]++;
			}
	}
	for (int i = 0; i < 3; i++)
	{
		if (has_vacuum_slab[i] == 1)
			imesh[i] = 1;
	}
	return imesh;
}

void w_kpt(char ikmesh, int imesh[3])
{
	FILE* fp = fopen("NEWKPT", "w");
	fprintf(fp, "Automatic mesh \n");
	fprintf(fp, "0\n");
	if (ikmesh == 'm' || ikmesh == 'M')
		fprintf(fp, "%s\n", "Monkhorst-pack");
	else
		fprintf(fp, "%s\n", "Gamma-centeredk");
	fprintf(fp, "%d %d %d \n", imesh[0], imesh[1], imesh[2]);
	fprintf(fp, "0 0 0 \n");
	fclose(fp);
	printf("Written NEWKPT file!\n");
}

void E_ikppra(const char* file, char ikmesh, int ikppra)
{
	_matrix vec = GetVec(file);
	_matrix rec = RecMat(vec);
	FILE* fp = fopen(file, "r");
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	bool has_vacuum_slab[3];
	has_vacuum_slab_or_not(pos, has_vacuum_slab);
	int* imesh = kpta(rec, ikppra, atom_num, ikmesh, has_vacuum_slab);
	w_kpt(ikmesh, imesh);
}

void E_ikspac(const char* file, char ikmesh, double kspac)
{
	_matrix vec = GetVec(file);
	_matrix rec = RecMat(vec);
	FILE* fp = fopen(file, "r");
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	bool has_vacuum_slab[3];
	has_vacuum_slab_or_not(pos, has_vacuum_slab);
	int* imesh = kpvf(rec, kspac, atom_num, ikmesh, has_vacuum_slab);
	w_kpt(ikmesh, imesh);
}
