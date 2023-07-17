#include"structure_operator.h"
#include"tools.h"
#include<algorithm>
#pragma warning(disable:4996)
#define ZERO_TOLERANCE 1e-4

using namespace std;
char   title[MAX_NCHAR];        // Title of system
double latt;                    // Scaling factor
int    ifix;                    // Select 0 or normal 1
int    iflg;                    // Direct 0 or Cartes 1
int    nant[2];                 // number of atoms and types
int    typenum[MAX_NELEM];      // number of atoms of each type
int    elemnum[MAX_NELEM];      // atomic number
char   elemsym[MAX_NELEM][3];   // atomic symbol
double vec[3][3];               // lattice vector
double xyz[MAX_NATOM][3];       // atomic coordinate
char   fix[MAX_NATOM][3];       // Selective fix on each atom

void EV_cross(const double a[], const  double b[], double c[])
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

void get_unitvector(double vec[3][3], double dx[3], double dy[3], double dz[3])
{
	double norm1, norm2, norm3;
	for (int i = 0; i < 3; i++)
	{
		norm1 += vec[0][i] * vec[0][i];
		norm2 += vec[1][i] * vec[1][i];
		norm3 += vec[2][i] * vec[2][i];
	}
	for (int i = 0; i < 3; i++)
	{
		dx[i] = vec[0][i] / sqrt(norm1);
		dy[i] = vec[1][i] / sqrt(norm2);
		dz[i] = vec[2][i] / sqrt(norm3);
	}
}

void transpose_matrix(double vec[3][3], double ntemp_vec[3][3])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ntemp_vec[i][j] = vec[j][i];
}

static void show_spacegroup_type(const SpglibSpacegroupType spgtype)
{
	FILE* fp = fopen("SYMMCAR", "at+");
	fprintf(fp, "Number: %d\n", spgtype.number);
	//	fprintf(fp,"International: %s\n", spgtype.international_short);
	fprintf(fp, "Long Name: %s\n", spgtype.international_full);
	//	fprintf(fp,"International: %s\n", spgtype.international);
	//	fprintf(fp,"Schoenflies: %s\n", spgtype.schoenflies);
	//	fprintf(fp,"Hall symbol: %s\n", spgtype.hall_symbol);
	//	fprintf(fp,"Point group intl: %s\n", spgtype.pointgroup_international);
	fprintf(fp, "Schoenflies Names: %s\n", spgtype.pointgroup_schoenflies);
	//	fprintf(fp,"Arithmetic cc num: %d\n", spgtype.arithmetic_crystal_class_number);
	//	fprintf(fp,"Arithmetic cc sym: %s\n", spgtype.arithmetic_crystal_class_symbol);
	fclose(fp);
}

static int show_spg_dataset(double lattice[3][3],
	const double origin_shift[3],
	double position[][3],
	const int num_atom,
	const int types[])
{
	SpglibDataset* dataset;
	char ptsymbol[6];
	int pt_trans_mat[3][3];
	int i, j;
	const char* wl = "abcdefghijklmnopqrstuvwxyz";

	for (i = 0; i < num_atom; i++) {
		for (j = 0; j < 3; j++) {
			position[i][j] += origin_shift[j];
		}
	}

	dataset = spg_get_dataset(lattice,
		position,
		types,
		num_atom,
		SYMPREAC);

	if (dataset == NULL) {
		return 1;
	}

	FILE* fp = fopen("SYMMCAR", "at+");
	fprintf(fp, "International Tables: %d\n", dataset->spacegroup_number);
	fprintf(fp, "International: %s \n", dataset->international_symbol);
	fprintf(fp, "Hall symbol: %s\n", dataset->hall_symbol);
	int spacegroup_number = dataset->spacegroup_number;
	/*char symbol[20];
	int space_group = spg_get_international(symbol, lattice, position, types, num_atom, SYMPREAC);*/
	if (spacegroup_number == 1 || spacegroup_number == 2)
		fprintf(fp, "Holohedry: %s\n", "Triclinic");
	if (spacegroup_number >= 3 && spacegroup_number < 16)
		fprintf(fp, "Holohedry: %s\n", "Monoclinic");
	if (spacegroup_number >= 16 && spacegroup_number < 75)
		fprintf(fp, "Holohedry: %s\n", "Orthogonal");
	if (spacegroup_number >= 75 && spacegroup_number < 143)
		fprintf(fp, "Holohedry: %s\n", "Tetragonal");
	if (spacegroup_number >= 143 && spacegroup_number < 168)
		fprintf(fp, "Holohedry: %s\n", "Trigonal");
	if (spacegroup_number >= 168 && spacegroup_number < 195)
		fprintf(fp, "Holohedry: %s\n", "Hexagonal");
	if (spacegroup_number >= 195 && spacegroup_number < 230)
		fprintf(fp, "Holohedry: %s\n", "Cubic");
	if (spg_get_pointgroup(ptsymbol,
		pt_trans_mat,
		dataset->rotations,
		dataset->n_operations))
	{
		fprintf(fp, "Crystal System: %s\n", ptsymbol);
		fclose(fp);
		SpglibSpacegroupType spgtype;
		spgtype = spg_get_spacegroup_type(dataset->hall_number);
		if (spgtype.number) {
			show_spacegroup_type(spgtype);
			spg_free_dataset(dataset);
			return 0;
		}
		else {
			spg_free_dataset(dataset);
			return 1;
		}
	}
	else {
		spg_free_dataset(dataset);
		return 1;
	}
}

static void show_cell(double lattice[3][3],
	double position[][3],
	const int types[],
	const int num_atom)
{
	int i;
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "Lattice parameter:\n");
	for (i = 0; i < 3; i++) {
		fprintf(fp, "%f %f %f\n", lattice[i][0], lattice[i][1], lattice[i][2]);
	}
	fprintf(fp, "Atomic positions:\n");
	for (i = 0; i < num_atom; i++) {
		fprintf(fp, "%d: %f %f %f\n",
			types[i], position[i][0], position[i][1], position[i][2]);
	}
	fclose(fp);
}

void translate_typenum_type(int* type, int atom_type, int typenum[])//typenum->type
{
	int k = 0;
	for (int i = 0; i < atom_type; i++)
		for (int j = 0; j < typenum[i]; j++)
			type[k++] = i + 1;
}

void translate_type_typenum(int* type, int nant[], int typenum[])//type->typenum
{
	for (int i = 0; i < nant[1]; i++)
		typenum[i] = 1;
	for (int i = 0, j = 0; i < nant[0] - 1; i++)
	{
		if (type[i] < 0)
			continue;
		if (type[i] == type[i + 1])
			typenum[j]++;
		else
			j++;
	}
}

int EV_primitive(const char file1[], const char file2[])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return 0;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	double ntemp_vec[3][3];
	transpose_matrix(vec, ntemp_vec);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "*** spg_find_primitive (unitcell --> primitive) ***:\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	int num_primitive_atom;
	num_primitive_atom = spg_find_primitive(ntemp_vec, xyz, types, nant[0], SYMPREAC);
	if (num_primitive_atom)
	{
		fclose(_fp);
		transpose_matrix(ntemp_vec, vec);
		show_cell(vec, xyz, types, num_primitive_atom);
		nant[0] = num_primitive_atom;
		int* primitive_typenum = (int*)malloc(sizeof(int) * nant[1]);
		translate_type_typenum(types, nant, primitive_typenum);
		FILE* fp2 = fopen(file2, "w");
		savposcar(fp2, title, latt, ifix, iflg, nant, primitive_typenum, elemnum, elemsym, vec, xyz, fix);
		fclose(fp2);
		FILE* fp_ = fopen("structure_operator", "at+");
		fprintf(fp_, "****************************************************\n");
		fclose(fp_);
		free(types);
		free(primitive_typenum);
		return 1;
	}
	else
	{
		fprintf(fp, "Primitive cell was not found.\n");
		fprintf(fp, "****************************************************\n");
		fclose(fp);
		free(types);
		return 0;
	}
}

int EV_unitcell(const char file1[], const char file2[])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return 0;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	double ntemp_vec[3][3];
	transpose_matrix(vec, ntemp_vec);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "*** spg_refine_cell ***:\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	/* 4 times larger memory space must be prepared. */
	int* n_types = (int*)malloc(sizeof(int) * nant[0] * 4);
	for (int i = 0; i < nant[0]; i++)
		n_types[i] = types[i];
	double(*n_xyz)[3];
	n_xyz = new double[4 * nant[0]][3];
	for (int i = 0; i < nant[0]; i++)
		for (int j = 0; j < 3; j++)
			n_xyz[i][j] = xyz[i][j];

	int num_atom_bravais;
	num_atom_bravais = spg_refine_cell(ntemp_vec, n_xyz, n_types, nant[0], SYMPREAC);
	if (num_atom_bravais)
	{
		fclose(_fp);
		transpose_matrix(ntemp_vec, vec);
		show_cell(vec, n_xyz, n_types, num_atom_bravais);
		nant[0] = num_atom_bravais;
		int* unit_typenum = (int*)malloc(sizeof(int) * nant[1]);
		for (int i = 0; i < nant[0] - 1; i++)
		{
			for (int j = 0; j < nant[0] - 1 - i; j++)
			{
				if (n_types[j] > n_types[j + 1])
				{
					swap(n_types[j], n_types[j + 1]);
					for (int k = 0; k < 3; k++)
					{
						swap(n_xyz[j][k], n_xyz[j + 1][k]);
					}
				}
			}
		}
		translate_type_typenum(n_types, nant, unit_typenum);
		FILE* fp2 = fopen(file2, "w");
		savposcar(fp2, title, latt, ifix, iflg, nant, unit_typenum, elemnum, elemsym, vec, n_xyz, fix);
		fclose(fp2);
		FILE* fp_ = fopen("structure_operator", "at+");
		fprintf(fp_, "****************************************************\n");
		fclose(fp_);
		return 1;
	}
	else
	{
		fprintf(fp, "Refine cell failed.\n");
		fprintf(fp, "****************************************************\n");
		fclose(fp);
		return 0;
	}
	free(types);
	free(typenum);
	free(n_types);
	for (int i = 0; i < nant[0] * 4; i++)
	{
		free(n_xyz[i]);
	}
	delete[] n_xyz;
}

void EV_get_symmetry(const char file1[])
{
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");

		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	double ntemp_vec[3][3];
	transpose_matrix(vec, ntemp_vec);
	FILE* _fp = fopen("SYMMCAR", "w+");
	fprintf(_fp, "*********************Find Symmetry******************\n");
	fclose(_fp);
	double origin_shift[3] = { 0.1, 0.1, 0 };
	show_spg_dataset(ntemp_vec, origin_shift, xyz, nant[0], types);
	FILE* fp_ = fopen("SYMMCAR", "at+");
	fprintf(fp_, "****************************************************\n");
	fclose(fp_);
	free(types);
}

void EV_supercell(const char file1[], const char file2[], int super[3])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "*********************SuperCell*********************\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	fclose(_fp);
	double ntemp_vec[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ntemp_vec[i][j] = vec[i][j] * super[i];
	double(*n1_xyz)[3];
	n1_xyz = new double[super[0] * nant[0]][3];
	int k1 = 0;
	for (int i = 0; i < super[0] * nant[0]; i++)
	{
		n1_xyz[i][0] = (xyz[i / super[0]][0] + k1) / super[0];
		n1_xyz[i][1] = xyz[i / super[0]][1];
		n1_xyz[i][2] = xyz[i / super[0]][2];
		if (k1 < super[0] - 1)
			k1++;
		else
			k1 = 0;
	}
	nant[0] *= super[0];
	double(*n2_xyz)[3];
	n2_xyz = new double[nant[0] * super[1]][3];
	int k2 = 0;
	for (int i = 0; i < super[1] * nant[0]; i++)
	{
		n2_xyz[i][1] = (n1_xyz[i / super[1]][1] + k2) / super[1];
		n2_xyz[i][0] = n1_xyz[i / super[1]][0];
		n2_xyz[i][2] = n1_xyz[i / super[1]][2];
		if (k2 < super[1] - 1)
			k2++;
		else
			k2 = 0;
	}
	nant[0] *= super[1];
	double(*n3_xyz)[3];
	n3_xyz = new double[nant[0] * super[2]][3];
	int k3 = 0;
	for (int i = 0; i < super[2] * nant[0]; i++)
	{
		n3_xyz[i][2] = (n2_xyz[i / super[2]][2] + k3) / super[2];
		n3_xyz[i][1] = n2_xyz[i / super[2]][1];
		n3_xyz[i][0] = n2_xyz[i / super[2]][0];
		if (k3 < super[2] - 1)
			k3++;
		else
			k3 = 0;
	}
	nant[0] *= super[2];
	int super_num = super[0] * super[1] * super[2];
	for (int i = 0; i < nant[1]; i++)
		typenum[i] *= super_num;
	int* n_types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(n_types, nant[1], typenum);
	show_cell(ntemp_vec, n3_xyz, n_types, nant[0]);
	FILE* fp2 = fopen(file2, "w");
	savposcar(fp2, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, ntemp_vec, n3_xyz, fix);
	fclose(fp2);
	FILE* fp_ = fopen("structure_operator", "at+");
	fprintf(fp_, "****************************************************\n");
	fclose(fp_);
	free(types);
	delete[] n1_xyz;
	delete[] n2_xyz;
	delete[] n3_xyz;
	free(n_types);
}

void get_defvect(char mode[], double current_strain, double defvect[3][3])
{
	//tensile strain along x axis
	if (!strcmp(mode, "xx"))
	{
		defvect[0][0] = 1.0000 + current_strain;
	}
	//tensile strain along y axis
	if (!strcmp(mode, "yy"))
	{
		defvect[1][1] = 1.0000 + current_strain;
	}
	//tensile strain along z axis
	if (!strcmp(mode, "zz"))
	{
		defvect[2][2] = 1.0000 + current_strain;
	}
	// shear strain in x-y plane
	if (!strcmp(mode, "xy"))
	{
		defvect[0][1] = current_strain / 2.0000;
		defvect[1][0] = current_strain / 2.0000;
	}
	//shear strain in x-z plane
	if (!strcmp(mode, "xz"))
	{
		defvect[0][2] = current_strain / 2.0000;
		defvect[2][0] = current_strain / 2.0000;
	}
	//shear strain in y-z plane
	if (!strcmp(mode, "yz"))
	{
		defvect[1][2] = current_strain / 2.0000;
		defvect[2][1] = current_strain / 2.0000;
	}
}

void EV_affine(const char file1[], char model[], char mode[], double init_value, double step_length, int step_num)
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	if (!strcmp(model, "-tension") && (!strcmp(mode, "xy") || !strcmp(mode, "yz") ||
		!strcmp(mode, "xz") || !strcmp(mode, "yx") || !strcmp(mode, "zy") || !strcmp(mode, "zx")))
	{
		fprintf(fp, "%s and %s do not match!\n", model, mode);
		return;
	}
	if ((!strcmp(model, "-simshear") || !strcmp(model, "-purshear")) && (!strcmp(mode, "xx") || !strcmp(mode, "yy") || !strcmp(mode, "zz")))
	{
		fprintf(fp, "%s and %s do not match!\n", model, mode);
		return;
	}
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "***********************AFFINE*********************\n");
	fclose(_fp);
	char(*file_name)[50];
	file_name = new char[step_num + 1][50];
	for (int i = 0; i < step_num + 1; i++)
	{
		FILE* fp = fopen("structure_operator", "at+");
		double xyz_e[MAX_NATOM][3] = { 0 };
		double vec_e[3][3] = { 0 };
		double defvect[3][3] = { {1.0000, 0.0000, 0.0000},{0.0000, 1.0000, 0.0000},{0.0000, 0.0000, 1.0000} };
		double current_strain = init_value + step_length * i;
		get_defvect(mode, current_strain, defvect);
		if (!strcmp("-simshear", model))
		{
			if (!strcmp(mode, "xy"))
			{
				defvect[1][0] = 0;
				defvect[0][1] *= 2;
			}
			if (!strcmp(mode, "yx"))
			{
				defvect[0][1] = 0;
				defvect[1][0] *= 2;
			}
			if (!strcmp(mode, "xz"))
			{
				defvect[2][0] = 0;
				defvect[0][2] *= 2;
			}
			if (!strcmp(mode, "zx"))
			{
				defvect[0][2] = 0;
				defvect[2][0] *= 2;
			}
			if (!strcmp(mode, "yz"))
			{
				defvect[2][1] = 0;
				defvect[1][2] *= 2;
			}
			if (!strcmp(mode, "zy"))
			{
				defvect[1][2] = 0;
				defvect[2][1] *= 2;
			}
		}
		sprintf(file_name[i], "AFFPOS_%.4lf", current_strain);
		fprintf(fp, "FILE NAME: %s\n", file_name[i]);
		fclose(fp);
		FILE* fp_aff = fopen(file_name[i], "w");
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				vec_e[i][j] = vec[i][0] * defvect[0][j] + vec[i][1] * defvect[1][j] + vec[i][2] * defvect[2][j];
		show_cell(vec_e, xyz, types, nant[0]);
		savposcar(fp_aff, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz, fix);
		fclose(fp_aff);
		FILE* fp_ = fopen("structure_operator", "at+");
		fprintf(fp_, "****************************************************\n");
		fclose(fp_);
	}
	free(types);
	delete[] file_name;
}

void EV_alias_tensile(const char file[], char mode[], double istart, double iend, double ispacing, double value)
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "FILE NAME: %s\n", file);
	fclose(fp);
	FILE* fp1 = fopen(file, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	double dx[3], dy[3], dz[3];
	get_unitvector(vec, dx, dy, dz);
	if (iflg == 0)
		direct_to_carts(iflg, vec, xyz, nant);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "***********************ALIAS*********************\n");
	fclose(_fp);
	int ax;
	double unitvec[3];
	if (!strcmp(mode, "xx"))
	{
		ax = 0;
		unitvec[0] = dx[0];
		unitvec[1] = dx[1];
		unitvec[2] = dx[2];
	}
	else if (!strcmp(mode, "yy"))
	{
		ax = 1;
		unitvec[0] = dy[0];
		unitvec[1] = dy[1];
		unitvec[2] = dz[2];
	}
	else if (!strcmp(mode, "zz"))
	{
		ax = 2;
		unitvec[0] = dz[0];
		unitvec[1] = dz[1];
		unitvec[2] = dz[2];
	}
	else
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is a wrong parameter!\n", mode);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	int step_num;
	if (ispacing == 0)
		step_num = 1;
	else
		step_num = (iend - istart) / ispacing + 1;
	char(*file_name)[50];
	file_name = new char[step_num][50];
	/*ifix = 0;
	for (int i = 0; i < nant[0]; i++)
	{
		fix[i][0] = 'F';
		fix[i][1] = 'F';
		fix[i][2] = 'T';
	}*/
	for (int i = 0; i < step_num; i++)
	{
		double xyz_e[MAX_NATOM][3] = { 0 };
		double vec_e[3][3] = { 0 };
		FILE* fp = fopen("structure_operator", "at+");
		double current = istart + ispacing * i;
		sprintf(file_name[i], "ALIPOS_%.4lf", current);
		fprintf(fp, "FILE NAME: %s\n", file_name[i]);
		fclose(fp);
		FILE* fp_ali = fopen(file_name[i], "w");
		for (int j = 0; j < nant[0]; j++)
		{
			if (xyz[j][ax] > value)
			{
				xyz_e[j][0] = xyz[j][0] + current * unitvec[0];
				xyz_e[j][1] = xyz[j][1] + current * unitvec[1];
				xyz_e[j][2] = xyz[j][2] + current * unitvec[2];
			}
			else
			{
				xyz_e[j][0] = xyz[j][0];
				xyz_e[j][1] = xyz[j][1];
				xyz_e[j][2] = xyz[j][2];
			}
		}
		for (int m = 0; m < 3; m++)
			for (int n = 0; n < 3; n++)
				vec_e[m][n] = vec[m][n];
		vec_e[ax][0] += current * unitvec[0];
		vec_e[ax][1] += current * unitvec[1];
		vec_e[ax][2] += current * unitvec[2];
		show_cell(vec_e, xyz_e, types, nant[0]);
		if (iflg != 1) iflg = 1;
		carts_to_direct(iflg, vec_e, xyz_e, nant);
		savposcar(fp_ali, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz_e, fix);
		fclose(fp_ali);
		FILE* fp_ = fopen("structure_operator", "at+");
		fprintf(fp_, "****************************************************\n");
		fclose(fp_);
	}
	free(types);
	delete[] file_name;
}

void EV_alias_shear(const char file[], char mode[], double istart1, double iend1, double ispacing1,
	double istart2, double iend2, double ispacing2, double value)
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "FILE NAME: %s\n", file);
	fclose(fp);
	FILE* fp1 = fopen(file, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	double dx[3], dy[3], dz[3];
	get_unitvector(vec, dx, dy, dz);
	if (iflg == 0)
		direct_to_carts(iflg, vec, xyz, nant);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "***********************ALIAS*********************\n");
	fclose(_fp);
	int step_num1, step_num2;
	if (ispacing1 == 0)
		step_num1 = 0;
	else
		step_num1 = (iend1 - istart1) / ispacing1 + 1;
	if (ispacing2 == 0)
		step_num2 = 0;
	else
		step_num2 = (iend2 - istart2) / ispacing2 + 1;
	char(*file_name)[50];
	file_name = new char[max(step_num1, step_num2)][50];
	/*ifix = 0;
	for (int i = 0; i < nant[0]; i++)
	{
		fix[i][0] = 'F';
		fix[i][1] = 'F';
		fix[i][2] = 'T';
	}*/
	for (int i = 0; i < max(step_num1, step_num2); i++)
	{
		double xyz_e[MAX_NATOM][3] = { 0 };
		double vec_e[3][3] = { 0 };
		FILE* fp = fopen("structure_operator", "at+");
		double current1 = min(iend1, istart1 + ispacing1 * i);
		double current2 = min(iend2, istart2 + ispacing2 * i);
		sprintf(file_name[i], "ALIPOS_%.4lf_%.4lf", current1, current2);
		fprintf(fp, "FILE NAME: %s\n", file_name[i]);
		fclose(fp);
		FILE* fp_ali = fopen(file_name[i], "w");
		if (!strcmp("xy", mode))
		{
			for (int j = 0; j < nant[0]; j++)
			{
				double distance = 0;
				double ab_cross[3];
				EV_cross(dx, dy, ab_cross);
				distance = fabs(xyz[j][0] * ab_cross[0] + xyz[j][1] * ab_cross[1] + xyz[j][2] * ab_cross[2]);
				if (distance > value)
				{
					xyz_e[j][0] = current1 * dx[0] + current2 * dy[0] + xyz[j][0];
					xyz_e[j][1] = current1 * dx[1] + current2 * dy[1] + xyz[j][1];
					xyz_e[j][2] = current1 * dx[2] + current2 * dy[2] + xyz[j][2];
				}
				else
				{
					xyz_e[j][0] = xyz[j][0];
					xyz_e[j][1] = xyz[j][1];
					xyz_e[j][2] = xyz[j][2];
				}
			}
			for (int m = 0; m < 3; m++)
				for (int n = 0; n < 3; n++)
					vec_e[m][n] = vec[m][n];
			vec_e[2][0] += current1 * dx[0] + current2 * dy[0];
			vec_e[2][1] += current1 * dx[1] + current2 * dy[1];
			vec_e[2][2] += current1 * dx[2] + current2 * dy[2];
			show_cell(vec_e, xyz_e, types, nant[0]);
			if (iflg != 1) iflg = 1;
			carts_to_direct(iflg, vec_e, xyz_e, nant);
			savposcar(fp_ali, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz_e, fix);
			fclose(fp_ali);
			FILE* fp_ = fopen("structure_operator", "at+");
			fprintf(fp_, "****************************************************\n");
			fclose(fp_);
		}
		if (!strcmp("xz", mode))
		{
			for (int j = 0; j < nant[0]; j++)
			{
				int distance = 0;
				double ab_cross[3];
				EV_cross(dx, dz, ab_cross);
				distance = fabs(xyz[j][0] * ab_cross[0] + xyz[j][1] * ab_cross[1] + xyz[j][2] * ab_cross[3]);
				if (distance > value)
				{
					xyz_e[j][0] = current1 * dx[0] + current2 * dz[0] + xyz[j][0];
					xyz_e[j][1] = current1 * dx[1] + current2 * dz[1] + xyz[j][1];
					xyz_e[j][2] = current1 * dx[2] + current2 * dz[2] + xyz[j][2];
				}
				else
				{
					xyz_e[j][0] = xyz[j][0];
					xyz_e[j][1] = xyz[j][1];
					xyz_e[j][2] = xyz[j][2];
				}
			}
			for (int m = 0; m < 3; m++)
				for (int n = 0; n < 3; n++)
					vec_e[m][n] = vec[m][n];
			vec_e[1][0] += current1 * dx[0] + current2 * dz[0];
			vec_e[1][1] += current1 * dx[1] + current2 * dz[1];
			vec_e[1][2] += current1 * dx[2] + current2 * dz[2];
			show_cell(vec_e, xyz_e, types, nant[0]);
			carts_to_direct(iflg, vec_e, xyz_e, nant);
			if (iflg != 1) iflg = 1;
			savposcar(fp_ali, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz_e, fix);
			fclose(fp_ali);
			FILE* fp_ = fopen("structure_operator", "at+");
			fprintf(fp_, "****************************************************\n");
			fclose(fp_);
		}
		if (!strcmp("yz", mode))
		{
			for (int j = 0; j < nant[0]; j++)
			{
				int distance = 0;
				double ab_cross[3];
				EV_cross(dy, dz, ab_cross);
				distance = fabs(xyz[j][0] * ab_cross[0] + xyz[j][1] * ab_cross[1] + xyz[j][2] * ab_cross[3]);
				if (distance > value)
				{
					xyz_e[j][0] = current1 * dy[0] + current2 * dz[0] + xyz[j][0];
					xyz_e[j][1] = current1 * dy[1] + current2 * dz[1] + xyz[j][1];
					xyz_e[j][2] = current1 * dy[2] + current2 * dz[2] + xyz[j][2];
				}
				else
				{
					xyz_e[j][0] = xyz[j][0];
					xyz_e[j][1] = xyz[j][1];
					xyz_e[j][2] = xyz[j][2];
				}
			}
			for (int m = 0; m < 3; m++)
				for (int n = 0; n < 3; n++)
					vec_e[m][n] = vec[m][n];
			vec_e[0][0] += current1 * dy[0] + current2 * dz[0];
			vec_e[0][1] += current1 * dy[1] + current2 * dz[1];
			vec_e[0][2] += current1 * dy[2] + current2 * dz[2];
			show_cell(vec_e, xyz_e, types, nant[0]);
			carts_to_direct(iflg, vec_e, xyz_e, nant);
			if (iflg != 1) iflg = 1;
			savposcar(fp_ali, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz_e, fix);
			fclose(fp_ali);
			FILE* fp_ = fopen("structure_operator", "at+");
			fprintf(fp_, "****************************************************\n");
			fclose(fp_);
		}
	}
	free(types);
	delete[] file_name;
}

void EV_rotproj(const char file1[], const char file2[], double rot[3])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "*********************PROJECTION*********************\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	fclose(_fp);
	for (int i = 0; i < 3; i++)
		rot[i] = rot[i] / 180 * Pi;
	//The rotation matrix rotxmat
	double rotxmat1[3][3] = { 1,0,0,0,cos(rot[0]),-sin(rot[0]),0,sin(rot[0]),cos(rot[0]) };
	double rotxmat[3][3];
	transpose_matrix(rotxmat1, rotxmat);
	//xx axis rotation
	double rotxvect[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rotxvect[i][j] = 0;
			for (int k = 0; k < 3; k++)
				rotxvect[i][j] += vec[i][k] * rotxmat[k][j];
		}
	}
	//The rotation matrix rotymat
	double rotymat1[3][3] = { cos(rot[1]),0,sin(rot[1]),0,1,0,-sin(rot[1]),0,cos(rot[1]) };
	double rotymat[3][3];
	transpose_matrix(rotymat1, rotymat);
	//yy axis rotation
	double rotyvect[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rotyvect[i][j] = 0;
			for (int k = 0; k < 3; k++)
				rotyvect[i][j] += rotxvect[i][k] * rotymat[k][j];
		}
	}
	//The rotation matrix rotzmat
	double rotzmat1[3][3] = { cos(rot[2]),-sin(rot[2]),0,sin(rot[2]),cos(rot[2]),0,0,0,1 };
	double rotzmat[3][3];
	transpose_matrix(rotzmat1, rotzmat);
	//zz axis rotation
	double rotzvect[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rotzvect[i][j] = 0;
			for (int k = 0; k < 3; k++)
				rotzvect[i][j] += rotyvect[i][k] * rotzmat[k][j];
		}
	}
	double vec_e[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			vec_e[i][j] = rotzvect[i][j];
	show_cell(vec_e, xyz, types, nant[0]);
	FILE* fp_proj = fopen(file2, "w");
	savposcar(fp_proj, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz, fix);
	fclose(fp_proj);
	FILE* fp_ = fopen("structure_operator", "at+");
	fprintf(fp_, "****************************************************\n");
	fclose(fp_);
	free(types);
}

void EV_indproj(const char file1[], const char file2[], int pvh, int pvk, int pvl, int uvu, int uvv, int uvw)
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "*********************PROJECTION*********************\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	fclose(_fp);
	// Define the oldx and oldy
	double oldx[3], oldy[3], oldz[3];
	oldx[0] = pvh; oldx[1] = pvk; oldx[2] = pvl;
	oldy[0] = uvu; oldy[1] = uvv; oldy[2] = uvw;
	oldz[0] = 0.0; oldz[1] = 0.0; oldz[2] = 0.0;
	if (oldx[0] * oldy[0] + oldx[1] * oldy[1] + oldx[2] * oldy[2] != 0)
	{
		FILE* fp_error = fopen("structure_operator", "at");
		fprintf(fp_error, "------------------ERROR!----------------\n");
		fprintf(fp_error, "Upward vector is not in projected plane!\n");
		fprintf(fp_error, "u, v, w, h, k and l must satisfy the condition : \n");
		fprintf(fp_error, "hu + kv + lw = 0\n");
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	//define the reciprocal vector of privect
	double tmpvect[5][5];
	double recipvect[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			tmpvect[i][j] = vec[i][j];
	for (int i = 0; i < 3; i++)
		for (int j = 3; j < 5; j++)
			tmpvect[i][j] = tmpvect[i][j - 3];
	for (int i = 3; i < 5; i++)
		for (int j = 0; j < 5; j++)
			tmpvect[i][j] = tmpvect[i - 3][j];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			recipvect[i][j] = tmpvect[i + 1][j + 1] * tmpvect[i + 2][j + 2] - tmpvect[i + 1][j + 2] * tmpvect[i + 2][j + 1];
	// define the newY and newX relative of XYZ
	double newx[5], newy[5], newz[5];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			newx[i] += oldx[j] * recipvect[j][i];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			newy[i] += oldy[j] * vec[j][i];
	// define the newz relative of xyz
	for (int i = 3; i < 5; i++)
	{
		newx[i] = newx[i - 3];
		newy[i] = newy[i - 3];
	}
	for (int i = 0; i < 3; i++)
		newz[i] = newx[i + 1] * newy[i + 2] - newy[i + 1] * newx[i + 2];
	// define the transformation matrix
	double projmat[3][3];
	for (int i = 0; i < 3; i++)
	{
		projmat[0][i] = newx[i] / sqrt(newx[0] * newx[0] + newx[1] * newx[1] + newx[2] * newx[2]);
		projmat[1][i] = newy[i] / sqrt(newy[0] * newy[0] + newy[1] * newy[1] + newy[2] * newy[2]);
		projmat[2][i] = newz[i] / sqrt(newz[0] * newz[0] + newz[1] * newz[1] + newz[2] * newz[2]);
	}
	double projmatI[3][3];
	transpose_matrix(projmat, projmatI);
	double vec_e[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				vec_e[i][j] += vec[i][k] * projmatI[k][j];
	show_cell(vec_e, xyz, types, nant[0]);
	FILE* fp_proj = fopen(file2, "w");
	savposcar(fp_proj, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz, fix);
	fclose(fp_proj);
	FILE* fp_ = fopen("structure_operator", "at+");
	fprintf(fp_, "****************************************************\n");
	fclose(fp_);
	free(types);
}

void EV_matproj(const char file1[], const char file2[], double projmat[3][3])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "*********************PROJECTION*********************\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	fclose(_fp);
	double projmatI[3][3];
	transpose_matrix(projmat, projmatI);
	double vec_e[3][3];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			vec_e[i][j] = 0;
			for (int k = 0; k < 3; k++)
				vec_e[i][j] += vec[i][k] * projmatI[k][j];
		}
	}
	show_cell(vec_e, xyz, types, nant[0]);
	FILE* fp_proj = fopen(file2, "w");
	savposcar(fp_proj, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz, fix);
	fclose(fp_proj);
	FILE* fp_ = fopen("structure_operator", "at+");
	fprintf(fp_, "****************************************************\n");
	fclose(fp_);
	free(types);
}

void EV_redefine(const char file1[], const char file2[], int rot[3][3], double symprec, int a1, int x1, int a2, int x2,bool fixdirect)
{
	FILE* fptran = fopen("redef.in", "w");
	fprintf(fptran, "VASPMATE will redefine cell from the following matrix.\n");
	for (int i = 0; i < 3; i++)
		fprintf(fptran, "  %d  %d  %d   # must be three integers\n", rot[i][0], rot[i][1], rot[i][2]);
	fclose(fptran);
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	if ((rot[0][0] * rot[1][1] * rot[2][2] + rot[0][1] * rot[1][2] * rot[2][0] + rot[0][2] * rot[1][0] * rot[2][1] -
		rot[0][0] * rot[1][2] * rot[2][1] - rot[0][2] * rot[1][1] * rot[2][0] - rot[0][1] * rot[1][0] * rot[2][2]) == 0)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR! The determinant of the matrix is 0!\n");
		fclose(fp_error);
	}
	int super[3] = { abs(rot[0][0]) + abs(rot[1][0]) + abs(rot[2][0]),
		abs(rot[0][1]) + abs(rot[1][1]) + abs(rot[2][1]),
		abs(rot[0][2]) + abs(rot[1][2]) + abs(rot[2][2]) };
	EV_supercell(file1, "temp", super);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "**********************REDEFINE*********************\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	fclose(_fp);
	char   title_e[MAX_NCHAR];        // Title of system
	double latt_e;                    // Scaling factor
	int    ifix_e;                    // Select 0 or normal 1
	int    iflg_e;                    // Direct 0 or Cartes 1
	int    nant_e[2];                 // number of atoms and types
	int    typenum_e[MAX_NELEM];      // number of atoms of each type
	int    elemnum_e[MAX_NELEM];      // atomic number
	char   elemsym_e[MAX_NELEM][3];   // atomic symbol
	double vec_e[3][3];               // lattice vector
	double xyz_e[MAX_NATOM][3];       // atomic coordinate
	char   fix_e[MAX_NATOM][3];       // Selective fix on each atom
	FILE* fp_super = fopen("temp", "r");
	readposcar(fp_super, title_e, latt_e, ifix_e, iflg_e, nant_e, typenum_e, elemnum_e, elemsym_e, vec_e, xyz_e, fix_e);
	fclose(fp_super);
	int length = nant_e[0];
	direct_to_carts(iflg_e, vec_e, xyz_e, nant_e);
	int* types_e = (int*)malloc(sizeof(int) * nant_e[0]);
	translate_typenum_type(types_e, nant_e[1], typenum_e);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			vec_e[i][j] = rot[i][0] * vec[0][j] + rot[i][1] * vec[1][j] + rot[i][2] * vec[2][j];
	carts_to_direct(iflg_e, vec_e, xyz_e, nant_e);
	double dis[3];
	for (int i = 0; i < nant_e[0] - 1; i++)
	{
		for (int j = i + 1; j < nant_e[0]; j++)
		{
			dis[0] = xyz_e[i][0] - xyz_e[j][0];
			dis[1] = xyz_e[i][1] - xyz_e[j][1];
			dis[2] = xyz_e[i][2] - xyz_e[j][2];
			for (int k = 0; k < 3; k++)
			{
				while (dis[k] < -0.5)dis[k] += 1;
				while (dis[k] > 0.5)dis[k] -= 1;
			}
			if (dis[0] * dis[0] + dis[1] * dis[1] + dis[2] * dis[2] < symprec * symprec)
			{
				types_e[i] = -1;
				break;
			}
		}
	}
	for (int i = 0; i < length; i++)
	{
		if (types_e[i] < 0)
			nant_e[0]--;
	}
	double(*_xyz)[3];
	_xyz = new double[nant_e[0]][3];
	for (int i = 0, j = 0; i < nant_e[0]; j++)
	{
		if (types_e[j] < 0)
			continue;
		_xyz[i][0] = xyz_e[j][0];
		_xyz[i][1] = xyz_e[j][1];
		_xyz[i][2] = xyz_e[j][2];
		i++;
	}
	for (int i = 0; i < nant_e[0]; i++)
		for (int j = 0; j < 3; j++)
		{
			while (_xyz[i][j] >= 0.9999)
				_xyz[i][j]--;
			while (_xyz[i][j] < 0)
				_xyz[i][j]++;
		}
	std::sort(types_e, types_e + length);
	int temp_nant[] = { length,nant_e[1] };
	translate_type_typenum(types_e, temp_nant, typenum_e);
	int* n_types_e = (int*)malloc(sizeof(int) * nant_e[0]);
	translate_typenum_type(n_types_e, nant_e[1], typenum_e);
	double _vec[3][3];
	//reshape box
	if (fixdirect)
	{
		int i, j;
		int a3, x3;
		double temp[3];
		double a, b, c, ab, bb, ac, bc, cc;
		double vb[3], vc[3];
		double temp_vec[3][3];
		transpose_matrix(vec_e, temp_vec);
		a = sqrt(temp_vec[0][0] * temp_vec[0][0] + temp_vec[1][0] * temp_vec[1][0] + temp_vec[2][0] * temp_vec[2][0]);
		b = sqrt(temp_vec[0][1] * temp_vec[0][1] + temp_vec[1][1] * temp_vec[1][1] + temp_vec[2][1] * temp_vec[2][1]);
		ab = (temp_vec[0][0] * temp_vec[0][1] + temp_vec[1][0] * temp_vec[1][1] + temp_vec[2][0] * temp_vec[2][1]) / a;
		for (i = 0; i < 3; i++)
			vb[i] = temp_vec[i][1] - temp_vec[i][0] * ab / a;
		bb = sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);

		ac = (temp_vec[0][0] * temp_vec[0][2] + temp_vec[1][0] * temp_vec[1][2] + temp_vec[2][0] * temp_vec[2][2]) / a;
		bc = (vb[0] * temp_vec[0][2] + vb[1] * temp_vec[1][2] + vb[2] * temp_vec[2][2]) / bb;

		for (i = 0; i < 3; i++)
			vc[i] = temp_vec[i][2] - temp_vec[i][0] * ac / a - vb[i] * bc / bb;

		if (temp_vec[0][2] * (temp_vec[1][0] * temp_vec[2][1] - temp_vec[2][0] * temp_vec[1][1]) +
			temp_vec[1][2] * (temp_vec[2][0] * temp_vec[0][1] - temp_vec[0][0] * temp_vec[2][1]) +
			temp_vec[2][2] * (temp_vec[0][0] * temp_vec[1][1] - temp_vec[1][0] * temp_vec[0][1]) > 0)
			cc = sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);
		else
			cc = -sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);
		for (i = 0; i < 3; i++)
		{
			if (i != a1 && i != a2)
				a3 = i;
			if (i != x1 && i != x2)
				x3 = i;
		}

		temp_vec[x1][a1] = a;
		temp_vec[x1][a2] = ab;
		temp_vec[x1][a3] = ac;
		temp_vec[x2][a1] = 0;
		temp_vec[x2][a2] = bb;
		temp_vec[x2][a3] = bc;
		temp_vec[x3][a1] = 0;
		temp_vec[x3][a2] = 0;
		temp_vec[x3][a3] = cc;

		for (i = 0; i < nant_e[0]; i++)
		{
			temp[0] = _xyz[i][0];
			temp[1] = _xyz[i][1];
			temp[2] = _xyz[i][2];
			_xyz[i][0] = temp[a1];
			_xyz[i][1] = temp[a2];
			_xyz[i][2] = temp[a3];
		}
		transpose_matrix(temp_vec, _vec);
	}
	else
	{
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				_vec[i][j] = vec_e[i][j];
	}
	show_cell(_vec, _xyz, n_types_e, nant_e[0]);
	FILE* fp_redef = fopen(file2, "w");
	savposcar(fp_redef, title_e, latt_e, ifix_e, iflg_e, nant_e, typenum_e, elemnum_e, elemsym_e, _vec, _xyz, fix_e);
	fclose(fp_redef);
	FILE* fp_ = fopen("structure_operator", "at+");
	fprintf(fp_, "****************************************************\n");
	fclose(fp_);
	free(types);
	free(types_e);
	free(n_types_e);
	remove("temp");
	delete[] _xyz;
	printf("Written redef.in file!\n");
	printf("Written REDPOS file!\n");
}

void EV_recell(const char file1[], const char file2[])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	show_cell(vec, xyz, types, nant[0]);
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "**********************RECELL***********************\n");
	fprintf(_fp, "FILE NAME: %s\n", file2);
	fclose(_fp);
	//get spacegroupnum
	char symbol[20];
	double ntemp_vec[3][3];
	transpose_matrix(vec, ntemp_vec);
	int spacegroup_number = spg_get_international(symbol, ntemp_vec, xyz, types, nant[0], SYMPREAC);
	if (spacegroup_number == 0)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp, "ERROR!Cannot get space group!\n");
		return;
	}
	//refine cell
	double n_vec[3][3];
	transpose_matrix(vec, n_vec);
	/* 4 times larger memory space must be prepared. */
	int* n_types = (int*)malloc(sizeof(int) * nant[0] * 4);
	for (int i = 0; i < nant[0]; i++)
		n_types[i] = types[i];
	double(*n_xyz)[3];
	n_xyz = new double[4 * nant[0]][3];
	for (int i = 0; i < nant[0]; i++)
		for (int j = 0; j < 3; j++)
			n_xyz[i][j] = xyz[i][j];

	int refine_cell = spg_refine_cell(n_vec, n_xyz, n_types, nant[0], SYMPREAC);
	transpose_matrix(n_vec, vec);
	if (refine_cell == 0)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp, "ERROR!Cannot recell the structure!\n");
		return;
	}
	nant[0] = refine_cell;
	int* unit_typenum = (int*)malloc(sizeof(int) * nant[1]);
	for (int i = 0; i < nant[0] - 1; i++)
	{
		for (int j = 0; j < nant[0] - 1 - i; j++)
		{
			if (n_types[j] > n_types[j + 1])
			{
				swap(n_types[j], n_types[j + 1]);
				for (int k = 0; k < 3; k++)
				{
					swap(n_xyz[j][k], n_xyz[j + 1][k]);
				}
			}
		}
	}
	translate_type_typenum(n_types, nant, unit_typenum);

	//_vec, n_xyz, n_types, unit_typenu
	double temp, lega, legb, legc, legar, legbr, legcr, cos_alpha, cos_bata, cos_gamma, cos_alphar,
		cos_batar, cos_gammar, tempvect[5][5], recipvect[3][3], cos_angle, temppos[MAX_NATOM][3];

	//Triclinic lattice : c < a < b
	if (spacegroup_number == 1 || spacegroup_number == 2)
	{
		lega = sqrt(vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]);
		legb = sqrt(vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]);
		legc = sqrt(vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]);
		cos_alpha = (vec[1][0] * vec[2][0] + vec[1][1] * vec[2][1] + vec[1][2] * vec[2][2]) / legb / legc;
		cos_bata = (vec[0][0] * vec[2][0] + vec[0][1] * vec[2][0] + vec[0][2] * vec[2][2]) / lega / legc;
		cos_gamma = (vec[0][0] * vec[1][0] + vec[0][1] * vec[1][1] + vec[0][2] * vec[1][2]) / lega / legb;
		//define the reciprocal matrix of privect
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				tempvect[i][j] = vec[i][j];
		for (int i = 0; i < 3; i++)
			for (int j = 3; j < 5; j++)
				tempvect[i][j] = tempvect[i][j - 3];
		for (int i = 3; i < 5; i++)
			for (int j = 0; j < 5; j++)
				tempvect[i][j] = tempvect[i - 3][j];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				recipvect[i][j] = tempvect[i + 1][j + 1] * tempvect[i + 2][j + 2] - tempvect[i + 1][j + 2] * tempvect[i + 2][j + 1];
		//end of defining the reciprocal matrix of privect
		legar = sqrt(recipvect[0][0] * recipvect[0][0] + recipvect[0][1] * recipvect[0][1] + recipvect[0][2] * recipvect[0][2]);
		legbr = sqrt(recipvect[1][0] * recipvect[1][0] + recipvect[1][1] * recipvect[1][1] + recipvect[1][2] * recipvect[1][2]);
		legcr = sqrt(recipvect[2][0] * recipvect[2][0] + recipvect[2][1] * recipvect[2][1] + recipvect[2][2] * recipvect[2][2]);
		cos_alphar = (vec[0][0] * recipvect[0][0] + vec[0][1] * recipvect[0][1] + vec[0][2] * recipvect[0][2]) / lega / legar;
		cos_batar = (vec[1][0] * recipvect[1][0] + vec[1][1] * recipvect[1][1] + vec[1][2] * recipvect[1][2]) / legb / legbr;
		cos_gammar = (vec[2][0] * recipvect[2][0] + vec[2][1] * recipvect[2][1] + vec[2][2] * recipvect[2][2]) / legc / legcr;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				tempvect[i][j] = 0.0000;
		// a <= b <= c
		if (lega <= legc && legb <= legc && lega <= legb)
		{
			tempvect[2][2] = lega;
			tempvect[0][2] = legb * cos_gamma;
			tempvect[0][0] = sqrt(legb * legb - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legc * cos_bata;
			tempvect[1][1] = legc * cos_gammar;
			tempvect[1][0] = sqrt(legc * legc - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			//To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / legb / legc;
			if (fabs(cos_angle - cos_alpha) >= ZERO_TOLERANCE)
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][0];
				temppos[j][0] = n_xyz[j][1];
				temppos[j][1] = n_xyz[j][2];
			}
		}
		// a <= c <= b
		if (lega <= legb && legc <= legb && lega <= legc)
		{
			tempvect[2][2] = lega;
			tempvect[0][2] = legc * cos_bata;
			tempvect[0][0] = sqrt(legc * legc - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legb * cos_gamma;
			tempvect[1][1] = legb * cos_batar;
			tempvect[1][0] = sqrt(legb * legb - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			//To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / legb / legc;
			if (fabs(cos_angle - cos_bata) >= ZERO_TOLERANCE)
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][0];
				temppos[j][0] = n_xyz[j][2];
				temppos[j][1] = n_xyz[j][1];
			}
		}
		// b <= a <= c
		if (legb <= legc && lega <= legc && legb <= lega)
		{
			tempvect[2][2] = legb;
			tempvect[0][2] = lega * cos_gamma;
			tempvect[0][0] = sqrt(lega * lega - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legc * cos_alpha;
			tempvect[1][1] = legb * cos_gamma;
			tempvect[1][0] = sqrt(legc * legc - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legc;
			if (fabs(cos_angle - cos_bata >= ZERO_TOLERANCE))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][1];
				temppos[j][0] = n_xyz[j][0];
				temppos[j][1] = n_xyz[j][2];
			}
		}
		// b <= c <= a
		if (legb <= lega && legc <= lega && legb <= legc)
		{
			tempvect[2][2] = legb;
			tempvect[0][2] = legc * cos_alpha;
			tempvect[0][0] = sqrt(legc * legc - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = lega * cos_gamma;
			tempvect[1][1] = legb * cos_alpha;
			tempvect[1][0] = sqrt(lega * lega - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legc;
			if (fabs(cos_angle - cos_bata >= ZERO_TOLERANCE))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][1];
				temppos[j][0] = n_xyz[j][2];
				temppos[j][1] = n_xyz[j][0];
			}
		}
		// c <= a <= b
		if (legc <= legb && lega <= legb && legc <= lega)
		{
			tempvect[2][2] = legc;
			tempvect[0][2] = lega * cos_bata;
			tempvect[0][0] = sqrt(lega * lega - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legb * cos_alpha;
			tempvect[1][1] = legb * cos_batar;
			tempvect[1][0] = sqrt(legb * legb - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legb;
			if (fabs(cos_angle - cos_gamma >= ZERO_TOLERANCE))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][2];
				temppos[j][0] = n_xyz[j][0];
				temppos[j][1] = n_xyz[j][1];
			}
		}
		// c <= b <= a
		if (legc <= lega && legb <= lega && legc <= legb)
		{
			tempvect[2][2] = legc;
			tempvect[0][2] = legb * cos_alpha;
			tempvect[0][0] = sqrt(legb * legb - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = lega * cos_bata;
			tempvect[1][1] = lega * cos_alphar;
			tempvect[1][0] = sqrt(lega * lega - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legb;
			if (fabs(cos_angle - cos_gamma >= ZERO_TOLERANCE))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][2];
				temppos[j][0] = n_xyz[j][1];
				temppos[j][1] = n_xyz[j][0];
			}
		}
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				vec[i][j] = tempvect[i][j];
		for (int i = 0; i < nant[0]; i++)
			for (int j = 0; j < 3; j++)
				n_xyz[i][j] = temppos[i][j];
		// alpha > 90 degree, bata > 90 degree
		if (vec[0][2] > 0)
		{
			vec[0][0] = -vec[0][0];
			vec[0][2] = -vec[0][2];
			for (int i = 0; i < nant[0]; i++)
				n_xyz[i][0] = 1 - n_xyz[i][0];
		}
		if (vec[1][2] > 0)
		{
			for (int i = 0; i < 3; i++)
				vec[1][i] = -vec[1][i];
			for (int i = 0; i < nant[0]; i++)
				n_xyz[i][1] = 1 - n_xyz[i][1];
		}
	}

	//Monoclinic lattice : c < a
	if (spacegroup_number >= 3 && spacegroup_number < 16)
	{
		lega = sqrt(vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]);
		legc = sqrt(vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]);
		if (lega < legc)
		{
			swap(vec[0][0], vec[2][2]);
			vec[0][2] = vec[2][0];
			vec[2][0] = 0.00;
			for (int j = 0; j < nant[0]; j++)
				swap(n_xyz[j][0], n_xyz[j][2]);
		}
		if (lega > legc)
		{
			temp = vec[2][2];
			vec[2][2] = legc;
			vec[0][0] = lega * temp / legc;
			vec[0][2] = lega * vec[2][0] / legc;
			vec[2][0] = 0.00;
		}
	}

	//Othorhombic lattice : c < a < b
	if (spacegroup_number >= 16 && spacegroup_number < 75)
	{
		for (int i = 0; i < 2; i++)
		{
			if (vec[2][2] > vec[i][i])
			{
				swap(vec[2][2], vec[i][i]);
				for (int j = 0; j < nant[0]; j++)
					swap(n_xyz[j][2], n_xyz[j][i]);
			}
		}
		if (vec[0][0] > vec[1][1])
		{
			swap(vec[0][0], vec[1][1]);
			for (int j = 0; j < nant[0]; j++)
				swap(n_xyz[j][0], n_xyz[j][1]);
		}
	}

	show_cell(vec, n_xyz, n_types, nant[0]);
	FILE* fp2 = fopen(file2, "w");
	savposcar(fp2, title, latt, ifix, iflg, nant, unit_typenum, elemnum, elemsym, vec, n_xyz, fix);
	fclose(fp2);
	FILE* fp_ = fopen("structure_operator", "at+");
	fprintf(fp_, "****************************************************\n");
	fclose(fp_);
}

void EV_fixatomcoor(const char file1[], const char file2[], char axis[], double start, double end, char addfix[3])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	direct_to_carts(iflg, vec, xyz, nant);
	ifix = 0;
	int ax;
	if (!strcmp(axis, "x") || !strcmp(axis, "y") || !strcmp(axis, "z"))
	{
		if (!strcmp(axis, "x")) ax = 0;
		if (!strcmp(axis, "y")) ax = 1;
		if (!strcmp(axis, "z")) ax = 2;
		for (int i = 0; i < nant[0]; i++)
		{
			if (xyz[i][ax] >= start && xyz[i][ax] <= end)
			{
				for (int j = 0; j < 3; j++)
					fix[i][j] = addfix[j];
			}
		}
	}
	double(*abc)[3];
	abc = new double[nant[0]][3];
	double norm1 = vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2];
	double norm2 = vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2];
	double norm3 = vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2];
	for (int i = 0; i < nant[0]; i++)
	{
		abc[i][0] = (xyz[i][0] * vec[0][0] + xyz[i][1] * vec[0][1] + xyz[i][2] * vec[0][2]) / norm1;
		abc[i][1] = (xyz[i][0] * vec[1][0] + xyz[i][1] * vec[1][1] + xyz[i][2] * vec[1][2]) / norm2;
		abc[i][2] = (xyz[i][0] * vec[2][0] + xyz[i][1] * vec[2][1] + xyz[i][2] * vec[2][2]) / norm3;
	}
	if (!strcmp(axis, "a") || !strcmp(axis, "b") || !strcmp(axis, "c"))
	{
		if (!strcmp(axis, "a")) ax = 0;
		if (!strcmp(axis, "b")) ax = 1;
		if (!strcmp(axis, "c")) ax = 2;
		for (int i = 0; i < nant[0]; i++)
		{
			if (abc[i][ax] >= start && abc[i][ax] <= end)
			{
				for (int j = 0; j < 3; j++)
					fix[i][j] = addfix[j];
			}
		}
	}
	FILE* fp2 = fopen("structure_operator", "at+");
	fprintf(fp2, "The atoms from %lf to %lf in %s axis\n are fixed as %c %c %c!\n", start, end, axis, addfix[0], addfix[1], addfix[2]);
	fprintf(fp2, "OUTPUT: %s\n", file2);
	fprintf(fp2, "****************************************************\n");
	fclose(fp2);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
	delete[] abc;
}

void EV_fixatomindex(const char file1[], const char file2[], int argc, char *argv[], char addfix[3])
{
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	vector<int> num;
	for (int i = 4; i < argc - 3; i++)
	{
		if (!strcmp(argv[i], "all"))
		{
			for (int i = 0; i < nant[0]; i++)
			{
				num.push_back(i + 1);
			}
			break;
		}
		else if (strstr(argv[i], "-") != NULL) // 
		{
			int start, end;
			sscanf(argv[i], "%d%*c%d", &start, &end);
			if (start > end || start<1 || end>nant[0])
			{
				printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
				return;
			}
			for (int i = start; i < end + 1; i++)
				num.push_back(i);
		}
		else if (IsNum(argv[i]))
		{
			if (atoi(argv[i]) <= 0 || atoi(argv[i]) > nant[0])
			{
				printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
				return;
			}
			num.push_back(atoi(argv[i]));
		}
		else if (!IsNum(argv[i]))
			printf("Please Input Atom_Index!\n");
	}
	//erase same value
	sort(num.begin(), num.end());
	num.erase(unique(num.begin(), num.end()), num.end());
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	ifix = 0;
	for (int i = 0; i < num.size(); i++)
	{
		for (int j = 0; j < 3; j++)
			fix[num[i] - 1][j] = addfix[j];
	}
	FILE* fp2 = fopen("structure_operator", "at+");
	fprintf(fp2, "The atom ");
	for (int i = 0; i < num.size(); i++)
		fprintf(fp2, "%d ", num[i]);
	fprintf(fp2, "are fixed as %c %c %c!\n", addfix[0], addfix[1], addfix[2]);
	fprintf(fp2, "OUTPUT: %s\n", file2);
	fprintf(fp2, "****************************************************\n");
	fclose(fp2);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
}

void EV_fixatomele(const char file1[], const char file2[], vector<string> elem, char addfix[3])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int cnt = 0;
	for (int i = 0; i < nant[1]; i++)
	{
		if (find(elem.begin(), elem.end(), elemsym[i]) != elem.end())
		{
			ifix = 0;
			for (int k = 0; k < typenum[i]; k++)
				for (int j = 0; j < 3; j++)
					fix[k + cnt][j] = addfix[j];
		}
		cnt += typenum[i];
	}
	FILE* fp2 = fopen("structure_operator", "at+");
	fprintf(fp2, "The element ");
	for (int i = 0; i < elem.size(); i++)
		fprintf(fp2, "%s ", elem[i].c_str());
	fprintf(fp2, "are fixed as %c %c %c!\n", addfix[0], addfix[1], addfix[2]);
	fprintf(fp2, "OUTPUT: %s\n", file2);
	fprintf(fp2, "****************************************************\n");
	fclose(fp2);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
}

void EV_cleanfix(const char file[])
{
	FILE* fp1 = fopen(file, "r");
	if (fp1 == NULL)
		return;
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	ifix = 1;
	FILE* fp2 = fopen(file, "w");
	savposcar(fp2, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp2);

}

void EV_carts_to_direct(const char file1[], const char file2[])
{
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		printf("ERROR!The %s is not exist!\n", file1);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix); //default direct
	fclose(fp1);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
}

void EV_direct_to_carts(const char file1[], const char file2[])
{
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		printf("ERROR!The %s is not exist!\n", file1);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);//default direct
	fclose(fp1);
	direct_to_carts(iflg, vec, xyz, nant);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
}

void EV_atomsort_coord(const char file1[], const char file2[], char mode[])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	direct_to_carts(iflg, vec, xyz, nant);
	int axis;
	if (!strcmp(mode, "x"))
		axis = 0;
	if (!strcmp(mode, "y"))
		axis = 1;
	if (!strcmp(mode, "z"))
		axis = 2;
	int count = 0;
	for (int i = 0; i < nant[1]; i++)
	{
		for (int j = 0; j < typenum[i]; j++)
		{
			for (int k = j + 1; k < typenum[i]; k++)
			{
				if (xyz[k + count][axis] < xyz[j + count][axis])
				{
					swap(xyz[k + count][0], xyz[j + count][0]);
					swap(xyz[k + count][1], xyz[j + count][1]);
					swap(xyz[k + count][2], xyz[j + count][2]);
					swap(fix[k + count][0], fix[j + count][0]);
					swap(fix[k + count][1], fix[j + count][1]);
					swap(fix[k + count][2], fix[j + count][2]);
				}
			}
		}
		count += typenum[i];
	}
	FILE* fp2 = fopen("structure_operator", "at+");
	fprintf(fp2, "Atoms are resorted according to the %s!\n", mode);
	fprintf(fp2, "OUTPUT: %s\n", file2);
	fprintf(fp2, "****************************************************\n");
	fclose(fp2);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
}

void EV_atomsort_element(const char file1[], const char file2[], char ele1[], char ele2[])
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	int flag1 = -1, flag2 = -1;
	for (int i = 0; i < nant[1]; i++)
	{
		if (!strcmp(ele1, elemsym[i]))
			flag1 = i;
		if (!strcmp(ele2, elemsym[i]))
			flag2 = i;
	}
	if (flag1 == flag2 || flag1 < 0 || flag2 < 0)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The element is wrong!\n");
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	if (flag1 > flag2)
		swap(flag1, flag2);
	vector < vector<double>> xyz1;
	vector < vector<char>> fix1;
	vector < vector<double>> xyz2;
	vector < vector<char>> fix2;
	int cnt1 = 0;
	for (int i = 0; i < nant[1]; i++)
	{
		if (i == flag1)
			for (int j = 0; j < typenum[i]; j++)
			{
				vector<double> coo = { xyz[cnt1 + j][0],xyz[cnt1 + j][1],xyz[cnt1 + j][2] };
				vector<char> str = { fix[cnt1 + j][0],fix[cnt1 + j][1],fix[cnt1 + j][2] };
				xyz1.push_back(coo);
				fix1.push_back(str);
			}
		if (i == flag2)
			for (int j = 0; j < typenum[i]; j++)
			{
				vector<double> coo = { xyz[cnt1 + j][0],xyz[cnt1 + j][1],xyz[cnt1 + j][2] };
				vector<char> str = { fix[cnt1 + j][0],fix[cnt1 + j][1],fix[cnt1 + j][2] };
				xyz2.push_back(coo);
				fix2.push_back(str);
			}
		cnt1 += typenum[i];
	}
	swap(typenum[flag1], typenum[flag2]);
	swap(elemnum[flag1], elemnum[flag2]);
	swap(elemsym[flag1], elemsym[flag2]);
	int cnt2 = 0;
	for (int i = 0; i < nant[1]; i++)
	{
		if (i == flag1)
			for (int j = 0; j < typenum[i]; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					xyz[cnt2 + j][k] = xyz2[j][k];
					fix[cnt2 + j][k] = fix2[j][k];
				}
			}
		if (i == flag2)
			for (int j = 0; j < typenum[i]; j++)
			{
				for (int k = 0; k < 3; k++)
				{
					xyz[cnt2 + j][k] = xyz1[j][k];
					fix[cnt2 + j][k] = fix1[j][k];
				}
			}
		cnt2 += typenum[i];
	}
	FILE* fp2 = fopen("structure_operator", "at+");
	fprintf(fp2, "Atoms are resorted according to the %s and %s!\n", ele1, ele2);
	fprintf(fp2, "OUTPUT: %s\n", file2);
	fprintf(fp2, "****************************************************\n");
	fclose(fp2);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
}

void EV_move(const char file1[], const char file2[], char label[], char axis[], double move[3], double mmin, double mmax)
{
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "****************************************************\n");
	fprintf(fp, "INPUT: %s\n", file1);
	fclose(fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	if (!strcmp(label, "-c"))
	{
		printf("Move in Carts coordinates!\n");
		direct_to_carts(iflg, vec, xyz, nant);
	}
	else if (!strcmp(label, "-d"))
	{
		printf("Move in Direct coordinates!\n");
		carts_to_direct(iflg, vec, xyz, nant);
	}
	else
		printf("Move in Initial coordinates!\n");
	if (mmin > mmax)
	{
		printf("There may be an error in the set range: [%lf %lf]!\n", mmin, mmax);
		return;
	}
	int ax;
	if (!strcmp("x", axis))
		ax = 0;
	else if (!strcmp("y", axis))
		ax = 1;
	else if (!strcmp("z", axis))
		ax = 2;
	else
		printf("axis is Error!\n");
	for (int i = 0; i < nant[0]; i++)
	{
		if (xyz[i][ax] >= mmin && xyz[i][ax] <= mmax)
			for (int j = 0; j < 3; j++)
				xyz[i][j] += move[j];
	}
	FILE* fp2 = fopen("structure_operator", "at+");
	fprintf(fp2, "Atoms between [%lf, %lf] on %s moved %lf on x, %lf on y, %lf on z!\n",
		mmin, mmax, axis, move[0], move[1], move[2]);
	fprintf(fp2, "OUTPUT: %s\n", file2);
	fprintf(fp2, "****************************************************\n");
	fclose(fp2);
	FILE* fp3 = fopen(file2, "w");
	savposcar(fp3, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp3);
}

void adsorbent(char basic[], char clus[], char destin[], double dis)
{
	POSCAR bas, clu, des;
	FILE* fp1 = fopen(basic, "r");
	FILE* fp2 = fopen(clus, "r");
	if (fp1 == NULL)
	{
		printf("%s IS NOT EXIST!\n", basic);
		return;
	}
	if (fp2 == NULL)
	{
		printf("%s IS NOT EXIST!\n", clus);
		return;
	}
	readposcar(fp1, bas);
	readposcar(fp2, clu);
	fclose(fp1); fclose(fp2);
	direct_to_carts(bas.iflg, bas.vec, bas.xyz, bas.nant);
	direct_to_carts(clu.iflg, clu.vec, clu.xyz, clu.nant);
	double bas_max = -1e10, clu_min = 1e10;
	for (int i = 0; i < bas.nant[0]; i++)
		bas_max = max(bas_max, bas.xyz[i][2]);
	for (int i = 0; i < clu.nant[0]; i++)
		clu_min = min(clu_min, clu.xyz[i][2]);
	double move = dis - (clu_min - bas_max);
	for (int i = 0; i < clu.nant[0]; i++)
		clu.xyz[i][2] += move;
	vector<vector<vector<double> > > pos;
	vector<vector<vector<char> > > fix;
	pos.resize(bas.nant[1]);
	fix.resize(bas.nant[1]);
	int cnt = 0;
	for (int i = 0; i < bas.nant[1]; i++)
	{
		pos[i].resize(bas.typenum[i]);
		fix[i].resize(bas.typenum[i]);
		for (int j = 0; j < bas.typenum[i]; j++)
		{
			pos[i][j].resize(3);
			fix[i][j].resize(3);
			for (int k = 0; k < 3; k++)
			{
				pos[i][j][k] = bas.xyz[cnt + j][k];
				fix[i][j][k] = bas.fix[cnt + j][k];
			}
		}
		//	pos[i].push_back(vector<double>{bas.xyz[cnt + j][0], bas.xyz[cnt + j][1], bas.xyz[cnt + j][2]});
		cnt += bas.typenum[i];
	}
	vector<string> elem;
	cnt = 0;
	for (int i = 0; i < clu.nant[1]; i++)
	{
		int flag = 0;
		for (int j = 0; j < bas.nant[1]; j++)
		{
			if (!strcmp(clu.elemsym[i], bas.elemsym[j]))
				flag = j;
		}
		if (flag != 0) // repeated
		{
			for (int j = 0; j < clu.typenum[i]; j++)
			{
				pos[flag].push_back(vector<double>{clu.xyz[cnt + j][0], clu.xyz[cnt + j][1], clu.xyz[cnt + j][2]});
				fix[flag].push_back(vector<char>{'T','T','T'});
			}
		}
		else
		{
			vector<vector<double> > tmp;
			vector<vector<char> > tmp1;
			elem.push_back(clu.elemsym[i]);
			for (int j = 0; j < clu.typenum[i]; j++)
			{
				tmp.push_back(vector<double>{clu.xyz[cnt + j][0], clu.xyz[cnt + j][1], clu.xyz[cnt + j][2]});
				tmp1.push_back(vector<char>{'T','T','T'});
			}
			pos.push_back(tmp);
			fix.push_back(tmp1);
		}
		cnt += clu.typenum[i];
	}
	for (int i = 0; i < bas.nant[1]; i++)
		strcpy(des.elemsym[i], bas.elemsym[i]);
	for (int i = 0; i < elem.size(); i++)
		strcpy(des.elemsym[i + bas.nant[1]], elem[i].c_str());
	des.nant[1] = pos.size();
	des.nant[0] = 0;
	for (int i = 0; i < pos.size(); i++)
	{
		des.typenum[i] = pos[i].size();
		for (int j = 0; j < pos[i].size(); j++)
		{
			for (int k = 0; k < 3; k++)
			{
				des.xyz[j + des.nant[0]][k] = pos[i][j][k];
				des.fix[j + des.nant[0]][k] = fix[i][j][k];
			}
		}
		des.nant[0] += des.typenum[i];
	}
	strcpy(des.title, "Generated by VASPMATE");
	des.latt = bas.latt;
	des.ifix = bas.ifix;
	des.iflg = bas.iflg;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			des.vec[i][j] = bas.vec[i][j];
	FILE* fp = fopen(destin, "w");
	savposcar(fp, des);
	fclose(fp);
	printf("Written %s file!\n", destin);
}