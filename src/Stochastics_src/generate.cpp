#include"../../include/Stochastics_include/variant.h"
#include"../../include/Stochastics_include/generate.h"
#define MAX_SITE_3D 4096
#define SYMPREC 1e-3
using namespace std;
const double MY_PI = acos(-1);
string GetValueAfterEqual(char buf[])
{
	int iat = 0;
	string value;
	for (int i = 0; i < strlen(buf); i++)
	{
		if (buf[i] == '#')
			break;
		if (buf[i] == '=')
			iat = 1;
		if (iat)
		{
			if (buf[i] != '=' && buf[i] != ' ' && buf[i] != '\t' && buf[i] != '\n' && buf[i] != '\r')
				value.push_back(buf[i]);
		}
	}
	return value;
}

void write_POTCAT(int atom[], int num)
{
	FILE* fp_POT = fopen("evops.in", "r");
	char route[80] = { 0 };
	while (!feof(fp_POT))
	{
		char temp[80];
		fgets(temp, 80, fp_POT);
		if (strstr(temp, "POT_DIR"))
			strcpy(route, GetValueAfterEqual(temp).c_str());
	}
	for (int i = 0; i < num; i++)
	{
		char arg[100];
		sprintf(arg, "cat  %s/POT_%d  >> POTCAR", route, atom[i]);
		system(arg);
	}
}

void read(const char poscar[], int atom[], int& num)
{
	FILE* fp_olog = fopen("olog", "at+");
	FILE* fp = fopen(poscar, "r");
	if (fp == NULL) fprintf(fp_olog, "cannot open %s \n", poscar);
	char   title[MAX_FNAME];                   // title of the system
	double latt;                     // Scaling Lattice parameter 
	int    ifix;                     // Select 0 or normal 1
	int    iflg;                     // Direct 0 or Cartes 1
	int    nant[2];                      // nant[0];nant[1] number of atoms and types
	int    typenum[MAX_NELEM];                 // number of atoms for each ion type
	int    elemnum[MAX_NELEM];                 // atomic number of each element
	char   elemsym[MAX_NELEM][3];              // atomic symbol of each element
	double vec[3][3];                 // lattice vectors
	double(*xyz)[3] = new double[MAX_NATOM][3];  // atomic coordinate
	char(*fix)[3] = new char[MAX_NATOM][3];  // Selective fix on each atom
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
	num = nant[1];
	for (int i = 0; i < num; i++)
		atom[i] = elemnum[i];
	fclose(fp_olog);
	delete[] xyz;
	delete[] fix;
}

void gene_first_population()  // , const char *ratio)
{
	FILE* fp_olog = fopen("../olog", "at+");
	static int counts = 0;
	counts += 1;
	fprintf(fp_olog, "first generation: pid:	%d counts %d \n", getpid(), counts);
	// cd prototype libiary 
	char lib_dir[MAX_FNAME];
	//		strcat( lib_dir,ratio);
	// access all structue and  substitute
	FILE* fp1;
	FILE* fp2;
	char   title[MAX_FNAME];                   // title of the system
	double latt;                     // Scaling Lattice parameter 
	int    ifix;                     // Select 0 or normal 1
	int    iflg;                     // Direct 0 or Cartes 1
	int    nant[2];                      // nant[0];nant[1] number of atoms and types
	int    typenum[MAX_NELEM];                 // number of atoms for each ion type
	int    elemnum[MAX_NELEM];                 // atomic number of each element
	char   elemsym[MAX_NELEM][3];              // atomic symbol of each element
	double vec[3][3];                 // lattice vectors
	double(*xyz)[3] = new double[MAX_NATOM][3];  // atomic coordinate
	char(*fix)[3] = new char[MAX_NATOM][3];  // Selective fix on each atom
	char name[9];
	int j = 0;
	DIR* dir_ptr;
	struct dirent* direntp;
	if ((dir_ptr = opendir("../lib/")) == NULL)
		fprintf(fp_olog, "lib cannot open\n");
	else
	{
		FILE* fp3 = fopen("../corre", "w");
		vector<string> filename;
		while ((direntp = readdir(dir_ptr)) != NULL)
		{
			string str1(direntp->d_name);
			if (str1 == "." || str1 == "..") continue;
			filename.push_back(str1);
		}
		sort(filename.begin(), filename.end(), [](string file1, string file2) {
			return atoi(file1.substr(0, file1.find('.')).c_str()) < atoi(file2.substr(0, file2.find('.')).c_str());
			});
		for (string name : filename)
		{
			string bin = "../lib/" + name;
			fp1 = fopen(bin.c_str(), "r");
			readposcar(fp1,
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
			fclose(fp1);
			fprintf(fp_olog, "reading:%s %d\n", name.c_str(), j);
			string savefile = "1st_" + to_string(j++);
			strcpy(title, " 1th generation");
			fprintf(fp3, "%s:	%s\n", savefile.c_str(), bin.c_str());
			//fp2=fopen(direntp->d_name,"w");
			fp2 = fopen(savefile.c_str(), "w");
			//fprintf(fp2,"%s %s\n",elemsym[0],elemsym[1] );
			savposcar(fp2,
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
				fix);   				// Selective fix on each atom
			fclose(fp2);
		}
		fclose(fp_olog);
		fclose(fp3);
		closedir(dir_ptr);
	}
	delete[] xyz;
	delete[] fix;
}

void getbox(int mode, int R_flag, double lattice[3][3], double ia, double ib, double ic, double ialpha, double ibeta, double igamma)
{
	double metric[3][3];
	metric[0][0] = ia * ia;
	metric[0][1] = ia * ib * cos(igamma / 180 * MY_PI);
	metric[0][2] = ia * ic * cos(ibeta / 180 * MY_PI);
	metric[1][0] = ia * ib * cos(igamma / 180 * MY_PI);
	metric[1][1] = ib * ib;
	metric[1][2] = ib * ic * cos(ialpha / 180 * MY_PI);
	metric[2][0] = ia * ic * cos(ibeta / 180 * MY_PI);
	metric[2][1] = ib * ic * cos(ialpha / 180 * MY_PI);
	metric[2][2] = ic * ic;

	if (mode == 1)//set_tricli(lattice, metric);
	{
		double a, b, c, alpha, beta, gamma, cg, cb, ca, sg;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		alpha = acos(metric[1][2] / b / c);
		beta = acos(metric[0][2] / a / c);
		gamma = acos(metric[0][1] / a / b);

		cg = cos(gamma);
		cb = cos(beta);
		ca = cos(alpha);
		sg = sin(gamma);

		lattice[0][0] = a;
		lattice[0][1] = b * cg;
		lattice[0][2] = c * cb;
		lattice[1][1] = b * sg;
		lattice[1][2] = c * (ca - cb * cg) / sg;
		lattice[2][2] = c * sqrt(1 - ca * ca - cb * cb - cg * cg +
			2 * ca * cb * cg) / sg;
	}

	if (mode == 2)//set_monocli(lattice, metric);
	{
		/* Lattice is expected to be C centring */
		double a, b, c, beta;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = a;
		lattice[1][1] = b;
		beta = acos(metric[0][2] / a / c);
		lattice[0][2] = c * cos(beta);
		lattice[2][2] = c * sin(beta);
	}

	if (mode == 3)//set_ortho(lattice, metric);
	{
		double a, b, c;
		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = a;
		lattice[1][1] = b;
		lattice[2][2] = c;
	}

	if (mode == 4)//set_tetra(lattice, metric);
	{
		double a, b, c;
		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = (a + b) / 2;
		lattice[1][1] = (a + b) / 2;
		lattice[2][2] = c;
	}

	if (mode == 5)//set_rhomb(lattice, metric);
	{
		if (R_flag)
		{
			double a, b, c, angle, ahex, chex;

			a = sqrt(metric[0][0]);
			b = sqrt(metric[1][1]);
			c = sqrt(metric[2][2]);
			angle = acos((metric[0][1] / a / b +
				metric[0][2] / a / c +
				metric[1][2] / b / c) / 3);

			/* Reference, http://cst-www.nrl.navy.mil/lattice/struk/rgr.html */
			ahex = 2 * (a + b + c) / 3 * sin(angle / 2);
			chex = (a + b + c) / 3 * sqrt(3 * (1 + 2 * cos(angle)));
			lattice[0][0] = ahex / 2;
			lattice[1][0] = -ahex / (2 * sqrt(3));
			lattice[2][0] = chex / 3;
			lattice[1][1] = ahex / sqrt(3);
			lattice[2][1] = chex / 3;
			lattice[0][2] = -ahex / 2;
			lattice[1][2] = -ahex / (2 * sqrt(3));
			lattice[2][2] = chex / 3;
		}
		else//set_trigo(lattice, metric);
		{
			double a, b, c;

			a = sqrt(metric[0][0]);
			b = sqrt(metric[1][1]);
			c = sqrt(metric[2][2]);
			lattice[0][0] = (a + b) / 2;
			lattice[0][1] = -(a + b) / 4;
			lattice[1][1] = (a + b) / 4 * sqrt(3);
			lattice[2][2] = c;
		}
	}
	if (mode == 6)//set_trigo(lattice, metric);
	{
		double a, b, c;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = (a + b) / 2;
		lattice[0][1] = -(a + b) / 4;
		lattice[1][1] = (a + b) / 4 * sqrt(3);
		lattice[2][2] = c;
	}

	if (mode == 7)//set_cubic(lattice, metric);
	{
		double a, b, c;

		a = sqrt(metric[0][0]);
		b = sqrt(metric[1][1]);
		c = sqrt(metric[2][2]);
		lattice[0][0] = (a + b + c) / 3;
		lattice[1][1] = (a + b + c) / 3;
		lattice[2][2] = (a + b + c) / 3;
	}
}

void SK_malloc_sitelist(int n, SK_site_list* list)
{
	if (n <= 0)
		return;
	if (list->num == 0)
	{
		list->pos = (double(*)[3])malloc(3 * n * sizeof(double));
		list->wyckoff = (int*)malloc(n * sizeof(int));
		list->equivalent_num = (int*)malloc(n * sizeof(int));
		list->equivalent_id = (int*)malloc(n * sizeof(int));
	}
	else
	{
		list->pos = (double(*)[3])realloc(list->pos, 3 * n * sizeof(double));
		list->wyckoff = (int*)realloc(list->wyckoff, n * sizeof(int));
		list->equivalent_num = (int*)realloc(list->equivalent_num, n * sizeof(int));
		list->equivalent_id = (int*)realloc(list->equivalent_id, n * sizeof(int));
	}
	list->num = n;
}

void SK_get_sitelist_3D(SK_site_list* list,
	const int sym_size,
	int(*rot)[3][3],
	const double(*trans)[3],
	SPGCONST double bravais_lattice[3][3],
	const int hall_number,
	double symprec)
{
	int i, j, k, n;
	char temp_symbol[7];
	double posible_pos[16] = {
		  0 ,1.0 / 12.0,1.0 / 8.0,
	1.0 / 6.0 ,1.0 / 4.0 ,1.0 / 3.0,
	3.0 / 8.0 ,5.0 / 12.0,1.0 / 2.0,
	7.0 / 12.0,5.0 / 8.0 ,2.0 / 3.0,
	3.0 / 4.0 ,5.0 / 6.0 ,7.0 / 8.0,
	11.0 / 12.0 };
	SK_malloc_sitelist(MAX_SITE_3D, list);
	n = 0;
	for (i = 0; i < 16; i++)
		for (j = 0; j < 16; j++)
			for (k = 0; k < 16; k++)
			{
				list->pos[n][0] = posible_pos[i];
				list->pos[n][1] = posible_pos[j];
				list->pos[n][2] = posible_pos[k];
				list->equivalent_id[n] = -1;
				n++;
			}
	int equivalent_id = 0;
	for (n = 0; n < list->num; n++)
	{
		if (list->equivalent_id[n] < 0)
		{
			double symprec2 = symprec * symprec;
			//get particle id from hallnumber
			int failed_flag;
			double vec[3];
			double dsp[3];
			int equivalent_num = 1;
			list->equivalent_id[n] = equivalent_id;
			for (i = 0; i < sym_size; i++)
			{
				failed_flag = 0;
				vec[0] = rot[i][0][0] * list->pos[n][0] + rot[i][0][1] * list->pos[n][1] + rot[i][0][2] * list->pos[n][2] + trans[i][0];
				vec[1] = rot[i][1][0] * list->pos[n][0] + rot[i][1][1] * list->pos[n][1] + rot[i][1][2] * list->pos[n][2] + trans[i][1];
				vec[2] = rot[i][2][0] * list->pos[n][0] + rot[i][2][1] * list->pos[n][1] + rot[i][2][2] * list->pos[n][2] + trans[i][2];
				while (vec[0] < -symprec)vec[0] += 1;
				while (vec[1] < -symprec)vec[1] += 1;
				while (vec[2] < -symprec)vec[2] += 1;
				while (vec[0] > 1 - symprec)vec[0] -= 1;
				while (vec[1] > 1 - symprec)vec[1] -= 1;
				while (vec[2] > 1 - symprec)vec[2] -= 1;
				for (j = n; j < list->num; j++)
				{
					dsp[0] = vec[0] - list->pos[j][0];
					dsp[1] = vec[1] - list->pos[j][1];
					dsp[2] = vec[2] - list->pos[j][2];
					while (dsp[0] < -0.5)dsp[0] += 1;
					while (dsp[1] < -0.5)dsp[1] += 1;
					while (dsp[2] < -0.5)dsp[2] += 1;
					while (dsp[0] > 0.5)dsp[0] -= 1;
					while (dsp[1] > 0.5)dsp[1] -= 1;
					while (dsp[2] > 0.5)dsp[2] -= 1;
					if ((dsp[0] * dsp[0] + dsp[1] * dsp[1] + dsp[2] * dsp[2]) < symprec2)
					{
						if (list->equivalent_id[j] != equivalent_id)
						{
							equivalent_num++;
							list->equivalent_id[j] = equivalent_id;
						}
						failed_flag = 1;
					}
				}
				if (!failed_flag)
				{
					equivalent_num++;
					SK_malloc_sitelist(list->num + 1, list);
					list->equivalent_id[list->num - 1] = equivalent_id;
					list->pos[list->num - 1][0] = vec[0];
					list->pos[list->num - 1][1] = vec[1];
					list->pos[list->num - 1][2] = vec[2];
				}
			}
			/*		list->wyckoff[n] = SK_get_Wyckoff_notation(temp_symbol,
						list->pos[n],
						sym_size,
						rot,
						trans,
						equivalent_num,
						bravais_lattice,
						hall_number,
						symprec);*/
			list->equivalent_num[n] = equivalent_num;
			for (j = n + 1; j < list->num; j++)
			{
				if (list->equivalent_id[n] == list->equivalent_id[j])//copy wyckoff data
				{
					list->wyckoff[j] = list->wyckoff[n];
					list->equivalent_num[j] = equivalent_num;
				}
			}
			equivalent_id++;
		}
	}
}

//vector<vector<double> > get_equiv_atoms(POSCAR& pos, int hall_number, int(*rot)[3][3], double(*tra)[3], double x, double y, double z, double symprec)
//{
//	int symmetry_num = spg_get_symmetry_from_database(rot, tra, hall_number);
//	pos.nant[0] = 0;
//	vector<vector<double> > xyz;
//	int i, j;
//	double symprec2 = symprec * symprec;
//	//get particle id from hallnumber
//
//	int failed_flag;
//	double vec[3];
//	double dsp[3];
//
//	for (i = 0; i < symmetry_num; i++)
//	{
//		failed_flag = 0;
//		vec[0] = rot[i][0][0] * x + rot[i][0][1] * y + rot[i][0][2] * z + tra[i][0];
//		vec[1] = rot[i][1][0] * x + rot[i][1][1] * y + rot[i][1][2] * z + tra[i][1];
//		vec[2] = rot[i][2][0] * x + rot[i][2][1] * y + rot[i][2][2] * z + tra[i][2];
//		for (j = 0; j < pos.nant[0]; j++)
//		{
//			dsp[0] = vec[0] - pos.xyz[j][0];
//			dsp[1] = vec[1] - pos.xyz[j][1];
//			dsp[2] = vec[2] - pos.xyz[j][2];
//			while (dsp[0] < -0.5)dsp[0] += 1;
//			while (dsp[1] < -0.5)dsp[1] += 1;
//			while (dsp[2] < -0.5)dsp[2] += 1;
//			while (dsp[0] > 0.5)dsp[0] -= 1;
//			while (dsp[1] > 0.5)dsp[1] -= 1;
//			while (dsp[2] > 0.5)dsp[2] -= 1;
//			if ((dsp[0] * dsp[0] + dsp[1] * dsp[1] + dsp[2] * dsp[2]) < symprec2)
//			{
//				failed_flag = 1;
//				break;
//			}
//		}
//		if (!failed_flag)
//		{
//			for (int i = 0; i < 3; i++)
//			{
//				while (vec[i] < 0) vec[i]++;
//				while (vec[i] >= 1) vec[i]--;
//			}
//			pos.nant[0]++;
//			xyz.push_back(vector<double>{vec[0], vec[1], vec[2]});
//		}
//	}
//	sort(xyz.begin(), xyz.end());
//	xyz.erase(unique(xyz.begin(), xyz.end()), xyz.end());
//	pos.nant[0] = xyz.size();
//	return xyz;
//}

void set_lattice(int type, double& a, double& b, double& c, double& alpha, double& beta, double& gamma)
{
	// rand_length: 2-12  rand_angle: 20 - 160
	if (type == 1)//Triclinic
	{
		a = rand() % 10000 / 1000.0 + 2;
		b = rand() % 10000 / 1000.0 + 2;
		c = rand() % 10000 / 1000.0 + 2;
		alpha = rand() % 14000 / 100.0 + 20;
		beta = rand() % 14000 / 100.0 + 20;
		gamma = rand() % 14000 / 100.0 + 20;
		while (gamma >= alpha + beta || alpha >= beta + gamma || beta >= alpha + gamma)
		{
			alpha = rand() % 14000 / 100.0 + 20;
			beta = rand() % 14000 / 100.0 + 20;
			gamma = rand() % 14000 / 100.0 + 20;
		}
	}
	if (type == 2)//Monoclinic
	{
		alpha = gamma = 90;
		a = rand() % 10000 / 1000.0 + 2;
		b = rand() % 10000 / 1000.0 + 2;
		c = rand() % 10000 / 1000.0 + 2;
		beta = rand() % 14000 / 100.0 + 20;
	}
	if (type == 3)//Orthorhombic
	{
		alpha = beta = gamma = 90;
		a = rand() % 10000 / 1000.0 + 2;
		b = rand() % 10000 / 1000.0 + 2;
		c = rand() % 10000 / 1000.0 + 2;
	}
	if (type == 4)//Tetragonal
	{
		alpha = beta = gamma = 90;
		a = b = rand() % 10000 / 1000.0 + 2;
		c = rand() % 10000 / 1000.0 + 2;
	}
	if (type == 5)//Trigonal
	{
		alpha = beta = 90;
		gamma = 120;
		a = b = rand() % 10000 / 1000.0 + 2;
		c = rand() % 10000 / 1000.0 + 2;
	}
	if (type == 6)//Hexagonal
	{
		alpha = beta = 90;
		gamma = 120;
		a = b = rand() % 10000 / 1000.0 + 2;
		c = rand() % 10000 / 1000.0 + 2;
	}
	if (type == 7)//Cubic
	{
		alpha = beta = gamma = 90;
		a = b = c = rand() % 10000 / 1000.0 + 2;
	}
}

void gene_Random_structrue(vector<int> atom, vector<string> element, int max_formula_atom, const char* file, double input_volume)
{
	if (atom.size() != element.size())
	{
		printf("The type of atoms is different from element!\n");
		return;
	}
	POSCAR pos;
	int hall_number = rand() % 423 + 108;// 108-530
	//printf("%d\n", hall_number);
	while (hall_number <= 107)
		hall_number = rand() % 530 + 1;
	int rot[192][3][3];
	double tra[192][3];
	// double x = rand() % 24 / 24.0, y = rand() % 24 / 24.0, z = rand() % 24 / 24.0;
	int mode;
	int R_flag = 1;
	double lattice[3][3] = { 0 }, ia, ib, ic, ialpha, ibeta, igamma;
	SpglibSpacegroupType SPG = spg_get_spacegroup_type(hall_number);
	int spacegroup_number = SPG.number;
	if (spacegroup_number == 1 || spacegroup_number == 2)
		mode = 1;
	if (spacegroup_number >= 3 && spacegroup_number < 16)
		mode = 2;
	if (spacegroup_number >= 16 && spacegroup_number < 75)
		mode = 3;
	if (spacegroup_number >= 75 && spacegroup_number < 143)
		mode = 4;
	if (spacegroup_number >= 143 && spacegroup_number < 168)
		mode = 5;
	if (spacegroup_number >= 168 && spacegroup_number < 195)
		mode = 6;
	if (spacegroup_number >= 195 && spacegroup_number < 230)
		mode = 7;
	//vector<vector<double> > xyz = get_equiv_atoms(pos, hall_number, rot, tra, x, y, z, SYMPREC);
	set_lattice(mode, ia, ib, ic, ialpha, ibeta, igamma);
	getbox(mode, R_flag, lattice, ia, ib, ic, ialpha, ibeta, igamma);
	SK_site_list* list = new SK_site_list;
	list->num = 0;
	int sym_size = spg_get_symmetry_from_database(rot, tra, hall_number);
	SK_get_sitelist_3D(list, sym_size, rot, tra, lattice, hall_number, SYMPREC);
	vector<int> norepeat_num;
	vector<int> norepeat_id;
	for (int i = 0; i < list->num; i++)
	{
		if (find(norepeat_id.begin(), norepeat_id.end(), list->equivalent_id[i]) == norepeat_id.end())
		{
			norepeat_num.push_back(list->equivalent_num[i]);
			norepeat_id.push_back(list->equivalent_id[i]);
		}
	}
	//sort(norepeat_num.begin(), norepeat_num.end());
	int total_num = 0;
	for (int i = 0; i < atom.size(); i++)
		total_num += atom[i];
	vector<pair<int, int> > v;
	for (int i = 0; i < norepeat_num.size(); i++)
		v.push_back(pair<int, int>(i, norepeat_num[i]));
	closest close;
	sort(v.begin(), v.end(), mysort);
	vector<int> subscript = close.getclosest(v, max_formula_atom);
	if (subscript.size() == 0)
		subscript.push_back(0);
	vector<int> real_subscript(subscript.size());
	for (int i = 0; i < subscript.size(); i++)
		real_subscript[i] = v[subscript[i]].first;
	vector<int> selectnum;
	vector<int> selectid;
	for (int i = 0; i < real_subscript.size(); i++)
	{
		selectnum.push_back(norepeat_num[real_subscript[i]]);
		selectid.push_back(norepeat_id[real_subscript[i]]);
	}
	vector<vector<double> > selectpos;
	for (int i = 0; i < list->num; i++)
	{
		if (find(selectnum.begin(), selectnum.end(), list->equivalent_num[i]) != selectnum.end() &&
			find(selectid.begin(), selectid.end(), list->equivalent_id[i]) != selectid.end())
			selectpos.push_back(vector<double>{list->pos[i][0], list->pos[i][1], list->pos[i][2]});
	}
	if (selectpos.size() >= max_formula_atom)
	{
		for (int i = 0; i < max_formula_atom; i++)
			for (int j = 0; j < 3; j++)
				pos.xyz[i][j] = selectpos[i][j];
	}
	else if (selectpos.size() < max_formula_atom)
	{
		for (int i = 0; i < selectpos.size(); i++)
			for (int j = 0; j < 3; j++)
				pos.xyz[i][j] = selectpos[i][j];
		int cnt = selectpos.size();
		int lack_num = max_formula_atom - selectpos.size();
		int extra = 0;
		while (lack_num > 0)
		{
			for (int i = extra; i < norepeat_num.size(); i++)
				if (find(real_subscript.begin(), real_subscript.end(), extra) == real_subscript.end())
				{
					extra = i;
					break;
				}
			for (int i = 0; i < list->num; i++)
			{
				if (list->equivalent_num[i] == norepeat_num[extra] && list->equivalent_id[i] == norepeat_id[extra])
				{
					for (int j = 0; j < 3; j++)
						pos.xyz[cnt][j] = list->pos[i][j];
					cnt++;
					lack_num--;
					if (lack_num == 0)
						break;
				}
			}
		}
	}
	/*if (total_num < pos.nant[0])
	{
		vector<int> random(total_num, pos.nant[0]);
		for (int i = 0; i < total_num; i++)
		{
			int position = rand() % pos.nant[0];
			while (find(random.begin(), random.end(), position) != random.end())
				position = rand() % pos.nant[0];
			random[i] = position;
		}
		int cnt = 0;
		for (int i = 0; i < pos.nant[0]; i++)
		{
			if (find(random.begin(), random.end(), i) != random.end())
				for (int j = 0; j < 3; j++)
					pos.xyz[i][j] = xyz[i][j];
		}
	}
	else if (total_num > pos.nant[0])
	{
		for (int i = 0; i < total_num; i++)
			for (int j = 0; j < 3; j++)
				pos.xyz[i][j] = (i < pos.nant[0]) ? xyz[i][j] : (rand() % 10000 / 10000.0);
	}*/
	pos.nant[0] = max_formula_atom;
	strcpy(pos.title, file);
	pos.latt = 1;
	pos.ifix = 1;
	pos.iflg = 0;
	pos.nant[1] = element.size();
	for (int i = 0; i < atom.size(); i++)
		pos.typenum[i] = atom[i];
	for (int i = 0; i < element.size(); i++)
		strcpy(pos.elemsym[i], element[i].c_str());
	for (int i = 0; i < pos.nant[1]; i++)
	{
		int an;
		char at[3];
		strcpy(at, pos.elemsym[i]);
		getAtomNum(at, an);
		pos.elemnum[i] = an;
	}
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			pos.vec[i][j] = lattice[i][j];
	if (pos.nant[0] <= total_num)
	{
		int times = total_num % pos.nant[0] == 0 ? total_num / pos.nant[0] : total_num / pos.nant[0] + 1;
		int min = 1e10;
		vector<int> super(3, 1);
		for (int i = 1; i <= (int)pow(times, 1.0 / 3); i++)
		{
			int left = 1, right = times;
			while (left <= right)
			{
				if (i * left * right == times)
				{
					if ((i + left + right) < min)
					{
						super[0] = i;
						super[1] = left;
						super[2] = right;
						min = super[0] + super[1] + super[2];
					}
					left++;
				}
				else if (i * left * right > times)
					right--;
				else if (i * left * right < times)
					left++;
			}
		}
		double ntemp_vec[3][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				ntemp_vec[i][j] = pos.vec[i][j] * super[i];
		double(*n1_xyz)[3];
		n1_xyz = new double[super[0] * pos.nant[0]][3];
		int k1 = 0;
		for (int i = 0; i < super[0] * pos.nant[0]; i++)
		{
			n1_xyz[i][0] = (pos.xyz[i / super[0]][0] + k1) / super[0];
			n1_xyz[i][1] = pos.xyz[i / super[0]][1];
			n1_xyz[i][2] = pos.xyz[i / super[0]][2];
			if (k1 < super[0] - 1)
				k1++;
			else
				k1 = 0;
		}
		pos.nant[0] *= super[0];
		double(*n2_xyz)[3];
		n2_xyz = new double[pos.nant[0] * super[1]][3];
		int k2 = 0;
		for (int i = 0; i < super[1] * pos.nant[0]; i++)
		{
			n2_xyz[i][1] = (n1_xyz[i / super[1]][1] + k2) / super[1];
			n2_xyz[i][0] = n1_xyz[i / super[1]][0];
			n2_xyz[i][2] = n1_xyz[i / super[1]][2];
			if (k2 < super[1] - 1)
				k2++;
			else
				k2 = 0;
		}
		pos.nant[0] *= super[1];
		double(*n3_xyz)[3];
		n3_xyz = new double[pos.nant[0] * super[2]][3];
		int k3 = 0;
		for (int i = 0; i < super[2] * pos.nant[0]; i++)
		{
			n3_xyz[i][2] = (n2_xyz[i / super[2]][2] + k3) / super[2];
			n3_xyz[i][1] = n2_xyz[i / super[2]][1];
			n3_xyz[i][0] = n2_xyz[i / super[2]][0];
			if (k3 < super[2] - 1)
				k3++;
			else
				k3 = 0;
		}
		pos.nant[0] *= super[2];
		for (int i = 0; i < pos.nant[0]; i++)
			for (int j = 0; j < 3; j++)
				pos.xyz[i][j] = n3_xyz[i][j];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				pos.vec[i][j] = ntemp_vec[i][j];
		delete[] n1_xyz;
		delete[] n2_xyz;
		delete[] n3_xyz;
	}
	// fix vec 
	double cur_volume = volume(pos.vec);
	if (input_volume == 0)
		input_volume = cur_volume;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			pos.vec[i][j] *= pow(input_volume / cur_volume, 1.0 / 3);
	pos.nant[0] = total_num;
	FILE* fp = fopen(file, "w");
	savposcar(fp, pos);
	fclose(fp);
}

vector<string> gene_first_population_Random(vector<int> atom, vector<string> element, int POPULATION_NUM, int formula_num, double volume)
{
	int total_num = 0;
	vector<string> ans(POPULATION_NUM);
	for (int i = 0; i < atom.size(); i++)
		total_num += atom[i];
	if (formula_num > total_num)
		printf("We suggest set FormulaAtomNum = %d\n", total_num);
	for (int i = 0; i < POPULATION_NUM; i++)
	{
		char file[50];
		sprintf(file, "POSCAR%d", i);
		gene_Random_structrue(atom, element, formula_num, file, volume);
		//check atom distance
		int MAX = 100;
		FILE* fp = fopen(file, "r");
		POSCAR pos;
		readposcar(fp, pos);
		fclose(fp);
		while (!DistAtom(pos.nant, pos.latt, pos.vec, pos.xyz) && MAX > 0)
		{
			gene_Random_structrue(atom, element, formula_num, file, volume);
			MAX--;
		}
		ans[i] = file;
		//force_fix(file, file, max_step, scale);
	}
	return ans;
}

