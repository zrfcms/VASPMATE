#include"../include/read_write.h"
#include<vector>
#include<algorithm>
#include<iostream>
using std::vector;

void has_vacuum_slab_or_not(POSCAR pos, bool has_vacuum_slab[3])
{
	for (int i = 0; i < 3; i++)
		has_vacuum_slab[i] = determine_has_vacuum_slab(pos, i);
}

bool determine_has_vacuum_slab(POSCAR pos, int direction) //direction = 0,1,2
{
	carts_to_direct(pos.iflg, pos.vec, pos.xyz, pos.nant);
	double threshold = 9;
	vector<vector<double> > frac_coor(pos.nant[0], vector<double>(3));
	for (int i = 0; i < pos.nant[0]; i++)
		for (int j = 0; j < 3; j++)
			frac_coor[i][j] = pos.xyz[i][j];
	sort(frac_coor.begin(), frac_coor.end());
	int n = frac_coor.size();
	double latt_length = 0;
	for (int i = 0; i < 3; i++)
		latt_length += pos.vec[direction][i] * pos.vec[direction][i];
	latt_length = sqrt(latt_length);
	if (frac_coor[0][direction] == frac_coor[n - 1][direction])
	{
		if (latt_length >= threshold)
			return true;
	}
	else
	{
		if ((latt_length - (frac_coor[n - 1][direction] - frac_coor[0][direction]) * latt_length) >= threshold)
			return true;
		else
		{
			for (int i = 1; i < n; i++)
			{
				if ((frac_coor[i][direction] - frac_coor[i - 1][direction]) * latt_length >= threshold)
					return true;
			}
		}
	}
	return false;
}

POSCAR delete_atom(POSCAR pos, int n)
{
	if (n >= pos.nant[0])
		return pos;
	int cnt = 0, cnt1 = 0;
	POSCAR pos1;
	char   title[MAX_NCHAR];        // Title of system
	double latt;                    // Scaling factor
	int    ifix;                    // Select 0 or normal 1
	int    iflg;                    // Direct 0 or Cartes 1
	int    nant[2];                 // number of atoms and types
	int    typenum[MAX_NELEM];      // number of atoms of each type
	int    elemnum[MAX_NELEM];      // atomic number
	char   elemsym[MAX_NELEM][3];   // atomic symbol
	double vec[3][3];               // lattice vector
	double xyz[MAX_NATOM][3];		// atomic coordinate
	char fix[MAX_NATOM][3];         // Selective fix on each atom
	strcpy(pos1.title, pos.title);
	pos1.latt = pos.latt;
	pos1.ifix = pos.ifix;
	pos1.iflg = pos.iflg;
	pos1.nant[0] = pos.nant[0];
	pos1.nant[1] = pos.nant[1];
	vector<int> n_typenum(pos.nant[1]);
	for (int i = 0; i < pos.nant[1]; i++)
		n_typenum[i] = pos.typenum[i];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			pos1.vec[i][j] = pos.vec[i][j];
	for (int i = 0; i < pos.nant[1]; i++)
	{
		for (int j = 0; j < pos.typenum[i]; j++)
		{
			if (cnt + j == n)
			{
				n_typenum[i]--;
				pos1.nant[0] --;
			}
			else
			{
				for (int k = 0; k < 3; k++)
				{
					pos1.xyz[cnt1][k] = pos.xyz[cnt + j][k];
					if (pos.ifix == 0)
						pos1.fix[cnt1][k] = pos.fix[cnt + j][k];
				}
				cnt1++;
			}
		}
		cnt += pos.typenum[i];
	}
	int new_ele = pos1.nant[1];
	int cnt2 = 0;
	for (int i = 0; i < pos1.nant[1]; i++)
	{
		if (n_typenum[i] == 0)
			new_ele--;
		else
		{
			pos1.typenum[cnt2] = n_typenum[i];
			pos1.elemnum[cnt2] = pos.elemnum[i];
			strcpy(pos1.elemsym[cnt2], pos.elemsym[i]);
			cnt2++;
		}
	}
	pos1.nant[1] = new_ele;
	return pos1;
}

void carts_to_direct(int& iflg, double vec[3][3], double xyz[][3], int nant[])
{
	if (iflg != 1)
		return;
	//translate Cartesian to Direct
	double vec_inverse[3][3];
	double vec_value = vec[0][0] * vec[1][1] * vec[2][2] + vec[0][1] * vec[1][2] * vec[2][0] + vec[1][0] * vec[2][1] * vec[0][2] -
		vec[0][2] * vec[1][1] * vec[2][0] - vec[0][1] * vec[1][0] * vec[2][2] - vec[0][0] * vec[1][2] * vec[2][1];
	if (vec_value == 0)
		return;
	vec_inverse[0][0] = (vec[1][1] * vec[2][2] - vec[1][2] * vec[2][1]) / vec_value;
	vec_inverse[0][1] = (-vec[1][0] * vec[2][2] + vec[1][2] * vec[2][0]) / vec_value;
	vec_inverse[0][2] = (vec[1][0] * vec[2][1] - vec[1][1] * vec[2][0]) / vec_value;
	vec_inverse[1][0] = (-vec[0][1] * vec[2][2] + vec[0][2] * vec[2][1]) / vec_value;
	vec_inverse[1][1] = (vec[0][0] * vec[2][2] - vec[0][2] * vec[2][0]) / vec_value;
	vec_inverse[1][2] = (-vec[0][0] * vec[2][1] + vec[0][1] * vec[2][0]) / vec_value;
	vec_inverse[2][0] = (vec[0][1] * vec[1][2] - vec[0][2] * vec[1][1]) / vec_value;
	vec_inverse[2][1] = (-vec[0][0] * vec[1][2] + vec[0][2] * vec[1][0]) / vec_value;
	vec_inverse[2][2] = (vec[0][0] * vec[1][1] - vec[0][1] * vec[1][0]) / vec_value;
	for (int i = 0; i < nant[0]; i++)
	{
		double p[3] = { xyz[i][0],xyz[i][1],xyz[i][2] };
		for (int j = 0; j < 3; j++)
		{
			xyz[i][j] = p[0] * vec_inverse[j][0] + p[1] * vec_inverse[j][1] + p[2] * vec_inverse[j][2];
			//xyz[i][j] = xyz[i][j]*vec_inverse[j][0] + xyz[i][j] * vec_inverse[j][1] + xyz[i][j] * vec_inverse[j][2];
		}
	}
	iflg = 0;
}

void direct_to_carts(int& iflg, double vec[3][3], double xyz[][3], int nant[])
{
	if (iflg != 0)
		return;
	for (int i = 0; i < nant[0]; i++)
	{
		double p[3] = { xyz[i][0],xyz[i][1],xyz[i][2] };
		for (int j = 0; j < 3; j++)
		{
			xyz[i][j] = p[0] * vec[0][j] + p[1] * vec[1][j] + p[2] * vec[2][j];
		}
	}
	iflg = 1;
}

void getAtomNum(char elemsym[3],    // String to store Element name
	int& elemnum        // atomic number
) {
	char atomstr[][3] = {
	"Nv", //   0
	"H",  //   1
	"He", //   2
	"Li", //   3
	"Be", //   4
	"B",  //   5
	"C",  //   6
	"N",  //   7
	"O",  //   8
	"F",  //   9
	"Ne", //  10
	"Na", //  11
	"Mg", //  12
	"Al", //  13
	"Si", //  14
	"P",  //  15
	"S",  //  16
	"Cl", //  17
	"Ar", //  18
	"K",  //  19
	"Ca", //  20
	"Sc", //  21
	"Ti", //  22
	"V",  //  23
	"Cr", //  24
	"Mn", //  25
	"Fe", //  26
	"Co", //  27
	"Ni", //  28
	"Cu", //  29
	"Zn", //  30
	"Ga", //  31
	"Ge", //  32
	"As", //  33
	"Se", //  34
	"Br", //  35
	"Kr", //  36
	"Rb", //  37
	"Sr", //  38
	"Y",  //  39
	"Zr", //  40
	"Nb", //  41
	"Mo", //  42
	"Tc", //  43
	"Ru", //  44
	"Rh", //  45
	"Pd", //  46
	"Ag", //  47
	"Cd", //  48
	"In", //  49
	"Sn", //  50
	"Sb", //  51
	"Te", //  52
	"I",  //  53
	"Xe", //  54
	"Cs", //  55
	"Ba", //  56
	"La", //  57
	"Ce", //  58
	"Pr", //  59
	"Nd", //  60
	"Pm", //  61
	"Sm", //  62
	"Eu", //  63
	"Gd", //  64
	"Tb", //  65
	"Dy", //  66
	"Ho", //  67
	"Er", //  68
	"Tm", //  69
	"Yb", //  70
	"Lu", //  71
	"Hf", //  72
	"Ta", //  73
	"W",  //  74
	"Re", //  75
	"Os", //  76
	"Ir", //  77
	"Pt", //  78
	"Au", //  79
	"Hg", //  80
	"Tl", //  81
	"Pb", //  82
	"Bi", //  83
	"Po", //  84
	"At", //  85
	"Rn", //  86
	"Fr", //  87
	"Ra", //  88
	"Ac", //  89
	"Th", //  90
	"Pa", //  91
	"U",  //  92
	"Np", //  93
	"Pu", //  94
	"Am", //  95
	"Cm", //  96
	"Bk", //  97
	"Cf", //  98
	"Es", //  99
	"Fm", // 100
	"Md", // 101
	"No", // 102
	"Lr", // 103
	"Rf", // 104
	"Db", // 105
	"Sg", // 106
	"Bh", // 107
	"Hs", // 108
	"Mt"  // 109
	};

	int i, j;
	for (j = 0; j < 109; j++) {
		if (strcmp(elemsym, atomstr[j]) == 0 || (atoi(elemsym) == j && j != 0)) {
			elemnum = j; break;
		}
	}
}

void getAtomSym(const int& elemnum,     // atomic number
	char elemsym[3]    // String to store Element name
) {
	char atomstr[][3] = {
	"Nv", //   0
	"H",  //   1
	"He", //   2
	"Li", //   3
	"Be", //   4
	"B",  //   5
	"C",  //   6
	"N",  //   7
	"O",  //   8
	"F",  //   9
	"Ne", //  10
	"Na", //  11
	"Mg", //  12
	"Al", //  13
	"Si", //  14
	"P",  //  15
	"S",  //  16
	"Cl", //  17
	"Ar", //  18
	"K",  //  19
	"Ca", //  20
	"Sc", //  21
	"Ti", //  22
	"V",  //  23
	"Cr", //  24
	"Mn", //  25
	"Fe", //  26
	"Co", //  27
	"Ni", //  28
	"Cu", //  29
	"Zn", //  30
	"Ga", //  31
	"Ge", //  32
	"As", //  33
	"Se", //  34
	"Br", //  35
	"Kr", //  36
	"Rb", //  37
	"Sr", //  38
	"Y",  //  39
	"Zr", //  40
	"Nb", //  41
	"Mo", //  42
	"Tc", //  43
	"Ru", //  44
	"Rh", //  45
	"Pd", //  46
	"Ag", //  47
	"Cd", //  48
	"In", //  49
	"Sn", //  50
	"Sb", //  51
	"Te", //  52
	"I",  //  53
	"Xe", //  54
	"Cs", //  55
	"Ba", //  56
	"La", //  57
	"Ce", //  58
	"Pr", //  59
	"Nd", //  60
	"Pm", //  61
	"Sm", //  62
	"Eu", //  63
	"Gd", //  64
	"Tb", //  65
	"Dy", //  66
	"Ho", //  67
	"Er", //  68
	"Tm", //  69
	"Yb", //  70
	"Lu", //  71
	"Hf", //  72
	"Ta", //  73
	"W",  //  74
	"Re", //  75
	"Os", //  76
	"Ir", //  77
	"Pt", //  78
	"Au", //  79
	"Hg", //  80
	"Tl", //  81
	"Pb", //  82
	"Bi", //  83
	"Po", //  84
	"At", //  85
	"Rn", //  86
	"Fr", //  87
	"Ra", //  88
	"Ac", //  89
	"Th", //  90
	"Pa", //  91
	"U",  //  92
	"Np", //  93
	"Pu", //  94
	"Am", //  95
	"Cm", //  96
	"Bk", //  97
	"Cf", //  98
	"Es", //  99
	"Fm", // 100
	"Md", // 101
	"No", // 102
	"Lr", // 103
	"Rf", // 104
	"Db", // 105
	"Sg", // 106
	"Bh", // 107
	"Hs", // 108
	"Mt"  // 109
	};

	strcpy(elemsym, atomstr[elemnum]);
}


int savposcar(FILE* fp,
	char   title[],                   // title of the system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                      // nant[0],nant[1] number of atoms and types
	int    typenum[],                 // number of atoms for each ion type
	int    elemnum[],                 // atomic number of each element
	char   elemsym[][3],              // atomic symbol of each element
	double vec[3][3],                 // lattice vectors
	double xyz[][3],                  // atomic coordinates
	char   fix[][3]                   // Selective fix on each atom
) {
	//	printf(" >>> Saving to ASCII file\n");

	int  i, j;
	//	printf(" Title= %s created by spl2gsg.x \n",title);
	fprintf(fp, "%s \n", title);

	//	printf("%16.9lf \n", latt);
	fprintf(fp, "%16.9lf \n", latt);
	fprintf(fp, "  %16.9lf %16.9lf %16.9lf \n  %16.9lf %16.9lf %16.9lf \n  %16.9lf %16.9lf %16.9lf \n",
		vec[0][0], vec[0][1], vec[0][2], vec[1][0], vec[1][1], vec[1][2], vec[2][0], vec[2][1], vec[2][2]);

	for (i = 0; i < nant[1]; i++)
	{
		//		printf("    %s  ", elemsym[i]);
		fprintf(fp, "    %s  ", elemsym[i]);
	}
	fprintf(fp, "\n");

	for (i = 0; i < nant[1]; i++)
	{
		//		printf(" %5d ", typenum[i]);   
		fprintf(fp, " %5d ", typenum[i]);
	}
	fprintf(fp, "\n");

	if (ifix == 0)
	{
		fprintf(fp, "Select \n");
	}

	if (iflg == 0)
	{
		fprintf(fp, "Direct \n");
	}
	else if (iflg == 1)
	{
		fprintf(fp, "Cartes \n");
	}
	else
	{
		fprintf(fp, " Wrong with flag read from bin2pos.x!! \n");
	}

	for (j = 0; j < nant[0]; j++)
	{
		if (ifix == 0)
		{
			//	printf("  %16.9lf %16.9lf %16.9lf %c %c %c \n",
			//		xyz[j][0],xyz[j][1],xyz[j][2],fix[j][0],fix[j][1],fix[j][2]);
			fprintf(fp, "  %16.9lf %16.9lf %16.9lf %c %c %c \n",
				xyz[j][0], xyz[j][1], xyz[j][2], fix[j][0], fix[j][1], fix[j][2]);
		}
		else
		{
			//printf("  %16.9lf %16.9lf %16.9lf \n",xyz[j][0],xyz[j][1],xyz[j][2]);
			fprintf(fp, "  %16.9lf %16.9lf %16.9lf \n", xyz[j][0], xyz[j][1], xyz[j][2]);
		}
	}

	return 0;
}

int readposcar(FILE* fp,
	char   title[],                   // Title of system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                    // Number of atoms and types
	int    typenum[],                 // number of atoms for each element
	int    elemnum[],                 // atomic number for each ion type
	char   elemsym[][3],              // atomic symbol for each ion type
	double vec[3][3],                 // Normal lattice vectors
	double xyz[][3],                  // Normal coordinates for all atoms
	char   fix[][3]                   // Selective fix for each atom
) {
	//	printf(" >>> Reading ASCII file\n");
	int    iline;
	int    flag;
	char   tline[MAX_NCOLU];
	double sf;                                     // scaling factor
	int    i, j;
	int    pos;
	char   elemnumstr[MAX_NELEM][5];
	char   str[10];
	int    iat;                                    // index of atoms
	double x, y, z;
	char   fx, fy, fz;

	iline = 0;
	nant[0] = 0;
	nant[1] = 0;
	flag = 2;
	ifix = 0;
	iflg = 0;
	iat = 0;
	while (1) {                                       //loop each line to the end
		fgets(tline, MAX_NCOLU, fp);             //Read one line
		iline++;                                  //move to next line 
		if (feof(fp) != 0) break;                 //Check end of the file: yes exit loop

		// read title from line #1
		if (iline == 1) {
			for (int i = 0; i < MAX_FNAME; i++)
			{
				if (tline[i] == '\n') { break; }
				title[i] = tline[i];
			}
			//				printf(" Title= %s\n", title);
			continue;
		}

		// read scaling factor from line #2
		if (iline == 2) {
			sscanf(tline, "%lf", &sf);
			latt = sf;
#if DEBUG
			//					printf("%16.9lf \n", latt);
#endif
			continue;
		}

		// read lattice vectors a from line #3
		if (iline == 3) {
			sscanf(tline, "%lf%lf%lf", &vec[0][0], &vec[0][1], &vec[0][2]);
#if DEBUG
			//					printf("  %16.9lf %16.9lf %16.9lf \n", vec[0][0],vec[0][1],vec[0][2]);
#endif
			continue;
		}

		// read lattice vectors b from line #4
		if (iline == 4) {
			sscanf(tline, "%lf%lf%lf", &vec[1][0], &vec[1][1], &vec[1][2]);
			//	    		printf("  %16.9lf %16.9lf %16.9lf \n", vec[1][0],vec[1][1],vec[1][2]);
			continue;
		}

		// read lattice vectors c from line #5
		if (iline == 5) {
			sscanf(tline, "%lf%lf%lf", &vec[2][0], &vec[2][1], &vec[2][2]);
			//		printf("  %16.9lf %16.9lf %16.9lf \n", vec[2][0],vec[2][1],vec[2][2]);
			continue;
		}

		// read elemental symbols from line #6
		if (iline == 6) {
			for (i = 0; i < MAX_NELEM; i++) {
				for (j = 0; j < strlen(tline); j++) { // get first non-blank position and string
					if (tline[j] != ' ' && tline[j] != '\t' && tline[j] != '\n' && tline[j] != '\r') {
						pos = j; break;      //exit the innermost enclosing loop
					}
					else if (tline[j] == '\n' || tline[j] == '\r') {
						pos = 9999; break;   //exit the innermost enclosing loop
					}
				}

				if (pos == 9999) break;
				sscanf(tline, "%s", elemsym[i]);
				for (j = 0; j < strlen(elemsym[i]); j++) {
					tline[pos + j] = ' ';
				}
			}
			nant[1] = i;               //Total number of elemental types

			for (i = 0; i < nant[1]; i++)
			{
				//    		printf("    %s  ", elemsym[i]);
			}
			//		printf("\n");

			continue;
		}

		// read number of atoms for each ion type
		if (iline == 7) {
			for (i = 0; i < nant[1]; i++) {
				for (j = 0; j < strlen(tline); j++) {  // get first non-blank position and number or string
					if (tline[j] != ' ' && tline[j] != '\t' && tline[j] != '\n' && tline[j] != '\r') {
						pos = j; break;
					}
				}
				sscanf(tline, "%s", elemnumstr[i]);
				for (j = 0; j < strlen(elemnumstr[i]); j++) {
					tline[pos + j] = ' ';
				}
				typenum[i] = atoi(elemnumstr[i]);
				if (typenum[i] <= 0) {
					printf("Error: found %i atom for ion type %i in CONTCAR (Line %i).\n",
						typenum[i], i + 1, iline);
					fclose(fp); exit(1);
				}
				nant[0] += typenum[i];  //Total number of atoms
			}

			continue;
		}

		// In case of Selective Dynamic
		if (iline == 8) { // In normal case, line 9 is atom position 
			sscanf(tline, "%s", str);
			if (str[0] == 'S' || str[0] == 's')
			{
				flag = 0; ifix = 0; iflg = 0; //Selective dynamics
			}
			else if (str[0] == 'D' || str[0] == 'd')
			{
				flag = 1; ifix = 1; iflg = 0;
			}
			else if (str[0] == 'C' || str[0] == 'c')
			{
				flag = 1; ifix = 1; iflg = 1;
			}
			else printf("Error with ifix and iflg!");

			continue;
		}

		// read coordinates or atomic postions 
		if (flag == 0) {         // In selective case, line 10 is atomic position 
			//if( iline == 9 ){  // In selective case, line 9  is coordinate system
			sscanf(tline, "%s", str);
			if (str[0] == 'D' || str[0] == 'd')
			{
				iflg = 0;
			}
			else if (str[0] == 'C' || str[0] == 'c')
			{
				iflg = 1;
			}
			else printf("Error with iflg!");

			flag = 1;     //Move to the next line 10 to read atomic position
			continue;
			//}
		}
		else if (flag == 1) {  // atomic positons from line 9 or line 10                   
			if (ifix == 0)
			{
				sscanf(tline, "%lf %lf %lf %c %c %c", &x, &y, &z, &fx, &fy, &fz);
				xyz[iat][0] = x; xyz[iat][1] = y; xyz[iat][2] = z;
				fix[iat][0] = fx; fix[iat][1] = fy; fix[iat][2] = fz;
			}
			else
			{
				sscanf(tline, "%lf %lf %lf ", &x, &y, &z);
				xyz[iat][0] = x; xyz[iat][1] = y; xyz[iat][2] = z;
				fix[iat][0] = 'T'; fix[iat][1] = 'T'; fix[iat][2] = 'T';
			}

			if (++iat >= nant[0])  break;

			continue;
		}
		else {
			printf("Error with flag!");
			continue;  //go to next line
		}
	} // End of while

/* To transform elemental symbol to atomic number */
	for (i = 0; i < nant[1]; i++)
	{
		int an;
		char at[3];
		strcpy(at, elemsym[i]);
		getAtomNum(at, an);
		elemnum[i] = an;
	}
	// default direct coordinates
	if (iflg == 1)
		carts_to_direct(iflg, vec, xyz, nant);
	return 0;
}

int savposcar(FILE* fp, POSCAR poscar)
{
	//	printf(" >>> Saving to ASCII file\n");

	int  i, j;
	//	printf(" poscar.title= %s created by spl2gsg.x \n",poscar.title);
	fprintf(fp, "%s \n", poscar.title);

	//	printf("%16.9lf \n", poscar.latt);
	fprintf(fp, "%16.9lf \n", poscar.latt);
	fprintf(fp, "  %16.9lf %16.9lf %16.9lf \n  %16.9lf %16.9lf %16.9lf \n  %16.9lf %16.9lf %16.9lf \n",
		poscar.vec[0][0], poscar.vec[0][1], poscar.vec[0][2], poscar.vec[1][0], poscar.vec[1][1], poscar.vec[1][2], poscar.vec[2][0], poscar.vec[2][1], poscar.vec[2][2]);

	for (i = 0; i < poscar.nant[1]; i++)
	{
		//		printf("    %s  ", poscar.elemsym[i]);
		fprintf(fp, "    %s  ", poscar.elemsym[i]);
	}
	fprintf(fp, "\n");

	for (i = 0; i < poscar.nant[1]; i++)
	{
		//		printf(" %5d ", poscar.typenum[i]);   
		fprintf(fp, " %5d ", poscar.typenum[i]);
	}
	fprintf(fp, "\n");

	if (poscar.ifix == 0)
	{
		fprintf(fp, "Select \n");
	}

	if (poscar.iflg == 0)
	{
		fprintf(fp, "Direct \n");
	}
	else if (poscar.iflg == 1)
	{
		fprintf(fp, "Cartes \n");
	}
	else
	{
		fprintf(fp, " Wrong with flag read from bin2pos.x!! \n");
	}

	for (j = 0; j < poscar.nant[0]; j++)
	{
		if (poscar.ifix == 0)
		{
			//	printf("  %16.9lf %16.9lf %16.9lf %c %c %c \n",
			//		poscar.xyz[j][0],poscar.xyz[j][1],poscar.xyz[j][2],poscar.fix[j][0],poscar.fix[j][1],poscar.fix[j][2]);
			fprintf(fp, "  %16.9lf %16.9lf %16.9lf %c %c %c \n",
				poscar.xyz[j][0], poscar.xyz[j][1], poscar.xyz[j][2], poscar.fix[j][0], poscar.fix[j][1], poscar.fix[j][2]);
		}
		else
		{
			//printf("  %16.9lf %16.9lf %16.9lf \n",poscar.xyz[j][0],poscar.xyz[j][1],poscar.xyz[j][2]);
			fprintf(fp, "  %16.9lf %16.9lf %16.9lf \n", poscar.xyz[j][0], poscar.xyz[j][1], poscar.xyz[j][2]);
		}
	}

	return 0;
}

int readposcar(FILE* fp, POSCAR& poscar)
{
	//	printf(" >>> Reading ASCII file\n");
	int    iline;
	int    flag;
	char   tline[MAX_NCOLU];
	double sf;                                     // scaling factor
	int    i, j;
	int    pos;
	char   elemnumstr[MAX_NELEM][5];
	char   str[10];
	int    iat;                                    // index of atoms
	double x, y, z;
	char   fx, fy, fz;

	iline = 0;
	poscar.nant[0] = 0;
	poscar.nant[1] = 0;
	flag = 2;
	poscar.ifix = 0;
	poscar.iflg = 0;
	iat = 0;
	while (1) {                                       //loop each line to the end
		fgets(tline, MAX_NCOLU, fp);             //Read one line
		iline++;                                  //move to next line 
		if (feof(fp) != 0) break;                 //Check end of the file: yes exit loop

		// read poscar.title from line #1
		if (iline == 1) {
			for (int i = 0; i < MAX_FNAME; i++)
			{
				if (tline[i] == '\n') { break; }
				poscar.title[i] = tline[i];
			}
			//				printf(" poscar.title= %s\n", poscar.title);
			continue;
		}

		// read scaling factor from line #2
		if (iline == 2) {
			sscanf(tline, "%lf", &sf);
			poscar.latt = sf;
#if DEBUG
			//					printf("%16.9lf \n", poscar.latt);
#endif
			continue;
		}

		// read poscar.lattice vectors a from line #3
		if (iline == 3) {
			sscanf(tline, "%lf%lf%lf", &poscar.vec[0][0], &poscar.vec[0][1], &poscar.vec[0][2]);
#if DEBUG
			//					printf("  %16.9lf %16.9lf %16.9lf \n", poscar.vec[0][0],poscar.vec[0][1],poscar.vec[0][2]);
#endif
			continue;
		}

		// read poscar.lattice vectors b from line #4
		if (iline == 4) {
			sscanf(tline, "%lf%lf%lf", &poscar.vec[1][0], &poscar.vec[1][1], &poscar.vec[1][2]);
			//	    		printf("  %16.9lf %16.9lf %16.9lf \n", poscar.vec[1][0],poscar.vec[1][1],poscar.vec[1][2]);
			continue;
		}

		// read poscar.lattice vectors c from line #5
		if (iline == 5) {
			sscanf(tline, "%lf%lf%lf", &poscar.vec[2][0], &poscar.vec[2][1], &poscar.vec[2][2]);
			//		printf("  %16.9lf %16.9lf %16.9lf \n", poscar.vec[2][0],poscar.vec[2][1],poscar.vec[2][2]);
			continue;
		}

		// read elemental symbols from line #6
		if (iline == 6) {
			for (i = 0; i < MAX_NELEM; i++) {
				for (j = 0; j < strlen(tline); j++) { // get first non-blank position and string
					if (tline[j] != ' ' && tline[j] != '\t' && tline[j] != '\n' && tline[j] != '\r') {
						pos = j; break;      //exit the innermost enclosing loop
					}
					else if (tline[j] == '\n' || tline[j] == '\r') {
						pos = 9999; break;   //exit the innermost enclosing loop
					}
				}

				if (pos == 9999) break;
				sscanf(tline, "%s", poscar.elemsym[i]);
				for (j = 0; j < strlen(poscar.elemsym[i]); j++) {
					tline[pos + j] = ' ';
				}
			}
			poscar.nant[1] = i;               //Total number of elemental types

			for (i = 0; i < poscar.nant[1]; i++)
			{
				//    		printf("    %s  ", poscar.elemsym[i]);
			}
			//		printf("\n");

			continue;
		}

		// read number of atoms for each ion type
		if (iline == 7) {
			for (i = 0; i < poscar.nant[1]; i++) {
				for (j = 0; j < strlen(tline); j++) {  // get first non-blank position and number or string
					if (tline[j] != ' ' && tline[j] != '\t' && tline[j] != '\n' && tline[j] != '\r') {
						pos = j; break;
					}
				}
				sscanf(tline, "%s", elemnumstr[i]);
				for (j = 0; j < strlen(elemnumstr[i]); j++) {
					tline[pos + j] = ' ';
				}
				poscar.typenum[i] = atoi(elemnumstr[i]);
				if (poscar.typenum[i] <= 0) {
					printf("Error: found %i atom for ion type %i in CONTCAR (Line %i).\n",
						poscar.typenum[i], i + 1, iline);
					fclose(fp); exit(1);
				}
				poscar.nant[0] += poscar.typenum[i];  //Total number of atoms
			}

			continue;
		}

		// In case of Selective Dynamic
		if (iline == 8) { // In normal case, line 9 is atom position 
			sscanf(tline, "%s", str);
			if (str[0] == 'S' || str[0] == 's')
			{
				flag = 0; poscar.ifix = 0; poscar.iflg = 0; //Selective dynamics
			}
			else if (str[0] == 'D' || str[0] == 'd')
			{
				flag = 1; poscar.ifix = 1; poscar.iflg = 0;
			}
			else if (str[0] == 'C' || str[0] == 'c')
			{
				flag = 1; poscar.ifix = 1; poscar.iflg = 1;
			}
			else printf("Error with poscar.ifix and poscar.iflg!");

			continue;
		}

		// read coordinates or atomic postions 
		if (flag == 0) {         // In selective case, line 10 is atomic position 
			//if( iline == 9 ){  // In selective case, line 9  is coordinate system
			sscanf(tline, "%s", str);
			if (str[0] == 'D' || str[0] == 'd')
			{
				poscar.iflg = 0;
			}
			else if (str[0] == 'C' || str[0] == 'c')
			{
				poscar.iflg = 1;
			}
			else printf("Error with poscar.iflg!");

			flag = 1;     //Move to the next line 10 to read atomic position
			continue;
			//}
		}
		else if (flag == 1) {  // atomic positons from line 9 or line 10                   
			if (poscar.ifix == 0)
			{
				sscanf(tline, "%lf %lf %lf %c %c %c", &x, &y, &z, &fx, &fy, &fz);
				poscar.xyz[iat][0] = x; poscar.xyz[iat][1] = y; poscar.xyz[iat][2] = z;
				poscar.fix[iat][0] = fx; poscar.fix[iat][1] = fy; poscar.fix[iat][2] = fz;
			}
			else
			{
				sscanf(tline, "%lf %lf %lf ", &x, &y, &z);
				poscar.xyz[iat][0] = x; poscar.xyz[iat][1] = y; poscar.xyz[iat][2] = z;
				poscar.fix[iat][0] = 'T'; poscar.fix[iat][1] = 'T'; poscar.fix[iat][2] = 'T';
			}

			if (++iat >= poscar.nant[0])  break;

			continue;
		}
		else {
			printf("Error with flag!");
			continue;  //go to next line
		}
	} // End of while

/* To transform elemental symbol to atomic number */
	for (i = 0; i < poscar.nant[1]; i++)
	{
		int an;
		char at[3];
		strcpy(at, poscar.elemsym[i]);
		getAtomNum(at, an);
		poscar.elemnum[i] = an;
	}
	// default direct coordinates
	if (poscar.iflg == 1)
		carts_to_direct(poscar.iflg, poscar.vec, poscar.xyz, poscar.nant);
	return 0;
}

int readposbin(FILE* fp,
	char   title[],                   // title of the system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                      // nant[0],nant[1] number of atoms and types
	int    typenum[],                 // number of atoms for each ion type
	int    elemnum[],                 // atomic number of each element
	char   elemsym[][3],              // atomic symbol of each element
	double vec[3][3],                 // lattice vectors
	double xyz[][3],                  // atomic coordinates
	char   fix[][3]                   // selective fix for atoms
) {
	//	printf(" >>> Reading BINARY file\n");
	int i, j;

	VEC pVEC[1];
	TYP pTYP[MAX_NELEM];
	POS pPOS[MAX_NATOM];

	int num_read = 0;
	num_read = fread(&pVEC, sizeof(VEC), 1, fp);
	if (0 == num_read)
	{
		perror("read file failed.\n");
		return 1;

	}
	else
	{
		strcpy(title, pVEC[0].st);
		//                printf(" Title= %s\n", title);

		nant[0] = pVEC[0].na;
		nant[1] = pVEC[0].nt;
		ifix = pVEC[0].ifix;
		iflg = pVEC[0].iflg;
		//            printf("Natm=%5d Ntyp=%5d ifix=%5d iflg=%5d\n", nant[0],nant[1],ifix,iflg);

		latt = pVEC[0].abc;
		//				printf("%16.9lf \n", latt);

		vec[0][0] = pVEC[0].ax; vec[0][1] = pVEC[0].ay; vec[0][2] = pVEC[0].az;
		vec[1][0] = pVEC[0].bx; vec[1][1] = pVEC[0].by; vec[1][2] = pVEC[0].bz;
		vec[2][0] = pVEC[0].cx; vec[2][1] = pVEC[0].cy; vec[2][2] = pVEC[0].cz;
		//				printf("  %16.9lf %16.9lf %16.9lf \n", vec[0][0],vec[0][1],vec[0][2]);
		//				printf("  %16.9lf %16.9lf %16.9lf \n", vec[1][0],vec[1][1],vec[1][2]);
		//				printf("  %16.9lf %16.9lf %16.9lf \n", vec[2][0],vec[2][1],vec[2][2]);
	}

	fread(&pTYP, sizeof(TYP), nant[1], fp);
	for (i = 0; i < nant[1]; i++)
	{
		int  an;
		char at[3];
		elemnum[i] = pTYP[i].id;
		an = elemnum[i];
		getAtomSym(an, at);
		strcpy(elemsym[i], at);
		//               printf("%5d %s", elemnum[i],elemsym[i]);
	}
	//      printf("\n");

	for (i = 0; i < nant[1]; i++)
	{
		typenum[i] = pTYP[i].nt;
		//             printf("%5d ", typenum[i]);
	}
	//    printf("\n");

	fread(&pPOS, sizeof(POS), nant[0], fp);
	for (j = 0; j < nant[0]; j++)
	{
		if (ifix == 0)
		{
			xyz[j][0] = pPOS[j].px; xyz[j][1] = pPOS[j].py; xyz[j][2] = pPOS[j].pz;
			fix[j][0] = pPOS[j].fx; fix[j][1] = pPOS[j].fy; fix[j][2] = pPOS[j].fz;
			//		printf("  %16.9lf %16.9lf %16.9lf %c %c %c \n", 
			//			xyz[j][0],xyz[j][1],xyz[j][2],fix[j][0],fix[j][1],fix[j][2]);
		}
		else
		{
			xyz[j][0] = pPOS[j].px; xyz[j][1] = pPOS[j].py; xyz[j][2] = pPOS[j].pz;
			fix[j][0] = 'T'; fix[j][1] = 'T'; fix[j][2] = 'T';
			//		printf("  %16.9lf %16.9lf %16.9lf \n", xyz[j][0],xyz[j][1],xyz[j][2]);
		}
	}

	return 0;
}



int savposbin(FILE* fp,
	char   title[],                   // title of the system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                    // nant[0],nant[1] number of atoms and types
	int    typenum[],                 // number of atoms for each ion type
	int    elemnum[],                 // atomic number of each element
	char   elemsym[][3],              // atomic symbol of each element
	double vec[3][3],                 // lattice vectors
	double xyz[][3],                  // atomic coordinates
	char   fix[][3]                   // Selective fix for each atom
) {
	//	printf(" >>> Saving to BINARY file\n");
	int i, j;

	VEC pVEC[1];
	TYP pTYP[MAX_NELEM];
	POS pPOS[MAX_NATOM];

	strcpy(pVEC[0].st, title);
	//   printf(" Title= %s\n", pVEC[0].st);

	pVEC[0].na = nant[0];
	pVEC[0].nt = nant[1];
	pVEC[0].ifix = ifix;
	pVEC[0].iflg = iflg;
	//  printf(" Natm=%5d Ntyp=%5d ifix=%5d iflg=%5d\n", 
  //	pVEC[0].na,pVEC[0].nt,pVEC[0].ifix,pVEC[0].iflg);

	pVEC[0].abc = latt;
	//	printf(" %16.9lf \n", pVEC[0].abc);

	pVEC[0].ax = vec[0][0]; pVEC[0].ay = vec[0][1]; pVEC[0].az = vec[0][2];
	pVEC[0].bx = vec[1][0]; pVEC[0].by = vec[1][1]; pVEC[0].bz = vec[1][2];
	pVEC[0].cx = vec[2][0]; pVEC[0].cy = vec[2][1]; pVEC[0].cz = vec[2][2];
	//	printf("  %16.9lf %16.9lf %16.9lf \n", pVEC[0].ax,pVEC[0].ay,pVEC[0].az);
	//	printf("  %16.9lf %16.9lf %16.9lf \n", pVEC[0].bx,pVEC[0].by,pVEC[0].bz);
	//	printf("  %16.9lf %16.9lf %16.9lf \n", pVEC[0].cx,pVEC[0].cy,pVEC[0].cz);

	fwrite(&pVEC, sizeof(VEC), 1, fp);

	for (i = 0; i < pVEC[0].nt; i++)
	{
		int an;
		char at[3];
		strcpy(at, elemsym[i]);
		getAtomNum(at, an);
		elemnum[i] = an;
		pTYP[i].id = elemnum[i];
		//               printf("%5d %s", pTYP[i].id, elemsym[i]);
	}
	//      printf("\n");


	for (i = 0; i < pVEC[0].nt; i++)
	{
		pTYP[i].nt = typenum[i];
		//             printf("%5d ", pTYP[i].nt);
	}
	//    printf("\n");

	fwrite(&pTYP, pVEC[0].nt * sizeof(TYP), 1, fp);

	for (j = 0; j < pVEC[0].na; j++)
	{
		if (ifix == 0)
		{
			pPOS[j].px = xyz[j][0];
			pPOS[j].py = xyz[j][1];
			pPOS[j].pz = xyz[j][2];
			pPOS[j].fx = fix[j][0];
			pPOS[j].fy = fix[j][1];
			pPOS[j].fz = fix[j][2];
			//                   printf("  %16.9lf %16.9lf %16.9lf %c %c %c \n", 
		   //			pPOS[j].px, pPOS[j].py, pPOS[j].pz, pPOS[j].fx, pPOS[j].fy, pPOS[j].fz);

		}
		else
		{
			pPOS[j].px = xyz[j][0];
			pPOS[j].py = xyz[j][1];
			pPOS[j].pz = xyz[j][2];
			//           	printf("  %16.9lf %16.9lf %16.9lf \n", 
		   //			pPOS[j].px, pPOS[j].py, pPOS[j].pz);
		}
	}

	fwrite(&pPOS, pVEC[0].na * sizeof(POS), 1, fp);

	return 0;
}







