#include"../../include/Stochastics_include/read_write.h"
#include"../../include/Stochastics_include/variant.h"
#include"../../include/Stochastics_include/tools.h"
#include<algorithm>
using namespace std;
void copy_file_sto(const  char* target, const  char* source)
{
	FILE* fp1 = fopen(source, "r");
	if (fp1 == NULL)
	{
		printf("%s  cannot open \n", source);
		return;
	}
	FILE* fp2 = fopen(target, "at+");
	char buff[1024] = " ";
	while (fgets(buff, 1024, fp1) != NULL)
		fputs(buff, fp2);
	fclose(fp1);
	fclose(fp2);
}

double  RandDigit()                                 // creating the random number from 0.0 to 1.0)
{
	double temp;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	srand(tv.tv_usec);                       // time seed for random digit
	temp = rand() % 10000 / 10000.1;
	return temp;
}

// subroutine "DistAtom": whether the interatomic distance of any two atoms >= 1.0000 or not
int DistAtom(
	int    nant[],                    // total number of atoms and types
	double& latt_o,			  // Scaling Lattice parameter
	double vec_o[3][3],               // lattice vectors
	double xyz_o[][3]                 // atomic coordinates
) {
	int i, j, k;
	int logic = 1;
	double d = 0.5 * pow(volume(vec_o) / nant[0], 1.0 / 3.0) * latt_o;
	double dist;                              // the distance between two atoms
	double xyz[MAX_NATOM][3];
	for (i = 0; i < nant[0]; i++) {
		for (j = 0; j < 3; j++) {
			xyz[i][j] = latt_o * (xyz_o[i][0] * vec_o[0][j] + xyz_o[i][1] * vec_o[1][j] + xyz_o[i][2] * vec_o[2][j]);
		}
	}                                         // "Dirct" >>> "Cartesian"
	for (i = 0; i < nant[0]; i++) {
		for (j = i + 1; j < nant[0]; j++) {
			dist = 0.0;
			for (k = 0; k < 3; k++) {
				dist = dist + (xyz[i][k] - xyz[j][k]) * (xyz[i][k] - xyz[j][k]);
			}
			dist = sqrt(dist);
			if (dist < d) 
				logic = 0; return logic;
		}
	}
	return logic;
}

void Exchange(
	int    nant_o[],                  // total number of atoms and types
	int    typenum_o[],               // number of each type 
	double xyz_o[][3],                // atomic coordinates
	double xyz_e[][3]
) {
	// Exchange operator which exchanges the coordinates of two atoms of different types
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>>Exchange operator \n");
	fprintf(fp_buf, ">>>The atoms selected for exchange are\n");
	int NumExg = min(5, nant_o[0] / 2);                             // Number of exchanges
	int Num1, Num2;
	double temp[3];
	int Nlower, Nuper;
	for (int i = 0; i < NumExg; i++)
	{
	angin:
		Nlower = 0;
		Nuper = typenum_o[0];
		Num1 = floor(nant_o[0] * RandDigit());
		Num2 = floor(nant_o[0] * RandDigit());

		for (int j = 0; j < nant_o[1]; j++)
		{
			if (Num1 >= Nlower && Num1 < Nuper &&
				Num2 >= Nlower && Num2 < Nuper && nant_o[1] != 1)
			{
				goto angin;
			}
			Nlower = Nlower + typenum_o[j];
			Nuper = Nuper + typenum_o[j + 1];
		}                                     // two atoms of different types

		for (int j = 0; j < 3; j++) {
			temp[j] = xyz_o[Num1][j];
			xyz_o[Num1][j] = xyz_o[Num2][j];
			xyz_o[Num2][j] = temp[j];
		}
		fprintf(fp_buf, "%d  %d\n", Num1, Num2);
	}
	for (int i = 0; i < nant_o[0]; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			xyz_e[i][j] = xyz_o[i][j];
		}
	}
	fclose(fp_buf);
}

void Strain(
	char   mode[],
	double vec_o[3][3],               // lattice vectors
	double vec_e[3][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>>Strain operator \n");
	int i, j;
	double strain[6];
	double defvect[3][3] = { {1.0000, 0.0000, 0.0000},{0.0000, 1.0000, 0.0000},{0.0000, 0.0000, 1.0000} };
	for (i = 0; i < 6; i++) {
		strain[i] = RandDigit() / 2.0000;         // the maximum strain: 0.5000 we should use normal 
	}
	// tensile strain along x axis
	if (!strcmp(mode, "xx"))
	{
		fprintf(fp_buf, ">>>the values of tensile strain along xx axis is \n%1.4f \n", strain[0]);
		defvect[0][0] = 1.0000 + strain[0];
	}
	// tensile strain along y axis
	if (!strcmp(mode, "yy"))
	{
		fprintf(fp_buf, ">>>the values of tensile strain along yy axis is \n%1.4f \n", strain[1]);
		defvect[1][1] = 1.0000 + strain[1];
	}
	// tensile strain along z axis
	if (!strcmp(mode, "zz"))
	{
		fprintf(fp_buf, ">>>the values of tensile strain along zz axis is \n%1.4f \n", strain[2]);
		defvect[2][2] = 1.0000 + strain[2];
	}
	// shear strain in x-y plane
	if (!strcmp(mode, "xy"))
	{
		fprintf(fp_buf, ">>>the values of the shear strain in xy plane is \n%1.4f \n", strain[3]);
		defvect[0][1] = strain[3] / 2.0000;
		defvect[1][0] = strain[3] / 2.0000;
	}
	// shear strain in x-z plane
	if (!strcmp(mode, "xz"))
	{
		fprintf(fp_buf, ">>>the values of the shear strain in xz plane is \n%1.4f \n", strain[4]);
		defvect[0][2] = strain[4] / 2.0000;
		defvect[2][0] = strain[4] / 2.0000;
	}
	// shear strain in y-z plane
	if (!strcmp(mode, "yz"))
	{
		fprintf(fp_buf, ">>>the values of the shear strain in yz plane is \n%1.4f \n", strain[5]);
		defvect[1][2] = strain[5] / 2.0000;
		defvect[2][1] = strain[5] / 2.0000;
	}
	// strain
	if (!strcmp(mode, "xyz"))
	{
		fprintf(fp_buf, ">>>the values of the pure strain is that: \n");
		fprintf(fp_buf, "%1.4f  %1.4f  %1.4f  %1.4f  %1.4f  %1.4f \n",
			strain[0], strain[1], strain[2], strain[3], strain[4], strain[5]);
		defvect[0][0] = 1.0000 + strain[0]; defvect[0][1] = strain[3] / 2.0000; defvect[0][2] = strain[4] / 2.0000;
		defvect[1][0] = strain[3] / 2.0000; defvect[1][1] = 1.0000 + strain[1]; defvect[1][2] = strain[5] / 2.0000;
		defvect[2][0] = strain[4] / 2.0000; defvect[2][1] = strain[5] / 2.0000; defvect[2][2] = 1.0000 + strain[2];
	}
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			vec_e[i][j] = vec_o[i][0] * defvect[0][j] + vec_o[i][1] * defvect[1][j] + vec_o[i][2] * defvect[2][j];
		}
	}
	fclose(fp_buf);
}

void Ripple(
	int    nant_o[],                  // total number of atoms and types
	double xyz_o[][3],                // atomic coordinates
	double xyz_e[][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	// Ripple opetator which shifts the coordinates of each atom long a random axis
	fprintf(fp_buf, ">>> Ripple operator \n");
	int i;
	int temp;
	char axis;
	int x, y, z;
	double MaxDisp = 0.5000;                     // the maximum possible displacement of any atom in the displaced(z) direction
	double theta_1, theta_2;                   // randomly chosen from the interval 0~2pi to vary the regions of the structure that are strongly displaced
	double angletemp_1, angletemp_2;
	temp = floor(6 * RandDigit());
	switch (temp)
	{
	case 0: x = 0, y = 1, z = 2; axis = 'z'; break;
	case 1: x = 0, y = 2, z = 1; axis = 'y'; break;
	case 2: x = 1, y = 0, z = 2; axis = 'z'; break;
	case 3: x = 1, y = 2, z = 0; axis = 'x'; break;
	case 4: x = 2, y = 0, z = 1; axis = 'y'; break;
	case 5: x = 2, y = 1, z = 0; axis = 'x'; break;
	}                                          // choosing a random axis
	fprintf(fp_buf, ">>>The ration to control axis and angle is temp = \n%d\n", temp);
	//fprintf(fp_buf, ">>>The displaced(z) direction is \n%c \n", axis);
	theta_1 = 2 * Pi * RandDigit();
	theta_2 = 2 * Pi * RandDigit();
	angletemp_1 = theta_1 / Pi * 180.00;
	angletemp_2 = theta_2 / Pi * 180.00;
	fprintf(fp_buf, ">>>angletemp_x, angletemp_y=\n%1.4f  %1.4f \n", angletemp_1, angletemp_2);
	for (i = 0; i < nant_o[0]; i++) {
		xyz_e[i][x] = xyz_o[i][x];
		xyz_e[i][y] = xyz_o[i][y];
		xyz_e[i][z] = xyz_o[i][z] + MaxDisp * cos(2 * Pi * xyz_o[i][x] + theta_1) * cos(2 * Pi * xyz_o[i][y] + theta_2);
		if (xyz_e[i][z] >= 1) { --xyz_e[i][z]; }
		if (xyz_e[i][z] < 0) { ++xyz_e[i][z]; }
	}
	fclose(fp_buf);
}

//The random move operator moves all atoms’ positions randomly, keeping the lattice vectors unchanged
void Randmove(
	int    nant_o[],
	double xyz_o[][3],
	double xyz_e[][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>> Random move perator \n >>>The distance each atom then moves is \n");
	int i, j;
	double dxyz;
	for (i = 0; i < nant_o[0]; i++)
	{
		for (j = 0; j < 3; j++)
		{
			dxyz = 2 * RandDigit() - 1;            // -1.0~1.0
			xyz_e[i][j] = xyz_o[i][j] + dxyz;
			if (xyz_e[i][j] >= 1) { --xyz_e[i][j]; }
			if (xyz_e[i][j] < 0) { ++xyz_e[i][j]; }
			fprintf(fp_buf, "%lf  ", dxyz);
		}
		fprintf(fp_buf, "\n");
	}
	fclose(fp_buf);
}

//This operator slips a group of atoms by a random distance along a random direction that is parallel to a specific plane
void Slip(
	int    nant_o[],                  // total number of atoms and types
	double xyz_o[][3],                // atomic coordinates
	double xyz_e[][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>> Slip operator \n");
	int i;
	double xlower = 0.25 + 0.5 * RandDigit();
	double ylower = 0.25 + 0.5 * RandDigit();
	double zlower = 0.25 + 0.5 * RandDigit();        // 0.25~0.75
	double xdisp = RandDigit();
	double ydisp = RandDigit();
	double zdisp = RandDigit();
	fprintf(fp_buf, ">>>slip on y-z plane: zlower,slip on z-x plane: xlower,slip on x-y plane: ylower=\n%1.2f %1.2f %1.2f\n", zlower, xlower, ylower);
	fprintf(fp_buf, ">>>displace along y axis,displace along z axis,displace along x axis=\n%1.4f %1.4f %1.4f\n", ydisp, zdisp, xdisp);
	fclose(fp_buf);
	/*printf(" slip on y-z plane: xlower=%1.2f, displace along y axis=%1.4f \n", xlower, ydisp);
	printf(" slip on z-x plane: ylower=%1.2f, displace along z axis=%1.4f \n", ylower, zdisp);
	printf(" slip on x-y plane: zlower=%1.2f, displace along x axis=%1.4f \n", zlower, xdisp);*/
	for (i = 0; i < nant_o[0]; i++) {
		// slip on y-z plane
		if (xyz_o[i][2] > zlower)
		{
			xyz_e[i][1] = xyz_o[i][1] + ydisp;
			if (xyz_e[i][1] > 1) { --xyz_e[i][1]; }//
		}
		else
		{
			xyz_e[i][1] = xyz_o[i][1];
		}

		// slip on z-x plane
		if (xyz_o[i][0] > xlower)
		{
			xyz_e[i][2] = xyz_o[i][2] + zdisp;
			if (xyz_e[i][2] > 1) { --xyz_e[i][2]; }
		}
		else
		{
			xyz_e[i][2] = xyz_o[i][2];
		}

		// slip on x-y plane
		if (xyz_o[i][1] > ylower)
		{
			xyz_e[i][0] = xyz_o[i][0] + xdisp;
			if (xyz_e[i][0] > 1) { --xyz_e[i][0]; }
		}
		else
		{
			xyz_e[i][0] = xyz_o[i][0];
		}
	}
}

void Twist(
	int    nant_o[],                  // total number of atoms and types
	double xyz_o[][3],                // atomic coordinates
	double xyz_e[][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>> Twist operator \n");
	int i, j, zz;
	char axis;
	double angletemp;
	double angle[3];
	zz = floor(3 * RandDigit());
	switch (zz)
	{
	case 0: axis = 'x'; break;
	case 1: axis = 'y'; break;
	case 2: axis = 'z'; break;
	}
	for (i = 0; i < 3; i++) {
		angle[i] = 2 * Pi * RandDigit();
	}                                          // select a random axis and a random angle
	angletemp = angle[zz] / Pi * 180.00;
	fprintf(fp_buf, ">>>Selecting a randomly chosen axis: \n%c\n >>>The rotation angle is \n%3.2f \n", axis, angletemp);
	fclose(fp_buf);
	// twist along x axis
	if (zz == 0)
	{
		for (i = 0; i < nant_o[0]; i++) {
			xyz_e[i][0] = xyz_o[i][0];
			xyz_e[i][1] = xyz_o[i][1] * cos(angle[0] * xyz_o[i][0]) + xyz_o[i][2] * sin(angle[0] * xyz_o[i][0]);
			xyz_e[i][2] = -xyz_o[i][1] * sin(angle[0] * xyz_o[i][0]) + xyz_o[i][2] * cos(angle[0] * xyz_o[i][0]);
		}
	}
	// twist along y axis
	if (zz == 1)
	{
		for (i = 0; i < nant_o[0]; i++) {
			xyz_e[i][0] = xyz_o[i][0] * cos(angle[1] * xyz_o[i][1]) - xyz_o[i][2] * sin(angle[1] * xyz_o[i][1]);
			xyz_e[i][1] = xyz_o[i][1];
			xyz_e[i][2] = xyz_o[i][0] * sin(angle[1] * xyz_o[i][1]) + xyz_o[i][2] * cos(angle[1] * xyz_o[i][1]);
		}
	}
	// twsit along z axis
	if (zz == 2)
	{
		for (i = 0; i < nant_o[0]; i++) {
			xyz_e[i][0] = xyz_o[i][0] * cos(angle[2] * xyz_o[i][2]) + xyz_o[i][1] * sin(angle[2] * xyz_o[i][2]);
			xyz_e[i][1] = -xyz_o[i][0] * sin(angle[2] * xyz_o[i][2]) + xyz_o[i][1] * cos(angle[2] * xyz_o[i][2]);
			xyz_e[i][2] = xyz_o[i][2];
		}
	}
	for (i = 0; i < nant_o[0]; i++) {
		for (j = 0; j < 3; j++) {
			if (xyz_e[i][j] >= 1) { --xyz_e[i][j]; }
			if (xyz_e[i][j] < 0) { ++xyz_e[i][j]; }
		}
	}
}

void Rotation(
	double vec_o[3][3],     // lattice vectors
	double vec_e[3][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>> Rotation operator \n");
	int i;
	double angletemp[3];
	double angle[3];
	double vec_t[3][3];
	double vec_m[3][3];
	for (i = 0; i < 3; i++) {
		angle[i] = 2 * Pi * RandDigit();
		angletemp[i] = angle[i] / Pi * 180.00;
	}
	fprintf(fp_buf, ">>> the rotation angles in x, y and z axis are \n%3.2f %3.2f %3.2f\n", angletemp[0], angletemp[1], angletemp[2]);
	//fprintf(fp_buf, ">>> the rotation angles in x, y and z axis are %3.2f, %3.2f and %3.2f , respectively \n", angle[0], angle[1], angle[2]);
	fclose(fp_buf);
	//printf("the rotation angles in x, y and z axis are %3.2f, %3.2f and %3.2f degree, respectively \n",angletemp[0], angletemp[1], angletemp[2]);
	// rotate along x axis
	for (i = 0; i < 3; i++) {
		vec_t[i][0] = vec_o[i][0];
		vec_t[i][1] = vec_o[i][1] * cos(angle[0]) + vec_o[i][2] * sin(angle[0]);
		vec_t[i][2] = -vec_o[i][1] * sin(angle[0]) + vec_o[i][2] * cos(angle[0]);
	}
	// rotate along y axis
	for (i = 0; i < 3; i++) {
		vec_m[i][0] = vec_t[i][0] * cos(angle[1]) - vec_t[i][2] * sin(angle[1]);
		vec_m[i][1] = vec_t[i][1];
		vec_m[i][2] = vec_t[i][0] * sin(angle[1]) + vec_t[i][2] * cos(angle[1]);
	}
	// rotate along z axis
	for (i = 0; i < 3; i++) {
		vec_e[i][0] = vec_m[i][0] * cos(angle[2]) + vec_m[i][1] * sin(angle[2]);
		vec_e[i][1] = -vec_m[i][0] * sin(angle[2]) + vec_m[i][1] * cos(angle[2]);
		vec_e[i][2] = vec_m[i][2];
	}
}

void Crossover(double& latt1_o,
	double& latt2_o,                  // Scaling Lattice parameter 
	double& latt_e,
	int    nant_o[],                    // total number of atoms and types
	int    typenum_o[],               // number of each type
	double vec1_o[3][3],
	double vec2_o[3][3],              // lattice vectors
	double vec_e[3][3],
	double xyz1_o[][3],
	double xyz2_o[][3],               // atomic coordinates
	double xyz_e[][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>>Crossover operator \n");
	int i, j, k;
	latt_e = 1.0;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			vec_e[i][j] = (latt1_o * vec1_o[i][j] + latt2_o * vec2_o[i][j]) / 2.0000;
		}
	}
	double segment = 0.25 + 0.5 * RandDigit();       // 0.25~0.75
	fprintf(fp_buf, "POSCAR_1: %1.2f; POSCAR_2: %1.2f \n", segment, 1.0000 - segment);
	int countnum;
	double temp[1000][3];
	int zz;
	char axis;
	int Nlower = 0;
	int Nuper = typenum_o[0];
	zz = floor(3 * RandDigit());
	switch (zz)
	{
	case 0: axis = 'x'; break;
	case 1: axis = 'y'; break;
	case 2: axis = 'z'; break;
	}
	fprintf(fp_buf, ">>>Selecting a random axis:\n%c \n", axis); // select a random axis
	fprintf(fp_buf, ">>>the add atom coordinate are\n");
	for (i = 0; i < nant_o[1]; i++) {
		countnum = 0;
		for (j = Nlower; j < Nuper; j++) {
			if (xyz1_o[j][zz] < segment)
			{
				for (k = 0; k < 3; k++) { temp[countnum][k] = xyz1_o[j][k]; }
				countnum++;
			}
		}
		for (j = Nlower; j < Nuper; j++) {
			if (xyz2_o[j][zz] >= segment)
			{
				for (k = 0; k < 3; k++) { temp[countnum][k] = xyz2_o[j][k]; }
				countnum++;
			}
		}
		for (countnum; countnum < typenum_o[i]; countnum++) {
			for (j = 0; j < 3; j++) {
				temp[countnum][j] = RandDigit();
				fprintf(fp_buf, "%lf ", temp[countnum][j]);
			}
			fprintf(fp_buf, "\n");
		}
		for (j = Nlower; j < Nuper; j++) {
			for (k = 0; k < 3; k++) {
				xyz_e[j][k] = temp[j - Nlower][k];
			}
		}
		Nlower = Nlower + typenum_o[i];
		Nuper = Nuper + typenum_o[i + 1];
	}
	fclose(fp_buf);
}

void Stripple(int    nant_o[],                    // total number of atoms and types
	double vec_o[3][3],               // lattice vectors
	double vec_e[3][3],
	double xyz_o[][3],                // atomic coordinates
	double xyz_e[][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>> Stripple operator which merges strain with ripple \n");
	fclose(fp_buf);
	char mode[] = { "xyz" };
	// Strain operator
	Strain(mode, vec_o, vec_e);
	// Ripple opetator which shifts the coordinates of each atom long a random axis
	Ripple(nant_o, xyz_o, xyz_e);
}

void Permustrain(int    nant_o[],                    // total number of atoms and types
	int    typenum_o[],
	double vec_o[3][3],               // lattice vectors
	double vec_e[3][3],
	double xyz_o[][3],                // atomic coordinates
	double xyz_e[][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>> Permustrain operator which combines strain with exchange \n");
	fclose(fp_buf);
	char mode[] = { "xyz" };
	// Strain operator
	Strain(mode, vec_o, vec_e);
	// Exchange operator which exchanges the coordinates of two atoms of different types
	Exchange(nant_o, typenum_o, xyz_o, xyz_e);
}

void Rotstrain(double vec_o[3][3],               // lattice vectors
	double vec_e[3][3]
) {
	FILE* fp_buf = fopen("../buf", "at+");
	fprintf(fp_buf, ">>> Rotstrain operator which combines rotation with strain \n");
	fclose(fp_buf);
	double vec_t[3][3];
	// Rotation operator
	Rotation(vec_o, vec_t);
	int temp;
	char mode[3];
	temp = floor(3 * RandDigit());
	switch (temp)
	{
	case 0: strcpy(mode, "xy"); break;
	case 1: strcpy(mode, "xz"); break;
	case 2: strcpy(mode, "yz"); break;
	}
	// Strain operator
	Strain(mode, vec_t, vec_e);
}

void print_Variant_Operator()
{
	printf("---------------------------------------------------------------------  \n");
	printf(" Evolutionary Variant Operator  \n");
	printf(" Choose Functionals: 0- Stripple                                       \n");
	printf("                  :  1- Permustrain                                    \n");
	printf("                  :  2- Random move                                    \n");
	printf("                  :  3- Slip                                           \n");
	printf("                  :  4- Twist                                          \n");
	printf("                  :  5- Crossover                                      \n");
	printf("                  :  6- Rotstrain                                      \n");
	printf("--------------------------------------------------------------------   \n");

}

int Variant_singer_parent(char POSCAR_Parent[], char POSCAR_Child[], int method, int Generation)
{
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
	double vec_e[3][3];               // lattice vector
	double xyz_e[MAX_NATOM][3];       // atomic coordinate
	char olog[1024] = { NULL };
	sprintf(olog, "../olog_%d", Generation);
	FILE* fp_parent;
	fp_parent = fopen(POSCAR_Parent, "r");
	char buf[1024];
	if (fgets(buf, 1024, fp_parent) == NULL)
	{
		fclose(fp_parent);
		return 0;
	}
	rewind(fp_parent);
	readposcar(fp_parent, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp_parent);
	int MAX = 10000;
	if (0 == method)
	{
		Stripple(nant, vec, vec_e, xyz, xyz_e);
		double v = volume(vec);
		double v_e = volume(vec_e);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				vec_e[i][j] *= pow(v / v_e, 1.0 / 3);
		while (!DistAtom(nant, latt, vec_e, xyz_e) && MAX > 0) // to limit interatomic distance >= 1.0000 
		{
			FILE* fclear = fopen("../buf", "w");
			fclose(fclear);
			Stripple(nant, vec, vec_e, xyz, xyz_e);
			double v = volume(vec);
			double v_e = volume(vec_e);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					vec_e[i][j] *= pow(v / v_e, 1.0 / 3);
			MAX--;
		}
		copy_file_sto(olog, "../buf");
		remove("../buf");
	}
	else if (1 == method)
	{
		Permustrain(nant, typenum, vec, vec_e, xyz, xyz_e);   // to limit interatomic distance >= 1.0000 
		double v = volume(vec);
		double v_e = volume(vec_e);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				vec_e[i][j] *= pow(v / v_e, 1.0 / 3);
		while (!DistAtom(nant, latt, vec_e, xyz_e) && MAX > 0)
		{
			FILE* fclear = fopen("../buf", "w");
			fclose(fclear);
			Permustrain(nant, typenum, vec, vec_e, xyz, xyz_e);
			double v = volume(vec);
			double v_e = volume(vec_e);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					vec_e[i][j] *= pow(v / v_e, 1.0 / 3);
			MAX--;
		}
		copy_file_sto(olog, "../buf");
		remove("../buf");
	}
	else if (2 == method)
	{
		Randmove(nant, xyz, xyz_e);
		while (!DistAtom(nant, latt, vec, xyz_e) && MAX > 0)   // to limit interatomic distance >= 1.0000 
		{
			FILE* fclear = fopen("../buf", "w");
			fclose(fclear);
			Randmove(nant, xyz, xyz_e);
			double v = volume(vec);
			double v_e = volume(vec_e);
			MAX--;
		}
		copy_file_sto(olog, "../buf");
		remove("../buf");
	}

	else if (3 == method)
	{
		Slip(nant, xyz, xyz_e);
		while (!DistAtom(nant, latt, vec, xyz_e) && MAX > 0)   // to limit interatomic distance >= 1.0000 
		{
			FILE* fclear = fopen("../buf", "w");
			fclose(fclear);
			Slip(nant, xyz, xyz_e);
			MAX--;
		}
		copy_file_sto(olog, "../buf");
		remove("../buf");
	}
	else if (4 == method)
	{
		Twist(nant, xyz, xyz_e);
		while (!DistAtom(nant, latt, vec, xyz_e) && MAX > 0)   // to limit interatomic distance >= 1.0000 
		{
			FILE* fclear = fopen("../buf", "w");
			fclose(fclear);
			Twist(nant, xyz, xyz_e);
			MAX--;
		}
		copy_file_sto(olog, "../buf");
		remove("../buf");
	}
	else if (5 == method)
	{
		Rotstrain(vec, vec_e);
		double v = volume(vec);
		double v_e = volume(vec_e);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				vec_e[i][j] *= pow(v / v_e, 1.0 / 3);
		while (!DistAtom(nant, latt, vec_e, xyz) && MAX > 0)   // to limit interatomic distance >= 1.0000 
		{
			FILE* fclear = fopen("../buf", "w");
			fclose(fclear);
			Rotstrain(vec, vec_e);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					vec_e[i][j] *= v / v_e;
			MAX--;
		}
		copy_file_sto(olog, "../buf");
		remove("../buf");
	}
	/* Write POSCAR file */
	FILE* fp_child;
	fp_child = fopen(POSCAR_Child, "w");
	strcpy(title, POSCAR_Child);
	if (2 == method || 3 == method || 4 == method)
	{
		savposcar(fp_child, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz_e, fix);
	}
	else if (5 == method)
	{
		savposcar(fp_child, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz, fix);
	}
	else if (1 == method || 0 == method)
	{
		savposcar(fp_child, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz_e, fix);
	}
	fclose(fp_child);
	return 1;
}

int Variant_two_parent(char POSCAR_Parent1[], char POSCAR_Parent2[], char POSCAR_Child[], int Generation)
{
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

	char   _title[MAX_NCHAR];        // Title of system
	double _latt;                    // Scaling factor
	int    _ifix;                    // Select 0 or normal 1
	int    _iflg;                    // Direct 0 or Cartes 1
	int    _nant[2];                 // number of atoms and types
	int    _typenum[MAX_NELEM];      // number of atoms of each type
	int    _elemnum[MAX_NELEM];      // atomic number
	char   _elemsym[MAX_NELEM][3];   // atomic symbol
	double _vec[3][3];               // lattice vector
	double _xyz[MAX_NATOM][3];       // atomic coordinate
	char   _fix[MAX_NATOM][3];       // Selective fix on each atom

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

	FILE* fp_parent;
	fp_parent = fopen(POSCAR_Parent1, "r");
	FILE* _fp_parent;
	_fp_parent = fopen(POSCAR_Parent2, "r");
	char buf[1024];
	if (fgets(buf, 1024, fp_parent) == NULL || fgets(buf, 1024, _fp_parent) == NULL)
		return 1;
	rewind(fp_parent);
	rewind(_fp_parent);
	readposcar(fp_parent, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp_parent);
	readposcar(_fp_parent, _title, _latt, _ifix, _iflg, _nant, _typenum, _elemnum, _elemsym, _vec, _xyz, _fix);
	fclose(_fp_parent);
	if (nant[0] != _nant[0])return 1;
	Crossover(latt, _latt, latt_e, nant, typenum, vec, _vec, vec_e, xyz, _xyz, xyz_e);
	char olog[1024] = { NULL };
	sprintf(olog, "../olog_%d", Generation);
	int MAX = 100;
	while (!DistAtom(nant, latt_e, vec_e, xyz_e) && MAX > 0)// to limit interatomic distance >= 1.0000 
	{
		FILE* fclear = fopen("../buf", "w");
		fclose(fclear);
		Crossover(latt, _latt, latt_e, nant, typenum, vec, _vec, vec_e, xyz, _xyz, xyz_e);
		MAX--;
	}
	copy_file_sto(olog, "../buf");
	remove("../buf");
	strcpy(title_e, POSCAR_Child);
	FILE* fp_child;
	fp_child = fopen(POSCAR_Child, "w");
	savposcar(fp_child, title_e, latt_e, ifix, iflg, nant, typenum, elemnum, elemsym, vec_e, xyz_e, fix);
	fclose(fp_child);
	return 0;
}
/*
int main()
{

Variant_singer_parent ("1st_0","POSA_0",0);
Variant_singer_parent ("1st_0","POSA_1",1);
Variant_singer_parent ("1st_0","POSA_2",2);
Variant_singer_parent ("1st_0","POSA_3",3);
Variant_singer_parent ("1st_0","POSA_4",4);
Variant_singer_parent ("1st_0","POSA_5",5);
Variant_two_parent("1st_11","1st_14","POS12" );
}

*/


