#include"acneb.h"
double neb_similarity(const char name1[50], const char name2[50])
{
	//read file
	POSCAR input;
	FILE* fp1 = fopen(name1, "r");
	if (fp1 == NULL)
	{
		printf("%s is open failed\n", name1);
		exit(1);
	}
	readposcar(fp1, input);
	direct_to_carts(input.iflg, input.vec, input.xyz, input.nant);
	fclose(fp1);

	POSCAR input2;
	FILE* fp2 = fopen(name2, "r");
	if (fp2 == NULL)
	{
		printf("%s is open failed\n", name2);
		exit(1);
	}
	readposcar(fp2, input2);
	direct_to_carts(input2.iflg, input2.vec, input2.xyz, input2.nant);
	fclose(fp2);

	if (input.nant[0] != input2.nant[0])
	{
		printf("The number of atoms in %s and %s are different\n", name1, name2);
		exit(1);
	}
	double L = 0.00;
	for (int i = 0; i < input.nant[0]; i++)
	{
		double x = input2.xyz[i][0] - input.xyz[i][0];
		double y = input2.xyz[i][1] - input.xyz[i][1];
		double z = input2.xyz[i][2] - input.xyz[i][2];
		double l = x * x + y * y + z * z;
		L = l + L;
	}
	double dis = sqrt(L);
	return dis;
}

void neb_Points(const char name1[50], const char name2[50], int N)
{
	POSCAR input1;
	FILE* fp1 = fopen(name1, "r");
	if (fp1 == NULL)
	{
		printf("%s is open failed\n", name1);
		exit(1);
	}
	readposcar(fp1, input1);
	direct_to_carts(input1.iflg, input1.vec, input1.xyz, input1.nant);
	fclose(fp1);

	POSCAR input2, input3;
	FILE* fp2 = fopen(name2, "r");
	if (fp2 == NULL)
	{
		printf("%s is open failed\n", name2);
		exit(1);
	}
	readposcar(fp2, input2);
	direct_to_carts(input2.iflg, input2.vec, input2.xyz, input2.nant);
	fseek(fp2, 0, SEEK_SET);
	readposcar(fp2, input3);
	direct_to_carts(input3.iflg, input3.vec, input3.xyz, input3.nant);
	fclose(fp2);

	if (input1.nant[0] != input2.nant[0])
	{
		printf("The number of atoms in %s and %s are different\n", name1, name2);
		return;
	}

	char name[] = { 0 };
	char str[256] = { 0 };
	for (int j = 1; j <= N; j++)
	{
		for (int i = 0; i < input1.nant[0]; i++)
		{
			if (input1.xyz[i][0] == input2.xyz[i][0] && input1.xyz[i][1] == input2.xyz[i][1] && input1.xyz[i][2] == input2.xyz[i][2])
			{
				input3.xyz[i][0] = input1.xyz[i][0];
				input3.xyz[i][1] = input1.xyz[i][1];
				input3.xyz[i][2] = input1.xyz[i][2];
			}
			else if (input1.xyz[i][0] != input2.xyz[i][0] && input1.xyz[i][1] == input2.xyz[i][1] && input1.xyz[i][2] == input2.xyz[i][2])
			{
				double d1 = input2.xyz[i][0] - input1.xyz[i][0];
				double n1 = d1 / (N + 1);
				input3.xyz[i][0] = input1.xyz[i][0] + j * n1;
				input3.xyz[i][1] = input1.xyz[i][1];
				input3.xyz[i][2] = input1.xyz[i][2];
			}
			else if (input1.xyz[i][0] == input2.xyz[i][0] && input1.xyz[i][1] != input2.xyz[i][1] && input1.xyz[i][2] == input2.xyz[i][2])
			{
				double d2 = input2.xyz[i][1] - input1.xyz[i][1];
				double n2 = d2 / (N + 1);
				input3.xyz[i][0] = input1.xyz[i][0];
				input3.xyz[i][1] = input1.xyz[i][1] + j * n2;
				input3.xyz[i][2] = input1.xyz[i][2];
			}
			else if (input1.xyz[i][0] == input2.xyz[i][0] && input1.xyz[i][1] == input2.xyz[i][1] && input1.xyz[i][2] != input2.xyz[i][2])
			{
				double d3 = input2.xyz[i][2] - input1.xyz[i][2];
				double n3 = d3 / (N + 1);
				input3.xyz[i][0] = input1.xyz[i][0];
				input3.xyz[i][1] = input1.xyz[i][1];
				input3.xyz[i][2] = input1.xyz[i][2] + j * n3;
			}
			else if (input1.xyz[i][0] != input2.xyz[i][0] && input1.xyz[i][1] != input2.xyz[i][1] && input1.xyz[i][2] == input2.xyz[i][2])
			{
				double k1 = (input1.xyz[i][1] - input2.xyz[i][1]) / (input1.xyz[i][0] - input2.xyz[i][0]);
				double b1 = input1.xyz[i][1] - k1 * input1.xyz[i][0];
				double d4 = input2.xyz[i][0] - input1.xyz[i][0];
				double n4 = d4 / (N + 1);
				input3.xyz[i][0] = input1.xyz[i][0] + j * n4;
				input3.xyz[i][1] = k1 * input3.xyz[i][0] + b1;
				input3.xyz[i][2] = input1.xyz[i][2];
			}
			else if (input1.xyz[i][0] != input2.xyz[i][0] && input1.xyz[i][1] == input2.xyz[i][1] && input1.xyz[i][2] != input2.xyz[i][2])
			{
				double k2 = (input1.xyz[i][2] - input2.xyz[i][2]) / (input1.xyz[i][0] - input2.xyz[i][0]);
				double b2 = input1.xyz[i][2] - k2 * input1.xyz[i][0];
				double d5 = input2.xyz[i][0] - input1.xyz[i][0];
				double n5 = d5 / (N + 1);
				input3.xyz[i][0] = input1.xyz[i][0] + j * n5;
				input3.xyz[i][1] = input1.xyz[i][1];
				input3.xyz[i][2] = k2 * input3.xyz[i][0] + b2;
			}
			else if (input1.xyz[i][0] == input2.xyz[i][0] && input1.xyz[i][1] != input2.xyz[i][1] && input1.xyz[i][2] != input2.xyz[i][2])
			{
				double k3 = (input1.xyz[i][2] - input2.xyz[i][2]) / (input1.xyz[i][1] - input2.xyz[i][1]);
				double b3 = input1.xyz[i][2] - k3 * input1.xyz[i][1];
				double d6 = input2.xyz[i][1] - input1.xyz[i][1];
				double n6 = d6 / (N + 1);
				input3.xyz[i][0] = input1.xyz[i][0];
				input3.xyz[i][1] = input1.xyz[i][1] + j * n6;
				input3.xyz[i][2] = k3 * input3.xyz[i][1] + b3;
			}
			else
			{
				double k = (input1.xyz[i][1] - input2.xyz[i][1]) / (input1.xyz[i][0] - input2.xyz[i][0]);
				double b = input1.xyz[i][1] - k * input1.xyz[i][0];
				double d = input2.xyz[i][0] - input1.xyz[i][0];
				double n = d / (N + 1);
				input3.xyz[i][0] = input1.xyz[i][0] + j * n;
				input3.xyz[i][1] = k * input3.xyz[i][0] + b;
				input3.xyz[i][2] = (input3.xyz[i][0] - input1.xyz[i][0]) / (input2.xyz[i][0] - input1.xyz[i][0]) * (input2.xyz[i][2] - input1.xyz[i][2]) + input1.xyz[i][2];
			}
		}
		FILE* fp3 = fopen("POSCAR", "w");
		savposcar(fp3, input3);
		fclose(fp3);
		sprintf(name, "0%d", j);
		mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXO);
		sprintf(str, "cp %s %s", "POSCAR", name);
		system(str);
	}
	mkdir("00", S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXO);
	sprintf(name, "0%d", N + 1);
	mkdir(name, S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXO);
	sprintf(str, "cp %s %s", name1, "00/POSCAR");
	system(str);
	sprintf(str, "cp %s %s", "ini.OUTCAR", "00/OUTCAR");
	system(str);
	char tmp1[50];
	sprintf(tmp1, "%s/OUTCAR", name);
	sprintf(str, "cp %s %s", "fin.OUTCAR", tmp1);
	system(str);
	char tmp2[50];
	sprintf(tmp2, "%s/POSCAR", name);
	sprintf(str, "cp %s %s", name2, tmp2);
	system(str);
	remove("POSCAR");
}

int neb_data(int N)
{
	char path[50], find_str1[100], find_str2[100], find_str3[100], save_str1[100], save_str2[100], save_str3[100], save[50], save1[50], save2[50], save3[50], save4[50], file_str[1024];
	int line = 0, l = 0, i = 0;
	double d = 0.000000;
	FILE* fp1, * fp2;
	sprintf(find_str1, "%s", "energy  without entropy");
	sprintf(find_str2, "%s", "projections on to tangent");
	sprintf(find_str3, "%s", "NEB: distance to prev, next image");
	sprintf(save3, "0.000000");
	fp2 = fopen("neb.dat", "w");
	fprintf(fp2, "Number    Distance      Energy     Relative Energy     Force\n");
	for (int j = 0; j <= N; j++)
	{
		sprintf(path, "0%d/OUTCAR", j);
		fp1 = fopen(path, "r");
		if (fp1 == NULL)
		{
			printf("%s is open fail!\n");
			return -1;
		}
		if (j == 0 || j == N)
		{
			sprintf(save2, "    0.000000");
		}
		while (fgets(file_str, sizeof(file_str), fp1))
		{
			line++;
			if (strstr(file_str, find_str1))
			{
				l = line;
				sprintf(save_str1, "%s", file_str);
				strncpy(save1, save_str1 + 29, 15);
			}
			if (strstr(file_str, find_str2))
			{
				i = line;
				sprintf(save_str2, "%s", file_str);
				strncpy(save2, save_str2 + 59, 12);
			}
			if (strstr(file_str, find_str3))
			{
				sprintf(save_str3, "%s", file_str);
				strncpy(save3, save_str3 + 53, 24);
				strncpy(save4, save_str3 + 65, 12);
			}
		}
		if (j == N)
		{
			d = d + atof(save4);
		}
		else
		{
			d = d + atof(save3);
		}
		if (j == 0)
		{
			sprintf(save, save1);
		}
		fprintf(fp2, "   %d      %.6lf   %.6lf      %.6lf        %.6lf\n", j, d, atof(save1), atof(save1) - atof(save), atof(save2));
	}
	fclose(fp1);
	fclose(fp2);
	return 0;
}

int neb_movie(int N)
{
	const char* name[] = { "00/POSCAR", "01/CONTCAR", "02/CONTCAR", "03/CONTCAR", "04/CONTCAR", "05/CONTCAR", "06/CONTCAR", "07/CONTCAR" };
	const char* name0[] = { "00/POSCAR", "01/POSCAR", "02/POSCAR", "03/POSCAR", "04/POSCAR", "05/POSCAR", "06/POSCAR", "07/POSCAR" };
	char path[50], find_str1[100], find_str2[100], save_str1[100], save_str2[100], save1[50], save2[50], file_str[1024];
	POSCAR input;
	FILE* fp1, * fp2;
	fp2 = fopen("movie.xyz", "w");
	if (fp2 == NULL)
	{
		printf("ERROR！\n");
		return -1;
	}
	int line = 0, l = 0, i = 0;
	sprintf(find_str1, "%s", "energy  without entropy");
	sprintf(find_str2, "%s", "FORCES: max atom");
	for (int j = 0; j <= N; j++)
	{
		sprintf(path, "0%d/OUTCAR", j);
		fp1 = fopen(path, "r");
		if (fp1 == NULL)
		{
			printf("ERROR！\n");
			return -1;
		}
		if (j != N)
		{
			FILE* fp = fopen(name[j], "r");
			if (fp == NULL)
			{
				printf("%s is open failed\n", name[j]);
				exit(1);
			}
			readposcar(fp, input);
			direct_to_carts(input.iflg, input.vec, input.xyz, input.nant);
			fclose(fp);
		}
		else
		{
			FILE* fp = fopen(name0[N], "r");
			if (fp == NULL)
			{
				printf("%s is open failed\n", name0[j]);
				exit(1);
			}
			readposcar(fp, input);
			direct_to_carts(input.iflg, input.vec, input.xyz, input.nant);
			fclose(fp);
		}
		fprintf(fp2, "%d\n", input.nant[0]);
		while (fgets(file_str, sizeof(file_str), fp1))
		{
			line++;
			if (strstr(file_str, find_str1))
			{
				l = line;
				sprintf(save_str1, "%s", file_str);
				strncpy(save1, save_str1 + 29, 15);
			}
			if (strstr(file_str, find_str2))
			{
				i = line;
				sprintf(save_str2, "%s", file_str);
				strncpy(save2, save_str2 + 25, 14);
			}
		}
		fprintf(fp2, "F=%lf    E=%lf\n", atof(save2), atof(save1));
		int cnt = 0;
		for (int i = 0; i < input.nant[1]; i++)
		{
			for (int j = 0; j < input.typenum[i]; j++)
			{
				fprintf(fp2, "%s\t%lf\t%lf\t%lf\n", input.elemsym[i], input.xyz[cnt + j][0], input.xyz[cnt + j][1], input.xyz[cnt + j][2]);
			}
			cnt += input.typenum[i];
		}
	}
	fclose(fp1);
	fclose(fp2);
	return 0;
}

void neb(int argc, char* argv[])
{
	if (!strcmp(argv[2], "-sim"))
	{
		double d;
		if (argc >= 5)
			d = neb_similarity(argv[3], argv[4]);
		else
			d = neb_similarity("ini", "fin");
		printf("The difference between ini and fin is %lf\n", d);
		double dis = d / 0.8;
		int N = (int)(dis + 0.5);
		printf("The suggest number of points %d\n", N);
	}
	if (!strcmp(argv[2], "-ins"))
	{
		double d = neb_similarity("ini", "fin");
		if (d <= 5)
		{
			double dis = d / 0.8;
			int N = (int)(dis + 0.5);
			printf("The suggest number of points %d\n", N);
			neb_Points("ini", "fin", N);
		}
		else
		{
			printf("Unreasonable structures\n");
		}
	}
	if (!strcmp(argv[2], "-out"))
	{
		double d = neb_similarity("ini", "fin");
		if (d <= 5)
		{
			double dis = d / 0.8;
			int N = (int)(dis + 0.5);
			printf("The suggest number of points %d\n", N);
			int n;
			n = N + 1;
			neb_data(n);
			neb_movie(n);
		}
	}
}


