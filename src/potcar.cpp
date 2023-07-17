#include"../include/potcar.h"
#include"../include/tools.h"
using namespace std;
string getpath(const char file[], const char key[])
{
	FILE* fp = fopen(file, "r");
	string value;
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return value;
	}
	char buf[1024];
	int iat = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		if (strstr(buf, key) != NULL)
		{
			for (int i = 0; i < strlen(buf); i++)
			{
				if (buf[i] == '#')
					break;
				if (buf[i] == '=')
					iat = 1;
				if (iat)
				{
					if (buf[i] != '=' && buf[i] != ' ')
						value.push_back(buf[i]);
				}
			}
		}
	}
	fclose(fp);
	//delet '\n'
	value.erase(value.end() - 1);
	return value;
}

static void cat_file(const char source_file[], const char desti_file[])
{
	FILE* fp1 = fopen(source_file, "r");
	if (fp1 == NULL)
	{
		printf("%s IS NOT EXIST!\n", source_file);
		return;
	}
	char buf[1024];
	FILE* fp2 = fopen(desti_file, "at+");
	while (fgets(buf, 1024, fp1) != NULL)
		fputs(buf, fp2);
	fclose(fp1);
	fclose(fp2);
}

void pot_merge(const char file[], char mode[], vector<string> label)
{
	FILE* fp1 = fopen(file, "r");
	if (fp1 == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return;
	}
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
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	remove("POTCAR");
	struct passwd* pw = getpwuid(getuid());
	char potpath[100];
	sprintf(potpath, "%s/%s", pw->pw_dir, ".potpath");
	string postfix;
	if (label.size() != 0)
	{
		for (int i = 0; i < label.size(); i++)
		{
			if (i != label.size() - 1)
				postfix += label[i] + "_";
			else
				postfix += label[i];
		}
	}
	for (int i = 0; i < nant[1]; i++)
	{
		if (!strcmp("-PBE", mode))
		{
			string PBE_PATH = getpath(potpath, "PBE");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", PBE_PATH.c_str(), elemnum[i]);
			else
				sprintf(potcar, "%s/POT_%d_%s", PBE_PATH.c_str(), elemnum[i], postfix.c_str());
			cat_file(potcar, "POTCAR");
		}
		else if (!strcmp("-LDA", mode))
		{
			string LDA_PATH = getpath(potpath, "LDA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", LDA_PATH.c_str(), elemnum[i]);
			else
				sprintf(potcar, "%s/POT_%d_%s", LDA_PATH.c_str(), elemnum[i], postfix.c_str());
			cat_file(potcar, "POTCAR");
		}
		else if (!strcmp("-GGA", mode))
		{
			string GGA_PATH = getpath(potpath, "GGA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", GGA_PATH.c_str(), elemnum[i]);
			else
				sprintf(potcar, "%s/POT_%d_%s", GGA_PATH.c_str(), elemnum[i], postfix.c_str());
			cat_file(potcar, "POTCAR");
		}
	}
}

void pot_merge_element(vector<string> element, char mode[], vector<string> label)
{
	remove("POTCAR");
	string path;
	char buf[1024];
	vector<string> v;
	char cur_dir[1024];
	getcwd(cur_dir, 1024);
	while (1)
	{
		getcwd(buf, 1024);
		path = buf;
		int num = count(path.begin(), path.end(), '/');
		string temp;
		for (int i = 0; i < path.size(); i++)
		{
			if (path[i] != '/' && i != path.size() - 1)
				temp += path[i];
			else
			{
				v.push_back(temp);
				temp.clear();
			}
		}
		if ((*(v.end() - 2) == "home"))
			break;
		if (num == 1)
		{
			printf("Can't find root dir!\n");
		}
		chdir("../");
	}
	chdir(cur_dir);
	char potpath[100];
	sprintf(potpath, "%s/%s", path.c_str(), ".potpath");
	string postfix;
	if (label.size() != 0)
	{
		for (int i = 0; i < label.size(); i++)
		{
			if (i != label.size() - 1)
				postfix += label[i] + "_";
			else
				postfix += label[i];
		}
	}
	vector<int> elemnum;
	for (int i = 0; i < element.size(); i++)
	{
		int index;
		getAtomNum(const_cast<char*>(element[i].c_str()), index);
		elemnum.push_back(index);
	}
	for (int i = 0; i < elemnum.size(); i++)
	{
		if (!strcmp("-PBE", mode))
		{
			string PBE_PATH = getpath(potpath, "PBE");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", PBE_PATH.c_str(), elemnum[i]);
			else
				sprintf(potcar, "%s/POT_%d_%s", PBE_PATH.c_str(), elemnum[i], postfix.c_str());
			cat_file(potcar, "POTCAR");
		}
		else if (!strcmp("-LDA", mode))
		{
			string LDA_PATH = getpath(potpath, "LDA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", LDA_PATH.c_str(), elemnum[i]);
			else
				sprintf(potcar, "%s/POT_%d_%s", LDA_PATH.c_str(), elemnum[i], postfix.c_str());
			cat_file(potcar, "POTCAR");
		}
		else if (!strcmp("-GGA", mode))
		{
			string GGA_PATH = getpath(potpath, "GGA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", GGA_PATH.c_str(), elemnum[i]);
			else
				sprintf(potcar, "%s/POT_%d_%s", GGA_PATH.c_str(), elemnum[i], postfix.c_str());
			cat_file(potcar, "POTCAR");
		}
	}
}

void check()
{
	int flag = 0;
	if (access("INCAR", 0))
	{
		printf("Error: The INCAR File NOT Found.\n");
		flag = 1;
	}
	if (access("POSCAR", 0))
	{
		printf("Error: The POSCAR File NOT Found.\n");
		flag = 1;
	}
	if (access("KPOINTS", 0))
	{
		printf("Error: The KPOINTS File NOT Found.\n");
		flag = 1;
	}
	if (access("POTCAR", 0))
	{
		printf("Error: The POTCAR File NOT Found.\n");
		flag = 1;
	}
	if (flag)
		return;
	else
		printf(" INCAR, POSCAR, KPOINTS and POTCAR seem to be OK.\n");
	FILE* fp1 = fopen("POSCAR", "r");
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
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp1);
	FILE* fp2 = fopen("POTCAR", "r");
	char buf[1024];
	vector<string> potelem;
	while (fgets(buf, 1024, fp2) != NULL)
	{
		if (strstr(buf, "VRHFIN") != NULL)
		{
			string elem;
			int iat = 0;
			for (int i = 0; i < strlen(buf); i++)
			{
				if (buf[i] == '=')
					iat = 1;
				if (buf[i] == ':')
					break;
				if (iat)
				{
					if (buf[i] != '=' && buf[i] != ' ')
						elem.push_back(buf[i]);
				}
			}
			potelem.push_back(elem);
		}
	}
	fclose(fp2);
	int potflag = 0;
	if (potelem.size() != nant[1])
		potflag = 1;
	else
	{
		for (int i = 0; i < potelem.size(); i++)
		{
			if (potelem[i] == elemsym[i])
				continue;
			else
			{
				potflag = 1;
				break;
			}
		}
	}
	if (!potflag)
		printf("Now you can submit VASP job.\n");
	else
	{
		printf("Element in POSCAR not corresponding to POTCAR.\n");
		printf("POTCAR: ");
		for (int i = 0; i < potelem.size(); i++)
			printf("%s ", potelem[i].c_str());
		printf("\n");
		printf("POSCAR: ");
		for (int i = 0; i < nant[1]; i++)
			printf("%s ", elemsym[i]);
		printf("\n");
	}
	string ISTART = GetInfoINCAR("ISTART");
	string ICHARG = GetInfoINCAR("ICHARG");
	if (ISTART.size() != 0 || ISTART == "1" || ISTART == "2" || ISTART == "3")
	{
		FILE* fp = fopen("WAVECAR", "r");
		if (fp == NULL)
			printf("Attention:ISTART = %d, But No WAVECAR is found!\n", atoi(ISTART.c_str()));
	}
	if (ICHARG.size() != 0 || ICHARG == "11" || ICHARG == "12" || ICHARG == "13")
	{
		FILE* fp = fopen("CHGCAR", "r");
		if (fp == NULL)
			printf("Attention:ICHARG = %d, But No CHGCAR is found!\n", atoi(ICHARG.c_str()));
	}
}