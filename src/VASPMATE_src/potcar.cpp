#include"../../include/VASPMATE_include/tools.h"
#include"../../include/VASPMATE_include/potcar.h"

using namespace std;
string getpath(const char file[], char key[])
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
	remove("NEWPOT");
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
		if (!strcmp("-PBE", mode) || !strcmp("PBE", mode))
		{
			string PBE_PATH = getpath(potpath, "PBE");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", PBE_PATH.c_str(), elemnum[i]);
			else
			{
				sprintf(potcar, "%s/POT_%d_%s", PBE_PATH.c_str(), elemnum[i], postfix.c_str());
				if (access(potcar, F_OK) != 0)
					sprintf(potcar, "%s/POT_%d", PBE_PATH.c_str(), elemnum[i]);
			}
			cat_file(potcar, "NEWPOT");
		}
		else if (!strcmp("-LDA", mode) || !strcmp("LDA", mode))
		{
			string LDA_PATH = getpath(potpath, "LDA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", LDA_PATH.c_str(), elemnum[i]);
			else
			{
				sprintf(potcar, "%s/POT_%d_%s", LDA_PATH.c_str(), elemnum[i], postfix.c_str());
				if (access(potcar, F_OK) != 0)
					sprintf(potcar, "%s/POT_%d", LDA_PATH.c_str(), elemnum[i]);
			}
			cat_file(potcar, "NEWPOT");
		}
		else if (!strcmp("-GGA", mode) || !strcmp("GGA", mode))
		{
			string GGA_PATH = getpath(potpath, "GGA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", GGA_PATH.c_str(), elemnum[i]);
			else
			{
				sprintf(potcar, "%s/POT_%d_%s", GGA_PATH.c_str(), elemnum[i], postfix.c_str());
				if (access(potcar, F_OK) != 0)
					sprintf(potcar, "%s/POT_%d", GGA_PATH.c_str(), elemnum[i]);
			}
			cat_file(potcar, "NEWPOT");
		}
	}
	printf("Written NEWPOT file!\n");
}

void pot_merge_element(vector<string> element, char mode[], vector<string> label)
{
	remove("NEWPOT");
	string path;
	char buf[1024];
	vector<string> v;
	/*char cur_dir[1024];
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
	chdir(cur_dir);*/
	struct passwd* pw = getpwuid(getuid());
	char potpath[100];
	sprintf(potpath, "%s/%s", pw->pw_dir, ".potpath");
	//char potpath[100];
	//sprintf(potpath, "%s/%s", path.c_str(), ".potpath");
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
		if (!strcmp("-PBE", mode) || !strcmp("PBE", mode))
		{
			string PBE_PATH = getpath(potpath, "PBE");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", PBE_PATH.c_str(), elemnum[i]);
			else
			{
				sprintf(potcar, "%s/POT_%d_%s", PBE_PATH.c_str(), elemnum[i], postfix.c_str());
				if (access(potcar, F_OK) != 0)
					sprintf(potcar, "%s/POT_%d", PBE_PATH.c_str(), elemnum[i]);
			}
			cat_file(potcar, "NEWPOT");
		}
		else if (!strcmp("-LDA", mode) || !strcmp("LDA", mode))
		{
			string LDA_PATH = getpath(potpath, "LDA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", LDA_PATH.c_str(), elemnum[i]);
			else
			{
				sprintf(potcar, "%s/POT_%d_%s", LDA_PATH.c_str(), elemnum[i], postfix.c_str());
				if (access(potcar, F_OK) != 0)
					sprintf(potcar, "%s/POT_%d", LDA_PATH.c_str(), elemnum[i]);
			}
			cat_file(potcar, "NEWPOT");
		}
		else if (!strcmp("-GGA", mode) || !strcmp("GGA", mode))
		{
			string GGA_PATH = getpath(potpath, "GGA");
			char potcar[100];
			if (label.size() == 0)
				sprintf(potcar, "%s/POT_%d", GGA_PATH.c_str(), elemnum[i]);
			else
			{
				sprintf(potcar, "%s/POT_%d_%s", GGA_PATH.c_str(), elemnum[i], postfix.c_str());
				if (access(potcar, F_OK) != 0)
					sprintf(potcar, "%s/POT_%d", GGA_PATH.c_str(), elemnum[i]);
			}
			cat_file(potcar, "NEWPOT");
		}
	}
	printf("Written NEWPOT file!\n");
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

void check_one(int filenumb)
{
	if(filenumb == 1)
	{
		int is_incar_ok = 0; //-1 error; 0 ok ; 1 attention
		if (access("INCAR", 0))
			{ printf("Error: The INCAR File NOT Found.\n"); is_incar_ok = -1; }
		else
		{
			string ISTART = GetInfoINCAR("ISTART");
			string ICHARG = GetInfoINCAR("ICHARG");
			string LDAU = GetInfoINCAR("LDAU");
			string LDAUTYPE = GetInfoINCAR("LDAUTYPE");
			vector<string> LDAUL = GetInfoINCAR("LDAUL", "");
			vector<string> LDAUU = GetInfoINCAR("LDAUU", "");
			vector<string> LDAUJ = GetInfoINCAR("LDAUJ", "");
			if (LDAU.size() != 0 && LDAU == ".TRUE.")
			{
				if (access("POSCAR", 0))
					printf("POSCAR file not found.\nVASPMATE will stop checking the consistency of LDAUL/LDAUU/LDAUJ with the number of elements in POSCAR.\n");
				else
				{
					POSCAR pos;
					if (LDAUTYPE.size() == 0)
					{
						printf("Attention: VASPATE noticed that LDAUTYPE is not set in INCAR.\n");
						is_incar_ok = 1;
					}
					if (LDAUL.size() == 0)
					{ 
						printf("Attention: VASPATE noticed that LDAUL is not set in INCAR.\n");  
						is_incar_ok = 1;
					}
					else if (LDAUL.size() != pos.nant[0])
					{
						printf("Attention: VASPMATE noticed that the parameters of LDAUL in INCAR are inconsistent with the number of atoms in POSCAR.\n");
						is_incar_ok = 1;
					}
					if (LDAUU.size() == 0)
					{
						printf("Attention: VASPATE noticed that LDAUU is not set in INCAR.\n");
						is_incar_ok = 1;
					}
					else if (LDAUU.size() != pos.nant[0])
					{
						printf("Attention: VASPMATE noticed that the parameters of LDAUU in INCAR are inconsistent with the number of atoms in POSCAR.\n");
						is_incar_ok = 1;
					}
					if (LDAUJ.size() == 0)
					{
						printf("Attention: VASPATE noticed that LDAUJ is not set in INCAR.\n");
						is_incar_ok = 1;
					}
					else if (LDAUJ.size() != pos.nant[0])
					{
						printf("Attention: VASPMATE noticed that the parameters of LDAUJ in INCAR are inconsistent with the number of atoms in POSCAR.\n");
						is_incar_ok = 1;
					}
				}
			}
			if (ISTART.size() != 0 || (ISTART == "1" || ISTART == "2" || ISTART == "3"))
			{
				FILE* fp = fopen("WAVECAR", "r");
				if (fp == NULL)
					printf("Attention:ISTART = %d, But No WAVECAR is found!\n", atoi(ISTART.c_str()));
			}
			if (ICHARG.size() != 0 || (ICHARG == "11" || ICHARG == "12" || ICHARG == "13"))
			{
				FILE* fp = fopen("CHGCAR", "r");
				if (fp == NULL)
					printf("Attention:ICHARG = %d, But No CHGCAR is found!\n", atoi(ICHARG.c_str()));
			}
			string ISPIN = GetInfoINCAR("ISPIN");
			vector<string> MAGMOM = GetInfoINCAR("MAGMOM", "");			
			if (ISPIN.size() != 0 || ISPIN == "2")
			{
				if (access("POSCAR", 0))
					printf("POSCAR file not found.\nVASPMATE will stop checking the consistency of MAGMOM with the number of elements in POSCAR.\n");
				else
				{
					POSCAR pos;
					if (MAGMOM.size() == 0)
					{
						printf("Attention: VASPATE noticed that MAGMOM is not set in INCAR.\n");
						is_incar_ok = 1;
					}
					else if (MAGMOM.size() != pos.nant[0])
					{
						printf("Attention: VASPMATE noticed that the parameters of MAGMOM in INCAR are inconsistent with the number of atoms in POSCAR.\n");
						is_incar_ok = 1;
					}
				}
			}
			string IVDW = GetInfoINCAR("IVDW");
			string LUSE_VDW = GetInfoINCAR("LUSE_VDW");
			if (IVDW.size() != 0 || LUSE_VDW.size() != 0)
			{
				if (LUSE_VDW == ".TRUE." && IVDW == "0")
				{
					printf("Attention: 'IVDW = 0' &&  'LUSE_VDW = .TRUE.' should not appear at the same time in the INCAR.\n");
					is_incar_ok = 1;
				}
				else if (LUSE_VDW == ".FALSE." && IVDW != "0")
				{
					printf("Attention: 'IVDW = %s' &&  'LUSE_VDW = .FALSE.' should not appear at the same time in the INCAR.\n", IVDW.c_str());
					is_incar_ok = 1;
				}
			}
			string ISIF = GetInfoINCAR("ISIF");
			string NSW = GetInfoINCAR("NSW");
			if (ISIF.size() != 0 || NSW.size() != 0)
			{
				if (ISIF == "3" && NSW == "0")
				{
					printf("Attention: 'ISIF = 3' &&  'NSW = 0' should not appear at the same time in the INCAR.\n");
					is_incar_ok = 1;
				}
			}
			if (is_incar_ok == 0)
				printf("The INCAR File seem to be OK.\n");
			else if (is_incar_ok == 1)
				printf("There may be some important points in the INCAR file that need attention.\n");
		}
	}
	if(filenumb == 2)
	{
		if (access("POSCAR", 0))
			printf("Error: The POSCAR File NOT Found.\n");
		else
			printf("The POSCAR File seem to be OK.\n");
	}
	if(filenumb == 3)
	{
		if (access("POTCAR", 0))
			printf("Error: The POTCAR File NOT Found.\n");
		else
		{
			FILE* fp1 = fopen("POSCAR", "r");
			POSCAR pos;
			readposcar(fp1,pos);
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
			if (potelem.size() != pos.nant[1])
				potflag = 1;
			else
			{
				for (int i = 0; i < potelem.size(); i++)
				{
					if (potelem[i] == pos.elemsym[i])
						continue;
					else
					{
						potflag = 1;
						break;
					}
				}
			}
			if (!potflag)
				printf("The POTCAR File seem to be OK.\n");
			else
			{
				printf("Attention：Element is not consistently matched!\n");
				printf("POTCAR: ");
				for (int i = 0; i < potelem.size(); i++)
					printf("%s ", potelem[i].c_str());
				printf("\n");
				printf("POSCAR: ");
				for (int i = 0; i < pos.nant[1]; i++)
					printf("%s ", pos.elemsym[i]);
				printf("\n");
			}
		}
	}
	if(filenumb == 4)
	{
		if (access("KPOINTS", 0))
			printf("Error: The KPOINTS File NOT Found.\n");
		else
			printf("The KPOINTS File seem to be OK.\n");
	}
}