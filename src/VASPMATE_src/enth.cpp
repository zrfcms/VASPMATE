#include"../../include/VASPMATE_include/enth.h"
#include"../../include/VASPMATE_include/mapdata.h"
#include"../../include/VASPMATE_include/outcar.h"
using namespace enth;

int enth::operator_enth(int argc, char* argv[])
{
	int p; int path;
	double en1 = 0; double en2 = 0;
	int check_par = check_para(argc, argv, "-ref", 1 , &p); //The fourth parameter is the minimum number of input parameters required.
	int check_path = check_para(argc, argv, "-path", 1 , &path); //The fourth parameter is the minimum number of input parameters required.
	string path_;
	vector<double> ref_en; ref_en.resize(0);
	bool is_muple = true;
    //if -ref exit, get number to table; else ulse enthalpies table
    if (check_par == 1)
	{
		for (int i = p + 1; i < argc; i++)
		{
			if (VM_isnumb(argv[i]))
				ref_en.push_back(atof(argv[i]));
			else
				break;
		}
	}
    //if -path exits, get hull; else signle enth
	if (check_path == 1)
	    path_ = string(argv[path + 1]);
	else if (check_path == 0)
        path_ = "./";
	else
		is_muple = false;
	if (is_muple)
		enth::micid_convexhull(path_, ref_en);
	else
		enth::Formation_Enthalpy(ref_en);
	return 0;
}

void enth::micid_convexhull(string path_, vector<double> ref_en)
{
	chdir(path_.c_str());
	int foldername = 1;
	vector<double> ratio;
	vector<double> result;
	ratio.push_back(0); result.push_back(0); //(0,0)
	ratio.push_back(1); result.push_back(0); //(1,0)
	while(foldername)
	{
		char pos_file[10];
		char eng_file[10];
		string index_name = string(4 - to_string(foldername).length(), '0') + to_string(foldername);
        if (!access(index_name.c_str(), 0))
			chdir(index_name.c_str());
		else
			break;
		//get energy
		double energy = 0;
		if(VM_fileexist("OSZICAR"))
			energy = get_energy_oszi();
		else
		{
			printf("There is no OSZICAR file in %s, VASPMATE skips this data collection!\n", index_name.c_str());
			foldername++;
			chdir("..");
			continue;
		}
		//get pos
		if(VM_fileexist("CONTCAR"))
			strcpy(pos_file, "CONTCAR");
		else if (VM_fileexist("POSCAR"))
			strcpy(pos_file, "POSCAR");
		else
		{
			printf("There is no CONTCAR/POSCAR file in %s, VASPMATE skips this data collection!\n", index_name.c_str());
			foldername++;
			chdir("..");
			continue;
		}
		FILE* fp = fopen(pos_file, "r");
		POSCAR pos_hull;
		readposcar(fp, pos_hull);
		if (pos_hull.nant[1] != 2)
		{
			printf("The number of structural elements in %s is not 2, VASPMATE skips this data collection!\n", index_name.c_str());
			foldername++;
			chdir("..");
			continue;
		}
		if (foldername == 1)
		{
			if (ref_en.size() < 2)
			{
				printf("Lack input reference values! Use the default elemental energy!\n");
			}
			for (int i = 0; i < ((ref_en.size() < pos_hull.nant[1]) ? ref_en.size() : pos_hull.nant[1]); i++)
				enthalpies[string(pos_hull.elemsym[i])] = ref_en[i];
			printf("The reference values for elements are: ");
			for (int i = 0; i < pos_hull.nant[1]; i++)
				printf("%2s %f;", pos_hull.elemsym[i], enthalpies[string(pos_hull.elemsym[i])]);
			printf("\n\n");
		}
		double en1 = enthalpies[string(pos_hull.elemsym[0])];
		double en2 = enthalpies[string(pos_hull.elemsym[1])];
		double rato = static_cast<double>(pos_hull.typenum[0]) / pos_hull.nant[0];
		ratio.push_back(rato);
		result.push_back((energy - en1*pos_hull.typenum[0] - en2*pos_hull.typenum[1]) / pos_hull.nant[0]);
		fclose(fp);
		chdir("..");
		foldername++;
	}
	FILE* fp_csv = fopen("convex_hull.csv", "w");
	fprintf(fp_csv,"ratio,formation_energy/(eV/atom)\n");
	for (int i = 0; i< ratio.size(); i++)
	{
		fprintf(fp_csv,"%lf,%lf\n", ratio[i], result[i]);
	}
	fclose(fp_csv);
	printf("Written convex_hull.csv file!\n");
}

void enth::Formation_Enthalpy(vector<double> ref_en)
{
	double energy = 0;
	if(VM_fileexist("OSZICAR"))
		energy = get_energy_oszi();
	else
	{
		printf("There is no OSZICAR file to read energy!\n");
		return;
	}
	char pos_file[10];
	if(VM_fileexist("CONTCAR"))
		strcpy(pos_file, "CONTCAR");
	else if (VM_fileexist("POSCAR"))
		strcpy(pos_file, "POSCAR");
	else
	{
		printf("There is no CONTCAR/POSCAR file to read POS!\n");
		return;
	}
	FILE* fp = fopen(pos_file, "r");
	POSCAR pos;
	readposcar(fp, pos);
	for (int i = 0; i < ((ref_en.size() < pos.nant[1]) ? ref_en.size() : pos.nant[1]); i++)
		enthalpies[string(pos.elemsym[i])] = ref_en[i];
	for (int i = 0; i < pos.nant[1]; i++)
		energy -= pos.typenum[i] * enthalpies[string(pos.elemsym[i])];
	energy /= pos.nant[0];
	FILE* fp_wri = fopen("ENTH_INFO.dat" , "w");
	printf("%lf\n", energy);
	fprintf(fp_wri, "%lf\n\n", energy);
	fprintf(fp_wri, "# The reference values for elements are:\n");
	for (int i = 0; i < pos.nant[1]; i++)
		fprintf(fp_wri, "# %2s %lf\n", pos.elemsym[i], enthalpies[string(pos.elemsym[i])]);
	fclose(fp_wri);
	printf("Written ENTH_INFO.dat file!\n");
}