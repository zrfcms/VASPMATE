#include<vector>
#include<string>
#include<iostream>
#include "../../include/Evolutionary_include/Evolutionary_main.h"
#include "../../include/Energy_Minimization/Energy_Minimization_main.h"
#include "../../include/Evolutionary_include/variant.h"
//#include<Windows.h>
using namespace std;
vector<string> Variant{ "Stripple","Permustrain","Randmove","Slip","Twist","Rotstrain","Crossover"};
class evo_Generate
{
public:
	evo_Generate(int argc, char* argv[]);
	int POPULATION_NUM = 0, formula_num = 0, max_step = 10000, md = 1;
	double scale = 0.1;
	char parent1[256], parent2[256];
	vector<int> var;
};

evo_Generate::evo_Generate(int argc, char* argv[])
{
	var = { 1, 1, 1, 1, 1, 1, 1};
	for (int i = 0; i < argc; i++)
		if (!strcmp("crossover", argv[i]))
			strcpy(parent2, argv[3]);
	strcpy(parent1, argv[2]);
	for (int i = 0; i < argc; i++)
	{
		if (!strcmp("Stripple", argv[i]))
			var[0] = 0;
		else if (!strcmp("Permustrain", argv[i]))
			var[1] = 0;
		else if (!strcmp("Randmove", argv[i]))
			var[2] = 0;
		else if (!strcmp("Slip", argv[i]))
			var[3] = 0;
		else if (!strcmp("Twist", argv[i]))
			var[4] = 0;
		else if (!strcmp("Rotstrain", argv[i]))
			var[5] = 0;
		else if (!strcmp("crossover", argv[i]))
			var[6] = 0;
		else if (!strcmp("-n", argv[i]))
			POPULATION_NUM = atoi(argv[i + 1]);
		else if (!strcmp("-min", argv[i]))
		{
			if (atoi(argv[i + 1]) == 0)
			{
				md = 0;
				max_step = atoi(argv[i + 2]);
				scale = atof(argv[i + 3]);
			}
			else if (atoi(argv[i + 1]) == 1)
				md = 1 ;
		}
	}
}

int Evolutionary_main(int argc, char* argv[])
{
	//srand(time(0));
	evo_Generate gen(argc , argv);
	string variant_log;
	for (int i = 0; i < 7; i++)
		if (gen.var[i] == 0)
			variant_log += (string)Variant[i] + " ";
	cout << variant_log << endl ;
	if (variant_log.empty())
	{
		variant_log = "No variant is chosen, please check!";
		return 0;
	}
	for (int i = 0; i < gen.POPULATION_NUM; i++)
	{
		char POSCAR_Child[50];
		sprintf(POSCAR_Child, "POSCAR%d", i);
		int method = rand() % 7;
		while (gen.var[method] == 1)  method = rand() % 7;
		if (method != 6)
		{
			if (access(gen.parent1, 0))
				printf("Error: The %s File NOT Found.\n", gen.parent1 );
			else
				Variant_singer_parent(gen.parent1, POSCAR_Child, method, 1);
		}
		else
		{
			if (access(gen.parent1, 0) || access(gen.parent2, 0))
				printf("Error: The %s or %s File NOT Found.\n", gen.parent1, gen.parent2);
			else
				Variant_two_parent(gen.parent1, gen.parent2, POSCAR_Child, 1);
		}
		if (!gen.md)
		{
			cout << "start md..." << endl ;
			int argc_md = 5;
			char** argv_md=(char**)malloc(5*sizeof(char*));
			argv_md[0] = 0; 
			for (int i = 0; i < argc; i++)
			{
				if (!strcmp("-min", argv[i]))
				{
					argv_md[3] = argv[i + 2]; //max_step
					argv_md[4] = argv[i + 3]; //scale
				}
			}
			argv_md[1] = (char*)malloc((strlen(POSCAR_Child)+1)*sizeof(char));
			strcpy(argv_md[1],POSCAR_Child);
			string f_md = POSCAR_Child + string("_md");
			argv_md[2] = (char*)malloc((strlen(f_md.c_str())+1)*sizeof(char));
			strcpy(argv_md[2],f_md.c_str());
			Energy_Minimization_main(argc_md , argv_md);
		}
	}
	return 0;
}
