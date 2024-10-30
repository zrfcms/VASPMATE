#include <iostream>
#include <sstream>
#include <string.h>
#include "../../include/Stochastics_include/Stochastics_main.h"
#include "../../include/Energy_Minimization/Energy_Minimization_main.h"
#include "../../include/Stochastics_include/generate.h"

class sto_Generate
{
public:
	sto_Generate(int argc, char* argv[]);
	std::vector<int> init_atom;
	std::vector<string> init_element;
	int POPULATION_NUM = 0, formula_num = 0, max_step = 10000, md = 1;
	double init_volume = 0, scale = 0.1;
};

template<typename T>
std::vector<T> Parse(const char temp[])
{
    vector<T> ans;
    stringstream ss((std::string)temp);
    string item;
    while (getline(ss, item, ',')) 
	{
        stringstream converter(item);
        T a;
        if (converter >> a)
            ans.push_back(a);
    }
    return ans;
}

sto_Generate::sto_Generate(int argc, char* argv[])
{
	int ato_numb = 2;
	int ele_numb = 0;
	int sto_num = 0;
	for (int i = 0; i < argc; i++)
	{
		if (!strcmp("-atom", argv[i]))
			ato_numb = i;
		else if (!strcmp("-elem", argv[i]))
			ele_numb = i;
		else if (!strcmp("-num", argv[i]))
			POPULATION_NUM = atoi(argv[i + 1]);
		else if (!strcmp("-fav", argv[i]))
		{
			formula_num = atoi(argv[i + 1]);
			init_volume = atof(argv[i + 2]);
		}
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
	sto_num = ele_numb - ato_numb - 1;
	if (sto_num == 0)
		printf("Please enter the number of atoms for each element!\n");
	init_atom.resize(sto_num);
	init_element.resize(sto_num);
	for (int i = 0; i < sto_num; i++)
	{
		init_atom[i] = atoi(argv[ato_numb + i + 1]);
		init_element[i]	= argv[ele_numb + i + 1];
	}
}

int Stochastics_main(int argc, char* argv[])
{
	srand(time(0));
	sto_Generate gen(argc, argv);
	std::vector<string> output = gene_first_population_Random(gen.init_atom, gen.init_element, gen.POPULATION_NUM, gen.formula_num, gen.init_volume);
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
		for (string f : output)
		{	
			string f_md = f + string("_md");
			argv_md[1] = (char*)malloc((strlen(f.c_str())+1)*sizeof(char));
			argv_md[2] = (char*)malloc((strlen(f_md.c_str())+1)*sizeof(char));
			strcpy(argv_md[1],f.c_str());
			strcpy(argv_md[2],(f_md.c_str()));
			Energy_Minimization_main(argc_md , argv_md);
		}
	}
	return 0;
}