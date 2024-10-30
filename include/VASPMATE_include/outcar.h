#pragma once
#ifndef _OUTCAR_
#define _OUTCAR_
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include <unordered_map>
#include <string>
#include <vector>
#include<algorithm>
#include <functional>
#include <cassert>
#include <iostream>
#include <string.h>
#include <tuple>
#include <math.h>
#if defined(__APPLE__) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#ifdef __APPLE__
#define fseeko fseek
#define ftello ftell
#endif
#include <fstream>
#include"cJSON.h"
#include"read_write.h"
#include"tinyxml2.h"
#include"spa_plot.h"
#include"tools.h"
using namespace std;
using namespace tinyxml2;
typedef vector<vector<vector<double> >  > DIELECT;
DIELECT get_dielectric(const char file[] = "vasprun.xml");
void get_linear_optical_spectrums_2d(DIELECT diel);
void get_linear_optical_spectrums_3d(DIELECT diel);
double get_energy();
double get_energy_oszi();
vector<double> get_stress();
namespace _outcar {
	vector<string> str_split(string s, vector<char> sp);
	struct Para_Info
	{
		Para_Info(const char key[], const char split[], int num, const char name[], function<void(FILE* fp, Para_Info*)> f)
		{
			strcpy(this->key, key);
			strcpy(this->split, split);
			this->name = name;
			this->num = num;
			this->print = f;
		}
		char key[100];
		char split[2];
		int num;
		string name;
		function<void(FILE* fp, Para_Info*)> print;
	};
	static unordered_map<string, string> m;
	void print_normal(FILE* fp, Para_Info* info);
	void print_null(FILE* fp, Para_Info* info);
	void print_GGA(FILE* fp, Para_Info* info);
	void print_ISPIN(FILE* fp, Para_Info* info);
	void print_ISIF(FILE* fp, Para_Info* info);
	void print_else(FILE* fp);
	void print_text(FILE* fp, const char name[], vector<string> v);
	static vector<Para_Info> parameter{
			{"POTCAR",":",1,"Potential_type",print_normal},
			{"LEXCH","=",1,"", print_null},
			{"GGA","=",1,"Functional_type", print_GGA},
			{"ISPIN","=",1,"Spin_polarization",print_ISPIN},
			{"LSORBIT", "=", 1, "Spin_orbit_coupling",print_normal},
			{"PREC","=",1,"Precise",print_normal},
			{"ENCUT","=",1,"Energy_cutoff",print_normal},
			{"ISIF","=",1,"Relaxation",print_ISIF},
			{"PSTRESS","=",1,"Pullay_stress",print_normal},
			{"EDIFFG","=",1,"Ionic_Convergence",print_normal},
			{"EDIFF","=",1,"Electronic_convergence",print_normal},
			{"ISMEAR","=",1,"Integration_scheme",print_normal},
			{"ALGO","=",1,"Minization_algorithm",print_normal},
			{"IBRION","=",1,"Ionic_update",print_normal},
			{"LSOL","=",1,"Solvation_Model",print_normal},
			//self-consistent
			{"ISTART","=",1,"",print_null},
			{"ICHARG","=",1,"",print_null},
			{"NSW","=",1,"",print_null},
			//metagga
			{"METAGGA","=",1,"MetaGGA_type",print_normal},
			//hse
			{"LHFCALC","=",1,"",print_null},
			{"HFSCREEN","=",1,"",print_null},
			{"AEXX","=",1,"",print_null},
			{"LIBXC1","=",1,"",print_null},
			//ldau
			{"LDAUTYPE","=",1,"",print_null},
			//vdw
			{"IVDW","=",1,"",print_null},
			{"LVDW","=",1,"",print_null},
			{"LUSE_VDW","=",1,"",print_null},
			{"LDIPOL","=",1,"Dipole_correction",print_normal},
	};
	class OUTCAR
	{
	public:
		OUTCAR() {
			for (int i = 0; i < parameter.size(); i++)
			{
				m[parameter[i].key] = "";
				DB_GetInfo(parameter[i]);
			}
		}
		void DB_GetInfo(Para_Info para, const char input_file[] = "OUTCAR");
		void output(const char output[] = "log.sdata", std::string flag = "-a", vector<string> file_name = {});
		vector<string> energy;
		vector<string> par_energy;
		vector<string> stress;
		vector<string> force;
		vector<string> charge;
		vector<string> magnitization_x, magnitization_y, magnitization_z;
	};
}
static void TranMapToJson(unordered_map<string, vector<string> >& data, cJSON* obj);
void TranLogdataToJson(const char inputfile[] = "log.sdata", const char outputfile[] = "log.json",  int MaxLogNum = 1024);
void TranJsonToLogdata(const char inputfile[] = "log.json", const char outputfile[] = "log.sdata",  int MaxLogNum = 1024);
void ParseJsonToMap(const cJSON* obj, unordered_map<string, vector<string>>& data);
void TranFileToCsv(const char inputfile[], const char outputfile[]);
void TranCsvToFile(const char inputfile[], const char outputfile[]);
tuple<vector<int>, vector<double>, int> EnergyConvergence(double EDIFFG);
tuple<vector<int>, vector<double>, int> ForceConvergence(double EDIFFG);
void RelaxJudgeConvergence();
tuple<vector<int>, vector<double>, int> StaticEnergyConvergence(double EDIFF);
tuple<vector<int>, vector<double>, int> mdConvergence(double EDIFF);
void StaticJudgeConvergence();
void mdJudgeConvergence();

void PlusLog(vector<string> filename = {}, const char outputfile[] = "newlog.sdata");
void PlusJson(vector<string> filename = {}, const char outputfile[] = "newlog.json");
#endif