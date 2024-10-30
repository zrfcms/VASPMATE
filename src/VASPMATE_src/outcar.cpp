#include"../../include/VASPMATE_include/outcar.h"
using namespace _outcar;
using namespace std;
#define NOEXIST(file) printf("%s IS NOT EXIST!\n",file);
#define blankline(buf) strspn(buf,"\t\n\r ") == strlen(buf)
#define WriteFile(file) printf("Written %s file!\n",file);

static const double epsilon0 = 8.8541878171E-12;
static const double sigma0 = 6.0853316295796981E-5;
static const double planck_constant = 4.13566733E-15;
static const double c0 = 2.99792458E8;
static const double pi = acos(-1);

DIELECT get_dielectric(const char file[])
{
	FILE* fp1 = fopen("REAL.IN", "w");
	FILE* fp2 = fopen("IMAG.IN", "w");
	fprintf(fp1, "   E(ev)      X         Y         Z        XY        YZ        ZX\n");
	fprintf(fp2, "   E(ev)      X         Y         Z        XY        YZ        ZX\n");
	DIELECT diel(2, vector<vector<double> >(7));
	tinyxml2::XMLDocument xml;
	if (xml.LoadFile(file) != XML_SUCCESS)
		return {};
	tinyxml2::XMLElement* rootNode = xml.RootElement();
	if (rootNode == NULL)
		return {};
	vector<string> imag_label{ "calculation" ,"dielectricfunction","imag","array","set","r" };
	for (auto val : imag_label)
		rootNode = rootNode->FirstChildElement(val.c_str());
	while (rootNode)
	{
		fprintf(fp2, "%s\n", rootNode->GetText());
		vector<double> tmp(7);
		sscanf(rootNode->GetText(), "%lf%lf%lf%lf%lf%lf%lf", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6]);
		for (int i = 0; i < 7; i++)
			diel[1][i].push_back(tmp[i]);
		rootNode = rootNode->NextSiblingElement();
	}
	rootNode = xml.RootElement();
	vector<string> real_label{ "calculation" ,"dielectricfunction","real","array","set","r" };
	for (auto val : real_label)
		rootNode = rootNode->FirstChildElement(val.c_str());
	while (rootNode)
	{
		fprintf(fp1, "%s\n", rootNode->GetText());
		vector<double> tmp(7);
		sscanf(rootNode->GetText(), "%lf%lf%lf%lf%lf%lf%lf", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6]);
		for (int i = 0; i < 7; i++)
			diel[0][i].push_back(tmp[i]);
		rootNode = rootNode->NextSiblingElement();
	}
	fclose(fp1);
	fclose(fp2);
	return diel;
}

void get_linear_optical_spectrums_2d(DIELECT diel)
{
	FILE* fp = fopen("POSCAR", "r");
	if (fp == nullptr)
	{
		NOEXIST("POSCAR");
		return;
	}
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	double lattice = sqrt(pos.vec[2][0] * pos.vec[2][0] + pos.vec[2][1] * pos.vec[2][1] + pos.vec[2][2] * pos.vec[2][2]);
	int n = diel[0][0].size();
	vector<double> omega(n);
	for (int i = 0; i < n; i++)
		omega[i] = 2 * pi * diel[0][0][i] / planck_constant;
	vector<vector<double> > real_polarizability_2D(2, vector<double>(n));
	vector<vector<double> > imag_polarizability_2D(2, vector<double>(n));
	vector<double> denominator(n);
	vector<vector<double> > transmission(2, vector<double>(n));
	vector<vector<double> > absorb(2, vector<double>(n));
	vector<vector<double> > reflection(2, vector<double>(n));
	vector<string> filename = { "real_optical_conductivity_2d","imag_optical_conductivity_2d","transmission_2d","absorption_2d","reflection_2d" };
	vector<FILE*> fps(filename.size());
	for (int i = 0; i < fps.size(); i++)
	{
		fps[i] = fopen(filename[i].c_str(), "w");
		i < 2 ? fprintf(fps[i], "   #Energy     XX     YY\n") : fprintf(fps[i], "   #Energy     XX(\%)     YY(\%)\n");
	}
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < n; j++)
		{
			real_polarizability_2D[i][j] = (diel[1][i + 1][j]) * lattice * epsilon0 * omega[j] * 1E-10;
			imag_polarizability_2D[i][j] = (1 - diel[0][i + 1][j]) * lattice * epsilon0 * omega[j] * 1E-10;
			denominator[j] = pow((1 + real_polarizability_2D[i][j] / epsilon0 / c0 * 0.5), 2) + pow(0.5 * imag_polarizability_2D[i][j] / epsilon0 / c0, 2);
			transmission[i][j] = 1.0 / denominator[j];
			absorb[i][j] = real_polarizability_2D[i][j] / epsilon0 / c0 / denominator[j];
			reflection[i][j] = 0.25 * (pow(real_polarizability_2D[i][j] / epsilon0 / c0, 2) + pow(imag_polarizability_2D[i][j] / epsilon0 / c0, 2)) / denominator[j];
		}
	}
	vector<vector<vector<double> > > ans{ real_polarizability_2D ,imag_polarizability_2D ,transmission ,absorb,reflection };
	for (int i = 0; i < fps.size(); i++)
	{
		for (int j = 0; j < n; j++)
			fprintf(fps[i], "   %.4lf     %.3lf     %.3lf\n", diel[0][0][j], i < 2 ? ans[i][0][j] / sigma0 : ans[i][0][j] * 100, i < 2 ? ans[i][1][j] / sigma0 : ans[i][1][j] * 100);
		WriteFile(filename[i].c_str());
		fclose(fps[i]);
	}
}

void get_linear_optical_spectrums_3d(DIELECT diel)
{
	int n = diel[0][0].size();
	vector<vector<double> >absorb(6, vector<double>(n));
	vector<vector<double> >refractive(6, vector<double>(n));
	vector<vector<double> >energylossspectrum(6, vector<double>(n));
	vector<vector<double> >extinction(6, vector<double>(n));
	vector<vector<double> >reflectivity(6, vector<double>(n));
	vector<double> freq(n);
	for (int i = 0; i < n; i++)
		freq[i] = diel[0][0][i] / planck_constant;
	vector<string> filename = { "absorption","refractive","energy_lossspectrum","extinction","reflectivity" };
	vector<FILE*> fps(filename.size());
	for (int i = 0; i < fps.size(); i++)
	{
		fps[i] = fopen(filename[i].c_str(), "w");
		i == 0 ? fprintf(fps[i], "	#Energy      xx(cm^-1)       yy(cm^-1)       zz(cm^-1)       xy(cm^-1)       yz(cm^-1)       zx(cm^-1)\n") :
			fprintf(fps[i], "	#Energy          xx             yy              zz              xy              yz              zx\n");
	}
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < n; j++)
		{
#define imag_epsilon diel[1][i + 1][j]
#define real_epsilon diel[0][i + 1][j]
			absorb[i][j] = 2 * pi * sqrt(2.0) * freq[j] * (sqrt(-real_epsilon + sqrt(imag_epsilon * imag_epsilon + real_epsilon * real_epsilon))) / (c0 * 100);
			refractive[i][j] = (sqrt(real_epsilon + sqrt(imag_epsilon * imag_epsilon + real_epsilon * real_epsilon))) / sqrt(2.0);
			energylossspectrum[i][j] = imag_epsilon / (imag_epsilon * imag_epsilon + real_epsilon * real_epsilon);
			extinction[i][j] = (sqrt(-real_epsilon + sqrt(imag_epsilon * imag_epsilon + real_epsilon * real_epsilon))) / sqrt(2.0);
			reflectivity[i][j] = ((refractive[i][j] - 1) * (refractive[i][j] - 1) + extinction[i][j] * extinction[i][j]) /
				((refractive[i][j] + 1) * (refractive[i][j] + 1) + extinction[i][j] * extinction[i][j]);
		}
	}
	vector<vector<vector<double> > > ans{ absorb,refractive ,energylossspectrum ,extinction,reflectivity };
	for (int i = 0; i < fps.size(); i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < 6; k++)
				fprintf(fps[i], "	%lf", ans[i][k][j]);
			fprintf(fps[i], "\n");
		}
		WriteFile(filename[i].c_str());
		fclose(fps[i]);
	}
}
double get_energy()
{
	double energy = 0;
	FILE* fp_ = fopen("OUTCAR", "r");
	if (fp_ == NULL)
	{
		printf("OUTCAR IS NOT EXIST!\n");
		return energy;
	}
	char buf[1024];
	while (fgets(buf, 1024, fp_) != NULL)
	{
		if (strstr(buf, "energy(sigma->0)") != NULL)
		{
			sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%lf", &energy);
		}
	}
	fclose(fp_);
	return energy;
}
double get_energy_oszi()
{
	double energy_oszi = 0;
	FILE* fp_ = fopen("OSZICAR", "r");
	if (fp_ == NULL)
	{
		printf("OSZICAR IS NOT EXIST!\n");
		return energy_oszi;
	}
	char buf[1024];
	while (fgets(buf, 1024, fp_) != NULL)
	{
		if (strstr(buf, "F") != NULL)
		{
			sscanf(buf, "%*s%*s%*s%*s%lf", &energy_oszi);
		}
	}
	fclose(fp_);
	return energy_oszi;
}
vector<double> get_stress()
{
	FILE* fp = fopen("OUTCAR", "r");
	if (fp == NULL)
	{
		printf("OUTCAR IS NOT EXIST!\n");
		return {};
	}
	vector<double> ans(6);//11          22          33          23          13          12
	char buf[1024];
	while (fgets(buf, 1024, fp) != NULL)
	{
		if (strstr(buf, "in kB") != NULL)
		{
			//XX          YY          ZZ          XY          YZ          ZX
			sscanf(buf, "%*s%*s%lf%lf%lf%lf%lf%lf", &ans[0], &ans[1], &ans[2], &ans[5], &ans[3], &ans[4]);
			break;
		}
	}
	fclose(fp);
	return ans;
}
vector<string> _outcar::str_split(string s, vector<char> sp)
{
	vector<string> ans;
	bool flag = 0;
	string str;
	for (int i = 0; i < s.size(); i++)
	{
		if (find(sp.begin(), sp.end(), s[i]) != sp.end())
		{
			flag = 1;
			continue;
		}
		else
		{
			if (flag == 1)
			{
				ans.push_back(str);
				str.clear();
				flag = 0;
			}
			str += s[i];
		}
	}
	if (!str.empty() && find(sp.begin(), sp.end(), str[0]) == sp.end())
		ans.push_back(str);
	return ans;
}
void _outcar::print_normal(FILE* fp, Para_Info* info) {
	{
		fprintf(fp, "#data	%s\n", info->name.c_str());
		fprintf(fp, "%s\n", m[info->key].c_str());
	}
};
void _outcar::print_null(FILE* fp, Para_Info* info) {};
void _outcar::print_GGA(FILE* fp, Para_Info* info) {
	{
		fprintf(fp, "#data	%s\n", info->name.c_str());
		if (m[info->key] == "--")
			m["GGA"] = m["LEXCH"];
		fprintf(fp, "%s\n", m["GGA"].c_str());
	}
}
void _outcar::print_ISPIN(FILE* fp, Para_Info* info) {
	{
		fprintf(fp, "#data	%s\n", info->name.c_str());
		if (m[info->key] == "1")
			fprintf(fp, "F\n");
		else if (m[info->key] == "2")
			fprintf(fp, "T\n");
	}
}
void _outcar::print_ISIF(FILE* fp, Para_Info* info)
{
	{
		fprintf(fp, "#data	%s\n", info->name.c_str());
		unordered_map<string, string > m_isif
		{
			{"0","atom position"},
			{"1","atom position"},
			{"2","atom position"},
			{"3","atom position&cell shape&cell volume"},
			{"4","atom position&cell shape"},
			{"5","cell shape"},
			{"6","cell shape&cell volume"},
			{"7","cell volume"},
		};
		fprintf(fp, "%s\n", m_isif[m[info->key]].c_str());
	}
}
void _outcar::print_else(FILE* fp)
{
	string val;
	if (m["ISTART"] == "1" && m["ICHARG"] == "11" && m["NSW"] == "0")
		val = "F";
	else
		val = "T";
	{
		fprintf(fp, "#data	Self_consistent\n");
		fprintf(fp, "%s\n", val.c_str());
	}
	val.clear();
	if (m["GGA"] == "--")
		m["GGA"] = m["LEXCH"];
	if (m["LHFCALC"] == "T" && m["GGA"] == "PE" && m["HFSCREEN"] == "0.2")
		val = "HSE06";
	else if (m["LHFCALC"] == "T" && m["GGA"] == "PE" && m["HFSCREEN"] == "0.3")
		val = "HSE03";
	else if (m["LHFCALC"] == "T" && m["GGA"] == "PS" && m["HFSCREEN"] == "0.2")
		val = "HSEsol";
	else if (m["LHFCALC"] == "T" && m["GGA"] == "PE")
		val = "PBE0";
	else if (m["LHFCALC"] == "T" && m["GGA"] == "B3")
		val = "B3LYP-VWN3";
	else if (m["LHFCALC"] == "T" && m["GGA"] == "B5")
		val = "B3LYP-VWN5";
	else if (m["LHFCALC"] == "T" && m["GGA"] == "LIBXC" && m["LIBXC1"] == "HYB_GGA_XC_B3PW91")
		val = "B3PW91";
	else if (m["LHFCALC"] == "T" && m["GGA"] == "LIBXC" && m["LIBXC1"] == "HYB_GGA_XC_B1WC")
		val = "B1-WC";
	else if (m["LHFCALC"] == "T" && m["METAGGA"] == "SCAN")
		val = "B3LYP-VWN5";
	else if (m["LHFCALC"] == "T" && m["AEXX"] == "1.0")
		val = "Hartree-Fock";
	{
		fprintf(fp, "#data	Hybrid_type\n");
		fprintf(fp, "%s\n", val.c_str());
	}
	val.clear();
	if (m["LDAUTYPE"] == "1")
		val = "Liechtenstein";
	else if (m["LDAUTYPE"] == "2")
		val = "Dudarev";
	else if (m["LDAUTYPE"] == "4")
		val = "Liechtenstein(no exchange splitting)";
	{
		fprintf(fp, "#data	DFT_U_type\n");
		fprintf(fp, "%s\n", val.c_str());
	}
	val.clear();
	if (m["IVDW"] == "1" || m["IVDW"] == "10")
		val = "DFT-D2";
	else if (m["IVDW"] == "11")
		val = "DFT-D3";
	else if (m["IVDW"] == "12")
		val = "DFT-D3";
	else if (m["IVDW"] == "13")
		val = "DFT-D4";
	else if (m["IVDW"] == "2" || m["IVDW"] == "20")
		val = "TS method";
	else if (m["IVDW"] == "21")
		val = "TS method(iterative Hirshfeld)";
	else if (m["IVDW"] == "202")
		val = "MBD-rSC";
	else if (m["IVDW"] == "263")
		val = "MBD-rSC/FI";
	else if (m["IVDW"] == "4")
		val = "dDsC";
	else if (m["IVDW"] == "3")
		val = "DFT-ulg";
	if (m["LUSE_VDW"] == "T" && m["GGA"] == "RE")
		val = "vdW-DF1";
	else if (m["LUSE_VDW"] == "T" && m["GGA"] == "ML")
		val = "vdW-DF2";
	else if (m["LUSE_VDW"] == "T" && m["GGA"] == "OR")
		val = "optPBE-vdW";
	else if (m["LUSE_VDW"] == "T" && m["GGA"] == "BO")
		val = "optB88-vdW";
	else if (m["LUSE_VDW"] == "T" && m["GGA"] == "MK")
		val = "optB86b-vdW";
	else if (m["LUSE_VDW"] == "T" && m["GGA"] == "MK" && m["ZAB_VDW"] == "-1.8867")
		val = "rev-vdW-DF2";
	else if (m["LUSE_VDW"] == "T" && m["GGA"] == "CX")
		val = "vdW-DF-cx";
	else if (m["LUSE_VDW"] == "T" && m["GGA"] == "ML")
		val = "rVV10";
	else if (m["LUSE_VDW"] == "T" && m["METAGGA"] == "SCAN")
		val = "SCAN+rVV10";
	else if (m["LUSE_VDW"] == "T" && m["METAGGA"] == "R2SCAN")
		val = "r2SCAN + rVV10";
	else if (m["LVDW"] == "T" && !m.count("IVDW"))
		val = "DFT-D2";
	{
		fprintf(fp, "#data	VDW_D_type\n");
		fprintf(fp, "%s\n", val.c_str());
	}
}
void _outcar::print_text(FILE* fp, const char name[], vector<string> v)
{
	if (v.empty())
		return;
	fprintf(fp, "#data	%s\n", name);
	for (int i = 0; i < v.size(); i++)
		fprintf(fp, "%s", v[i].c_str());
}
void _outcar::OUTCAR::DB_GetInfo(Para_Info para, const char input_file[])
{
	FILE* fp = fopen(input_file, "r");
	string value;
	if (fp == NULL)
	{
		printf("%s in not exist!\n", input_file);
		return;
	}
	char buf[1024];
	int iat = 0;
	int cnt = 0;
	while (fgets(buf, 1024, fp))
	{
		if (strstr(buf, "free  energy   TOTEN"))
		{
			this->energy.clear();
			cnt = 3;
			while (cnt--)
			{
				this->energy.push_back(buf);
				fgets(buf, 1024, fp);
			}
			continue;
		}
		if (strstr(buf, "Free energy of the ion-electron system (eV)"))
		{
			this->par_energy.clear();
			while (!strstr(buf, "free energy"))
			{
				this->par_energy.push_back(buf);
				fgets(buf, 1024, fp);
			}
			continue;
		}
		if (strstr(buf, "in kB"))
		{
			this->stress.clear();
			cnt = 2;
			while (cnt--)
			{
				this->stress.push_back(buf);
				fgets(buf, 1024, fp);
			}
			continue;
		}
		if (strstr(buf, "TOTAL-FORCE"))
		{
			this->force.clear();
			while (1)
			{
				this->force.push_back(buf);
				if (strstr(buf, "total drift"))
					break;
				fgets(buf, 1024, fp);
			}
			continue;
		}
		if (strstr(buf, "total charge") && !strstr(buf, "total charge-"))
		{
			this->charge.clear();
			while (1)
			{
				this->charge.push_back(buf);
				if (strstr(buf, "tot") && !strstr(buf, "#"))
					break;
				fgets(buf, 1024, fp);
			}
			continue;
		}
		if (strstr(buf, "magnetization (x)"))
		{
			this->magnitization_x.clear();
			while (1)
			{
				this->magnitization_x.push_back(buf);
				if (strstr(buf, "tot") && !strstr(buf, "#"))
					break;
				fgets(buf, 1024, fp);
			}
			continue;
		}
		if (strstr(buf, "magnetization (y)"))
		{
			this->magnitization_y.clear();
			while (1)
			{
				this->magnitization_y.push_back(buf);
				if (strstr(buf, "tot") && !strstr(buf, "#"))
					break;
				fgets(buf, 1024, fp);
			}
			continue;
		}
		if (strstr(buf, "magnetization (z)"))
		{
			this->magnitization_z.clear();
			while (1)
			{
				this->magnitization_z.push_back(buf);
				if (strstr(buf, "tot") && !strstr(buf, "#"))
					break;
				fgets(buf, 1024, fp);
			}
			continue;
		}
		string tmp = buf;
		string s = tmp.substr(0, find(tmp.begin(), tmp.end(), para.split[0]) - tmp.begin());
		vector<string> t = str_split(s, { ' ','\t','\r','\n' });
		if (t.empty())
			continue;
		else if (t.back() == para.key)
		{
			int cnt = 0;
			char* p = strrchr(buf, para.split[0]) + 1;
			char* token = strtok(p, " ");
			while (token != NULL && cnt < para.num)
			{
				value += token;
				token = strtok(NULL, " ");
				cnt++;
			}
			if (value.back() == '\n')
				value.pop_back();
			m[t.back()] = value;
			break;
		}
	}
	fclose(fp);
}
void _outcar::OUTCAR::output(const char output[], std::string flag, vector<string> file_name)
{
	FILE* fp = fopen(output, "at+");
	//basic parameter from input file
	auto cpfile = [fp](const char source[]) {
		FILE* fps = fopen(source, "r");
		if (fps == nullptr)
			return;
		char buf[1024];
		while (fgets(buf, 1024, fps))
			fputs(buf, fp);
		fclose(fps);
	};
	fprintf(fp, "#start\n");
	fprintf(fp, "#data	STRUCT\n");
	cpfile("CONTCAR");
	fprintf(fp, "#data	CONTCAR\n");
	cpfile("CONTCAR");
	fprintf(fp, "#data	POSCAR\n");
	cpfile("POSCAR");
	fprintf(fp, "#data	INCAR\n");
	cpfile("INCAR");
	fprintf(fp, "#data	KPOINTS\n");
	cpfile("KPOINTS");
	if (flag == "-a")
	{
		fprintf(fp, "#data	OSZICAR\n");
		cpfile("OSZICAR");
		fprintf(fp, "#data	OUTCAR\n");
		cpfile("OUTCAR");
	}
	if (file_name.size() > 0)
	{
		for (int i = 0; i < file_name.size(); i++)
		{
			ifstream file(file_name[i].c_str());
			if (!file) 
			{
        		cout << file_name[i] << " does not exist or failed to open!" << endl;
			}
    		else if (file.peek() == std::ifstream::traits_type::eof()) 
			{
        		cout << file_name[i] << " exists but is empty!" << endl;
    		}
			else
			{
				fprintf(fp, "#data	%s\n", file_name[i].c_str());
				cpfile(file_name[i].c_str());
			}
		}
	}
	for (int i = 0; i < parameter.size(); i++)
		parameter[i].print(fp, &parameter[i]);
	print_else(fp);
	//add par_energy to energy
	energy.push_back("\n");
	for (auto &str : par_energy)
		energy.emplace_back(str);
	print_text(fp, "Energy_in_OUTCAR", energy);
	print_text(fp, "Stress_in_OUTCAR", stress);
	print_text(fp, "Force_in_OUTCAR", force);
	print_text(fp, "Charge_in_OUTCAR", charge);
	vector<string> magnitization;
	for (int i = 0; i < magnitization_x.size(); i++)
		magnitization.push_back(magnitization_x[i]);
	for (int i = 0; i < magnitization_y.size(); i++)
		magnitization.push_back(magnitization_y[i]);
	for (int i = 0; i < magnitization_z.size(); i++)
		magnitization.push_back(magnitization_z[i]);
	print_text(fp, "Magnitization_in_OUTCAR", magnitization);
	fprintf(fp, "#end\n");
	fclose(fp);
	printf("Written %s file!\n", output);
}
static void TranMapToJson(unordered_map<string, vector<string> >& data, cJSON* obj)
{
	for (auto it = data.begin(); it != data.end(); it++)
	{
		string key = it->first;
		vector<string> val = it->second;
		if (val.empty())
		{
			cJSON_AddStringToObject(obj, key.c_str(), "");
			continue;
		}
		string str;
		for (string s : val)
			str += s;
		if (val.size() > 1)
			cJSON_AddStringToObject(obj, key.c_str(), str.c_str());
		else
		{
			vector<string> sp = str_split(val[0], { ' ','\t','\r','\n' });
			assert(!sp.empty());
			if (sp.size() > 1)
				cJSON_AddStringToObject(obj, key.c_str(), str.c_str());
			else
			{
				string strf = sp[0];
				/*if (strf.size() == 1 && strf[0] == 'T')
					cJSON_AddTrueToObject(obj, key.c_str());
				else if (strf.size() == 1 && strf[0] == 'F')
					cJSON_AddFalseToObject(obj, key.c_str());
				else*/
				{
					bool flag = 0;
					for (int i = 0; i < strf.size(); i++)
					{
						if (strf[i] != '.' && !isdigit(strf[i]))
						{
							cJSON_AddStringToObject(obj, key.c_str(), str.c_str());
							flag = 1;
							break;
						}
					}
					if (!flag)
						cJSON_AddNumberToObject(obj, key.c_str(), atof(str.c_str()));
				}
			}
		}
	}
};
void TranLogdataToJson(const char inputfile[], const char outputfile[], int MaxLogNum)
{
	FILE* fp = fopen(inputfile, "r");
	assert(fp != 0);
	char buf[1024];
	cJSON* root = cJSON_CreateObject();
	int curlog = 0;
	unordered_map<string, vector<string> > data;
	string keyword;
	vector<string> value;
	while (fgets(buf, 1024, fp))
	{
		if (strstr(buf, "#end"))
		{
			data[keyword] = value;
			cJSON* obj = cJSON_CreateObject();
			TranMapToJson(data, obj);
			string name = "Log" + to_string(++curlog);
			cJSON_AddItemToObject(root, name.c_str(), obj);
			data.clear();
			if (curlog > MaxLogNum)
				break;
			continue;
		}
		if (strstr(buf, "#data"))
		{
			if (!keyword.empty())
				data[keyword] = value;
			value.clear();
			keyword = str_split(buf, { ' ','\t','\r','\n' }).back();
			continue;
		}
		if (strspn(buf, "\t\n\r") == strlen(buf))
			continue;
		value.push_back(buf);
	}
	fclose(fp);
	auto output = cJSON_Print(root);
	FILE* fp_ = fopen(outputfile, "w");
	assert(fp_ != 0);
	fprintf(fp_, output);
	fclose(fp_);
	cJSON_Delete(root);
	printf("Written %s file!\n", outputfile);
};

void ParseJsonToMap(const cJSON* logItem, unordered_map<string, vector<string>>& data) 
{
    cJSON* child = logItem->child;
    while (child) 
	{
        if (cJSON_IsArray(child)) 
		{
            vector<string> values;
            cJSON* arrayItem = child->child;
            while (arrayItem) 
			{
                if (cJSON_IsString(arrayItem))
				{
                    values.push_back(arrayItem->valuestring);
                }
                arrayItem = arrayItem->next;
            }
            data[child->string] = values;
        } 
		else if (cJSON_IsString(child)) 
		{
            data[child->string] = {child->valuestring};
        }
        child = child->next;
    }
}

void TranJsonToLogdata(const char inputfile[], const char outputfile[], int MaxLogNum) 
{
	if(VM_fileexist(inputfile))
	{
		FILE* fp = fopen(inputfile, "r");
		FILE* fp_out = fopen(outputfile, "w");
		fseek(fp, 0, SEEK_END);
		long fsize = ftell(fp);
		fseek(fp, 0, SEEK_SET);
		char* buffer = new char[fsize + 1];
		fread(buffer, 1, fsize, fp);
		fclose(fp);
		buffer[fsize] = '\0';
		cJSON* root = cJSON_Parse(buffer);
		delete[] buffer;
		int curLog = 0;
        cJSON* logItem = nullptr;
        cJSON_ArrayForEach(logItem, root) 
		{
            if (++curLog > MaxLogNum)
                break;
            fprintf(fp_out, "#start\n");
            unordered_map<string, vector<string>> data;
            ParseJsonToMap(logItem, data);
            fprintf(fp_out, "#data %s\n", logItem->string);
            for (const auto& entry : data) 
			{
                fprintf(fp_out, "#data %s\n", entry.first.c_str());
                for (const auto& line : entry.second)
                    fprintf(fp_out, "%s\n", line.c_str());
            }
            fprintf(fp_out, "#end\n");
        }
		fclose(fp_out);
		cJSON_Delete(root);
	}
	else
	{
		printf("No %s file is found!", inputfile);
		return;
	}
	printf("Written %s file!\n", outputfile);
}

void TranFileToCsv(const char inputfile[], const char outputfile[])
{
	if(VM_fileexist(inputfile))
	{
		FILE* fp_r = fopen(inputfile, "r");
		FILE* fp_w = fopen(outputfile, "w");
		char buf[1024];
		while (fgets(buf, 1024, fp_r) != NULL)
		{
			bool wasSpace = false;
			string result;
			string buf_str = string(buf);
			buf_str.erase(0, buf_str.find_first_not_of(' '));
    		buf_str.erase(buf_str.find_last_not_of(' ') + 1);
			buf_str.erase(0, buf_str.find_first_not_of('	'));
    		buf_str.erase(buf_str.find_last_not_of('	') + 1);
			for (char &ch : buf_str) 
			{
				if (ch == ' ' || ch == '	') 
				{
					if (!wasSpace) 
					{
						result += ',';
						wasSpace = true;
					}
				} 
				else 
				{
					result += ch;
					wasSpace = false;
				}
			}
			fprintf(fp_w, "%s", result.c_str());
		}
		fclose(fp_r);
		fclose(fp_w);
		WriteFile(outputfile);
	}
	else
	{
		printf("No %s file is found!", inputfile);
		return;
	}
}

void TranCsvToFile(const char inputfile[], const char outputfile[])
{
	if(VM_fileexist(inputfile))
	{
		FILE* fp_r = fopen(inputfile, "r");
		FILE* fp_w = fopen(outputfile, "w");
		char buf[1024];
		while (fgets(buf, 1024, fp_r) != NULL)
		{
			for (char &ch : buf) 
				if (ch == ',')
					ch = ' ';
			fprintf(fp_w, "%s", buf);
		}
		fclose(fp_r);
		fclose(fp_w);
		WriteFile(outputfile);
	}
	else
	{
		printf("No %s file is found!", inputfile);
		return;
	}
}

void PlusLog(vector<string> filename, const char outputfile[])
{
	FILE* fp_w = fopen(outputfile, "a");
	for(int i = 0; i < filename.size(); i++)
	{
		if(VM_fileexist(filename[i].c_str()))
		{
			FILE* fp_r = fopen(filename[i].c_str(), "r");
			char buf[1024];
			while (fgets(buf, 1024, fp_r) != NULL)
			{
				fprintf(fp_w, "%s", buf);
			}
			fclose(fp_r);
		}
		else
		{
			printf("No %s file is found!", filename[i].c_str());
			continue;
		}
	}
	fclose(fp_w);
	printf("Add to %s file!\n", outputfile);
}

void PlusJson(vector<string> filename, const char outputfile[])
{
    FILE* fp_w = fopen(outputfile, "w");
	cJSON* merged_json = cJSON_CreateObject();
    for (int i = 0; i < filename.size(); i++) 
	{
		if(VM_fileexist(filename[i].c_str()))
		{
			FILE* file = fopen(filename[i].c_str(), "r");
			fseek(file, 0, SEEK_END);
			long length = ftell(file);
			fseek(file, 0, SEEK_SET);
			char* data = (char*)malloc(length + 1);
			fread(data, 1, length, file);
			data[length] = '\0';
			fclose(file);
			cJSON* json = cJSON_Parse(data);
			free(data);
			cJSON* current_element = NULL;
			cJSON_ArrayForEach(current_element, json) 
			{
				cJSON_AddItemToObject(merged_json, current_element->string, cJSON_Duplicate(current_element, 1));
			}
			cJSON_Delete(json);
		}
		else
		{
			printf("No %s file is found!\n", filename[i].c_str());
			continue;
		}
    }
    char *string = cJSON_Print(merged_json);
    fprintf(fp_w, "%s\n", string);
    free(string);
	cJSON_Delete(merged_json);
	printf("Written %s file!\n", outputfile);
}

tuple<vector<int>, vector<double>, int> EnergyConvergence(double EDIFFG)
{
    if (EDIFFG <= 0)
        return {};
    FILE *fp = fopen("OSZICAR", "r");
    vector<double> energy;
    vector<int> idx;
    char buf[1024];
    if (fp == nullptr)
    {
        printf("No OSZICAR! VASPMATE will read energy from OUTCAR!\n");
        fp = fopen("OUTCAR", "r");
        if (fp == nullptr)
        {
            printf("No OUTCAR! VASPMATE can't judge energy convergence!\n");
            return {};
        }
        double f;
        int i = 1;
        while (fgets(buf, 1024, fp) != NULL)
        {
            if (strstr(buf, "energy(sigma->0)") != NULL)
            {
                sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%lf", &f);
                energy.push_back(f);
                idx.push_back(i++);
            }
        }
    }
    else
    {
        int pos = 4;
        while (fgets(buf, 1024, fp))
        {
            if (strstr(buf, "E0"))
            {
                double f;
                int i;
                sscanf(buf, "%d%*s%*s%*s%lf", &i, &f);
                idx.push_back(i);
                energy.push_back(f);
            }
        }
        fclose(fp);
    }
    int flag = 0;
	FILE *fp_con = fopen("Converge.data", "a");
    if (fabs(energy[energy.size() - 1] - energy[energy.size() - 2]) <= EDIFFG)
    {
        flag = 1;
        printf("The calculation energy result has successfully converged!\n");
		fprintf(fp_con,"1 >> Energy Converged!\n");
    }
    else
    {
        flag = -1;
        printf("The calculation energy result did not converge successfully!\n");
		fprintf(fp_con,"3 >> Energy not Converged!\n");
    }
	fclose(fp_con);
    return make_tuple(idx, energy, flag);
}

tuple<vector<int>, vector<double>, int> ForceConvergence(double EDIFFG)
{
    if (EDIFFG >= 0)
        return {};
    FILE *fp = fopen("OUTCAR", "r");
	FILE *fp_ = fopen("POSCAR", "r");
	POSCAR pos;
    char buf[1024];
    vector<double> fmax_list;
    vector<int> idx;
    int cnt = 1;
    while (fgets(buf, 1024, fp))
    {
        if (strstr(buf, "TOTAL-FORCE"))
        {
            vector<vector<double>> f;
            fgets(buf, 1024, fp);
            fgets(buf, 1024, fp);
            double a, b, c;
            sscanf(buf, "%*lf%*lf%*lf%lf%lf%lf", &a, &b, &c);
			readposcar(fp_, pos);
			if (pos.ifix == 0)
			{
					f.push_back({0, 0, 0});
			}
			else if (pos.ifix == 1)
			{
				if ( pos.fix[cnt - 1][0]=='T' && pos.fix[cnt - 1][1]=='T' && pos.fix[cnt - 1][2]=='T') 
					f.push_back({a, b, c});
				else
					f.push_back({0, 0, 0});
			}
            double fmax = 0;
            for (auto i : f)
                fmax = max(fmax, sqrt((i[0] * i[0] + i[1] * i[1] + i[2] * i[2])));
            fmax_list.push_back(fmax);
            idx.push_back(cnt++);
        }
    }
    int flag = 0;
	FILE *fp_con = fopen("Converge.data", "a");
    if (fmax_list.back() <= fabs(EDIFFG))
    {
        flag = 1;
        printf("The calculation force result has successfully converged!\n");
		fprintf(fp_con,"2 >>Force Converged!\n");
    }
    else
    {
        flag = -1;
        printf("The calculation force result did not converge successfully!\n");
		fprintf(fp_con,"4 >> Force not Converged!\n");
    }
	fclose(fp_con);
    return make_tuple(idx, fmax_list, flag);
}

void RelaxJudgeConvergence()
{
    FILE *fp = fopen("OUTCAR", "r");
    if (fp == nullptr)
    {
        printf("No OUTCAR! VASPMATE can't read EDIFFG parameter!\n");
        return;
    }
    char buf[1024];
    double EDIFFG;
    while (fgets(buf, 1024, fp))
    {
        if (strstr(buf, "EDIFFG"))
        {
            sscanf(buf, "%*s%*s%lf", &EDIFFG);
            break;
        }
    }
    fclose(fp);
    printf("EDIFFG = %.5e in this calculation!\n", EDIFFG);
    if(EDIFFG > 0)
    {
        printf("EDIFFG > 0, VASPMATE will calculate energy convergence!\n");
        auto res = EnergyConvergence(EDIFFG);
        Plotcurve(get<0> (res), get<1> (res), string("Step"), string("Energy"), string("Energy_Convergence"), string("Energy_Convergence.spco"));
        printf("Written \"Energy_Convergence.spco\" file!\n");
    }
    else
    {
        printf("EDIFFG < 0, VASPMATE will calculate force convergence!\n");
        auto res = ForceConvergence(EDIFFG);
        Plotcurve(get<0> (res), get<1> (res), "Step", "Force", "Force_Convergence", "Force_Convergence.spco");
        printf("Written \"Force_Convergence.spco\" file!\n");
    }
}

tuple<vector<int>, vector<double>, int> StaticEnergyConvergence(double EDIFF)
{
    char buf[1024];
    vector<double> energy;
    vector<int> idx;
    int i = 1;
    FILE* fp = fopen("OSZICAR","r");
    if(fp == nullptr)
    {
        printf("No OSZICAR! VASPMATE will read energy from OUTCAR!\n");
        fp = fopen("OUTCAR", "r");
        if (fp == nullptr)
        {
            printf("No OUTCAR! VASPMATE can't judge energy convergence!\n");
            return {};
        }
        double f;
        int i = 1;
        while (fgets(buf, 1024, fp) != NULL)
        {
            if (strstr(buf, "energy(sigma->0)") != NULL)
            {
                sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%lf", &f);
                energy.push_back(f);
                idx.push_back(i++);
            }
        }
    }
	fgets(buf, 1024, fp);
    while(fgets(buf, 1024, fp))
    {
        if(strstr("E0", buf)) 
			break;
		else
        {
            double f;
            sscanf(buf, "%*s%*s%lf", &f);
            energy.push_back(f);
            idx.push_back(i++);
        }
    }
    int flag = 0;
	FILE *fp_con = fopen("Converge.data", "a");
    if (fabs(energy[energy.size() - 1] - energy[energy.size() - 2]) <= EDIFF)
    {
        flag = 1;
        printf("The calculation energy result has successfully converged!\n");
		fprintf(fp_con,"The calculation energy result has successfully converged!\n");
    }
    else
    {
        flag = -1;
        printf("The calculation energy result did not converge successfully!\n");
		fprintf(fp_con,"The calculation energy result did not converge successfully!\n");
    }
    return make_tuple(idx, energy, flag);
	fclose(fp_con);
}

tuple<vector<int>, vector<double>, int> mdConvergence(double EDIFF)
{
    char buf[1024];
    vector<double> energy;
    vector<int> idx;
    int i = 1;
    FILE* fp = fopen("OSZICAR","r");
    if(fp == nullptr)
    {
        printf("No OSZICAR! VASPMATE will read energy from OUTCAR!\n");
        fp = fopen("OUTCAR", "r");
        if (fp == nullptr)
        {
            printf("No OUTCAR! VASPMATE can't judge MD convergence!\n");
            return {};
        }
        double f;
        int i = 1;
        while (fgets(buf, 1024, fp) != NULL)
        {
            if (strstr(buf, "energy(sigma->0)") != NULL)
            {
                sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%lf", &f);
                energy.push_back(f);
                idx.push_back(i++);
            }
        }
    }
	fgets(buf, 1024, fp);
    while(fgets(buf, 1024, fp))
    {
        if(strstr("E0", buf)) 
			break;
		else
        {
            double f;
            sscanf(buf, "%*s%*s%lf", &f);
            energy.push_back(f);
            idx.push_back(i++);
        }
    }
    int flag = 0;
    return make_tuple(idx, energy, flag);
}

void StaticJudgeConvergence()
{
    FILE *fp = fopen("OUTCAR", "r");
    if (fp == nullptr)
    {
        printf("No OUTCAR! VASPMATE can't read EDIFF parameter!\n");
        return;
    }
    char buf[1024];
    double EDIFF;
    while (fgets(buf, 1024, fp))
    {
        if (strstr(buf, "EDIFF"))
        {
            sscanf(buf, "%*s%*s%lf", &EDIFF);
            break;
        }
    }
    fclose(fp);
	FILE *fp_con = fopen("Converge.data", "a");
    printf("EDIFF = %.5e in this calculation!\n", EDIFF);
	fprintf(fp_con,"EDIFF = %.5e in this calculation!\n", EDIFF);
    auto res = StaticEnergyConvergence(EDIFF);
    Plotcurve(get<0> (res), get<1> (res), "Step", "Energy", "Energy_Convergence", "Energy_Convergence.spco");
    printf("Written \"Energy_Convergence.spco\" file!\n");
	fclose(fp_con);
}

void mdJudgeConvergence()
{
    FILE *fp = fopen("OUTCAR", "r");
    if (fp == nullptr)
    {
        printf("No OUTCAR! VASPMATE can't read EDIFF parameter!\n");
        return;
    }
    char buf[1024];
    double EDIFF;
    while (fgets(buf, 1024, fp))
    {
        if (strstr(buf, "EDIFF"))
        {
            sscanf(buf, "%*s%*s%lf", &EDIFF);
            break;
        }
    }
    fclose(fp);
	FILE *fp_con = fopen("Converge.data", "a");
    printf("EDIFF = %.5e in this calculation!\n", EDIFF);
	fprintf(fp_con,"EDIFF = %.5e in this calculation!\n", EDIFF);
    auto res = mdConvergence(EDIFF);
    Plotcurve(get<0> (res), get<1> (res), "Step", "Energy", "MD_Convergence", "MD_Convergence.spco");
    printf("Written \"MD_Convergence.spco\" file!\n");
	fclose(fp_con);
}

