#include"../include/cif.h"
#include"../include/tools.h"

static const double PI = 3.14159265358979;

using namespace std;
void CIF_init(CIF& cif)
{
	cif.length_a = cif.length_b = cif.length_c = 0;
	cif.alpha = cif.beta = cif.gamma = 0;
	cif.pos_atom_num = 0;
}

void TranCellToVec(CIF& cif)
{
	//  a // x ,b in xy
	cif.vec[0][0] = cif.length_a;  // a1
	cif.vec[0][1] = cif.vec[0][2] = 0; // a2
	cif.vec[1][0] = cif.length_b * cos(cif.gamma * PI / 180); // b1
	cif.vec[2][0] = cif.length_c * cos(cif.beta * PI / 180); // c1
	cif.vec[1][2] = 0; // b3
	cif.vec[1][1] = sqrt(cif.length_b * cif.length_b - cif.vec[1][0] * cif.vec[1][0]);// b2
	cif.vec[2][1] = (cif.length_b * cif.length_c * cos(cif.alpha * PI / 180) - cif.vec[1][0] * cif.vec[2][0]) / cif.vec[1][1]; //c2
	cif.vec[2][2] = sqrt(cif.length_c * cif.length_c - cif.vec[2][0] * cif.vec[2][0] - cif.vec[2][1] * cif.vec[2][1]);
}

template<typename T>
T TranEquivToCal(vector<string> input, string equiv)
{
	for (int i = 0; i < equiv.size(); i++)
	{
		if (equiv[i] == 'x')
			equiv.replace(i, 1, "(" + TranDecimalTofraction(input[0]) + ")");
		if (equiv[i] == 'y')
			equiv.replace(i, 1, "(" + TranDecimalTofraction(input[1]) + ")");
		if (equiv[i] == 'z')
			equiv.replace(i, 1, "(" + TranDecimalTofraction(input[2]) + ")");
	}
	auto back = midToBack(equiv);
	double res = getValue(back);
	return res;
}

void getpos_xyz(CIF& cif)
{
	for (int i = 0; i < cif.elemsymbol.size(); i++)
	{
		int current_atom = cif.equiv_pos.size() * cif.cif_typenum[i];
		vector<vector<double> >tmp(current_atom, vector<double>(3));
		int cur = 0;
		for (int j = 0; j < cif.cif_typenum[i]; j++)
		{
			for (int k = 0; k < cif.equiv_pos.size(); k++)
			{
				for (int m = 0; m < 3; m++)
				{
					double value = TranEquivToCal<double>(vector<string>{cif.frac_xyz[i][j][0],
						cif.frac_xyz[i][j][1], cif.frac_xyz[i][j][2]}, cif.equiv_pos[k][m]);
					while (value < 0) value++;
					while (value >= 1) value--;
					tmp[cif.equiv_pos.size() * j + k][m] = value;
				}
			}
		}
		vector<vector<double> >::iterator it = tmp.begin(), _it = tmp.begin();
		vector < vector<double>> rep;
		for (; it != tmp.end(); it++)
		{
			vector < vector<double>> tp(tmp.begin(), it);
			for (int i = 0; i < tp.size(); i++)
			{
				double dis = 0;
				for (int j = 0; j < 3; j++)
					dis += (tp[i][j] - (*it)[j]) * (tp[i][j] - (*it)[j]);
				if (fabs(dis) <= 1e-4)
				{
					rep.push_back(*it);
					break;
				}
			}
		}
		if (!rep.empty())
		{
			int k = 0;
			for (_it; _it != tmp.end();)
			{
				if (*_it == rep[k]) {
					_it = tmp.erase(_it);
					k++;
				}
				else {
					++_it;
				}
			}
		}
		cif.pos_xyz.push_back(tmp);
	}
	cif.pos_typenum.resize(cif.pos_xyz.size());
	for (int i = 0; i < cif.pos_xyz.size(); i++)
	{
		cif.pos_typenum[i] = cif.pos_xyz[i].size();
		cif.pos_atom_num += cif.pos_typenum[i];
	}
}

void readcif(FILE* fp, CIF& cif)
{
	int start = -1;
	int pos_elem = 0, pos_x = 0, pos_y = 0, pos_z = 0;
	char buf[1024];
	int flag = 0;
	int flag_num = 0;
	int cnt = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		if (start != -1 && (strstr(buf, "loop_") != NULL || strspn(buf, " \t\n\r") == strlen(buf) || strstr(buf, "#") != NULL|| strstr(buf, "bond") != NULL))
		{
			start = -1;
			break;
		}
		if (strstr(buf, "_cell_length_a") != NULL)
		{
			sscanf(buf, "%*s%lf", &cif.length_a);
			continue;
		}
		if (strstr(buf, "_cell_length_b") != NULL)
		{
			sscanf(buf, "%*s%lf", &cif.length_b);
			continue;
		}
		if (strstr(buf, "_cell_length_c") != NULL)
		{
			sscanf(buf, "%*s%lf", &cif.length_c);
			continue;
		}
		if (strstr(buf, "_cell_angle_alpha") != NULL)
		{
			sscanf(buf, "%*s%lf", &cif.alpha);
			continue;
		}
		if (strstr(buf, "_cell_angle_beta") != NULL)
		{
			sscanf(buf, "%*s%lf", &cif.beta);
			continue;
		}
		if (strstr(buf, "_cell_angle_gamma") != NULL)
		{
			sscanf(buf, "%*s%lf", &cif.gamma);
			continue;
		}
		if (strstr(buf, "loop_") != NULL)
		{
			flag = 0;
			continue;
		}
		if (strstr(buf, "_symmetry_equiv_pos_site_id") != NULL)
		{
			flag_num = 1;
			continue;
		}
		if (strstr(buf, "_symmetry_equiv_pos_as_xyz") != NULL)
		{
			flag = 1;
			continue;
		}
		if (flag)
		{
			string tmp = buf;
			string::iterator _it, it_;
			if (flag_num)
			{
				string::iterator it = tmp.begin();
				if (!isdigit(*it))
				{
					flag = 1;
					continue;
				}
				int _f = 0;
				for (; it != tmp.end(); it++)
				{
					if (!_f && isdigit(*it))
					{
						_it = it;
						_f = 1;
						continue;
					}
					if (_f && !isdigit(*it))
					{
						it_ = it;
						break;
					}
				}
				tmp.erase(_it, it_);
			}
			vector<string> s(3);
			for (int i = 0, j = 0; i < tmp.size(); i++)
			{
				if (tmp[i] == ' ' || tmp[i] == '\'' || tmp[i] == '\t' || tmp[i] == '\n' || tmp[i] == '\r')
					continue;
				else if (tmp[i] == ',')
					j++;
				else
					s[j].push_back(tmp[i]);
			}
			cif.equiv_pos.push_back(vector<string>{s[0], s[1], s[2]});
		}
		if (strstr(buf, "_atom_site") != NULL)
		{
			start++;
			if (strstr(buf, "_atom_site_label") != NULL)
				pos_elem = start;
			if (strstr(buf, "_atom_site_fract_x") != NULL)
				pos_x = start;
			if (strstr(buf, "_atom_site_fract_y") != NULL)
				pos_y = start;
			if (strstr(buf, "_atom_site_fract_z") != NULL)
				pos_z = start;
		}
		else if (start != -1)
		{
			vector<string> tmp;
			char* token = NULL;
			token = strtok(buf, " \t\n\r\f");
			while (token != NULL)
			{
				tmp.push_back(token);
				token = strtok(NULL, " \t\n\r\f");
			}
		/*	string s;
			for (int i = 0; i < strlen(buf); i++)
			{
				if (buf[i] != ' ')
					s.push_back(buf[i]);
				else if (!s.empty())
				{
					tmp.push_back(s);
					s.clear();
				}
			}*/
			string elem = RemovePostfixNum(tmp[pos_elem]);
			vector<string>::iterator _elem = find(cif.elemsymbol.begin(), cif.elemsymbol.end(), elem);
			if (_elem != cif.elemsymbol.end())
			{
				cif.cif_typenum[_elem - cif.elemsymbol.begin()]++;
				cif.frac_xyz[_elem - cif.elemsymbol.begin()].push_back(vector<string> {to_string(atof(tmp[pos_x].c_str())),
					to_string(atof(tmp[pos_y].c_str())), to_string(atof(tmp[pos_z].c_str()))});
			}
			else
			{
				cif.cif_typenum.push_back(1);
				cif.frac_xyz.resize(cif.frac_xyz.size() + 1);
				cif.frac_xyz[cnt++].push_back(vector<string> {to_string(atof(tmp[pos_x].c_str())),
					to_string(atof(tmp[pos_y].c_str())), to_string(atof(tmp[pos_z].c_str()))});
			}
			cif.elemsymbol.push_back(elem);
			vector<string>::iterator it, _it;
			for (it = ++cif.elemsymbol.begin(); it != cif.elemsymbol.end();)
			{
				_it = find(cif.elemsymbol.begin(), it, *it);
				if (_it != it)
				{
					it = cif.elemsymbol.erase(it);
				}
				else
					it++;
			}
		}
	}
	if(cif.equiv_pos.empty())
		cif.equiv_pos.push_back(vector<string>{"x","y","z"});
}

void printinfo(CIF cif)
{
	cout << "a= " << cif.length_a << endl;
	cout << "b= " << cif.length_b << endl;
	cout << "c= " << cif.length_c << endl;
	cout << "alpha= " << cif.alpha << endl;
	cout << "beta= " << cif.beta << endl;
	cout << "gamma= " << cif.gamma << endl;
	for (int i = 0; i < 3; i++)
		cout << cif.vec[i][0] << " " << cif.vec[i][1] << " " << cif.vec[i][2] << " " << endl;
	for (int i = 0; i < cif.elemsymbol.size(); i++)
		cout << cif.elemsymbol[i] << " ";
	cout << endl;
	for (int i = 0; i < cif.frac_xyz.size(); i++)
		for (int j = 0; j < cif.frac_xyz[i].size(); j++)
			cout << cif.frac_xyz[i][j][0] << " " << cif.frac_xyz[i][j][1] << " " << cif.frac_xyz[i][j][2] << " " << endl;
	for (int i = 0; i < cif.cif_typenum.size(); i++)
		cout << cif.cif_typenum[i] << " ";
	cout << endl;
	for (int i = 0; i < cif.equiv_pos.size(); i++)
		cout << cif.equiv_pos[i][0] << " " << cif.equiv_pos[i][1] << " " << cif.equiv_pos[i][2] << " " << endl;
	for (int i = 0; i < cif.pos_typenum.size(); i++)
		cout << cif.pos_typenum[i] << " ";
	cout << endl;
	for (int i = 0; i < cif.pos_xyz.size(); i++)
		for (int j = 0; j < cif.pos_xyz[i].size(); j++)
			cout << cif.pos_xyz[i][j][0] << " " << cif.pos_xyz[i][j][1] << " " << cif.pos_xyz[i][j][2] << " " << endl;
}

void TranCifToPOS(CIF cif, POSCAR& pos)
{
	strcpy(pos.title, "Generated by VASPMATE");
	pos.latt = 1.0;
	pos.ifix = 1;
	pos.iflg = 0;
	pos.nant[0] = cif.pos_atom_num;
	pos.nant[1] = cif.pos_typenum.size();
	for (int i = 0; i < cif.pos_typenum.size(); i++)
		pos.typenum[i] = cif.pos_typenum[i];
	for (int i = 0; i < cif.elemsymbol.size(); i++)
		strcpy(pos.elemsym[i], cif.elemsymbol[i].c_str());
	for (int i = 0; i < pos.nant[1]; i++)
	{
		int an;
		char at[3];
		strcpy(at, pos.elemsym[i]);
		getAtomNum(at, an);
		pos.elemnum[i] = an;
	}
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			pos.vec[i][j] = cif.vec[i][j];
	int cnt = 0;
	for (int i = 0; i < pos.nant[1]; i++)
	{
		for (int j = 0; j < pos.typenum[i]; j++)
			for (int k = 0; k < 3; k++)
				pos.xyz[cnt + j][k] = cif.pos_xyz[i][j][k];
		cnt += pos.typenum[i];
	}
}

void TranCIFToPOSCAR(const char _cif[],const char _poscar[])
{
	CIF cif;
	CIF_init(cif);
	FILE* fp = fopen(_cif, "r");
	if (fp == NULL)
	{
		printf("%s is not exist!\n", _cif);
		return;
	}
	readcif(fp, cif);
	//printinfo(cif);
	getpos_xyz(cif);
	TranCellToVec(cif);
	POSCAR pos;
	TranCifToPOS(cif, pos);
	FILE* fp2 = fopen(_poscar, "w");
	savposcar(fp2, pos);
	printf("Written %s file!\n", _poscar);
	fclose(fp);
	fclose(fp2);
}

void TranVecToBox(double vec[3][3], double info[6])
{
	info[0] = sqrt(vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]); // a
	info[1] = sqrt(vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]); // b
	info[2] = sqrt(vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]); // c
	info[3] = acos((vec[1][0] * vec[2][0] + vec[1][1] * vec[2][1] + vec[1][2] * vec[2][2]) / info[1] / info[2]) * 180 / PI; //alpha
	info[4] = acos((vec[0][0] * vec[2][0] + vec[0][1] * vec[2][1] + vec[0][2] * vec[2][2]) / info[0] / info[2]) * 180 / PI; //beta
	info[5] = acos((vec[1][0] * vec[0][0] + vec[1][1] * vec[0][1] + vec[1][2] * vec[0][2]) / info[0] / info[1]) * 180 / PI; //gamma
}

void TranPOSCARToCIF(const char _poscar[], const char _cif[])
{
	FILE* _fp = fopen(_poscar, "r");
	POSCAR pos;
	readposcar(_fp, pos);
	fclose(_fp);
	double info[6];
	TranVecToBox(pos.vec, info);
	FILE* fp = fopen(_cif, "w");
	fprintf(fp, "data_created_by_VASPMATE\n");
	fprintf(fp, "_audit_creation_date          %s %s\n",__DATE__,__TIME__);
	fprintf(fp, "_pd_phase_name                     'CIF files'\n");
	fprintf(fp, "_cell_length_a                     %lf\n", info[0]);
	fprintf(fp, "_cell_length_b                     %lf\n", info[1]);
	fprintf(fp, "_cell_length_c                     %lf\n", info[2]);
	fprintf(fp, "_cell_angle_alpha                  %lf\n", info[3]);
	fprintf(fp, "_cell_angle_beta                   %lf\n", info[4]);
	fprintf(fp, "_cell_angle_gamma                  %lf\n", info[5]);
	fprintf(fp, "_symmetry_space_group_name_H-M       'P 1'\n");
	fprintf(fp, "_symmetry_Int_Tables_number            1\n");
	fprintf(fp, "loop_\n");
	fprintf(fp, "_symmetry_equiv_pos_as_xyz\n");
	fprintf(fp, "    'x, y, z'\n");
	fprintf(fp, "loop_\n");
	fprintf(fp,	"   _atom_site_label\n");
	fprintf(fp,	"   _atom_site_occupancy\n");
	fprintf(fp,	"   _atom_site_fract_x\n");
	fprintf(fp,	"   _atom_site_fract_y\n");
	fprintf(fp,	"   _atom_site_fract_z\n");
	fprintf(fp,	"   _atom_site_type_symbol\n");
	int cnt = 0;
	for (int i = 0; i < pos.nant[1]; i++)
	{
		for (int j = 0; j < pos.typenum[i]; j++)
		{
			fprintf(fp, "	%2s%03d  1.0    %lf    %lf    %lf    %s\n", pos.elemsym[i], j, pos.xyz[cnt+j][0], pos.xyz[cnt + j][1], pos.xyz[cnt + j][2], pos.elemsym[i]);
		}
		cnt += pos.typenum[i];
	}
	fclose(fp);
	printf("Written %s file!\n", _cif);
}