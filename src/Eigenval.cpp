#include"../include/Eigenval.h"
using namespace std;
EIGENVAL::EIGENVAL(const char file[], double efermi)
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		exit(1);
	}
	char buf[1024];
	int line = 0;
	int kpt = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		line++;
		if (line == 1)
		{
			sscanf(buf, "%*d%*d%*d%d", &ISPIN);
			eigenvalue.resize(ISPIN);
			continue;
		}
		if ((line >= 2 && line <= 5) || line == 7)
			continue;
		if (line == 6)
		{
			sscanf(buf, "%d%d%d", &value_electron, &num_kpts, &num_bands);
			for (int i = 0; i < eigenvalue.size(); i++)
			{
				eigenvalue[i].resize(num_kpts);
				for (int j = 0; j < eigenvalue[i].size(); j++)
					eigenvalue[i][j].resize(num_bands);
			}
			continue;
		}
		if (line >= 8)
		{
			if (strspn(buf, " \t\n\r") == strlen(buf))
			{
				kpt++;
				continue;
			}
			else if (len(buf) == 4)
			{
				double temp[4];
				sscanf(buf, "%lf%lf%lf%lf", &temp[0], &temp[1], &temp[2], &temp[3]);
				position.push_back(vector<double>{temp[0], temp[1], temp[2]});
				weight.push_back(temp[3]);
			}
			else
			{
				int number;
				double value[2];
				sscanf(buf, "%d%lf%lf", &number, &value[0], &value[1]);
				for (int i = 0; i < ISPIN; i++)
					eigenvalue[i][kpt][number - 1] = value[i];
			}
		}
	}
	for (int k = 0; k < ISPIN; k++)
	{
		vector<double> maxvalue(num_bands, -1e10), minvalue(num_bands, 1e10);
		double lumo_value = 1e10, homo_value = -1e10;
		for (int i = 0; i < num_kpts; i++)
			for (int j = 0; j < num_bands; j++)
			{
				maxvalue[j] = max(maxvalue[j], eigenvalue[k][i][j]);
				minvalue[j] = min(minvalue[j], eigenvalue[k][i][j]);
			}
		for (int i = 0; i < num_bands; i++)
		{
			if (maxvalue[i] >= efermi && minvalue[i] <= efermi) // fermi surface
			{
				if (k == 0)
					band_index_up.insert(i + 1);
				else if (k == 1)
					band_index_dw.insert(i + 1);
			}
			else if (minvalue[i] >= efermi)
			{
				if (minvalue[i] <= lumo_value)
				{
					cbm = i + 1;
					lumo_value = minvalue[i];
				}
			}
			else if (maxvalue[i] <= efermi)
			{
				if (maxvalue[i] >= homo_value)
				{
					vbm = i + 1;
					homo_value = maxvalue[i];
				}
			}
		}
	}
	if (band_index_up.empty())
		band_index_up.insert(vbm);
	if (ISPIN == 2 && band_index_dw.empty())
		band_index_dw.insert(vbm);
}

vector<double> EIGENVAL::getkpointsline(double vec[3][3])
{
	vector<double> kline;
	double rec_vec[3][3];
	RecMat(vec, rec_vec);
	vector<vector<double> > carts_kpt(this->num_kpts, vector<double>(3));
	for (int i = 0; i < position.size(); i++)
	{
		double p[3] = { position[i][0],position[i][1],position[i][2] };
		for (int j = 0; j < 3; j++)
		{
			carts_kpt[i][j] = p[0] * rec_vec[0][j] + p[1] * rec_vec[1][j] + p[2] * rec_vec[2][j];
		}
	}
	kline.resize(this->num_kpts);
	kline[0] = 0;
	vector<double> delta;
	for (int i = 0; i < this->num_kpts - 1; i++)
	{
		double d = sqrt((carts_kpt[i][0] - carts_kpt[i + 1][0]) * (carts_kpt[i][0] - carts_kpt[i + 1][0]) +
			(carts_kpt[i][1] - carts_kpt[i + 1][1]) * (carts_kpt[i][1] - carts_kpt[i + 1][1]) +
			(carts_kpt[i][2] - carts_kpt[i + 1][2]) * (carts_kpt[i][2] - carts_kpt[i + 1][2]));
		delta.push_back(d);
	}
	for (int i = 1; i < kline.size(); i++)
		kline[i] = kline[i - 1] + delta[i - 1];
	return kline;
}

void EIGENVAL::TranEigenToXcrysden(set<int> select_band_index, double efermi)
{
	FILE* fp1 = fopen("KPOINTS", "r");
	char buf[1024];
	if (fp1 == NULL)
	{
		printf("KPOINTS IS NOT EXIST!\n");
		return;
	}
	fgets(buf, 1024, fp1);
	fclose(fp1);
	int mesh[3];
	char ikmesh;
	sscanf(buf, "%c%d%d%d", &ikmesh, &mesh[0], &mesh[1], &mesh[2]);
	FILE* fp2 = fopen("POSCAR", "r");
	if (fp2 == NULL)
	{
		printf("POSCAR IS NOT EXIST!\n");
		return;
	}
	POSCAR pos;
	readposcar(fp2, pos);
	fclose(fp2);
	vector<double> eigenval(mesh[0] * mesh[1] * mesh[2]);
	int grid_address[mesh[0] * mesh[1] * mesh[2]][3];
	int grid_mapping_table[mesh[0] * mesh[1] * mesh[2]];
	int is_shift[3];
	if (ikmesh == 'G')
	{
		is_shift[0] = 0;
		is_shift[1] = 0;
		is_shift[2] = 0;
	}
	else
	{
		is_shift[0] = 1;
		is_shift[1] = 1;
		is_shift[2] = 1;
	}
	int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
	translate_typenum_type(types, pos.nant[1], pos.typenum);
	double n_vec[3][3];
	transpose_matrix(pos.vec, n_vec);
	SpglibDataset* dataset;
	dataset = spg_get_dataset(n_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
	int num_ir = spg_get_ir_reciprocal_mesh(grid_address,
		grid_mapping_table,
		mesh,
		is_shift,
		1,
		n_vec,
		pos.xyz,
		types,
		pos.nant[0],
		SYMPREAC);
	double rec_vec[3][3];
	RecMat(pos.vec, rec_vec);
	if (ISPIN == 1)
	{
		const char XcrysFile[] = "FERMISURFACE.bxsf";
		FILE* fp3 = fopen(XcrysFile, "w");
		fprintf(fp3, " BEGIN_INFO\n");
		fprintf(fp3, "   Fermi Energy:	%lf\n", efermi);
		fprintf(fp3, " END_INFO\n");
		fprintf(fp3, "\n");
		fprintf(fp3, " BEGIN_BLOCK_BANDGRID_3D\n");
		fprintf(fp3, "   THIS_FILE_WAS_GENERATED_BY_VASPMATE\n");
		fprintf(fp3, "   BEGIN_BANDGRID_3D\n");
		fprintf(fp3, "   %d\n", (int)band_index_up.size());
		fprintf(fp3, "   %d  %d  %d\n", mesh[0], mesh[1], mesh[2]);
		fprintf(fp3, "   0.0 0.0 0.0\n");
		for (int i = 0; i < 3; i++)
			fprintf(fp3, "   %lf  %lf  %f\n", rec_vec[i][0], rec_vec[i][1], rec_vec[i][2]);
		set<int>::iterator s_it, _s_it;
		if (select_band_index.empty())
		{
			s_it = band_index_up.begin();
			_s_it = band_index_up.end();
		}
		else
		{
			s_it = select_band_index.begin();
			_s_it = select_band_index.end();
		}
		for (; s_it != _s_it; s_it++)
		{
			fprintf(fp3, "   band:  %d\n", *s_it);
			vector<int> nonrepet;
			for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
			{
				if (grid_mapping_table[i] == i)
					nonrepet.push_back(i);
			}
			vector<int>::iterator it;
			for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
			{
				it = find(nonrepet.begin(), nonrepet.end(), i);
				if (it != nonrepet.end())
					eigenval[i] = eigenvalue[0][it - nonrepet.begin()][*s_it - 1];
				else
					eigenval[i] = eigenval[grid_mapping_table[i]];
				fprintf(fp3, "    %5lf  ", eigenval[i]);
				if ((i + 1) % 10 == 0 || i == mesh[0] * mesh[1] * mesh[2] - 1)
					fprintf(fp3, "\n");
			}
		}
		fprintf(fp3, "   END_BANDGRID_3D\n");
		fprintf(fp3, " END_BLOCK_BANDGRID_3D\n");
		fclose(fp3);
		printf("Written FEMISURFACE.bxsf file!\n");
	}
	else if (ISPIN == 2)
	{
		vector<const char*> filename = { "FERMISURFACE_UP.bxsf" ,"FERMISURFACE_DW.bxsf" };
		vector<set<int> > band_index = { band_index_up,band_index_dw };
		for (int i = 0; i < 2; i++)
		{
			FILE* _fp = fopen(filename[i], "w");
			fprintf(_fp, " BEGIN_INFO\n");
			fprintf(_fp, "   Fermi Energy:	%lf\n", efermi);
			fprintf(_fp, " END_INFO\n");
			fprintf(_fp, "\n");
			fprintf(_fp, " BEGIN_BLOCK_BANDGRID_3D\n");
			fprintf(_fp, "   THIS_FILE_WAS_GENERATED_BY_VASPMATE\n");
			fprintf(_fp, "   BEGIN_BANDGRID_3D\n");
			fprintf(_fp, "   %d\n", (int)band_index[i].size());
			fprintf(_fp, "   %d  %d  %d\n", mesh[0], mesh[1], mesh[2]);
			fprintf(_fp, "   0.0 0.0 0.0\n");
			for (int i = 0; i < 3; i++)
				fprintf(_fp, "   %lf  %lf  %f\n", rec_vec[i][0], rec_vec[i][1], rec_vec[i][2]);
			set<int>::iterator s_it, _s_it;
			if (select_band_index.empty())
			{
				s_it = band_index[i].begin();
				_s_it = band_index[i].end();
			}
			else
			{
				s_it = select_band_index.begin();
				_s_it = select_band_index.end();
			}
			for (; s_it != _s_it; s_it++)
			{
				fprintf(_fp, "   band:  %d\n", *s_it);
				vector<int> nonrepet;
				for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
				{
					if (grid_mapping_table[i] == i)
						nonrepet.push_back(i);
				}
				vector<int>::iterator it;
				for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
				{
					it = find(nonrepet.begin(), nonrepet.end(), i);
					if (it != nonrepet.end())
						eigenval[i] = eigenvalue[0][it - nonrepet.begin()][*s_it - 1];
					else
						eigenval[i] = eigenval[grid_mapping_table[i]];
					fprintf(_fp, "    %5lf  ", eigenval[i]);
					if ((i + 1) % 10 == 0 || i == mesh[0] * mesh[1] * mesh[2] - 1)
						fprintf(_fp, "\n");
				}
			}
			fprintf(_fp, "   END_BANDGRID_3D\n");
			fprintf(_fp, " END_BLOCK_BANDGRID_3D\n");
			fclose(_fp);
			printf("Written %s file!\n", filename[i]);
		}
	}
	delete[] types;
}

void EIGENVAL::TranEigenToFermiSurface(double efermi, string LORBIT, vector<vector<vector<vector<double> > > > ion_dos
	, vector<vector<vector<vector<double> > >> ion_dos_up, vector<vector<vector<vector<double> > >> ion_dos_dw, int argc, char* argv[])
{
	FILE* fp1 = fopen("KPOINTS", "r");
	char buf[1024];
	if (fp1 == NULL)
	{
		printf("KPOINTS IS NOT EXIST!\n");
		return;
	}
	fgets(buf, 1024, fp1);
	fclose(fp1);
	int mesh[3];
	char ikmesh;
	sscanf(buf, "%c%d%d%d", &ikmesh, &mesh[0], &mesh[1], &mesh[2]);
	FILE* fp2 = fopen("POSCAR", "r");
	if (fp2 == NULL)
	{
		printf("POSCAR IS NOT EXIST!\n");
		return;
	}
	POSCAR pos;
	readposcar(fp2, pos);
	fclose(fp2);
	int max_element = 0;
	for (int i = 0; i < pos.nant[1]; i++)
		max_element = max(max_element, pos.elemnum[i]);
	vector<string> ind;
	vector<string> orb;
	set<int> select_band_index;
	int start_orbit_argv = argc;
	if (argc > 2 && !strcmp("-b", argv[2])) //select band
	{
		for (int i = 3; i < argc; i++)
		{
			if (strcmp(argv[i], "-o"))
				select_band_index.insert(atoi(argv[i]));
			else
			{
				start_orbit_argv = i + 1;
				break;
			}
		}
	}
	if (argc > 2 && !strcmp("-o", argv[2])) //select orbit
		start_orbit_argv = 3;
	for (int i = start_orbit_argv; i < argc; i++)
	{
		if (atoi(LORBIT.c_str()) == 11 && !strcmp(argv[i], "p"))
		{
			orb.push_back("py");
			orb.push_back("pz");
			orb.push_back("px");
		}
		else if (atoi(LORBIT.c_str()) == 11 && !strcmp(argv[i], "d"))
		{
			orb.push_back("dxy");
			orb.push_back("dyz");
			orb.push_back("dz2");
			orb.push_back("dxz");
			orb.push_back("dx2");
		}
		else if (atoi(LORBIT.c_str()) == 11 && !strcmp(argv[i], "f"))
		{
			orb.push_back("f1");
			orb.push_back("f2");
			orb.push_back("f3");
			orb.push_back("f4");
			orb.push_back("f5");
			orb.push_back("f6");
			orb.push_back("f7");
		}
		else if (!IsRightOrbit(argv[i], atoi(LORBIT.c_str()), max_element))
			ind.push_back(argv[i]);
		else
			orb.push_back(argv[i]);
	}
	vector<string> orbit;
	if (atoi(LORBIT.c_str()) == 10 && max_element <= 57)
		orbit = { "s","p","d" };
	if (atoi(LORBIT.c_str()) == 10 && max_element > 57)
		orbit = { "s","p","d","f" };
	if (atoi(LORBIT.c_str()) == 11 && max_element <= 57)
		orbit = { "s","py","pz","px","dxy","dyz","dz2","dxz","dx2" };
	if (atoi(LORBIT.c_str()) == 11 && max_element > 57)
		orbit = { "s","py","pz","px","dxy","dyz","dz2","dxz","dx2","f1","f2","f3","f4","f5","f6","f7" };
	vector<int> num, orbit_ind;
	if (ind.empty() && !orb.empty())
	{
		for (int i = 0; i < orbit.size(); i++)
		{
			if (find(orb.begin(), orb.end(), orbit[i]) != orb.end())
				orbit_ind.push_back(i);
		}
		for (int i = 0; i < pos.nant[0]; i++)
			num.push_back(i + 1);
	}
	else if (orb.empty() && !ind.empty())
	{
		for (int i = 0; i < ind.size(); i++)
		{
			if (IsNum(ind[i]))
			{
				if (atoi(ind[i].c_str()) <= 0 || atoi(ind[i].c_str()) > pos.nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", pos.nant[0]);
					return;
				}
				num.push_back(atoi(ind[i].c_str()));
			}
			else if (strstr(ind[i].c_str(), "-") != NULL)
			{
				int start, end;
				sscanf(ind[i].c_str(), "%d%*c%d", &start, &end);
				if (start > end || start < 1 || end > pos.nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", pos.nant[0]);
					return;
				}
				for (int i = start; i < end + 1; i++)
					num.push_back(i);
			}
			else
			{
				//check symbol is right element
				int IsRightelem = 0;
				for (int j = 0; j < pos.nant[1]; j++)
				{
					if (ind[i] == pos.elemsym[j])
						IsRightelem = 1;
				}
				if (IsRightelem == 0)
				{
					printf("%s are not in element types!\n", ind[i].c_str());
					return;
				}
				//tran element to atom index
				int cnt = 0;
				for (int k = 0; k < pos.nant[1]; k++)
				{
					if (ind[i] == pos.elemsym[k])
					{
						for (int j = 0; j < pos.typenum[k]; j++)
							num.push_back(j + cnt + 1);
					}
					cnt += pos.typenum[k];
				}
			}
		}
		for (int i = 0; i < orbit.size(); i++)
			orbit_ind.push_back(i);
	}
	else if (!orb.empty() && !ind.empty())
	{
		for (int i = 0; i < ind.size(); i++)
		{
			if (IsNum(ind[i]))
			{
				if (atoi(ind[i].c_str()) <= 0 || atoi(ind[i].c_str()) > pos.nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", pos.nant[0]);
					return;
				}
				num.push_back(atoi(ind[i].c_str()));
			}
			else if (strstr(ind[i].c_str(), "-") != NULL)
			{
				int start, end;
				sscanf(ind[i].c_str(), "%d%*c%d", &start, &end);
				if (start > end || start < 1 || end > pos.nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", pos.nant[0]);
					return;
				}
				for (int i = start; i < end + 1; i++)
					num.push_back(i);
			}
			else
			{
				//check symbol is right element
				int IsRightelem = 0;
				for (int j = 0; j < pos.nant[1]; j++)
				{
					if (ind[i] == pos.elemsym[j])
						IsRightelem = 1;
				}
				if (IsRightelem == 0)
				{
					printf("%s are not in element types!\n", ind[i].c_str());
					return;
				}
				//tran element to atom index
				int cnt = 0;
				for (int k = 0; k < pos.nant[1]; k++)
				{
					if (ind[i] == pos.elemsym[k])
					{
						for (int j = 0; j < pos.typenum[k]; j++)
							num.push_back(j + cnt + 1);
					}
					cnt += pos.typenum[k];
				}
			}
		}
		orb = orbit;
		for (int i = 0; i < orbit.size(); i++)
		{
			if (find(orb.begin(), orb.end(), orbit[i]) != orb.end())
				orbit_ind.push_back(i);
		}
	}
	vector<string> select(pos.nant[0], "F");
	for (int i = 0; i < num.size(); i++)
		select[num[i] - 1] = "T";
	vector<string> elemlabel;
	for (int i = 0; i < pos.nant[1]; i++)
	{
		for (int j = 0; j < pos.typenum[i]; j++)
		{
			char elem[20];
			sprintf(elem, "%s%d", pos.elemsym[i], j + 1);
			elemlabel.push_back(elem);
		}
	}
	char select_file[] = "SELECT_ATOMS_LIST";
	FILE* fp3 = fopen(select_file, "w");
	fprintf(fp3, "ATOMS_ID ATOM_LABEL  X_POSITION   Y_POSITION   Z_POSITION   SELECTED?\n");
	for (int i = 0; i < pos.nant[0]; i++)
		fprintf(fp3, "	%d	%s	%lf	%lf	%lf	%s\n", i + 1, elemlabel[i].c_str(), pos.xyz[i][0], pos.xyz[i][1], pos.xyz[i][2], select[i].c_str());
	fclose(fp3);
	printf("Written %s file!\n", select_file);
	char Select_file[] = "SELECT_ORBITS_LIST";
	FILE* fp4 = fopen(Select_file, "w");
	fprintf(fp4, "ORBITALS_ID  ORBITAL_LABEL  SELECTED?\n");
	for (int i = 0; i < orbit.size(); i++)
		fprintf(fp4, "   %2d            %3s            %s\n", i + 1, orbit[i].c_str(), (find(orb.begin(), orb.end(), orbit[i]) == orb.end() ? "F" : "T"));
	fclose(fp4);
	printf("Written %s file!\n", Select_file);
	vector<double> eigenval(mesh[0] * mesh[1] * mesh[2]);
	int grid_address[mesh[0] * mesh[1] * mesh[2]][3];
	int grid_mapping_table[mesh[0] * mesh[1] * mesh[2]];
	int is_shift[3];
	if (ikmesh == 'G')
	{
		is_shift[0] = 0;
		is_shift[1] = 0;
		is_shift[2] = 0;
	}
	else
	{
		is_shift[0] = 1;
		is_shift[1] = 1;
		is_shift[2] = 1;
	}
	int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
	translate_typenum_type(types, pos.nant[1], pos.typenum);
	double n_vec[3][3];
	transpose_matrix(pos.vec, n_vec);
	SpglibDataset* dataset;
	dataset = spg_get_dataset(n_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
	int num_ir = spg_get_ir_reciprocal_mesh(grid_address,
		grid_mapping_table,
		mesh,
		is_shift,
		1,
		n_vec,
		pos.xyz,
		types,
		pos.nant[0],
		SYMPREAC);
	double rec_vec[3][3];
	RecMat(pos.vec, rec_vec);
	if (ISPIN == 1)
	{
		set<int> cal_band;
		if (select_band_index.empty())
			cal_band = band_index_up;
		else
			cal_band = select_band_index;
		vector<vector<double> > dos_weight(cal_band.size(), vector<double>(num_kpts));
		vector<vector<double> > dos_weight_tot(cal_band.size(), vector<double>(mesh[0] * mesh[1] * mesh[2]));
		set<int>::iterator dos_it = cal_band.begin();
		int cnt = 0;
		for (; dos_it != cal_band.end(); dos_it++)
		{
			for (int j = 0; j < num_kpts; j++)
			{
				double temp = 0;
				for (int m = 0; m < num.size(); m++)
				{
					for (int n = 0; n < orbit_ind.size(); n++)
					{
						temp += ion_dos[j][*dos_it - 1][num[m] - 1][orbit_ind[n]];
					}
				}
				dos_weight[cnt][j] = temp;
			}
			cnt++;
		}
		const char FermiSurfFile[] = "FERMISURFACE.frmsf";
		FILE* fp3 = fopen(FermiSurfFile, "w");
		fprintf(fp3, "%d  %d  %d\n", mesh[0], mesh[1], mesh[2]);
		fprintf(fp3, "1\n");
		fprintf(fp3, "%d\n", (int)cal_band.size());
		for (int i = 0; i < 3; i++)
			fprintf(fp3, "   %lf  %lf  %f\n", rec_vec[i][0], rec_vec[i][1], rec_vec[i][2]);
		set<int>::iterator s_it = cal_band.begin(), _s_it = cal_band.end();
		for (int cou = 0; s_it != _s_it; s_it++, cou++)
		{
			vector<int> nonrepet;
			for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
			{
				if (grid_mapping_table[i] == i)
					nonrepet.push_back(i);
			}
			vector<int>::iterator it;
			for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
			{
				it = find(nonrepet.begin(), nonrepet.end(), i);
				if (it != nonrepet.end())
				{
					eigenval[i] = eigenvalue[0][it - nonrepet.begin()][*s_it - 1];
					dos_weight_tot[cou][i] = dos_weight[cou][it - nonrepet.begin()];
				}
				else
				{
					eigenval[i] = eigenval[grid_mapping_table[i]];
					dos_weight_tot[cou][i] = dos_weight_tot[cou][grid_mapping_table[i]];
				}
				fprintf(fp3, "%5lf\n", eigenval[i] - efermi);
			}
		}
		if (!orb.empty() || !ind.empty())
		{
			for (int k = 0; k < dos_weight_tot.size(); k++)
				for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
					fprintf(fp3, "%lf\n", dos_weight_tot[k][i]);
		}
		fclose(fp3);
		printf("Written FERMISURFACE.frmsf file!\n");
	}
	else if (ISPIN == 2)
	{
		set<int> cal_band_up, cal_band_dw;
		if (select_band_index.empty())
		{
			cal_band_up = band_index_up;
			cal_band_dw = band_index_dw;
		}
		else
		{
			cal_band_up = select_band_index;
			cal_band_dw = select_band_index;
		}
		vector<vector<double> > dos_weight_up(cal_band_up.size(), vector<double>(num_kpts));
		vector<vector<double> > dos_weight_dw(cal_band_dw.size(), vector<double>(num_kpts));
		vector<vector<double> > dos_weight_up_tot(cal_band_up.size(), vector<double>(mesh[0] * mesh[1] * mesh[2]));
		vector<vector<double> > dos_weight_dw_tot(cal_band_dw.size(), vector<double>(mesh[0] * mesh[1] * mesh[2]));
		set<int>::iterator dos_it_up = cal_band_up.begin();
		set<int>::iterator dos_it_dw = cal_band_dw.begin();
		int cnt1 = 0, cnt2 = 0;
		for (; dos_it_up != cal_band_up.end(); dos_it_up++)
		{
			for (int j = 0; j < num_kpts; j++)
			{
				double temp_up = 0;
				for (int m = 0; m < num.size(); m++)
					for (int n = 0; n < orbit_ind.size(); n++)
						temp_up += ion_dos_up[j][*dos_it_up - 1][num[m] - 1][orbit_ind[n]];
				dos_weight_up[cnt1][j] = temp_up;
			}
			cnt1++;
		}
		for (; dos_it_dw != cal_band_dw.end(); dos_it_dw++)
		{
			for (int j = 0; j < num_kpts; j++)
			{
				double temp_dw = 0;
				for (int m = 0; m < num.size(); m++)
					for (int n = 0; n < orbit_ind.size(); n++)
						temp_dw += ion_dos_dw[j][*dos_it_dw - 1][num[m] - 1][orbit_ind[n]];
				dos_weight_dw[cnt2][j] = temp_dw;
			}
			cnt2++;
		}
		vector<vector<vector<double> >> dos_weight = { dos_weight_up, dos_weight_dw };
		vector<vector<vector<double> >> dos_weight_tot = { dos_weight_up_tot, dos_weight_dw_tot };
		vector<const char*> filename = { "FERMISURFACE_UP.frmsf" ,"FERMISURFACE_DW.frmsf" };
		vector<set<int> > cal_band_tot = { cal_band_up,cal_band_dw };
		for (int k = 0; k < 2; k++)
		{
			FILE* _fp = fopen(filename[k], "w");
			fprintf(_fp, "%d  %d  %d\n", mesh[0], mesh[1], mesh[2]);
			fprintf(_fp, "1\n");
			fprintf(_fp, "%d\n", (int)cal_band_tot[k].size());
			for (int i = 0; i < 3; i++)
				fprintf(_fp, "   %lf  %lf  %f\n", rec_vec[i][0], rec_vec[i][1], rec_vec[i][2]);
			set<int>::iterator s_it = cal_band_tot[k].begin(), _s_it = cal_band_tot[k].end();
			for (int cou = 0; s_it != _s_it; s_it++, cou++)
			{
				vector<int> nonrepet;
				for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
				{
					if (grid_mapping_table[i] == i)
						nonrepet.push_back(i);
				}
				vector<int>::iterator it;
				for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
				{
					it = find(nonrepet.begin(), nonrepet.end(), i);
					if (it != nonrepet.end())
					{
						eigenval[i] = eigenvalue[0][it - nonrepet.begin()][*s_it - 1];
						dos_weight_tot[k][cou][i] = dos_weight[k][cou][it - nonrepet.begin()];
					}
					else
					{
						eigenval[i] = eigenval[grid_mapping_table[i]];
						dos_weight_tot[k][cou][i] = dos_weight_tot[k][cou][grid_mapping_table[i]];
					}
					fprintf(_fp, "%5lf\n", eigenval[i] - efermi);
				}
			}
			if (!orb.empty() || !ind.empty())
			{
				for (int l = 0; l < dos_weight_tot[k].size(); l++)
					for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
						fprintf(_fp, "%5lf\n", dos_weight_tot[k][l][i]);
			}
			fclose(_fp);
			printf("Written %s file!\n", filename[k]);
		}
	}
	delete[] types;
}

void EIGENVAL::get3dband(vector<int> band_index, double efermi)
{
	// output KX.grd KY.grd
	FILE* fp1 = fopen("POSCAR", "r");
	if (fp1 == NULL)
	{
		printf("POSCAR IS NOT EXIST!\n");
		return;
	}
	POSCAR pos;
	readposcar(fp1, pos);
	fclose(fp1);
	double rec[3][3];
	RecMat(pos.vec, rec);
	FILE* fp2 = fopen("KPOINTS", "r");
	if (fp2 == NULL)
	{
		printf("KPOINTS IS NOT EXIST!\n");
		return;
	}
	char buf[1024];
	int line = 0;
	vector<vector<double> >kpos;
	int cnt = 0;
	int kmesh[3];
	while (fgets(buf, 1024, fp2) != NULL)
	{
		line++;
		if (line == 1 || line == 3)
			continue;
		else if (line == 2)
			kpos.resize(atoi(buf));
		else
		{
			double x, y, z;
			int weight;
			sscanf(buf, "%lf%lf%lf%d", &x, &y, &z, &weight);
			kpos[cnt++] = vector<double>{ x, y, z };
		}
	}
	fclose(fp2);
	vector<vector<double> >kpos_cartes(kpos.size(), vector<double>(3));
	for (int i = 0; i < kpos.size(); i++)
	{
		kpos_cartes[i][0] = kpos[i][0] * rec[0][0] + kpos[i][1] * rec[1][0] + kpos[i][2] * rec[2][0];
		kpos_cartes[i][1] = kpos[i][0] * rec[0][1] + kpos[i][1] * rec[1][1] + kpos[i][2] * rec[2][1];
	}
	FILE* fpkx = fopen("KX.grd", "w");
	FILE* fpky = fopen("KY.grd", "w");
	for (int i = 0; i < kpos.size(); i++)
	{
		fprintf(fpkx, " %lf\n", kpos_cartes[i][0]);
		fprintf(fpky, " %lf\n", kpos_cartes[i][1]);
	}
	fclose(fpkx);
	fclose(fpky);
	printf("Written KX.grd file!\n");
	printf("Written KY.grd file!\n");
	// Output BAND_*.grd
	if (band_index.empty())
	{
		if (ISPIN == 1)
		{
			vector<int> ho_lu = { *band_index_up.begin() ,*band_index_up.begin() + 1 };
			vector<const char*> filename = { "BAND_HOMO.grd","BAND_LUMO.grd" };
			for (int i = 0; i < 2; i++)
			{
				FILE* fp = fopen(filename[i], "w");
				for (int j = 0; j < kpos.size(); j++)
					fprintf(fp, "    %5lf\n", eigenvalue[0][j][ho_lu[i] - 1] - efermi);
				fclose(fp);
				printf("Written %s file!\n", filename[i]);
			}
		}
		else
		{
			vector<vector<int> > ho_lu = { {*band_index_up.begin() ,*band_index_up.begin() + 1},{*band_index_dw.begin() ,*band_index_dw.begin() + 1 } };
			vector<vector<const char*> > filename = { { "BAND_HOMO_UP.grd","BAND_LUMO_UP.grd" } ,{ "BAND_HOMO_DW.grd","BAND_LUMO_DW.grd" } };
			for (int k = 0; k < ISPIN; k++)
			{
				for (int i = 0; i < 2; i++)
				{
					FILE* fp = fopen(filename[k][i], "w");
					for (int j = 0; j < kpos.size(); j++)
						fprintf(fp, "    %5lf\n", eigenvalue[k][j][ho_lu[k][i] - 1] - efermi);
					fclose(fp);
					printf("Written %s file!\n", filename[k][i]);
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < band_index.size(); i++)
		{
			for (int k = 0; k < ISPIN; k++)
			{
				vector<string> tot_file(ISPIN, "BAND_B");
				if (ISPIN == 1)
					tot_file[0] += to_string(band_index[i]) + ".grd";
				else
				{
					tot_file[0] += to_string(band_index[i]) + "_UP.grd";
					tot_file[1] += to_string(band_index[i]) + "_DW.grd";
				}
				FILE* fp = fopen(tot_file[k].c_str(), "w");
				for (int j = 0; j < kpos.size(); j++)
					fprintf(fp, "    %5lf\n", eigenvalue[k][j][band_index[i] - 1] - efermi);
				fclose(fp);
				printf("Written %s file!\n", tot_file[k].c_str());
			}
		}
	}
}

void EIGENVAL::read_kpts_mapping_table(vector<vector<double> >& kmesh_primcell, vector<vector<double> >& kmesh_supercell, vector<vector<int> >& gvector_supercell)
{
	FILE* fp = fopen("KPOINTS_MAPPING_TABLE", "r");
	if (fp == NULL)
	{
		printf("KPOINTS_MAPPING_TABLE IS NOT EXIST!\n");
		return;
	}
	char buf[1024];
	fgets(buf, 1024, fp);
	while (fgets(buf, 1024, fp))
	{
		vector<double> a(3), b(3);
		vector<int> c(3);
		sscanf(buf, "%lf%lf%lf%lf%lf%lf%d%d%d", &a[0], &a[1], &a[2], &b[0], &b[1], &b[2], &c[0], &c[1], &c[2]);
		kmesh_primcell.push_back(a);
		kmesh_supercell.push_back(b);
		gvector_supercell.push_back(c);
	}
	fclose(fp);
}

void EIGENVAL::getunfold()
{
	ionizing::WAVEHIGH wave{ "WAVECAR" };
	int npws_max = wave.getHeader().npws_max;
	vector<int> index_G_primcell(npws_max);
	int nbands = wave.getHeader()._nBands;
	int nkpts = wave.getHeader()._nKpoints;
	int nspin = wave.getInfo()._nSpin;
	vector<vector<vector<double> > > E_km(nbands, vector<vector<double> >(nkpts, vector<double>(nspin)));
	vector<vector<vector<double> > > P_km(nbands, vector<vector<double> >(nkpts, vector<double>(nspin)));
	Eigen::VectorXi npws = wave.getNPlaneWaves();
	ionizing::MatX3d kmesh = wave.getKVectors();
	vector<Eigen::MatrixX3i> gmesh(nkpts);
	for (int i = 0; i < nkpts; i++)
		gmesh[i] = wave.getGVectors(i, false);
	double tran[3][3];
	get_tran_matrix("redef.in", tran);
	FILE* fp = fopen("POSCAR", "r");
	if (fp == NULL)
	{
		printf("POSCAR IS NOT EXIST!\n");
		return;
	}
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	double inv_tran[3][3];
	brinv(tran, inv_tran);
	double prim_vec[3][3];
	multi_mat(inv_tran, pos.vec, prim_vec);
	double rec_vec[3][3];
	RecMat(prim_vec, rec_vec);
	string LORBIT = GetInfoINCAR("LORBIT");
	if (LORBIT.length() == 0)
		LORBIT = "0";
	DOS dos;
	dos.readDoscar(nspin, atoi(LORBIT.c_str()));
	double efermi = dos.fermi_energy();
	EIGENVAL eigen("EIGENVAL", efermi);
	vector<vector<vector<double> > > eigenvalue = eigen.get_eigenvalue(); //ispin*num_kpts*num_bands
	BAND band;
	vector<double> kline = band.hse_getkpointsline(nkpts,this->position, prim_vec);
	vector<vector<double> > kmesh_primcell, kmesh_supercell;
	vector<vector<int> > gvector_supercell;
	read_kpts_mapping_table(kmesh_primcell, kmesh_supercell, gvector_supercell);
	int ibz_kpt = band.getibz_kpt();
	int ng = 0;
	vector<string> file(nspin);
	if (nspin == 1)
	{
		file[0] = "EBS.dat";
	}
	else if (nspin == 2)
	{
		file[0] = "EBS_UP.dat";
		file[1] = "EBS_DW.dat";
	}
	for (int is = 0; is < nspin; is++)
	{
		FILE* fp1 = fopen(file[is].c_str(), "w");
		for (int ik = 0; ik < nkpts - ibz_kpt; ik++)
		{
			get_overlap_G_single_kpt(npws(ik + ibz_kpt), gmesh[ik + ibz_kpt], index_G_primcell, ng, inv_tran);
			fprintf(fp1, "#Kpoint-index: %d\n", ik + 1);
			for (int ig = 0; ig < ng - 1; ig++)
			{
				vector<double> g_plus(3);
				for (int j = 0; j < 3; j++)
					g_plus[j] = gmesh[ik + ibz_kpt](index_G_primcell[ig], j) + gvector_supercell[ik][j];
				for (int ipw = 0; ipw < npws(ik + ibz_kpt); ipw++)
				{
					if (g_plus[0] == gmesh[ik + ibz_kpt](ipw, 0) && g_plus[1] == gmesh[ik + ibz_kpt](ipw, 1) && g_plus[2] == gmesh[ik + ibz_kpt](ipw, 2))
					{
						for (int ib = 0; ib < nbands; ib++)
						{
							const ionizing::Veccd& coeff_vec = wave.getBandCoeff(is, ik + ibz_kpt, ib, true);
							//P_km[ib][ik][is] += fabs(coeff_vec(ipw) * std::conj(coeff_vec(ipw)));
						}
					}
				}
			}
			for (int ib = 0; ib < nbands; ib++)
			{
				E_km[ib][ik][is] = eigenvalue[is][ik + ibz_kpt][ib] - efermi;
				if (P_km[ib][ik][is] > 0.05)
					fprintf(fp1, "%lf	%lf	%lf\n", kline[ik], E_km[ib][ik][is], P_km[ib][ik][is]);
			}
		}
		fclose(fp1);
		printf("Written %s file!\n", file[is].c_str());
	}
}

void get_overlap_G_single_kpt(int npws_single_kpt, Eigen::MatrixX3i G_single_kpt, vector<int>& index_G_primcell, int& ng, double tran[3][3])
{
	ng = 0;
	double eps = 1e-5;
	for (int ipw = 0; ipw < npws_single_kpt; ipw++)
	{
		double a = tran[0][0] * G_single_kpt(ipw, 0) + tran[0][1] * G_single_kpt(ipw, 1) + tran[0][2] * G_single_kpt(ipw, 2);
		double b = tran[1][0] * G_single_kpt(ipw, 0) + tran[1][1] * G_single_kpt(ipw, 1) + tran[1][2] * G_single_kpt(ipw, 2);
		double c = tran[2][0] * G_single_kpt(ipw, 0) + tran[2][1] * G_single_kpt(ipw, 1) + tran[2][2] * G_single_kpt(ipw, 2);
		double delta = sqrt(fabs(a) - (int)a) * (fabs(a) - (int)a) + (fabs(b) - (int)b) * (fabs(b) - (int)b) + (fabs(c) - (int)c) * (fabs(c) - (int)c);
		if (delta < eps)
			index_G_primcell[ng++] = ipw;
	}
}
