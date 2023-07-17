#include"dos.h"
#include"read_write.h"
#include"tools.h"
#include<spa_plot.h>
//
//static int IsOrbit(char symbol[])
//{
//	if (!strcmp(symbol, "s") || !strcmp(symbol, "p") || !strcmp(symbol, "d") || !strcmp(symbol, "f") || !strcmp(symbol, "px")
//		|| !strcmp(symbol, "py") || !strcmp(symbol, "pz") || !strcmp(symbol, "dxy") || !strcmp(symbol, "dyz") || !strcmp(symbol, "dxz")
//		|| !strcmp(symbol, "dz2") || !strcmp(symbol, "dx2") || !strcmp(symbol, "f1") || !strcmp(symbol, "f2") || !strcmp(symbol, "f3")
//		|| !strcmp(symbol, "f4") || !strcmp(symbol, "f5") || !strcmp(symbol, "f6") || !strcmp(symbol, "f7"))
//		return 1;
//	else
//		return 0;
//}

int IsRightOrbit(string orbit, int LORBIT, int max_element)
{
	if (LORBIT == 10)
	{
		if ((orbit == "s" || orbit == "p" || orbit == "d") && (max_element <= 57))
			return 1;
		if ((orbit == "s" || orbit == "p" || orbit == "d" || orbit == "f") && (max_element > 57))
			return 1;
	}
	if (LORBIT == 11)
	{
		if ((orbit == "s" || orbit == "px" || orbit == "py" || orbit == "pz" || orbit == "dxy"
			|| orbit == "dyz" || orbit == "dxz" || orbit == "dz2" || orbit == "dx2" || orbit == "s"
			|| orbit == "p" || orbit == "d") && (max_element <= 57))
			return 1;
		if ((orbit == "s" || orbit == "px" || orbit == "py" || orbit == "pz" || orbit == "dxy"
			|| orbit == "dyz" || orbit == "dxz" || orbit == "dz2" || orbit == "dx2" || orbit == "f1"
			|| orbit == "f2" || orbit == "f3" || orbit == "f4" || orbit == "f5" || orbit == "f6"
			|| orbit == "f7" || orbit == "s" || orbit == "p" || orbit == "d" || orbit == "f") && (max_element > 57))
			return 1;
	}
	return 0;
}

vector<int> TranArgvToAtomIndex(int argc, char* argv[], int nant[], int typenum[], char elemsym[][3], double xyz[][3], char label[])
{
	vector<int> num;
	vector<string> select(nant[0], "F");
	vector<string> elemlabel;
	for (int i = 0; i < nant[1]; i++)
	{
		for (int j = 0; j < typenum[i]; j++)
		{
			char elem[20];
			sprintf(elem, "%s%d", elemsym[i], j + 1);
			elemlabel.push_back(elem);
		}
	}
	for (int i = 3; i < argc; i++)
	{
		if (!strcmp(argv[i], "all"))
		{
			for (int i = 0; i < nant[0]; i++)
			{
				num.push_back(i + 1);
				select[i] = "T";
			}
			break;
		}
		else if (strstr(argv[i], "-") != NULL) // 
		{
			int start, end;
			sscanf(argv[i], "%d%*c%d", &start, &end);
			if (start > end || start<1 || end>nant[0])
			{
				printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
				return num;
			}
			for (int i = start; i < end + 1; i++)
			{
				num.push_back(i);
				select[i - 1] = "T";
			}
		}
		else if (IsNum(argv[i]))
		{
			if (atoi(argv[i]) <= 0 || atoi(argv[i]) > nant[0])
			{
				printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
				return num;
			}
			num.push_back(atoi(argv[i]));
			select[atoi(argv[i]) - 1] = "T";
		}
		else if (!IsNum(argv[i]))
		{
			int IsRightelem = 0;
			for (int j = 0; j < nant[1]; j++)
			{
				if (!strcmp(argv[i], elemsym[j]))
					IsRightelem = 1;
			}
			if (IsRightelem == 0)
			{
				printf("%s is not in element types.\n", argv[i]);
				continue;
			}
			else
			{
				int cnt = 0;
				for (int k = 0; k < nant[1]; k++)
				{
					if (!strcmp(argv[i], elemsym[k]))
					{
						for (int j = 0; j < typenum[k]; j++)
						{
							num.push_back(j + cnt + 1);
							select[j + cnt] = "T";
						}
					}
					cnt += typenum[k];
				}
			}
		}
	}
	//erase same value
	sort(num.begin(), num.end());
	num.erase(unique(num.begin(), num.end()), num.end());
	if (label == NULL)
		return num;
	char select_file[] = "SELECT_ATOMS_LIST";
	FILE* fp = fopen(select_file, "w");
	fprintf(fp, "ATOMS_ID ATOM_LABEL  X_POSITION   Y_POSITION   Z_POSITION   SELECTED?\n");
	for (int i = 0; i < nant[0]; i++)
		fprintf(fp, "	%d	%s	%lf	%lf	%lf	%s\n", i + 1, elemlabel[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], select[i].c_str());
	fclose(fp);
	printf("Written %s file!\n", select_file);
	return num;
}

void DOS::readDoscar(int ISPIN, int LORBIT)
{
	/*if (ISPIN != 1 && ISPIN != 2)
	{
		printf("ISPIN IS ERROR IN INCAR!\n");
		return;
	}
	if (LORBIT != 10 && LORBIT != 11)
	{
		printf("LORBIT IS ERROR IN INCAR!\n");
		return;
	}*/
	FILE* fp = fopen("DOSCAR", "r");
	if (fp == NULL)
	{
		perror("DOSCAR IS NOT EXIST!\n");
		exit(1);
	}
	char buf[1024];
	int line = 0;
	int iat = 0;
	int resize_flag = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		line++;
		if (line == 1)
		{
			sscanf(buf, "%d%d%d%d", &IonWithSphere, &Ion, &Pdos, &NCDIJ);
			continue;
		}
		if (line == 2)
		{
			sscanf(buf, "%lf%lf%lf%lf%lf", &volume, &vec[0], &vec[1], &vec[2], &POTIM);
			continue;
		}
		if (line == 3)
		{
			sscanf(buf, "%lf", &init_T);
			continue;
		}
		if (line == 4)
		{
			car = buf;
			continue;
		}
		if (line == 5)
		{
			system = buf;
			continue;
		}
		if (line == 6)
		{
			sscanf(buf, "%lf%lf%d%lf%lf", &Emax, &Emin, &NEDOS, &efermi, &dosnum);
			//read fermi energy from self-consistent calculation
			FILE* fp1 = fopen("FERMI_LEVEL", "r");
			if (fp1 == NULL)
			{
				printf("You should use the self-consistent calculation fermi_level!\n");
				continue;
			}
			else
			{
				fscanf(fp1, "%lf", &efermi);
				fclose(fp1);
			}
			continue;
		}
		if (line >= 7 && line < 7 + NEDOS)
		{
			if (ISPIN == 1)
			{
				double temp_energy = 0;
				double temp_dos = 0;
				double temp_integ_dos = 0;
				sscanf(buf, "%lf%lf%lf", &temp_energy, &temp_dos, &temp_integ_dos);
				energy.push_back(temp_energy);
				dos_nospin.push_back(temp_dos);
				integ_dos_nospin.push_back(temp_integ_dos);
			}
			else if (ISPIN == 2)
			{
				double temp_energy = 0;
				double temp_dos_up = 0;
				double temp_dos_dw = 0;
				double temp_integ_dos_up = 0;
				double temp_integ_dos_dw = 0;
				sscanf(buf, "%lf%lf%lf%lf%lf", &temp_energy, &temp_dos_up, &temp_dos_dw,
					&temp_integ_dos_up, &temp_integ_dos_dw);
				energy.push_back(temp_energy);
				dos_up.push_back(temp_dos_up);
				dos_dw.push_back(temp_dos_dw);
				integ_dos_up.push_back(temp_integ_dos_up);
				integ_dos_dw.push_back(temp_integ_dos_dw);
			}
			continue;
		}
		if (line == 6 + (NEDOS + 1) * (iat + 1))
		{
			iat++;
			continue; //the info is repeated with line 6
		}
		if (LORBIT == 10 && ISPIN == 1)
		{
			double temp_energy_ion = 0;
			double temp_s_dos = 0;
			double temp_p_dos = 0;
			double temp_d_dos = 0;
			double temp_f_dos = 0;
			sscanf(buf, "%lf%lf%lf%lf%lf", &temp_energy_ion, &temp_s_dos, &temp_p_dos, &temp_d_dos, &temp_f_dos);
			if (!resize_flag)
			{
				energy_ion.resize(Ion);
				s_dos.resize(Ion);
				p_dos.resize(Ion);
				d_dos.resize(Ion);
				f_dos.resize(Ion);
				resize_flag = 1;
			}
			energy_ion[iat - 1].push_back(temp_energy_ion);
			s_dos[iat - 1].push_back(temp_s_dos);
			p_dos[iat - 1].push_back(temp_p_dos);
			d_dos[iat - 1].push_back(temp_d_dos);
			f_dos[iat - 1].push_back(temp_f_dos);
		}
		if (LORBIT == 11 && ISPIN == 1)
		{
			double temp_energy_ion = 0;
			double temp_s0_dos = 0;
			double temp_py_dos = 0;
			double temp_pz_dos = 0;
			double temp_px_dos = 0;
			double temp_dxy_dos = 0;
			double temp_dyz_dos = 0;
			double temp_dz2_dos = 0;
			double temp_dxz_dos = 0;
			double temp_dx2_dos = 0;
			double temp_f1_dos = 0;
			double temp_f2_dos = 0;
			double temp_f3_dos = 0;
			double temp_f4_dos = 0;
			double temp_f5_dos = 0;
			double temp_f6_dos = 0;
			double temp_f7_dos = 0;
			sscanf(buf, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp_energy_ion,
				&temp_s0_dos,
				&temp_py_dos,
				&temp_pz_dos,
				&temp_px_dos,
				&temp_dxy_dos,
				&temp_dyz_dos,
				&temp_dz2_dos,
				&temp_dxz_dos,
				&temp_dx2_dos,
				&temp_f1_dos,
				&temp_f2_dos,
				&temp_f3_dos,
				&temp_f4_dos,
				&temp_f5_dos,
				&temp_f6_dos,
				&temp_f7_dos
			);
			if (!resize_flag)
			{
				energy_ion.resize(Ion);
				s0_dos.resize(Ion);
				py_dos.resize(Ion);
				pz_dos.resize(Ion);
				px_dos.resize(Ion);
				dxy_dos.resize(Ion);
				dyz_dos.resize(Ion);
				dz2_dos.resize(Ion);
				dxz_dos.resize(Ion);
				dx2_dos.resize(Ion);
				f1_dos.resize(Ion);
				f2_dos.resize(Ion);
				f3_dos.resize(Ion);
				f4_dos.resize(Ion);
				f5_dos.resize(Ion);
				f6_dos.resize(Ion);
				f7_dos.resize(Ion);
				resize_flag = 1;
			}
			energy_ion[iat - 1].push_back(temp_energy_ion);
			s0_dos[iat - 1].push_back(temp_s0_dos);
			py_dos[iat - 1].push_back(temp_py_dos);
			pz_dos[iat - 1].push_back(temp_pz_dos);
			px_dos[iat - 1].push_back(temp_px_dos);
			dxy_dos[iat - 1].push_back(temp_dxy_dos);
			dyz_dos[iat - 1].push_back(temp_dyz_dos);
			dz2_dos[iat - 1].push_back(temp_dz2_dos);
			dxz_dos[iat - 1].push_back(temp_dxz_dos);
			dx2_dos[iat - 1].push_back(temp_dx2_dos);
			f1_dos[iat - 1].push_back(temp_f1_dos);
			f2_dos[iat - 1].push_back(temp_f2_dos);
			f3_dos[iat - 1].push_back(temp_f3_dos);
			f4_dos[iat - 1].push_back(temp_f4_dos);
			f5_dos[iat - 1].push_back(temp_f5_dos);
			f6_dos[iat - 1].push_back(temp_f6_dos);
			f7_dos[iat - 1].push_back(temp_f7_dos);
		}
		if (LORBIT == 10 && ISPIN == 2)
		{
			double temp_energy_ion = 0;
			double temp_s_dos_up = 0;
			double temp_p_dos_up = 0;
			double temp_d_dos_up = 0;
			double temp_s_dos_dw = 0;
			double temp_p_dos_dw = 0;
			double temp_d_dos_dw = 0;
			double temp_f_dos_up = 0;
			double temp_f_dos_dw = 0;
			sscanf(buf, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &temp_energy_ion, &temp_s_dos_up, &temp_s_dos_dw, &temp_p_dos_up,
				&temp_p_dos_dw, &temp_d_dos_up, &temp_d_dos_dw, &temp_f_dos_up, &temp_f_dos_dw);
			energy_ion.resize(Ion);
			s_dos_up.resize(Ion);
			p_dos_up.resize(Ion);
			d_dos_up.resize(Ion);
			f_dos_up.resize(Ion);
			s_dos_dw.resize(Ion);
			p_dos_dw.resize(Ion);
			d_dos_dw.resize(Ion);
			f_dos_dw.resize(Ion);
			energy_ion[iat - 1].push_back(temp_energy_ion);
			s_dos_up[iat - 1].push_back(temp_s_dos_up);
			p_dos_up[iat - 1].push_back(temp_p_dos_up);
			d_dos_up[iat - 1].push_back(temp_d_dos_up);
			f_dos_up[iat - 1].push_back(temp_d_dos_up);
			s_dos_dw[iat - 1].push_back(temp_s_dos_dw);
			p_dos_dw[iat - 1].push_back(temp_p_dos_dw);
			d_dos_dw[iat - 1].push_back(temp_d_dos_dw);
			f_dos_dw[iat - 1].push_back(temp_d_dos_up);
		}
		if (LORBIT == 11 && ISPIN == 2)
		{
			double temp_energy_ion = 0;
			double temp_s0_dos_up = 0;
			double temp_py_dos_up = 0;
			double temp_pz_dos_up = 0;
			double temp_px_dos_up = 0;
			double temp_dxy_dos_up = 0;
			double temp_dyz_dos_up = 0;
			double temp_dz2_dos_up = 0;
			double temp_dxz_dos_up = 0;
			double temp_dx2_dos_up = 0;
			double temp_f1_dos_up = 0;
			double temp_f2_dos_up = 0;
			double temp_f3_dos_up = 0;
			double temp_f4_dos_up = 0;
			double temp_f5_dos_up = 0;
			double temp_f6_dos_up = 0;
			double temp_f7_dos_up = 0;
			double temp_s0_dos_dw = 0;
			double temp_py_dos_dw = 0;
			double temp_pz_dos_dw = 0;
			double temp_px_dos_dw = 0;
			double temp_dxy_dos_dw = 0;
			double temp_dyz_dos_dw = 0;
			double temp_dz2_dos_dw = 0;
			double temp_dxz_dos_dw = 0;
			double temp_dx2_dos_dw = 0;
			double temp_f1_dos_dw = 0;
			double temp_f2_dos_dw = 0;
			double temp_f3_dos_dw = 0;
			double temp_f4_dos_dw = 0;
			double temp_f5_dos_dw = 0;
			double temp_f6_dos_dw = 0;
			double temp_f7_dos_dw = 0;
			sscanf(buf, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
				&temp_energy_ion,
				&temp_s0_dos_up,
				&temp_s0_dos_dw,
				&temp_py_dos_up,
				&temp_py_dos_dw,
				&temp_pz_dos_up,
				&temp_pz_dos_dw,
				&temp_px_dos_up,
				&temp_px_dos_dw,
				&temp_dxy_dos_up,
				&temp_dxy_dos_dw,
				&temp_dyz_dos_up,
				&temp_dyz_dos_dw,
				&temp_dz2_dos_up,
				&temp_dz2_dos_dw,
				&temp_dxz_dos_up,
				&temp_dxz_dos_dw,
				&temp_dx2_dos_up,
				&temp_dx2_dos_dw,
				&temp_f1_dos_up,
				&temp_f1_dos_dw,
				&temp_f2_dos_up,
				&temp_f2_dos_dw,
				&temp_f3_dos_up,
				&temp_f3_dos_dw,
				&temp_f4_dos_up,
				&temp_f4_dos_dw,
				&temp_f5_dos_up,
				&temp_f5_dos_dw,
				&temp_f6_dos_up,
				&temp_f6_dos_dw,
				&temp_f7_dos_up,
				&temp_f7_dos_dw);
			energy_ion.resize(Ion);
			s0_dos_up.resize(Ion);
			py_dos_up.resize(Ion);
			pz_dos_up.resize(Ion);
			px_dos_up.resize(Ion);
			dxy_dos_up.resize(Ion);
			dyz_dos_up.resize(Ion);
			dz2_dos_up.resize(Ion);
			dxz_dos_up.resize(Ion);
			dx2_dos_up.resize(Ion);
			s0_dos_dw.resize(Ion);
			py_dos_dw.resize(Ion);
			pz_dos_dw.resize(Ion);
			px_dos_dw.resize(Ion);
			dxy_dos_dw.resize(Ion);
			dyz_dos_dw.resize(Ion);
			dz2_dos_dw.resize(Ion);
			dxz_dos_dw.resize(Ion);
			dx2_dos_dw.resize(Ion);
			f1_dos_up.resize(Ion);
			f1_dos_dw.resize(Ion);
			f2_dos_up.resize(Ion);
			f2_dos_dw.resize(Ion);
			f3_dos_up.resize(Ion);
			f3_dos_dw.resize(Ion);
			f4_dos_up.resize(Ion);
			f4_dos_dw.resize(Ion);
			f5_dos_up.resize(Ion);
			f5_dos_dw.resize(Ion);
			f6_dos_up.resize(Ion);
			f6_dos_dw.resize(Ion);
			f7_dos_up.resize(Ion);
			f7_dos_dw.resize(Ion);
			energy_ion[iat - 1].push_back(temp_energy_ion);
			s0_dos_up[iat - 1].push_back(temp_s0_dos_up);
			py_dos_up[iat - 1].push_back(temp_py_dos_up);
			pz_dos_up[iat - 1].push_back(temp_pz_dos_up);
			px_dos_up[iat - 1].push_back(temp_px_dos_up);
			dxy_dos_up[iat - 1].push_back(temp_dxy_dos_up);
			dyz_dos_up[iat - 1].push_back(temp_dyz_dos_up);
			dz2_dos_up[iat - 1].push_back(temp_dz2_dos_up);
			dxz_dos_up[iat - 1].push_back(temp_dxz_dos_up);
			dx2_dos_up[iat - 1].push_back(temp_dx2_dos_up);
			s0_dos_dw[iat - 1].push_back(temp_s0_dos_dw);
			py_dos_dw[iat - 1].push_back(temp_py_dos_dw);
			pz_dos_dw[iat - 1].push_back(temp_pz_dos_dw);
			px_dos_dw[iat - 1].push_back(temp_px_dos_dw);
			dxy_dos_dw[iat - 1].push_back(temp_dxy_dos_dw);
			dyz_dos_dw[iat - 1].push_back(temp_dyz_dos_dw);
			dz2_dos_dw[iat - 1].push_back(temp_dz2_dos_dw);
			dxz_dos_dw[iat - 1].push_back(temp_dxz_dos_dw);
			dx2_dos_dw[iat - 1].push_back(temp_dx2_dos_dw);
			f1_dos_up[iat - 1].push_back(temp_f1_dos_up);
			f1_dos_dw[iat - 1].push_back(temp_f1_dos_dw);
			f2_dos_up[iat - 1].push_back(temp_f2_dos_up);
			f2_dos_dw[iat - 1].push_back(temp_f2_dos_dw);
			f3_dos_up[iat - 1].push_back(temp_f3_dos_up);
			f3_dos_dw[iat - 1].push_back(temp_f3_dos_dw);
			f4_dos_up[iat - 1].push_back(temp_f4_dos_up);
			f4_dos_dw[iat - 1].push_back(temp_f4_dos_dw);
			f5_dos_up[iat - 1].push_back(temp_f5_dos_up);
			f5_dos_dw[iat - 1].push_back(temp_f5_dos_dw);
			f6_dos_up[iat - 1].push_back(temp_f6_dos_up);
			f6_dos_dw[iat - 1].push_back(temp_f6_dos_dw);
			f7_dos_up[iat - 1].push_back(temp_f7_dos_up);
			f7_dos_dw[iat - 1].push_back(temp_f7_dos_dw);
		}
	}
	orbit_dos = {
		//0-3
		s_dos, p_dos, d_dos, f_dos,
		//4-19
		s0_dos, py_dos, pz_dos, px_dos, dxy_dos, dyz_dos, dz2_dos, dxz_dos, dx2_dos,
		f1_dos, f2_dos, f3_dos, f4_dos, f5_dos, f6_dos, f7_dos,
		//20-27
		s_dos_up, p_dos_up, d_dos_up, f_dos_up, s_dos_dw,p_dos_dw, d_dos_dw, f_dos_dw,
		//28-59
		s0_dos_up, py_dos_up, pz_dos_up, px_dos_up,dxy_dos_up, dyz_dos_up, dz2_dos_up, dxz_dos_up, dx2_dos_up,
		f1_dos_up, f2_dos_up, f3_dos_up, f4_dos_up,f5_dos_up, f6_dos_up, f7_dos_up,
		s0_dos_dw, py_dos_dw, pz_dos_dw, px_dos_dw,dxy_dos_dw, dyz_dos_dw, dz2_dos_dw, dxz_dos_dw, dx2_dos_dw,
		f1_dos_dw, f2_dos_dw, f3_dos_dw, f4_dos_dw,f5_dos_dw, f6_dos_dw, f7_dos_dw
	};
	Idos.resize(orbit_dos.size());
	for (int i = 0; i < orbit_dos.size(); i++)
	{
		Idos[i].resize(Ion);
		for (int j = 0; j < Ion; j++)
		{
			if (orbit_dos[i].size() != 0)
				Idos[i][j] = Intergral(energy, orbit_dos[i][j]);
		}
	}
	fclose(fp);
	return;
}

double DOS::fermi_energy()
{
	FILE* fp1 = fopen("Fermi_Energy", "w");
	fprintf(fp1, "		#The Fermi_Energy is read from DOSCAR.The data in TODS.dat has shifted efermi to 0 ev.\n");
	fprintf(fp1, "		%lf", this->efermi);
	fclose(fp1);
	printf("Written Fermi_Energy File!\n");
	return this->efermi;
}

void DOS::TDos()
{
	FILE* fp = fopen("TDOS.dat", "w");
	FILE* _fp = fopen("ITDOS.dat", "w");
	string ISPIN = GetInfoINCAR("ISPIN");
	if (ISPIN.length() == 0 || atoi(ISPIN.c_str()) == 1)
	{
		vector<double> IDOS = Intergral(this->energy, this->dos_nospin);
		fprintf(fp, "  #Energy        TDOS\n");
		fprintf(_fp, "  #Energy        TDOS\n");
		for (int i = 0; i < this->NEDOS; i++)
		{
			fprintf(fp, "%lf        %lf\n", this->energy[i] - this->efermi, this->dos_nospin[i]);
			fprintf(_fp, "%lf        %lf\n", this->energy[i] - this->efermi, IDOS[i]);
		}
	}
	else if (atoi(ISPIN.c_str()) == 2)
	{
		vector<double> IDOS_up = Intergral(this->energy, this->dos_up);
		vector<double> IDOS_dw = Intergral(this->energy, this->dos_dw);
		fprintf(fp, "  #Energy        TDOS\n");
		fprintf(_fp, "  #Energy        TDOS\n");
		fprintf(fp, "  #Energy        TDOS-UP        TDOS-DOWN\n");
		fprintf(_fp, "  #Energy        TDOS-UP        TDOS-DOWN\n");
		for (int i = 0; i < this->NEDOS; i++)
		{
			fprintf(fp, "%lf        %lf        %lf\n", this->energy[i] - this->efermi, this->dos_up[i], this->dos_dw[i]);
			fprintf(_fp, "%lf        %lf        %lf\n", this->energy[i] - this->efermi, IDOS_up[i], IDOS_dw[i]);
		}
	}
	fclose(fp);
	fclose(_fp);
	string spa = convert_TDOS("TDOS.dat", atoi(ISPIN.c_str()));
	string Ispa = convert_TDOS("ITDOS.dat", atoi(ISPIN.c_str()));
	printf("Written TDOS.dat File!\n");
	printf("Written ITDOS.dat File!\n");
	printf("Written %s File!\n", spa.c_str());
	printf("Written %s File!\n", Ispa.c_str());
	fermi_energy();
}

void DOS::ADos(int num, const char elemlabel[], int headind, int ISPIN, int LORBIT, int max_element)
{
	if (num <= 0 || num > this->Ion)
	{
		printf("THE ATOM_INDEX IS ERROR! Please input atom_index in 1 <= & <= %d\n", this->Ion);
		return;
	}
	if (ISPIN == 1 && LORBIT == 10)
	{
		char file[20];
		char _file[20];
		if (elemlabel == NULL)
		{
			sprintf(file, "PDOS_A%d.dat", num);
			sprintf(_file, "IPDOS_A%d.dat", num);
		}
		else
		{
			sprintf(file, "PDOS_%s%d.dat", elemlabel, num - headind);
			sprintf(_file, "IPDOS_%s%d.dat", elemlabel, num - headind);
		}
		FILE* fp = fopen(file, "w");
		FILE* _fp = fopen(_file, "w");
		if (max_element <= 57)
		{
			fprintf(fp, "     #Energy           s          p          d          tot\n");
			fprintf(_fp, "     #Energy           s          p          d          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos = 0;
				double tot_idos = 0;
				fprintf(fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 0; k < 3; k++)
				{
					fprintf(fp, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos += this->orbit_dos[k][num - 1][i];
					tot_idos += this->Idos[k][num - 1][i];
				}
				fprintf(fp, "%lf     ", tot_dos);
				fprintf(_fp, "%lf     ", tot_idos);
				fprintf(fp, "\n");
				fprintf(_fp, "\n");
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp, "     #Energy           s          p          d          f          tot\n");
			fprintf(_fp, "     #Energy           s          p          d          f          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos = 0;
				double tot_idos = 0;
				fprintf(fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 0; k < 4; k++)
				{
					fprintf(fp, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos += this->orbit_dos[k][num - 1][i];
					tot_idos += this->Idos[k][num - 1][i];
				}
				fprintf(fp, "%lf     ", tot_dos);
				fprintf(_fp, "%lf     ", tot_idos);
				fprintf(fp, "\n");
				fprintf(_fp, "\n");
			}
		}
		fclose(fp);
		fclose(_fp);
		printf("Written %s FILE!\n", file);
		printf("Written %s FILE!\n", _file);
		string spa = convert_ADOS(file);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_file);
		printf("Written %s File!\n", _spa.c_str());
	}
	if (ISPIN == 1 && LORBIT == 11)
	{
		char file[20];
		char _file[20];
		if (elemlabel == NULL)
		{
			sprintf(file, "PDOS_A%d.dat", num);
			sprintf(_file, "IPDOS_A%d.dat", num);
		}
		else
		{
			sprintf(file, "PDOS_%s%d.dat", elemlabel, num - headind);
			sprintf(_file, "IPDOS_%s%d.dat", elemlabel, num - headind);
		}
		FILE* fp = fopen(file, "w");
		FILE* _fp = fopen(_file, "w");
		if (max_element <= 57)
		{
			fprintf(fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot\n");
			fprintf(_fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos = 0;
				double tot_idos = 0;
				fprintf(fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 4; k < 13; k++)
				{
					fprintf(fp, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos += this->orbit_dos[k][num - 1][i];
					tot_idos += this->Idos[k][num - 1][i];
				}
				fprintf(fp, "%lf     ", tot_dos);
				fprintf(_fp, "%lf     ", tot_idos);
				fprintf(fp, "\n");
				fprintf(_fp, "\n");
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total\n");
			fprintf(_fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos = 0;
				double tot_idos = 0;
				fprintf(fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 4; k < 20; k++)
				{
					fprintf(fp, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos += this->orbit_dos[k][num - 1][i];
					tot_idos += this->Idos[k][num - 1][i];
				}
				fprintf(fp, "%lf     ", tot_dos);
				fprintf(_fp, "%lf     ", tot_idos);
				fprintf(fp, "\n");
				fprintf(_fp, "\n");
			}
		}
		fclose(fp);
		fclose(_fp);
		printf("Written %s FILE!\n", file);
		printf("Written %s FILE!\n", _file);
		string spa = convert_ADOS(file);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_file);
		printf("Written %s File!\n", _spa.c_str());
	}
	if (ISPIN == 2 && LORBIT == 10)
	{
		char file1[20], file2[20], _file1[20], _file2[20];
		if (elemlabel == NULL)
		{
			sprintf(file1, "PDOS_A%d_UP.dat", num);
			sprintf(file2, "PDOS_A%d_DW.dat", num);
			sprintf(_file1, "IPDOS_A%d_UP.dat", num);
			sprintf(_file2, "IPDOS_A%d_DW.dat", num);
		}
		else
		{
			sprintf(file1, "PDOS_%s%d_UP.dat", elemlabel, num - headind);
			sprintf(file2, "PDOS_%s%d_DW.dat", elemlabel, num - headind);
			sprintf(_file1, "IPDOS_%s%d_UP.dat", elemlabel, num - headind);
			sprintf(_file2, "IPDOS_%s%d_DW.dat", elemlabel, num - headind);
		}
		FILE* fp1 = fopen(file1, "w");
		FILE* fp2 = fopen(file2, "w");
		FILE* _fp1 = fopen(_file1, "w");
		FILE* _fp2 = fopen(_file2, "w");
		if (max_element <= 57)
		{
			fprintf(fp1, "     #Energy           s          p          d          tot\n");
			fprintf(_fp1, "     #Energy           s          p          d          tot\n");
			fprintf(fp2, "     #Energy           s          p          d          tot\n");
			fprintf(_fp2, "     #Energy           s          p          d          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos1 = 0;
				double tot_idos1 = 0;
				double tot_dos2 = 0;
				double tot_idos2 = 0;
				fprintf(fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 20; k < 23; k++)
				{
					fprintf(fp1, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp1, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos1 += this->orbit_dos[k][num - 1][i];
					tot_idos1 += this->Idos[k][num - 1][i];
				}
				fprintf(fp1, "%lf     ", tot_dos1);
				fprintf(_fp1, "%lf     ", tot_idos1);
				fprintf(fp1, "\n");
				fprintf(_fp1, "\n");
				for (int k = 24; k < 27; k++)
				{
					fprintf(fp2, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp2, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos2 += this->orbit_dos[k][num - 1][i];
					tot_idos2 += this->Idos[k][num - 1][i];
				}
				fprintf(fp2, "%lf     ", tot_dos2);
				fprintf(_fp2, "%lf     ", tot_idos2);
				fprintf(fp2, "\n");
				fprintf(_fp2, "\n");
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp1, "     #Energy           s          p          d          f          tot\n");
			fprintf(_fp1, "     #Energy           s          p          d          f          tot\n");
			fprintf(fp2, "     #Energy           s          p          d          f          tot\n");
			fprintf(_fp2, "     #Energy           s          p          d          f          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos1 = 0;
				double tot_idos1 = 0;
				double tot_dos2 = 0;
				double tot_idos2 = 0;
				fprintf(fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 20; k < 24; k++)
				{
					fprintf(fp1, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp1, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos1 += this->orbit_dos[k][num - 1][i];
					tot_idos1 += this->Idos[k][num - 1][i];
				}
				fprintf(fp1, "%lf     ", tot_dos1);
				fprintf(_fp1, "%lf     ", tot_idos1);
				fprintf(fp1, "\n");
				fprintf(_fp1, "\n");
				for (int k = 24; k < 28; k++)
				{
					fprintf(fp2, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp2, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos2 += this->orbit_dos[k][num - 1][i];
					tot_idos2 += this->Idos[k][num - 1][i];
				}
				fprintf(fp2, "%lf     ", tot_dos2);
				fprintf(_fp2, "%lf     ", tot_idos2);
				fprintf(fp2, "\n");
				fprintf(_fp2, "\n");
			}
		}
		fclose(fp1);
		fclose(fp2);
		fclose(_fp1);
		fclose(_fp2);
		printf("Written %s FILE!\n", file1);
		printf("Written %s FILE!\n", file2);
		printf("Written %s FILE!\n", _file1);
		printf("Written %s FILE!\n", _file2);
		string spa = convert_ADOS(file1,file2);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_file1, _file2);
		printf("Written %s File!\n", _spa.c_str());
	}
	if (ISPIN == 2 && LORBIT == 11)
	{
		char file1[20], file2[20], _file1[20], _file2[20];
		if (elemlabel == NULL)
		{
			sprintf(file1, "PDOS_A%d_UP.dat", num);
			sprintf(file2, "PDOS_A%d_DW.dat", num);
			sprintf(_file1, "IPDOS_A%d_UP.dat", num);
			sprintf(_file2, "IPDOS_A%d_DW.dat", num);
		}
		else
		{
			sprintf(file1, "PDOS_%s%d_UP.dat", elemlabel, num - headind);
			sprintf(file2, "PDOS_%s%d_DW.dat", elemlabel, num - headind);
			sprintf(_file1, "IPDOS_%s%d_UP.dat", elemlabel, num - headind);
			sprintf(_file2, "IPDOS_%s%d_DW.dat", elemlabel, num - headind);
		}
		FILE* fp1 = fopen(file1, "w");
		FILE* fp2 = fopen(file2, "w");
		FILE* _fp1 = fopen(_file1, "w");
		FILE* _fp2 = fopen(_file2, "w");
		if (max_element <= 57)
		{
			fprintf(fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			fprintf(_fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			fprintf(fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			fprintf(_fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos1 = 0;
				double tot_idos1 = 0;
				double tot_dos2 = 0;
				double tot_idos2 = 0;
				fprintf(fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 28; k < 37; k++)
				{
					fprintf(fp1, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp1, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos1 += this->orbit_dos[k][num - 1][i];
					tot_idos1 += this->Idos[k][num - 1][i];
				}
				fprintf(fp1, "%lf     ", tot_dos1);
				fprintf(_fp1, "%lf     ", tot_idos1);
				fprintf(fp1, "\n");
				fprintf(_fp1, "\n");
				for (int k = 44; k < 53; k++)
				{
					fprintf(fp2, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp2, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos2 += this->orbit_dos[k][num - 1][i];
					tot_idos2 += this->Idos[k][num - 1][i];
				}
				fprintf(fp2, "%lf     ", tot_dos2);
				fprintf(_fp2, "%lf     ", tot_idos2);
				fprintf(fp2, "\n");
				fprintf(_fp2, "\n");
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			fprintf(_fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			fprintf(fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			fprintf(_fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				double tot_dos1 = 0;
				double tot_idos1 = 0;
				double tot_dos2 = 0;
				double tot_idos2 = 0;
				fprintf(fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[num - 1][i] - this->efermi);
				for (int k = 28; k < 44; k++)
				{
					fprintf(fp1, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp1, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos1 += this->orbit_dos[k][num - 1][i];
					tot_idos1 += this->Idos[k][num - 1][i];
				}
				fprintf(fp1, "%lf     ", tot_dos1);
				fprintf(_fp1, "%lf     ", tot_idos1);
				fprintf(fp1, "\n");
				fprintf(_fp1, "\n");
				for (int k = 44; k < 60; k++)
				{
					fprintf(fp2, "%lf     ", this->orbit_dos[k][num - 1][i]);
					fprintf(_fp2, "%lf     ", this->Idos[k][num - 1][i]);
					tot_dos2 += this->orbit_dos[k][num - 1][i];
					tot_idos2 += this->Idos[k][num - 1][i];
				}
				fprintf(fp2, "%lf     ", tot_dos2);
				fprintf(_fp2, "%lf     ", tot_idos2);
				fprintf(fp2, "\n");
				fprintf(_fp2, "\n");
			}
		}
		fclose(fp1);
		fclose(fp2);
		fclose(_fp1);
		fclose(_fp2);
		printf("Written %s FILE!\n", file1);
		printf("Written %s FILE!\n", file2);
		printf("Written %s FILE!\n", _file1);
		printf("Written %s FILE!\n", _file2);
		string spa = convert_ADOS(file1, file2);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_file1, _file2);
		printf("Written %s File!\n", _spa.c_str());
	}
	fermi_energy();
}

void DOS::EDos(vector<string> element, int typenum[], char elemsym[][3], int nant[], int max_element, int ISPIN, int LORBIT)
{
	int cnt = 0;
	char(*file)[20];
	file = new char[nant[0]][20];
	char(*_file)[20];
	_file = new char[nant[0]][20];
	char(*file1)[20];
	file1 = new char[nant[0]][20];
	char(*_file1)[20];
	_file1 = new char[nant[0]][20];
	char(*file2)[20];
	file2 = new char[nant[0]][20];
	char(*_file2)[20];
	_file2 = new char[nant[0]][20];
	for (int i = 0; i < nant[1]; i++)
	{
		if (find(element.begin(), element.end(), elemsym[i]) == element.end())
		{
			cnt += typenum[i];
			continue;
		}
		else
		{
			if (ISPIN == 1 && LORBIT == 10)
			{
				sprintf(file[i], "PDOS_%s.dat", elemsym[i]);
				FILE* fp = fopen(file[i], "w");
				sprintf(_file[i], "IPDOS_%s.dat", elemsym[i]);
				FILE* _fp = fopen(_file[i], "w");
				if (max_element <= 57)
				{
					fprintf(fp, "     #Energy           s          p          d          tot\n");
					fprintf(_fp, "     #Energy           s          p          d          tot\n");
					for (int k = 0; k < this->NEDOS; k++)
					{
						double tot_dos = 0;
						double tot_idos = 0;
						fprintf(fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						for (int m = 0; m < 3; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp, "%lf     ", temp_dos);
							fprintf(_fp, "%lf     ", temp_idos);
						}
						fprintf(fp, "%lf     \n", tot_dos);
						fprintf(_fp, "%lf     \n", tot_idos);
					}
				}
				else if (max_element > 57)
				{
					fprintf(fp, "     #Energy           s          p          d          f          tot\n");
					fprintf(_fp, "     #Energy           s          p          d          f          tot\n");
					for (int k = 0; k < this->NEDOS; k++)
					{
						double tot_dos = 0;
						double tot_idos = 0;
						fprintf(fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						for (int m = 0; m < 4; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp, "%lf     ", temp_dos);
							fprintf(_fp, "%lf     ", temp_idos);
						}
						fprintf(fp, "%lf     \n", tot_dos);
						fprintf(_fp, "%lf     \n", tot_idos);
					}
				}
				fclose(fp);
				fclose(_fp);
				printf("Written %s file!\n", file[i]);
				printf("Written %s file!\n", _file[i]);
				string spa = convert_ADOS(file[i]);
				printf("Written %s File!\n", spa.c_str());
				string _spa = convert_ADOS(_file[i]);
				printf("Written %s File!\n", _spa.c_str());
			}
			if (ISPIN == 1 && LORBIT == 11)
			{
				sprintf(file[i], "PDOS_%s.dat", elemsym[i]);
				FILE* fp = fopen(file[i], "w");
				sprintf(_file[i], "IPDOS_%s.dat", elemsym[i]);
				FILE* _fp = fopen(_file[i], "w");
				if (max_element <= 57)
				{
					fprintf(fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
					fprintf(_fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
					for (int k = 0; k < this->NEDOS; k++)
					{
						double tot_dos = 0;
						double tot_idos = 0;
						fprintf(fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						for (int m = 4; m < 13; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp, "%lf     ", temp_dos);
							fprintf(_fp, "%lf     ", temp_idos);
						}
						fprintf(fp, "%lf     \n", tot_dos);
						fprintf(_fp, "%lf     \n", tot_idos);
					}
				}
				else if (max_element > 57)
				{
					fprintf(fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
					fprintf(_fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 4; m < 20; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp, "%lf     ", temp_dos);
							fprintf(_fp, "%lf     ", temp_idos);
						}
						fprintf(fp, "%lf     \n", tot_dos);
						fprintf(_fp, "%lf     \n", tot_idos);
					}
				}
				fclose(fp);
				fclose(_fp);
				printf("Written %s file!\n", file[i]);
				printf("Written %s file!\n", _file[i]);
				string spa = convert_ADOS(file[i]);
				printf("Written %s File!\n", spa.c_str());
				string _spa = convert_ADOS(_file[i]);
				printf("Written %s File!\n", _spa.c_str());
			}
			if (ISPIN == 2 && LORBIT == 10)
			{
				sprintf(file1[i], "PDOS_%s_UP.dat", elemsym[i]);
				FILE* fp1 = fopen(file1[i], "w");
				sprintf(file2[i], "PDOS_%s_DW.dat", elemsym[i]);
				FILE* fp2 = fopen(file2[i], "w");
				sprintf(_file1[i], "IPDOS_%s_UP.dat", elemsym[i]);
				FILE* _fp1 = fopen(_file1[i], "w");
				sprintf(_file2[i], "IPDOS_%s_DW.dat", elemsym[i]);
				FILE* _fp2 = fopen(_file2[i], "w");
				if (max_element <= 57)
				{
					fprintf(fp1, "     #Energy           s          p          d          tot\n");
					fprintf(_fp1, "     #Energy           s          p          d          tot\n");
					fprintf(fp2, "     #Energy           s          p          d          tot\n");
					fprintf(_fp2, "     #Energy           s          p          d          tot\n");
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 20; m < 23; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp1, "%lf     ", temp_dos);
							fprintf(_fp1, "%lf     ", temp_idos);
						}
						fprintf(fp1, "%lf     \n", tot_dos);
						fprintf(_fp1, "%lf     \n", tot_idos);
					}
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 24; m < 27; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp2, "%lf     ", temp_dos);
							fprintf(_fp2, "%lf     ", temp_idos);
						}
						fprintf(fp2, "%lf     \n", tot_dos);
						fprintf(_fp2, "%lf     \n", tot_idos);
					}
				}
				else if (max_element > 57)
				{
					fprintf(fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
					fprintf(_fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
					fprintf(fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
					fprintf(_fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 20; m < 24; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp1, "%lf     ", temp_dos);
							fprintf(_fp1, "%lf     ", temp_idos);
						}
						fprintf(fp1, "%lf     \n", tot_dos);
						fprintf(_fp1, "%lf     \n", tot_idos);
					}
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 24; m < 28; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp2, "%lf     ", temp_dos);
							fprintf(_fp2, "%lf     ", temp_idos);
						}
						fprintf(fp2, "%lf     \n", tot_dos);
						fprintf(_fp2, "%lf     \n", tot_idos);
					}
				}
				fclose(fp1);
				fclose(_fp1);
				fclose(fp2);
				fclose(_fp2);
				printf("Written %s file!\n", file1[i]);
				printf("Written %s file!\n", file2[i]);
				printf("Written %s file!\n", _file1[i]);
				printf("Written %s file!\n", _file2[i]);
				string spa = convert_ADOS(file1[i], file2[i]);
				printf("Written %s File!\n", spa.c_str());
				string _spa = convert_ADOS(_file1[i], _file2[i]);
				printf("Written %s File!\n", _spa.c_str());
			}
			if (ISPIN == 2 && LORBIT == 11)
			{
				sprintf(file1[i], "PDOS_%s_UP.dat", elemsym[i]);
				FILE* fp1 = fopen(file1[i], "w");
				sprintf(file2[i], "PDOS_%s_DW.dat", elemsym[i]);
				FILE* fp2 = fopen(file2[i], "w");
				sprintf(_file1[i], "IPDOS_%s_UP.dat", elemsym[i]);
				FILE* _fp1 = fopen(_file1[i], "w");
				sprintf(_file2[i], "IPDOS_%s_DW.dat", elemsym[i]);
				FILE* _fp2 = fopen(_file2[i], "w");
				fprintf(fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
				fprintf(_fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
				fprintf(fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
				fprintf(_fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
				if (max_element <= 57)
				{
					for (int k = 0; k < this->NEDOS; k++)
					{
						double tot_dos = 0;
						double tot_idos = 0;
						fprintf(fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						for (int m = 28; m < 37; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp1, "%lf     ", temp_dos);
							fprintf(_fp1, "%lf     ", temp_idos);
						}
						fprintf(fp1, "%lf     \n", tot_dos);
						fprintf(_fp1, "%lf     \n", tot_idos);
					}
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 44; m < 53; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp2, "%lf     ", temp_dos);
							fprintf(_fp2, "%lf     ", temp_idos);
						}
						fprintf(fp2, "%lf     \n", tot_dos);
						fprintf(_fp2, "%lf     \n", tot_idos);
					}
				}
				else if (max_element > 57)
				{
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp1, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 28; m < 44; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp1, "%lf     ", temp_dos);
							fprintf(_fp1, "%lf     ", temp_idos);
						}
						fprintf(fp1, "%lf     \n", tot_dos);
						fprintf(_fp1, "%lf     \n", tot_idos);
					}
					for (int k = 0; k < this->NEDOS; k++)
					{
						fprintf(fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						fprintf(_fp2, "%lf     ", this->energy_ion[0][k] - this->efermi);
						double tot_dos = 0;
						double tot_idos = 0;
						for (int m = 44; m < 60; m++)
						{
							double temp_dos = 0;
							double temp_idos = 0;
							for (int j = 0; j < typenum[i]; j++)
							{
								temp_dos += this->orbit_dos[m][j + cnt][k];
								temp_idos += this->Idos[m][j + cnt][k];
							}
							tot_dos += temp_dos;
							tot_idos += temp_idos;
							fprintf(fp2, "%lf     ", temp_dos);
							fprintf(_fp2, "%lf     ", temp_idos);
						}
						fprintf(fp2, "%lf     \n", tot_dos);
						fprintf(_fp2, "%lf     \n", tot_idos);
					}
				}
				fclose(fp1);
				fclose(_fp1);
				fclose(fp2);
				fclose(_fp2);
				printf("Written %s file!\n", file1[i]);
				printf("Written %s file!\n", file2[i]);
				printf("Written %s file!\n", _file1[i]);
				printf("Written %s file!\n", _file2[i]);
				string spa = convert_ADOS(file1[i], file2[i]);
				printf("Written %s File!\n", spa.c_str());
				string _spa = convert_ADOS(_file1[i], _file2[i]);
				printf("Written %s File!\n", _spa.c_str());
			}
		}
		cnt += typenum[i];
	}
	fermi_energy();
	delete[] file;
	delete[] _file;
	delete[] file1;
	delete[] _file1;
	delete[] file2;
	delete[] _file2;
}

void DOS::SDos(vector<int> num, string postfix, int max_element, int ISPIN, int LORBIT)
{
	if (ISPIN == 1 && LORBIT == 10)
	{
		char sum_file[50];
		sprintf(sum_file, "PDOS_SUM_%s.dat", postfix.c_str());
		FILE* fp = fopen(sum_file, "w");
		char _sum_file[50];
		sprintf(_sum_file, "IPDOS_SUM_%s.dat", postfix.c_str());
		FILE* _fp = fopen(_sum_file, "w");
		if (max_element <= 57)
		{
			fprintf(fp, "     #Energy           s          p          d          tot\n");
			fprintf(_fp, "     #Energy           s          p          d          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 0; k < 3; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp, "%lf     ", temp_dos);
					fprintf(_fp, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp, "%lf     \n", tot_dos);
				fprintf(_fp, "%lf     \n", tot_idos);
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp, "     #Energy           s          p          d          f          tot\n");
			fprintf(_fp, "     #Energy           s          p          d          f          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 0; k < 4; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp, "%lf     ", temp_dos);
					fprintf(_fp, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp, "%lf     \n", tot_dos);
				fprintf(_fp, "%lf     \n", tot_idos);
			}
		}
		printf("Written %s file!\n", sum_file);
		printf("Written %s file!\n", _sum_file);
		string spa = convert_ADOS(sum_file);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_sum_file);
		printf("Written %s File!\n", _spa.c_str());
		fclose(fp);
		fclose(_fp);
	}
	if (ISPIN == 1 && LORBIT == 11)
	{
		char sum_file[50];
		sprintf(sum_file, "PDOS_SUM_%s.dat", postfix.c_str());
		FILE* fp = fopen(sum_file, "w");
		char _sum_file[50];
		sprintf(_sum_file, "IPDOS_SUM_%s.dat", postfix.c_str());
		FILE* _fp = fopen(_sum_file, "w");
		if (max_element <= 57)
		{
			fprintf(fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			fprintf(_fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 4; k < 13; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp, "%lf     ", temp_dos);
					fprintf(_fp, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp, "%lf     \n", tot_dos);
				fprintf(_fp, "%lf     \n", tot_idos);
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			fprintf(_fp, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 4; k < 19; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp, "%lf     ", temp_dos);
					fprintf(_fp, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp, "%lf     \n", tot_dos);
				fprintf(_fp, "%lf     \n", tot_idos);
			}
		}
		printf("Written %s file!\n", sum_file);
		printf("Written %s file!\n", _sum_file);
		string spa = convert_ADOS(sum_file);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_sum_file);
		printf("Written %s File!\n", _spa.c_str());
		fclose(fp);
		fclose(_fp);
	}
	if (ISPIN == 2 && LORBIT == 10)
	{
		char sum_file1[50], sum_file2[50], _sum_file1[50], _sum_file2[50];
		sprintf(sum_file1, "PDOS_SUM_UP_%s.dat", postfix.c_str());
		FILE* fp1 = fopen(sum_file1, "w");
		sprintf(sum_file2, "PDOS_SUM_DW_%s.dat", postfix.c_str());
		FILE* fp2 = fopen(sum_file2, "w");
		sprintf(_sum_file1, "IPDOS_SUM_UP_%s.dat", postfix.c_str());
		FILE* _fp1 = fopen(_sum_file1, "w");
		sprintf(_sum_file2, "IPDOS_SUM_UP_%s.dat", postfix.c_str());
		FILE* _fp2 = fopen(_sum_file2, "w");
		if (max_element <= 57)
		{
			fprintf(fp1, "     #Energy           s          p          d          tot\n");
			fprintf(_fp1, "     #Energy           s          p          d          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 20; k < 23; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp1, "%lf     ", temp_dos);
					fprintf(_fp1, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp1, "%lf     \n", tot_dos);
				fprintf(_fp1, "%lf     \n", tot_idos);
			}
			fprintf(fp2, "     #Energy           s          p          d          tot\n");
			fprintf(_fp2, "     #Energy           s          p          d          tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 24; k < 27; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp2, "%lf     ", temp_dos);
					fprintf(_fp2, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp2, "%lf     \n", tot_dos);
				fprintf(_fp2, "%lf     \n", tot_idos);
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp1, "     #Energy           s          p          d           f         tot\n");
			fprintf(_fp1, "     #Energy           s          p          d           f         tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 20; k < 24; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp1, "%lf     ", temp_dos);
					fprintf(_fp1, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp1, "%lf     \n", tot_dos);
				fprintf(_fp1, "%lf     \n", tot_idos);
			}
			fprintf(fp2, "     #Energy           s          p          d           f         tot\n");
			fprintf(_fp2, "     #Energy           s          p          d           f         tot\n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 24; k < 28; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp2, "%lf     ", temp_dos);
					fprintf(_fp2, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp2, "%lf     \n", tot_dos);
				fprintf(_fp2, "%lf     \n", tot_idos);
			}
		}
		printf("Written %s file!\n", sum_file1);
		printf("Written %s file!\n", sum_file2);
		printf("Written %s file!\n", _sum_file1);
		printf("Written %s file!\n", _sum_file2);
		string spa = convert_ADOS(sum_file1, sum_file2);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_sum_file1, _sum_file2);
		printf("Written %s File!\n", _spa.c_str());
		fclose(fp1);
		fclose(fp2);
		fclose(_fp1);
		fclose(_fp2);
	}
	if (ISPIN == 2 && LORBIT == 11)
	{
		char sum_file1[50], sum_file2[50], _sum_file1[50], _sum_file2[50];
		sprintf(sum_file1, "PDOS_SUM_UP_%s.dat", postfix.c_str());
		FILE* fp1 = fopen(sum_file1, "w");
		sprintf(sum_file2, "PDOS_SUM_DW_%s.dat", postfix.c_str());
		FILE* fp2 = fopen(sum_file2, "w");
		sprintf(_sum_file1, "IPDOS_SUM_UP_%s.dat", postfix.c_str());
		FILE* _fp1 = fopen(_sum_file1, "w");
		sprintf(_sum_file2, "IPDOS_SUM_UP_%s.dat", postfix.c_str());
		FILE* _fp2 = fopen(_sum_file2, "w");
		if (max_element <= 57)
		{
			fprintf(fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			fprintf(_fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 28; k < 37; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp1, "%lf     ", temp_dos);
					fprintf(_fp1, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp1, "%lf     \n", tot_dos);
				fprintf(_fp1, "%lf     \n", tot_idos);
			}
			fprintf(fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			fprintf(_fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       tot  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 44; k < 53; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp2, "%lf     ", temp_dos);
					fprintf(_fp2, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp2, "%lf     \n", tot_dos);
				fprintf(_fp2, "%lf     \n", tot_idos);
			}
		}
		else if (max_element > 57)
		{
			fprintf(fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			fprintf(_fp1, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp1, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 28; k < 44; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp1, "%lf     ", temp_dos);
					fprintf(_fp1, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp1, "%lf     \n", tot_dos);
				fprintf(_fp1, "%lf     \n", tot_idos);
			}
			fprintf(fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			fprintf(_fp2, "     #Energy           s          py          pz          px         dxy        dyz         dz2         dxz         dx2       f1        f2         f3         f4         f5         f6         f7         total  \n");
			for (int i = 0; i < this->NEDOS; i++)
			{
				fprintf(fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				fprintf(_fp2, "%lf     ", this->energy_ion[0][i] - this->efermi);
				double tot_dos = 0;
				double tot_idos = 0;
				for (int k = 44; k < 60; k++)
				{
					double temp_dos = 0;
					double temp_idos = 0;
					for (int j = 0; j < num.size(); j++)
					{
						temp_dos += this->orbit_dos[k][num[j] - 1][i];
						temp_idos += this->Idos[k][num[j] - 1][i];
					}
					fprintf(fp2, "%lf     ", temp_dos);
					fprintf(_fp2, "%lf     ", temp_idos);
					tot_dos += temp_dos;
					tot_idos += temp_idos;
				}
				fprintf(fp2, "%lf     \n", tot_dos);
				fprintf(_fp2, "%lf     \n", tot_idos);
			}
		}
		printf("Written %s file!\n", sum_file1);
		printf("Written %s file!\n", sum_file2);
		printf("Written %s file!\n", _sum_file1);
		printf("Written %s file!\n", _sum_file2);
		string spa = convert_ADOS(sum_file1, sum_file2);
		printf("Written %s File!\n", spa.c_str());
		string _spa = convert_ADOS(_sum_file1, _sum_file2);
		printf("Written %s File!\n", _spa.c_str());
		fclose(fp1);
		fclose(fp2);
		fclose(_fp1);
		fclose(_fp2);
	}
	fermi_energy();
}

void DOS::ODos(vector<vector<string> > ind, vector<vector<string> > orb, vector<COMB> orbit, int typenum[], char elemsym[][3], int nant[],
	int max_element, int ISPIN, int LORBIT)
{
	map<string, vector<vector<double> > > data1 =
	{
		{"s", this->s_dos},{"p", this->p_dos},{"d", this->d_dos},{"f",this->f_dos}
	};
	map<string, vector<vector<double> > > data2 =
	{
		{"s", this->s0_dos},{"px", this->px_dos},{"py", this->py_dos},{"pz", this->pz_dos},{"dxy",this->dxy_dos},
		{"dxz",this->dxz_dos},{"dyz",this->dyz_dos},{"dx2",this->dx2_dos},{"dz2",this->dz2_dos},{"f1",this->f1_dos},
		{"f2",this->f2_dos},{"f3",this->f3_dos},{"f4",this->f4_dos},{"f5",this->f5_dos},{"f6",this->f6_dos},{"f7",this->f7_dos}
	};
	map<string, vector<vector<double> > > data3_up =
	{
		{"s", this->s_dos_up},{"p", this->p_dos_up},{"d", this->d_dos_up},{"f",this->f_dos_up}
	};
	map<string, vector<vector<double> > > data3_dw =
	{
		{"s", this->s_dos_dw},{"p", this->p_dos_dw},{"d", this->d_dos_dw},{"f",this->f_dos_dw}
	};
	map<string, vector<vector<double> > > data4_up =
	{
		{"s", this->s0_dos_up},{"px", this->px_dos_up},{"py", this->py_dos_up},{"pz", this->pz_dos_up},{"dxy",this->dxy_dos_up},
		{"dxz",this->dxz_dos_up},{"dyz",this->dyz_dos_up},{"dx2",this->dx2_dos_up},{"dz2",this->dz2_dos_up},{"f1",this->f1_dos_up},
		{"f2",this->f2_dos_up},{"f3",this->f3_dos_up},{"f4",this->f4_dos_up},{"f5",this->f5_dos_up},{"f6",this->f6_dos_up},{"f7",this->f7_dos_up}
	};
	map<string, vector<vector<double> > > data4_dw =
	{
		{"s", this->s0_dos_dw},{"px", this->px_dos_dw},{"py", this->py_dos_dw},{"pz", this->pz_dos_dw},{"dxy",this->dxy_dos_dw},
		{"dxz",this->dxz_dos_dw},{"dyz",this->dyz_dos_dw},{"dx2",this->dx2_dos_dw},{"dz2",this->dz2_dos_dw},{"f1",this->f1_dos_dw},
		{"f2",this->f2_dos_dw},{"f3",this->f3_dos_dw},{"f4",this->f4_dos_dw},{"f5",this->f5_dos_dw},{"f6",this->f6_dos_dw},{"f7",this->f7_dos_dw}
	};
	vector<vector<double> > data;	//data[orbit.size()][NEDOS]
	if (ISPIN == 1)
		data.resize(orbit.size());
	else
		data.resize(orbit.size() * 2);
	for (int i = 0; i < data.size(); i++)
		data[i].resize(this->NEDOS);
	vector<string> title;
	FILE* fp = fopen("PDOS_USER.dat", "w");
	fprintf(fp, "	#Energy");
	//Symbol [1 2 s][3-5 s px py][6-10 all]... 
	// 1&2_s 3-5_s&px&py 6-10_all
	if (ISPIN == 1)
	{
		for (int i = 0; i < orbit.size(); i++)
		{
			string str;
			for (int k = 0; k < ind[i].size(); k++)
				str += ind[i][k] + (k == ind[i].size() - 1 ? "" : "&");
			str += "_";
			for (int j = 0; j < orb[i].size(); j++)
				str += orb[i][j] + (j == orb[i].size() - 1 ? "" : "&");
			title.push_back(str);
		}
		for (int i = 0; i < title.size(); i++)
			fprintf(fp, "	%s", title[i].c_str());
	}
	else
	{
		for (int i = 0; i < orbit.size(); i++)
		{
			string str1 = "up_";
			string str2 = "dw_";
			for (int k = 0; k < ind[i].size(); k++)
			{
				str1 += ind[i][k] + (k == ind[i].size() - 1 ? "" : "&");
				str2 += ind[i][k] + (k == ind[i].size() - 1 ? "" : "&");
			}
			str1 += "_";
			str2 += "_";
			for (int j = 0; j < orb[i].size(); j++)
			{
				str1 += orb[i][j] + (j == orb[i].size() - 1 ? "" : "&");
				str2 += orb[i][j] + (j == orb[i].size() - 1 ? "" : "&");
			}
			title.push_back(str1);
			title.push_back(str2);
		}
		for (int i = 0; i < title.size(); i++)
			fprintf(fp, "	%s", title[i].c_str());
	}
	fprintf(fp, "\n");
	for (int i = 0; i < orbit.size(); i++)
	{
		if (ISPIN == 1 && LORBIT == 10)
			for (int j = 0; j < orbit[i].orb.size(); j++)
				for (int k = 0; k < orbit[i].index.size(); k++)
					data[i] = AddVec(data[i], data1[orbit[i].orb[j]][orbit[i].index[k] - 1]);
		if (ISPIN == 1 && LORBIT == 11)
			for (int j = 0; j < orbit[i].orb.size(); j++)
				for (int k = 0; k < orbit[i].index.size(); k++)
					data[i] = AddVec(data[i], data2[orbit[i].orb[j]][orbit[i].index[k] - 1]);
		if (ISPIN == 2 && LORBIT == 10)
			for (int j = 0; j < orbit[i].orb.size(); j++)
				for (int k = 0; k < orbit[i].index.size(); k++)
				{
					data[i] = AddVec(data[i], data3_up[orbit[i].orb[j]][orbit[i].index[k] - 1]);
					data[i + orbit.size()] = AddVec(data[i + orbit.size()], data3_dw[orbit[i].orb[j]][orbit[i].index[k] - 1]);
				}
		if (ISPIN == 2 && LORBIT == 11)
			for (int j = 0; j < orbit[i].orb.size(); j++)
				for (int k = 0; k < orbit[i].index.size(); k++)
				{
					data[i] = AddVec(data[i], data4_up[orbit[i].orb[j]][orbit[i].index[k] - 1]);
					data[i + orbit.size()] = AddVec(data[i + orbit.size()], data4_dw[orbit[i].orb[j]][orbit[i].index[k] - 1]);
				}
	}
	if (ISPIN == 1)
	{
		for (int i = 0; i < this->NEDOS; i++)
		{
			fprintf(fp, "	%lf", this->energy_ion[0][i] - this->efermi);
			for (int j = 0; j < orbit.size(); j++)
			{
				fprintf(fp, "	%lf", data[j][i]);
			}
			fprintf(fp, "\n");
		}
	}
	else
	{
		for (int i = 0; i < this->NEDOS; i++)
		{
			fprintf(fp, "	%lf", this->energy_ion[0][i] - this->efermi);
			for (int j = 0; j < orbit.size(); j++)
			{
				fprintf(fp, "	%lf", data[j][i]);
				fprintf(fp, "	%lf", data[j + orbit.size()][i]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	printf("Written PDOS_USER.dat file!\n");
	string spa = convert_ADOS("PDOS_USER.dat");
	printf("Written %s File!\n", spa.c_str());
	fermi_energy();
}

void DOS::BandCenter(vector<int> num, string atomid, int max_element, int ISPIN, int LORBIT)
{
	FILE* fp = fopen("Band_Center", "w");
	if (ISPIN == 1 && LORBIT == 10)
	{
		if (max_element <= 57)
		{
			vector<vector<double> > dos_sum(3, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum(3, vector<double>(this->NEDOS));
			double center[3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < num.size(); j++)
				{
					dos_sum[i] = this->orbit_dos[i][num[j] - 1] + dos_sum[i];
					Idos_sum[i] = this->Idos[i][num[j] - 1] + Idos_sum[i];
				}
			for (int i = 0; i < 3; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum[i]);
				center[i] = temp[this->NEDOS - 1] / Idos_sum[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	Atom_ID   s-band-center    p-band-center    d-band-center  (in units of eV)\n");
			fprintf(fp, "	%s	%lf	%lf	%lf\n", atomid.c_str(), center[0], center[1], center[2]);
		}
		else if (max_element > 57)
		{
			vector<vector<double> > dos_sum(4, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum(4, vector<double>(this->NEDOS));
			double center[4];
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < num.size(); j++)
				{
					dos_sum[i] = this->orbit_dos[i][num[j] - 1] + dos_sum[i];
					Idos_sum[i] = this->Idos[i][num[j] - 1] + Idos_sum[i];
				}
			for (int i = 0; i < 4; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum[i]);
				center[i] = temp[this->NEDOS - 1] / Idos_sum[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	Atom_ID   s-band-center    p-band-center    d-band-center    f-band-center  (in units of eV)\n");
			fprintf(fp, "	%s	%lf	%lf	%lf %lf\n", atomid.c_str(), center[0], center[1], center[2], center[3]);
		}
	}
	if (ISPIN == 1 && LORBIT == 11)
	{
		if (max_element <= 57)
		{
			vector<vector<double> > dos_sum(3, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum(3, vector<double>(this->NEDOS));
			vector<vector<double> > sp_dos_sum(9, vector<double>(this->NEDOS));
			vector<vector<double> > sp_Idos_sum(9, vector<double>(this->NEDOS));
			double center[3], sp_center[9];
			for (int i = 0; i < num.size(); i++)
			{
				dos_sum[0] = this->orbit_dos[4][num[i] - 1] + dos_sum[0]; 		//s
				sp_dos_sum[0] = this->orbit_dos[4][num[i] - 1] + sp_dos_sum[0];	//s0
				Idos_sum[0] = this->Idos[4][num[i] - 1] + Idos_sum[0];
				sp_Idos_sum[0] = this->Idos[4][num[i] - 1] + sp_Idos_sum[0];
				for (int j = 5, k = 1; j < 8; j++, k++)
				{
					dos_sum[1] = this->orbit_dos[j][num[i] - 1] + dos_sum[1];	//p
					Idos_sum[1] = this->Idos[j][num[i] - 1] + Idos_sum[1];
					sp_dos_sum[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum[k]; //py pz px
					sp_Idos_sum[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum[k];
				}
				for (int j = 8, k = 4; j < 13; j++, k++)
				{
					dos_sum[2] = this->orbit_dos[j][num[i] - 1] + dos_sum[2];	//d
					Idos_sum[2] = this->Idos[j][num[i] - 1] + Idos_sum[2];
					sp_dos_sum[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum[k]; //dxy, dyz, dz2, dxz, dx2,
					sp_Idos_sum[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum[k];
				}
			}
			for (int i = 0; i < 3; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum[i]);
				center[i] = temp[this->NEDOS - 1] / Idos_sum[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	Atom_ID   s-band-center    p-band-center    d-band-center  (in units of eV)\n");
			fprintf(fp, "	%s	%lf	%lf	%lf\n\n", atomid.c_str(), center[0], center[1], center[2]);
			fprintf(fp, "#	Atom_ID   s   py   pz   px   dxy   dyz   dz2   dxz   dx2   (in units of eV)\n");
			fprintf(fp, "	%s", atomid.c_str());
			for (int i = 0; i < 9; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * sp_dos_sum[i]);
				sp_center[i] = temp[this->NEDOS - 1] / sp_Idos_sum[i][this->NEDOS - 1];
				fprintf(fp, "	%lf", sp_center[i]);
			}
			fprintf(fp, "\n");
		}
		else if (max_element > 57)
		{
			vector<vector<double> > dos_sum(4, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum(4, vector<double>(this->NEDOS));
			vector<vector<double> > sp_dos_sum(16, vector<double>(this->NEDOS));
			vector<vector<double> > sp_Idos_sum(16, vector<double>(this->NEDOS));
			double center[4], sp_center[11];
			for (int i = 0; i < num.size(); i++)
			{
				dos_sum[0] = this->orbit_dos[4][num[i] - 1] + dos_sum[0];
				Idos_sum[0] = this->Idos[4][num[i] - 1] + Idos_sum[0];
				sp_dos_sum[0] = this->orbit_dos[4][num[i] - 1] + sp_dos_sum[0];
				sp_Idos_sum[0] = this->Idos[4][num[i] - 1] + sp_Idos_sum[0];
				for (int j = 5, k = 1; j < 8; j++, k++)
				{
					dos_sum[1] = this->orbit_dos[j][num[i] - 1] + dos_sum[1];
					Idos_sum[1] = this->Idos[j][num[i] - 1] + Idos_sum[1];
					sp_dos_sum[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum[k]; //py pz px
					sp_Idos_sum[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum[k];
				}
				for (int j = 8, k = 4; j < 13; j++, k++)
				{
					dos_sum[2] = this->orbit_dos[j][num[i] - 1] + dos_sum[2];
					Idos_sum[2] = this->Idos[j][num[i] - 1] + Idos_sum[2];
					sp_dos_sum[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum[k]; //py pz px
					sp_Idos_sum[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum[k];
				}
				for (int j = 13, k = 9; j < 20; j++, k++)
				{
					dos_sum[3] = this->orbit_dos[j][num[i] - 1] + dos_sum[3];
					Idos_sum[3] = this->Idos[j][num[i] - 1] + Idos_sum[3];
					sp_dos_sum[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum[k]; //py pz px
					sp_Idos_sum[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum[k];
				}
			}
			for (int i = 0; i < 4; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum[i]);
				center[i] = temp[this->NEDOS - 1] / Idos_sum[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	s-band-center    p-band-center    d-band-center    f-band-center  (in units of eV)\n\n");
			fprintf(fp, "	%s	%lf	%lf	%lf %lf\n\n", atomid.c_str(), center[0], center[1], center[2], center[3]);
			fprintf(fp, "#	Atom_ID   s   py   pz   px   dxy   dyz   dz2   dxz   dx2   f1   f2   f3   f4   f5	f6	f7(in units of eV)\n");
			fprintf(fp, "	%s", atomid.c_str());
			for (int i = 0; i < 16; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * sp_dos_sum[i]);
				sp_center[i] = temp[this->NEDOS - 1] / sp_Idos_sum[i][this->NEDOS - 1];
				fprintf(fp, "	%lf", sp_center[i]);
			}
			fprintf(fp, "\n");
		}
	}
	if (ISPIN == 2 && LORBIT == 10)
	{
		if (max_element <= 57)
		{
			vector<vector<double> > dos_sum_up(3, vector<double>(this->NEDOS)), dos_sum_dw(3, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum_up(3, vector<double>(this->NEDOS)), Idos_sum_dw(3, vector<double>(this->NEDOS));
			double center_up[3];
			double center_dw[3];
			for (int i = 20; i < 23; i++)
				for (int j = 0; j < num.size(); j++)
				{
					dos_sum_up[i] = this->orbit_dos[i][num[j] - 1] + dos_sum_up[i];
					Idos_sum_up[i] = this->Idos[i][num[j] - 1] + Idos_sum_up[i];
				}
			for (int i = 24; i < 27; i++)
				for (int j = 0; j < num.size(); j++)
				{
					dos_sum_dw[i] = this->orbit_dos[i][num[j] - 1] + dos_sum_dw[i];
					Idos_sum_dw[i] = this->Idos[i][num[j] - 1] + Idos_sum_dw[i];
				}
			for (int i = 0; i < 3; i++)
			{
				vector<double> temp_up = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_up[i]);
				center_up[i] = temp_up[this->NEDOS - 1] / Idos_sum_up[i][this->NEDOS - 1];
				vector<double> temp_dw = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_dw[i]);
				center_dw[i] = temp_dw[this->NEDOS - 1] / Idos_sum_dw[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	Atom_ID   s-band-center    p-band-center    d-band-center  (in units of eV)\n");
			fprintf(fp, "	SPIN-UP	%s	%lf	%lf	%lf\n", atomid.c_str(), center_up[0], center_up[1], center_up[2]);
			fprintf(fp, "	SPIN-DW	%s	%lf	%lf	%lf\n", atomid.c_str(), center_dw[0], center_dw[1], center_dw[2]);
		}
		else if (max_element > 57)
		{
			vector<vector<double> > dos_sum_up(4, vector<double>(this->NEDOS)), dos_sum_dw(4, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum_up(4, vector<double>(this->NEDOS)), Idos_sum_dw(4, vector<double>(this->NEDOS));
			double center_up[4];
			double center_dw[4];
			for (int i = 20; i < 24; i++)
				for (int j = 0; j < num.size(); j++)
				{
					dos_sum_up[i] = this->orbit_dos[i][num[j] - 1] + dos_sum_up[i];
					Idos_sum_up[i] = this->Idos[i][num[j] - 1] + Idos_sum_up[i];
				}
			for (int i = 24; i < 28; i++)
				for (int j = 0; j < num.size(); j++)
				{
					dos_sum_dw[i] = this->orbit_dos[i][num[j] - 1] + dos_sum_dw[i];
					Idos_sum_dw[i] = this->Idos[i][num[j] - 1] + Idos_sum_dw[i];
				}
			for (int i = 0; i < 4; i++)
			{
				vector<double> temp_up = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_up[i]);
				center_up[i] = temp_up[this->NEDOS - 1] / Idos_sum_up[i][this->NEDOS - 1];
				vector<double> temp_dw = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_dw[i]);
				center_dw[i] = temp_dw[this->NEDOS - 1] / Idos_sum_dw[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	Atom_ID   s-band-center    p-band-center    d-band-center      f-band-center	(in units of eV)\n");
			fprintf(fp, "	SPIN-UP	%s	%lf	%lf	%lf	%lf\n", atomid.c_str(), center_up[0], center_up[1], center_up[2], center_up[3]);
			fprintf(fp, "	SPIN-DW	%s	%lf	%lf	%lf	%lf\n", atomid.c_str(), center_dw[0], center_dw[1], center_dw[2], center_dw[3]);
		}
	}
	if (ISPIN == 2 && LORBIT == 11)
	{
		if (max_element <= 57)
		{
			vector<vector<double> > dos_sum_up(3, vector<double>(this->NEDOS)), dos_sum_dw(3, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum_up(3, vector<double>(this->NEDOS)), Idos_sum_dw(3, vector<double>(this->NEDOS));
			vector<vector<double> > sp_dos_sum_up(9, vector<double>(this->NEDOS)), sp_dos_sum_dw(9, vector<double>(this->NEDOS));
			vector<vector<double> > sp_Idos_sum_up(9, vector<double>(this->NEDOS)), sp_Idos_sum_dw(9, vector<double>(this->NEDOS));
			double center_up[3], center_dw[3], sp_center_up[9], sp_center_dw[9];
			for (int i = 0; i < num.size(); i++)
			{
				dos_sum_up[0] = this->orbit_dos[28][num[i] - 1] + dos_sum_up[0];
				Idos_sum_up[0] = this->Idos[28][num[i] - 1] + Idos_sum_up[0];
				sp_dos_sum_up[0] = this->orbit_dos[28][num[i] - 1] + sp_dos_sum_up[0];
				sp_Idos_sum_up[0] = this->Idos[28][num[i] - 1] + sp_Idos_sum_up[0];
				for (int j = 29, k = 1; j < 32; j++, k++)
				{
					dos_sum_up[1] = this->orbit_dos[j][num[i] - 1] + dos_sum_up[1];
					Idos_sum_up[1] = this->Idos[j][num[i] - 1] + Idos_sum_up[1];
					sp_dos_sum_up[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_up[k];
					sp_Idos_sum_up[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_up[k];
				}
				for (int j = 32, k = 4; j < 37; j++, k++)
				{
					dos_sum_up[2] = this->orbit_dos[j][num[i] - 1] + dos_sum_up[2];
					Idos_sum_up[2] = this->Idos[j][num[i] - 1] + Idos_sum_up[2];
					sp_dos_sum_up[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_up[k];
					sp_Idos_sum_up[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_up[k];
				}
			}
			for (int i = 0; i < num.size(); i++)
			{
				dos_sum_dw[0] = this->orbit_dos[44][num[i] - 1] + dos_sum_dw[0];
				Idos_sum_dw[0] = this->Idos[44][num[i] - 1] + Idos_sum_dw[0];
				sp_dos_sum_dw[0] = this->orbit_dos[28][num[i] - 1] + sp_dos_sum_dw[0];
				sp_Idos_sum_dw[0] = this->Idos[28][num[i] - 1] + sp_Idos_sum_dw[0];
				for (int j = 45, k = 1; j < 48; j++, k++)
				{
					dos_sum_dw[1] = this->orbit_dos[j][num[i] - 1] + dos_sum_dw[1];
					Idos_sum_dw[1] = this->Idos[j][num[i] - 1] + Idos_sum_dw[1];
					sp_dos_sum_dw[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_dw[k];
					sp_Idos_sum_dw[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_dw[k];
				}
				for (int j = 48, k = 4; j < 53; j++, k++)
				{
					dos_sum_dw[2] = this->orbit_dos[j][num[i] - 1] + dos_sum_dw[2];
					Idos_sum_dw[2] = this->Idos[j][num[i] - 1] + Idos_sum_dw[2];
					sp_dos_sum_dw[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_dw[k];
					sp_Idos_sum_dw[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_dw[k];
				}
			}
			for (int i = 0; i < 3; i++)
			{
				vector<double> temp_up = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_up[i]);
				center_up[i] = temp_up[this->NEDOS - 1] / Idos_sum_up[i][this->NEDOS - 1];
				vector<double> temp_dw = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_dw[i]);
				center_dw[i] = temp_dw[this->NEDOS - 1] / Idos_sum_dw[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	Atom_ID   s-band-center    p-band-center    d-band-center  (in units of eV)\n");
			fprintf(fp, "SPIN-UP	%s	%lf	%lf	%lf\n", atomid.c_str(), center_up[0], center_up[1], center_up[2]);
			fprintf(fp, "SPIN-DW	%s	%lf	%lf	%lf\n\n", atomid.c_str(), center_dw[0], center_dw[1], center_dw[2]);
			fprintf(fp, "#	Atom_ID   s   py   pz   px   dxy   dyz   dz2   dxz   dx2   (in units of eV)\n");
			fprintf(fp, "SPIN-UP	%s", atomid.c_str());
			for (int i = 0; i < 9; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * sp_dos_sum_up[i]);
				sp_center_up[i] = temp[this->NEDOS - 1] / sp_Idos_sum_up[i][this->NEDOS - 1];
				fprintf(fp, "	%lf", sp_center_up[i]);
			}
			fprintf(fp, "\n");
			fprintf(fp, "SPIN-DW	%s", atomid.c_str());
			for (int i = 0; i < 9; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * sp_dos_sum_dw[i]);
				sp_center_dw[i] = temp[this->NEDOS - 1] / sp_Idos_sum_dw[i][this->NEDOS - 1];
				fprintf(fp, "	%lf", sp_center_dw[i]);
			}
			fprintf(fp, "\n");
		}
		else if (max_element > 57)
		{
			vector<vector<double> > dos_sum_up(4, vector<double>(this->NEDOS)), dos_sum_dw(4, vector<double>(this->NEDOS));
			vector<vector<double> > Idos_sum_up(4, vector<double>(this->NEDOS)), Idos_sum_dw(4, vector<double>(this->NEDOS));
			vector<vector<double> > sp_dos_sum_up(16, vector<double>(this->NEDOS)), sp_dos_sum_dw(16, vector<double>(this->NEDOS));
			vector<vector<double> > sp_Idos_sum_up(16, vector<double>(this->NEDOS)), sp_Idos_sum_dw(16, vector<double>(this->NEDOS));
			double center_up[4], center_dw[4], sp_center_up[16], sp_center_dw[16];
			for (int i = 0; i < num.size(); i++)
			{
				dos_sum_up[0] = this->orbit_dos[28][num[i] - 1] + dos_sum_up[0];
				Idos_sum_up[0] = this->Idos[28][num[i] - 1] + Idos_sum_up[0];
				sp_dos_sum_up[0] = this->orbit_dos[28][num[i] - 1] + sp_dos_sum_up[0];
				sp_Idos_sum_up[0] = this->Idos[28][num[i] - 1] + sp_Idos_sum_up[0];
				for (int j = 29, k = 1; j < 32; j++, k++)
				{
					dos_sum_up[1] = this->orbit_dos[j][num[i] - 1] + dos_sum_up[1];
					Idos_sum_up[1] = this->Idos[j][num[i] - 1] + Idos_sum_up[1];
					sp_dos_sum_up[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_up[k];
					sp_Idos_sum_up[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_up[k];
				}
				for (int j = 32, k = 4; j < 37; j++, k++)
				{
					dos_sum_up[2] = this->orbit_dos[j][num[i] - 1] + dos_sum_up[2];
					Idos_sum_up[2] = this->Idos[j][num[i] - 1] + Idos_sum_up[2];
					sp_dos_sum_up[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_up[k];
					sp_Idos_sum_up[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_up[k];
				}
				for (int j = 37, k = 9; j < 44; j++, k++)
				{
					dos_sum_up[3] = this->orbit_dos[j][num[i] - 1] + dos_sum_up[3];
					Idos_sum_up[3] = this->Idos[j][num[i] - 1] + Idos_sum_up[3];
					sp_dos_sum_up[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_up[k];
					sp_Idos_sum_up[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_up[k];
				}
			}
			for (int i = 0; i < num.size(); i++)
			{
				dos_sum_dw[0] = this->orbit_dos[28][num[i] - 1] + dos_sum_dw[0];
				Idos_sum_dw[0] = this->Idos[28][num[i] - 1] + Idos_sum_dw[0];
				sp_dos_sum_dw[0] = this->orbit_dos[28][num[i] - 1] + sp_dos_sum_dw[0];
				sp_Idos_sum_dw[0] = this->Idos[28][num[i] - 1] + sp_Idos_sum_dw[0];
				for (int j = 45, k = 1; j < 48; j++, k++)
				{
					dos_sum_dw[1] = this->orbit_dos[j][num[i] - 1] + dos_sum_dw[1];
					Idos_sum_dw[1] = this->Idos[j][num[i] - 1] + Idos_sum_dw[1];
					sp_dos_sum_dw[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_dw[k];
					sp_Idos_sum_dw[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_dw[k];
				}
				for (int j = 48, k = 4; j < 53; j++, k++)
				{
					dos_sum_dw[2] = this->orbit_dos[j][num[i] - 1] + dos_sum_dw[2];
					Idos_sum_dw[2] = this->Idos[j][num[i] - 1] + Idos_sum_dw[2];
					sp_dos_sum_dw[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_dw[k];
					sp_Idos_sum_dw[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_dw[k];
				}
				for (int j = 53, k = 9; j < 60; j++, k++)
				{
					dos_sum_dw[3] = this->orbit_dos[j][num[i] - 1] + dos_sum_dw[3];
					Idos_sum_dw[3] = this->Idos[j][num[i] - 1] + Idos_sum_dw[3];
					sp_dos_sum_dw[k] = this->orbit_dos[j][num[i] - 1] + sp_dos_sum_dw[k];
					sp_Idos_sum_dw[k] = this->Idos[j][num[i] - 1] + sp_Idos_sum_dw[k];
				}
			}
			for (int i = 0; i < 4; i++)
			{
				vector<double> temp_up = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_up[i]);
				center_up[i] = temp_up[this->NEDOS - 1] / Idos_sum_up[i][this->NEDOS - 1];
				vector<double> temp_dw = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * dos_sum_dw[i]);
				center_dw[i] = temp_dw[this->NEDOS - 1] / Idos_sum_dw[i][this->NEDOS - 1];
			}
			fprintf(fp, "#	Atom_ID   s-band-center    p-band-center    d-band-center  f-band-center(in units of eV)\n");
			fprintf(fp, "SPIN-UP	%s	%lf	%lf	%lf	%lf\n", atomid.c_str(), center_up[0], center_up[1], center_up[2], center_up[3]);
			fprintf(fp, "SPIN-DW	%s	%lf	%lf	%lf	%lf\n\n", atomid.c_str(), center_dw[0], center_dw[1], center_dw[2], center_dw[3]);
			fprintf(fp, "#	Atom_ID   s   py   pz   px   dxy   dyz   dz2   dxz   dx2   (in units of eV)\n");
			fprintf(fp, "SPIN-UP	%s", atomid.c_str());
			for (int i = 0; i < 16; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * sp_dos_sum_up[i]);
				sp_center_up[i] = temp[this->NEDOS - 1] / sp_Idos_sum_up[i][this->NEDOS - 1];
				fprintf(fp, "	%lf", sp_center_up[i]);
			}
			fprintf(fp, "\n");
			fprintf(fp, "SPIN-DW	%s", atomid.c_str());
			for (int i = 0; i < 16; i++)
			{
				vector<double> temp = Intergral(this->energy - this->efermi, (this->energy - this->efermi) * sp_dos_sum_dw[i]);
				sp_center_dw[i] = temp[this->NEDOS - 1] / sp_Idos_sum_dw[i][this->NEDOS - 1];
				fprintf(fp, "	%lf", sp_center_dw[i]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n\n");
	fprintf(fp, "Remark:\n");
	fprintf(fp, "Energy window for integration is from %lf to %lf.\n", this->energy[0] - this->efermi, this->energy[NEDOS - 1] - this->efermi);
	fprintf(fp, "Band Center is respect to Fermi level, i.e., E_F = 0 eV.\n");
	fclose(fp);
	printf("Written %s file!\n", "Band_Center");
}

typedef struct INPUT
{
	vector<int> num;
	int head_index;
	string elemlabel;
}INPUT;

void DOS::GetDos(int argc, char* argv[])
{
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
	FILE* fp = fopen("POSCAR", "r");
	if (fp == NULL)
	{
		perror("POSCAR IS NOT EXIST!\n");
		exit(1);
	}
	readposcar(fp, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	int max_element = 0;
	for (int i = 0; i < nant[1]; i++)
		max_element = max(max_element, elemnum[i]);
	if (strcmp(argv[1], "--dos"))
		return;
	DOS dos;
	string ISPIN = GetInfoINCAR("ISPIN");
	if (ISPIN.length() == 0)
		ISPIN = "1";
	string LORBIT = GetInfoINCAR("LORBIT");
	if (LORBIT.length() == 0)
		LORBIT = "0";
	dos.readDoscar(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	if (!strcmp(argv[2], "-efermi"))
		printf("%lf\n", dos.efermi);
	if (!strcmp(argv[2], "-t") || !strcmp(argv[2], "-T")) //total VASPMATE --dos -t/-T
		dos.TDos();
	if (!strcmp(argv[2], "-a") || !strcmp(argv[2], "-A")) //atom VASPMATE --dos -a/-A num(1,2,3...)
	{
		for (int i = 0; i < nant[0]; i++)
			dos.ADos(i + 1, NULL, 0, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element);
	}
	if (!strcmp(argv[2], "-e") || !strcmp(argv[2], "-E")) //element VASPMATE --dos -e/-E 
	{
		vector<string> elem;
		elem.resize(nant[1]);
		for (int i = 0; i < nant[1]; i++)
			elem.push_back(elemsym[i]);
		dos.EDos(elem, typenum, elemsym, nant, max_element, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	}
	if (!strcmp(argv[2], "-s") || !strcmp(argv[2], "-S"))
	{
		vector<int> num;
		vector<INPUT> data;
		for (int i = 3; i < argc; i++)
		{
			INPUT input;
			if (!strcmp(argv[i], "all"))
			{
				for (int i = 0; i < nant[0]; i++)
					num.push_back(i + 1);
			}
			else if (strstr(argv[i], "-") != NULL) // 
			{
				int start, end;
				sscanf(argv[i], "%d%*c%d", &start, &end);
				if (start > end || start<1 || end>nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				for (int i = start; i < end + 1; i++)
					num.push_back(i);
			}
			else if (IsNum(argv[i]))
			{
				if (atoi(argv[i]) <= 0 || atoi(argv[i]) > nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				num.push_back(atoi(argv[i]));
			}
			else if (!IsNum(argv[i]))
			{
				int IsRightelem = 0;
				for (int j = 0; j < nant[1]; j++)
				{
					if (!strcmp(argv[i], elemsym[j]))
						IsRightelem = 1;
				}
				if (IsRightelem == 0)
				{
					printf("%s is not in element types.\n", argv[i]);
					continue;
				}
				else
				{
					input.elemlabel = argv[i];
					int cnt = 0;
					for (int k = 0; k < nant[1]; k++)
					{
						if (!strcmp(argv[i], elemsym[k]))
						{
							input.head_index = cnt;
							for (int j = 0; j < typenum[k]; j++)
								input.num.push_back(j + cnt + 1);
							data.push_back(input);
						}
						cnt += typenum[k];
					}
				}
			}
		}
		sort(num.begin(), num.end());
		num.erase(unique(num.begin(), num.end()), num.end());
		for (int i = 0; i < num.size(); i++)
			dos.ADos(num[i], NULL, 0, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element);
		for (int i = 0; i < data.size(); i++)
			for (int j = 0; j < data[i].num.size(); j++)
				dos.ADos(data[i].num[j], data[i].elemlabel.c_str(), data[i].head_index, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element);
	}
	if (!strcmp(argv[2], "-sa") || !strcmp(argv[2], "-SA"))	//select atom VASPMATE --dos --sa/-SA (a b c...)/(a-b)
	{
		vector<int> num;
		for (int i = 3; i < argc; i++)
		{
			if (!strcmp(argv[i], "all"))
			{
				for (int i = 0; i < nant[0]; i++)
					num.push_back(i + 1);
				break;
			}
			else if (strstr(argv[i], "-") != NULL) // 
			{
				int start, end;
				sscanf(argv[i], "%d%*c%d", &start, &end);
				if (start > end || start<1 || end>nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				for (int i = start; i < end + 1; i++)
					num.push_back(i);
			}
			else if (IsNum(argv[i]))
			{
				if (atoi(argv[i]) <= 0 || atoi(argv[i]) > nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				num.push_back(atoi(argv[i]));
			}
			else if (!IsNum(argv[i]))
			{
				printf("Please Input Atom Index!\n");
				return;
			}
		}
		//erase same value
		sort(num.begin(), num.end());
		num.erase(unique(num.begin(), num.end()), num.end());
		for (int i = 0; i < num.size(); i++)
			dos.ADos(num[i], NULL, 0, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element);
	}
	if (!strcmp(argv[2], "-se") || !strcmp(argv[2], "-SE") || !strcmp(argv[2], "-select"))	//select atom VASPMATE --dos --sa/-SA (a b c...)/(a-b)
	{
		vector<INPUT> data;
		for (int i = 3; i < argc; i++)
		{
			INPUT input;
			if (IsNum(argv[i]))
			{
				printf("Please Input Element Symbol!\n");
				return;
			}
			else if (!IsNum(argv[i]))
			{
				int IsRightelem = 0;
				for (int j = 0; j < nant[1]; j++)
				{
					if (!strcmp(argv[i], elemsym[j]))
						IsRightelem = 1;
				}
				if (IsRightelem == 0)
				{
					printf("%s is not in element types.\n", argv[i]);
					continue;
				}
				else
				{
					input.elemlabel = argv[i];
					int cnt = 0;
					for (int k = 0; k < nant[1]; k++)
					{
						if (!strcmp(argv[i], elemsym[k]))
						{
							input.head_index = cnt;
							for (int j = 0; j < typenum[k]; j++)
								input.num.push_back(j + cnt + 1);
							data.push_back(input);
						}
						cnt += typenum[k];
					}
				}
			}
		}
		for (int i = 0; i < data.size(); i++)
			for (int j = 0; j < data[i].num.size(); j++)
				dos.ADos(data[i].num[j], data[i].elemlabel.c_str(), data[i].head_index, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element);
	}
	if (!strcmp(argv[2], "-m") || !strcmp(argv[2], "-M") || !strcmp(argv[2], "-ma") || !strcmp(argv[2], "-me"))	//Multiple Atoms or Elements
	{
		vector<int> num = TranArgvToAtomIndex(argc, argv, nant, typenum, elemsym, xyz, const_cast<char*>("w"));
		string postfix;
		for (int i = 3; i < argc; i++)
		{
			if (i != argc - 1)
				postfix = postfix + argv[i] + "_";
			else
				postfix += argv[i];
		}
		dos.SDos(num, postfix, max_element, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	}
	if (!strcmp(argv[2], "-o") || !strcmp(argv[2], "-O") || !strcmp(argv[2], "-oa") || !strcmp(argv[2], "-oe"))	//orbit VASPMATE --dos -o element(index) s p d ...
	{
		vector<COMB> orbit;
		vector<vector<string> > ind(argc / 2);
		vector<vector<string> > orb(argc / 2);
		COMB comb;
		int cnt = 0;
		for (int i = 3; i < argc; i++)
		{
			//save atom index
			// 1 2 3 
			if (IsNum(argv[i]))
			{
				ind[cnt].push_back(argv[i]);
				//check the index is right
				if (atoi(argv[i]) < 0 || atoi(argv[i]) > nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				comb.index.push_back(atoi(argv[i]));
			}
			//1-3
			else if (strstr(argv[i], "-") != NULL)
			{
				ind[cnt].push_back(argv[i]);
				int start, end;
				sscanf(argv[i], "%d%*c%d", &start, &end);
				if (start > end || start<1 || end>nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				for (int i = start; i < end + 1; i++)
					comb.index.push_back(i);
			}
			//B N
			else if (strcmp("all", argv[i]) && !IsRightOrbit(argv[i], atoi(LORBIT.c_str()), max_element))
			{
				ind[cnt].push_back(argv[i]);
				//check symbol is right element
				int IsRightelem = 0;
				for (int j = 0; j < nant[1]; j++)
				{
					if (!strcmp(argv[i], elemsym[j]) || argv[i] == "all")
						IsRightelem = 1;
				}
				if (IsRightelem == 0)
				{
					printf("%s is not in element types.\n", argv[i]);
					return;
				}
				//tran element to atom index
				int cnt = 0;
				for (int k = 0; k < nant[1]; k++)
				{
					if (!strcmp(argv[i], elemsym[k]))
					{
						for (int j = 0; j < typenum[k]; j++)
							comb.index.push_back(j + cnt + 1);
					}
					cnt += typenum[k];
				}
			}
			//save orbit index
			else if (!strcmp("all", argv[i]))
			{
				orb[cnt].push_back(argv[i]);
				comb.orb.resize(20);
				if (atoi(LORBIT.c_str()) == 10 && max_element <= 57)
				{
					comb.orb = { "s","p","d" };
				}
				if (atoi(LORBIT.c_str()) == 10 && max_element > 57)
				{
					comb.orb = { "s","p","d","f" };
				}
				if (atoi(LORBIT.c_str()) == 11 && max_element <= 57)
				{
					comb.orb = { "s","px","py","pz","dxy","dxz","dyz","dx2","dz2" };
				}
				if (atoi(LORBIT.c_str()) == 11 && max_element > 57)
				{
					comb.orb = { "s","px","py","pz","dxy","dxz","dyz","dx2","dz2","f1","f2","f3","f4","f5","f6","f7" };
				}
			}
			else
			{
				if (atoi(LORBIT.c_str()) == 11 && !strcmp(argv[i], "p"))
				{
					orb[cnt].push_back("p");
					comb.orb.push_back("py");
					comb.orb.push_back("pz");
					comb.orb.push_back("px");
				}
				else if (atoi(LORBIT.c_str()) == 11 && !strcmp(argv[i], "d"))
				{
					orb[cnt].push_back("d");
					comb.orb.push_back("dxy");
					comb.orb.push_back("dyz");
					comb.orb.push_back("dz2");
					comb.orb.push_back("dxz");
					comb.orb.push_back("dx2");
				}
				else if (atoi(LORBIT.c_str()) == 11 && !strcmp(argv[i], "f"))
				{
					orb[cnt].push_back("f");
					comb.orb.push_back("f1");
					comb.orb.push_back("f2");
					comb.orb.push_back("f3");
					comb.orb.push_back("f4");
					comb.orb.push_back("f5");
					comb.orb.push_back("f6");
					comb.orb.push_back("f7");
				}
				else
				{
					orb[cnt].push_back(argv[i]);
					comb.orb.push_back(argv[i]);
				}
			}
			if (i == argc - 1 || (i < argc - 1 && (IsRightOrbit(argv[i], atoi(LORBIT.c_str()), max_element) || !strcmp("all", argv[i])) &&
				!IsRightOrbit(argv[i + 1], atoi(LORBIT.c_str()), max_element)))
			{
				orbit.push_back(comb);
				comb.index.clear();
				comb.orb.clear();
				cnt++;
			}
		}
		dos.ODos(ind, orb, orbit, typenum, elemsym, nant, max_element, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	}
	if (!strcmp(argv[2], "-bc") || !strcmp(argv[2], "-BC")) //BandCenter VASPMATE --dos -bc element(atom_index,all)
	{
		string atomid;
		for (int i = 3; i < argc; i++)
		{
			if (i != argc - 1)
				atomid = atomid + argv[i] + "&";
			else
				atomid += argv[i];
		}
		vector<int> num = TranArgvToAtomIndex(argc, argv, nant, typenum, elemsym, xyz, NULL);
		dos.BandCenter(num, atomid, max_element, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	}
	fclose(fp);
}
