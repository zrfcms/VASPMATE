#include"../include/band.h"

void BAND::readPROCAR(int ISPIN, int LORBIT)
{
	if (ISPIN != 1 && ISPIN != 2)
	{
		printf("ISPIN IS ERROR IN INCAR!\n");
		return;
	}
	if (LORBIT != 10 && LORBIT != 11)
	{
		printf("LORBIT IS ERROR IN INCAR!\n");
		return;
	}
	FILE* fp = fopen("PROCAR", "r");
	if (fp == NULL)
	{
		printf("PROCAR IS NOT EXIST!\n");
		return;
	}
	char buf[1024];
	int line = 0;
	int kpt = 0;
	int bd = 0;
	int ion = 0;
	int flag = 0;
	while (fgets(buf, 1024, fp) != NULL)
	{
		line++;
		if (ISPIN == 1)
		{
			if (line == 1)
			{
				this->title = buf;
				continue;
			}
			if (line == 2)
			{
				sscanf(buf, "%*[^:]%*[:]%d%*[^:]%*[:]%d%*[^:]%*[:]%d", &this->num_kpts, &this->num_bands, &this->num_ions);
				// initial size of vector
				this->k_coor.resize(this->num_kpts);
				for (int i = 0; i < this->k_coor.size(); i++) this->k_coor[i].resize(3);
				this->weight.resize(this->num_kpts);
				this->energy.resize(this->num_kpts);
				for (int i = 0; i < this->energy.size(); i++) this->energy[i].resize(this->num_bands);
				this->occ.resize(this->num_kpts);
				for (int i = 0; i < this->energy.size(); i++) this->occ[i].resize(this->num_bands);
				this->ion_dos.resize(this->num_kpts);
				for (int i = 0; i < this->ion_dos.size(); i++)
				{
					this->ion_dos[i].resize(this->num_bands);
					for (int j = 0; j < this->ion_dos[i].size(); j++)
					{
						this->ion_dos[i][j].resize(this->num_ions + 1);
						for (int k = 0; k < this->ion_dos[i][j].size(); k++)
						{
							if (LORBIT == 10)
								this->ion_dos[i][j][k].resize(4);
							else if (LORBIT == 11)
								this->ion_dos[i][j][k].resize(16);//1+3+5+7
						}
					}
				}
				continue;
			}
			if (strstr(buf, "k-point") != NULL)
			{
				sscanf(buf, "%*[^:]%*[:]%lf%lf%lf%*[^=]%*[=]%lf", &this->k_coor[kpt][0], &this->k_coor[kpt][1], &this->k_coor[kpt][2], &this->weight[kpt]);
				bd = 0;
				kpt++;
				continue;
			}
			if (strstr(buf, "band") != NULL)
			{
				sscanf(buf, "%*[^#]%*[#]%*s%lf%*s%*s%lf", &this->energy[kpt - 1][bd], &this->occ[kpt - 1][bd]);
				ion = 0;
				bd++;
				continue;
			}
			if (strstr(buf, "ion") != NULL)
			{
				flag = 1;
				continue;
			}
			if (flag == 1)
			{
				if (LORBIT == 10)
					sscanf(buf, "%*s%lf%lf%lf%lf", &this->ion_dos[kpt - 1][bd - 1][ion][0], &this->ion_dos[kpt - 1][bd - 1][ion][1],
						&this->ion_dos[kpt - 1][bd - 1][ion][2], &this->ion_dos[kpt - 1][bd - 1][ion][3]);
				else if (LORBIT == 11)
					sscanf(buf, "%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &this->ion_dos[kpt - 1][bd - 1][ion][0],
						&this->ion_dos[kpt - 1][bd - 1][ion][1], &this->ion_dos[kpt - 1][bd - 1][ion][2], &this->ion_dos[kpt - 1][bd - 1][ion][3],
						&this->ion_dos[kpt - 1][bd - 1][ion][4], &this->ion_dos[kpt - 1][bd - 1][ion][5], &this->ion_dos[kpt - 1][bd - 1][ion][6],
						&this->ion_dos[kpt - 1][bd - 1][ion][7], &this->ion_dos[kpt - 1][bd - 1][ion][8], &this->ion_dos[kpt - 1][bd - 1][ion][9],
						&this->ion_dos[kpt - 1][bd - 1][ion][10], &this->ion_dos[kpt - 1][bd - 1][ion][11], &this->ion_dos[kpt - 1][bd - 1][ion][12],
						&this->ion_dos[kpt - 1][bd - 1][ion][13], &this->ion_dos[kpt - 1][bd - 1][ion][14], &this->ion_dos[kpt - 1][bd - 1][ion][15]);
				ion++;
				if ((ion == this->num_ions + 1 && this->num_ions != 1) || (ion == this->num_ions && this->num_ions == 1))
					flag = 0;
			}
		}
		if (ISPIN == 2)
		{
			if (line == 1)
			{
				this->title = buf;
				continue;
			}
			if (buf[0] == '#')
			{
				if (kpt == 0)
				{
					sscanf(buf, "%*[^:]%*[:]%d%*[^:]%*[:]%d%*[^:]%*[:]%d", &this->num_kpts, &this->num_bands, &this->num_ions);
					// initial size of vector
					this->k_coor.resize(this->num_kpts);
					for (int i = 0; i < this->k_coor.size(); i++) this->k_coor[i].resize(3);
					this->weight.resize(this->num_kpts);
					this->energy_up.resize(this->num_kpts);
					this->energy_dw.resize(this->num_kpts);
					for (int i = 0; i < this->energy_up.size(); i++) this->energy_up[i].resize(this->num_bands);
					for (int i = 0; i < this->energy_dw.size(); i++) this->energy_dw[i].resize(this->num_bands);
					this->occ_up.resize(this->num_kpts);
					this->occ_dw.resize(this->num_kpts);
					for (int i = 0; i < this->occ_up.size(); i++) this->occ_up[i].resize(this->num_bands);
					for (int i = 0; i < this->occ_dw.size(); i++) this->occ_dw[i].resize(this->num_bands);
					this->ion_dos_up.resize(this->num_kpts);
					this->ion_dos_dw.resize(this->num_kpts);
					for (int i = 0; i < this->ion_dos_up.size(); i++)
					{
						this->ion_dos_up[i].resize(this->num_bands);
						this->ion_dos_dw[i].resize(this->num_bands);
						for (int j = 0; j < this->ion_dos_up[i].size(); j++)
						{
							this->ion_dos_up[i][j].resize(this->num_ions + 1);
							this->ion_dos_dw[i][j].resize(this->num_ions + 1);
							for (int k = 0; k < this->ion_dos_up[i][j].size(); k++)
							{
								if (LORBIT == 10)
								{
									this->ion_dos_up[i][j][k].resize(4);
									this->ion_dos_dw[i][j][k].resize(4);
								}
								else if (LORBIT == 11)
								{
									this->ion_dos_up[i][j][k].resize(16);//1+3+5+7
									this->ion_dos_dw[i][j][k].resize(16);
								}
							}
						}
					}
					continue;
				}
				else
					continue;
			}
			if (strstr(buf, "k-point") != NULL)
			{
				if (kpt < this->num_kpts)
					sscanf(buf, "%*[^:]%*[:]%lf%lf%lf%*[^=]%*[=]%lf", &this->k_coor[kpt][0], &this->k_coor[kpt][1],
						&this->k_coor[kpt][2], &this->weight[kpt]);
				bd = 0;
				kpt++;
				continue;
			}
			if (strstr(buf, "band") != NULL)
			{
				if (kpt - 1 < this->num_kpts)
					sscanf(buf, "%*[^#]%*[#]%*s%lf%*s%*s%lf", &this->energy_up[kpt - 1][bd], &this->occ_up[kpt - 1][bd]);
				else
					sscanf(buf, "%*[^#]%*[#]%*s%lf%*s%*s%lf", &this->energy_dw[kpt - 1 - this->num_kpts][bd], &this->occ_dw[kpt - 1 - this->num_kpts][bd]);
				ion = 0;
				bd++;
				continue;
			}
			if (strstr(buf, "ion") != NULL)
			{
				flag = 1;
				continue;
			}
			if (flag == 1)
			{
				if (LORBIT == 10)
				{
					if (kpt - 1 < this->num_kpts)
						sscanf(buf, "%*s%lf%lf%lf%lf", &this->ion_dos_up[kpt - 1][bd - 1][ion][0], &this->ion_dos_up[kpt - 1][bd - 1][ion][1],
							&this->ion_dos_up[kpt - 1][bd - 1][ion][2], &this->ion_dos_up[kpt - 1][bd - 1][ion][3]);
					else
						sscanf(buf, "%*s%lf%lf%lf%lf", &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][0], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][1],
							&this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][2], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][3]);
				}
				else if (LORBIT == 11)
				{
					if (kpt - 1 < this->num_kpts)
						sscanf(buf, "%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &this->ion_dos_up[kpt - 1][bd - 1][ion][0],
							&this->ion_dos_up[kpt - 1][bd - 1][ion][1], &this->ion_dos_up[kpt - 1][bd - 1][ion][2], &this->ion_dos_up[kpt - 1][bd - 1][ion][3],
							&this->ion_dos_up[kpt - 1][bd - 1][ion][4], &this->ion_dos_up[kpt - 1][bd - 1][ion][5], &this->ion_dos_up[kpt - 1][bd - 1][ion][6],
							&this->ion_dos_up[kpt - 1][bd - 1][ion][7], &this->ion_dos_up[kpt - 1][bd - 1][ion][8], &this->ion_dos_up[kpt - 1][bd - 1][ion][9],
							&this->ion_dos_up[kpt - 1][bd - 1][ion][10], &this->ion_dos_up[kpt - 1][bd - 1][ion][11], &this->ion_dos_up[kpt - 1][bd - 1][ion][12],
							&this->ion_dos_up[kpt - 1][bd - 1][ion][13], &this->ion_dos_up[kpt - 1][bd - 1][ion][14], &this->ion_dos_up[kpt - 1][bd - 1][ion][15]);
					else
						sscanf(buf, "%*s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][0],
							&this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][1], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][2], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][3],
							&this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][4], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][5], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][6],
							&this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][7], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][8], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][9],
							&this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][10], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][11], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][12],
							&this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][13], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][14], &this->ion_dos_dw[kpt - 1 - this->num_kpts][bd - 1][ion][15]);
				}
				ion++;
				if ((ion == this->num_ions + 1 && this->num_ions != 1) || (ion == this->num_ions && this->num_ions == 1))
					flag = 0;
			}
		}
	}
}

vector<double> BAND::getkpointsline(double vec[3][3])
{
	this->start_point = 0;
	vector<double> kline;
	FILE* fp = fopen("NEWKPATH", "r");
	if (fp == NULL)
	{
		printf("NEWKPATH IS NOT EXIST!\n");
		return kline;
	}
	char buf[1024];
	int line = 0;
	int point = 0;
	vector<string> label;
	vector<vector<double> > high_symmetry;
	while (fgets(buf, 1024, fp) != NULL)
	{
		line++;
		if (line == 2)
		{
			sscanf(buf, "%d", &point);
			continue;
		}
		if (strspn(buf, " \t\n") == strlen(buf))
			continue;
		if (line > 4)
		{
			char la[] = "XX";
			vector<double> hs(3, 0);
			sscanf(buf, "%lf%lf%lf%s", &hs[0], &hs[1], &hs[2], la);
			high_symmetry.push_back(hs);
			label.push_back(la);
		}
	}
	fclose(fp);
	double rec_vec[3][3];
	RecMat(vec, rec_vec);
	vector<vector<double> > carts_kpt(this->num_kpts, vector<double>(3));
	for (int i = 0; i < this->k_coor.size(); i++)
	{
		double p[3] = { k_coor[i][0],k_coor[i][1],k_coor[i][2] };
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
		double d = (i + 1) % point == 0 ? 0 : sqrt((carts_kpt[i][0] - carts_kpt[i + 1][0]) * (carts_kpt[i][0] - carts_kpt[i + 1][0]) +
			(carts_kpt[i][1] - carts_kpt[i + 1][1]) * (carts_kpt[i][1] - carts_kpt[i + 1][1]) +
			(carts_kpt[i][2] - carts_kpt[i + 1][2]) * (carts_kpt[i][2] - carts_kpt[i + 1][2]));
		delta.push_back(d);
	}
	for (int i = 1; i < kline.size(); i++)
		kline[i] = kline[i - 1] + delta[i - 1];
	FILE* fp1 = fopen("KLABELS", "w");
	fprintf(fp1, "K-Label    Coordinate of high-symmetry k-point in band-structure plots\n");
	fprintf(fp1, "	%s	%lf\n", label[0].c_str(), kline[0]);
	for (int i = 1; i < this->num_kpts / point; i++)
	{
		if (high_symmetry[2 * i - 1][0] == high_symmetry[2 * i][0] && high_symmetry[2 * i - 1][1] == high_symmetry[2 * i][1] && high_symmetry[2 * i - 1][2] == high_symmetry[2 * i][2])
			fprintf(fp1, "	%s	%lf\n", label[2 * i].c_str(), kline[i * point]);
		else
			fprintf(fp1, "	%s	%lf\n", (label[2 * i - 1] + '|' + label[2 * i]).c_str(), kline[i * point]);
	}
	fprintf(fp1, "	%s	%lf\n", label[label.size() - 1].c_str(), kline[kline.size() - 1]);
	fprintf(fp1, "\n\n");
	fprintf(fp1, "*Tips: Tip: Please label each high symmetry point in NEWKPATH file. Otherwise, they will be identified as 'XX' in KLABELS file\n");
	fclose(fp1);
	printf("Written KLABELS file!\n");
	return kline;
}

vector<double> BAND::hse_getkpointsline(int nkpts, vector<vector<double> > kcoor, double vec[3][3])
{
	vector<double> kline;
	FILE* fp = fopen("NEWKPATH", "r");
	if (fp == NULL)
	{
		printf("NEWKPATH IS NOT EXIST!\n");
		return kline;
	}
	char buf[1024];
	int line = 0;
	vector<string> label;
	vector<vector<double> > high_symmetry;
	while (fgets(buf, 1024, fp) != NULL)
	{
		line++;
		if (strspn(buf, " \t\n") == strlen(buf))
			continue;
		if (line > 4)
		{
			char la[] = "XX";
			vector<double> hs(3, 0);
			sscanf(buf, "%lf%lf%lf%s", &hs[0], &hs[1], &hs[2], la);
			high_symmetry.push_back(hs);
			label.push_back(la);
		}
	}
	fclose(fp);
	FILE* _fp = fopen("KPOINTS", "r");
	char kbuf[1024];
	fgets(kbuf, 1024, _fp);
	fclose(_fp);
	vector<int> temp = ExtractNumbersFromString(kbuf, sizeof(kbuf));
	this->start_point = temp[3];
	vector<int> point(temp.begin() + 6, temp.begin() + 6 + temp[5]);
	double rec_vec[3][3];
	RecMat(vec, rec_vec);
	vector<vector<double> > carts_kpt(nkpts - this->start_point, vector<double>(3));
	for (int i = this->start_point; i < kcoor.size(); i++)
	{
		double p[3] = { kcoor[i][0],kcoor[i][1],kcoor[i][2] };
		for (int j = 0; j < 3; j++)
		{
			carts_kpt[i - this->start_point][j] = p[0] * rec_vec[0][j] + p[1] * rec_vec[1][j] + p[2] * rec_vec[2][j];
		}
	}
	kline.resize(nkpts - this->start_point);
	kline[0] = 0;
	vector<double> delta;
	int cnt1 = 0;
	for (int i = 0; i < point.size(); i++)
	{
		for (int j = 0; j < point[i]; j++)
		{
			if (cnt1 + j == nkpts - this->start_point)
				break;
			double d = (j + 1) % point[i] == 0 ? 0 : sqrt((carts_kpt[cnt1 + j][0] - carts_kpt[cnt1 + j + 1][0]) * (carts_kpt[cnt1 + j][0] - carts_kpt[cnt1 + j + 1][0]) +
				(carts_kpt[cnt1 + j][1] - carts_kpt[cnt1 + j + 1][1]) * (carts_kpt[cnt1 + j][1] - carts_kpt[cnt1 + j + 1][1]) +
				(carts_kpt[cnt1 + j][2] - carts_kpt[cnt1 + j + 1][2]) * (carts_kpt[cnt1 + j][2] - carts_kpt[cnt1 + j + 1][2]));
			delta.push_back(d);
		}
		cnt1 += point[i];
	}
	for (int i = 1; i < kline.size(); i++)
		kline[i] = kline[i - 1] + delta[i - 1];
	FILE* fp1 = fopen("KLABELS", "w");
	fprintf(fp1, "K-Label    Coordinate of high-symmetry k-point in band-structure plots\n");
	fprintf(fp1, "	%s	%lf\n", label[0].c_str(), kline[0]);
	int cnt2 = 0;
	for (int i = 1; i < point.size(); i++)
	{
		cnt2 += point[i - 1];
		if (high_symmetry[2 * i - 1][0] == high_symmetry[2 * i][0] && high_symmetry[2 * i - 1][1] == high_symmetry[2 * i][1] && high_symmetry[2 * i - 1][2] == high_symmetry[2 * i][2])
			fprintf(fp1, "	%s	%lf\n", label[2 * i].c_str(), kline[cnt2]);
		else
			fprintf(fp1, "	%s	%lf\n", (label[2 * i - 1] + '|' + label[2 * i]).c_str(), kline[cnt2]);
	}
	fprintf(fp1, "	%s	%lf\n", label[label.size() - 1].c_str(), kline[kline.size() - 1]);
	fprintf(fp1, "\n\n");
	fprintf(fp1, "*Tips: Tip: Please label each high symmetry point in NEWKPATH file. Otherwise, they will be identified as 'XX' in KLABELS file\n");
	fclose(fp1);
	printf("Written KLABELS file!\n");
	return kline;
}

void BAND::basicband(int ISPIN, vector<double> kline, double efermi)
{
	if (ISPIN == 1)
	{
		FILE* fp1 = fopen("BAND.dat", "w");
		//label
		FILE* fpk = fopen("KLABELS", "r");
		int line = 0;
		char buf[1024];
		vector<double> dis;
		while (fgets(buf, 1024, fpk) != NULL)
		{
			line++;
			if (line == 1)
				continue;
			if (strspn(buf, " \t\n") == strlen(buf))
				break;
			string s;
			double d;
			sscanf(buf, "%s%lf", s.c_str(), &d);
			dis.push_back(d);
		}
		fclose(fpk);
		// normal
		fprintf(fp1, "#K-Path(1/A) Energy-Level(eV)\n");
		fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
		for (int i = 0; i < this->num_bands; i++)
		{
			fprintf(fp1, "# Band-Index   %d\n", i + 1);
			if (i % 2 == 0)
			{
				for (int j = start_point; j < this->num_kpts; j++)
				{
					fprintf(fp1, " %lf	%lf\n", kline[j - this->start_point], this->energy[j][i] - efermi);
				}
			}
			else
			{
				for (int j = this->num_kpts - 1; j >= start_point; j--)
				{
					fprintf(fp1, " %lf	%lf\n", kline[j - this->start_point], this->energy[j][i] - efermi);
				}
			}
			fprintf(fp1, "\n");
		}
		fclose(fp1);
		printf("Written BAND.dat file!\n");
		FILE* fp2 = fopen("REFORMATTED_BAND.dat", "w");
		fprintf(fp2, "#K-Path     Band-1    Band-2    Band-3    Band-4      ...     Band-n\n");
		for (int i = 0; i < this->num_kpts; i++)
		{
			fprintf(fp2, "  %lf  ", kline[i - this->start_point]);
			for (int j = 0; j < this->num_bands; j++)
				fprintf(fp2, "%lf  ", this->energy[i][j] - efermi);
			fprintf(fp2, "\n");
		}
		fclose(fp2);
		printf("Written REFORMATTED_BAND.dat file!\n");
	}
	else if (ISPIN == 2)
	{
		FILE* fp1 = fopen("BAND.dat", "w");
		fprintf(fp1, "#K-Path(1/A)         Spin-up(eV)   Spin-down(eV)\n");
		fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
		//label
		FILE* fpk = fopen("KLABELS", "r");
		int line = 0;
		char buf[1024];
		vector<double> dis;
		while (fgets(buf, 1024, fpk) != NULL)
		{
			line++;
			if (line == 1)
				continue;
			if (strspn(buf, " \t\n") == strlen(buf))
				break;
			string s;
			double d;
			sscanf(buf, "%s%lf", s.c_str(), &d);
			dis.push_back(d);
		}
		fclose(fpk);
		for (int i = 0; i < this->num_bands; i++)
		{
			fprintf(fp1, "# Band-Index   %d\n", i + 1);
			if (i % 2 == 0)
			{
				for (int j = start_point; j < this->num_kpts; j++)
				{
					fprintf(fp1, " %lf	%lf	%lf\n", kline[j - this->start_point], this->energy_up[j][i] - efermi, this->energy_dw[j][i] - efermi);
				}
			}
			else
			{
				for (int j = this->num_kpts - 1; j >= start_point; j--)
				{
					fprintf(fp1, " %lf	%lf	%lf\n", kline[j - this->start_point], this->energy_up[j][i] - efermi, this->energy_dw[j][i] - efermi);
				}
			}
			fprintf(fp1, "\n");
		}
		fclose(fp1);
		printf("Written BAND.dat file!\n");
		FILE* fp2_up = fopen("REFORMATTED_BAND_UP.dat", "w");
		fprintf(fp2_up, "#K-Path     Band-1    Band-2    Band-3    Band-4      ...     Band-n\n");
		for (int i = 0; i < this->num_kpts; i++)
		{
			fprintf(fp2_up, "  %lf  ", kline[i - this->start_point]);
			for (int j = 0; j < this->num_bands; j++)
				fprintf(fp2_up, "%lf  ", this->energy_up[i][j] - efermi);
			fprintf(fp2_up, "\n");
		}
		fclose(fp2_up);
		printf("Written REFORMATTED_BAND_UP.dat file!\n");
		FILE* fp2_dw = fopen("REFORMATTED_BAND_DW.dat", "w");
		fprintf(fp2_dw, "#K-Path     Band-1    Band-2    Band-3    Band-4      ...     Band-n\n");
		for (int i = 0; i < this->num_kpts; i++)
		{
			fprintf(fp2_dw, "  %lf  ", kline[i - this->start_point]);
			for (int j = 0; j < this->num_bands; j++)
				fprintf(fp2_dw, "%lf  ", this->energy_dw[i][j] - efermi);
			fprintf(fp2_dw, "\n");
		}
		fclose(fp2_dw);
	}
	convert_basic_band("BAND.dat", ISPIN);
	printf("Written BAND.spa file!\n");
}

void BAND::Aband(int num, const char elemlabel[], int headind, int ISPIN, int LORBIT, int max_element, vector<double> kline, double efermi)
{
	if (num <= 0 || num > this->num_ions)
	{
		printf("THE ATOM_INDEX IS ERROR! Please input atom_index in 1 <= & <= %d\n", this->num_ions);
		return;
	}
	if (ISPIN == 1)
	{
		char file[20];
		if (elemlabel == NULL)
			sprintf(file, "PBAND_A%d.dat", num);
		else
			sprintf(file, "PBAND_%s%d.dat", elemlabel, num - headind);
		FILE* fp = fopen(file, "w");
		if (LORBIT == 10)
		{
			if (max_element <= 57)
			{
				fprintf(fp, "#K-Path          Energy     s     p     d     tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 4; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 4; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
			else
			{
				fprintf(fp, "#K-Path          Energy     s     p     d     f     tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 5; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 5; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
		}
		else if (LORBIT == 11)
		{
			if (max_element <= 57)
			{
				fprintf(fp, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 10; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 10; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
			else
			{
				fprintf(fp, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 17; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							for (int k = 0; k < 17; k++)
								fprintf(fp, "	%lf", this->ion_dos[j][i][num - 1][k]);
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
		}
		fclose(fp);
		printf("Written %s file!\n", file);
		string spa = convert_fat_band(file);
		printf("Written %s File!\n", spa.c_str());
	}
	else if (ISPIN == 2)
	{
		char file1[20];
		char file2[20];
		if (elemlabel == NULL)
		{
			sprintf(file1, "PBAND_A%d_UP.dat", num);
			sprintf(file2, "PBAND_A%d_DW.dat", num);
		}
		else
		{
			sprintf(file1, "PBAND_%s%d_UP.dat", elemlabel, num - headind);
			sprintf(file2, "PBAND_%s%d_DW.dat", elemlabel, num - headind);
		}
		FILE* fp1 = fopen(file1, "w");
		FILE* fp2 = fopen(file2, "w");
		if (LORBIT == 10)
		{
			if (max_element <= 57)
			{
				fprintf(fp1, "#K-Path          Energy     s     p     d     tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     p     d     tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 4; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 4; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
			else
			{
				fprintf(fp1, "#K-Path          Energy     s     p     d     f     tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     p     d     f     tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 5; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 5; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
		}
		else if (LORBIT == 11)
		{
			if (max_element <= 57)
			{
				fprintf(fp1, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 10; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 10; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
			else
			{
				fprintf(fp1, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 17; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							for (int k = 0; k < 17; k++)
							{
								fprintf(fp1, "	%lf", this->ion_dos_up[j][i][num - 1][k]);
								fprintf(fp2, "	%lf", this->ion_dos_dw[j][i][num - 1][k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
		}
		fclose(fp1);
		fclose(fp2);
		printf("Written %s file!\n", file1);
		printf("Written %s file!\n", file2);
		string spa1 = convert_fat_band(file1);
		printf("Written %s File!\n", spa1.c_str());
		string spa2 = convert_fat_band(file2);
		printf("Written %s File!\n", spa2.c_str());
	}
}

void BAND::Eband(vector<string> element, int typenum[], char elemsym[][3], int nant[], int max_element,
	int ISPIN, int LORBIT, vector<double> kline, double efermi)
{
	int cnt = 0;
	char(*file1)[20];
	file1 = new char[nant[0]][20];
	char(*file2)[20];
	file2 = new char[nant[0]][20];
	for (int m = 0; m < nant[1]; m++)
	{
		if (find(element.begin(), element.end(), elemsym[m]) == element.end())
		{
			cnt += typenum[m];
			continue;
		}
		else
		{
			if (ISPIN == 1)
			{
				char(*file)[20];
				file = new char[nant[0]][20];
				sprintf(file[m], "PBAND_%s.dat", elemsym[m]);
				FILE* fp = fopen(file[m], "w");
				if (LORBIT == 10)
				{
					if (max_element <= 57)
					{
						fprintf(fp, "#K-Path          Energy     s     p     d     tot\n");
						fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(4, 0);
									for (int k = 0; k < 4; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(4, 0);
									for (int k = 0; k < 4; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							fprintf(fp, "\n");
						}
					}
					else if (max_element > 57)
					{
						fprintf(fp, "#K-Path          Energy     s     p     d     f     tot\n");
						fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(5, 0);
									for (int k = 0; k < 5; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(5, 0);
									for (int k = 0; k < 5; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							fprintf(fp, "\n");
						}
					}
				}
				else if (LORBIT == 11)
				{
					if (max_element <= 57)
					{
						fprintf(fp, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
						fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(10, 0);
									for (int k = 0; k < 10; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(10, 0);
									for (int k = 0; k < 10; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							fprintf(fp, "\n");
						}
					}
					else if (max_element > 57)
					{
						fprintf(fp, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
						fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(17, 0);
									for (int k = 0; k < 17; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
									vector<double> temp_dos(17, 0);
									for (int k = 0; k < 17; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
											temp_dos[k] += this->ion_dos[j][i][n + cnt][k];
										fprintf(fp, "	%lf", temp_dos[k]);
									}
									fprintf(fp, "\n");
								}
							}
							fprintf(fp, "\n");
						}
					}
				}
				fclose(fp);
				printf("Written %s file!\n", file[m]);
				string spa = convert_fat_band(file[m]);
				printf("Written %s File!\n", spa.c_str());
			}
			else if (ISPIN == 2)
			{
				sprintf(file1[m], "PBAND_%s_UP.dat", elemsym[m]);
				FILE* fp1 = fopen(file1[m], "w");
				sprintf(file2[m], "PBAND_%s_DW.dat", elemsym[m]);
				FILE* fp2 = fopen(file2[m], "w");
				if (LORBIT == 10)
				{
					if (max_element <= 57)
					{
						fprintf(fp1, "#K-Path          Energy     s     p     d     tot\n");
						fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						fprintf(fp2, "#K-Path          Energy     s     p     d     tot\n");
						fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp1, "# Band-Index   %d\n", i + 1);
							fprintf(fp2, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(4, 0);
									vector<double> temp_dos_dw(4, 0);
									for (int k = 0; k < 4; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(4, 0);
									vector<double> temp_dos_dw(4, 0);
									for (int k = 0; k < 4; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else if (max_element > 57)
					{
						fprintf(fp1, "#K-Path          Energy     s     p     d     f     tot\n");
						fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						fprintf(fp2, "#K-Path          Energy     s     p     d     f     tot\n");
						fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp1, "# Band-Index   %d\n", i + 1);
							fprintf(fp2, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(5, 0);
									vector<double> temp_dos_dw(5, 0);
									for (int k = 0; k < 5; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(5, 0);
									vector<double> temp_dos_dw(5, 0);
									for (int k = 0; k < 5; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
				}
				else if (LORBIT == 11)
				{
					if (max_element <= 57)
					{
						fprintf(fp1, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
						fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						fprintf(fp2, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
						fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp1, "# Band-Index   %d\n", i + 1);
							fprintf(fp2, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(10, 0);
									vector<double> temp_dos_dw(10, 0);
									for (int k = 0; k < 10; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(10, 0);
									vector<double> temp_dos_dw(10, 0);
									for (int k = 0; k < 10; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else if (max_element > 57)
					{
						fprintf(fp1, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
						fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						fprintf(fp2, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
						fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
						for (int i = 0; i < this->num_bands; i++)
						{
							fprintf(fp1, "# Band-Index   %d\n", i + 1);
							fprintf(fp2, "# Band-Index   %d\n", i + 1);
							if (i % 2 == 0)
							{
								for (int j = this->start_point; j < this->num_kpts; j++)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(10, 0);
									vector<double> temp_dos_dw(10, 0);
									for (int k = 0; k < 17; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							else
							{
								for (int j = this->num_kpts - 1; j >= this->start_point; j--)
								{
									fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
									fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
									vector<double> temp_dos_up(17, 0);
									vector<double> temp_dos_dw(17, 0);
									for (int k = 0; k < 17; k++)
									{
										for (int n = 0; n < typenum[m]; n++)
										{
											temp_dos_up[k] += this->ion_dos_up[j][i][n + cnt][k];
											temp_dos_dw[k] += this->ion_dos_dw[j][i][n + cnt][k];
										}
										fprintf(fp1, "	%lf", temp_dos_up[k]);
										fprintf(fp2, "	%lf", temp_dos_dw[k]);
									}
									fprintf(fp1, "\n");
									fprintf(fp2, "\n");
								}
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
				}
				fclose(fp1);
				fclose(fp2);
				printf("Written %s file!\n", file1[m]);
				printf("Written %s file!\n", file2[m]);
				string spa1 = convert_fat_band(file1[m]);
				printf("Written %s File!\n", spa1.c_str());
				string spa2 = convert_fat_band(file2[m]);
				printf("Written %s File!\n", spa2.c_str());
			}
			cnt += typenum[m];
		}
	}
	delete[] file1;
	delete[] file2;
}

void BAND::Sband(vector<int> num, string postfix, int max_element, int ISPIN, int LORBIT, vector<double> kline, double efermi)
{
	if (ISPIN == 1)
	{
		char sum_file[50];
		sprintf(sum_file, "PBAND_SUM_%s.dat", postfix.c_str());
		FILE* fp = fopen(sum_file, "w");
		if (LORBIT == 10)
		{
			if (max_element <= 57)
			{
				fprintf(fp, "#K-Path          Energy     s     p     d     tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(4, 0);
							for (int k = 0; k < 4; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(4, 0);
							for (int k = 0; k < 4; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
			else if (max_element > 57)
			{
				fprintf(fp, "#K-Path          Energy     s     p     d     f     tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(5, 0);
							for (int k = 0; k < 5; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(5, 0);
							for (int k = 0; k < 5; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
		}
		else if (LORBIT == 11)
		{
			if (max_element <= 57)
			{
				fprintf(fp, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(10, 0);
							for (int k = 0; k < 10; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(10, 0);
							for (int k = 0; k < 10; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
			else if (max_element > 57)
			{
				fprintf(fp, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
				fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(17, 0);
							for (int k = 0; k < 17; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp, " %lf	%lf", kline[j - this->start_point], this->energy[j][i] - efermi);
							vector<double> temp_dos(17, 0);
							for (int k = 0; k < 17; k++)
							{
								for (int n = 0; n < num.size(); n++)
									temp_dos[k] += this->ion_dos[j][i][num[n] - 1][k];
								fprintf(fp, "	%lf", temp_dos[k]);
							}
							fprintf(fp, "\n");
						}
					}
					fprintf(fp, "\n");
				}
			}
		}
		fclose(fp);
		printf("Written %s file!\n", sum_file);
		string spa = convert_fat_band(sum_file);
		printf("Written %s File!\n", spa.c_str());
	}
	else if (ISPIN == 2)
	{
		char sum_file1[50];
		char sum_file2[50];
		sprintf(sum_file1, "PBAND_SUM_UP_%s.dat", postfix.c_str());
		FILE* fp1 = fopen(sum_file1, "w");
		sprintf(sum_file2, "PBAND_SUM_DW_%s.dat", postfix.c_str());
		FILE* fp2 = fopen(sum_file2, "w");
		if (LORBIT == 10)
		{
			if (max_element <= 57)
			{
				fprintf(fp1, "#K-Path          Energy     s     p     d     tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     p     d     tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(4, 0);
							vector<double> temp_dos_dw(4, 0);
							for (int k = 0; k < 4; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(4, 0);
							vector<double> temp_dos_dw(4, 0);
							for (int k = 0; k < 4; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
			else if (max_element > 57)
			{
				fprintf(fp1, "#K-Path          Energy     s     p     d     f     tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     p     d     f     tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(5, 0);
							vector<double> temp_dos_dw(5, 0);
							for (int k = 0; k < 5; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(5, 0);
							vector<double> temp_dos_dw(5, 0);
							for (int k = 0; k < 5; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
		}
		else if (LORBIT == 11)
		{
			if (max_element <= 57)
			{
				fprintf(fp1, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(10, 0);
							vector<double> temp_dos_dw(10, 0);
							for (int k = 0; k < 10; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(10, 0);
							vector<double> temp_dos_dw(10, 0);
							for (int k = 0; k < 10; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
			else if (max_element > 57)
			{
				fprintf(fp1, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
				fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				fprintf(fp2, "#K-Path          Energy     s     py     pz     px    dxy    dyz    dz2    dxz    dx2    f1    f2    f3    f4    f5    f6    f7    tot\n");
				fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
				for (int i = 0; i < this->num_bands; i++)
				{
					fprintf(fp1, "# Band-Index   %d\n", i + 1);
					fprintf(fp2, "# Band-Index   %d\n", i + 1);
					if (i % 2 == 0)
					{
						for (int j = this->start_point; j < this->num_kpts; j++)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(10, 0);
							vector<double> temp_dos_dw(10, 0);
							for (int k = 0; k < 17; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					else
					{
						for (int j = this->num_kpts - 1; j >= this->start_point; j--)
						{
							fprintf(fp1, " %lf	%lf", kline[j - this->start_point], this->energy_up[j][i] - efermi);
							fprintf(fp2, " %lf	%lf", kline[j - this->start_point], this->energy_dw[j][i] - efermi);
							vector<double> temp_dos_up(17, 0);
							vector<double> temp_dos_dw(17, 0);
							for (int k = 0; k < 17; k++)
							{
								for (int n = 0; n < num.size(); n++)
								{
									temp_dos_up[k] += this->ion_dos_up[j][i][num[n] - 1][k];
									temp_dos_dw[k] += this->ion_dos_dw[j][i][num[n] - 1][k];
								}
								fprintf(fp1, "	%lf", temp_dos_up[k]);
								fprintf(fp2, "	%lf", temp_dos_dw[k]);
							}
							fprintf(fp1, "\n");
							fprintf(fp2, "\n");
						}
					}
					fprintf(fp1, "\n");
					fprintf(fp2, "\n");
				}
			}
		}
		fclose(fp1);
		fclose(fp2);
		printf("Written %s file!\n", sum_file1);
		printf("Written %s file!\n", sum_file2);
		string spa1 = convert_fat_band(sum_file1);
		printf("Written %s File!\n", spa1.c_str());
		string spa2 = convert_fat_band(sum_file2);
		printf("Written %s File!\n", spa2.c_str());
	}
}

void BAND::Oband(vector<int> num, vector<int> orbit_ind, vector<double> kline, int ISPIN, int LORBIT, double efermi)
{
	if (ISPIN == 1)
	{
		char file[] = "PBAND_SUM.dat";
		FILE* fp = fopen(file, "w");
		fprintf(fp, "#K-Path          Energy      Orbital-Weight\n");
		fprintf(fp, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
		for (int i = 0; i < this->num_bands; i++)
		{
			fprintf(fp, "# Band-Index    %d\n", i + 1);
			if (i % 2 == 0)
			{
				for (int j = this->start_point; j < this->num_kpts; j++)
				{
					double temp = 0;
					for (int m = 0; m < num.size(); m++)
					{
						for (int n = 0; n < orbit_ind.size(); n++)
						{
							temp += this->ion_dos[j][i][num[m] - 1][orbit_ind[n]];
						}
					}
					fprintf(fp, "   %lf    %lf     %.3lf\n", kline[j - this->start_point], this->energy[j][i] - efermi, temp);
				}
			}
			else
			{
				for (int j = this->num_kpts - 1; j >= this->start_point; j--)
				{
					double temp = 0;
					for (int m = 0; m < num.size(); m++)
					{
						for (int n = 0; n < orbit_ind.size(); n++)
						{
							temp += this->ion_dos[j][i][num[m] - 1][orbit_ind[n]];
						}
					}
					fprintf(fp, "   %lf    %lf     %.3lf\n", kline[j - this->start_point], this->energy[j][i] - efermi, temp);
				}
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		printf("Written %s file!\n", file);
		string spa = convert_fat_band(file);
		printf("Written %s File!\n", spa.c_str());
	}
	else
	{
		char file1[] = "PBAND_SUM_UP.dat";
		char file2[] = "PBAND_SUM_DW.dat";
		FILE* fp1 = fopen(file1, "w");
		fprintf(fp1, "#K-Path          Energy      Orbital-Weight\n");
		fprintf(fp1, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
		for (int i = 0; i < this->num_bands; i++)
		{
			fprintf(fp1, "# Band-Index    %d\n", i + 1);
			if (i % 2 == 0)
			{
				for (int j = this->start_point; j < this->num_kpts; j++)
				{
					double temp = 0;
					for (int m = 0; m < num.size(); m++)
					{
						for (int n = 0; n < orbit_ind.size(); n++)
						{
							temp += this->ion_dos_up[j][i][num[m] - 1][orbit_ind[n]];
						}
					}
					fprintf(fp1, "   %lf    %lf     %.3lf\n", kline[j - this->start_point], this->energy_up[j][i] - efermi, temp);
				}
			}
			else
			{
				for (int j = this->num_kpts - 1; j >= this->start_point; j--)
				{
					double temp = 0;
					for (int m = 0; m < num.size(); m++)
					{
						for (int n = 0; n < orbit_ind.size(); n++)
						{
							temp += this->ion_dos_up[j][i][num[m] - 1][orbit_ind[n]];
						}
					}
					fprintf(fp1, "   %lf    %lf     %.3lf\n", kline[j - this->start_point], this->energy_up[j][i] - efermi, temp);
				}
			}
			fprintf(fp1, "\n");
		}
		fclose(fp1);
		printf("Written %s file!\n", file1);
		string spa1 = convert_fat_band(file1);
		printf("Written %s File!\n", spa1.c_str());
		FILE* fp2 = fopen(file2, "w");
		fprintf(fp2, "#K-Path          Energy      Orbital-Weight\n");
		fprintf(fp2, "# NKPTS & NBANDS:  %d  %d\n", this->num_kpts - this->start_point, this->num_bands);
		for (int i = 0; i < this->num_bands; i++)
		{
			fprintf(fp2, "# Band-Index    %d\n", i + 1);
			if (i % 2 == 0)
			{
				for (int j = this->start_point; j < this->num_kpts; j++)
				{
					double temp = 0;
					for (int m = 0; m < num.size(); m++)
					{
						for (int n = 0; n < orbit_ind.size(); n++)
						{
							temp += this->ion_dos_dw[j][i][num[m] - 1][orbit_ind[n]];
						}
					}
					fprintf(fp1, "   %lf    %lf     %.3lf\n", kline[j - this->start_point], this->energy_dw[j][i] - efermi, temp);
				}
			}
			else
			{
				for (int j = this->num_kpts - 1; j >= this->start_point; j--)
				{
					double temp = 0;
					for (int m = 0; m < num.size(); m++)
					{
						for (int n = 0; n < orbit_ind.size(); n++)
						{
							temp += this->ion_dos_dw[j][i][num[m] - 1][orbit_ind[n]];
						}
					}
					fprintf(fp1, "   %lf    %lf     %.3lf\n", kline[j - this->start_point], this->energy_dw[j][i] - efermi, temp);
				}
			}
			fprintf(fp2, "\n");
		}
		fclose(fp2);
		printf("Written %s file!\n", file2);
		string spa2 = convert_fat_band(file2);
		printf("Written %s File!\n", spa2.c_str());
	}
}

int Ifconduction_band(vector<double> band, double efermi)
{
	int flag = 0;
	if (band[0] >= efermi)
		flag = 1; // conduction
	else
		flag = -1; // Valence
	for (int i = 1; i < band.size(); i++)
	{
		if ((flag == 1 && band[i] <= efermi) || (flag == -1 && band[i] >= efermi))
		{
			flag = 0;
			break;
		}
	}
	return flag;
}

vector<vector<int> > BAND::BandGap(int ISPIN, vector<double> kline, double efermi, bool write)
{
	int metal = 0;
	vector<vector<int> > ret;
	if (ISPIN == 1)
	{
		double VBM = -1e5, CBM = 1e5;
		vector<int> v, c;
		double V_location[3], C_location[3];
		int lumo = 0, homo = 0;
		for (int i = 0; i < this->energy[0].size() - 1; i++)
		{
			if (this->energy[0][i]<efermi && this->energy[0][i + 1]>efermi)
			{
				homo = i + 1;
				lumo = i + 2;
			}
		}
		vector <vector<double>> band_col(this->num_bands);
		for (int j = 0; j < this->num_bands; j++)
			for (int i = this->start_point; i < this->num_kpts; i++)
				band_col[j].push_back(this->energy[i][j]);
		for (int i = 0; i < this->num_bands; i++)
		{
			if (Ifconduction_band(band_col[i], efermi) == -1) // Valence
			{
				for (int j = 0; j < band_col[i].size(); j++)
				{
					if (band_col[i][j] > VBM)
					{
						v.clear(); //update
						VBM = band_col[i][j];
						v.push_back(j);
					}
					else if (band_col[i][j] == VBM)
						v.push_back(j);
				}
			}
			else if (Ifconduction_band(band_col[i], efermi) == 1) // conduction
			{
				for (int j = 0; j < band_col[i].size(); j++)
				{
					if (band_col[i][j] < CBM)
					{
						c.clear(); //update
						CBM = band_col[i][j];
						c.push_back(j);
					}
					else if (band_col[i][j] == CBM)
						c.push_back(j);
				}
			}
			else
				metal = 1;
		}
		double band_gap = CBM - VBM;
		ret.push_back(vector<int> {homo, lumo});
		ret.push_back(c);
		ret.push_back(v);
		if (!write)
			return ret;
		FILE* fp = fopen("BAND_GAP", "w");
		fprintf(fp, " ------------------------ Band Information -------------------------\n");
		if (metal == 1)
		{
			fprintf(fp, "          Band Character:    Metallic\n");
			fprintf(fp, " -------------------------------------------------------------------\n");
			fclose(fp);
		}
		else
		{
			int direct = 0;
			for (int i = 0; i < v.size(); i++)
				for (int j = 0; j < c.size(); j++)
				{
					if (v[i] == c[j])
					{
						direct = 1;
						for (int k = 0; k < 3; k++)
						{
							V_location[k] = this->k_coor[v[i]+ start_point][k];
							C_location[k] = this->k_coor[c[j]+ start_point][k];
						}
						break;
					}
				}
			if (!direct)
			{
				for (int i = 0; i < 3; i++)
				{
					V_location[i] = this->k_coor[v[0]+ start_point][i];
					C_location[i] = this->k_coor[c[0]+ start_point][i];
				}
			}
			fprintf(fp, "          Band Gap (eV):  %lf\n", band_gap);
			fprintf(fp, "         Band Character:  %s\n", direct ? "Direct" : "Indirect");
			fprintf(fp, "      Fermi Energy (eV):  %lf\n", efermi);
			fprintf(fp, "     VBM Eigenvalue(eV): % lf\n", VBM);
			fprintf(fp, "     CBM Eigenvalue(eV):  %lf\n", CBM);
			fprintf(fp, "  Highest-Occupied Band:         %d\n", homo);
			fprintf(fp, "   Lowest-Occupied Band:         %d\n", lumo);
			fprintf(fp, "        Location of VBM:  %lf  %lf  %lf\n", V_location[0], V_location[1], V_location[2]);
			fprintf(fp, "        Location of CBM:  %lf  %lf  %lf\n", C_location[0], C_location[1], C_location[2]);
			fprintf(fp, " -------------------------------------------------------------------\n");
			fclose(fp);
		}
	}
	else
	{
		int metal_up = 0, metal_dw = 0;
		double VBM_up = -1e5, CBM_up = 1e5, VBM_dw = -1e5, CBM_dw = 1e5;
		vector<int> v_up, c_up, v_dw, c_dw;
		double V_location_up[3] , C_location_up[3] , V_location_dw[3] , C_location_dw[3] ;
		vector <vector<double>> band_col_up(this->num_bands);
		vector <vector<double>> band_col_dw(this->num_bands);
		int lumo_up = 0, homo_up = 0;
		for (int i = 0; i < this->energy_up[0].size() - 1; i++)
		{
			if (this->energy_up[0][i]<efermi && this->energy_up[0][i + 1]>efermi)
			{
				homo_up = i + 1;
				lumo_up = i + 2;
			}
		}
		int lumo_dw = 0, homo_dw = 0;
		for (int i = 0; i < this->energy_dw[0].size() - 1; i++)
		{
			if (this->energy_dw[0][i]<efermi && this->energy_dw[0][i + 1]>efermi)
			{
				homo_dw = i + 1;
				lumo_dw = i + 2;
			}
		}
		for (int j = 0; j < this->num_bands; j++)
			for (int i = this->start_point; i < this->num_kpts; i++)
			{
				band_col_up[j].push_back(this->energy_up[i][j]);
				band_col_dw[j].push_back(this->energy_dw[i][j]);
			}
		for (int i = 0; i < this->num_bands; i++)
		{
			if (Ifconduction_band(band_col_up[i], efermi) == -1) // Valence
			{
				for (int j = 0; j < band_col_up[i].size(); j++)
				{
					if (band_col_up[i][j] > VBM_up)
					{
						v_up.clear(); //update
						VBM_up = band_col_up[i][j];
						v_up.push_back(j);
					}
					else if (band_col_up[i][j] == VBM_up)
						v_up.push_back(j);
				}
			}
			else if (Ifconduction_band(band_col_up[i], efermi) == 1) // conduction
			{
				for (int j = 0; j < band_col_up[i].size(); j++)
				{
					if (band_col_up[i][j] < CBM_up)
					{
						c_up.clear(); //update
						CBM_up = band_col_up[i][j];
						c_up.push_back(j);
					}
					else if (band_col_up[i][j] == CBM_up)
						c_up.push_back(j);
				}
			}
			else
				metal_up = 1;
		}
		for (int i = 0; i < this->num_bands; i++)
		{
			if (Ifconduction_band(band_col_dw[i], efermi) == -1) // Valence
			{
				for (int j = 0; j < band_col_dw[i].size(); j++)
				{
					if (band_col_dw[i][j] > VBM_dw)
					{
						v_dw.clear(); //dwdate
						VBM_dw = band_col_dw[i][j];
						v_dw.push_back(j);
					}
					else if (band_col_dw[i][j] == VBM_dw)
						v_dw.push_back(j);
				}
			}
			else if (Ifconduction_band(band_col_dw[i], efermi) == 1) // conduction
			{
				for (int j = 0; j < band_col_dw[i].size(); j++)
				{
					if (band_col_dw[i][j] < CBM_dw)
					{
						c_dw.clear(); //dwdate
						CBM_dw = band_col_dw[i][j];
						c_dw.push_back(j);
					}
					else if (band_col_dw[i][j] == CBM_dw)
						c_dw.push_back(j);
				}
			}
			else
				metal_dw = 1;
		}
		ret.push_back(vector<int> {homo_up, homo_dw, lumo_up, lumo_dw});
		ret.push_back(v_up); ret.push_back(v_dw);
		ret.push_back(c_up); ret.push_back(c_dw);
		if (!write)
			return ret;
		double band_gap_up = CBM_up - VBM_up;
		double band_gap_dw = CBM_dw - VBM_dw;
		FILE* fp = fopen("BAND_GAP", "w");
		fprintf(fp, " ------------------------ Band Information -------------------------\n");
		if (metal_up == 1 || metal_dw == 1)
		{
			fprintf(fp, "          Band Character:    Metallic\n");
			fprintf(fp, " -------------------------------------------------------------------\n");
			fclose(fp);
		}
		else
		{
			int direct_up = 0;
			for (int i = 0; i < v_up.size(); i++)
				for (int j = 0; j < c_up.size(); j++)
				{
					if (v_up[i] == c_up[j])
					{
						direct_up = 1;
						for (int k = 0; k < 3; k++)
						{
							V_location_up[k] = this->k_coor[i][k];
							C_location_up[k] = this->k_coor[j][k];
						}
						break;
					}
				}
			if (!direct_up)
			{
				for (int i = 0; i < 3; i++)
				{
					V_location_up[i] = this->k_coor[v_up[0]][i];
					C_location_up[i] = this->k_coor[c_up[0]][i];
				}
			}
			int direct_dw = 0;
			for (int i = 0; i < v_dw.size(); i++)
				for (int j = 0; j < c_dw.size(); j++)
				{
					if (v_dw[i] == c_dw[j])
					{
						direct_dw = 1;
						for (int k = 0; k < 3; k++)
						{
							V_location_dw[k] = this->k_coor[i][k];
							C_location_dw[k] = this->k_coor[j][k];
						}
						break;
					}
				}
			if (!direct_dw)
			{
				for (int i = 0; i < 3; i++)
				{
					V_location_dw[i] = this->k_coor[v_dw[0]][i];
					C_location_dw[i] = this->k_coor[c_dw[0]][i];
				}
			}
			double CBM_t, VBM_t;
			int lumo_t = 0, homo_t = 0;
			double V_location_t[3] , C_location_t[3] ;
			if (CBM_up > CBM_dw)
			{
				CBM_t = CBM_dw;
				lumo_t = lumo_dw;
				for (int i = 0; i < 3; i++)
					C_location_t[i] = C_location_dw[i];
			}
			else
			{
				CBM_t = CBM_up;
				lumo_t = lumo_up;
				for (int i = 0; i < 3; i++)
					C_location_t[i] = C_location_up[i];
			}
			if (VBM_up < VBM_dw)
			{
				VBM_t = VBM_dw;
				homo_t = homo_dw;
				for (int i = 0; i < 3; i++)
					V_location_t[i] = V_location_dw[i];
			}
			else
			{
				VBM_t = VBM_up;
				homo_t = homo_up;
				for (int i = 0; i < 3; i++)
					V_location_t[i] = V_location_up[i];
			}
			double band_gap_t = CBM_t - VBM_t;
			fprintf(fp, "            Spin Channel:     <UP>       <DOWN>      <TOTAL>\n");
			fprintf(fp, "           Band Gap (eV):  %lf       %lf       %lf\n", band_gap_up, band_gap_dw, band_gap_t);
			fprintf(fp, "          Band Character:  %s       %s\n", direct_up ? "Direct" : "Indirect", direct_dw ? "Direct" : "Indirect");
			fprintf(fp, "       Fermi Energy (eV):  %lf       %lf       %lf\n", efermi, efermi, efermi);
			fprintf(fp, "     VBM Eigenvalue (eV):  %lf       %lf       %lf\n", VBM_up, VBM_dw, VBM_t);
			fprintf(fp, "     CBM Eigenvalue (eV):  %lf       %lf       %lf\n", CBM_up, CBM_dw, CBM_t);
			fprintf(fp, "   Highest-Occupied Band:        %d        %d          %d\n", homo_up, homo_dw, homo_t);
			fprintf(fp, "  Lowest-Unoccupied Band:        %d        %d          %d\n", lumo_up, lumo_dw, lumo_t);
			fprintf(fp, "    Location of VBM (UP):  %lf  %lf  %lf\n", V_location_up[0], V_location_up[1], V_location_up[2]);
			fprintf(fp, "    Location of CBM (UP):  %lf  %lf  %lf\n", C_location_up[0], C_location_up[1], C_location_up[2]);
			fprintf(fp, "    Location of VBM (DW):  %lf  %lf  %lf\n", V_location_dw[0], V_location_dw[1], V_location_dw[2]);
			fprintf(fp, "    Location of CBM (DW):  %lf  %lf  %lf\n", C_location_dw[0], C_location_dw[1], C_location_dw[2]);
			fprintf(fp, "    Location of VBM (TO):  %lf  %lf  %lf\n", V_location_t[0], V_location_t[1], V_location_t[2]);
			fprintf(fp, "    Location of CBM (TO):  %lf  %lf  %lf\n", C_location_t[0], C_location_t[1], C_location_t[2]);
			fprintf(fp, " ------------------------------------- -----------------------------\n");
			fclose(fp);
		}
	}
	printf("Written BAND_GAP file!\n");
	return ret;
}

double BAND::Effective_Mass(int fitting_point, int position, int band_index, bool direct, vector<double> kline, double efermi)
{
	vector<double> energy(fitting_point), kpath(fitting_point);
	if (direct) // right
	{
		int cnt = 0;
		for (; cnt < fitting_point && position < this->num_kpts; )
		{
			energy[cnt] = (this->energy[position][band_index] - efermi) * Hartree;
			kpath[cnt] = kline[position] * Bohr;
			position++; cnt++;
		}
	}
	else
	{
		int cnt = 0;
		for (; cnt < fitting_point && position > 0;)
		{
			energy[cnt] = (this->energy[position][band_index] - efermi) * Hartree;
			kpath[cnt] = kline[position] * Bohr;
			position--; cnt++;
		}
	}
	double coefficient[fitting_point];
	EMatrix(kpath, energy, fitting_point, 3, coefficient);
	double em = 1 / (2 * coefficient[3]);
	return em;
}

typedef struct INPUT
{
	vector<int> num;
	int head_index;
	string elemlabel;
}INPUT;

void BAND::getband(int argc, char* argv[])
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
	fclose(fp);
	double n_vec[3][3];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			n_vec[i][j] = vec[i][j] * latt;
	int max_element = 0;
	for (int i = 0; i < nant[1]; i++)
		max_element = max(max_element, elemnum[i]);
	if (strcmp(argv[1], "--band"))
		return;
	string ISPIN = GetInfoINCAR("ISPIN");
	if (ISPIN.length() == 0)
		ISPIN = "1";
	string LORBIT = GetInfoINCAR("LORBIT");
	if (LORBIT.length() == 0)
		LORBIT = "0";
	readPROCAR(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	DOS dos;
	dos.readDoscar(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	double efermi = dos.fermi_energy();
	vector<double> kline;
	if (argv[2][1] == 'h' || argv[2][1] == 'H')
		kline = hse_getkpointsline(this->num_kpts, this->k_coor, n_vec);
	else
		kline = getkpointsline(n_vec);
	if (!strcmp(argv[2], "-b") || !strcmp(argv[2], "-B") || !strcmp(argv[2], "-hb") || !strcmp(argv[2], "-HB"))
	{
		basicband(atoi(ISPIN.c_str()), kline, efermi);
	}
	if (!strcmp(argv[2], "-a") || !strcmp(argv[2], "-A") || !strcmp(argv[2], "-ha") || !strcmp(argv[2], "-HA")) // output each atom
	{
		for (int i = 0; i < nant[0]; i++)
			Aband(i + 1, NULL, 0, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element, kline, efermi);
	}
	if (!strcmp(argv[2], "-e") || !strcmp(argv[2], "-E") || !strcmp(argv[2], "-he") || !strcmp(argv[2], "-HE")) // output each element
	{
		vector<string> element;
		element.resize(nant[1]);
		for (int i = 0; i < nant[1]; i++)
			element.push_back(elemsym[i]);
		Eband(element, typenum, elemsym, nant, max_element, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), kline, efermi);
	}
	if (!strcmp(argv[2], "-s") || !strcmp(argv[2], "-S") || !strcmp(argv[2], "-hs") || !strcmp(argv[2], "-HS")) // VASPMATE --band -s <atom-index> <element-index>
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
			Aband(num[i], NULL, 0, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element, kline, efermi);
		for (int i = 0; i < data.size(); i++)
			for (int j = 0; j < data[i].num.size(); j++)
				Aband(data[i].num[j], data[i].elemlabel.c_str(), data[i].head_index, atoi(ISPIN.c_str()),
					atoi(LORBIT.c_str()), max_element, kline, efermi);
	}
	if (!strcmp(argv[2], "-sa") || !strcmp(argv[2], "-SA") || !strcmp(argv[2], "-hsa") || !strcmp(argv[2], "-HSA"))	// output select atoms
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
			Aband(num[i], NULL, 0, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), max_element, kline, efermi);
	}
	if (!strcmp(argv[2], "-se") || !strcmp(argv[2], "-SE") || !strcmp(argv[2], "-hsa") || !strcmp(argv[2], "-HSA"))	// output select elements
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
			{
				if (!strcmp(argv[2], "-se") || !strcmp(argv[2], "-SE"))
					Aband(data[i].num[j], data[i].elemlabel.c_str(), data[i].head_index, atoi(ISPIN.c_str()),
						atoi(LORBIT.c_str()), max_element, kline, efermi);
			}
	}
	if (!strcmp(argv[2], "-m") || !strcmp(argv[2], "-M") || !strcmp(argv[2], "-ma") || !strcmp(argv[2], "-me")
		|| !strcmp(argv[2], "-hm") || !strcmp(argv[2], "-HM") || !strcmp(argv[2], "-hma") || !strcmp(argv[2], "-hme"))	//Multiple Atoms or Elements
	{
		vector<int> num = TranArgvToAtomIndex(argc, argv, nant, typenum, elemsym, xyz, NULL);
		string postfix;
		for (int i = 3; i < argc; i++)
		{
			if (i != argc - 1)
				postfix = postfix + argv[i] + "_";
			else
				postfix += argv[i];
		}
		Sband(num, postfix, max_element, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), kline, efermi);
	}
	if (!strcmp(argv[2], "-o") || !strcmp(argv[2], "-O") || !strcmp(argv[2], "-oa") || !strcmp(argv[2], "-oe")
		|| !strcmp(argv[2], "-ho") || !strcmp(argv[2], "-HO") || !strcmp(argv[2], "-hoa") || !strcmp(argv[2], "-hoe")) //select atom and orbit VASPMATE --band -o (a b c...)/(a-b) s p d f(s px ...)
	{
		vector<string> ind;
		vector<string> orb;
		for (int i = 3; i < argc; i++)
		{
			if (!IsRightOrbit(argv[i], atoi(LORBIT.c_str()), max_element) && strcmp(argv[i], "all"))
				ind.push_back(argv[i]);
			else if (atoi(LORBIT.c_str()) == 11 && !strcmp(argv[i], "p"))
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
			else
				orb.push_back(argv[i]);
		}
		vector<int> num;
		for (int i = 0; i < ind.size(); i++)
		{
			if (IsNum(ind[i]))
			{
				if (atoi(ind[i].c_str()) <= 0 || atoi(ind[i].c_str()) > nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				num.push_back(atoi(ind[i].c_str()));
			}
			else if (strstr(ind[i].c_str(), "-") != NULL)
			{
				int start, end;
				sscanf(ind[i].c_str(), "%d%*c%d", &start, &end);
				if (start > end || start < 1 || end > nant[0])
				{
					printf("THE ATOM_INDEX IS ERROR!Please input atoms in 1 <= & <= %d\n", nant[0]);
					return;
				}
				for (int i = start; i < end + 1; i++)
					num.push_back(i);
			}
			else
			{
				//check symbol is right element
				int IsRightelem = 0;
				for (int j = 0; j < nant[1]; j++)
				{
					if (ind[i] == elemsym[j])
						IsRightelem = 1;
				}
				if (IsRightelem == 0)
				{
					printf("%s are not in element types!\n", ind[i].c_str());
					return;
				}
				//tran element to atom index
				int cnt = 0;
				for (int k = 0; k < nant[1]; k++)
				{
					if (ind[i] == elemsym[k])
					{
						for (int j = 0; j < typenum[k]; j++)
							num.push_back(j + cnt + 1);
					}
					cnt += typenum[k];
				}
			}
		}
		if (find(orb.begin(), orb.end(), "all") != orb.end())
		{
			if (atoi(LORBIT.c_str()) == 10 && max_element <= 57)
			{
				orb = { "s","p","d" };
			}
			if (atoi(LORBIT.c_str()) == 10 && max_element > 57)
			{
				orb = { "s","p","d","f" };
			}
			if (atoi(LORBIT.c_str()) == 11 && max_element <= 57)
			{
				orb = { "s","px","py","pz","dxy","dxz","dyz","dx2","dz2" };
			}
			if (atoi(LORBIT.c_str()) == 11 && max_element > 57)
			{
				orb = { "s","px","py","pz","dxy","dxz","dyz","dx2","dz2","f1","f2","f3","f4","f5","f6","f7" };
			}
		}
		vector<string> select(nant[0], "F");
		for (int i = 0; i < num.size(); i++)
			select[num[i] - 1] = "T";
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
		char select_file[] = "SELECT_ATOMS_LIST";
		FILE* fp1 = fopen(select_file, "w");
		fprintf(fp1, "ATOMS_ID ATOM_LABEL  X_POSITION   Y_POSITION   Z_POSITION   SELECTED?\n");
		for (int i = 0; i < nant[0]; i++)
			fprintf(fp1, "	%d	%s	%lf	%lf	%lf	%s\n", i + 1, elemlabel[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], select[i].c_str());
		fclose(fp1);
		printf("Written %s file!\n", select_file);
		char Select_file[] = "SELECT_ORBITS_LIST";
		FILE* fp2 = fopen(Select_file, "w");
		fprintf(fp2, "ORBITALS_ID  ORBITAL_LABEL  SELECTED?\n");
		vector<int> orbit_ind;
		vector<string> orbit;
		if (atoi(LORBIT.c_str()) == 10 && max_element <= 57)
			orbit = { "s","p","d" };
		if (atoi(LORBIT.c_str()) == 10 && max_element > 57)
			orbit = { "s","p","d","f" };
		if (atoi(LORBIT.c_str()) == 11 && max_element <= 57)
			orbit = { "s","py","pz","px","dxy","dyz","dz2","dxz","dx2" };
		if (atoi(LORBIT.c_str()) == 11 && max_element > 57)
			orbit = { "s","py","pz","px","dxy","dyz","dz2","dxz","dx2","f1","f2","f3","f4","f5","f6","f7" };
		for (int i = 0; i < orbit.size(); i++)
		{
			fprintf(fp2, "   %2d            %3s            %s\n", i + 1, orbit[i].c_str(), (find(orb.begin(), orb.end(), orbit[i]) == orb.end() ? "F" : "T"));
			if (find(orb.begin(), orb.end(), orbit[i]) != orb.end())
				orbit_ind.push_back(i);
		}
		fclose(fp2);
		printf("Written %s file!\n", Select_file);
		Oband(num, orbit_ind, kline, atoi(ISPIN.c_str()), atoi(LORBIT.c_str()), efermi);
	}
	if (!strcmp("-bg", argv[2]) || !strcmp(argv[2], "-BG") || !strcmp("-hbg", argv[2]) || !strcmp(argv[2], "-HBG"))
	{
		BandGap(atoi(ISPIN.c_str()), kline, efermi);
	}
	if (!strcmp("-em", argv[2]) || !strcmp("-hem", argv[2]))	//effective mass
	{
		if (atoi(ISPIN.c_str()) == 2)
		{
			printf("VASPMATE Currently only supported effective mass calculation for Non-Charged & Non-Magnetic Semiconductor!\n");
			return;
		}
		vector<vector<int> >info = BandGap(atoi(ISPIN.c_str()), kline, efermi, 0);
		int fittint_point = 6;
		/*	info: (vector<int> {homo, lumo}); (c); (v);	 */
		if (info.size() == 0)
		{
			printf("The Band Character is Metallic,Band_Gap < 0\n");
			return;
		}
		vector<vector<string>> pd;//point+direct
		int homo = 0, lumo = 0, select_flag = 0;
		for (int i = 3; i < argc; i += 2)
		{
			if (!strcmp(argv[i], "-np"))
			{
				fittint_point = atoi(argv[i + 1]);
				continue;
			}
			if (!strcmp(argv[i], "-nb"))
			{
				select_flag = 1;
				homo = atoi(argv[i + 1]) - 1;
				lumo = homo + 1;
				continue;
			}
			pd.push_back(vector<string>{argv[i], argv[i + 1]});
		}
		char file[50];
		if (!strcmp("-em", argv[2]))
			strcpy(file, "KPOINTS");
		else
			strcpy(file, "NEWKPATH");
		FILE* fp = fopen(file, "r");
		if (fp == NULL)
		{
			printf("%s IS NOT EXIST!\n", file);
			return;
		}
		char buf[1024];
		int line = 0;
		vector<point> kp;
		while (fgets(buf, 1024, fp) != NULL)
		{
			line++;
			if (line > 4 && strspn(buf, " \t\n\r") != strlen(buf))
			{
				double pos[3];
				char label[10];
				sscanf(buf, "%lf%lf%lf%s", &pos[0], &pos[1], &pos[2], label);
				point p(pos, label);
				kp.push_back(p);
			}
		}
		fclose(fp);
		vector<vector<vector<double> > > select_position(pd.size());
		for (int i = 0; i < pd.size(); i++)
		{
			int flag = 0;
			if (strcmp(pd[i][1].c_str(), "l") && strcmp(pd[i][1].c_str(), "r"))
			{
				printf("The second parameter must be ``l`` or ``r``!\n");
				return;
			}
			for (int j = 0; j < kp.size(); j++)
			{
				if (kp[j].label == pd[i][0])
				{
					select_position[i].push_back(kp[j].pos);
					flag = 1;
					break;
				}
			}
			if (!flag)
			{
				printf("The first parameter %s is not High-Symm-Point!\n", pd[i][0].c_str());
				return;
			}
		}
		/*	info: (vector<int> {homo, lumo}); (c); (v);	 */
		if (!select_flag)
		{
			homo = info[0][0] - 1;
			lumo = info[0][1] - 1;
		}
		FILE* fp1 = fopen("Effective_Mass", "at+");
		fprintf(fp1, "+ -------------------Summary------------------- +\n");
		fprintf(fp1, "Fitting Point: %d\n", fittint_point);
		fprintf(fp1, "Band Index:            LUMO = %d	HOMO = %d\n", lumo + 1, homo + 1);
		for (int i = 0; i < select_position.size(); i++)
		{
			int band_index = -1, _band_index = -1;
			int position = 0, _position = 0;
			int flag = 0, _flag = 0;
			bool direct = strcmp(pd[i][1].c_str(), "l") ? 1 : 0;
			for (int j = 0; j < select_position[i].size(); j++)
			{
				for (int k = 0; k < info[1].size(); k++)
					if (IsSameVector<double>(select_position[i][j], k_coor[info[1][k]]))
					{
						band_index = lumo;
						position = info[1][k];
						flag = 1;
						break;
					}
				for (int k = 0; k < info[2].size(); k++)
					if (IsSameVector<double>(select_position[i][j], k_coor[info[2][k]]))
					{
						if (!flag)
						{
							band_index = homo;
							position = info[2][k];
						}
						// VBM and CBM is same position
						else
						{
							_band_index = homo;
							_position = info[2][k];
							_flag = 1;
						}
						break;
					}
				if (band_index == -1)
				{
					printf("The %s is not VBM or CBM, please select again!\n", pd[i][0].c_str());
					fprintf(fp1, "+ -------------------Summary------------------- +\n");
					fclose(fp1);
					return;
				}
			}
			if ((position == 0 && direct == 0) || (position == this->num_kpts - 1))
			{
				fprintf(fp1, "Please check your direct.There is no data in this direct!\n");
				continue;
			}
			fprintf(fp1, "LABEL: %s Direct: %s ", pd[i][0].c_str(), direct ? "Right" : "Left");
			double em = Effective_Mass(fittint_point, position, band_index, direct, kline, efermi);
			double _em;
			if (_flag)
			{
				_em = Effective_Mass(fittint_point, _position, _band_index, direct, kline, efermi);
				fprintf(fp1, "Electron: %lf	Hole: %lf\n", em, _em);
			}
			else
			{
				if (band_index == lumo)
					fprintf(fp1, "Electron: %lf\n", em);
				else
					fprintf(fp1, "Hole: %lf\n", em);
			}
		}
		fprintf(fp1, "+ -------------------Summary------------------- +\n");
		fprintf(fp1, "\n");
		fclose(fp1);
		printf("Written Effective_Mass file!\n");
	}
}