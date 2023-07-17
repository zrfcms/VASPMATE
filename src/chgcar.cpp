#include"../include/chgcar.h"
void CHGCAR::readchgcar(const char file[])
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		printf("%s is not exist!\n", file);
		fclose(fp);
		return;
	}
	readposcar(fp, this->pos);
	char buf[1024];
	int aug_num;
	fscanf(fp, "%d %d %d", &this->chg_size[0], &this->chg_size[1], &this->chg_size[2]);
	this->chg.resize(this->chg_size[0] * this->chg_size[1] * this->chg_size[2]);
	for (int i = 0; i < this->chg_size[0] * this->chg_size[1] * this->chg_size[2]; i++)
		fscanf(fp, "%lf ", &this->chg[i]);
	this->Chg.resize(this->chg_size[0]);
	for (int i = 0; i < this->chg_size[0]; i++)
	{
		this->Chg[i].resize(this->chg_size[1]);
		for (int j = 0; j < this->chg_size[1]; j++)
			this->Chg[i][j].resize(this->chg_size[2]);
	}
	int cnt = 0;
	for (int i = 0; i < chg_size[2]; i++)
		for (int j = 0; j < chg_size[1]; j++)
			for (int k = 0; k < chg_size[0]; k++)
			{
				this->Chg[k][j][i] = this->chg[cnt];
				cnt++;
			}
	while (1)
	{
		if (fgets(buf, 1024, fp) == NULL)
		{
			fclose(fp);
			return;
		}
		else if (strstr(buf, "augmentation") == NULL)
		{
			strcpy(this->Separator, buf);
			break;
		}
		else
			sscanf(buf, "%*s%*s%*d%d", &aug_num);
		vector<double> tmp(aug_num);
		for (int i = 0; i < aug_num; i++)
			fscanf(fp, "%lf ", &tmp[i]);
		this->chg_aug.push_back(tmp);
	}
	fscanf(fp, "%d %d %d", &this->mag_size[0], &this->mag_size[1], &this->mag_size[2]);
	this->mag.resize(this->mag_size[0] * this->mag_size[1] * this->mag_size[2]);
	for (int i = 0; i < this->mag_size[0] * this->mag_size[1] * this->mag_size[2]; i++)
		fscanf(fp, "%lf ", &this->mag[i]);
	while (1)
	{
		if (fgets(buf, 1024, fp) == NULL)
		{
			fclose(fp);
			return;
		}
		else
			sscanf(buf, "%*s%*s%*d%d", &aug_num);
		vector<double> tmp(aug_num);
		for (int i = 0; i < aug_num; i++)
			fscanf(fp, "%lf ", &tmp[i]);
		this->mag_aug.push_back(tmp);
	}
	fclose(fp);
}

void CHGCAR::operator_bader(int argc, char* argv[])
{
	/*CHGCAR chg;
	chg.readchgcar(argv[2]);*/
	//writechgcar(argv[3], chg.pos, chg.chg_size, chg.chg, chg.chg_aug,
	//	chg.Separator, chg.mag_size, chg.mag, chg.mag_aug, 1);
	if (!strcmp(argv[2], "-comb"))
	{
		if (argc < 5)
			return;
		CHGCAR chg1, chg2;
		chg1.readchgcar(argv[3]);
		chg2.readchgcar(argv[4]);
		printf("Atoms in file1: %d, Atoms in file2: %d\n", chg1.pos.nant[0], chg2.pos.nant[0]);
		if (chg1.pos.nant[0] != chg2.pos.nant[0])
			printf("Atoms pos.nant[0] not same in two files!\n");
		int point1 = chg1.chg_size[0] * chg1.chg_size[1] * chg1.chg_size[2];
		int point2 = chg2.chg_size[0] * chg2.chg_size[1] * chg2.chg_size[2];
		printf("Points in file1: %d, Points in file2: %d\n", point1, point2);
		if (point1 != point2)
		{
			printf("Number of points not same in two files!\n");
			return;
		}
		CHGCAR n_chg(chg1);
		if (argc == 5)
			for (int i = 0; i < point1; i++)
				n_chg.chg[i] = chg1.chg[i] + chg2.chg[i];
		else if (argc == 6)
			for (int i = 0; i < point1; i++)
				n_chg.chg[i] = atoi(argv[5]) * chg1.chg[i] + chg2.chg[i];
		else if (argc == 7)
			for (int i = 0; i < point1; i++)
				n_chg.chg[i] = atoi(argv[5]) * chg1.chg[i] + atoi(argv[6]) * chg2.chg[i];
		n_chg.writechgcar("CHGCAR_SUM.vasp", 0);
	}
}

CHGCAR CHGCAR::operator_cdd(CHGCAR chg1, CHGCAR chg2, bool oper)
{
	CHGCAR n_chg(chg1);
	if (chg1.pos.nant[0] != chg2.pos.nant[0])
		printf("Atoms pos.nant[0] not same in two files!\n");
	int point1 = chg1.chg_size[0] * chg1.chg_size[1] * chg1.chg_size[2];
	int point2 = chg2.chg_size[0] * chg2.chg_size[1] * chg2.chg_size[2];
	printf("Points in file1: %d, Points in file2: %d\n", point1, point2);
	if (point1 != point2)
	{
		printf("Number of points not same in two files!\n");
		return chg1;
	}
	for (int i = 0; i < point1; i++)
		n_chg.chg[i] = chg1.chg[i] + (oper ? chg2.chg[i] : -chg2.chg[i]);
	return n_chg;
}

void CHGCAR::series_oper(int argc, char* argv[], bool oper)
{
	CHGCAR chg;
	chg.readchgcar(argv[3]);
	for (int i = 4; i < argc; i++)
	{
		CHGCAR chg1;
		chg1.readchgcar(argv[i]);
		chg = operator_cdd(chg, chg1, oper);
	}
	if (oper)
	{
		chg.writechgcar("CHGSUMM.vasp", 0);
		printf("Writteb CHGSUMM.vasp file!\n");
	}
	else
	{
		chg.writechgcar("CHGDIFF.vasp", 0);
		printf("Writteb CHGDIFF.vasp file!\n");
	}
}

void CHGCAR::writechgcar(const char file[], bool label)
{
	FILE* fp = fopen(file, "w");
	savposcar(fp, this->pos);
	fprintf(fp, "\n");
	fprintf(fp, "   %d   %d   %d\n", this->chg_size[0], this->chg_size[1], this->chg_size[2]);
	int cnt = 0;
	for (int i = 0; i < this->chg.size(); i++)
	{
		cnt++;
		fprintf(fp, " % E ", this->chg[i]);
		if (cnt == 5)
		{
			fprintf(fp, "\n");
			cnt = 0;
		}
	}
	if (cnt)
		fprintf(fp, "\n");
	if (!label)
	{
		fclose(fp);
		return;
	}
	else
	{
		for (int i = 0; i < this->chg_aug.size(); i++)
		{
			int aug_cnt = 0;
			fprintf(fp, "augmentation occupancies   %d %d\n", i + 1, this->chg_aug[i].size());
			for (int j = 0; j < this->chg_aug[i].size(); j++)
			{
				aug_cnt++;
				fprintf(fp, "% E ", this->chg_aug[i][j]);
				if (aug_cnt == 5)
				{
					fprintf(fp, "\n"); aug_cnt = 0;
				}
			}
			fprintf(fp, "\n");
		}
		//mag
		if (this->Separator == NULL)
		{
			fclose(fp);
			return;
		}
		else
		{
			fprintf(fp, "   %d   %d   %d\n", this->mag_size[0], this->mag_size[1], this->mag_size[2]);
			for (int i = 0; i < this->mag.size(); i++)
			{
				cnt++;
				fprintf(fp, "% E ", this->mag[i]);
				if (cnt == 5)
				{
					fprintf(fp, "\n");
					cnt = 0;
				}
			}
			if (cnt)
				fprintf(fp, "\n");
			for (int i = 0; i < this->mag_aug.size(); i++)
			{
				int aug_cnt = 0;
				fprintf(fp, "augmentation occupancies   %d %d\n", i + 1, this->mag_aug[i].size());
				for (int j = 0; j < this->mag_aug[i].size(); j++)
				{
					aug_cnt++;
					fprintf(fp, "% E ", this->mag_aug[i][j]);
					if (aug_cnt == 5)
					{
						fprintf(fp, "\n"); aug_cnt = 0;
					}
				}
				fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
}

void writechgcar(const char file[], POSCAR pos, int chg_size[3], vector<double> chg, const vector<vector<double> >& chg_aug,
	char Separator[100], int mag_size[3], const vector<double>& mag , const vector<vector<double> >& mag_aug, bool label)
{
	FILE* fp = fopen(file, "w");
	savposcar(fp, pos);
	fprintf(fp, "\n");
	fprintf(fp, "   %d   %d   %d\n", chg_size[0], chg_size[1], chg_size[2]);
	int cnt = 0;
	for (int i = 0; i < chg.size(); i++)
	{
		cnt++;
		fprintf(fp, " % E ", chg[i]);
		if (cnt == 5)
		{
			fprintf(fp, "\n");
			cnt = 0;
		}
	}
	if (cnt)
		fprintf(fp, "\n");
	if (!label)
	{
		fclose(fp);
		return;
	}
	else
	{
		for (int i = 0; i < chg_aug.size(); i++)
		{
			int aug_cnt = 0;
			fprintf(fp, "augmentation occupancies   %d %d\n", i + 1, chg_aug[i].size());
			for (int j = 0; j < chg_aug[i].size(); j++)
			{
				aug_cnt++;
				fprintf(fp, "% E ", chg_aug[i][j]);
				if (aug_cnt == 5)
				{
					fprintf(fp, "\n"); aug_cnt = 0;
				}
			}
			fprintf(fp, "\n");
		}
		//mag
		if (Separator == NULL)
		{
			fclose(fp);
			return;
		}
		else
		{
			fprintf(fp, "   %d   %d   %d\n", mag_size[0], mag_size[1], mag_size[2]);
			for (int i = 0; i < mag.size(); i++)
			{
				cnt++;
				fprintf(fp, "% E ", mag[i]);
				if (cnt == 5)
				{
					fprintf(fp, "\n");
					cnt = 0;
				}
			}
			if (cnt)
				fprintf(fp, "\n");
			for (int i = 0; i < mag_aug.size(); i++)
			{
				int aug_cnt = 0;
				fprintf(fp, "augmentation occupancies   %d %d\n", i + 1, mag_aug[i].size());
				for (int j = 0; j < mag_aug[i].size(); j++)
				{
					aug_cnt++;
					fprintf(fp, "% E ", mag_aug[i][j]);
					if (aug_cnt == 5)
					{
						fprintf(fp, "\n"); aug_cnt = 0;
					}
				}
				fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
}

void CHGCAR::writespincar(const char file[])
{
	if (GetInfoINCAR("ISPIN") != "2")
		return;
	FILE* fp = fopen(file, "w");
	savposcar(fp, this->pos);
	fprintf(fp, "\n");
	fprintf(fp, "   %d   %d   %d\n", this->mag_size[0], this->mag_size[1], this->mag_size[2]);
	int cnt = 0;
	for (int i = 0; i < this->mag.size(); i++)
	{
		cnt++;
		fprintf(fp, " % E ", this->mag[i]);
		if (cnt == 5)
		{
			fprintf(fp, "\n");
			cnt = 0;
		}
	}
	if (cnt)
		fprintf(fp, "\n");
	fclose(fp);
	printf("Written %s file!\n", file);
}

void CHGCAR::writespinUp_Dwcar(const char up_file[], const char dw_file[])
{
	if (GetInfoINCAR("ISPIN") != "2")
		return;
	vector<double> chg_up(this->chg.size()), chg_dw(this->chg.size());
	for (int i = 0; i < this->chg.size(); i++)
	{
		chg_up[i] = (this->chg[i] + this->mag[i]) / 2;
		chg_dw[i] = (this->chg[i] - this->mag[i]) / 2;
	}
	FILE* fp1 = fopen(up_file, "w");
	FILE* fp2 = fopen(dw_file, "w");
	savposcar(fp1, this->pos);
	fprintf(fp1, "\n");
	savposcar(fp2, this->pos);
	fprintf(fp2, "\n");
	fprintf(fp1, "   %d   %d   %d\n", this->chg_size[0], this->chg_size[1], this->chg_size[2]);
	fprintf(fp2, "   %d   %d   %d\n", this->chg_size[0], this->chg_size[1], this->chg_size[2]);
	int cnt = 0;
	for (int i = 0; i < this->chg.size(); i++)
	{
		cnt++;
		fprintf(fp1, " % E ", chg_up[i]);
		fprintf(fp2, " % E ", chg_dw[i]);
		if (cnt == 5)
		{
			fprintf(fp1, "\n");
			fprintf(fp2, "\n");
			cnt = 0;
		}
	}
	if (cnt)
	{
		fprintf(fp1, "\n");
		fprintf(fp2, "\n");
	}
	fclose(fp1);
	fclose(fp2);
	printf("Written %s file!\n", up_file);
	printf("Written %s file!\n", dw_file);
}

void CHGCAR::ChgcarToCube(const char file1[], const char file2[])
{
	readchgcar(file1);
	FILE* fp = fopen(file2, "w");
	fprintf(fp, "Cube File generated by VASPMATE\n");
	fprintf(fp, "Total:  %d grid points\n", this->chg_size[0] * this->chg_size[1] * this->chg_size[2]);
	fprintf(fp, "  %3d  %10lf  %10lf  %10lf\n", this->pos.nant[0], 0, 0, 0);
	for (int i = 0; i < 3; i++)
		fprintf(fp, "  %3d  %10lf  %10lf  %10lf\n", this->chg_size[i], this->pos.vec[i][0] * this->pos.latt / b2a / this->chg_size[0],
			this->pos.vec[i][1] * this->pos.latt / b2a / this->chg_size[1], this->pos.vec[i][2] * this->pos.latt / b2a / this->chg_size[2]);
	direct_to_carts(this->pos.iflg, this->pos.vec, this->pos.xyz, this->pos.nant);
	int cnt = 0;
	for (int i = 0; i < this->pos.nant[1]; i++)
	{
		for (int j = 0; j < this->pos.typenum[i]; j++)
			fprintf(fp, "  %3d  %10lf  %10lf  %10lf  %lf\n", this->pos.elemnum[i], 0.0, this->pos.xyz[cnt + j][0] / b2a,
				this->pos.xyz[cnt + j][1] / b2a, this->pos.xyz[cnt + j][2] / b2a);
		cnt += this->pos.typenum[i];
	}
	double vol = volume(this->pos.vec);
	cnt = 0;
	for (int i = 0; i < chg_size[0]; i++)
		for (int j = 0; j < chg_size[1]; j++)
			for (int k = 0; k < chg_size[2]; k++)
			{
				fprintf(fp, " % E ", this->Chg[i][j][k] / au2eV / vol);
				cnt++;
				if (cnt == 6)
				{
					fprintf(fp, "\n");
					cnt = 0;
				}
			}
	fclose(fp);
	printf("Written %s file!\n", file2);
}

Cube::Cube(const char file[])
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return;
	}
	char buf[1024];
	int line = -1, cnt = 0;
	nmo = 0;
	while (fgets(buf, 1024, fp))
	{
		line++;
		if (line <= 1)
			continue;
		if (line == 2)
		{
			sscanf(buf, "%d%lf%lf%lf", &atom_num, &org[0], &org[1], &org[2]);
			if (atom_num < 0)
				nmo = 1;
			atom_num = abs(atom_num);
			elemnum.resize(atom_num);
			xyz.resize(atom_num);
			for (int i = 0; i < atom_num; i++)
				xyz[i].resize(3);
			charge.resize(atom_num);
			continue;
		}
		if (line == 3)
		{
			sscanf(buf, "%d%lf%lf%lf", &cub_size[0], &vec[0][0], &vec[0][1], &vec[0][2]);
			continue;
		}
		if (line == 4)
		{
			sscanf(buf, "%d%lf%lf%lf", &cub_size[1], &vec[1][0], &vec[1][1], &vec[1][2]);
			continue;
		}
		if (line == 5)
		{
			sscanf(buf, "%d%lf%lf%lf", &cub_size[2], &vec[2][0], &vec[2][1], &vec[2][2]);
			continue;
		}
		if (line <= 5 + atom_num)
		{
			sscanf(buf, "%d%lf%lf%lf%lf", &elemnum[cnt], &charge[cnt], &xyz[cnt][0], &xyz[cnt][1], &xyz[cnt][2]);
			cnt++;
			if (cnt == atom_num)
				break;
			continue;
		}
	}
	if (nmo == 1) //molecular orbit
	{
		fscanf(fp, "%d", &nmo);
		orbit_index.resize(nmo);
		for (int i = 0; i < nmo; i++)
			fscanf(fp, "%d", &orbit_index[i]);
	}
	for (int i = 0; i < this->elemnum.size(); i++)
	{
		int flag = 0;
		for (int j = 0; j < typenum.size(); j++)
			if (elemnum[i] == typenum[j][0])
			{
				typenum[j].push_back(elemnum[i]);
				flag = 1;
			}
		if (!flag)
			typenum.push_back(vector<int>{elemnum[i]});
	}
	cub.resize(max(1, nmo));
	for (int n = 0; n < max(1, nmo); n++)
	{
		cub[n].resize(cub_size[0]);
		for (int i = 0; i < cub_size[0]; i++) //x
		{
			cub[n][i].resize(cub_size[1]);
			for (int j = 0; j < cub_size[1]; j++) //y
			{
				cub[n][i][j].resize(cub_size[2]);
				for (int k = 0; k < cub_size[2]; k++) //z
					fscanf(fp, "%lf", &cub[n][i][j][k]);
			}
		}
	}
	fclose(fp);
}

void Cube::CubeToChgcar(const char file[], int orbit_num)
{
	if (orbit_num != 0 && nmo == 0)
	{
		printf("No molecular orbital exists in the current file!\n");
		return;
	}
	else if (orbit_num != 0)
	{
		vector<int> ::iterator it = find(orbit_index.begin(), orbit_index.end(), orbit_num);
		if (it == orbit_index.end())
		{
			printf("The index of orbit is error!Please select again\n");
			printf("The index of orbit in the current file is\n");
			for (int i = 0; i < nmo; i++)
				printf("%d ", orbit_index[i]);
			printf("\n");
		}
		else
			orbit_num = it - orbit_index.begin();
	}
	POSCAR pos;
	strcpy(pos.title, "CHGCAR generated by VASPMATE");
	pos.latt = 1;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			pos.vec[i][j] = this->vec[i][j] * b2a * this->cub_size[i] - this->org[i];
	pos.ifix = 1;
	pos.iflg = 1;
	pos.nant[0] = this->elemnum.size();
	pos.nant[1] = this->typenum.size();
	for (int i = 0; i < pos.nant[1]; i++)
		pos.typenum[i] = this->typenum[i].size();
	for (int i = 0; i < pos.nant[0]; i++)
		for (int j = 0; j < 3; j++)
			pos.xyz[i][j] = this->xyz[i][j];
	for (int i = 0; i < pos.nant[1]; i++)
	{
		pos.elemnum[i] = this->typenum[i][0];
		getAtomSym(pos.elemnum[i], pos.elemsym[i]);
	}
	double vol = volume(pos.vec);
	FILE* fp = fopen(file, "w");
	savposcar(fp, pos);
	fprintf(fp,"\n");
	fprintf(fp, "   %d   %d   %d\n", this->cub_size[0], this->cub_size[1], this->cub_size[2]);
	int cnt = 0;
	for (int i = 0; i < cub_size[2]; i++)
		for (int j = 0; j < cub_size[1]; j++)
			for (int k = 0; k < cub_size[0]; k++)
			{
				fprintf(fp, " % E ", this->cub[orbit_num][k][j][i] * au2eV * vol);
				cnt++;
				if (cnt == 5)
				{
					fprintf(fp, "\n");
					cnt = 0;
				}
			}
	fclose(fp);
	printf("Written %s file!\n", file);
}

void Chgcarinterpolation(const char file1[], const char file2[],int scale)
{
	CHGCAR chgcar;
	chgcar.readchgcar(file1);
	POSCAR pos = chgcar.getchgcarpos();
	int init_size[3];
	chgcar.getchg_size(init_size);
	vector<double> init_pot = chgcar.getchg();
	int out_size[3];
	vector<double> out_pot;
	ITP_french(init_size, init_pot, out_size, out_pot, scale);
	writechgcar(file2,pos, out_size, out_pot);
}
