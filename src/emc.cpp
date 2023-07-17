#include"../include/emc.h"
using std::vector;
using std::string;

EMC::EMC(const char file[])
{
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		exit(1);
	}
	char buf[1024];
	fgets(buf, 1024, fp);
	sscanf(buf, "%d", &mode);
	fgets(buf, 1024, fp);
	sscanf(buf, "%c", &ikmesh);
	fgets(buf, 1024, fp);
	sscanf(buf, "%lf", &kspacing);
	fgets(buf, 1024, fp);
	sscanf(buf, "%d", &fit_point);
	fgets(buf, 1024, fp);
	sscanf(buf, "%lf", &cutoff);
	fgets(buf, 1024, fp);
	sscanf(buf, "%d", &task);
	int cnt = task;
	while (cnt--)
	{
		path p;
		fgets(buf, 1024, fp);
		sscanf(buf, "%lf%lf%lf%lf%lf%lf%s", &p.start[0], &p.start[1], &p.start[2], &p.end[0], &p.end[1], &p.end[2], p.direct);
		pa.push_back(p);
	}
	fclose(fp);
	nkpts_emc.resize(task, vector<vector<double> >(fit_point, vector<double>(3)));
	for (int i = 0; i < task; i++)
	{
		double kvector[3] = { 0 };
		double length = 0;
		for (int j = 0; j < 3; j++)
		{
			kvector[j] = pa[i].end[j] - pa[i].start[j];
			length += kvector[j] * kvector[j];
		}
		length = sqrt(length);
		double cos[3] = { 0 };
		double delta[3] = { 0 };
		for (int j = 0; j < 3; j++)
		{
			cos[j] = kvector[j] / length;
			delta[j] = cutoff * cos[j];
		}
		FILE* fp1 = fopen("POSCAR", "r");
		POSCAR pos;
		if (fp1 == NULL)
		{
			printf("POSCAR IS NOT EXIST!\n");
			exit(1);
		}
		readposcar(fp1, pos);//vec : lattice of convention cell
		fclose(fp1);
		double ret[3][3];
		RecMat(pos.vec, ret);
		double ret_inv[3][3];
		brinv(ret, ret_inv);
		double p[3] = { delta[0],delta[1],delta[2] };
		delta[0] = p[0] * ret_inv[0][0] + p[1] * ret_inv[1][0] + p[2] * ret_inv[2][0];
		delta[1] = p[0] * ret_inv[0][1] + p[1] * ret_inv[1][1] + p[2] * ret_inv[2][1];
		delta[2] = p[0] * ret_inv[0][2] + p[1] * ret_inv[1][2] + p[2] * ret_inv[2][2];
		for (int j = 0; j < fit_point; j++)
			for (int k = 0; k < 3; k++)
				nkpts_emc[i][j][k] = pa[i].start[k] + delta[k] * j;
	}
}

void EMC::get_emc_option(int band_index)
{
	if (mode == 1)
		get_emc_kpath("POSCAR");
	else if (mode == 2)
		get_emc(band_index);
}

void EMC::get_emc_kpath(const char file[])
{
	vector<int> weight;
	vector<vector<double> > irr_coor;
	int* imesh = ir_reciprocal_mesh(file, 0, kspacing, ikmesh, NULL, const_cast<char*>("kv"), weight, irr_coor);
	FILE* fp1 = fopen("EMCKPT", "w");
	fprintf(fp1, "%c	%lf	%d	%d	# Parameters to Generate KPOINTS (Don't Edit This Line)\n", ikmesh, kspacing, irr_coor.size(), fit_point, task);
	fprintf(fp1, "	%d\n", irr_coor.size() + fit_point * task);
	fprintf(fp1, "Reciprocal lattice\n");
	for (int i = 0; i < irr_coor.size(); i++)
		fprintf(fp1, "	%lf	%lf	%lf	%d\n", irr_coor[i][0], irr_coor[i][1], irr_coor[i][2], weight[i]);
	for (int i = 0; i < task; i++)
		for (int j = 0; j < fit_point; j++)
			fprintf(fp1, "	%lf	%lf	%lf	0\n", nkpts_emc[i][j][0], nkpts_emc[i][j][1], nkpts_emc[i][j][2]);
	fclose(fp1);
	printf("Written EMCKPT file!\n");
}

void EMC::get_emc(int band_index)
{
	int cbm = 0, vbm = 0;
	DOS dos;
	string ISPIN = GetInfoINCAR("ISPIN");
	if (ISPIN.length() == 0)
		ISPIN = "1";
	string LORBIT = GetInfoINCAR("LORBIT");
	if (LORBIT.length() == 0)
		LORBIT = "0";
	dos.readDoscar(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
	double efermi = dos.fermi_energy();
	FILE* fp = fopen("POSCAR", "r");
	if (fp == NULL)
	{
		printf("POSCAR IS NOT EXIST!\n");
		return;
	}
	POSCAR pos;
	readposcar(fp, pos);
	read_level_emc(efermi, pos.vec, cbm, vbm, band_index);
	fitting(efermi);
	write_emc(cbm, vbm);
}

void EMC::read_level_emc(double efermi, double vec[3][3], int& cbm, int& vbm, int band_index)
{
	EIGENVAL eigen("EIGENVAL", efermi);
	if (band_index != -1)
	{
		vbm = band_index;
		cbm = vbm + 1;
	}
	else
	{
		cbm = eigen.get_cbm();
		vbm = eigen.get_vbm();
	}
	int nkpts = task * fit_point;
	this->ispin = eigen.get_ispin();
	int nband = eigen.get_bands();
	vector<vector<vector<double> > > level(ispin, vector<vector<double> >(nkpts, vector<double>(nband)));
	vector<vector<vector<double> > > eigenvalue = eigen.get_eigenvalue();
	level_emc.resize(2, vector<vector<vector<double> > >(ispin, vector<vector<double> >(task, vector<double>(fit_point)))); //2*ispin*task*fit_points
	kstep_emc.resize(task, vector<double>(fit_point));
	vector<int> weight;
	vector<vector<double> > irr_coor;
	ir_reciprocal_mesh("POSCAR", 0, kspacing, ikmesh, NULL, const_cast<char*>("kv"), weight, irr_coor);
	int ibz = irr_coor.size();
	for (int i = 0; i < ispin; i++)
		for (int j = 0; j < nkpts; j++)
			for (int k = 0; k < nband; k++)
				level[i][j][k] = eigenvalue[i][j + ibz][k];
	vector<double> kline = eigen.getkpointsline(vec);
	vector<double> kstep(nkpts);
	for (int i = 0; i < nkpts; i++)
		kstep[i] = kline[i + ibz] - kline[ibz];
	int cnt = 0;
	for (int i = 0; i < ispin; i++)
		for (int j = 0; j < task; j++)
			for (int k = 0; k < fit_point; k++)
			{
				level_emc[0][i][j][k] = level[i][cnt][cbm - 1];
				level_emc[1][i][j][k] = level[i][cnt][vbm - 1];
				kstep_emc[j][k] = kstep[cnt] - kstep[j * fit_point];
				cnt++;
			}
}

vector<vector<double> > EMC::fitting(double efermi)
{
	mass.resize(task, vector<double>(2));
	r2.resize(task, vector<double>(2));
	for (int i = 0; i < ispin; i++)
		for (int j = 0; j < task; j++)
			for (int k = 0; k < 2; k++)
			{
				vector<double> energy(fit_point), kpath(fit_point);
				for (int l = 0; l < fit_point; l++)
				{
					energy[l] = (level_emc[k][i][j][l] - efermi) * Hartree;
					kpath[l] = kstep_emc[j][l] * Bohr;
				}
				for (int l = 4; l <= fit_point; l++)
				{
					vector<double> vx(kpath.begin(), kpath.begin() + l);
					vector<double> vy(energy.begin(), energy.begin() + l);
					double* coefficient = (double*)malloc(sizeof(double) * l);
					EMatrix(vx, vy, l, 3, coefficient);
					double R2 = cal_R2(vx, vy, coefficient);
					if (R2 > r2[j][k])
					{
						r2[j][k] = R2;
						mass[j][k] = 1 / (2 * coefficient[3]);
					}
					free(coefficient);
				}
			}
	return mass;
}

void EMC::write_emc(int cbm, int vbm)
{
	FILE* fp = fopen("Effective_Mass", "w");
	fprintf(fp, "+ -------------------Summary------------------- +\n");
	fprintf(fp, "Fitting Point: %d\n", fit_point);
	fprintf(fp, "Band Index:            LUMO = %d	HOMO = %d\n", cbm, vbm);
	fprintf(fp, "Effective-Mass (in m0)	Electron (R2)	Hole (R2)\n");
	for (int i = 0; i < task; i++)
		fprintf(fp, "K-Path Index %d: %s	%lf (%lf)	%lf (%lf)\n", i + 1, pa[i].direct, mass[i][0],r2[i][0], mass[i][1], r2[i][1]);
	fprintf(fp, "+ -------------------Summary------------------- +\n");
	fprintf(fp, "\n");
	fclose(fp);
	printf("Written Effective_Mass file!\n");
}
