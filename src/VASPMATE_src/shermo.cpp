#include"shermo.h"

int TranposcarTomsym_element_t(POSCAR& pos, msym_element_t** element)
{
	for (int i = 0; i < pos.nant[0]; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			double tmp = pos.xyz[i][j] - pos.xyz[0][j];
			if (tmp >= 0.5)
				pos.xyz[i][j] -= 1;
			if (tmp <= -0.5)
				pos.xyz[i][j] += 1;
		}
	}
	direct_to_carts(pos.iflg, pos.vec, pos.xyz, pos.nant);
	msym_element_t* atom = (msym_element_t*)malloc(pos.nant[0] * sizeof(msym_element_t));
	memset(atom, 0, pos.nant[0] * sizeof(msym_element_t));
	int cnt = 0;
	for (int i = 0; i < pos.nant[1]; i++)
	{
		for (int j = 0; j < pos.typenum[i]; j++)
		{
			strcpy(atom[cnt + j].name, pos.elemsym[i]);
			for (int k = 0; k < 3; k++)
				atom[cnt + j].v[k] = pos.xyz[cnt + j][k];
		}
		cnt += pos.typenum[i];
	}
	*element = atom;
	return pos.nant[0];
}

init::init(const char poscar[], const char outcar[])
{
	vector<double> elemass = initmass();
	FILE* fp = fopen(poscar, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", poscar);
		return;
	}
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	//get pointgroup
	strcpy(pointgroup, "NULL!");
	msym_element_t* elements = NULL;
	int length = TranposcarTomsym_element_t(pos, &elements);
	msym_context ctx = (msym_context)msymCreateContext();
	msymSetElements(ctx, length, elements);
	free(elements);  elements = NULL;
	msymFindSymmetry(ctx);
	msymGetPointGroupName(ctx, sizeof(char[6]), pointgroup);
	printf("  Molecular Symmetry:");
	if (!strcmp(pointgroup, "NULL"))
		printf("can't find molecular pointgroup!\n");
	else
		printf("%s\n", pointgroup);
	int cnt = 0;
	for (int i = 0; i < pos.nant[1]; i++)
	{
		for (int j = 0; j < pos.typenum[i]; j++)
		{
			mass.push_back(elemass[pos.elemnum[i]]);
			xyz.push_back(vector<double>{pos.xyz[cnt + j][0], pos.xyz[cnt + j][1], pos.xyz[cnt + j][2]});
		}
		cnt += pos.typenum[i];
	}
	// get wavenum
	FILE* fp_ = fopen(outcar, "r");
	if (fp_ == NULL)
	{
		printf("%s IS NOT EXIST!\n", outcar);
		return;
	}
	char buf[1024];
	while (fgets(buf, 1024, fp_) != NULL)
	{
		double temp = 0;
		if (strstr(buf, "f/i=") != NULL && strstr(buf, "cm-1") != NULL)
		{
			sscanf(buf, "%*d%*s%*lf%*s%*lf%*s%lf%*s", &temp);
			temp *= -1;
			tot_wavenum.push_back(temp);
		}
		else if (strstr(buf, "f/i=") == NULL && strstr(buf, "cm-1") != NULL)
		{
			sscanf(buf, "%*d%*s%*s%*lf%*s%*lf%*s%lf%*s", &temp);
			tot_wavenum.push_back(temp);
		}
		else if (strstr(buf, "energy(sigma->0)") != NULL)
		{
			sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%lf", &E);
		}
	}
	fclose(fp_);
}

InertiaValue::InertiaValue(vector<double> mass, vector<vector<double>> xyz)
{
	totmass = 0;
	for (int i = 0; i < mass.size(); i++)
		totmass += mass[i];
	double cenmassx = 0, cenmassy = 0, cenmassz = 0;
	for (int i = 0; i < mass.size(); i++)
	{
		cenmassx += xyz[i][0] * mass[i] / totmass;
		cenmassy += xyz[i][1] * mass[i] / totmass;
		cenmassz += xyz[i][2] * mass[i] / totmass;
	}
	double v11 = 0, v22 = 0, v33 = 0, v12 = 0, v13 = 0, v23 = 0;
	for (int i = 0; i < mass.size(); i++)
	{
		v11 += mass[i] * ((xyz[i][1] - cenmassy) * (xyz[i][1] - cenmassy) + (xyz[i][2] - cenmassz) * (xyz[i][2] - cenmassz)) / b2a / b2a;
		v22 += mass[i] * ((xyz[i][0] - cenmassx) * (xyz[i][0] - cenmassx) + (xyz[i][2] - cenmassz) * (xyz[i][2] - cenmassz)) / b2a / b2a;
		v33 += mass[i] * ((xyz[i][0] - cenmassx) * (xyz[i][0] - cenmassx) + (xyz[i][1] - cenmassy) * (xyz[i][1] - cenmassy)) / b2a / b2a;
		v12 -= mass[i] * (xyz[i][0] - cenmassx) * (xyz[i][1] - cenmassy) / b2a / b2a;
		v13 -= mass[i] * (xyz[i][0] - cenmassx) * (xyz[i][2] - cenmassz) / b2a / b2a;
		v23 -= mass[i] * (xyz[i][1] - cenmassy) * (xyz[i][2] - cenmassz) / b2a / b2a;
	}
	inertmat << v11, v12, v13, v12, v22, v23, v13, v23, v33;
	EigenSolver<MatrixXd> es(inertmat);
	inert = es.eigenvalues();
}

int thermol::getrotsym(char pointgroup[6])
{
	int rotsym = 0;
	if (pointgroup == "C1" || pointgroup == "Ci" || pointgroup == "Cs" || pointgroup == "Civ")
		rotsym = 1;
	else if (pointgroup == "Dih")
		rotsym = 2;
	else if (pointgroup == "Ih")
		rotsym = 60;
	else if (pointgroup[0] == 'S')
		rotsym = (pointgroup[1] - -'0') / 2;
	else if (pointgroup == "T" || pointgroup == "Td" || pointgroup == "Th")
		rotsym = 12;
	else if (pointgroup == "Oh")
		rotsym = 24;
	else if (pointgroup[0] == 'C')
		rotsym = pointgroup[1] - '0';
	else if (pointgroup[0] == 'D')
		rotsym = (pointgroup[1] - '0') * 2;
	else
	{
		printf(" Warning: Cannot identify point group! Assume rotational symmetry number to be 1!\n");
		rotsym = 1;
	}
	return rotsym;
}

void getthermol_data(thermol& th, vector<double> tot_wavenum, double E)
{
	FILE* fp = fopen("Thermol_Info.dat", "at+");
	//get private info
	vector<thermol> separate(2); //trans rot
	double p_pa = th.getp() * atm2Pa;
	double T = th.getT();
	double totmass = th.gettotmass();
	int imode = th.getimode();
	int ilinear = th.getilinear();
	VectorXcd inert = th.getinert();
	int rotsym = th.getrotsym();
	int ilowfreq = th.getilowfreq();
	double ravib = th.getravib();
	int nSpinmulti = th.getnSpinmulti();
	int sclZPE = th.getsclZPE();
	int sclheat = th.getsclheat();
	int sclS = th.getsclS();
	int sclCV = th.getsclCV();
	vector<double> wavenum;
	if (imode == 0)
	{
		if (ilinear)
			for (int i = 0; i < tot_wavenum.size() - 5; i++)
				wavenum.push_back(tot_wavenum[i]);
		else
			for (int i = 0; i < tot_wavenum.size() - 6; i++)
				wavenum.push_back(tot_wavenum[i]);
	}
	else
	{
		for (int i = 0; i < tot_wavenum.size(); i++)
			wavenum.push_back(tot_wavenum[i]);
	}
	vector<double> freq(wavenum.size());
	for (int i = 0; i < wavenum.size(); i++)
		freq[i] = wavenum[i] * wave2freq;
	// calculate contribution
	fprintf(fp, " Note: Only for translation, U is different to H, and CV is different to CP\n\n");
	if (imode == 0)
	{
		fprintf(fp, "							======== Translation ======== \n");
		separate[0].Q = pow((2 * pi * (totmass * amu2kg) * kb * T / h / h), (3E0 / 2E0)) * R * T / p_pa;
		separate[0].Cv = 3E0 / 2E0 * R;
		separate[0].Cp = 5E0 / 2E0 * R;
		separate[0].U = 3E0 / 2E0 * R * T / 1000;
		separate[0].S = R * (log(separate[0].Q / NA) + 5E0 / 2E0);
		separate[0].H = 5E0 / 2E0 * R * T / 1000;
		if ((inert[0].real() + inert[1].real() + inert[2].real()) < 1E-10) //single atom
		{
			separate[1].Q = 1;
			separate[1].U = 1;
			separate[1].Cv = 1;
			separate[1].S = 1;
		}
		else
		{
			vector<double>inertkg(3);
			for (int i = 0; i < 3; i++)
				inertkg[i] = inert[i].real() * amu2kg * pow((b2a * 1E-10), 2.0);
			if (ilinear)
			{
				separate[1].Q = 8 * pi * pi * inertkg[2] * kb * T / rotsym / h / h;
				separate[1].U = R * T / 1000;
				separate[1].Cv = R;
				separate[1].S = R * (log(separate[1].Q) + 1);
			}
			else
			{
				separate[1].Q = 8 * pi * pi / rotsym / h / h / h * pow((2 * pi * kb * T), (3E0 / 2E0)) * sqrt(inertkg[0] * inertkg[1] * inertkg[2]);
				separate[1].U = 3E0 * R * T / 2E0 / 1000;
				separate[1].Cv = 3E0 * R / 2E0;
				separate[1].S = R * (log(separate[1].Q) + 3E0 / 2E0);
			}
		}
		fprintf(fp, "  Translational q:	%E	q/NA: %E\n", separate[0].Q, separate[0].Q / NA);
		fprintf(fp, "  Translational U:	%lf kJ/mol	%lf kcal/mol\n", separate[0].U, separate[0].U / cal2J);
		fprintf(fp, "  Translational H:	%lf kJ/mol	%lf kcal/mol\n", separate[0].H, separate[0].H / cal2J);
		fprintf(fp, "  Translational S:	%lf J/mol/K	%lf cal/mol/K  -TS:%lf kcal/mol\n", separate[0].S, separate[0].S / cal2J, -separate[0].S / cal2J / 1000 * T);
		fprintf(fp, "  Translational CV:	%lf J/mol/K	%lf cal/mol/K\n", separate[0].Cv, separate[0].Cv / cal2J);
		fprintf(fp, "  Translational CP:	%lf J/mol/K	%lf cal/mol/K\n", separate[0].Cp, separate[0].Cp / cal2J);
		fprintf(fp, "\n");
		fprintf(fp, "  							========= Rotation ========\n");
		fprintf(fp, "  Rotational q:	%E\n", separate[1].Q);
		fprintf(fp, "  Rotational U:	%lf kJ/mol	%lf kcal/mol	=H\n", separate[1].U, separate[1].U / cal2J);
		fprintf(fp, "  Rotational S:	%lf J/mol/K	%lf cal/mol/K   -TS:	%lf kcal/mol\n", separate[1].S, separate[1].S / cal2J, -separate[1].S / cal2J / 1000 * T);
		fprintf(fp, "  Rotational CV:	%lf J/mol/K	%lf cal/mol/K	=CP\n", separate[1].Cv, separate[1].Cv / cal2J);
	}
	else if (imode == 1)
	{
		fprintf(fp, "  Translation contribution is ignored since imode=1\n");
		fprintf(fp, "\n");
		fprintf(fp, "  Rotation contribution is ignored since imode = 1\n");
		separate[0].Cv = 0;
		separate[0].Cp = 0;
		separate[0].U = 0;
		separate[0].S = 0;
		separate[0].H = 0;
		separate[0].Q = 1;
		separate[1].Cv = 0;
		separate[1].U = 0;
		separate[1].S = 0;
		separate[1].Q = 1;
	}
	//Calculating vibration contribution
	double qvib_v0 = 1, qvib_bot = 1;
	for (int i = 0; i < freq.size(); i++)
	{
		if (freq[i] <= 0) continue;
		double freqtmp = freq[i];
		if (ilowfreq == 1 && wavenum[i] < ravib) freqtmp = ravib * wave2freq;// artificially raises all frequencies
		double tmpv0 = 1 / (1 - exp(-h * freqtmp / (kb * T)));
		double tmpbot = exp(-h * freqtmp / (kb * 2 * T)) / (1 - exp(-h * freqtmp / (kb * T)));
		qvib_v0 *= tmpv0;
		qvib_bot *= tmpbot;
	}
	double U_vib_heat = 0, CV_vib = 0, S_vib = 0, ZPE = 0;
	for (int i = 0; i < freq.size(); i++)
	{
		double tmpZPE = 0, tmpheat = 0, tmpCV = 0, tmpS = 0;
		double freqtrunc, prefac_trunc, term_trunc, prefac, term;
		if (freq[i] <= 0) continue;
		if (ilowfreq == 1)
		{
			freqtrunc = ravib * wave2freq;
			prefac_trunc = h * freqtrunc / (kb * T);
			term_trunc = exp(-h * freqtrunc / (kb * T));
		}
		//ZPE
		tmpZPE = wavenum[i] * sclZPE / 2 / au2cm_1 * au2kJ_mol;
		//Heating contribution to U
		if (T > 0)
		{
			prefac = h * freq[i] * sclheat / (kb * T);
			term = exp(-h * freq[i] * sclheat / (kb * T));
			if (ilowfreq == 1 && wavenum[i] < ravib)
			{
				prefac = prefac_trunc;
				term = term_trunc;
			}
			tmpheat = R * T * prefac * term / (1 - term) / 1000;
		}
		//Cv
		prefac = h * freq[i] * sclCV / (kb * T);
		term = exp(-h * freq[i] * sclCV / (kb * T));
		if (ilowfreq == 1 && wavenum[i] < ravib)
		{
			prefac = prefac_trunc;
			term = term_trunc;
		}
		tmpCV = R * prefac * prefac * term / (1 - term) / (1 - term);
		//S
		prefac = h * freq[i] * sclS / (kb * T);
		term = exp(-h * freq[i] * sclS / (kb * T));
		if (ilowfreq == 1 && wavenum[i] < ravib)
		{
			prefac = prefac_trunc;
			term = term_trunc;
		}
		tmpS = R * (prefac * term / (1 - term) - log(1 - term)); //RRHO
		if (ilowfreq == 2)	//Grimmes interpolation
		{
			double miu = h / (8 * pi * pi * freq[i]);
			double Bav = 1E-44;//kg * m ^ 2
			double miup = miu * Bav / (miu + Bav);
			double Sfree = R * (0.5 + log(sqrt(8 * pi * pi * pi * miup * kb * T / h / h)));
			double wei = 1 / (1 + pow((100 / wavenum[i]), 4));
			tmpS = wei * tmpS + (1 - wei) * Sfree;
		}
		ZPE += tmpZPE;
		U_vib_heat += tmpheat;
		CV_vib += tmpCV;
		S_vib += tmpS;
	}
	double U_vib = U_vib_heat + ZPE;
	fprintf(fp, "\n");
	fprintf(fp, "							======== Vibration ========\n");
	fprintf(fp, "  Vibrational q(V=0):	%lf\n", qvib_v0);
	fprintf(fp, "  Vibrational q(bot):	%lf\n", qvib_bot);
	fprintf(fp, "  Vibrational U(T)-U(0):	%lf kJ/mol	%lf kcal/mol   =H(T)-H(0)\n", U_vib_heat, U_vib_heat / cal2J);
	fprintf(fp, "  Vibrational U:	%lf kJ/mol	%lf kcal/mol    =H\n", U_vib, U_vib / cal2J);
	fprintf(fp, "  Vibrational S: 	%lf J/mol/K	%lf cal/mol/K   -TS:	%lf kcal/mol\n", S_vib, S_vib / cal2J, -S_vib / cal2J / 1000 * T);
	fprintf(fp, "  Vibrational CV:	%lf J/mol/K	%lf cal/mol/K   =CP\n", CV_vib, CV_vib / cal2J);
	fprintf(fp, "  Zero-point energy (ZPE):	%lf kJ/mol	%lf kcal/mol	%lf a.u.\n", ZPE, ZPE / cal2J, ZPE / au2kJ_mol);
	//Calculating electron contribution
	fprintf(fp, "\n");
	double q_ele = nSpinmulti, U_ele = 0, CV_ele = 0, S_ele = R * log(q_ele);
	if (imode == 1)
	{
		q_ele = 0;
		S_ele = 0;
	}
	fprintf(fp, "							========= Electron ========\n");
	fprintf(fp, "  Electronic q:	%lf\n", q_ele);
	fprintf(fp, "  Electronic U:	%lf kJ/mol	%lf kcal/mol    =H\n", U_ele, U_ele / cal2J);
	fprintf(fp, "  Electronic S:	%lf J/mol/K	%lf cal/mol/K   -TS:	%lf kcal/mol\n", S_ele, S_ele / cal2J, -S_ele / cal2J / 1000 * T);
	fprintf(fp, "  Electronic CV:	%lf J/mol/K	%lf cal/mol/K   =CP\n", CV_ele, CV_ele / cal2J);
	//Outputting total result
	fprintf(fp, "\n");
	fprintf(fp, "							=========== Total =========\n");
	fprintf(fp, "  Total q(V=0):	%E\n", separate[0].Q * separate[1].Q * qvib_v0 * q_ele);
	fprintf(fp, "  Total q(bot):	%E\n", separate[0].Q * separate[1].Q * qvib_bot * q_ele);
	fprintf(fp, "  Total q(V=0)/NA:	%E\n", separate[0].Q * separate[1].Q * qvib_v0 * q_ele / NA);
	fprintf(fp, "  Total q(bot)/NA:	%E\n", separate[0].Q * separate[1].Q * qvib_bot * q_ele / NA);
	th.Cv = separate[0].Cv + separate[1].Cv + CV_vib + CV_ele;
	th.Cp = separate[0].Cp + separate[1].Cv + CV_vib + CV_ele;
	th.S = separate[0].S + separate[1].S + S_vib + S_ele;
	fprintf(fp, "  Total CV:	%lf J/mol/K	%lf cal/mol/K\n", th.Cv, th.Cv / cal2J);
	fprintf(fp, "  Total CP:	%lf J/mol/K	%lf cal/mol/K\n", th.Cp, th.Cp / cal2J);
	fprintf(fp, "  Total S:	%lf J/mol/K	%lf cal/mol/K	-TS:	%lf J/mol	%lf ev\n", th.S, th.S / cal2J, -th.S * T, -th.S * T / au2kJ_mol * au2eV / 1000);
	if (T == 0)
	{
		th.U = ZPE;
		th.H = ZPE;
		th.G = ZPE;
	}
	else
	{
		th.U = separate[0].U + separate[1].U + U_vib + U_ele;
		th.H = separate[0].H + separate[1].U + U_vib + U_ele;
		th.G = th.H - T * th.S / 1000;
	}
	fprintf(fp, "  Zero point energy (ZPE):	%lf kJ/mol/K	%lf kcal/mol	%lf a.u.	%lf ev\n", ZPE, ZPE / cal2J, ZPE / au2kJ_mol, ZPE / au2kJ_mol * au2eV);
	fprintf(fp, "  Thermal correction to U:	%lf kJ/mol/K	%lf kcal/mol	%lf a.u.	%lf ev\n", th.U, th.U / cal2J, th.U / au2kJ_mol, th.U / au2kJ_mol * au2eV);
	fprintf(fp, "  Thermal correction to H:	%lf kJ/mol/K	%lf kcal/mol	%lf a.u.	%lf ev\n", th.H, th.H / cal2J, th.H / au2kJ_mol, th.H / au2kJ_mol * au2eV);
	fprintf(fp, "  Thermal correction to G:	%lf kJ/mol/K	%lf kcal/mol	%lf a.u.	%lf ev\n", th.G, th.G / cal2J, th.G / au2kJ_mol, th.G / au2kJ_mol * au2eV);
	double U0 = E + ZPE / au2kJ_mol;
	double U_final = E + th.U / au2kJ_mol;
	double H_final = E + th.H / au2kJ_mol;
	double G_final = E + th.G / au2kJ_mol;
	fprintf(fp, "  Electronic energy:	%lf a.u.	%lf ev\n", E, E * au2eV);
	fprintf(fp, "  Sum of electronic energy and ZPE, namely U/H/G at 0 K:	%lf a.u.	%lf ev\n", U0, U0 * au2eV);
	fprintf(fp, "  Sum of electronic energy and thermal correction to U:	%lf a.u.	%lf ev\n", U_final, U_final * au2eV);
	fprintf(fp, "  Sum of electronic energy and thermal correction to H:	%lf a.u.	%lf ev\n", H_final, H_final * au2eV);
	fprintf(fp, "  Sum of electronic energy and thermal correction to G:	%lf a.u.	%lf ev\n", G_final, G_final * au2eV);
	fprintf(fp, "\n");
	fprintf(fp, "  Thanks to Sobereva!(sobereva@sina.com)\n");
	fprintf(fp, "  More information about Thermol and Gaussian please refer to http://sobereva.com/552 and  https://gaussian.com/thermo/)\n");
	fclose(fp);
	printf("  Zero point energy (ZPE): %lf kcal/mol %lf ev\n", ZPE / cal2J, ZPE / au2kJ_mol * au2eV);
	printf("  Thermal correction to U: %lf kcal/mol %lf ev\n", th.U / cal2J, th.U / au2kJ_mol * au2eV);
	printf("  Thermal correction to H: %lf kcal/mol %lf ev\n", th.H / cal2J, th.H / au2kJ_mol * au2eV);
	printf("  Thermal correction to G: %lf kcal/mol %lf ev\n", th.G / cal2J, th.G / au2kJ_mol * au2eV);
	printf("  Total S                : %lf J/mol/K %lf ev/K\n", th.S, th.S / au2kJ_mol * au2eV / 1000);
	printf("  TS                    : %lf J/mol %lf ev\n", th.S * T, th.S * T / au2kJ_mol * au2eV / 1000);
}

void shermo(int argc, char* argv[])
{
	printf("================Thermo Energy Calculation!================\n");
	int imode = 0, p = 1, ilowfreq = 0, nSpinmulti = 1;
	double T = 298.15, ravib = 100;
	for (int i = 2; i < argc; i++)
	{
		if (!strcmp("-i", argv[i]) || !strcmp("-im", argv[i]))
			imode = atoi(argv[i + 1]);
		if (!strcmp("-t", argv[i]) || !strcmp("-T", argv[i]))
			T = atof(argv[i + 1]);
		if (!strcmp("-p", argv[i]) || !strcmp("-P", argv[i]))
			p = atoi(argv[i + 1]);
		if (!strcmp("-s", argv[i]) || !strcmp("-sm", argv[i]))
			nSpinmulti = atof(argv[i + 1]);
		if (!strcmp("-l", argv[i]) || !strcmp("-lf", argv[i]))
			ilowfreq = atoi(argv[i + 1]);
		if (!strcmp("-v", argv[i]) || !strcmp("-cv", argv[i]))
			ravib = atof(argv[i + 1]);
	}
	printf("  Temperature(K): %lf\n", T);
	printf("  Pressure(Atm): %d\n", p);
	printf("  Spin multiplicity(Number of Unpaired electron + 1): %d\n", nSpinmulti);
	printf("  The treatment of low frequencies: ");
	if (ilowfreq == 0)
		printf("Harmonic\n");
	if (ilowfreq == 1)
	{
		printf("Raising low frequencies.\n");
		printf("  Raising lower frequencies to  %lf (cm^-1)!\n", ravib);
	}
	if (ilowfreq == 2)
		printf("Grimme's entropy interpolation\n");
	FILE* fp = fopen("Thermol_Info.dat", "w+");
	fprintf(fp, "			 ======== Thermo Energy Calculated by VASPMATE ========\n");
	fprintf(fp, "				  =========== Calculate Inparameter =========\n");
	fprintf(fp, "  Temperature(K): %lf\n", T);
	fprintf(fp, "  Pressure(Atm): %d\n", p);
	fprintf(fp, "  Spin multiplicity(Number of Unpaired electron + 1): %d\n", nSpinmulti);
	fprintf(fp, "  The treatment of low frequencies: ");
	if (ilowfreq == 0)
		fprintf(fp, "Harmonic\n");
	if (ilowfreq == 1)
	{
		fprintf(fp, "  Raising low frequencies.\n");
		fprintf(fp, "  Raising lower frequencies to  %lf (cm^-1)!\n", ravib);
	}
	if (ilowfreq == 2)
		fprintf(fp, "Grimme's entropy interpolation\n");
	fclose(fp);
	init _in("POSCAR", "OUTCAR");
	double E = _in.getE() / au2eV;
	thermol th(_in, imode, p, T, ilowfreq, ravib, nSpinmulti);
	getthermol_data(th, _in.gettot_wavenum(), E);
	printf("  More calculate details please read Thermol_Info.dat file!\n");
	printf("===========================================================\n");
}
