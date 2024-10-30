#include"../../include/VASPMATE_include/elastic.h"
using std::cout;
using std::endl;
using std::vector;
using std::string;
using namespace elas;
static void TranStrainToMatrix(vector<double> strain, vector<vector<double> >& matrix)
{
	matrix.resize(3, vector<double>(3));
	matrix[0][0] = strain[0] + 1;
	matrix[1][1] = strain[1] + 1;
	matrix[2][2] = strain[2] + 1;
	matrix[1][2] = matrix[2][1] = strain[3] / 2;
	matrix[0][2] = matrix[2][0] = strain[4] / 2;
	matrix[0][1] = matrix[1][0] = strain[5] / 2;
}
static void cal_VHR(Eigen::Matrix<double, 6, 6> Sti_matrix,vector<vector<double> > *VHR)
{
	//%%%%%%%%%Voigt
	double Pv = (*VHR)[0][0] = Sti_matrix(0,0) + Sti_matrix(1,1) + Sti_matrix(2,2);
	double Qv = (*VHR)[0][1] = Sti_matrix(0,1) + Sti_matrix(0,2) + Sti_matrix(1,2);
	double Rv = (*VHR)[0][2] = Sti_matrix(3,3) + Sti_matrix(4,4) + Sti_matrix(5,5);

	double Ev = (*VHR)[0][3] = ((Pv+2.00*Qv)*(Pv-Qv+3.00*Rv))/(3.00*(2.00*Pv+3.00*Qv+Rv));
	double Gv = (*VHR)[0][4] = (Pv-Qv+3.00*Rv)/15.00;
	double Kv = (*VHR)[0][5] = Ev*Gv/(3.00*(3.00*Gv-Ev));
    double Muv = (*VHR)[0][6]=((*VHR)[0][3]/(2.00*(*VHR)[0][4]))-1.00;
	
	//%%%%%%%%%Reuss
	Eigen::Matrix<double, 6, 6> Fle_matrix(Sti_matrix.inverse());
	double Pr = (*VHR)[1][0] = Fle_matrix(0,0) + Fle_matrix(1,1) + Fle_matrix(2,2);
	double Qr = (*VHR)[1][1] = Fle_matrix(0,1) + Fle_matrix(0,2) + Fle_matrix(1,2);
	double Rr = (*VHR)[1][2] = Fle_matrix(3,3) + Fle_matrix(4,4) + Fle_matrix(5,5);

	double Er = (*VHR)[1][3] = 15.00/(3.00*Pr+2.00*Qr+Rr);
	double Gr = (*VHR)[1][4] = 15.00/(4.00*(Pr-Qr)+3.00*Rr);
	double Kr = (*VHR)[1][5] = Er*Gr/(3.00*(3.00*Gr-Er));
    double Mur = (*VHR)[1][6]= (Er/(2.00*Gr))-1.00;

	//%%%%%%%%%Hill
	double Eh = (*VHR)[2][0] = (Ev+Er)/2.00;
	double Gh = (*VHR)[2][1] = (Gv+Gr)/2.00;
	double Kh = (*VHR)[2][2] = (Kv+Kr)/2.00;
    double Muh = (*VHR)[2][3]= (Muv+Mur)/2.00;

	//Elastic Anisotropy Index
	double Ac = (*VHR)[2][4] = (Gv-Gr)/(Gv+Gr);
	double Au = (*VHR)[2][5] = 5*Gv/Gr+Kv/Kr-6;

	//Stable: Solve eigenvalues and determine if they are less than 0
	Eigen::EigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> solver(Sti_matrix);
	Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> eigenvalues = solver.eigenvalues();
	Eigen::Matrix<double, Eigen::Dynamic, 1> realEigenvalues = eigenvalues.real();
	for (int i = 0; i < realEigenvalues.rows(); i++) 
	{
    	if (realEigenvalues(i)< 0) (*VHR)[2][6] = 1;
		else (*VHR)[2][6] = 0;
	}
}
static void PrintElaDt(FILE* fp,Eigen::Matrix<double, 6, 6> Sti_matrix)
{
	fprintf(fp, " Elastic tensor: \n\n");
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			fprintf(fp, "%10lf ", Sti_matrix(i,j));
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n Compliance tensor: \n\n");
	Eigen::Matrix<double, 6, 6> Compliance(Sti_matrix.inverse());
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
			fprintf(fp, "%10lf ", Compliance(i, j));
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n Young's, shear and bulk moduli and Poisson ratio \n\n");
	//ubroutine "VHR" is used to calculation the elasticity modulus according elastic constants using the Voigt-Reuss-Hill approximations
	/*To ensure code cleanliness, design a VHR matrix to store datas:
	Pv Qv Rv Ev Gv Kv Muv
	Pr Qr Rr Er Gr Kr Mur
	Eh Gh Kh Muh Ac Au Stable
	*/
	vector<vector<double> > VHR;
	VHR.resize(3, vector<double>(7));
	cal_VHR(Sti_matrix, &VHR);
	fprintf(fp, "    Voigt approximate:   %10lf   %10lf   %10lf   %10lf\n", VHR[0][3],VHR[0][4],VHR[0][5],VHR[0][6]);
	fprintf(fp, "    Reuss approximate:   %10lf   %10lf   %10lf   %10lf\n", VHR[1][3],VHR[1][4],VHR[1][5],VHR[1][6]);
	fprintf(fp, "    Hill approximate :   %10lf   %10lf   %10lf   %10lf\n", VHR[2][0],VHR[2][1],VHR[2][2],VHR[2][3]);
	fprintf(fp, "\n Pugh ratio (G/K):   %10lf\n", VHR[2][1]/VHR[2][2]);
	fprintf(fp, " Cauchy pressure (Pc=C12-C44):   %10lf\n", Sti_matrix(0,1)-Sti_matrix(3,3));
	fprintf(fp, "\n Elastic Anisotropy Index \n\n");
	fprintf(fp, "    Chung-Buessem Anisotropy Index (Ac=(Gv-Gr)/(Gv+Gr)):   %10lf\n", VHR[2][4]);
	fprintf(fp, "    Universal Elastic Anisotropy Index (Au=5*Gv/Gr+Kv/Kr-6):  %10lf\n", VHR[2][5]);
	fprintf(fp, "\n Elastic Stability Conditions:");
	if( std::fabs(VHR[2][6]) < 1e-10 ) fprintf(fp, " Stable\n");
	else fprintf(fp, " Unstable\n");
}
static void Print_csv(FILE* fp,Eigen::Matrix<double, 6, 6> Sti_matrix)
{
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			fprintf(fp, "%lf", Sti_matrix(i,j));
			if( j < 5)
				fprintf(fp, ",");
		}
		fprintf(fp, "\n");
	}
}
double elastic_en::ela_get_energy_oszi()
{
	double energy_oszi = 0;
	FILE* fp_ = fopen("OSZICAR", "r");
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
vector<vector<int> > deformation_type{
		{1,2,3,4,5,6},
		{-2,1,4,-3,6,-5},
		{3,-5,-1,6,2,-4},
		{-4,-6,5,1,-3,2},
		{5,4,6,-2,-1,-3},
		{-6,3,-2,5,-4,1}
};
vector<vector<int> > choice{
	{1},{1,3},{1,3},{1,3},{1,3},{1,3},{1,3,5},{1,2,3,4,5},{1,2,3,4,5,6}
};
vector<string> Roman{ "I","II","III","IV","V","VI" };
template<typename T>
vector<double> operator*(vector<T> v, double u)
{
	vector<double> ans;
	for (auto val : v)
		ans.push_back(val * u);
	return ans;
}
vector<string> vasp_file{ "INCAR","KPOINTS","POTCAR" };
elastic::elastic(vector<double> strain, const char file[])
{
	this->strain = strain;
	cal_method = { &elastic::cal_cub,&elastic::cal_hexa,&elastic::cal_trig1,&elastic::cal_trig2,
		&elastic::cal_tetra1,&elastic::cal_tetra2,&elastic::cal_ortho,&elastic::cal_mono,&elastic::cal_tric };
	Celas.resize(6, vector<double>(6));
	FILE* fp = fopen(file, "r");
	if (fp == nullptr)
	{
		cout << file << " IS NOT EXIST!" << endl;
		return;
	}
	string touch_RECELL = "touch RECELL";
	string RECELL = "RECELL";
	if (!access(RECELL.c_str(), 0))
		cout << "RECELL is exist! && skip this step!" << endl;
	else
		system(touch_RECELL.c_str());
	FILE* fp2 = fopen("RECELL", "at+");
	EV_recell(file,"RECELL");
	readposcar(fp2, pos);
	fclose(fp2);
	double ntemp_vec[3][3];
	transpose_matrix(pos.vec, ntemp_vec);
	int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
	translate_typenum_type(types, pos.nant[1], pos.typenum);
	SpglibDataset* dataset = spg_get_dataset(ntemp_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
	int spacegroup_number = dataset->spacegroup_number;
	if (spacegroup_number == 1 || spacegroup_number == 2)
		Laue = N; 
	if (spacegroup_number >= 3 && spacegroup_number < 16)
		Laue = M;
	if (spacegroup_number >= 16 && spacegroup_number < 75)
		Laue = O;
	if (spacegroup_number >= 75 && spacegroup_number < 89)
		Laue = T2;
	if (spacegroup_number >= 89 && spacegroup_number < 143)
		Laue = T1;
	if (spacegroup_number >= 143 && spacegroup_number < 149)
		Laue = R2;
	if (spacegroup_number >= 149 && spacegroup_number < 168)
		Laue = R1;
	if (spacegroup_number >= 168 && spacegroup_number < 195)
		Laue = H;
	if (spacegroup_number >= 195 && spacegroup_number < 230)
		Laue = C;
	lists = choice[Laue];
	stress.resize(lists.size(), vector<vector<double> >(strain.size(), vector<double>(6)));
	coeff.resize(lists.size(), vector<double>(6));
	spg_free_dataset(dataset);
	free(types);
}
void elastic::generatefile(const char filename[], vector<vector<double> > matrix)
{
	double n_vec[3][3];
	for (int i = 0; i < 3; i++)
	{
		double p[3] = { pos.vec[i][0],pos.vec[i][1],pos.vec[i][2] };
		for (int j = 0; j < 3; j++)
			n_vec[i][j] = p[0] * matrix[j][0] + p[1] * matrix[j][1] + p[2] * matrix[j][2];
	}
	POSCAR p(pos);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			p.vec[i][j] = n_vec[i][j];
	FILE* fp = fopen(filename, "w");
	savposcar(fp, p);
	fclose(fp);
}
void elastic::getstress(const char file[])
{
	vector<vector<vector<double> > > tmp(lists.size(), vector<vector<double> >(6, vector<double>(strain.size())));
	for (int i = 0; i < lists.size(); i++)
	{
		chdir(Roman[i].c_str());
		FILE* fp = fopen(file, "w");
		fprintf(fp, "	Strain       XX(GPa)       YY(GPa)       ZZ(GPa)       XY(GPa)       YZ(GPa)       ZX(GPa)\n");
		for (int j = 0; j < strain.size(); j++)
		{
			string foldname = "strain_" + to_string(strain[j]);
			chdir(foldname.c_str());
			stress[i][j] = get_stress() * (-0.1);
			fprintf(fp, "	%lf       %lf       %lf       %lf       %lf       %lf       %lf\n", strain[j], stress[i][j][0], stress[i][j][1],
				stress[i][j][2], stress[i][j][5], stress[i][j][3], stress[i][j][4]);
			chdir("..");
		}
		fclose(fp);
		chdir("..");
	}
	for (int i = 0; i < lists.size(); i++)
	{
		for (int j = 0; j < 6; j++)
		{
			for (int k = 0; k < stress[i].size(); k++)
				tmp[i][j][k] = stress[i][k][j];
			Eigen::VectorXd Coefficient(FitterLeastSquareMethod(strain, tmp[i][j], 1));
			coeff[i][j] = Coefficient[1];
		}
	}
}
void elastic::generate()
{
	for (int i = 0; i < lists.size(); i++)
	{
		string mkdir = "mkdir " + Roman[i];
		if (!access(Roman[i].c_str(), 0))
			cout << Roman[i] << " is exist! && skip this step!" << endl;
		else
			system(mkdir.c_str());
		cout << "Created " << Roman[i] << " folder!" << endl;
		chdir(Roman[i].c_str());
		for (auto val : strain)
		{
			vector<vector<double> > matrix;
			string foldname = "strain_" + to_string(val);
			string mkdir1 = "mkdir " + foldname;
			if (!access(foldname.c_str(), 0))
				cout << foldname << " is exist! && skip this step!" << endl;
			else
				system(mkdir1.c_str());
			cout << "Created " << foldname << " folder!" << endl;
			chdir(foldname.c_str());
			TranStrainToMatrix(deformation_type[lists[i] - 1] * val, matrix);
			generatefile("POSCAR", matrix);
			for (string file : vasp_file)
			{
				string source_file = "../../" + file;
				string dest_file = "./" + file;
				if (!access(source_file.c_str(), 0))
					copy_file(dest_file.c_str(), source_file.c_str());
			}
			chdir("..");
		}
		chdir("..");
	}
}
void elastic::calculate()
{
	cal_method[Laue](this);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			if (i > j)
				Celas[i][j] = Celas[j][i];
	PrintElasticProperty();
}
void elastic::PrintElasticProperty(const char file[])
{
	FILE* fp = fopen(file, "w");
	Eigen::Matrix<double, 6, 6> _Celas;
	for (int i = 0; i < 6; i++) 
	{
        for (int j = 0; j < 6; j++) 
		{
            _Celas(i, j) = Celas[i][j];
        }
    }
	PrintElaDt(fp,_Celas);
	fclose(fp);
	FILE* fp_csv = fopen("Elastic_tensor.csv", "w");
	Print_csv(fp_csv,_Celas);
	fclose(fp_csv);	
	cout << "Create " << file << " !" << endl;
	cout << "Create " << "Elastic_tensor.csv" << " !" << endl;
}
void elastic::cal_cub()
{
	//%%%%%			Cubic system
	//%%%%% (3 independent elastic constants)
	//%%%%% c11  c12  c12    0    0    0
	//%%%%% c12  c11  c12    0    0    0
	//%%%%% c12  c12  c11    0    0    0
	//%%%%% 0    0    0		c44   0    0
	//%%%%% 0    0    0		 0   c44   0
	//%%%%% 0    0    0		 0    0   c44
	Eigen::Matrix<double, 3, 3> A;
	//c11 c12 c44
	A << 1, 5, 0,
		2, 4, 0,
		0, 0, 4;
	Eigen::Matrix<double, 3, 1> B;
	B << coeff[0][0], coeff[0][1], coeff[0][3];
	Eigen::Matrix<double, 3, 1> C(A.inverse() * B);

	Celas[0][0] = C[0]; //c11
	Celas[0][1] = C[1]; //c12
	Celas[0][2] = Celas[0][1]; //c13 = c12

	Celas[1][1] = Celas[0][0]; //c22 = c11
	Celas[1][2] = Celas[0][1]; //c23 = c12

	Celas[2][2] = Celas[0][0]; //c33 = c11

	Celas[3][3] = Celas[4][4] = Celas[5][5] = C[2]; //c44 = c55 = c66
}
void elastic::cal_hexa()
{
	//%%%%% Hexagonal system
	//%%%%% (5 independent elastic constants)
	//%%%%% c11  c12  c13    0    0    0
	//%%%%% c12  c11  c13    0    0    0
	//%%%%% c13  c13  c33    0    0    0
	//%%%%% 0    0    0    c44    0    0
	//%%%%% 0    0    0      0  c44    0
	//%%%%% 0    0    0      0    0  c66 = (c11 - c12) / 2
	Eigen::Matrix<double, 5, 5> A;
	//c11 c12 c13 c33 c44
	A << 1, 2, 3, 0, 0,
		2, 1, 3, 0, 0,
		0, 0, 3, 3, 0,
		0, 0, 0, 0, 4,
		3, -5, -1, 0, 0;
	Eigen::Matrix<double, 5, 1> B;
	B << coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], coeff[1][0];
			cout << B[0] << endl;
	Eigen::Matrix<double, 5, 1> C(A.inverse() * B);
	Celas[0][0] = C[0]; //c11
	Celas[0][1] = C[1]; //c12
	Celas[0][2] = C[2]; //c13

	Celas[1][1] = Celas[0][0]; //c22 = c12
	Celas[1][2] = Celas[0][2]; //c23 = c13

	Celas[2][2] = C[3]; //c33

	Celas[3][3] = Celas[4][4] = C[4]; //c44 = c55

	Celas[5][5] = (Celas[0][0] - Celas[0][1]) / 2; //c66
}
void elastic::cal_trig1()
{
	//%%%%%		Rhombohedral I system
	//%%%%% (6 independent elastic constants)
	//%%%%% c11  c12  c13  c14    0    0
	//%%%%% c12  c11  c13 -c14    0    0
	//%%%%% c13  c13  c33    0    0    0
	//%%%%% c14 -c14  0    c44    0    0
	//%%%%% 0    0    0    0     c44  c14
	//%%%%% 0    0    0    0     c14  c66 = (c11 - c12) / 2
	Eigen::Matrix<double, 6, 6> A;
	//c11 c12 c13 c14 c33 c44
	A << 1, 2, 3, 4, 0, 0,
		2, 1, 3, -4, 0, 0,
		0, 0, 3, 0, 3, 0,
		0, 0, 0, -1, 0, 4,
		0, 0, 0, 6, 0, 5,
		3, -5, -1, 6, 0, 0;
	Eigen::Matrix<double, 6, 1> B;
	B << coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], coeff[0][4], coeff[1][0];
	Eigen::Matrix<double, 6, 1> C(A.inverse() * B);
	Celas[0][0] = C[0]; //c11
	Celas[0][1] = C[1]; //c12
	Celas[0][2] = C[2]; //c13
	Celas[0][3] = C[3]; //c14

	Celas[1][1] = Celas[0][0]; //c22 = c11
	Celas[1][2] = Celas[0][2]; //c23 = c13
	Celas[1][3] = -Celas[0][3]; //c24 = -c14

	Celas[2][2] = C[4]; //c33

	Celas[3][3] = C[5]; //c44

	Celas[4][4] = Celas[3][3]; //c55 = c44
	Celas[4][5] = Celas[0][3]; //c56 = c14

	Celas[5][5] = (Celas[0][0] - Celas[0][1]) / 2; //c66
}
void elastic::cal_trig2()
{
	//%%%%%		Rhombohedral II system
	//%%%%% (7 independent elastic constants)
	//%%%%% c11  c12  c13  c14  c15    0
	//%%%%% c12  c11  c13 -c14 -c15    0
	//%%%%% c13  c13  c33    0    0    0
	//%%%%% c14 -c14    0  c44    0 -c15
	//%%%%% c15 -c15    0    0  c44  c14
	//%%%%% 0    0      0 -c15  c14  c66 = (c11 - c12) / 2
	Eigen::Matrix<double, 7, 7> A;
	//c11 c12 c13 c14 c15 c33 c44
	A << 1, 2, 3, 4, 5, 0, 0,
		2, 1, 3, -4, -5, 0, 0,
		0, 0, 3, 0, 0, 3, 0,
		0, 0, 0, -1, -6, 0, 4,
		0, 0, 0, 0, 6, -1, 5,
		3, -5, -1, 6, 2, 0, 0,
		-5, 3, -1, -6, -2, 0, 0;
	Eigen::Matrix<double, 7, 1> B;
	B << coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], coeff[0][4], coeff[1][0], coeff[1][1];
	Eigen::Matrix<double, 7, 1> C(A.inverse() * B);
	Celas[0][0] = C[0]; //c11
	Celas[0][1] = C[1]; //c12
	Celas[0][2] = C[2]; //c13
	Celas[0][3] = C[3]; //c14
	Celas[0][4] = C[4]; //c15

	Celas[1][1] = Celas[0][0]; //c22 = c11
	Celas[1][2] = Celas[0][2]; //c23 = c13
	Celas[1][3] = -Celas[0][3]; //c24 = -c14
	Celas[1][4] = -Celas[0][4]; //c25 = -c15

	Celas[2][2] = C[4]; //c33

	Celas[3][3] = C[5]; //c44
	Celas[3][5] = -Celas[0][4]; //c46 = -c15

	Celas[4][4] = Celas[3][3]; //c55 = c44
	Celas[4][5] = Celas[0][3]; //c56 = c14

	Celas[5][5] = (Celas[0][0] - Celas[0][1]) / 2; //c66
}
void elastic::cal_tetra1()
{
	//%%%%%		Tetragonal I system
	//%%%%% (6 independent elastic constants)
	//%%%%% c11  c12  c13    0    0    0
	//%%%%% c12  c11  c13    0    0    0
	//%%%%% c13  c13  c33    0    0    0
	//%%%%% 0    0    0     c44   0    0
	//%%%%% 0    0    0      0   c44   0
	//%%%%% 0    0    0      0    0   c66
	Eigen::Matrix<double, 6, 6> A;
	A << 1, 2, 3, 0, 0, 0,
		2, 1, 3, 0, 0, 0,
		0, 0, 3, 3, 0, 0,
		0, 0, 0, 0, 4, 0,
		0, 0, 0, 0, 0, 6,
		3, -5, -1, 0, 0, 0;
	Eigen::Matrix<double, 6, 1> B;
	B << coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], coeff[0][5], coeff[1][0];
	Eigen::Matrix<double, 6, 1> C(A.inverse() * B);
	Celas[0][0] = C[0]; //c11
	Celas[0][1] = C[1]; //c12
	Celas[0][2] = C[2]; //c13

	Celas[1][1] = Celas[0][0]; //c22 = c11
	Celas[1][2] = Celas[0][2]; //c23 = c13

	Celas[2][2] = C[3]; //c33

	Celas[3][3] = Celas[4][4] = C[4]; //c44 = c55

	Celas[5][5] = C[5]; //c66
}
void elastic::cal_tetra2()
{
	//%%%%%		Tetragonal II system
	//%%%%% (7 independent elastic constants)
	//%%%%% c11  c12  c13    0    0   c16
	//%%%%% c12  c11  c13    0    0  -c16
	//%%%%% c13  c13  c33    0    0    0
	//%%%%% 0    0    0     c44   0    0
	//%%%%% 0    0    0      0   c44   0
	//%%%%% c16 -c16  0      0    0   c66
	Eigen::Matrix<double, 7, 7> A;
	A << 1, 2, 3, 0, 6, 0, 0,
		2, 1, 3, 0, -6, 0, 0,
		0, 0, 3, 3, 0, 0, 0,
		0, 0, 0, 0, 0, 4, 0,
		0, 0, 0, 0, 0, 0, 6,
		3, -5, -1, 0, -4, 0, 0,
		-5, 3, -1, 0, 4, 0, 0;
	Eigen::Matrix<double, 7, 1> B;
	B << coeff[0][0], coeff[0][1], coeff[0][2], coeff[0][3], coeff[0][5], coeff[1][0], coeff[1][1];
	Eigen::Matrix<double, 7, 1> C(A.inverse() * B);
	Celas[0][0] = C[0]; //c11
	Celas[0][1] = C[1]; //c12
	Celas[0][2] = C[2]; //c13
	Celas[0][5] = C[4]; //c16

	Celas[1][1] = Celas[0][0]; //c22 = c12
	Celas[1][2] = Celas[0][2]; //c23 = c13
	Celas[1][5] = -Celas[0][5]; //c26 = -c16

	Celas[2][2] = C[3]; //c33

	Celas[3][3] = Celas[4][4] = C[5]; //c44 = c55

	Celas[5][5] = C[6]; //c66
}
void elastic::cal_ortho()
{
	//%%%%%		Orthorhombic system
	//%%%%% (9 independent elastic constants)
	//%%%%% c11  c12  c13    0    0    0
	//%%%%%	c12  c22  c23    0    0    0
	//%%%%%	c13  c23  c33    0    0    0
	//%%%%%	0    0    0		c44   0    0
	//%%%%%	0    0    0		 0  c55    0
	//%%%%%	0    0    0		 0    0  c66
	Eigen::Matrix<double, 6, 6> _strain;
	Eigen::Matrix<double, 6, 6> _stress;
	_strain.fill(0);
	_stress.fill(0);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 3; j++)
		{
			_strain(i, j) = deformation_type[lists[j] - 1][i];
			_stress(i, j) = coeff[j][i];
		}
	_strain(3, 3) = _strain(4, 4) = _strain(5, 5) = 1;
	_stress(3, 3) = coeff[0][3] / _strain(3,0);
	_stress(4, 4) = coeff[0][4] / _strain(4,0);
	_stress(5, 5) = coeff[0][5] / _strain(5,0);
	auto ans = _stress * (_strain.inverse());
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			Celas[i][j] = ans(i,j);
}
void elastic::cal_mono()
{
	//[e14 e16 * [c44 c46 = [s14 s16
	//e24 e26]	 c46 c66]    s24 s26]
	Eigen::Matrix<double, 2, 2> tmp1; 
	tmp1 << deformation_type[lists[0] - 1][3], deformation_type[lists[0] - 1][5],
		deformation_type[lists[1] - 1][3], deformation_type[lists[1] - 1][5];
	//c66 c46
	Eigen::Matrix<double, 2, 2> tmp2;
	tmp2 << coeff[0][3], coeff[0][5],
		coeff[1][3], coeff[1][5];
	auto tmp3 = tmp1.inverse() * tmp2;
	double c44 = tmp3(0, 0);
	double c66 = tmp3(1, 1);
	Eigen::Matrix<double, 6, 6> _strain;
	Eigen::Matrix<double, 6, 6> _stress;
	_strain.fill(0);
	_stress.fill(0);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 5; j++)
		{
			_strain(i, j) = deformation_type[lists[j] - 1][i];
			_stress(i, j) = coeff[j][i];
		}
	_strain(5, 5) = 1;
	_stress(5, 3) = c44;
	_stress(5, 5) = c66;
	auto ans = _stress * (_strain.inverse());
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			Celas[i][j] = ans(i, j);
}
void elastic::cal_tric()
{
	Eigen::Matrix<double, 6, 6> _strain;
	Eigen::Matrix<double, 6, 6> _stress;
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
		{
			_strain(i, j) = deformation_type[lists[j] - 1][i];
			_stress(i, j) = coeff[j][i];
		}
	auto ans = _stress * (_strain.inverse());
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			Celas[i][j] = ans(i, j);
}


//Calculate using the energy-strain
vector<string> elsa_number{ "01","02","03","04","05","06","07","08","09","10","11",
							"12","13","14","15","16","17","18","19","20","21"};//Like Roman
elastic_en::elastic_en(vector<double> strain_energy, const char file[])
{
	this->strain_energy = strain_energy;
	cal_method_energy = { &elastic_en::cal_cub,&elastic_en::cal_hexa,&elastic_en::cal_trig1,&elastic_en::cal_trig2,
		&elastic_en::cal_tetra1,&elastic_en::cal_tetra2,&elastic_en::cal_ortho,&elastic_en::cal_mono,&elastic_en::cal_tric};
	Celas.resize(6, vector<double>(6));
	FILE* fp = fopen(file, "r");
	if (fp == nullptr)
	{
		cout << file << " IS NOT EXIST!" << endl;
		fclose(fp);
		return;
	}
	string touch_RECELL = "touch RECELL";
	string RECELL = "RECELL";
	if (!access(RECELL.c_str(), 0))
		cout << "RECELL is exist! && skip this step!" << endl;
	else
		system(touch_RECELL.c_str());
	FILE* fp2 = fopen("RECELL", "at+");
	EV_recell(file,"RECELL");
	readposcar(fp2, pos);
	fclose(fp);
	fclose(fp2);
	vol = volume(pos.vec);
	double ntemp_vec[3][3];
	transpose_matrix(pos.vec, ntemp_vec);
	int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
	translate_typenum_type(types, pos.nant[1], pos.typenum);
	SpglibDataset* dataset = spg_get_dataset(ntemp_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
	int spacegroup_number = dataset->spacegroup_number;
	if (spacegroup_number == 1 || spacegroup_number == 2)
		{nelastic = 21; Laue = N;}
	if (spacegroup_number >= 3 && spacegroup_number < 16)
		{nelastic = 13; Laue = M;}
	if (spacegroup_number >= 16 && spacegroup_number < 75)
		{nelastic = 9; Laue = O;}
	if (spacegroup_number >= 75 && spacegroup_number < 89)
		{nelastic = 7; Laue = T2;}
	if (spacegroup_number >= 89 && spacegroup_number < 143)
		{nelastic = 6; Laue = T1;}
	if (spacegroup_number >= 143 && spacegroup_number < 149)
		{nelastic = 7; Laue = R2;}
	if (spacegroup_number >= 149 && spacegroup_number < 168)
		{nelastic = 6; Laue = R1;}
	if (spacegroup_number >= 168 && spacegroup_number < 195)
		{nelastic = 5; Laue = H;}
	if (spacegroup_number >= 195 && spacegroup_number < 230)
		{nelastic = 3; Laue = C;}
	//lists.size()= nelastic;
	defMat.resize(nelastic, vector<double>(6));
	defVect.resize(6);
	conf_e2s.resize(nelastic);
	engxy.resize(nelastic, vector<double>(strain_energy.size()));
	dataxy.resize(nelastic, vector<double>(strain_energy.size()));
	spg_free_dataset(dataset);
	free(types);
}
void elastic_en::generatefile(const char filename[], vector<vector<double> > matrix)
{
	double n_vec[3][3];
	for (int i = 0; i < 3; i++)
	{
		double p[3] = { pos.vec[i][0],pos.vec[i][1],pos.vec[i][2] };
		for (int j = 0; j < 3; j++)
			n_vec[i][j] = p[0] * matrix[j][0] + p[1] * matrix[j][1] + p[2] * matrix[j][2];
	}
	POSCAR p(pos);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			p.vec[i][j] = n_vec[i][j];
	FILE* fp = fopen(filename, "w");
	savposcar(fp, p);
	fclose(fp);
}
void elastic_en::getenergy(const char file[])
{
	FILE* fp = fopen(file, "w");
	fprintf(fp, "Strain ");
	for (int s = 0; s < strain_energy.size(); s++)
		fprintf(fp, "%8.4lf ", strain_energy[s]);
	double E_min = 0 ;
	bool not_none_empty = true;
	for (int i = 0; i < nelastic; i++)
	{
		fprintf(fp, "\n");
		chdir(elsa_number[i].c_str());
		fprintf(fp, "%6s " , elsa_number[i].c_str());
		for (int j = 0; j < strain_energy.size(); j++)
		{
			engxy[i][j] = 0;
			string foldname = "strain_" + to_string(strain_energy[j]);
			chdir(foldname.c_str());
			FILE* fp_ = fopen("OSZICAR", "r");
			if (fp_ == NULL)
			{
				printf("ERROR!!! OSZICAR in %d %s IS NOT EXIST!\n", i, foldname.c_str());
				not_none_empty = false;
			}
			else
			{
				std::ifstream file("OSZICAR", std::ios::binary | std::ios::ate);
				std::streamsize fileSize = file.tellg(); // Get the size of the file
    			file.close();
				if(fileSize == 0)
				{
					printf("ERROR!!! OSZICAR in %d %s IS EMPTY!\n", i, foldname.c_str());
					not_none_empty = false;
				}
			}
			engxy[i][j] = ela_get_energy_oszi();
			if(engxy[i][j] == 0)
				if (not_none_empty == true)
					printf("WARNING!!! ENERGY in %d %s IS ZERO! Please check the OSZICAR.\n", i, foldname.c_str());
			if (engxy[i][j] < E_min)
            	E_min = engxy[i][j];
			fprintf(fp, "%8.4lf " , engxy[i][j]);
			chdir("..");
		}
		chdir("..");
	}
	for (int i = 0; i < nelastic; i++)
	{
		int k = 0;
		vector<double> data_x;
		vector<double> data_y;
		data_x.resize(strain_energy.size());
		data_y.resize(strain_energy.size());
		for (int j = 0; j < strain_energy.size(); j++)
		{
			dataxy[i][j] = (engxy[i][j] - E_min) / vol * 160.2000;
			//To delete the bad points
			if (dataxy[i][j] >= 0.0 && dataxy[i][j] <= 2.0) 
			{
				data_y[k] = dataxy[i][j];
        		data_x[k] = strain_energy[j];
				k++;		
    		}
		}
		//polynomial fitting according the least square method
		Eigen::VectorXd Coefficient(FitterLeastSquareMethod(data_x, data_y, 2));	 
		conf_e2s[i] = Coefficient[2];
	}
}
void elastic_en::calculate()
{
	cal_method_energy[Laue](this);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			if (i > j)
				Celas[i][j] = Celas[j][i];
	PrintElasticProperty();
}
void elastic_en::PrintElasticProperty(const char file[])
{
	FILE* fp = fopen(file, "w");
	Eigen::Matrix<double, 6, 6> _Celas;
	for (int i = 0; i < 6; i++) 
	{
        for (int j = 0; j < 6; j++) 
		{
            _Celas(i, j) = Celas[i][j];
        }
    }
	PrintElaDt(fp,_Celas);
	fclose(fp);
	FILE* fp_csv = fopen("Elastic_tensor.csv", "w");
	Print_csv(fp_csv,_Celas);
	fclose(fp_csv);	
	cout << "Create " << file << " !" << endl;
	cout << "Create " << "Elastic_tensor.csv" << " !" << endl;
}
void elastic_en::generate()
{
	for (int i = 0; i < nelastic; i++)
	{
		string mkdir = "mkdir " + elsa_number[i];
		if (!access(elsa_number[i].c_str(), 0))
			cout << elsa_number[i] << " is exist! && skip this step!" << endl;
		else
			system(mkdir.c_str());
		cout << "Created " << elsa_number[i] << " folder!" << endl;
		chdir(elsa_number[i].c_str());
		for (auto val : strain_energy)
		{
			string foldname = "strain_" + to_string(val);
			string mkdir1 = "mkdir " + foldname;
			if (!access(foldname.c_str(), 0))
				cout << foldname << " is exist! && skip this step!" << endl;
			else
				system(mkdir1.c_str());
			cout << "Created " << foldname << " folder!" << endl;
			chdir(foldname.c_str());
			switch (Laue)
			{
				case C:{matrix_cub(i, val);break;}
				case H:{matrix_hexa(i, val);break;}
				case R1:{matrix_trig1(i, val);break;}
				case R2:{matrix_trig2(i, val);break;}
				case T1:{matrix_tetra1(i, val);break;}
				case T2:{matrix_tetra2(i, val);break;}
				case O:{matrix_ortho(i, val);break;}
				case M:{matrix_mono(i, val);break;}
				case N:{matrix_tric(i, val);break;}
			}
			for (string file : vasp_file)
			{
				string source_file = "../../" + file;
				string dest_file = "./" + file;
				if (!access(source_file.c_str(), 0))
					copy_file(dest_file.c_str(), source_file.c_str());
			}
			chdir("..");
		}
		chdir("..");
	}
}
void elastic_en::matrix_cub(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{0.0, 0.0, 0.0, val, val, val},
              {val, val, 0.0, 0.0, 0.0, 0.0},
              {val, val, val, 0.0, 0.0, 0.0}};
	for (int k = 0; k < 6; k++) 
		defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_hexa(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0, 0.0, 0.0, val},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, val, 0.0},
			  {val, val, val, 0.0, 0.0, 0.0}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_trig1(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0, 0.0, 0.0, val},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, val, 0.0},
			  {val, val, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, 0.0, val, val}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_trig2(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0, 0.0, 0.0, val},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, val, 0.0},
			  {val, val, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, 0.0, val, val},
			  {0.0, val, 0.0, 0.0, 0.0, val}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_tetra1(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0, 0.0, 0.0, val},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, val, 0.0},
			  {val, val, val, 0.0, 0.0, 0.0},
			  {0.0, val, val, 0.0, 0.0, 0.0}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_tetra2(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, 0.0, 0.0, 0.0, val},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, val, 0.0},
			  {val, val, val, 0.0, 0.0, 0.0},
			  {0.0, val, val, 0.0, 0.0, 0.0},
			  {val, 0.0, 0.0, 0.0, 0.0, val}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_ortho(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, 0.0, 0.0, 0.0, 0.0, 0.0},
              {0.0, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, 0.0, 0.0},
			  {0.0, 0.0, 0.0, 0.0, val, 0.0},
			  {0.0, 0.0, 0.0, 0.0, 0.0, val},
			  {val, val, 0.0, 0.0, 0.0, 0.0},
			  {0.0, val, val, 0.0, 0.0, 0.0},
			  {val, 0.0, val, 0.0, 0.0, 0.0}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_mono(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, 0.0, 0.0, 0.0, 0.0, 0.0},
              {0.0, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, 0.0, 0.0},
			  {0.0, 0.0, 0.0, 0.0, val, 0.0},
			  {0.0, 0.0, 0.0, 0.0, 0.0, val},
			  {val, val, 0.0, 0.0, 0.0, 0.0},
			  {0.0, val, val, 0.0, 0.0, 0.0},
			  {val, 0.0, val, 0.0, 0.0, 0.0},
			  {val, 0.0, 0.0, 0.0, val, 0.0},
			  {0.0, val, 0.0, 0.0, val, 0.0},
			  {0.0, 0.0, val, 0.0, val, 0.0},
			  {0.0, 0.0, 0.0, val, 0.0, val}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_en::matrix_tric(int j, double val)
{
	vector<vector<double> > matrix;
	defMat = {{val, 0.0, 0.0, 0.0, 0.0, 0.0},
              {0.0, val, 0.0, 0.0, 0.0, 0.0},
              {0.0, 0.0, val, 0.0, 0.0, 0.0},
			  {0.0, 0.0, 0.0, val, 0.0, 0.0},
			  {0.0, 0.0, 0.0, 0.0, val, 0.0},
			  {0.0, 0.0, 0.0, 0.0, 0.0, val},
			  {val, val, 0.0, 0.0, 0.0, 0.0},
			  {val, 0.0, val, 0.0, 0.0, 0.0},
			  {val, 0.0, 0.0, val, 0.0, 0.0},
			  {val, 0.0, 0.0, 0.0, val, 0.0},
			  {val, 0.0, 0.0, 0.0, 0.0, val},
			  {0.0, val, val, 0.0, 0.0, 0.0},
			  {0.0, val, 0.0, val, 0.0, 0.0},
			  {0.0, val, 0.0, 0.0, val, 0.0},
			  {0.0, val, 0.0, 0.0, 0.0, val},
			  {0.0, 0.0, val, val, 0.0, 0.0},
			  {0.0, 0.0, val, 0.0, val, 0.0},
			  {0.0, 0.0, val, 0.0, 0.0, val},
			  {0.0, 0.0, 0.0, val, val, 0.0},
			  {0.0, 0.0, 0.0, val, 0.0, val},
			  {0.0, 0.0, 0.0, 0.0, val, val}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}

void elastic_en::cal_cub()
{
    Celas[0][0] = 2.00 * conf_e2s[1] - 2.00 / 3.00 * conf_e2s[2];
    Celas[1][1] = Celas[0][0];
    Celas[2][2] = Celas[0][0];
    Celas[3][3] = 2.00 / 3.00 * conf_e2s[0];
    Celas[4][4] = Celas[3][3];
    Celas[5][5] = Celas[3][3];

    Celas[0][1] = 2.00 / 3.00 * conf_e2s[2] - conf_e2s[1];
    Celas[1][2] = Celas[0][1];
    Celas[0][2] = Celas[0][1];

    Celas[1][0] = Celas[0][1];
    Celas[2][1] = Celas[1][2];
    Celas[3][0] = Celas[0][2];
}
void elastic_en::cal_hexa()
{
    Celas[0][0] = (conf_e2s[0] + 4.00 * conf_e2s[1]) / 2.00;
    Celas[1][1] = Celas[0][0];
    Celas[2][2] = conf_e2s[2] * 2.00;
    Celas[3][3] = conf_e2s[3];
    Celas[4][4] = Celas[3][3];
    Celas[5][5] = 2.00 * conf_e2s[1];

    Celas[0][1] = (conf_e2s[0] - 4.00 * conf_e2s[1]) / 2.00;
    Celas[0][2] = (conf_e2s[4] - conf_e2s[0] - conf_e2s[2]) / 2.00;
    Celas[1][2] = Celas[0][2];

    Celas[1][0] = Celas[0][1];
    Celas[2][1] = Celas[1][2];
    Celas[2][0] = Celas[0][2];
}
void elastic_en::cal_trig1()
{
	//Rhombohedral I system
    Celas[0][0] = (conf_e2s[0] + 4.00 * conf_e2s[1]) / 2.00;
    Celas[1][1] = Celas[0][0];
    Celas[2][2] = conf_e2s[2] * 2.00;
    Celas[3][3] = conf_e2s[3];
    Celas[4][4] = Celas[3][3];
    Celas[5][5] = 2.00 * conf_e2s[1];

    Celas[0][1] = (conf_e2s[0] - 4.00 * conf_e2s[1]) / 2.00;
    Celas[0][2] = (conf_e2s[4] - conf_e2s[0] - conf_e2s[2]) / 2.00;
    Celas[0][3] = conf_e2s[5] - (conf_e2s[3] / 2.00) - conf_e2s[1];
    Celas[1][2] = Celas[0][2];
    Celas[1][3] = -Celas[0][3];
    Celas[4][5] = Celas[0][3];

	for (int i = 1; i < 6; i++) 
	{
    	for (int j = 0; j < i; j++) {
        	Celas[i][j] = Celas[j][i];
    	}
	}
}
void elastic_en::cal_trig2()
{
	//Rhombohedral II system
	Celas[0][0] = (conf_e2s[0] + 4 * conf_e2s[1]) / 2.00;
    Celas[1][1] = Celas[0][0];
    Celas[2][2] = conf_e2s[2] * 2;
    Celas[3][3] = conf_e2s[3];
    Celas[4][4] = Celas[3][3];
    Celas[5][5] = 2.00 * conf_e2s[1];

    Celas[0][1] = (conf_e2s[0] - 4 * conf_e2s[1]) / 2.00;
    Celas[0][2] = (conf_e2s[4] - conf_e2s[0] - conf_e2s[2]) / 2.00;
    Celas[1][3] = Celas[0][1];
    Celas[0][3] = conf_e2s[5] - conf_e2s[3]/2.00 - conf_e2s[1];
    Celas[0][4] = -(conf_e2s[6] - Celas[0][0]/2.00 - Celas[3][3]/2.00);
    Celas[3][5] = -Celas[0][4];
    Celas[1][3] = -Celas[0][3];
    Celas[1][4] = -Celas[0][4];
    Celas[4][5] = Celas[0][3];
	Celas[1][2] = Celas[0][2];

	for (int i = 1; i < 6; i++) 
	{
    	for (int j = 0; j < i; j++) {
        	Celas[i][j] = Celas[j][i];
    	}
	}
}
void elastic_en::cal_tetra1()
{
	//Tetragonal I system
    Celas[0][0] = conf_e2s[0] - (conf_e2s[4] - 2.00 * conf_e2s[5] + conf_e2s[2]);
    Celas[1][1] = Celas[0][0];
    Celas[2][2] = conf_e2s[2] * 2.00;
    Celas[3][3] = conf_e2s[3];
    Celas[4][4] = Celas[3][3];
    Celas[5][5] = conf_e2s[1] * 2.00;

    Celas[0][1] = conf_e2s[4] - 2.00 * conf_e2s[5] + conf_e2s[2];
    Celas[0][2] = (conf_e2s[4] - conf_e2s[2] - conf_e2s[0]) / 2.00;
    Celas[1][2] = Celas[0][2];
}
void elastic_en::cal_tetra2()
{
	//Tetragonal II system
    Celas[0][0] = conf_e2s[0] - (conf_e2s[4] - 2.00 * conf_e2s[5] + conf_e2s[2]);
    Celas[1][1] = Celas[0][0];
    Celas[2][2] = conf_e2s[2] * 2.00;
    Celas[3][3] = conf_e2s[3];
    Celas[4][4] = Celas[3][3];
    Celas[5][5] = conf_e2s[1] * 2.00;

    Celas[0][1] = conf_e2s[4] - 2.00 * conf_e2s[5] + conf_e2s[2];
    Celas[0][2] = (conf_e2s[4] - conf_e2s[2] - conf_e2s[0]) / 2.00;
    Celas[0][5] = conf_e2s[6] - Celas[0][0] / 2.00 - Celas[5][5] / 2.00;
    Celas[1][2] = Celas[0][2];
    Celas[1][5] = -Celas[0][5];
}
void elastic_en::cal_ortho()
{
	//Orthorhombic system
    Celas[0][0] = 2 * conf_e2s[0];
    Celas[1][1] = 2 * conf_e2s[1];
    Celas[2][2] = 2 * conf_e2s[2];
    Celas[3][3] = 2 * conf_e2s[3];
    Celas[4][4] = 2 * conf_e2s[4];
    Celas[5][5] = 2 * conf_e2s[5];

    Celas[0][1] = conf_e2s[6] - conf_e2s[0] - conf_e2s[1];
    Celas[1][2] = conf_e2s[7] - conf_e2s[1] - conf_e2s[2];
    Celas[0][2] = conf_e2s[8] - conf_e2s[0] - conf_e2s[2];
}
void elastic_en::cal_mono()
{
	//Monoclinic system
    Celas[0][0] = 2 * conf_e2s[0];
    Celas[1][1] = 2 * conf_e2s[1];
    Celas[2][2] = 2 * conf_e2s[2];
    Celas[3][3] = 2 * conf_e2s[3];
    Celas[4][4] = 2 * conf_e2s[4];
    Celas[5][5] = 2 * conf_e2s[5];

    Celas[0][1] = conf_e2s[6] - conf_e2s[0] - conf_e2s[1];
    Celas[1][2] = conf_e2s[7] - conf_e2s[1] - conf_e2s[2];
    Celas[0][2] = conf_e2s[8] - conf_e2s[0] - conf_e2s[2];
    Celas[0][4] = conf_e2s[9] - conf_e2s[0] - conf_e2s[4];
    Celas[1][4] = conf_e2s[10] - conf_e2s[1] - conf_e2s[4];
    Celas[2][4] = conf_e2s[11] - conf_e2s[2] - conf_e2s[4];
    Celas[3][5] = conf_e2s[12] - conf_e2s[3] - conf_e2s[5];

	for (int i = 1; i < 6; i++) 
	{
    	for (int j = 0; j < i; j++) {
        	Celas[i][j] = Celas[j][i];
    	}
	}
}
void elastic_en::cal_tric()
{
	//Triclinic system
	Celas[0][0] = 2 * conf_e2s[0];
	Celas[1][1] = 2 * conf_e2s[1];
	Celas[2][2] = 2 * conf_e2s[2];
	Celas[3][3] = 2 * conf_e2s[3];
	Celas[4][4] = 2 * conf_e2s[4];
	Celas[5][5] = 2 * conf_e2s[5];

	Celas[0][1] = conf_e2s[6] - conf_e2s[0] - conf_e2s[1];
	Celas[0][2] = conf_e2s[7] - conf_e2s[0] - conf_e2s[2];
	Celas[0][3] = conf_e2s[8] - conf_e2s[0] - conf_e2s[3];
	Celas[0][4] = conf_e2s[9] - conf_e2s[0] - conf_e2s[4];
	Celas[0][5] = conf_e2s[10] - conf_e2s[0] - conf_e2s[5];
	Celas[1][2] = conf_e2s[11] - conf_e2s[1] - conf_e2s[2];
	Celas[1][3] = conf_e2s[12] - conf_e2s[1] - conf_e2s[3];
	Celas[1][4] = conf_e2s[13] - conf_e2s[1] - conf_e2s[4];
	Celas[1][5] = conf_e2s[14] - conf_e2s[1] - conf_e2s[5];
	Celas[2][3] = conf_e2s[15] - conf_e2s[2] - conf_e2s[3];
	Celas[2][4] = conf_e2s[16] - conf_e2s[2] - conf_e2s[4];
	Celas[2][5] = conf_e2s[17] - conf_e2s[2] - conf_e2s[5];
	Celas[3][4] = conf_e2s[18] - conf_e2s[3] - conf_e2s[4];
	Celas[3][5] = conf_e2s[19] - conf_e2s[3] - conf_e2s[5];
	Celas[4][5] = conf_e2s[20] - conf_e2s[4] - conf_e2s[5];

	for (int i = 1; i < 6; i++) 
	{
    	for (int j = 0; j < i; j++) {
        	Celas[i][j] = Celas[j][i];
    	}
	}
}