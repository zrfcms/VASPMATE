#include"../../include/VASPMATE_include/elastic_DL.h"
using std::cout;
using std::endl;
using std::vector;
using std::string;
using namespace elas_DL;
#define ZERO_TOLERANCE1 1e-4
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
void elastic_DL::Show_cell(double lattice[3][3], double position[][3], const int types[], const int num_atom)
{
	int i;
	FILE* fp = fopen("structure_operator", "at+");
	fprintf(fp, "Lattice parameter:\n");
	for (i = 0; i < 3; i++) {
		fprintf(fp, "%f %f %f\n", lattice[i][0], lattice[i][1], lattice[i][2]);
	}
	fprintf(fp, "Atomic positions:\n");
	for (i = 0; i < num_atom; i++) {
		fprintf(fp, "%d: %f %f %f\n",
			types[i], position[i][0], position[i][1], position[i][2]);
	}
	fclose(fp);
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
	if( std::fabs(VHR[2][6]) < 1e-10 ) fprintf(fp, " Stable");
	else fprintf(fp, " Unstable");
}
template<typename T>
vector<double> operator*(vector<T> v, double u)
{
	vector<double> ans;
	for (auto val : v)
		ans.push_back(val * u);
	return ans;
}
static vector<string> vasp_file{ "INCAR","KPOINTS","POTCAR" };
static vector<string> elsa_number{ "01","02","03","04","05","06","07","08","09","10","11",
							"12","13","14","15","16","17","18","19","20","21"};//Like Roman
elastic_DL::elastic_DL(int numDL , vector<double> strain_energyDL, const char file[])
{
	this->strain_energyDL = strain_energyDL;
	cal_method_energy = { &elastic_DL::cal_cub,&elastic_DL::cal_hexa,&elastic_DL::cal_trig,&elastic_DL::cal_tetra,&elastic_DL::cal_ortho,&elastic_DL::cal_mono};
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
	Recell(numDL,file,"RECELL");
	readposcar(fp2, pos);
	fclose(fp);
	fclose(fp2);
	vol = volume(pos.vec);
	if (numDL==1)
		{nelastic = 15; Laue = R;}
	if (numDL==2)
		{nelastic = 13; Laue = M;}
	if (numDL==3)
		{nelastic = 9; Laue = O;}
	if (numDL==4)
		{nelastic = 11; Laue = T;}
	if (numDL==5)
		{nelastic = 9; Laue = H;}
	if (numDL==6)
		{nelastic = 9; Laue = C;}
	defMat.resize(nelastic, vector<double>(6));
	defVect.resize(6);
	conf_e2s.resize(nelastic);
	engxy.resize(nelastic, vector<double>(strain_energyDL.size()));
	dataxy.resize(nelastic, vector<double>(strain_energyDL.size()));
}
void elastic_DL::generatefile(const char filename[], vector<vector<double> > matrix)
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
void elastic_DL::getenergy(const char file[])
{
	FILE* fp = fopen(file, "w");
	fprintf(fp, "    Strain       ");
	for (int s = 0; s < strain_energyDL.size(); s++)
		fprintf(fp, "%.6lf       ", strain_energyDL[s]);
	double E_min = 0 ;
	for (int i = 0; i < nelastic; i++)
	{
		fprintf(fp, "\n");
		chdir(elsa_number[i].c_str());
		fprintf(fp, "    %s           " , elsa_number[i].c_str());
		for (int j = 0; j < strain_energyDL.size(); j++)
		{
			engxy[i][j] = 0;
			string foldname = "strain_" + to_string(strain_energyDL[j]);
			chdir(foldname.c_str());
			engxy[i][j] = get_energy_oszi();
			if (engxy[i][j] < E_min)
            	E_min = engxy[i][j];
			fprintf(fp, "%.6lf       " , engxy[i][j]);
			chdir("..");
		}
		chdir("..");
	}
	for (int i = 0; i < nelastic; i++)
	{
		int k = 0;
		vector<double> data_x;
		vector<double> data_y;
		data_x.resize(strain_energyDL.size());
		data_y.resize(strain_energyDL.size());
		for (int j = 0; j < strain_energyDL.size(); j++)
		{
			dataxy[i][j] = (engxy[i][j] - E_min) / vol * 160.2000;
			//To delete the bad points
			if (dataxy[i][j] >= 0.0 && dataxy[i][j] <= 2.0) 
			{
				data_y[k] = dataxy[i][j];
        		data_x[k] = strain_energyDL[j];
				k++;		
    		}
		}
		//polynomial fitting according the least square method
		Eigen::VectorXd Coefficient(FitterLeastSquareMethod(data_x, data_y, 2));	 
		conf_e2s[i] = Coefficient[2];
	}
}
void elastic_DL::calculate()
{
	cal_method_energy[Laue](this);
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			if (i > j)
				Celas[i][j] = Celas[j][i];
	PrintElasticProperty();
}
void elastic_DL::PrintElasticProperty(const char file[])
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
	cout << "Create " << file << " !" << endl;
}
void elastic_DL::generate()
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
		for (auto val : strain_energyDL)
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
				case T:{matrix_tetra(i, val);break;}
				case O:{matrix_ortho(i, val);break;}
				case M:{matrix_mono(i, val);break;}
				case R:{matrix_trig(i, val);break;}
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
void elastic_DL::matrix_cub(int j, double val)
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
void elastic_DL::matrix_hexa(int j, double val)
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
void elastic_DL::matrix_trig(int j, double val)
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
			  {0.0, val, val, 0.0, 0.0, 0.0},
			  {0.0, val, 0.0, val, 0.0, 0.0},
			  {0.0, val, 0.0, 0.0, val, 0.0},
			  {0.0, 0.0, 0.0, val, 0.0, val},
			  {0.0, 0.0, 0.0, 0.0, val, val}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_DL::matrix_tetra(int j, double val)
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
			  {val, 0.0, 0.0, 0.0, 0.0, val},
			  {0.0, val, 0.0, 0.0, 0.0, val}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}
void elastic_DL::matrix_ortho(int j, double val)
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
void elastic_DL::matrix_mono(int j, double val)
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
			  {val, 0.0, 0.0, 0.0, 0.0, val},
			  {0.0, val, 0.0, 0.0, 0.0, val},
			  {0.0, 0.0, val, 0.0, 0.0, val},
			  {0.0, 0.0, 0.0, val, val, 0.0}};
	for (int k = 0; k < 6; k++) 
			defVect[k] = defMat[j][k];
	TranStrainToMatrix(defVect, matrix);
	generatefile("POSCAR", matrix);
}

void elastic_DL::cal_cub()
{
    Celas[0][0] = 2.00 * conf_e2s[0];
    Celas[1][1] = 2.00 * conf_e2s[1];
    Celas[2][2] = 2.00 * conf_e2s[2];
    Celas[3][3] = 2.00 * conf_e2s[3];
    Celas[4][4] = 2.00 * conf_e2s[4];
    Celas[5][5] = 2.00 * conf_e2s[5];

    Celas[0][1] = conf_e2s[6] - conf_e2s[0] - conf_e2s[1];
    Celas[1][2] = conf_e2s[7] - conf_e2s[1] - conf_e2s[2];
    Celas[0][2] = conf_e2s[8] - conf_e2s[0] - conf_e2s[2];
}
void elastic_DL::cal_hexa()
{
    Celas[0][0] = 2.00 * conf_e2s[0];
    Celas[1][1] = 2.00 * conf_e2s[1];
    Celas[2][2] = 2.00 * conf_e2s[2];
    Celas[3][3] = 2.00 * conf_e2s[3];
    Celas[4][4] = 2.00 * conf_e2s[4];
    Celas[5][5] = 2.00 * conf_e2s[5];

    Celas[0][1] = conf_e2s[6] - conf_e2s[0] - conf_e2s[1];
    Celas[1][2] = conf_e2s[7] - conf_e2s[1] - conf_e2s[2];
    Celas[0][2] = conf_e2s[8] - conf_e2s[0] - conf_e2s[2];
}
void elastic_DL::cal_trig()
{
	//Rhombohedral system
    Celas[0][0] = 2.00 * conf_e2s[0];
    Celas[1][1] = 2.00 * conf_e2s[1];
    Celas[2][2] = 2.00 * conf_e2s[2];
    Celas[3][3] = 2.00 * conf_e2s[3];
    Celas[4][4] = 2.00 * conf_e2s[4];
    Celas[5][5] = 2.00 * conf_e2s[5];

    Celas[0][1] = conf_e2s[6] - conf_e2s[0] - conf_e2s[1];
	Celas[0][2] = conf_e2s[7] - conf_e2s[0] - conf_e2s[2];
	Celas[0][3] = conf_e2s[8] - conf_e2s[0] - conf_e2s[3];
	Celas[0][4] = conf_e2s[9] - conf_e2s[0] - conf_e2s[4];
    Celas[1][2] = conf_e2s[10] - conf_e2s[1] - conf_e2s[2];
    Celas[1][3] = conf_e2s[11] - conf_e2s[1] - conf_e2s[3];
    Celas[1][4] = conf_e2s[12] - conf_e2s[1] - conf_e2s[4];
    Celas[3][5] = conf_e2s[13] - conf_e2s[3] - conf_e2s[5];
    Celas[4][5] = conf_e2s[14] - conf_e2s[4] - conf_e2s[5];
}
void elastic_DL::cal_tetra()
{
	//Tetragonal system
    Celas[0][0] = 2.00 * conf_e2s[0];
    Celas[1][1] = 2.00 * conf_e2s[1];
    Celas[2][2] = 2.00 * conf_e2s[2];
    Celas[3][3] = 2.00 * conf_e2s[3];
    Celas[4][4] = 2.00 * conf_e2s[4];
    Celas[5][5] = 2.00 * conf_e2s[5];

	Celas[0][1] = conf_e2s[6] - conf_e2s[0] - conf_e2s[1];
	Celas[1][2] = conf_e2s[7] - conf_e2s[1] - conf_e2s[2];
	Celas[0][2] = conf_e2s[8] - conf_e2s[0] - conf_e2s[2];
	Celas[0][5] = conf_e2s[9] - conf_e2s[0] - conf_e2s[5];
	Celas[1][5] = conf_e2s[10] - conf_e2s[1] - conf_e2s[5];
}
void elastic_DL::cal_ortho()
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
void elastic_DL::cal_mono()
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
    Celas[0][5] = conf_e2s[9] - conf_e2s[0] - conf_e2s[5];
    Celas[1][5] = conf_e2s[10] - conf_e2s[1] - conf_e2s[5];
    Celas[2][5] = conf_e2s[11] - conf_e2s[2] - conf_e2s[5];
    Celas[3][4] = conf_e2s[12] - conf_e2s[3] - conf_e2s[4];
}

void elastic_DL::Recell(int numDL , const char file1[], const char file2[])
{
	FILE* _fp = fopen("structure_operator", "at+");
	fprintf(_fp, "****************************************************\n");
	fprintf(_fp, "INPUT: %s\n", file1);
	fclose(_fp);
	FILE* fp1 = fopen(file1, "r");
	if (fp1 == NULL)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(fp_error, "ERROR!The %s is not exist!\n", file1);
		fprintf(fp_error, "****************************************************\n");
		fclose(fp_error);
		return;
	}
	readposcar(fp1, pos.title, pos.latt, pos.ifix, pos.iflg, pos.nant, pos.typenum, pos.elemnum, pos.elemsym, pos.vec, pos.xyz, pos.fix);
	fclose(fp1);
	int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
	translate_typenum_type(types, pos.nant[1], pos.typenum);
	Show_cell(pos.vec, pos.xyz, types, pos.nant[0]);
	FILE* fp_ = fopen("structure_operator", "at+");
	fprintf(fp_, "**********************RECELL***********************\n");
	fprintf(fp_, "FILE NAME: %s\n", file2);
	fclose(fp_);
	
	//get spacegroupnum
	char symbol[20];
	double ntemp_vec[3][3];
	transpose_matrix(pos.vec, ntemp_vec);
	int spacegroup_number = spg_get_international(symbol, ntemp_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
	if (spacegroup_number == 0)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(_fp, "ERROR!Cannot get space group!\n");
		return;
	}
	spacegroup_number == 1;
	//refine cell
	double n_vec[3][3];
	transpose_matrix(pos.vec, n_vec);
	
	/* 4 times larger memory space must be prepared. */
	int* n_types = (int*)malloc(sizeof(int) * pos.nant[0] * 4);
	for (int i = 0; i < pos.nant[0]; i++)
		n_types[i] = types[i];
	double(*n_xyz)[3];
	n_xyz = new double[4 * pos.nant[0]][3];
	for (int i = 0; i < pos.nant[0]; i++)
		for (int j = 0; j < 3; j++)
			n_xyz[i][j] = pos.xyz[i][j];

	/*
	int refine_cell = spg_refine_cell(n_vec, n_xyz, n_types, pos.nant[0], SYMPREAC);
	transpose_matrix(n_vec, pos.vec);
	if (refine_cell == 0)
	{
		FILE* fp_error = fopen("structure_operator", "at+");
		fprintf(_fp, "ERROR!Cannot recell the structure!\n");
		return;
	}
	pos.nant[0] = refine_cell;*/
	int* unit_typenum = (int*)malloc(sizeof(int) * pos.nant[1]);
	for (int i = 0; i < pos.nant[0] - 1; i++)
	{
		for (int j = 0; j < pos.nant[0] - 1 - i; j++)
		{
			if (n_types[j] > n_types[j + 1])
			{
				swap(n_types[j], n_types[j + 1]);
				for (int k = 0; k < 3; k++)
				{
					swap(n_xyz[j][k], n_xyz[j + 1][k]);
				}
			}
		}
	}
	translate_type_typenum(n_types, pos.nant, unit_typenum);

	//_vec, n_xyz, n_types, unit_typenu
	double temp, lega, legb, legc, legar, legbr, legcr, cos_alpha, cos_bata, cos_gamma, cos_alphar,
		cos_batar, cos_gammar, tempvect[5][5], recipvect[3][3], cos_angle, temppos[MAX_NATOM][3];

	//Triclinic lattice : c < a < b
	if (spacegroup_number == 1 || spacegroup_number == 2)
	{
		lega = sqrt(pos.vec[0][0] * pos.vec[0][0] + pos.vec[0][1] * pos.vec[0][1] + pos.vec[0][2] * pos.vec[0][2]);
		legb = sqrt(pos.vec[1][0] * pos.vec[1][0] + pos.vec[1][1] * pos.vec[1][1] + pos.vec[1][2] * pos.vec[1][2]);
		legc = sqrt(pos.vec[2][0] * pos.vec[2][0] + pos.vec[2][1] * pos.vec[2][1] + pos.vec[2][2] * pos.vec[2][2]);
		cos_alpha = (pos.vec[1][0] * pos.vec[2][0] + pos.vec[1][1] * pos.vec[2][1] + pos.vec[1][2] * pos.vec[2][2]) / legb / legc;
		cos_bata = (pos.vec[0][0] * pos.vec[2][0] + pos.vec[0][1] * pos.vec[2][0] + pos.vec[0][2] * pos.vec[2][2]) / lega / legc;
		cos_gamma = (pos.vec[0][0] * pos.vec[1][0] + pos.vec[0][1] * pos.vec[1][1] + pos.vec[0][2] * pos.vec[1][2]) / lega / legb;
		//define the reciprocal matrix of privect
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				tempvect[i][j] = pos.vec[i][j];
		for (int i = 0; i < 3; i++)
			for (int j = 3; j < 5; j++)
				tempvect[i][j] = tempvect[i][j - 3];
		for (int i = 3; i < 5; i++)
			for (int j = 0; j < 5; j++)
				tempvect[i][j] = tempvect[i - 3][j];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				recipvect[i][j] = tempvect[i + 1][j + 1] * tempvect[i + 2][j + 2] - tempvect[i + 1][j + 2] * tempvect[i + 2][j + 1];
		//end of defining the reciprocal matrix of privect
		legar = sqrt(recipvect[0][0] * recipvect[0][0] + recipvect[0][1] * recipvect[0][1] + recipvect[0][2] * recipvect[0][2]);
		legbr = sqrt(recipvect[1][0] * recipvect[1][0] + recipvect[1][1] * recipvect[1][1] + recipvect[1][2] * recipvect[1][2]);
		legcr = sqrt(recipvect[2][0] * recipvect[2][0] + recipvect[2][1] * recipvect[2][1] + recipvect[2][2] * recipvect[2][2]);
		cos_alphar = (pos.vec[0][0] * recipvect[0][0] + pos.vec[0][1] * recipvect[0][1] + pos.vec[0][2] * recipvect[0][2]) / lega / legar;
		cos_batar = (pos.vec[1][0] * recipvect[1][0] + pos.vec[1][1] * recipvect[1][1] + pos.vec[1][2] * recipvect[1][2]) / legb / legbr;
		cos_gammar = (pos.vec[2][0] * recipvect[2][0] + pos.vec[2][1] * recipvect[2][1] + pos.vec[2][2] * recipvect[2][2]) / legc / legcr;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				tempvect[i][j] = 0.0000;
		// a <= b <= c
		if (lega <= legc && legb <= legc && lega <= legb)
		{
			tempvect[2][2] = lega;
			tempvect[0][2] = legb * cos_gamma;
			tempvect[0][0] = sqrt(legb * legb - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legc * cos_bata;
			tempvect[1][1] = legc * cos_gammar;
			tempvect[1][0] = sqrt(legc * legc - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			//To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / legb / legc;
			if (fabs(cos_angle - cos_alpha) >= ZERO_TOLERANCE1)
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < pos.nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][0];
				temppos[j][0] = n_xyz[j][1];
				temppos[j][1] = n_xyz[j][2];
			}
		}
		// a <= c <= b
		if (lega <= legb && legc <= legb && lega <= legc)
		{
			tempvect[2][2] = lega;
			tempvect[0][2] = legc * cos_bata;
			tempvect[0][0] = sqrt(legc * legc - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legb * cos_gamma;
			tempvect[1][1] = legb * cos_batar;
			tempvect[1][0] = sqrt(legb * legb - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			//To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / legb / legc;
			if (fabs(cos_angle - cos_bata) >= ZERO_TOLERANCE1)
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < pos.nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][0];
				temppos[j][0] = n_xyz[j][2];
				temppos[j][1] = n_xyz[j][1];
			}
		}
		// b <= a <= c
		if (legb <= legc && lega <= legc && legb <= lega)
		{
			tempvect[2][2] = legb;
			tempvect[0][2] = lega * cos_gamma;
			tempvect[0][0] = sqrt(lega * lega - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legc * cos_alpha;
			tempvect[1][1] = legb * cos_gamma;
			tempvect[1][0] = sqrt(legc * legc - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legc;
			if (fabs(cos_angle - cos_bata >= ZERO_TOLERANCE1))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < pos.nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][1];
				temppos[j][0] = n_xyz[j][0];
				temppos[j][1] = n_xyz[j][2];
			}
		}
		// b <= c <= a
		if (legb <= lega && legc <= lega && legb <= legc)
		{
			tempvect[2][2] = legb;
			tempvect[0][2] = legc * cos_alpha;
			tempvect[0][0] = sqrt(legc * legc - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = lega * cos_gamma;
			tempvect[1][1] = legb * cos_alpha;
			tempvect[1][0] = sqrt(lega * lega - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legc;
			if (fabs(cos_angle - cos_bata >= ZERO_TOLERANCE1))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < pos.nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][1];
				temppos[j][0] = n_xyz[j][2];
				temppos[j][1] = n_xyz[j][0];
			}
		}
		// c <= a <= b
		if (legc <= legb && lega <= legb && legc <= lega)
		{
			tempvect[2][2] = legc;
			tempvect[0][2] = lega * cos_bata;
			tempvect[0][0] = sqrt(lega * lega - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = legb * cos_alpha;
			tempvect[1][1] = legb * cos_batar;
			tempvect[1][0] = sqrt(legb * legb - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legb;
			if (fabs(cos_angle - cos_gamma >= ZERO_TOLERANCE1))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < pos.nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][2];
				temppos[j][0] = n_xyz[j][0];
				temppos[j][1] = n_xyz[j][1];
			}
		}
		// c <= b <= a
		if (legc <= lega && legb <= lega && legc <= legb)
		{
			tempvect[2][2] = legc;
			tempvect[0][2] = legb * cos_alpha;
			tempvect[0][0] = sqrt(legb * legb - tempvect[0][2] * tempvect[0][2]);
			tempvect[1][2] = lega * cos_bata;
			tempvect[1][1] = lega * cos_alphar;
			tempvect[1][0] = sqrt(lega * lega - tempvect[1][2] * tempvect[1][2] - tempvect[1][1] * tempvect[1][1]);
			// To make the Angle unchanged
			cos_angle = (tempvect[0][0] * tempvect[1][0] + tempvect[0][1] * tempvect[1][1] + tempvect[0][2] * tempvect[1][2]) / lega / legb;
			if (fabs(cos_angle - cos_gamma >= ZERO_TOLERANCE1))
				tempvect[1][0] = -tempvect[1][0];
			for (int j = 0; j < pos.nant[0]; j++)
			{
				temppos[j][2] = n_xyz[j][2];
				temppos[j][0] = n_xyz[j][1];
				temppos[j][1] = n_xyz[j][0];
			}
		}
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				pos.vec[i][j] = tempvect[i][j];
		for (int i = 0; i < pos.nant[0]; i++)
			for (int j = 0; j < 3; j++)
				n_xyz[i][j] = temppos[i][j];
		// alpha > 90 degree, bata > 90 degree
		if (pos.vec[0][2] > 0)
		{
			pos.vec[0][0] = -pos.vec[0][0];
			pos.vec[0][2] = -pos.vec[0][2];
			for (int i = 0; i < pos.nant[0]; i++)
				n_xyz[i][0] = 1 - n_xyz[i][0];
		}
		if (pos.vec[1][2] > 0)
		{
			for (int i = 0; i < 3; i++)
				pos.vec[1][i] = -pos.vec[1][i];
			for (int i = 0; i < pos.nant[0]; i++)
				n_xyz[i][1] = 1 - n_xyz[i][1];
		}
	}

	Show_cell(pos.vec, n_xyz, n_types, pos.nant[0]);
	FILE* fp2 = fopen(file2, "w");
	savposcar(fp2, pos.title, pos.latt, pos.ifix, pos.iflg, pos.nant, unit_typenum, pos.elemnum, pos.elemsym, pos.vec, n_xyz, pos.fix);
	fclose(fp2);
	FILE* _fp_ = fopen("structure_operator", "at+");
	fprintf(_fp_, "****************************************************\n");
	fclose(_fp_);
}
