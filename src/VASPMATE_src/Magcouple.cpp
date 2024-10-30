#include"../../include/VASPMATE_include/Magcouple.h"
using std::cout;
using std::endl;
using std::vector;
using std::string;

void generate_mag_couple_comb_Helper(vector<int>& arr, int current, int n, vector<vector<int>>& result) 
{
    if (current == n) 
        result.push_back(arr);
    else 
    {
        for (int value : {1, -1})
        {
            arr[current] = value;
            generate_mag_couple_comb_Helper(arr, current + 1, n, result);
        }
    }
}
vector<vector<int>> generate_mag_couple_comb(int n) 
{
    vector<int> arr(n); 
    vector<vector<int>> result; 
    generate_mag_couple_comb_Helper(arr, 0, n, result);
    return result; 
}
vector<string> magcouple_vasp_file{ "INCAR","KPOINTS","POTCAR","SUPERPOS"};
map<string, double> couple_default_mag = {
    {"H",0},
    {"He",0},
    {"Li",0},
    {"Be",0},
    {"B",0},
    {"C",0},
    {"N",0},
    {"O",0},
    {"F",0},
    {"Ne",0},
    {"Na",0},
    {"Mg",0},
    {"Al",0},
    {"Si",0},
    {"P",0},
    {"S",0},
    {"Cl",0},
    {"Ar",0},
    {"K",0},
    {"Ca",0},
    {"Sc",0},
    {"Ti",0},
    {"V",5},
    {"Cr",5},
    {"Mn",5},
    {"Fe",5},
    {"Co",5},
    {"Ni",5},
    {"Cu",1.73},
    {"Zn",0},
    {"Ga",0},
    {"Ge",0},
    {"As",0},
    {"Se",0},
    {"Br",0},
    {"Kr",0},
    {"Rb",0},
    {"Sr",0},
    {"Y",0},
    {"Zr",0},
    {"Nb",0},
    {"Mo",5},
    {"Tc",0},
    {"Ru",2.2},
    {"Rh",0},
    {"Pd",0},
    {"Ag",05},
    {"Cd",0},
    {"In",0},
    {"Sn",0},
    {"Sb",0},
    {"Te",0},
    {"I",0},
    {"Xe",0},
    {"Cs",0},
    {"Ba",0},
    {"La",0},
    {"Ce",5},
    {"Pr",3.58},
    {"Nd",3.62},
    {"Pm",2.68},
    {"Sm",0.85},
    {"Eu",10},
    {"Gd",7.94},
    {"Tb",9.72},
    {"Dy",10.65},
    {"Ho",10.6},
    {"Er",9.58},
    {"Tm",7.56},
    {"Yb",4.54},
    {"Lu",0},
    {"Hf",0},
    {"Ta",0},
    {"W",5},
    {"Re",0},
    {"Os",2.2},
    {"Ir",0},
    {"Pt",0},
    {"Au",0},
    {"Hg",0},
    {"Tl",0},
    {"Pb",0},
    {"Bi",0},
    {"Po",0},
    {"At",0},
    {"Rn",0},
    {"Fr",0},
    {"Ra",0},
    {"Ac",0},
    {"Th",0},
    {"Pa",0},
    {"U",0},
    {"Np",0},
    {"Pu",0},
    {"Am",0},
    {"Cm",0},
    {"Bk",0},
    {"Cf",0},
    {"Es",0},
    {"Fm",0},
    {"Md",0},
    {"No",0},
    {"Lr",0},
    {"Rf",0},
    {"Db",0},
    {"Sg",0},
    {"Bh",0},
    {"Hs",0},
    {"Mt",0}
};

vector<int> Magcouple::cal_couple(vector<int> coup_vec)
{
    int n = coup_vec.size();
    vector<int> result(n*(n-1)/2, 1);
    for (int i = 0; i < n; i++) 
        for (int j = i + 1; j < n; j++) 
            result[i] *= coup_vec[j];
    return result;
}

void Magcouple::find_operation_atom(int symmetry_num, int (*rot)[3][3], double (*tra)[3], double x, double y, double z, double symprec, int flag)
{
    int i,j;
    double symprec2=symprec*symprec;
    //get particle id from hallnumber
    //printf("N_op=%d\n",symmetry_num);
    double vec[3];
    double dsp[3];
  
    for(i=0;i<symmetry_num;i++)
    {
        vec[0]=rot[i][0][0]*x+rot[i][0][1]*y+rot[i][0][2]*z+tra[i][0];
        vec[1]=rot[i][1][0]*x+rot[i][1][1]*y+rot[i][1][2]*z+tra[i][1];
        vec[2]=rot[i][2][0]*x+rot[i][2][1]*y+rot[i][2][2]*z+tra[i][2];
        for(j=0;j<pos.nant[0];j++)
        {
            dsp[0]=vec[0]-pos.xyz[j][0];
            dsp[1]=vec[1]-pos.xyz[j][1];
            dsp[2]=vec[2]-pos.xyz[j][2];
            while(dsp[0]<-0.5)dsp[0]+=1;
            while(dsp[1]<-0.5)dsp[1]+=1;
            while(dsp[2]<-0.5)dsp[2]+=1;
            while(dsp[0]>0.5)dsp[0]-=1;
            while(dsp[1]>0.5)dsp[1]-=1;
            while(dsp[2]>0.5)dsp[2]-=1;
            if((dsp[0]*dsp[0]+dsp[1]*dsp[1]+dsp[2]*dsp[2])<symprec2)
            {
                mag_super_numb[j] = flag;
                break;
            }
        }
    }
}

void Magcouple::generatefile_INCAR(vector<int> mag_couple_comb_i)
{
    vector<string> value_str;
    vector<double> value_dou;
    value_str.resize(pos.nant[0]);
    value_dou.resize(pos.nant[0], 1);
    for (int i = 1; i <= diff_mag_position; i++)
    {
        for (int j = 0; j < mag_super_numb.size(); j++)
        {
            if(mag_super_numb[j] == i)
            {
                value_dou[j] = mag_couple_comb_i[i - 1];
                break;
            }
            else
                continue;
        }
    }
    int numb = 0;
    for (int i = 0; i < pos.nant[1]; i++) 
    {
        for (int j = 0; j < pos.typenum[i]; j++)
        {
            value_dou[numb + j] = value_dou[numb + j] * couple_default_mag[pos.elemsym[i]];
            //cout <<  value_dou[numb + j] << endl;
        }
        numb = numb + pos.typenum[i];
    }
    INCAR_fix("ISPIN", { "2" });
    for (size_t i = 0; i < mag_super_numb.size(); i++)
        value_str[i] = to_string(value_dou[i]);
    INCAR_fix("MAGMOM", value_str);
}

Magcouple::Magcouple(const char file[], const char table[])
{
    FILE* fp = fopen(file, "r");
	if (fp == nullptr)
	{
		cout << file << " IS NOT EXIST!" << endl;
		return;
	}
    EV_unitcell("INPOS", "UNITPOS");
    FILE* fp_unit = fopen("UNITPOS", "r");
    readposcar(fp_unit, pos_u);
    fclose(fp_unit);
    int super[3] = {3, 3, 1};
    EV_supercell("UNITPOS", "SUPERPOS", super);
    FILE* fp_super = fopen("SUPERPOS", "r");
	readposcar(fp_super, pos);
	fclose(fp_super);

    //check table
	if (table != nullptr)
	{
		FILE* fpt = fopen(table, "r");
		if (fpt == nullptr)
		{
			printf("%s IS NOT EXIST!\n", table);
			printf("Use default MAG value!\n");
		}
		else
		{
			char buf[1024];
			while (fgets(buf, 1024, fpt))
			{
				char elem[10];
				double mag;
				if (sscanf(buf, "%s%lf", elem, &mag) != 2)
				{
					remove(std::begin(buf), std::end(buf), '\n');
					printf("[%s] in [%s] file may be wrong format!\n", buf, table);
					continue;
				}
				if (!couple_default_mag.count(elem))
				{
					printf("elem [%s] in [%s] file may be wrong elemnt!\n", elem, table);
					continue;
				}
				else
					couple_default_mag[elem] = mag;
			}
		}
		fclose(fpt);
	}

    //use Spglib to get the wyckoff of each mag atom
    double ntemp_vec[3][3];
    int* types = (int*)malloc(sizeof(int) * pos_u.nant[0]);
    transpose_matrix(pos_u.vec, ntemp_vec);
    SpglibDataset* dataset = spg_get_dataset(ntemp_vec, pos_u.xyz, types, pos_u.nant[0], SYMPREAC);
    mag_wyckoffs.resize(pos_u.nant[0], "n/a");
    mag_ele.resize(pos_u.nant[0], "n/a");
    int numb = 0;
	for (int i = 0; i < pos_u.nant[1]; i++)
	{
        if (couple_default_mag.count(pos_u.elemsym[i]))
        {
            if ( couple_default_mag[pos_u.elemsym[i]] >= 1e-3)
            {
                for (int j = 0; j < pos_u.typenum[i]; j++)
                {
                    mag_wyckoffs[numb + j] = to_string(dataset->wyckoffs[numb + j]);
                }
            }
        }
		else
		{
			printf("Error! %s in your INPOS files is not an element type!\n", pos_u.elemsym[i]);
			return;
		}
        numb = numb + pos_u.typenum[i];
	}
    /*for (int i = 0; i < mag_wyckoffs.size(); i++)
    {
        cout << mag_wyckoffs[i] <<endl;
    }*/
    int* types2 = (int*)malloc(sizeof(int) * pos.nant[0]);
    translate_typenum_type(types2, pos.nant[1], pos.typenum);
    transpose_matrix(pos.vec, ntemp_vec);
    SpglibDataset* dataset2 = spg_get_dataset(ntemp_vec, pos.xyz, types2, pos.nant[0], SYMPREAC);
    magsymmetry_num = dataset2->n_operations;
    magrot = dataset2->rotations;
    magtra = dataset2->translations;
    //spg_free_dataset(dataset);
    free(types);

    //Next, let’s find the position (only wyckoff) that stands out from the rest.
    mag_position.resize(0, vector<double>(3));
    int numb2 = 0;
    for (int i = 0; i < pos_u.nant[1]; i++) //Temporarily input all MAGMOM values in ferromagnetic form.
    {
        for (int j = 0; j < pos_u.typenum[i]; j++)
        {
            vector<double> posxyz = {0, 0, 0};
            mag_ele[numb2 + j] = pos_u.elemsym[i];
            if((numb2 + j == 0) && (mag_wyckoffs[numb2 + j] != "n/a"))
            {
                diff_mag_position++;
                posxyz = {pos_u.xyz[numb2 + j][0], pos_u.xyz[numb2 + j][1], pos_u.xyz[numb2 + j][2]};
                mag_position.push_back(posxyz);
            }
            else
            {
                if((mag_wyckoffs[numb2 + j] != mag_wyckoffs[numb2 + j - 1]) && (mag_wyckoffs[numb2 + j] != "n/a"))
                {
                    diff_mag_position++;
                    posxyz = {pos_u.xyz[numb2 + j][0], pos_u.xyz[numb2 + j][1], pos_u.xyz[numb2 + j][2]};
                    mag_position.push_back(posxyz);
                }
            }
            //add_order_magm[numb + j] = default_mag[pos.elemsym[i]];
        }
        numb2 = numb2 + pos_u.typenum[i];
    }

    mag_super_numb.resize(pos.nant[0]);
    for (int i = 0; i < mag_position.size(); i++)
    {
        find_operation_atom(magsymmetry_num, magrot, magtra, mag_position[i][0], mag_position[i][1], mag_position[i][2], SYM_PREC, i + 1);
        //cout << mag_position[i][0] << " " << mag_position[i][1] << " " << mag_position[i][2] << endl;
    }

    //generate_mag_couple_comb
    mag_couple_comb = generate_mag_couple_comb(diff_mag_position);
}

void Magcouple::generate()
{
    //cout << diff_mag_position << endl;
    FILE* fp_INCAR = fopen("INCAR", "r");
	if (fp_INCAR == nullptr)
	{
		cout << "INCAR IS NOT EXIST!" << endl;
		return;
	}
    for (int i = 0; i < mag_couple_comb.size(); i++)
    {
        string couple_comb = "";
        for (int j = 0; j < diff_mag_position; j++)
        {
            couple_comb = couple_comb + "_" + to_string(mag_couple_comb[i][j]);
        }
        string mkdir = "mkdir " + couple_comb;
        if (!access(couple_comb.c_str(), 0))
            cout << couple_comb << " is exist! && skip this step!" << endl;
        else
            system(mkdir.c_str());
        cout << "Created " << couple_comb << " folder!" << endl;
        chdir(couple_comb.c_str());
        for (string file : magcouple_vasp_file)
        {
            string source_file = "../" + file;
            string dest_file = "./" + file;
            if (!access(source_file.c_str(), 0))
                copy_file(dest_file.c_str(), source_file.c_str());
        }
        string source_POS = "./SUPERPOS";
		string dest_POS = "./POSCAR";
		if (!access(source_POS.c_str(), 0))
			copy_file(dest_POS.c_str(), source_POS.c_str());
        generatefile_INCAR(mag_couple_comb[i]);
        //go on
        chdir("..");
    }
}

void Magcouple::derive()
{
    int number_row = diff_mag_position*(diff_mag_position + 1)/2 + 1;
    int number_line = pow(2,diff_mag_position);
    couple_matrix.resize(number_line, number_row);
    solve_matrix.resize(number_line, 1);
    energy_matrix.resize(number_line, 1);
    for (int i = 0; i < mag_couple_comb.size(); i++)
    {
        string couple_comb = "";
        for (int j = 0; j < diff_mag_position; j++)
        {
            couple_comb = couple_comb + "_" + to_string(mag_couple_comb[i][j]);
        }
        string mkdir = "mkdir " + couple_comb;
        chdir(couple_comb.c_str());
        //go on
        energy_matrix(i, 0) = get_energy();
        chdir("..");
    }
    for (int i = 0; i < number_line; i++)
    {
        vector<int> couple_line = cal_couple(mag_couple_comb[i]);
        for (int j = 0; j < diff_mag_position*(diff_mag_position - 1)/2; j++)
        {
            couple_matrix(i, j) = couple_line[j];
        }
        for (int j = diff_mag_position*(diff_mag_position - 1)/2; j < number_row - 1; j++)
        {
            couple_matrix(i, j) = mag_couple_comb[i][j- diff_mag_position*(diff_mag_position - 1)/2];
        }
        for (int j = number_row - 1; j < number_row; j++)
        {
            couple_matrix(i, j) = 1;
        }        
    }
    solve_matrix=couple_matrix.ldlt().solve(energy_matrix);
    cout << solve_matrix <<endl;
}