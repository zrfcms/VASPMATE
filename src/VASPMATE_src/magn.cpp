#include"../../include/VASPMATE_include/magn.h"
#include"../../include/VASPMATE_include/tools.h"
#include"../../include/VASPMATE_include/outcar.h"
using std::cout;
using std::map;
using std::endl;
using std::vector;
using std::string;

static void copy_newinp(const char source_file[], const char desti_file[])
{
	FILE* fp1 = fopen(source_file, "r");
	if (fp1 == NULL)
	{
		printf("%s IS NOT EXIST!\n", source_file);
		return;
	}
	char buf[1024];
	FILE* fp2 = fopen(desti_file, "w");
	while (fgets(buf, 1024, fp1) != NULL)
		fputs(buf, fp2);
	fclose(fp1);
	fclose(fp2);
}
vector<vector<int>> findTriplets(int n) 
{
    vector<vector<int>> triplets;
    triplets.resize(0, vector<int>(3));
    for (int i = 1; i <= n ; i++)
    {
        for (int a = 1; a <= i; a++) 
        {
            if (i % a == 0) 
            {
                int i_div_a = i / a;
                for (int b = a; b <= sqrt(i_div_a); b++) 
                {
                    if (i_div_a % b == 0) 
                    {
                        int c = i_div_a / b;
                        vector<int> abc = {a, b, c};
                        if (find(triplets.begin(), triplets.end(), abc) == triplets.end())
                            triplets.push_back(abc);
                    }
                }
            }
        }
    }
    return triplets;
}
vector<string> mag_vasp_file{ "INCAR","KPOINTS","POTCAR","POSCAR"};
map<string, double> default_mag = {
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

bool areRowsSimilar(const vector<double>& row1, const vector<double>& row2) 
{
    bool all_Same = true;
    bool all_Opposite = true;
    for (size_t i = 0; i < row1.size(); i++) 
    {
        if (!(abs(row1[i] - row2[i]) < MIN_EPSILON))
            all_Same = false;
        if (!(abs(row1[i] + row2[i]) < MIN_EPSILON))
            all_Opposite = false;
        if (!all_Same && !all_Opposite) 
            return false;
    }   
    return all_Same || all_Opposite;
}
void generate_combination_Helper(vector<int>& arr, int current, int n, vector<vector<int>>& result) 
{
    if (current == n) 
        result.push_back(arr);
    else 
    {
        for (int value : {1, 0,-1}) 
        {
            arr[current] = value;
            generate_combination_Helper(arr, current + 1, n, result);
        }
    }
}
vector<vector<int>> generate_combination(int n) 
{
    vector<int> arr(n); 
    vector<vector<int>> result; 
    generate_combination_Helper(arr, 0, n, result);
    return result; 
}


magorder::magorder(const char file1[], const char file2[], const char table[])
{
    FILE* fp_pos = fopen(file1, "r");
	if (fp_pos == nullptr)
	{
		cout << file1 << " IS NOT EXIST!" << endl;
		return;
	}
    FILE* fp_incar = fopen(file2, "r");
	if (fp_incar == nullptr)
	{
		cout << file2 << " IS NOT EXIST!" << endl;
		return;
	}
    //test spg_get_dataset Wyckoff and Equivalent atom
    /*if (dataset != NULL){
    printf("Wyckoff positions:\n");
    for (int i = 0; i < dataset->n_atoms; i++)
        printf("Atom %d: Wyckoff position %d\n", i, dataset->wyckoffs[i]);
    printf("\nEquivalent atoms:\n");
    for (int i = 0; i < dataset->n_atoms; i++)
        printf("Atom %d: Equivalent atom %d\n", i, dataset->equivalent_atoms[i]);
    }
    else
        printf("Failed to obtain the dataset.\n");*/
    	//fix magmom value according to table
	if (table != nullptr)
	{
		FILE* fpt = fopen(table, "r");
		if (fpt == nullptr)
		{
			printf("%s IS NOT EXIST!\n", table);
			printf("Use default MAG value!\n");
            return;
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
				if (!default_mag.count(elem))
				{
					printf("elem [%s] in [%s] file may be wrong elemnt!\n", elem, table);
					continue;
				}
				else
					default_mag[elem] = mag;
			}
		}
		fclose(fpt);
	}
    fclose(fp_pos);
}
void magorder::read_pos_mag(const char file[], char mode[])
{
    POSCAR pos;
    FILE* fp_pos = fopen(file , "r");
    readposcar(fp_pos, pos);
    fclose(fp_pos);
    if (pos.nant[0] > MAX_ATOM_MAG)
	{
		cout << "The number of atoms in the cell is too large. If you wish to use this function, please try to reduce the size of the cell as much as possible." << endl;
		return;
	} 
    double ntemp_vec[3][3];
    transpose_matrix(pos.vec, ntemp_vec);
    int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
    translate_typenum_type(types, pos.nant[1], pos.typenum);
    SpglibDataset* dataset = spg_get_dataset(ntemp_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
    amag_wyckoffs.resize(pos.nant[0], "n/a");
    amag_ele.resize(pos.nant[0], "n/a");
    amag_equivalent_atoms.resize(pos.nant[0]);
    magtype.order_mag.resize(0, vector<double>(pos.nant[0]));
    magtype.order_mode.resize(0, vector<string>(1));
    int numb = 0;
	for (int i = 0; i < pos.nant[1]; i++)
	{
        if (default_mag.count(pos.elemsym[i]))
        {
            if ( default_mag[pos.elemsym[i]] >= 1e-3)
            {
                for (int j = 0; j < pos.typenum[i]; j++)
                {
                    amag_wyckoffs[numb + j] = to_string(dataset->wyckoffs[numb + j]);
                    //amag_equivalent_atoms[i + j] = dataset->equivalent_atoms[i + j];
                }
            }
        }
		else
		{
			printf("Error! %s in your INPOS files is not an element type!\n", pos.elemsym[i]);
			return;
		}
        numb = numb + pos.typenum[i];
	}
    /*for (int i = 0 ; i < amag_wyckoffs.size(); i++)
        cout << amag_wyckoffs[i] << endl;*/
    spg_free_dataset(dataset);
    free(types);
    //Next, let’s find the position that stands out from the rest.
    int numb2 = 0;
    for (int i = 0; i < pos.nant[1]; i++) //Temporarily input all MAGMOM values in ferromagnetic form.
    {
        for (int j = 0; j < pos.typenum[i]; j++)
        {
            amag_ele[numb2 + j] = pos.elemsym[i];
            if(numb2 + j == 0)
                diffposition++;
            else
            {
                if((amag_ele[numb2 + j] != amag_ele[numb2 + j - 1]) || (amag_wyckoffs[numb2 + j] != amag_wyckoffs[numb2 + j - 1]))
                    diffposition++;
            }
            //add_order_magm[numb + j] = default_mag[pos.elemsym[i]];
        }
        numb2 = numb2 + pos.typenum[i];
    }
    combination = generate_combination(diffposition); //{{1, 1, 1, 1, 0}, {1, 1, 1, -1, 0}, ......}
    /*for (int i = 0; i < combination.size(); i++){
        for (int j = 0; j < combination[i].size(); j++)
            cout << combination[i][j] << " ";
        cout << endl;
    }*/
    cal_order_mag(pos);
    cal_order_mode(pos);
    vector<int> erase_numb = Similarity_check(magtype.order_mag); //Remove identical rows and rows that differ by exactly -1 times.
    sort(erase_numb.begin(), erase_numb.end(), greater<int>());
    for (int i = 0; i < erase_numb.size(); i++)
    {
        magtype.order_mode.erase(magtype.order_mode.begin() + erase_numb[i]);
        magtype.order_mag.erase(magtype.order_mag.begin() + erase_numb[i]);
        //cout << " " << erase_numb[i] << endl;
    }
    /*for (int i = 0; i < magtype.order_mode.size(); i++)
    {
        cout << magtype.order_mode[i][0] << endl;
        for (int j = 0; j < pos.nant[0]; j++)
            cout << magtype.order_mag[i][j] << endl;
    }*/ 
    int j = 1;
    char file_list[1024] = "MAGMOM_list";
    FILE* fp_wri = fopen(file_list, "w");
    fprintf(fp_wri, " number:                  mode:                  MAGMOM:\n\n");
    for (int i = 0; i < magtype.order_mode.size(); i++)
	{ 
        if(strcmp("a", mode))
        {  
            bool is_mode = false;
            //for (int k = 0; mode[k] != NULL; k++) 
            if (strcmp(magtype.order_mode[i][0].c_str(), mode) == 0)
                is_mode = true;
            if(is_mode == false)
                continue;
        }
        string j_formatted = to_string(j);
        string j_number = string(3 - j_formatted.length(), '0') + j_formatted;
        j_formatted = j_number + "_" + magtype.order_mode[i][0];
        string mkdir = "mkdir " + j_formatted;
		if (!access(j_formatted.c_str(), 0))
			cout << j_formatted << " is exist! && skip this step!" << endl;
		else
			system(mkdir.c_str());
		cout << "Created " << j_formatted << " folder!" << endl;
		chdir(j_formatted.c_str());
		for (string file : mag_vasp_file)
		{
			string source_file = "../" + file;
			string dest_file = "./" + file;
			if (!access(source_file.c_str(), 0))
			    copy_file(dest_file.c_str(), source_file.c_str());
		}
        generatefile_INCAR(magtype.order_mode[i][0], magtype.order_mag[i], pos.nant[0]);
		chdir("..");
        generatefile_MAGMON_list(fp_wri, j_number, magtype.order_mode[i][0], magtype.order_mag[i], pos.nant[0]);
        j++;
	}
    fclose(fp_wri);
}
void magorder::generate(char mode[], int cell_supernumb)
{
    //cout << cell_supernumb << endl;
    if(cell_supernumb == 0)
    {
        cout << "It's detected that you need to revert to the primitive cell." << endl;
        EV_unitcell("INPOS", "UNITPOS");
        copy_file("POSCAR","UNITPOS");
        read_pos_mag("UNITPOS", mode);
    }
    else if(cell_supernumb == 1)
    {
        copy_file("POSCAR","INPOS");
        read_pos_mag("INPOS", mode);
    }
    else if(cell_supernumb < 0)
    {
        cout << "Error: The input for the cell size is a negative number, please re-enter a positive integer." << endl;
        return;
    }
    else
    {
        cout << "When it's detected that you need to expand the cell size max to "<< cell_supernumb <<", we will generate a separate folder for each super cell." << endl;
        auto triplets = findTriplets(cell_supernumb);
        for (int i = 0; i < triplets.size(); i++)
        {
            for(int j = 0; j < triplets[i].size(); j++)
            {
                cout << triplets[i][j]<< " ";
            }
            cout << endl;
        }
        for (int i = 0; i < triplets.size(); i++)
        {
            string supercell = "supercell_" + to_string(triplets[i][0]) + "_" + to_string(triplets[i][1]) + "_" + to_string(triplets[i][2]);
            string mkdir = "mkdir " + supercell;
            if (!access(supercell.c_str(), 0))
                cout << supercell << " is exist! && skip this step!" << endl;
            else
                system(mkdir.c_str());
            cout << "Created " << supercell << " folder!" << endl;
            chdir(supercell.c_str());
            for (string file : mag_vasp_file)
            {
                string source_file = "../" + file;
                string dest_file = "./" + file;
                if (!access(source_file.c_str(), 0))
                    copy_file(dest_file.c_str(), source_file.c_str());
            }
            string source_POS = "../INPOS";
			string dest_POS = "./INPOS";
			if (!access(source_POS.c_str(), 0))
			    copy_file(dest_POS.c_str(), source_POS.c_str());
            if(triplets[i][0] * triplets[i][1] * triplets[i][2] == 1)
            {
                copy_file("POSCAR","INPOS");
                read_pos_mag("INPOS", mode);
            }
            else
            {
                int super[3] = {triplets[i][0], triplets[i][1], triplets[i][2]};
                EV_supercell("INPOS", "SUPERPOS", super);
                copy_file("POSCAR","SUPERPOS");
                read_pos_mag("SUPERPOS", mode);
            }
            chdir("..");
        }
    }
}

void magorder::cal_order_mag(POSCAR pos)
{
    //Now, let's generate the magnitude of the magnetic moment for each position, which will be stored in the vector "order_mag".
    vector<double> add_order_magm;
    add_order_magm.resize(pos.nant[0]);
    for (int k = 0; k < combination.size(); k++)
    {
        int numb = 0;
        int differ = 0;
        for (int i = 0; i < pos.nant[1]; i++) 
        {
            for (int j = 0; j < pos.typenum[i]; j++)
            {
                add_order_magm[numb + j] = default_mag[pos.elemsym[i]] * combination[k][differ];
                if(numb + j == 0)
                    differ++;
                else
                {
                    if((amag_ele[numb + j] != amag_ele[numb + j - 1]) || (amag_wyckoffs[numb + j] != amag_wyckoffs[numb + j - 1]))
                        differ++;
                }
            }
            numb = numb + pos.typenum[i];
        }
        magtype.order_mag.push_back(add_order_magm);
    }
}
void magorder::cal_order_mode(POSCAR pos)
{
    vector<string> add_order_mode;
    for (int i = 0; i < magtype.order_mag.size(); i++)
    {
        bool all_zeros = true; //Check for non-magnetism ------ ISPIN 1
        bool all_positive = true; ////Check for fm-magnetism
        double plus_zeros = 0; //Check for anti/fim-magnetism
        for (double num : magtype.order_mag[i]) 
        {
            if (num != 0) 
            {
                all_zeros = false;
                break;
            }
        }
        if(all_zeros == true)
            add_order_mode = {"nfm"};
        else
        {
            for (double num : magtype.order_mag[i]) 
            {
                if (num < 0) 
                {
                    all_positive = false;
                    break;
                }
            }
            if(all_positive == true)
                add_order_mode = {"fm"};
            else
            {
                for (double num : magtype.order_mag[i]) 
                    plus_zeros += num;
                if(abs(plus_zeros) < 1e-6)
                    add_order_mode = {"afm"};
                else
                    add_order_mode = {"sfm"};   
            }
        }
        magtype.order_mode.push_back(add_order_mode);
    }
}
void magorder::generatefile_MAGMON_list(FILE* fp_wri, string j_number, string mode, vector<double> magmon, int nant_0)
{
    fprintf(fp_wri, "%7s:  %4s  " , j_number.c_str(), mode.c_str());
    if (!strcmp("clm", mode.c_str()))
        for (int j = 0; j < nant_0; j++)
            fprintf(fp_wri, " 3*%6.3f " , magmon[j]);
    else
    {
        for (int j = 0; j < nant_0; j++)
            fprintf(fp_wri, " %6.3f " , magmon[j]);
    }
    fprintf(fp_wri, " \n");
}
void magorder::generatefile_INCAR(string mode , vector<double> magmon, int nant_0)
{
    vector<string> value;
    value.resize(nant_0);
    copy_newinp("INCAR","temINP");
    if (!strcmp("nfm", mode.c_str()))
        INCAR_fix("ISPIN", { "1" });
    else   
        INCAR_fix("ISPIN", { "2" });
    for (size_t i = 0; i < magmon.size(); i++)
        value[i] = to_string(magmon[i]);
    if (strcmp("nfm", mode.c_str()))
        INCAR_fix("MAGMOM", value);
	remove("temINP");
}
void magorder::generatefile_INCAR_imag(string mode, string j_number, vector<double> magmon, int nant_0)
{
    vector<string> value;
    value.resize(nant_0);
    copy_newinp("INCAR","temINP");
    if (!strcmp("nfm", mode.c_str()))
        INCAR_fix("ISPIN", { "1" });
    else if (!strcmp("clm", mode.c_str()))
    {
        INCAR_fix("ISPIN", { "2" });
        INCAR_fix("LNONCOLLINE", { ".TRUE." });
    }
    else   
        INCAR_fix("ISPIN", { "2" });
    if (!strcmp("clm", mode.c_str()))
        for (size_t i = 0; i < magmon.size(); i++)
            value[i] = "3*" + to_string(magmon[i]);
    else
        for (size_t i = 0; i < magmon.size(); i++)
            value[i] = to_string(magmon[i]);
    if (strcmp("nfm", mode.c_str()))
        INCAR_fix("MAGMOM", value);
    copy_newinp("INCAR",j_number.c_str());
    copy_newinp("temINP","INCAR");
	remove("temINP");
}
vector<int> magorder::Similarity_check(vector<vector<double>> check_vector)
{
    for (int i = 0; i < check_vector.size(); i++)
    {
        bool shouldInvert = false;
        for (int j = 0; j < check_vector[i].size(); j++) 
        {
            if (check_vector[i][j] != 0) 
            {
                if (check_vector[i][j] < 0)
                    shouldInvert = true;
                break;
            }
        }
        if (shouldInvert) 
        {
            for (double& num : check_vector[i])
                num *= -1;
        }
    }
    vector<int> erase_numb;
    for (int i = 0; i < check_vector.size(); i++) 
        for (int j = i + 1; j <  check_vector.size(); j++) 
            if (areRowsSimilar( check_vector[i],  check_vector[j])) 
                 if (find(erase_numb.begin(), erase_numb.end(), j) == erase_numb.end())
                    erase_numb.push_back(j);
    return erase_numb;
}

void magorder::derive()
{
    bool found_super = false;
    string command = "find . -maxdepth 1 -type d -name 'supercell_*'";
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe)
    {
        found_super = true;
        char buffer_super[128];
        if (fgets(buffer_super, 128, pipe) == nullptr)
        {
            string folderName = ".";
            derive_POS(folderName);
        }
        else
        {
            while (!feof(pipe)) 
            {
                if (fgets(buffer_super, 128, pipe) != nullptr)
                {
                    string folderName(buffer_super);
                    string folderName_super;
                    regex pattern("[a-zA-Z0-9]+_[0-9]_[0-9]_[0-9]+");
                    smatch match;
                    if (regex_search(folderName, match, pattern) && !match.empty())
                        folderName_super = match[0];
                    chdir(folderName_super.c_str());
                    cout << "Enter the " << folderName_super << " folder and look for the most stable structure." << endl;
                    derive_POS(folderName_super);
                    chdir("../");
                }
            }
        }
        pclose(pipe);
    }
}

void magorder::derive_POS(string folderName)
{
    bool found_mag = false;
    string command = "find . -maxdepth 2 -type d -name '0*'";
    FILE* pipe_mag = popen(command.c_str(), "r");
    if (pipe_mag)
    {
        found_mag = true;
        char buffer_mag[128];
        string stable_dir = "none";
        double stable_eng = 1000;
        while (!feof(pipe_mag)) 
        {
            if (fgets(buffer_mag, 128, pipe_mag) != nullptr)
            {
                string folderName(buffer_mag);
                string folderName_mag;
                regex pattern("[0-9]+_[a-zA-Z0-9]+");
                smatch match;
                if (regex_search(folderName, match, pattern) && !match.empty())
                    folderName_mag = match[0];
                cout << folderName_mag <<endl;
                chdir(folderName_mag.c_str());
                double eng_each_mag = get_energy();
                if (eng_each_mag < stable_eng)
                {
                    stable_eng = eng_each_mag;
                    stable_dir = folderName_mag;
                }
                chdir("../");
            }
        }
        cout << stable_eng << "    ";
        cout << stable_dir <<endl;
        pclose(pipe_mag);
        chdir(stable_dir.c_str());
        string source_file = "./CONTCAR";
		string dest_file = "../MAGPOS";
		if (!access(source_file.c_str(), 0))
			copy_file(dest_file.c_str(), source_file.c_str());
        chdir("../");
    }
    else
    {
        cout << "Error: There is no corresponding folder here. Please double-check your command." << endl;
    }
}


//for INCAR
void magorder::incar_generate(char mode[])
{
    FILE* fp_pos_tem = fopen("INPOS", "r");
	if (fp_pos_tem == nullptr)
	{
		return;
	}
    FILE* fp_incar_tem = fopen("INCAR", "r");
	if (fp_incar_tem == nullptr)
	{
		return;
	}
    fclose(fp_pos_tem);
    fclose(fp_incar_tem);
    POSCAR pos;
    FILE* fp_pos = fopen("INPOS" , "r");
    readposcar(fp_pos, pos);
    fclose(fp_pos);
    if (pos.nant[0] > MAX_ATOM_MAG)
	{
		cout << "The number of atoms in the cell is too large. If you wish to use this function, please try to reduce the size of the cell as much as possible." << endl;
		return;
	}
    double ntemp_vec[3][3];
    transpose_matrix(pos.vec, ntemp_vec);
    int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
    translate_typenum_type(types, pos.nant[1], pos.typenum);
    SpglibDataset* dataset = spg_get_dataset(ntemp_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
    amag_wyckoffs.resize(pos.nant[0], "n/a");
    amag_ele.resize(pos.nant[0], "n/a");
    amag_equivalent_atoms.resize(pos.nant[0]);
    magtype.order_mag.resize(0, vector<double>(pos.nant[0]));
    magtype.order_mode.resize(0, vector<string>(1));
	for (int i = 0; i < pos.nant[0]; i++)
	{
        if (default_mag.count(pos.elemeach[i]))
        {
            if ( default_mag[pos.elemeach[i]] >= 1e-3)
            {
                amag_wyckoffs[i] = to_string(dataset->wyckoffs[i]);
                //amag_equivalent_atoms[i + j] = dataset->equivalent_atoms[i + j];
            }
        }
		else
		{
			printf("Error! %s in your INPOS files is not an element type!\n", pos.elemsym[i]);
			return;
		}
	}
    /*for (int i = 0 ; i < amag_wyckoffs.size(); i++)
        cout << amag_wyckoffs[i] << endl;*/
    spg_free_dataset(dataset);
    free(types);
    //Next, let’s find the position that stands out from the rest.
    for (int i = 0; i < pos.nant[0]; i++) //Temporarily input all MAGMOM values in ferromagnetic form.
    {
        amag_ele[i] = pos.elemeach[i];
        if(i == 0)
            diffposition++;
        else
        {
            if((amag_ele[i] != amag_ele[i - 1]) || (amag_wyckoffs[i] != amag_wyckoffs[i - 1]))
                diffposition++;
         }
        //add_order_magm[i] = default_mag[pos.elemsym[i]];
    }
    combination = generate_combination(diffposition); //{{1, 1, 1, 1, 0}, {1, 1, 1, -1, 0}, ......}
    /*for (int i = 0; i < combination.size(); i++){
        for (int j = 0; j < combination[i].size(); j++)
            cout << combination[i][j] << " ";
        cout << endl;
    }*/
    cal_order_mag(pos);
    cal_order_mode(pos);
    vector<int> erase_numb = Similarity_check(magtype.order_mag); //Remove identical rows and rows that differ by exactly -1 times.
    sort(erase_numb.begin(), erase_numb.end(), greater<int>());
    for (int i = 0; i < erase_numb.size(); i++)
    {
        magtype.order_mode.erase(magtype.order_mode.begin() + erase_numb[i]);
        magtype.order_mag.erase(magtype.order_mag.begin() + erase_numb[i]);
        //cout << " " << erase_numb[i] << endl;
    }
    /*for (int i = 0; i < magtype.order_mode.size(); i++)
    {
        cout << magtype.order_mode[i][0] << endl;
        for (int j = 0; j < pos.nant[0]; j++)
            cout << magtype.order_mag[i][j] << " ";
        cout << endl;
    }*/
    int j = 1;
    char file_list[1024] = "MAGMOM_list";
    FILE* fp_wri = fopen(file_list, "w");
    fprintf(fp_wri, " number:  mode:  MAGMOM:\n\n");
    for (int i = 0; i < magtype.order_mode.size(); i++)
	{ 
        if(strcmp("a", mode))
        {  
            bool is_mode = false;
            for (int k = 0; mode[k] != NULL; k++) 
                if (strcmp(magtype.order_mode[i][0].c_str(), mode) == 0)
                    is_mode = true;
            if(is_mode == false)
                continue;
        }
        string j_formatted = to_string(j);
        string j_forimag = to_string(j);
        string j_number = string(3 - j_formatted.length(), '0') + j_formatted;
        j_formatted = j_number + "_" + magtype.order_mode[i][0];
        j_forimag = "incar_" + magtype.order_mode[i][0] + "_" + j_number;
        generatefile_INCAR_imag(magtype.order_mode[i][0], j_forimag, magtype.order_mag[i], pos.nant[0]);
        generatefile_MAGMON_list(fp_wri, j_number, magtype.order_mode[i][0], magtype.order_mag[i], pos.nant[0]);
        j++;
	}
    if(!strcmp("clm", mode))
    {
        string j_forclm = "incar_clm_" + string(3 - to_string(j).length(), '0') + to_string(j);
        generatefile_INCAR_imag("clm", j_forclm, magtype.order_mag[0], pos.nant[0]);
        generatefile_MAGMON_list(fp_wri, string(3 - to_string(j).length(), '0') + to_string(j), "clm", magtype.order_mag[0], pos.nant[0]);
    }
    fclose(fp_wri);
}