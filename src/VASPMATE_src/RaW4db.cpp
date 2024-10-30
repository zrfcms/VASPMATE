#include"../../include/VASPMATE_include/RaW4db.h"
#include"../../include/VASPMATE_include/mapdata.h"
using namespace _RaW4db;

const char _RaW4db::sdui_content[] = "Auto_Symmetry_Analysis: On\n" \
"Name: Default\n" \ 
"Tab: Basic_Info\n" \ 
"Color: 225 225 225\n" \ 
"Textcolor: 0 0 25\n" \ 
"Group: NULL 10 10 385 560\n" \ 
"POSVIEW: NULL STRUCT 10 10 385 560\n" \ 
"GroupEnd\n" \ 
"Group: NULL 405 10 385 410\n" \ 
"Color: 220 220 255\n" \ 
"Output_Int: NULL ID 645 20 135 30\n" \ 
"Label: MaterialsID: 415 20 185 30\n" \ 
"Output_Char: NULL NAME 645 60 135 30\n" \ 
"Label: Formular: 415 60 185 30\n" \ 
"Choice: NULL Classification 645 100 135 30 Unary Binary Ternary Quaternary Others\n" \ 
"Label: Classification 415 100 185 30\n" \ 
"Output_Choice: NULL Crystal_system 645 140 135 30 Triclinic Monoclinic Orthorhombic Tetragonal Trigonal Hexagonal Cubic\n" \ 
"Label: Crystal_system: 415 140 185 30\n" \ 
"Output_Char: NULL Pearson_symbol 645 180 135 30\n" \ 
"Label: Pearson_symbol: 415 180 185 30\n" \ 
"Output_Char: NULL International_symbol 645 220 135 30\n" \ 
"Label: Space_Group 415 220 185 30\n" \ 
"Output_Float: NULL Band_gap 645 260 135 30\n" \ 
"Label: Band_Gap(eV) 415 260 185 30\n" \ 
"Output_Float: NULL Formation_energy 645 300 135 30\n" \ 
"Label: Formation_Energy(eV/Atom): 415 300 185 30\n" \ 
"Output_Choice: NULL Elastic_stability_conditions 645 340 135 30 Stable Unstable\n" \ 
"Label: Elastic_Stability: 415 340 185 30\n" \ 
"Output_Float: NULL Total_magnetization 645 380 135 30 \n" \ 
"Label: Magnetization(µB/f.u.): 415 380 185 30\n" \ 
"GroupEnd\n" \ 
"Message: NULL Description 405 430 385 140\n" \ 
"Textsize: 17\n" \ 
"TabEnd\n" \ 
"Tab: Crystal_Str.\n" \ 
"Group: Lattice: 10 30 385 290\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Float: NULL Lattice_a 255 40 135 30\n" \ 
"Label: a: 20 40 185 30\n" \ 
"Float: NULL Lattice_b 255 80 135 30\n" \ 
"Label: b: 20 80 185 30\n" \ 
"Float: NULL Lattice_c 255 120 135 30\n" \ 
"Label: c: 20 120 185 30\n" \ 
"Float: NULL Lattice_alpha 255 160 135 30\n" \ 
"Label: α: 20 160 185 30\n" \ 
"Float: NULL Lattice_beta 255 200 135 30\n" \ 
"Label: β: 20 200 185 30\n" \ 
"Float: NULL Lattice_gamma 255 240 135 30\n" \ 
"Label: γ: 20 240 185 30\n" \ 
"Float: NULL Lattice_volume 255 280 135 30\n" \ 
"Label: Volume: 20 280 185 30\n" \ 
"GroupEnd\n" \ 
"Group: Symmerty: 10 360 385 290\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Int: NULL System_size 255 370 135 30 \n" \ 
"Label: System_size 20 370 185 30\n" \ 
"Choice: NULL Crystal_system 255 410 135 30 Triclinic Monoclinic Orthorhombic Tetragonal Trigonal Hexagonal Cubic\n" \ 
"Label: Crystal_system: 20 410 185 30\n" \ 
"Char: NULL Pearson_symbol 255 450 135 30\n" \ 
"Label: Pearson_symbol: 20 450 185 30\n" \ 
"Char: NULL Hall_symbol 255 490 135 30\n" \ 
"Label: Hall_symbol: 20 490 185 30\n" \ 
"Char: NULL International_symbol 255 530 135 30\n" \ 
"Label: International_symbol: 20 530 185 30\n" \ 
"Int: NULL International_number 255 570 135 30\n" \ 
"Label: International_number: 20 570 185 30\n" \ 
"Char: NULL Pointgroup_symbol 255 610 135 30\n" \ 
"Label: Pointgroup_symbol: 20 610 185 30\n" \ 
"GroupEnd\n" \ 
"Text: Crystal_Structure: STRUCT 405 30 385 620\n" \ 
"Textsize: 17\n" \ 
"TabEnd\n" \ 
"Tab:  Electronic_Str.\n" \ 
"Group: Electronic_Property: 10 30 780 130\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Float: NULL Band_gap 255 40 135 30\n" \ 
"Label: Band_gap: 20 40 185 30\n" \ 
"Choice: NULL Bandgap_type 645 40 135 30 Direct Indirect\n" \ 
"Label: Bandgap_type: 415 40 185 30\n" \ 
"Choice: NULL Bond_type 255 80 135 30 Metallic Covalent Ionic\n" \ 
"Label: Bond_type: 20 80 185 30\n" \ 
"Float: NULL Fermi_energy 645 80 135 30\n" \ 
"Label: Fermi_energy: 415 80 185 30\n" \ 
"Float: NULL VBM_location 255 120 135 30 \n" \ 
"Label: VBM_location: 20 120 185 30\n" \ 
"Float: NULL CBM_location 645 120 135 30\n" \ 
"Label: CBM_location: 415 120 185 30\n" \ 
"GroupEnd\n" \ 
"Text: Band_Structure: Band_structure 10 200 385 370\n" \ 
"Textsize: 17\n" \ 
"Text: Density_of_State: Density_of_state 405 200 385 370\n" \ 
"Textsize: 17\n" \ 
"TabEnd\n" \ 
"Tab:  Thermodynamic_Pro.\n" \ 
"Group: Thermodynamic_Stablilty: 10 30 780 90\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Choice: NULL Thermodynamic_Stability 255 40 135 30 Unstable Stable\n" \ 
"Label: Thermodynamic_Stability: 20 40 185 30\n" \ 
"Float: NULL Formation_energy 645 40 135 30 \n" \ 
"Label: Formation_Energy(eV/Atom): 415 40 185 30\n" \ 
"Float: NULL Total_energy 255 80 135 30\n" \ 
"Label: Total_Energy(eV/Atom): 20 80 185 30\n" \ 
"Float: NULL Cohesive_energy 645 80 135 30 \n" \ 
"Label: Cohesive_Energy(eV/Atom): 415 80 185 30\n" \ 
"GroupEnd\n" \ 
"TabEnd\n" \ 
"Tab:  Mechanical_Pro.\n" \ 
"Group: Stiffness_Tensor: 10 30 385 225\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Table: NULL Stiffness_tensor 6 6 20 40 365 205\n" \ 
"GroupEnd\n" \ 
"Group: Compliance_Tensor: 10 295 385 225\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Table: NULL Compliance_tensor 6 6 20 300 365 205\n" \ 
"GroupEnd\n" \ 
"Group: Elastic_Properties: 405 30 385 490\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Float: NULL Voigt_Youngs_modulus 645 40 135 30\n" \ 
"Label: Voigt_Youngs_modulus: 415 40 185 30\n" \ 
"Float: NULL Voigt_shear_modulus 645 80 135 30\n" \ 
"Label: Voigt_shear_modulus: 415 80 185 30\n" \ 
"Float: NULL Voigt_bulk_modulus 645 120 135 30\n" \ 
"Label: Voigt_bulk_modulus: 415 120 185 30\n" \ 
"Float: NULL Voigt_Poisson_ratio 645 160 135 30\n" \ 
"Label: Voigt_Poisson_ratio: 415 160 185 30\n" \ 
"Float: NULL Reuss_Youngs_modulus 645 200 135 30\n" \ 
"Label: Reuss_Youngs_modulus: 415 200 185 30\n" \ 
"Float: NULL Reuss_shear_modulus 645 240 135 30\n" \ 
"Label: Reuss_shear_modulus: 415 240 185 30\n" \ 
"Float: NULL Reuss_bulk_modulus 645 280 135 30\n" \ 
"Label: Reuss_bulk_modulus: 415 280 185 30\n" \ 
"Float: NULL Reuss_Poisson_ratio 645 320 135 30\n" \ 
"Label: Reuss_Poisson_ratio: 415 320 185 30\n" \ 
"Float: NULL Hill_Youngs_modulus 645 360 135 30\n" \ 
"Label: Hill_Youngs_modulus: 415 360 185 30\n" \ 
"Float: NULL Hill_shear_modulus 645 400 135 30\n" \ 
"Label: Hill_shear_modulus: 415 400 185 30\n" \ 
"Float: NULL Hill_bulk_modulus 645 440 135 30\n" \ 
"Label: Hill_bulk_modulus: 415 440 185 30\n" \ 
"Float: NULL Hill_Poisson_ratio 645 480 135 30\n" \ 
"Label: Hill_Poisson_ratio: 415 480 185 30\n" \ 
"GroupEnd\n" \ 
"Group: Other_Properties: 10 560 780 130\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Float: NULL Pugh_ratio 645 570 135 30\n" \ 
"Label: Pugh_ratio: 415 570 185 30\n" \ 
"Float: NULL Cauchy_pressure 255 570 135 30\n" \ 
"Label: Cauchy_pressure: 20 570 185 30\n" \ 
"Float: NULL Chung_Buessem_anisotropy_index 255 610 135 30\n" \ 
"Label: Chung_Buessem_anisotropy_index: 20 610 185 30\n" \ 
"Float: NULL Universal_elastic_anisotropy_index  645 610 135 30\n" \ 
"Label: Universal_elastic_anisotropy_index: 415 610 185 30\n" \ 
"Choice: NULL Elastic_stability_conditions 255 650 135 30 Stable Unstable\n" \ 
"Label: Elastic_stability_conditions: 20 650 185 30\n" \ 
"GroupEnd\n" \ 
"TabEnd\n" \ 
"Tab:  Physchemical_Pro.\n" \ 
"Group: Magnetic_Properties: 10 30 385 90\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Choice: NULL Magnetism_ordering 255 40 135 30 Non-magnetic Ferromagnetic Antiferromagnetic Ferrimagnetic\n" \ 
"Label: Ordering: 20 40 185 30\n" \ 
"Float: NULL Total_magnetization 255 80 135 30 \n" \ 
"Label: Total_Magnetization(µB/f.u.): 20 80 185 30\n" \ 
"GroupEnd\n" \ 
"Text: Atomic_Magnetic_moment(µB): Atomic_magnetic_moment 405 30 385 90\n" \ 
"Textsize: 17\n" \ 
"TabEnd\n";

const char _RaW4db::sdui_Suppl[] = "Tab: Suppl.\n" \ 
"Group: Calculation_detail: 10 30 780 410\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Char: NULL Potential_type 255 40 135 30\n" \ 
"Label: Potential_type: 20 40 185 30\n" \ 
"Char: NULL Functional_type 645 40 135 30 \n" \ 
"Label: Functional_type: 415 40 185 30\n" \ 
"Char: NULL Precise 255 80 135 30\n" \ 
"Label: Precise: 20 80 185 30\n" \ 
"Float: NULL Energy_cutoff 645 80 135 30 \n" \ 
"Label: Energy_cutoff(eV): 415 80 185 30\n" \ 
"Char: NULL Minimization_algorithm 255 120 135 30\n" \ 
"Label: Minimization_algorithm 20 120 185 30\n" \ 
"Char: NULL Electronic_convergence 645 120 135 30 \n" \ 
"Label: Electronic_convergence 415 120 185 30\n" \ 
"Char: NULL Self_consistent 255 160 135 30\n" \ 
"Label: Self-consistent: 20 160 185 30\n" \ 
"Char: NULL Intergration_scheme 645 160 135 30 \n" \ 
"Label: Intergration_scheme: 415 160 185 30\n" \ 
"Char: NULL Spin_polarization 255 200 135 30\n" \ 
"Label: Spin_polarization: 20 200 185 30\n" \ 
"Char: NULL Spin_orbit_coupling 645 200 135 30 \n" \ 
"Label: Spin-orbit_coupling: 415 200 185 30\n" \ 
"Char: NULL Relaxation 255 240 135 30\n" \ 
"Label: Relaxation: 20 240 185 30\n" \ 
"Float: NULL Pullay_stress 645 240 135 30 \n" \ 
"Label: Pullay_stress(kB): 415 240 185 30\n" \ 
"Char: NULL Ionic_update 255 280 135 30\n" \ 
"Label: Ionic_update: 20 280 185 30\n" \ 
"Char: NULL Ionic_convergence 645 280 135 30\n" \ 
"Label: Ionic_Convergence: 415 280 185 30\n" \ 
"Char: NULL MetaGGA_type 255 320 135 30\n" \ 
"Label: MetaGGA_type: 20 320 185 30\n" \ 
"Char: NULL Hybrid_type 645 320 135 30 \n" \ 
"Label: Hybrid_type: 415 320 185 30\n" \ 
"Char: NULL DFT_U_type 255 360 135 30 \n" \ 
"Label: DFT+U_type: 20 360 185 30\n" \ 
"Char: NULL VDW_D_type 645 360 135 30 \n" \ 
"Label: VDW+D_type: 415 360 185 30\n" \ 
"Char: NULL Solvation_model 255 400 135 30 \n" \ 
"Label: Solvation_Model: 20 400 185 30\n" \ 
"Char: NULL Dipole_correction 645 400 135 30 \n" \ 
"Label: Dipole_Correction: 415 400 185 30\n" \ 
"GroupEnd\n" \ 
"Text: INCAR INCAR 10 480 385 225\n" \ 
"Textsize: 17\n" \ 
"Text: KPOINTS KPOINTS 405 480 385 225\n" \ 
"Textsize: 17\n" \ 
"Text: Energy Energy_in_OUTCAR 10 745 385 225\n" \ 
"Textsize: 17\n" \ 
"Text: Stress Stress_in_OUTCAR 405 745 385 225\n" \ 
"Textsize: 17\n" \ 
"Text: Force Force_in_OUTCAR 10 1010 385 225\n" \ 
"Textsize: 17\n" \ 
"Text: Charge Charge_in_OUTCAR 405 1010 385 225\n" \ 
"Textsize: 17\n" \ 
"Text: Magnitization Magnitization_in_OUTCAR 10 1275 385 225\n" \ 
"Textsize: 17\n" \ 
"Text: Elasticity Elasticity_in_OUTCAR 405 1275 385 225\n" \ 
"Textsize: 17\n";

const char _RaW4db::sdui_Ref[] = "Tab: Ref.\n" \ 
"Group: Publisher: 10 30 780 130\n" \ 
"Color: 220 220 255\n" \ 
"Textsize: 17\n" \ 
"Char:  NULL Author 255 40 135 30\n" \ 
"Label: Author: 20 40 185 30\n" \ 
"Char:  NULL Affiliation 645 40 135 30 \n" \ 
"Label: Affiliation: 415 40 185 30\n" \ 
"Char:  NULL Email 255 80 135 30\n" \ 
"Label: Email: 20 80 185 30\n" \ 
"Char:  NULL Date 645 80 135 30 \n" \ 
"Label: Date: 415 80 185 30\n" \ 
"Char:  NULL Source 255 120 135 30\n" \ 
"Label: Source: 20 120 185 30\n" \ 
"Char:  NULL Source_ID 645 120 135 30 \n" \ 
"Label: Source_ID: 415 120 185 30\n" \ 
"GroupEnd\n" \ 
"Text: NULL Reference 10 170 780 400\n" \ 
"TabEnd \n";

RaW4db::RaW4db()
{
    if(VM_fileexist("CONTCAR") && VM_fileexist("POSCAR"))
        RaW4db::read_CONTCAR("CONTCAR");
    else if(VM_fileexist("CONTCAR") && !VM_fileexist("POSCAR"))
        RaW4db::read_CONTCAR("CONTCAR");
    else if(!VM_fileexist("CONTCAR") && VM_fileexist("POSCAR"))
        RaW4db::read_CONTCAR("POSCAR");
    else
        throw logic_error("No CONTAR/POSCAR is found!"); //catch (const logic_error& e); cerr << e.what() << endl;
    if(VM_fileexist("OUTCAR"))
        RaW4db::read_OUTCAR("OUTCAR");
    if(VM_fileexist("BAND_GAP"))
        RaW4db::read_BANDGAP("BAND_GAP");
    if(VM_fileexist("PDOS_SUM.csv"))
        RaW4db::read_DOS();
    if(VM_fileexist("ELAS_INFO.dat"))
        RaW4db::read_ELAS();
}

int _RaW4db::RaW4db::operator_db_collect(int argc, char* argv[])
{
	int inc; int p_vec; int _table;
    int check_table = check_para(argc, argv, "-t", "-table", 1 , &_table);
	int check_inc = check_para(argc, argv, "-inc", "-include", 1 , &inc);
	vector<const char*> par_vec = {"-a", "-b"};
	int check_vec = check_mulpara(argc, argv, par_vec, 0 , &p_vec);
    vector<string> file_name; string db_name; string sdui_name; string table_name = "vaspinfo";
    file_name.push_back("INCAR"); file_name.push_back("KPOINTS"); file_name.push_back("CONTCAR"); file_name.push_back("POSCAR");
    if (check_vec == 1 && !strcmp("-a", argv[p_vec]))
    {
        if(!strcmp("-a", argv[p_vec]))
        {
            file_name.push_back("OSZICAR");
            file_name.push_back("OUTCAR");
        }
    }
	if(check_inc == 1)
	{
		for (int i = inc + 1; i < argc; i++)
        {
			bool add_file = true;
            if (!VM_fileexist(argv[i]))
            {
                printf("The file %s does not exist or is empty!\n", argv[i]);
                continue;
            }
            for (int j = 0; j < file_name.size(); j++)
                if (string(argv[i]) == file_name[j])
                    add_file = false;
            if (add_file)
                file_name.push_back(string(argv[i]));
        }
	}
    if (argc == 2)
        db_name = "vasp.db";
    else if (argc == 3)
    {
		if (check_vec == 1) //VASPMATE --test_db -a/-b
            db_name = "vasp.db";
		else
			db_name = string(argv[2]);
    }
    else if (argc > 3)
    {
		if (p_vec == 2 || inc == 2 || _table == 2) //VASPMATE --test_db -a/-b -include file1 file2; VASPMATE --test_db -include file1 file2
            db_name = "vasp.db";
        else
            db_name = string(argv[2]); //VASPMATE --test_db other.db -a/-b -include file1 file2; VASPMATE --test_db other.db -include file1 file2
    }
    if (check_table == 1)
        table_name = string(argv[_table + 1]);
    if (write_db(db_name.c_str(), table_name.c_str(), file_name) == 0 && write_sdui(this -> Name_of_sdui.c_str(), table_name.c_str(), file_name) == 0)
        printf("Data collection successful, stored in file: %s, table: %s, ID= %d!\n", db_name.c_str(), table_name.c_str(), this->insert_id);
    return 0;
}

int _RaW4db::RaW4db::read_CONTCAR(const char *READPOS)
{
    #define X(ident, name, type, value) contcar.ident = {name, type, value};
    CONTCAR_MEMBERS
    #undef X
    FILE* fpr = fopen(READPOS, "r");
    readposcar(fpr, this -> pos);
    fclose(fpr);
	double ntemp_vec[3][3];
	transpose_matrix(pos.vec, ntemp_vec);
	int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
	translate_typenum_type(types, pos.nant[1], pos.typenum);
	SpglibDataset* dataset = spg_get_dataset(ntemp_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
    int spacegroup_number = dataset->spacegroup_number;
    string crystal_system; int crystal_system_number; string pearson;
	if (spacegroup_number == 1 || spacegroup_number == 2)
		{crystal_system = "Triclinic"; crystal_system_number = 1; pearson += "a";}
	if (spacegroup_number >= 3 && spacegroup_number < 16)
		{crystal_system = "Monoclinic"; crystal_system_number = 2; pearson += "m";}
	if (spacegroup_number >= 16 && spacegroup_number < 75)
	    {crystal_system = "Orthorhombic"; crystal_system_number = 3; pearson += "o";}
	if (spacegroup_number >= 75 && spacegroup_number < 89)
	    {crystal_system = "Tetragonal"; crystal_system_number = 4; pearson += "t";}
	if (spacegroup_number >= 89 && spacegroup_number < 143)
		{crystal_system = "Tetragonal"; crystal_system_number = 4; pearson += "t";}
	if (spacegroup_number >= 143 && spacegroup_number < 149)
		{crystal_system = "Trigonal"; crystal_system_number = 5; pearson += "h";}
	if (spacegroup_number >= 149 && spacegroup_number < 168)
		{crystal_system = "Trigonal"; crystal_system_number = 5; pearson += "h";}
	if (spacegroup_number >= 168 && spacegroup_number < 195)
		{crystal_system = "Hexagonal"; crystal_system_number = 6; pearson += "h";}
	if (spacegroup_number >= 195 && spacegroup_number < 230)
		{crystal_system = "Cubic"; crystal_system_number = 7; pearson += "c";}
    //NAME
    string _NAME;
    for(int i = 0; i < pos.nant[1]; i++)
    {
        _NAME = _NAME + string(pos.elemsym[i]) + " " + to_string(pos.typenum[i]);
        if (i != pos.nant[1] - 1)
            _NAME = _NAME + " ";
    }
    contcar.NAME.value = _NAME;
    //PATH
    string _PATH; int Classification_number;
    if (pos.nant[1] == 1)
        {_PATH += "Unary/"; Classification_number = 1;}
    else if (pos.nant[1] == 2)
        {_PATH += "Binary/"; Classification_number = 2;}
    else if (pos.nant[1] == 3)
        {_PATH += "Ternary/"; Classification_number = 3;}
    else if (pos.nant[1] == 4)
        {_PATH += "Quaternary/"; Classification_number = 4;}
    else
        {_PATH += "Others/"; Classification_number = 5;}
    _PATH += crystal_system;
    contcar.PATH.value = _PATH;
    //STRUCT
    string temp = RaW4db::readFile("CONTCAR");
    contcar.STRUCT.value = temp;
    //Classification
    contcar.Classification.value = to_string(Classification_number);
    //Crystal_system
    contcar.Crystal_system.value = to_string(crystal_system_number);
    //Pearson_symbol
	SpglibSpacegroupType spgtype;
	spgtype = spg_get_spacegroup_type(dataset->hall_number);
	SpacegroupType spacegrouptype = spgdb_get_spacegroup_type(dataset->hall_number);
	Centering centering = spacegrouptype.centering;
	if (centering == PRIMITIVE)
		pearson += "P";
	else if (centering == BODY)
		pearson += "I";
	else if (centering == FACE)
		pearson += "F";
	else if (centering == A_FACE || centering == B_FACE || centering == C_FACE)
		pearson += "C";
	else if (centering == R_CENTER)
		pearson += "R";
    pearson += to_string(pos.nant[0]);
    contcar.Pearson_symbol.value = pearson;
    //International_symbol
    contcar.International_symbol.value = string(dataset->international_symbol);
    //Lattice_a Lattice_b Lattice_c
    double* cell_vectors = getcevec(pos.vec);
    contcar.Lattice_a.value = to_string(cell_vectors[0]);
    contcar.Lattice_b.value = to_string(cell_vectors[1]);
    contcar.Lattice_c.value = to_string(cell_vectors[2]);
    //Lattice_alpha Lattice_beta Lattice_gamma
    double* cell_angle = getangle(pos.vec);
    contcar.Lattice_alpha.value = to_string(cell_angle[0]);
    contcar.Lattice_beta.value = to_string(cell_angle[1]);
    contcar.Lattice_gamma.value = to_string(cell_angle[2]);
    //Lattice_volume
    contcar.Lattice_volume.value = to_string(volume(pos.vec));
    //System_size
    contcar.System_size.value = to_string(pos.nant[0]);
    //Hall_symbol
    contcar.Hall_symbol.value = string(dataset->hall_symbol);
    //International_number
    contcar.International_number.value = to_string(spacegroup_number);
    //Pointgroup_symbol
    char ptsymbol[6];int pt_trans_mat[3][3];
    spg_get_pointgroup(ptsymbol, pt_trans_mat, dataset->rotations, dataset->n_operations);
    contcar.Pointgroup_symbol.value = string(ptsymbol);
    spg_free_dataset(dataset);
	free(types);
    return 0;
}

int _RaW4db::RaW4db::read_OUTCAR(const char *OUTCAR)
{
    #define X(ident, name, type, value) outcar.ident = {name, type, value};
    OUTCAR_MEMBERS
    #undef X
    FILE* fpr = fopen(OUTCAR, "r");
	char buf[1024]; char temp[1024];
    fseek(fpr, 0, SEEK_SET);
    double energy_total; string mag_x, mag_y, mag_z; double mag_total_x = 0; double mag_total_y = 0; double mag_total_z = 0;
    vector <string> potcar_type; bool ferrima = false;
    int max_element = 0; int cnt = 0; bool self_con = false; bool NSW_relax = false;
	for (int i = 0; i < pos.nant[1]; i++)
		max_element = max(max_element, pos.elemnum[i]);
	while (fgets(buf, 1024, fpr))
	{
        if (strstr(buf, "magnetization (x)") != NULL)
		{
            mag_x.clear();
            int line_cnt = 0;
            mag_x += string(buf);
            while (line_cnt < pos.nant[0] && fgets(buf, 1024, fpr) != NULL)
            {
                mag_x += string(buf);
                if (buf[0] != '\n' && buf[0] != '\0' && buf[0] != '-' && buf[0] != '#')
                {
                    if(max_element <= 57)
                        sscanf(buf, "%*s%*s%*s%*s%s", &temp);
                    else
                        sscanf(buf, "%*s%*s%*s%*s%*s%s", &temp);
                    if(atof(temp) < 0)
                        ferrima = true;
                    line_cnt ++;
                    continue;
                }
            }
			while (1)
			{
				mag_x += string(buf);
				if (strstr(buf, "tot") && !strstr(buf, "#"))
                {
                    sscanf(buf, "%*s%*s%*s%*s%lf", &mag_total_x);
					break;
                }
				fgets(buf, 1024, fpr);
			}
			continue;
		}
        if (strstr(buf, "magnetization (y)") != NULL)
		{
			mag_y.clear();
			while (1)
			{
				mag_y += string(buf);
				if (strstr(buf, "tot") && !strstr(buf, "#"))
                {
                    sscanf(buf, "%*s%*s%*s%*s%lf", &mag_total_y);
					break;
                }
				fgets(buf, 1024, fpr);
			}
			continue;
		}
        if (strstr(buf, "magnetization (z)") != NULL)
		{
			mag_z.clear();
			while (1)
			{
				mag_z += string(buf);
				if (strstr(buf, "tot") && !strstr(buf, "#"))
                {
                    sscanf(buf, "%*s%*s%*s%*s%lf", &mag_total_z);
					break;
                }
				fgets(buf, 1024, fpr);
			}
			continue;
		}
        //Potential_type
        if (strstr(buf, "POTCAR") != NULL)
		{
            sscanf(buf, "%*s%s", &temp);
            outcar.POTCAR.value = string(temp);
            continue;
        }
        //Functional_type
        if (strstr(buf, "GGA type") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.Functional_type.value = string(temp);
            continue;
        }
        //Precise
        if (strstr(buf, "PREC") != NULL && !strstr(buf, "rms") && !strstr(buf, "preconditioning"))
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.PREC.value = string(temp);
            continue;
        }
        //Energy_cutoff
        if (strstr(buf, "ENCUT") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.ENCUT.value = string(temp);
            continue;
        }
        //Minimization_algorithm
        if (strstr(buf, "IALGO") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.IALGO.value = string(temp);
            continue;
        }
        //Electronic_convergence
        if (strstr(buf, "EDIFF") != NULL && !strstr(buf, "EDIFFG") && !strstr(buf, "--"))
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.EDIFF.value = string(temp);
            continue;
        }
        //Self_consistent
        if (strstr(buf, "EDIFF is reached") != NULL)
        {
            self_con = true;
            continue;
        }
        //Intergration_scheme
        if (strstr(buf, "IBRION") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.IBRION.value = string(temp);
            continue;
        }
        //Spin_polarization
        if (strstr(buf, "ISPIN") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.ISPIN.value = string(temp);
            continue;
        }
        //Spin_orbit_coupling
        if (strstr(buf, "LSORBIT") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
            outcar.LSORBIT.value = string(temp);
            continue;
        }
        //Spin_orbit_coupling
        if (strstr(buf, "NSW") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
            if(atoi(temp) > 0)
                NSW_relax = true;
            continue;
        }
        //Pullay_stress
        if (strstr(buf, "Pullay stress") != NULL)
        {
            sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%*s%s", &temp);
                outcar.Pullay_stress.value = string(temp);
            continue;
        }
        //Ionic_update
        if (strstr(buf, "POTIM") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
                outcar.POTIM.value = string(temp);
            continue;
        }
        //Ionic_convergence
        if (strstr(buf, "EDIFFG") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
                outcar.EDIFFG.value = string(temp);
            continue;
        }
        //MetaGGA_type
        if (strstr(buf, "MetaGGA") != NULL)
        {
            sscanf(buf, "%*s%s", &temp);
                outcar.METAGGA.value = string(temp);
            continue;
        }
        //Hybrid_type
        if (strstr(buf, "LHFCALC") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
                outcar.LHFCALC.value = string(temp);
            continue;
        }
        //DFT_U_type
        if (strstr(buf, "LDA+U") != NULL)
        {
            sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%*s%*s%s", &temp);
                outcar.LDAU.value = string(temp);
            continue;
        }
        //VDW_D_type
        if (strstr(buf, "IVDW") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
                outcar.VDW.value = string(temp);
            continue;
        }
        //Solvation_model
        if (strstr(buf, "LSOL") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
                outcar.LSOL.value = string(temp);
            continue;
        }
        //Dipole_correction
        if (strstr(buf, "IDIPOL") != NULL)
        {
            sscanf(buf, "%*s%*s%s", &temp);
                outcar.DIPOL.value = string(temp);
            continue;
        }
        if (strstr(buf, "free energy    TOTEN") != NULL)
        {
            int line_cnt = 0;
            outcar.Energy_in_OUTCAR.value.clear();
            while (line_cnt < 2)
			{
                if (buf[0] != '\n' && buf[0] != '\0')
                {
                    outcar.Energy_in_OUTCAR.value += string(buf);
                    line_cnt++;
                }
                if (strstr(buf, "energy(sigma->0)") != NULL)
                {
                    sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%lf", &energy_total);
                }
                fgets(buf, 1024, fpr);
			}
            continue;
        }
        if (strstr(buf, "in kB") != NULL)
        {
            outcar.Stress_in_OUTCAR.value = string(buf);
            continue;
        }
        if (strstr(buf, "TOTAL-FORCE") != NULL)
        {
			outcar.Force_in_OUTCAR.value.clear();
            while (1)
			{
				outcar.Force_in_OUTCAR.value += string(buf);
				if (strstr(buf, "total drift"))
					break;
				fgets(buf, 1024, fpr);
			}
			continue;
        }
        if (strstr(buf, "TOTAL ELASTIC MODULI (kBar)") != NULL)
        {
			outcar.Elasticity_in_OUTCAR.value.clear();
            int line_cnt = 0;
            while (line_cnt < 10)
			{
                if (buf[0] != '\n' && buf[0] != '\0')
                {
                    outcar.Elasticity_in_OUTCAR.value += string(buf);
                    line_cnt++;
                }
                fgets(buf, 1024, fpr);
			}
            continue;
        }
        if (strstr(buf, "total charge") != NULL)
        {
			if (strstr(buf, "total charge-density along one line") != NULL)
                continue;
            else
            {
                outcar.Charge_in_OUTCAR.value.clear();
                int line_cnt = 0;
                while (line_cnt < pos.nant[0] + 6)
                {
                    if (buf[0] != '\n' && buf[0] != '\0')
                    {
                        outcar.Charge_in_OUTCAR.value += string(buf);
                        line_cnt++;
                    }
                    fgets(buf, 1024, fpr);
                }
                continue;
            }
        }
    }
    //Total_energy
    outcar.Total_energy.value = to_string(energy_total);
    //Atomic_magnetic_moment
    outcar.Atomic_magnetic_moment.value = mag_x + mag_y + mag_z;
    //Formation_energy
    double energy_format = energy_total;
    for (int i = 0; i < pos.nant[1]; i++)
		energy_format -= pos.typenum[i] * enthalpies[string(pos.elemsym[i])];
	energy_format /= pos.nant[0];
    outcar.Formation_energy.value = to_string(energy_format);
    //Cohesive_energy
    //no data of single atom energy
    //Magnetism_ordering
    if(mag_total_x > 0 && ferrima == false)
        outcar.Magnetism_ordering.value = "2";
    else if(mag_total_x > 0 && ferrima == true)
        outcar.Magnetism_ordering.value = "4";
    else if(abs(mag_total_x) < 0.E-3 && ferrima == true)
        outcar.Magnetism_ordering.value = "3";
    else if(abs(mag_total_x) < 0.E-3 && ferrima == false)
        outcar.Magnetism_ordering.value = "1";
    else
        outcar.Magnetism_ordering.value = "0";
    //Total_magnetization
    if(outcar.Magnetism_ordering.value != "0")
        outcar.Total_magnetization.value = to_string(sqrt(mag_total_x*mag_total_x + mag_total_y*mag_total_y + mag_total_z*mag_total_z));
    //
    if(self_con)
        outcar.Self_consistent.value = "Reached";
    if(NSW_relax)
        outcar.Relaxation.value = "Relaxation calculation";
    else
        outcar.Relaxation.value = "Self-consistent calculation";
    if(outcar.LDAU.value == "")
        outcar.LDAU.value = "F";
    if(outcar.VDW.value == "")
        outcar.VDW.value = "F";
    if(outcar.LSOL.value == "")
        outcar.LSOL.value = "F";
    outcar.Magnitization_in_OUTCAR.value = outcar.Atomic_magnetic_moment.value;
    fclose(fpr);
    return 0;
}

int _RaW4db::RaW4db::read_BANDGAP(const char *BAND_GAP)
{
    #define X(ident, name, type, value) band.ident = {name, type, value};
    BAND_MEMBERS
    #undef X
    FILE* fpr = fopen(BAND_GAP, "r");
	char buf[1024]; char temp[1024];
    fseek(fpr, 0, SEEK_SET);
    double _Band_gap; double _Fermi_energy;
    double _VBM_location[3]; double _CBM_location[3];
	while (fgets(buf, 1024, fpr))
	{
        if (strstr(buf, "Band Gap (eV)") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%lf", &_Band_gap);
            continue;
		}
        if (strstr(buf, "Band Character") != NULL)
		{
            if (strstr(buf, "Indirect") != NULL)
                band.Bandgap_type.value = "2";
            else if (strstr(buf, "Direct") != NULL)
                band.Bandgap_type.value = "1";
            else if (strstr(buf, "Metallic") != NULL)
            {
                band.Bond_type.value = "1";
                break;                
            }
            continue;
		}
        if (strstr(buf, "Fermi Energy (eV)") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%lf", &_Fermi_energy);
            continue;
		}
        if (strstr(buf, "Location of VBM") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%lf%lf%lf", &_VBM_location[0], &_VBM_location[1], &_VBM_location[2]);
            continue;
		}
        if (strstr(buf, "Location of CBM") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%lf%lf%lf", &_CBM_location[0], &_CBM_location[1], &_CBM_location[2]);
            continue;
		}
    }
    if (band.Bond_type.value != "1")
    {
        //Band_gap
        band.Band_gap.value = to_string(_Band_gap);
        //VBM_location
        band.VBM_location.value = to_string(_VBM_location[0]) + " " + to_string(_VBM_location[1]) + to_string(_VBM_location[2]);
        //CBM_location
        band.CBM_location.value = to_string(_CBM_location[0]) + " " + to_string(_CBM_location[1]) + to_string(_CBM_location[2]);
    }
    else
    {
        FILE* fpr2 = fopen("Fermi_Energy", "r");
        fseek(fpr2, 0, SEEK_SET);
        while (fgets(buf, 1024, fpr2))
            if (strstr(buf, "#") == NULL)
                sscanf(buf, "%lf", &_Fermi_energy);
        fclose(fpr2);
    }
    //Fermi_energy
    band.Fermi_energy.value = to_string(_Fermi_energy);
    //Band_structure
    band.Band_structure.value = RaW4db::readFile("BAND_REFORMATTED.csv");
    fclose(fpr);
    return 0;
}

int _RaW4db::RaW4db::read_DOS()
{
    #define X(ident, name, type, value) dos.ident = {name, type, value};
    DOS_MEMBERS
    #undef X
    dos.Density_of_state.value = RaW4db::readFile("PDOS_SUM.csv");
    return 0;
}

int _RaW4db::RaW4db::read_ELAS(const char *ELAS_INFO)
{
    #define X(ident, name, type, value) elas.ident = {name, type, value};
    ELAS_MEMBERS
    #undef X
	char buf[1024]; char temp[1024];
    FILE* fpr = fopen(ELAS_INFO, "r");
    fseek(fpr, 0, SEEK_SET);
    double _Elastic_tensor[6][6]; double _Compliance_tensor[6][6];
    double _Voigt_approximate[4]; double _Reuss_approximate[4]; double _Hill_approximate[4];
    double _Pugh_ratio; double _Cauchy_pressure; double _Anisotropy[2];
	while (fgets(buf, 1024, fpr))
	{
        if (strstr(buf, "Elastic tensor:") != NULL)
		{
			int line_cnt = 0;
            while (line_cnt < 6 && fgets(buf, 1024, fpr) != NULL)
			{
                if (buf[0] != '\n' && buf[0] != '\0')
                {
                    sscanf(buf, "%lf%lf%lf%lf%lf%lf", &_Elastic_tensor[line_cnt][0], &_Elastic_tensor[line_cnt][1], &_Elastic_tensor[line_cnt][2], 
                                                    &_Elastic_tensor[line_cnt][3], &_Elastic_tensor[line_cnt][4], &_Elastic_tensor[line_cnt][5]);
                    line_cnt++;
                }
			}
            continue;
		}
        if (strstr(buf, "Compliance tensor:") != NULL)
		{
			int line_cnt = 0;
            while (line_cnt < 6 && fgets(buf, 1024, fpr) != NULL)
			{
                if (buf[0] != '\n' && buf[0] != '\0')
                {
                    sscanf(buf, "%lf%lf%lf%lf%lf%lf", &_Compliance_tensor[line_cnt][0], &_Compliance_tensor[line_cnt][1], &_Compliance_tensor[line_cnt][2], 
                                                    &_Compliance_tensor[line_cnt][3], &_Compliance_tensor[line_cnt][4], &_Compliance_tensor[line_cnt][5]);
                    line_cnt++;
                }
			}
            continue;
		}
        if (strstr(buf, "Voigt approximate:") != NULL)
		{
            sscanf(buf, "%*s%*s%lf%lf%lf%lf", &_Voigt_approximate[0], &_Voigt_approximate[1], &_Voigt_approximate[2], &_Voigt_approximate[3]);
            continue;
		}
        if (strstr(buf, "Reuss approximate:") != NULL)
		{
            sscanf(buf, "%*s%*s%lf%lf%lf%lf", &_Reuss_approximate[0], &_Reuss_approximate[1], &_Reuss_approximate[2], &_Reuss_approximate[3]);
            continue;
		}
        if (strstr(buf, "Hill approximate :") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%lf%lf%lf%lf", &_Hill_approximate[0], &_Hill_approximate[1], &_Hill_approximate[2], &_Hill_approximate[3]);
            continue;
		}
        if (strstr(buf, "Pugh ratio") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%lf", &_Pugh_ratio);
            continue;
		}
        if (strstr(buf, "Cauchy pressure") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%lf", &_Cauchy_pressure);
            continue;
		}
        if (strstr(buf, "Chung-Buessem") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%*s%lf", &_Anisotropy[0]);
            continue;
		}
        if (strstr(buf, "Universal Elastic") != NULL)
		{
            sscanf(buf, "%*s%*s%*s%*s%*s%lf", &_Anisotropy[1]);
            continue;
		}
        if (strstr(buf, "Conditions") != NULL)
		{
            if (strstr(buf, "Unstable") != NULL)
                elas.Elastic_stability_conditions.value = "2";
            else if (strstr(buf, "Stable") != NULL)
                elas.Elastic_stability_conditions.value = "1";
            continue;
		}
    }
    //Stiffness_tensor
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            elas.Stiffness_tensor.value = elas.Stiffness_tensor.value + to_string(_Elastic_tensor[i][j]) + "\n";
    //Compliance_tensor
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            elas.Compliance_tensor.value = elas.Compliance_tensor.value + to_string(_Compliance_tensor[i][j]) + "\n";
    //Voigt
    elas.Voigt_Youngs_modulus.value = to_string(_Voigt_approximate[0]);
    elas.Voigt_shear_modulus.value = to_string(_Voigt_approximate[1]);
    elas.Voigt_bulk_modulus.value = to_string(_Voigt_approximate[2]);
    elas.Voigt_Poisson_ratio.value = to_string(_Voigt_approximate[3]);
    //Reuss
    elas.Reuss_Youngs_modulus.value = to_string(_Reuss_approximate[0]);
    elas.Reuss_shear_modulus.value = to_string(_Reuss_approximate[1]);
    elas.Reuss_bulk_modulus.value = to_string(_Reuss_approximate[2]);
    elas.Reuss_Poisson_ratio.value = to_string(_Reuss_approximate[3]);
    //Reuss
    elas.Hill_Youngs_modulus.value = to_string(_Hill_approximate[0]);
    elas.Hill_shear_modulus.value = to_string(_Hill_approximate[1]);
    elas.Hill_bulk_modulus.value = to_string(_Hill_approximate[2]);
    elas.Hill_Poisson_ratio.value = to_string(_Hill_approximate[3]);
    //Pugh_ratio
    elas.Pugh_ratio.value = to_string(_Pugh_ratio);
    //Cauchy_pressure
    elas.Cauchy_pressure.value = to_string(_Cauchy_pressure);
    //Chung_Buessem_anisotropy_index
    elas.Chung_Buessem_anisotropy_index.value = to_string(_Anisotropy[0]);
    //Universal_elastic_anisotropy_index
    elas.Universal_elastic_anisotropy_index.value = to_string(_Anisotropy[1]);
    fclose(fpr);
    return 0;
}

int _RaW4db::RaW4db::write_db(const char* dbname, const char* tablename, vector<string> file_name)
{
    sqlite3 *db; int db_return = 0;
    db_return = VMdb::open_database(&db, dbname);
    db_return = VMdb::create_table(&db, tablename, VMdb::db_table_name);
    string Path_of_db = string(sqlite3_db_filename(db, "main"));
    size_t last_Slash = Path_of_db.rfind('/');
    this -> Name_of_db = Path_of_db.substr(last_Slash + 1);
    if (Name_of_db.size() >= 3 && Name_of_db.substr(Name_of_db.size() - 3) == ".db")
        this -> Name_of_sdui = Name_of_db.substr(0, Name_of_db.size() - 3) + ".sdui";
    //CONTCAR
    Member* contcar_members[] = 
    {
        #define X(ident, name, type, value) &contcar.ident,
        CONTCAR_MEMBERS
        #undef X
    };
    string contcar_attr = ""; vector<string> contcar_msg;
    for (int i = 0; i < sizeof(contcar_members) / sizeof(contcar_members[0]); i++)
    {
        if( i == sizeof(contcar_members) / sizeof(contcar_members[0]) - 1)
            contcar_attr = contcar_attr + contcar_members[i]->name;
        else
            contcar_attr = contcar_attr + contcar_members[i]->name + ",";
        contcar_msg.push_back(contcar_members[i]->value);
    }
    db_return = VMdb::insert_data_no_print(&db, tablename, contcar_attr.c_str(), contcar_msg);
    this -> insert_id = VMdb::get_last_id(&db);
    //OUTCAR
    Member* outcar_members[] = 
    {
        #define X(ident, name, type, value) &outcar.ident,
        OUTCAR_MEMBERS
        #undef X
    };
    string outcar_attr = ""; string outcar_condition = "ID=" + to_string(insert_id);
    for (int i = 0; i < sizeof(outcar_members) / sizeof(outcar_members[0]); i++)
    {
        if( i == sizeof(outcar_members) / sizeof(outcar_members[0]) - 1)
            outcar_attr = outcar_attr + "'" + outcar_members[i]->name + "'='" + outcar_members[i]->value + "'";
        else
            outcar_attr = outcar_attr + "'" + outcar_members[i]->name + "'='" + outcar_members[i]->value + "'" + ",";
    }
    db_return = VMdb::update_table_no_print(&db, tablename, outcar_attr.c_str(), outcar_condition.c_str());
    //BAND
    Member* band_members[] = 
    {
        #define X(ident, name, type, value) &band.ident,
        BAND_MEMBERS
        #undef X
    };
    string band_attr = ""; string band_condition = "ID=" + to_string(insert_id);
    for (int i = 0; i < sizeof(band_members) / sizeof(band_members[0]); i++)
    {
        if( i == sizeof(band_members) / sizeof(band_members[0]) - 1)
            band_attr = band_attr + "'" + band_members[i]->name + "'='" + band_members[i]->value + "'";
        else
            band_attr = band_attr + "'" + band_members[i]->name + "'='" + band_members[i]->value + "'" + ",";
    }
    db_return = VMdb::update_table_no_print(&db, tablename, band_attr.c_str(), band_condition.c_str());
    //DOS
    Member* dos_members[] = 
    {
        #define X(ident, name, type, value) &dos.ident,
        DOS_MEMBERS
        #undef X
    };
    string dos_attr = ""; string dos_condition = "ID=" + to_string(insert_id);
    for (int i = 0; i < sizeof(dos_members) / sizeof(dos_members[0]); i++)
    {
        if( i == sizeof(dos_members) / sizeof(dos_members[0]) - 1)
            dos_attr = dos_attr + "'" + dos_members[i]->name + "'='" + dos_members[i]->value + "'";
        else
            dos_attr = dos_attr + "'" + dos_members[i]->name + "'='" + dos_members[i]->value + "'" + ",";
    }
    db_return = VMdb::update_table_no_print(&db, tablename, dos_attr.c_str(), dos_condition.c_str());
    //ELAS
    Member* elas_members[] = 
    {
        #define X(ident, name, type, value) &elas.ident,
        ELAS_MEMBERS
        #undef X
    };
    string elas_attr = ""; string elas_condition = "ID=" + to_string(insert_id);
    for (int i = 0; i < sizeof(elas_members) / sizeof(elas_members[0]); i++)
    {
        if( i == sizeof(elas_members) / sizeof(elas_members[0]) - 1)
            elas_attr = elas_attr + "'" + elas_members[i]->name + "'='" + elas_members[i]->value + "'";
        else
            elas_attr = elas_attr + "'" + elas_members[i]->name + "'='" + elas_members[i]->value + "'" + ",";
    }
    db_return = VMdb::update_table_no_print(&db, tablename, elas_attr.c_str(), elas_condition.c_str());
    //File
    for (int i = 0; i < file_name.size(); i++)
    {
        //I think it’s not a good idea to put the reading process in the write section,
        //but I’m worried that the memory won't be able to hold such a large file and will crash.
        string file_temp = RaW4db::readFile(file_name[i].c_str());
        string file_attr = ""; string file_condition = "ID=" + to_string(insert_id);
        db_return = VMdb::alter_column_no_print(&db, tablename, file_name[i].c_str(), "TEXT DEFAULT NULL");
        file_attr = file_attr + "'" + file_name[i] + "'='" + file_temp + "'";
        db_return = VMdb::update_table_no_print(&db, tablename, file_attr.c_str(), file_condition.c_str());
    }
    sqlite3_close(db);
    return db_return;
}

int _RaW4db::RaW4db::write_sdui(const char* sduiname, const char* tablename, vector<string> file_name)
{
    FILE* f_sdui = fopen(sduiname, "w");
    fprintf(f_sdui, "DataBase: %s\n", Name_of_db.c_str());
    fprintf(f_sdui, "MainTable: %s\n", tablename);
    fprintf(f_sdui, "%s", sdui_content);
    fprintf(f_sdui, "%s", sdui_Suppl);
    int hight = 1275;
    for (int i = 2; i < file_name.size(); i++)
    {
        if(i % 2 == 0)
        {
            hight += 265;
            fprintf(f_sdui, "Text: %s %s 10 %d 385 225\n", file_name[i].c_str(), file_name[i].c_str(), hight);
        }
        else
        {
            fprintf(f_sdui, "Text: %s %s 405 %d 385 225\n", file_name[i].c_str(), file_name[i].c_str(), hight);
        }
        fprintf(f_sdui, "Textsize: 17\n");
    }
    fprintf(f_sdui, "TabEnd\n");
    fprintf(f_sdui, "%s", sdui_Ref);
    fclose(f_sdui);
    return 0;
}

string _RaW4db::RaW4db::readFile(const char* filename) 
{
    FILE* file = fopen(filename, "r");
    if (!file)
        return nullptr;
    fseek(file, 0, SEEK_END);
    long fsize = ftell(file);
    fseek(file, 0, SEEK_SET);
    char* content = new char[fsize + 1];
    if (!content) 
    {
        fclose(file);
        return nullptr;
    }
    fread(content, 1, fsize, file);
    fclose(file);
    content[fsize] = '\0';
    string original(content);
    delete[] content;
    return standard_content(original);
}

string _RaW4db::RaW4db::standard_content(const string& input) 
{
    string output;
    for (char c : input) 
    {
        if (c == '\'') 
        {
            output += "''";
        } 
        else 
        {
            output += c;
        }
    }
    return output;
}