#include"../../include/VASPMATE_include/aimd.h"
using namespace aimd;

int aimd::operator_aimd(int argc, char* argv[])
{
    if (argc == 2)
        return -1;
    if (!strcmp(argv[2], "-nve")) //VASPMATE --mds -nve -ts [steps]
    {
        int par_ts;
        int check_ts = check_para(argc, argv, "-ts", 1, &par_ts);
        if(check_ts == 0)
            return -1;
        else if (check_ts == 2)
            write_nve(atoi(argv[3]));
        else if (check_ts == 1)
            write_nve(atoi(argv[par_ts + 1]));
        else
            return -1;
    }
    else if (!strcmp(argv[2], "-nvt")) //VASPMATE --mds -nvt -T [temperature] -ts [steps]
    {
        int par_ts; int par_T;
        int check_ts = check_para(argc, argv, "-ts", 1, &par_ts);
        int check_T = check_para(argc, argv, "-t", "-T", 1, &par_T);
        if((check_ts || check_T) == 0)
            return -1;
        else if ((check_ts && check_T) == 2)
            write_nvt(atof(argv[3]), atoi(argv[4]));
        else if ((check_ts && check_T) == 1)
            write_nvt(atof(argv[par_T + 1]), atoi(argv[par_ts + 1]));
        else
            return -1;
    }
    else if (!strcmp(argv[2], "-npt")) //VASPMATE --mds -npt -T [temperature] -p [pressure] -ts [steps]
    {
        int par_ts; int par_T; int par_pre;
        int check_ts = check_para(argc, argv, "-ts", 1, &par_ts);
        int check_T = check_para(argc, argv, "-t", "-T", 1, &par_T);
        int check_pre = check_para(argc, argv, "-p", 1, &par_pre);
        if((check_ts || check_T || check_pre) == 0)
            return -1;
        else if ((check_ts && check_T && check_pre) == 2)
            write_npt(atof(argv[3]), atof(argv[4]), atoi(argv[5]));
        else if ((check_ts && check_T && check_pre) == 1)
            write_npt(atof(argv[par_T + 1]), atof(argv[par_pre + 1]), atoi(argv[par_ts + 1]));
        else
            return -1;
    }
    else if (!strcmp(argv[2], "-et")) //VASPMATE --mds -et/mt 
    {
        FILE* fp = fopen("OSZICAR", "r");
		if (fp == NULL)
		{
			printf("No OSZICAR file is found!");
			return 0;
		}
        get_et(fp);
        fclose(fp);
    }
    else if (!strcmp(argv[2], "-mt")) //VASPMATE --mds -et/mt 
    {
        string ISPIN = "1";
        if(VM_fileexist("INCAR"))
            ISPIN = GetInfoINCAR("ISPIN");
        else if (VM_fileexist("OUTCAR"))
            ISPIN = GetInfoOUTCAR("ISPIN");
        if (ISPIN == "1")
        {
            printf("ISPIN is 1, unable to obtain magnetic information!");
            return 0;
        }
        FILE* fp = fopen("OSZICAR", "r");
		if (fp == NULL)
		{
			printf("No OSZICAR file is found!");
			return 0;
		}
        get_mt(fp);
        fclose(fp);
    }
    else if (!strcmp(argv[2], "-ns") && argc > 3) //VASPMATE --mds -ns [time step numbers] 
    {
        if(VM_fileexist("XDATCAR"))
        {
            vector<int> time_step;
            for (int i = 3; i< argc; i++)
            {
                if (strstr(argv[i], "-") != NULL) // 
                {
                    int start, end;
                    sscanf(argv[i], "%d%*c%d", &start, &end);
                    if (start > end || start<1)
                    {
                        printf("The Time Steps %d-%d is Error!\n", start, end);
                        return 0;
                    }
                    for (int i = start; i < end + 1; i++)
                        time_step.push_back(i);
                }
                else if (IsNum(argv[i]))
                {
                    if (atoi(argv[i]) <= 0)
                    {
                        printf("The Time Steps %d is Error!\n", atoi(argv[i]));
                        return 0;
                    }
                    time_step.push_back(atoi(argv[i]));
                }
                else
                    return -1;
            }
            FILE* fp = fopen("XDATCAR", "r");
            get_postep(fp, time_step);  
            fclose(fp);
        }
        else
        {
			printf("No XDATCAR file is found!");
			return 0;            
        }
    }
    else
        return -1;
    return 0;
}

void aimd::write_nve(int steps)
{
    FILE* fp = fopen("AIMDINP", "w");
	fprintf(fp, "! Global Parameters\n");
	fprintf(fp, "  ISTART    =  0           # Read existing wavefunction; if there                            \n");
	fprintf(fp, "# ISPIN     =  1           # Spin polarised DFT                                              \n");
	fprintf(fp, "# ICHARG    =  2           # Non-self-consistent: GGA/LDA band structures                    \n");
	fprintf(fp, "  LREAL     =  Auto        # Projection operators: automatic                                 \n");
	fprintf(fp, "  ENCUT     =  520         # Cut-off energy for plane wave basis set, in eV                  \n");
	fprintf(fp, "  PREC      =  Accurate    # Precision level                                                 \n");
	fprintf(fp, "  ADDGRID   = .TRUE.       # Increase grid; helps GGA convergence                            \n");
    fprintf(fp, "! Static State Calculation\n");
	fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
	fprintf(fp, "  SIGMA     =  0.1         # Please check the width of the smearing                          \n");
	fprintf(fp, "  IBRION    =  0           # Choose molecular-dynamics                                       \n");
	fprintf(fp, "# LORBIT    =  11          # PAW radii for projected DOS                                     \n");
	fprintf(fp, "# NEDOS     =  2001        # DOSCAR points                                                   \n");
	fprintf(fp, "  NELM      =  60          # Max electronic SCF steps                                        \n");
	fprintf(fp, "! SCF energy convergence; in eV                                                              \n");
	fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence, in eV                                   \n");
	fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
    fprintf(fp, "  ISIF      =  2           # Compute stress tensor but do not change box volume/shape        \n");
	fprintf(fp, "! MD                                                                                         \n");
	fprintf(fp, "  NSW       = %5d        # Number of time steps                                        \n", steps);
	fprintf(fp, "  POTIM     =  2           # Time step in femto seconds                                      \n");
    fprintf(fp, "  SMASS     = -3                                                                             \n");
    fprintf(fp, "  MDALGO    =  1           # Using Andersen thermostat                                       \n");
    fprintf(fp, "  TEBEG     =  300         # Set the starting temperature (in K) for AIMD                    \n");
    fprintf(fp, "  ANDERSEN_PROB = 0.0      # Set Andersen collision probability to zero to get NVE enseble   \n");
    fprintf(fp, "  NBLOCK    =  1                                                                             \n");
	fprintf(fp, "! Write flags                                                                                \n");
	fprintf(fp, "  LWAVE     = .FALSE.      # Donot Write WAVECAR                                             \n");
	fprintf(fp, "  LCHARG    = .FALSE.      # Donot Write CHGCAR                                              \n");
    fclose(fp);
    printf("Write AIMDINP file!\n");
}

void aimd::write_nvt(double temp, int steps)
{
    FILE* fp = fopen("AIMDINP", "w");
	fprintf(fp, "! Global Parameters\n");
	fprintf(fp, "  ISTART    =  0           # Read existing wavefunction; if there                            \n");
	fprintf(fp, "# ISPIN     =  1           # Spin polarised DFT                                              \n");
	fprintf(fp, "# ICHARG    =  2           # Non-self-consistent: GGA/LDA band structures                    \n");
	fprintf(fp, "  LREAL     =  Auto        # Projection operators: automatic                                 \n");
	fprintf(fp, "  ENCUT     =  520         # Cut-off energy for plane wave basis set, in eV                  \n");
	fprintf(fp, "  PREC      =  Accurate    # Precision level                                                 \n");
	fprintf(fp, "  ADDGRID   = .TRUE.       # Increase grid; helps GGA convergence                            \n");
    fprintf(fp, "! Static State Calculation\n");
	fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
	fprintf(fp, "  SIGMA     =  0.1         # Please check the width of the smearing                          \n");
	fprintf(fp, "  IBRION    =  0           # Choose molecular-dynamics                                       \n");
	fprintf(fp, "# LORBIT    =  11          # PAW radii for projected DOS                                     \n");
	fprintf(fp, "# NEDOS     =  2001        # DOSCAR points                                                   \n");
	fprintf(fp, "  NELM      =  60          # Max electronic SCF steps                                        \n");
	fprintf(fp, "! SCF energy convergence; in eV                                                              \n");
	fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence, in eV                                   \n");
	fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
    fprintf(fp, "  ISIF      =  2           # Compute stress tensor but do not change box volume/shape        \n");
	fprintf(fp, "! MD                                                                                         \n");
	fprintf(fp, "  NSW       = %5d        # Number of time steps                                       \n", steps);
	fprintf(fp, "  POTIM     =  2           # Time step in femto seconds                                      \n");
    fprintf(fp, "  SMASS     =  1                                                                             \n");
    fprintf(fp, "  MDALGO    =  2           # Using Langevin thermostat                                       \n");
    fprintf(fp, "  TEBEG     = %.2f       # Set the starting temperature (in K) for AIMD                \n", temp);
    fprintf(fp, "  NBLOCK    =  1                                                                             \n");
	fprintf(fp, "! Write flags                                                                                \n");
	fprintf(fp, "  LWAVE     = .FALSE.      # Donot Write WAVECAR                                             \n");
	fprintf(fp, "  LCHARG    = .FALSE.      # Donot Write CHGCAR                                              \n");
    fclose(fp);
    printf("Write AIMDINP file!\n");
}

void aimd::write_npt(double temp, double press, int steps)
{
    FILE* fp = fopen("AIMDINP", "w");
	fprintf(fp, "! Global Parameters\n");
	fprintf(fp, "  ISTART    =  0           # Read existing wavefunction; if there                            \n");
	fprintf(fp, "# ISPIN     =  1           # Spin polarised DFT                                              \n");
	fprintf(fp, "# ICHARG    =  2           # Non-self-consistent: GGA/LDA band structures                    \n");
	fprintf(fp, "  LREAL     =  Auto        # Projection operators: automatic                                 \n");
	fprintf(fp, "  ENCUT     =  520         # Cut-off energy for plane wave basis set, in eV                  \n");
	fprintf(fp, "  PREC      =  Accurate    # Precision level                                                 \n");
	fprintf(fp, "  ADDGRID   = .TRUE.       # Increase grid; helps GGA convergence                            \n");
    fprintf(fp, "! Static State Calculation\n");
	fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
	fprintf(fp, "  SIGMA     =  0.1         # Please check the width of the smearing                          \n");
	fprintf(fp, "  IBRION    =  0           # Choose molecular-dynamics                                       \n");
	fprintf(fp, "# LORBIT    =  11          # PAW radii for projected DOS                                     \n");
	fprintf(fp, "# NEDOS     =  2001        # DOSCAR points                                                   \n");
	fprintf(fp, "  NELM      =  60          # Max electronic SCF steps                                        \n");
	fprintf(fp, "! SCF energy convergence; in eV                                                              \n");
	fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence, in eV                                   \n");
	fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
    fprintf(fp, "  ISIF      =  3           # Compute stress tensor and change box volume/shape               \n");
	fprintf(fp, "! MD                                                                                         \n");
	fprintf(fp, "  NSW       = %5d        # Number of time steps                                       \n", steps);
	fprintf(fp, "  POTIM     =  2           # Time step in femto seconds                                      \n");
    fprintf(fp, "  SMASS     =  1                                                                             \n");
    fprintf(fp, "  MDALGO    =  3           # Using Langevin thermostat                                       \n");
    fprintf(fp, "  TEBEG     =  %.2f       # Set the starting temperature (in K) for AIMD               \n", temp);
    fprintf(fp, "  TEEND     =  %.2f       # Set the end temperature (in K) for AIMD                    \n", temp);
    fprintf(fp, "  PSTRESS   =  %.2f       # Sets the external pressure in kB                          \n", press);
    fprintf(fp, "# LANGEVIN_GAMMA =10 10 10 # Langevin friction coefficient for three atomic species         \n");
    fprintf(fp, "# LANGEVIN_GAMMA_L =10.0   # Langevin friction coefficient for lattice degrees of freedom    \n");
    fprintf(fp, "  NBLOCK    =  1                                                                             \n");
	fprintf(fp, "! Write flags                                                                                \n");
	fprintf(fp, "  LWAVE     = .FALSE.      # Donot Write WAVECAR                                             \n");
	fprintf(fp, "  LCHARG    = .FALSE.      # Donot Write CHGCAR                                              \n");
    fclose(fp);
    printf("Write AIMDINP file!\n");
}

void aimd::get_et(FILE* fp)
{
    vector<double> eng_step;
	char buf[1024];
    bool contain_E0 = false;
    fseek(fp, 0, SEEK_SET);
	while (fgets(buf, 1024, fp) != NULL)
	{
		double energy_oszi = 0;
        if (strstr(buf, "E0") != NULL)
		{
			contain_E0 = true;
            sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%*s%lf", &energy_oszi);
            eng_step.push_back(energy_oszi);
		}
	}
    if (contain_E0)
    {
        FILE* fp_write = fopen("Energy_steps.dat", "w");
        fprintf(fp_write, "steps energy\n");
        for(int i = 0; i < eng_step.size(); i++)
        {
            fprintf(fp_write, "%5d %.3lf\n", i+1, eng_step[i]);
        }
        fclose(fp_write);
        printf("Write Energy_steps.dat file!\n");
    }
    else
       printf("There is no energy information included in this OSIZICAR, please check it!\n");
}

void aimd::get_mt(FILE* fp)
{
    vector<double> mag_step;
	char buf[1024];
    bool contain_mag = false;
    fseek(fp, 0, SEEK_SET);
	while (fgets(buf, 1024, fp) != NULL)
	{
		double mag_oszi = 0;
        if (strstr(buf, "mag") != NULL)
		{
			contain_mag = true;
            sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lf", &mag_oszi);
            mag_step.push_back(mag_oszi);
		}
	}
    if (contain_mag)
    {
        FILE* fp_write = fopen("Magmom.dat", "w");
        fprintf(fp_write, "steps magmom\n");
        for(int i = 0; i < mag_step.size(); i++)
        {
            fprintf(fp_write, "%5d %7.3lf\n", i+1, mag_step[i]);
        }
        fclose(fp_write);
        printf("Write Magmom.dat file!\n");
    }
    else
       printf("There is no magnetic information included in this OSIZICAR, please check it!\n");
}

void aimd::get_postep(FILE* fp, vector<int> time_step)
{
    POSCAR pos;
    readxdatcar(fp , pos);
    char buf_max[1024];
    int MAX_step = 0;
    fseek(fp, 0, SEEK_SET);
    while (fgets(buf_max, 1024, fp) != NULL)
        if (strstr(buf_max, "configuration=") != NULL)
            sscanf(buf_max, "%*s%*s%d", &MAX_step);
    for (int i = 0; i< time_step.size(); i++)
    {
        if (time_step[i] > MAX_step)
        {
            printf("Too large number: \"%d\" ! Skip it!\n", time_step[i]);
            continue;
        }
        int line = 0;
        char buf[1024];
        int find_step = 0;
        int iat = 0;
        fseek(fp, 0, SEEK_SET);
        while (fgets(buf, 1024, fp) != NULL)
        {
            double x = 0, y = 0, z = 0;
            if (line > find_step + pos.nant[0] + 1 && find_step != 0)
                break;
            else
            {
                int line_step;
                if (strstr(buf, "configuration=") != NULL)
                {
                    sscanf(buf, "%*s%*s%d", &line_step);
                    if (line_step == time_step[i])
                        find_step = line;
                }
                else if (line > find_step && find_step != 0)
                {
                    sscanf(buf, "%lf %lf %lf ", &x, &y, &z);
                    pos.xyz[iat][0] = x; pos.xyz[iat][1] = y; pos.xyz[iat][2] = z;
                    pos.fix[iat][0] = 'T'; pos.fix[iat][1] = 'T'; pos.fix[iat][2] = 'T';
                    if (++iat >= pos.nant[0])
                    {
                        line ++;
                        continue;
                    }
                }
            }
            line ++;
        }
        string pos_title = string(pos.title);
        pos_title.erase(std::remove(pos_title.begin(), pos_title.end(), ' '), pos_title.end());
        string time_step_str = to_string(time_step[i]);
        string padded_time_step = string(4 - time_step_str.length(), '0') + time_step_str;
        string postep_index = "POSTEP_" + padded_time_step; // POSTEP_0001
        string postep_title = pos_title + "_" + padded_time_step; // Fe2O3_0001
        char tem_title[80];
        strcpy(tem_title, pos.title);
        strcpy(pos.title, postep_title.c_str());
        FILE* fp_write = fopen(postep_index.c_str(), "w");
        savposcar(fp_write, pos);
        strcpy(pos.title, tem_title);
        printf("Written %s file!\n", postep_index.c_str());
    }
}