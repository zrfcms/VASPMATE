#pragma once
#ifndef _RaW4db_
#define _RaW4db_

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<stdexcept>
#include<cmath>
#include"tools.h"
#include"read_write.h"
#include"sqlitedb.h"
#include"structure_operator.h"
using namespace std;

#define CONTCAR_MEMBERS \
X(NAME, "NAME", "string", "") \
X(PATH, "PATH", "string", "") \
X(STRUCT, "STRUCT", "string", "") \
X(Classification, "Classification", "int", "") \
X(Crystal_system, "Crystal_system", "int", "") \
X(Pearson_symbol, "Pearson_symbol", "string", "") \
X(International_symbol, "International_symbol", "string", "") \
X(Lattice_a, "Lattice_a", "double", "") \
X(Lattice_b, "Lattice_b", "double", "") \
X(Lattice_c, "Lattice_c", "double", "") \
X(Lattice_alpha, "Lattice_alpha", "double", "") \
X(Lattice_beta, "Lattice_beta", "double", "") \
X(Lattice_gamma, "Lattice_gamma", "double", "") \
X(Lattice_volume, "Lattice_volume", "double", "") \
X(System_size, "System_size", "int", "") \
X(Hall_symbol, "Hall_symbol", "string", "") \
X(International_number, "International_number", "int", "") \
X(Pointgroup_symbol, "Pointgroup_symbol", "string", "") \

#define OUTCAR_MEMBERS \
X(Formation_energy, "Formation_energy", "double", "") \
X(Total_magnetization, "Total_magnetization", "double", "") \
X(Total_energy, "Total_energy", "double", "") \
X(Cohesive_energy, "Cohesive_energy", "double", "") \
X(Magnetism_ordering, "Magnetism_ordering", "int", "") \
X(Atomic_magnetic_moment, "Atomic_magnetic_moment", "string", "") \
X(POTCAR, "Potential_type", "string", "") \
X(Functional_type, "Functional_type", "string", "") \
X(PREC, "Precise", "string", "") \
X(ENCUT, "Energy_cutoff", "double", "") \
X(IALGO, "Minimization_algorithm", "string", "") \
X(EDIFF, "Electronic_convergence", "string", "") \
X(Self_consistent, "Self_consistent", "string", "") \
X(IBRION, "Intergration_scheme", "string", "") \
X(ISPIN, "Spin_polarization", "string", "") \
X(LSORBIT, "Spin_orbit_coupling", "string", "") \
X(Relaxation, "Relaxation", "string", "") \
X(Pullay_stress, "Pullay_stress", "string", "") \
X(POTIM, "Ionic_update", "string", "") \
X(EDIFFG, "Ionic_convergence", "string", "") \
X(METAGGA, "MetaGGA_type", "string", "") \
X(LHFCALC, "Hybrid_type", "string", "") \
X(LDAU, "DFT_U_type", "string", "") \
X(VDW, "VDW_D_type", "string", "") \
X(LSOL, "Solvation_model", "string", "") \
X(DIPOL, "Dipole_correction", "string", "") \
X(Energy_in_OUTCAR, "Energy_in_OUTCAR", "string", "") \
X(Stress_in_OUTCAR, "Stress_in_OUTCAR", "string", "") \
X(Force_in_OUTCAR, "Force_in_OUTCAR", "string", "") \
X(Charge_in_OUTCAR, "Charge_in_OUTCAR", "string", "") \
X(Magnitization_in_OUTCAR, "Magnitization_in_OUTCAR", "string", "") \
X(Elasticity_in_OUTCAR, "Elasticity_in_OUTCAR", "string", "") \

#define BAND_MEMBERS \
X(Band_gap, "Band_gap", "double", "") \
X(Bandgap_type, "Bandgap_type", "int", "") \
X(Bond_type, "Bond_type", "int", "") \
X(Fermi_energy, "Fermi_energy", "double", "") \
X(VBM_location, "VBM_location", "double", "") \
X(CBM_location, "CBM_location", "double", "") \
X(Band_structure, "Band_structure", "string", "") \

#define DOS_MEMBERS \
X(Density_of_state, "Density_of_state", "string", "") \

#define ELAS_MEMBERS \
X(Stiffness_tensor, "Stiffness_tensor", "string", "") \
X(Compliance_tensor, "Compliance_tensor", "string", "") \
X(Voigt_Youngs_modulus, "Voigt_Youngs_modulus", "double", "") \
X(Voigt_shear_modulus, "Voigt_shear_modulus", "double", "") \
X(Voigt_bulk_modulus, "Voigt_bulk_modulus", "double", "") \
X(Voigt_Poisson_ratio, "Voigt_Poisson_ratio", "double", "") \
X(Reuss_Youngs_modulus, "Reuss_Youngs_modulus", "double", "") \
X(Reuss_shear_modulus, "Reuss_shear_modulus", "double", "") \
X(Reuss_bulk_modulus, "Reuss_bulk_modulus", "double", "") \
X(Reuss_Poisson_ratio, "Reuss_Poisson_ratio", "double", "") \
X(Hill_Youngs_modulus, "Hill_Youngs_modulus", "double", "") \
X(Hill_shear_modulus, "Hill_shear_modulus", "double", "") \
X(Hill_bulk_modulus, "Hill_bulk_modulus", "double", "") \
X(Hill_Poisson_ratio, "Hill_Poisson_ratio", "double", "") \
X(Pugh_ratio, "Pugh_ratio", "double", "") \
X(Cauchy_pressure, "Cauchy_pressure", "double", "") \
X(Chung_Buessem_anisotropy_index, "Chung_Buessem_anisotropy_index", "double", "") \
X(Universal_elastic_anisotropy_index, "Universal_elastic_anisotropy_index", "double", "") \
X(Elastic_stability_conditions, "Elastic_stability_conditions", "double", "") \

namespace _RaW4db
{
    typedef struct Member {
        string name;
        string type;
        string value;
    } Member;

    typedef struct {
        #define X(ident, name, type, value) Member ident;
        CONTCAR_MEMBERS
        #undef X
    } CONTCAR_db;
    typedef struct {
        #define X(ident, name, type, value) Member ident;
        OUTCAR_MEMBERS
        #undef X
    } OUTCAR_db;
    typedef struct {
        #define X(ident, name, type, value) Member ident;
        BAND_MEMBERS
        #undef X
    } BAND_db;
    typedef struct {
        #define X(ident, name, type, value) Member ident;
        DOS_MEMBERS
        #undef X
    } DOS_db;
    typedef struct {
        #define X(ident, name, type, value) Member ident;
        ELAS_MEMBERS
        #undef X
    } ELAS_db;

    extern const char sdui_content[];
    extern const char sdui_Suppl[];
    extern const char sdui_Ref[];

    class RaW4db
    {
        public:
            RaW4db();
            int read_CONTCAR(const char *READPOS = "CONTCAR"); //use CONTCAR_db contcar
            int read_OUTCAR(const char *OUTCAR = "OUTCAR"); //use OUTCAR_db outcar
            int read_BANDGAP(const char *BAND_GAP = "BAND_GAP"); //use BAND_db band
            int read_DOS(); //use DOS_db dos
            int read_ELAS(const char *ELAS_INFO = "ELAS_INFO.dat"); //use ELAS_db elas
            int write_db(const char* dbname, const char* tablename, vector<string> file_name);
            int write_sdui(const char* sduiname, const char* tablename, vector<string> file_name);
            int operator_db_collect(int argc, char* argv[]);
        private:
            string readFile(const char* filename);
            string standard_content(const string& input);
            CONTCAR_db contcar;
            OUTCAR_db outcar;
            BAND_db band;
            DOS_db dos;
            ELAS_db elas;
            POSCAR pos;
            int insert_id;
            string Name_of_db; // /home/erye/test/db/db_DFT+U/vasp.db -> vasp.db
            string Name_of_sdui; // vasp.db -> vasp.sdui
    };
}

#endif

