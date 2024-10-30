#pragma once
#define _POSIX_C_SOURCE 2
#ifndef _INCAR_MAGN_
#define _INCAR_MAGN_

#define MAX_ATOM_MAG_    100   // max number of atoms in amag
#define MIN_EPSILON_    1e-6

#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<cmath>
#include<regex>
#include<unordered_map>
#include"../spglib/spglib.h"
#include"../include/VASPMATE_include/read_write.h"
#include"../include/VASPMATE_include/structure_operator.h"
#include"../include/VASPMATE_include/incar.h"
#include"../include/VASPMATE_include/magn.h"

namespace incar_mag
{
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

class incar_mag
{
    public:
        incar_mag(const char file[] = "INPOS",const char file2[] = "INCAR", const char table[] = nullptr);
        void generate(char* mode[] = nullptr);
        void read_pos_mag(const char file[] = "INPOS", char* mode[] = nullptr);
        //void write_file_mag();
        void generatefile_MAGMON_list(FILE* fp_wri, string j_number, string mode, vector<double> magmon, int nant_0);
        void generatefile_INCAR(string mode , vector<double> magmon, int nant_0);
        void cal_order_mag(POSCAR pos);
        void cal_order_mode(POSCAR pos);
        vector<int> Similarity_check(vector<vector<double>> check_vector);

    private:
        MAGTYPE magtype;
        int diffposition = 0; //Whenever the Wyckoff position and the element are not the same, one will be added.
        static const int MY_num_orderings = 64;
        static const int MY_max_unique_sites = 8;
        vector<string> amag_wyckoffs;
        vector<string> amag_ele;
        vector<int> amag_equivalent_atoms;
        vector<vector<int>> combination;
    };
};

#endif