#pragma once
#ifndef _MACP_
#define _MACP_

#define SYM_PREC 1e-1

#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<cmath>
#include<regex>
#include<unordered_map>
#include"../include/VASPMATE_include/tools.h"
#include"../include/VASPMATE_include/outcar.h"
#include"../include/VASPMATE_include/read_write.h"
#include"../include/VASPMATE_include/structure_operator.h"
#include"../include/VASPMATE_include/incar.h"

class Magcouple
{
public:
    Magcouple(const char file[] = "INPOS", const char table[] = nullptr);
    void generate();
    void derive();
    void generatefile_INCAR(vector<int> mag_couple_comb_i);
    void find_operation_atom(int symmetry_num, int (*rot)[3][3], double (*tra)[3], double x, double y, double z, double symprec, int flag);
    vector<int> cal_couple(vector<int> coup_vec);
private:
    POSCAR pos; //supercell POSCAR
    POSCAR pos_u; //unitcell POSCAR
    int diff_mag_position = 0;
    vector<string> mag_wyckoffs;
	vector<string> mag_ele;
    int magsymmetry_num;
    int (*magrot)[3][3];
    double (*magtra)[3];
    vector<vector<double>> mag_position; //unitcell
    vector<int> mag_super_numb; //To search for atoms that are symmetric to the original magnetic site positions after supercell, 
                        //symmetrical ones are designated as 1, 2, 3, 4.
    vector<vector<int>> mag_couple_comb;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> couple_matrix;
    Eigen::Matrix<double, Eigen::Dynamic, 1> energy_matrix;
    Eigen::Matrix<double, Eigen::Dynamic, 1> solve_matrix;
};

#endif