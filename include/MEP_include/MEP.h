#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<sys/stat.h>
#include"../Energy_Minimization/QB/QB.h"
#include"../Energy_Minimization/QSPG/QSPG.h"
#define RATE 0.3
void init_distance_matrix(int num,double***matrix);
void get_distance_matrix(SPG_tools spg,double**matrix);
void get_distance_matrix_correction(SPG_tools spg,double**matrix);
void free_distance_matrix(int num,double**matrix);
double energy_distance_matrix(int num,double**m_lin,double**m_cur);
double energydif_distance_matrix(int num,double**m_lin,double**m_cur,double**m_cur2);
extern int MEP_main(int argc,char **argv);