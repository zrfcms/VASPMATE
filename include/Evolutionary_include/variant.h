#pragma once
#ifndef _VARIANT_
#define _VARIANT_
#include"read_write.h"
#define BUFFSIZE   200

void copy_file(const  char* target, const  char* source);

double RandDigit();                                 // creating the random number from 0.0 to 1.0);L

// subroutine "DistAtom": whether the interatomic distance of any two atoms >= 1.0000 or not
int DistAtom(   int    nant[],                    // total number of atoms and types
		double &latt_o,			  // Scaling Lattice parameter
		double vec_o[3][3],               // lattice vectors
		double xyz_o[][3]);                 // atomic coordinates
void Exchange(	 int    nant_o[],                  // total number of atoms and types
		 int    typenum_o[],               // number of each type 
		 double xyz_o[][3],                // atomic coordinates
		 double xyz_e[][3]
		 );
void Strain(   char   mode[],
			   double vec_o[3][3],               // lattice vectors
               double vec_e[3][3]);


void Ripple(	 int    nant_o[],                  // total number of atoms and types
		 double xyz_o[][3],                // atomic coordinates
		 double xyz_e[][3]);



//The random move operator moves all atoms’ positions randomly, keeping the lattice vectors unchanged
void Randmove(int    nant_o[],
	double xyz_o[][3],
	double xyz_e[][3]);

//This operator slips a group of atoms by a random distance along a random direction that is parallel to a specific plane


void Slip(	 int    nant_o[],                  // total number of atoms and types
						 double xyz_o[][3],                // atomic coordinates
						 double xyz_e[][3]);

void Twist(	 int    nant_o[],                  // total number of atoms and types
						 double xyz_o[][3],                // atomic coordinates
						 double xyz_e[][3]);



void Rotation(   double vec_o[3][3],               // lattice vectors
                 double vec_e[3][3]);


void Crossover(  double &latt1_o,
								 double &latt2_o,                  // Scaling Lattice parameter 
							 	 double &latt_e,
                 int    nant_o[],                    // total number of atoms and types
                 int    typenum_o[],               // number of each type
        				 double vec1_o[3][3],
								 double vec2_o[3][3],              // lattice vectors
                 double vec_e[3][3],
					 	     double xyz1_o[][3],
		 	    			 double xyz2_o[][3],               // atomic coordinates
						 	 	 double xyz_e[][3]); 

void Stripple(   int    nant_o[],                    // total number of atoms and types
                 double vec_o[3][3],               // lattice vectors
                 double vec_e[3][3],
								 double xyz_o[][3],                // atomic coordinates
								 double xyz_e[][3]);

void Permustrain(int    nant_o[],                    // total number of atoms and types
								 int    typenum_o[],
                 double vec_o[3][3],               // lattice vectors
                 double vec_e[3][3],
								 double xyz_o[][3],                // atomic coordinates
								 double xyz_e[][3]);

void Rotstrain(  double vec_o[3][3],               // lattice vectors
                 double vec_e[3][3]);
void print_Variant_Operator();



int Variant_singer_parent(char POSCAR_Parent[], char POSCAR_Child[], int method, int Generation);



int Variant_two_parent(char POSCAR_Parent1[], char POSCAR_Parent2[],char POSCAR_Child[], int Generation );

#endif

