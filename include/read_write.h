#pragma once
#ifndef _POSCAR_
#define _POSCAR_

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/time.h>
#include <fcntl.h>
#include <sys/types.h>

#define MAX_FNAME       128   // max length of file name
#define MAX_NCOLU       1024    // max number of columns of the input file
#define MAX_NTEXT       50    // max number of text fields on each line of the input file
#define MAX_NLINE		100	  // max number lines of file
#define MAX_NELEM       50    // max number of elements, 64
#define MAX_NATOM       5000   // max number of atoms, 1024
#define MAX_NREAD       10000 // max number of string for reading
#define MAX_NCHAR       80    // max number of char for title

const double Pi = acos(-1);

typedef struct
{
	char   st[MAX_NCHAR];//the name of title
	int    na;  // nant[0],nant[1] number of atoms and types
	int    nt;
	int    ifix; // Select 0 or normal 1
	int    iflg;       // Direct 0 or Cartes 1
	double abc;// Scaling Lattice parameter 
	double ax; double ay; double az;//lattice vectors 
	double bx; double by; double bz;
	double cx; double cy; double cz;
} VEC;

typedef struct
{
	int id;  //atomic number
	int nt;//number of atoms 
} TYP;

typedef struct
{
	double px;
	double py;
	double pz;
	char   fx;
	char   fy;
	char   fz;
} POS;

typedef struct POSCAR
{
	char   title[MAX_NCHAR];        // Title of system
	double latt;                    // Scaling factor
	int    ifix;                    // Select 0 or normal 1
	int    iflg;                    // Direct 0 or Cartes 1
	int    nant[2];                 // number of atoms and types
	int    typenum[MAX_NELEM];      // number of atoms of each type
	int    elemnum[MAX_NELEM];      // atomic number
	char   elemsym[MAX_NELEM][3];   // atomic symbol
	double vec[3][3];               // lattice vector
	double xyz[MAX_NATOM][3];		// atomic coordinate
	char fix[MAX_NATOM][3];         // Selective fix on each atom
}POSCAR;

void has_vacuum_slab_or_not(POSCAR pos, bool has_vacuum_slab[3]);

bool determine_has_vacuum_slab(POSCAR pos, int direction);

POSCAR delete_atom(POSCAR pos, int n);

void carts_to_direct(int& iflg, double vec[3][3], double xyz[][3], int nant[]);

void direct_to_carts(int& iflg, double vec[3][3], double xyz[][3], int nant[]);

void getAtomNum(char elemsym[3],    // String to store Element name
	int& elemnum);        // atomic number

void getAtomSym(const int& elemnum,     // atomic number
	char elemsym[3]);    // String to store Element name

int savposcar(FILE* fp,
	char   title[],                   // title of the system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                      // nant[0],nant[1] number of atoms and types
	int    typenum[],                 // number of atoms for each ion type
	int    elemnum[],                 // atomic number of each element
	char   elemsym[][3],              // atomic symbol of each element
	double vec[3][3],                 // lattice vectors
	double xyz[][3],                  // atomic coordinates
	char   fix[][3]);                   // Selective fix on each atom

int savposcar(FILE* fp, POSCAR poscar);

int readposcar(FILE* fp, POSCAR& poscar);

int readposcar(FILE* fp,
	char   title[],                   // title of the system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                      // nant[0],nant[1] number of atoms and types
	int    typenum[],                 // number of atoms for each ion type
	int    elemnum[],                 // atomic number of each element
	char   elemsym[][3],              // atomic symbol of each element
	double vec[3][3],                 // lattice vectors
	double xyz[][3],                  // atomic coordinates
	char   fix[][3]);                   // Selective fix on each atom

int savposbin(FILE* fp,
	char   title[],                   // title of the system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                      // nant[0],nant[1] number of atoms and types
	int    typenum[],                 // number of atoms for each ion type
	int    elemnum[],                 // atomic number of each element
	char   elemsym[][3],              // atomic symbol of each element
	double vec[3][3],                 // lattice vectors
	double xyz[][3],                  // atomic coordinates
	char   fix[][3]);                   // Selective fix on each atom

int readposbin(FILE* fp,
	char   title[],                   // title of the system
	double& latt,                     // Scaling Lattice parameter 
	int& ifix,                     // Select 0 or normal 1
	int& iflg,                     // Direct 0 or Cartes 1
	int    nant[],                      // nant[0],nant[1] number of atoms and types
	int    typenum[],                 // number of atoms for each ion type
	int    elemnum[],                 // atomic number of each element
	char   elemsym[][3],              // atomic symbol of each element
	double vec[3][3],                 // lattice vectors
	double xyz[][3],                  // atomic coordinates
	char   fix[][3]);                   // Selective fix on each atom

#endif


