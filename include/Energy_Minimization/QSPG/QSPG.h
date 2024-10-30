#include"./../../../spglib/spglib.h"
typedef struct
{
	double (*pos)[3];
	int *type;
	int num;
	double mat[3][3];
	int ele_n=0;
	char (*elename)[32];
}SPG_tools;

/*
typedef struct {
    int spacegroup_number;
    int hall_number;
    char international_symbol[11];
    char hall_symbol[17];
    char choice[6];
    double transformation_matrix[3][3];
    double origin_shift[3];
    int n_operations;
    int (*rotations)[3][3];
    double (*translations)[3];
    int n_atoms;
    int *wyckoffs;
    char (*site_symmetry_symbols)[7];
    int *equivalent_atoms;
    int *crystallographic_orbits;
    double primitive_lattice[3][3];
    int *mapping_to_primitive;
    int n_std_atoms;
    double std_lattice[3][3];
    int *std_types;
    double (*std_positions)[3];
    double std_rotation_matrix[3][3];
    int *std_mapping_to_primitive;
    char pointgroup_symbol[6];
} SpglibDataset;
*/
#define QSPGDataset SpglibDataset
#define QSPGDataset_free_dataset spg_free_dataset
extern void SPG_free(SPG_tools* SPG);
extern void SPG2QB(SPG_tools* SPG,QB_tools *QB);
extern void QB2SPG(QB_tools *QB,SPG_tools* SPG);

extern void QSPG_refined(QB_tools* output,QB_tools* input,double symprec);
extern void QSPG_primitive(QB_tools* output,QB_tools* input,double symprec);
extern void QSPG_get_symmetry(QB_tools*input,QSPGDataset**data,double symprec);

extern void QB_f2c(QB_tools *QB,double input[3],double output[3]);
extern void QB_c2f(QB_tools *QB,double input[3],double output[3]);