#include"../../include/VASPMATE_include/tools.h"
#include"../../include/VASPMATE_include/incar.h"
#include <iostream>

using namespace std;

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

void print_global(FILE* fp)
{
	fprintf(fp, "! Global Parameters\n");
	fprintf(fp, "  ISTART    =  0           # Read existing wavefunction; if there                            \n");
	fprintf(fp, "# ISPIN     =  1           # Spin polarised DFT                                              \n");
	fprintf(fp, "# ICHARG    =  2           # Non-self-consistent: GGA/LDA band structures                    \n");
	fprintf(fp, "  LREAL     = .FALSE.      # Projection operators: automatic                                 \n");
	fprintf(fp, "  ENCUT     =  520         # Cut-off energy for plane wave basis set, in eV                  \n");
	fprintf(fp, "  PREC      =  Accurate    # Precision level                                                 \n");
	fprintf(fp, "# LWAVE     = .TRUE.       # Write WAVECAR or not                                            \n");
	fprintf(fp, "# LCHARG    = .TRUE.       # Write CHGCAR or not                                             \n");
	fprintf(fp, "  ADDGRID   = .TRUE.       # Increase grid; helps GGA convergence                            \n");
	fprintf(fp, "# LVTOT     = .TRUE.       # Write total electrostatic potential into LOCPOT or not          \n");
	fprintf(fp, "# LVHAR     = .TRUE.       # Write ionic + Hartree electrostatic potential into LOCPOT or not\n");
	fprintf(fp, "# NELECT    =              # No. of electrons: charged cells; be careful                     \n");
	fprintf(fp, "# LPLANE    = .TRUE.       # Real space distribution; supercells                             \n");
	fprintf(fp, "# NPAR      =  4           # Max is no. nodes; don't set for hybrids                         \n");
	fprintf(fp, "# NWRITE    =  2           # Medium-level output                                             \n");
	fprintf(fp, "# KPAR      =  2           # Divides k-grid into separate groups                             \n");
	fprintf(fp, "# NGX       =  500         # FFT grid mesh density for nice charge/potential plots           \n");
	fprintf(fp, "# NGY       =  500         # FFT grid mesh density for nice charge/potential plots           \n");
	fprintf(fp, "# NGZ       =  500         # FFT grid mesh density for nice charge/potential plots           \n");
	fprintf(fp, "\n");
}
string write_INCAR(int argc, char* argv[])
{
	std::string filename = "incar_";
	for (int i = 2; i < argc; i++)
	{
		filename += argv[i];
		if (i != argc - 1)
			filename += "_";
	}
	FILE* fp = fopen(filename.c_str(), "w");
	for (int i = 2; i < argc; i++)
	{
		if (!strcmp(argv[i], "stc") || !strcmp(argv[i], "ssc")) {
			print_global(fp);
			fprintf(fp, "! Static State Calculation\n");
			fprintf(fp, "  NSW       =  0           # Number of ionic steps                                           \n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
			fprintf(fp, "  SIGMA     =  0.1         # Please check the width of the smearing                          \n");
			fprintf(fp, "  IBRION    = -1           # Algorithm: 0-MD; 1-Quasi-New; 2-CG                              \n");
			fprintf(fp, "# LORBIT    =  11          # PAW radii for projected DOS                                     \n");
			fprintf(fp, "# NEDOS     =  2001        # DOSCAR points                                                   \n");
			fprintf(fp, "  NELM      =  60          # Max electronic SCF steps                                        \n");
			fprintf(fp, "\n");
			fprintf(fp, "! SCF energy convergence; in eV                                                              \n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence, in eV                                   \n");
			fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
		}
		if (!strcmp(argv[i], "nsc")) {
			print_global(fp);
			fprintf(fp, "! Non-Selfconsistent Calculation\n");
			fprintf(fp, "  NSW       =  0           # Number of ionic steps                                           \n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
			fprintf(fp, "  SIGMA     =  0.1         # Please check the width of the smearing                          \n");
			fprintf(fp, "  IBRION    = -1           # Algorithm: 0-MD; 1-Quasi-New; 2-CG                              \n");
			fprintf(fp, "# LORBIT    =  11          # PAW radii for projected DOS                                     \n");
			fprintf(fp, "# NEDOS     =  2001        # DOSCAR points                                                   \n");
			fprintf(fp, "  ICHARG    =  12          # Non-selfconsistent calculations; superposition charge densities \n");
			fprintf(fp, "  NELM      =  0           # Max electronic SCF steps                                        \n");
		}
		if (!strcmp(argv[i], "rlx") || !strcmp(argv[i], "lar")) {
			print_global(fp);
			fprintf(fp, "! Lattice Atomic Relaxtion\n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
			fprintf(fp, "  NSW       =  300         # Number of ionic steps                                           \n");
			fprintf(fp, "  SIGMA     =  0.05        # Please check the width of the smearing                          \n");
			fprintf(fp, "  IBRION    =  2           # Algorithm: 0-MD; 1-Quasi-New; 2-CG                              \n");
			fprintf(fp, "  ISIF      =  3           # Optimize atomic coordinates and lattice parameters              \n");
			fprintf(fp, "\n");
			fprintf(fp, "! SCF energy convergence; in eV\n");
			fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence, in eV                                   \n");
		}
		if (!strcmp(argv[i], "mds")) {
			print_global(fp);
			fprintf(fp, "! Electronic Relaxation\n");
			fprintf(fp, "  ISMEAR    =  0                                                                             \n");
			fprintf(fp, "  SIGMA     =  0.05                                                                          \n");
			fprintf(fp, "  EDIFF     =  1E-06                                                                         \n");
			fprintf(fp, "\n");
			fprintf(fp, "! Molecular Dynamics\n");
			fprintf(fp, "  IBRION    =  0           # Activate MD                                                     \n");
			fprintf(fp, "  NSW       =  100         # Max ionic steps                                                 \n");
			fprintf(fp, "  EDIFFG    = -1E-02       # Ionic convergence; eV/A                                         \n");
			fprintf(fp, "  POTIM     =  1           # Timestep in fs                                                  \n");
			fprintf(fp, "  SMASS     =  0           # MD Algorithm: -3-microcanonical ensemble; 0-canonical ensemble  \n");
			fprintf(fp, "# TEBEG     =  100         # Start temperature K                                             \n");
			fprintf(fp, "# TEEND     =  100         # Final temperature K                                             \n");
			fprintf(fp, "# MDALGO    =  1           # Andersen Thermostat                                             \n");
			fprintf(fp, "# ISYM      =  0           # Symmetry: 0=none; 2=GGA; 3=hybrids                              \n");
		}
		if (!strcmp(argv[i], "sar")) {
			print_global(fp);
			fprintf(fp, "! Electronic Relaxation\n");
			fprintf(fp, "  NELM      =  60          # Max electronic SCF steps                                        \n");
			fprintf(fp, "  NELMIN    =  6           # Min electronic SCF steps                                        \n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence; in eV                                   \n");
			fprintf(fp, "# GGA       =  PS          # PBEsol exchange-correlation                                     \n");
			fprintf(fp, "\n");
			fprintf(fp, "! Ionic Relaxation\n");
			fprintf(fp, "  NSW       =  100         # Max ionic steps                                                 \n");
			fprintf(fp, "  IBRION    =  2           # Algorithm: 0-MD; 1-Quasi-New; 2-CG                              \n");
			fprintf(fp, "  ISIF      =  2           # Stress/relaxation: 2-Ions, 3-Shape/Ions/V, 4-Shape/Ions         \n");
			fprintf(fp, "  EDIFFG    = -2E-02       # Ionic convergence; eV/AA                                        \n");
			fprintf(fp, "# ISYM      =  2           # Symmetry: 0=none; 2=GGA; 3=hybrids                              \n");
		}
		if (!strcmp(argv[i], "mpc") || !strcmp(argv[i], "mag")) {
			fprintf(fp, "! Magnetic Proproty Calculation\n");
			fprintf(fp, "  ISPIN     =  2           # Spin polarised DFT                                              \n");
			fprintf(fp, "# MAGMOM    =              # Set this parameters manually                                    \n");
			fprintf(fp, "  LASPH     = .TRUE.       # Non-spherical elements; d/f convergence                         \n");
			fprintf(fp, "  GGA_COMPAT= .FALSE.      # Apply spherical cutoff on gradient field                        \n");
			fprintf(fp, "  VOSKOWN   =  1           # Enhances the magnetic moments and the magnetic energies         \n");
			fprintf(fp, "  LMAXMIX   =  4           # For d elements increase LMAXMIX to 4, f: LMAXMIX = 6            \n");
			fprintf(fp, "# AMIX      =  0.2         # Mixing parameter to control SCF convergence                     \n");
			fprintf(fp, "# BMIX      =  0.0001      # Mixing parameter to control SCF convergence                     \n");
			fprintf(fp, "# AMIX_MAG  =  0.4         # Mixing parameter to control SCF convergence                     \n");
			fprintf(fp, "# BMIX_MAG  =  0.0001      # Mixing parameter to control SCF convergence                     \n");
		}
		if (!strcmp(argv[i], "opc") || !strcmp(argv[i], "opt")) {
			fprintf(fp, "! Optical Property Calculation\n");
			fprintf(fp, "  ALGO      =  Exact                                                                         \n");
			fprintf(fp, "  LOPTICS   = .TRUE.       # Switch on optic calculations                                    \n");
			fprintf(fp, "  CSHIFT    =  0.100       # Complex shift in Kramers-Kronig transformation                  \n");
			fprintf(fp, "  NEDOS     =  2000        # Number of grid points in DOS                                    \n");
			fprintf(fp, "# NBANDS    =              # Set this parameters manually, number of bands to include        \n");
		}
		if (!strcmp(argv[i], "soc")) {
			fprintf(fp, "! Spin-Orbit Coupling Calculation\n");
			fprintf(fp, "  LSORBIT   = .TRUE.       # Activate SOC                                                    \n");
			fprintf(fp, "  GGA_COMPAT= .FALSE.      # Apply spherical cutoff on gradient field                        \n");
			fprintf(fp, "  VOSKOWN   =  1           # Enhances the magnetic moments and the magnetic energies         \n");
			fprintf(fp, "  LMAXMIX   =  4           # For d elements increase LMAXMIX to 4, f: LMAXMIX = 6            \n");
			fprintf(fp, "# SAXIS     =  0 0 1       # Direction of the magnetic field                                 \n");
			fprintf(fp, "# MAGMOM    =  0 0 3       # Set this parameters manually, Local magnetic moment parallel to SAXIS, 3*NIONS*1.0 for non-collinear magnetic systems\n");
			fprintf(fp, "# NBANDS    =              # Set this parameters manually, 2 * number of bands of collinear-run\n");
		}
		if (!strcmp(argv[i], "hse")) {
			fprintf(fp, "! Hybrid function HSE06\n");
			fprintf(fp, "  LHFCALC   = .TRUE.       # Activate HF                                                     \n");
			fprintf(fp, "  AEXX      =  0.25        # 25% HF exact exchange, adjusted this value to reproduce experimental band gap\n");
			fprintf(fp, "  HFSCREEN  =  0.2         # Switch to screened exchange; e.g. HSE06                         \n");
			fprintf(fp, "  ALGO      =  ALL         # Electronic Minimisation Algorithm; ALGO=58                      \n");
			fprintf(fp, "  TIME      =  0.4         # Timestep for IALGO5X                                            \n");
			fprintf(fp, "  PRECFOCK  =  N           # HF FFT grid                                                     \n");
			fprintf(fp, "# NKRED     =  2           # Reduce k-grid-even only, see also NKREDX, NKREDY and NKREDZ     \n");
			fprintf(fp, "# HFLMAX    =  4           # HF cut-off: 4d, 6f                                              \n");
			fprintf(fp, "# LDIAG     = .TRUE.       # Diagnolise Eigenvalues                                          \n");
		}
		if (!strcmp(argv[i], "vdw")) {
			fprintf(fp, "! DFT-D3 Correction\n");
			fprintf(fp, "  IVDW      =  11          # DFT-D3 method of method with no damping                         \n");
		}
		if (!strcmp(argv[i], "sic")) {
			fprintf(fp, "! DFT+U Calculation\n");
			fprintf(fp, "  LDAU      = .TRUE.       # Activate DFT+U                                                  \n");
			fprintf(fp, "  LDATYPE   =  2           # Dudarev; only U-J matters                                       \n");
			fprintf(fp, "  LDAUL     =  2 -1        # Orbitals for each species                                       \n");
			fprintf(fp, "  LDAUU     =  2  0        # U for each species 		                                      \n");
			fprintf(fp, "  LDAUJ     =  0  0        # J for each species 		                                      \n");
			fprintf(fp, "  LMAXMIX   =  4           # Mixing cut-off; 4-d, 6-f                                        \n");
		}
		if (!strcmp(argv[i], "ecc")) {
			fprintf(fp, "! Elastic constants Calculation					                                          \n");
			fprintf(fp, "# NGZ       =  500         # FFT grid mesh density for nice charge/potential plots           \n");
			fprintf(fp, "  IBRION    =  6           # Determine the Hessian matrix                                    \n");
			fprintf(fp, "  NFREE     =  4           # How many displacements are used for each direction; 2-4         \n");
			fprintf(fp, "  ISIF      =  3           # Stress/relaxation: 3-Shape/Ions/V                               \n");
			fprintf(fp, "  NSW       =  1           # Max ionic steps                                                 \n");
			fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence, in eV                                   \n");
			//fprintf(fp, "PREC      =  Accurate      # High level                                                    \n");
		}
		if (!strcmp(argv[i], "bca")) {
			fprintf(fp, "! Bader Charge Analysis\n");
			fprintf(fp, "  LAECHG    = .TRUE.       # Write core charge into CHGCAR file                              \n");
			fprintf(fp, "  LCHARG    = .TRUE.       # Write CHGCAR file                                               \n");
		}
		if (!strcmp(argv[i], "elf")) {
			fprintf(fp, "! Electron Localization Function\n");
			fprintf(fp, "  ISTART    =  1           # Read existing wavefunction; if there                            \n");
			fprintf(fp, "  LELF      = .TRUE.       # Activate ELF                                                    \n");
		}
		if (!strcmp(argv[i], "fpm")) {
			fprintf(fp, "! Frozen Phonon Method\n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing                                               \n");
			fprintf(fp, "  SIGMA     =  0.01        # Smearing value in eV                                            \n");
			fprintf(fp, "  IBRION    = -1           # Ions are not moved                                              \n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence; in eV                                   \n");
			fprintf(fp, "  PREC      =  Accurate    # Precision level                                                 \n");
			fprintf(fp, "# ENCUT     =  520         # Cut-off energy for plane wave basis set, in eV                  \n");
			fprintf(fp, "  IALGO     =  38          # Davidson block iteration scheme                                 \n");
			fprintf(fp, "  LREAL     = .FALSE.      # Projection operators: false                                     \n");
			fprintf(fp, "  LWAVE     = .FLASE.      # Write WAVECAR or not                                            \n");
			fprintf(fp, "  LCHARG    = .FLASE.      # Write CHGCAR or not                                             \n");
			fprintf(fp, "  ADDGRID   = .TRUE.       # Increase grid; helps GGA convergence                            \n");
		}
		if (!strcmp(argv[i], "dfp")) {
			fprintf(fp, "! Density functional Perturbation\n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing                					              \n");
			fprintf(fp, "  SIGMA     =  0.05        # Smearing value in eV                                            \n");
			fprintf(fp, "  IBRION    =  8           # determines the Hessian matrix using DFPT                        \n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence; in eV                                   \n");
			fprintf(fp, "  PREC      =  Accurate    # Precision level                                                 \n");
			fprintf(fp, "# ENCUT     =  520         # Cut-off energy for plane wave basis set, in eV                  \n");
			fprintf(fp, "  IALGO     =  38          # Davidson block iteration scheme                                 \n");
			fprintf(fp, "  LREAL     = .FALSE.      # Projection operators: false                                     \n");
			fprintf(fp, "  LWAVE     = .FLASE.      # Write WAVECAR or not                                            \n");
			fprintf(fp, "  LCHARG    = .FLASE.      # Write CHGCAR or not                                             \n");
			fprintf(fp, "  ADDGRID   = .TRUE.       # Increase grid; helps GGA convergence                            \n");
		}
		if (!strcmp(argv[i], "neb")) {
			fprintf(fp, "! Budged Elastic Band (NEB)\n");
			fprintf(fp, "  IMAGES    =  5           # Defines the number of VASP calculations in separate directories \n");
			fprintf(fp, "  NSW       =  500         # Number of ionic steps                                           \n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
			fprintf(fp, "  SIGMA     =  0.05        # Please check the width of the smearing                          \n");
			fprintf(fp, "  IBRION    =  3           # Do MD with a zero time step                                     \n");
			fprintf(fp, "  POTIM     =  0           # Zero time step so that VASP does not move the ions              \n");
			fprintf(fp, "  SPRING    = -5.0         # Spring force (eV/A2) between images                             \n");
			fprintf(fp, "  LCLIMB    = .TRUE.       # Turn on the climbing image algorithm                            \n");
			fprintf(fp, "  ICHAIN    =  0           # Indicates which method to run. NEB (ICHAIN=0) is the default    \n");
			fprintf(fp, "  IOPT      =  1           # LBFGS = Limited-memory Broyden-Fletcher-Goldfarb-Shanno         \n");
		}
		if (!strcmp(argv[i], "tss") || !strcmp(argv[i], "dim")) {
			fprintf(fp, "! The dimer\n");
			fprintf(fp, "  NSW       =  500         # Number of ionic steps                                           \n");
			fprintf(fp, "  IBRION    =  3           # Do MD with a zero time step                                     \n");
			fprintf(fp, "  POTIM     =  0           # Zero time step so that VASP does not move the ions              \n");
			fprintf(fp, "  ICHAIN    =  2           # Use the dimer method required for the latest code               \n");
			fprintf(fp, "  DdR       =  0.005       # The dimer separation (twice the distance between images)        \n");
			fprintf(fp, "  DRotMax   =  1           # Maximum number of rotation steps per translation step           \n");
			fprintf(fp, "  DFNMin    =  0.01        # Magnitude of the rotational force below which the dimer is not rotated \n");
			fprintf(fp, "  DFNMax    =  1.0         # Magnitude of the rotational force below which dimer rotation stops     \n");
			fprintf(fp, "  IOPT      =  2           # CG = Conjugate Gradient                                         \n");
		}
		if (!strcmp(argv[i], "pbs")) {
			fprintf(fp, "! Static State Calculation\n");
			fprintf(fp, "  NSW       =  0           # Number of ionic steps                                           \n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
			fprintf(fp, "  SIGMA     =  0.1         # Please check the width of the smearing                          \n");
			fprintf(fp, "  IBRION    = -1           # Algorithm: 0-MD; 1-Quasi-New; 2-CG                              \n");
			fprintf(fp, "  LORBIT    =  11          # PAW radii for projected DOS                                     \n");
			fprintf(fp, "# NEDOS     =  2001        # DOSCAR points                                                   \n");
			fprintf(fp, "  NELM      =  60          # Max electronic SCF steps                                        \n");
			fprintf(fp, "\n");
			fprintf(fp, "! SCF energy convergence; in eV\n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence; in eV                                   \n");
			fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
		}
		if (!strcmp(argv[i], "bbs")) {
			fprintf(fp, "! Static State Calculation\n");
			fprintf(fp, "  NSW       =  0           # Number of ionic steps                                           \n");
			fprintf(fp, "  ISMEAR    =  0           # Gaussian smearing method                                        \n");
			fprintf(fp, "  SIGMA     =  0.1         # Please check the width of the smearing                          \n");
			fprintf(fp, "  IBRION    = -1           # Algorithm: 0-MD; 1-Quasi-New; 2-CG                              \n");
			fprintf(fp, "  LORBIT    =  11          # PAW radii for projected DOS                                     \n");
			fprintf(fp, "# NEDOS     =  2001        # DOSCAR points                                                   \n");
			fprintf(fp, "  NELM      =  60          # Max electronic SCF steps                                        \n");
			fprintf(fp, "\n");
			fprintf(fp, "! SCF energy convergence; in eV\n");
			fprintf(fp, "  EDIFF     =  1E-06       # SCF energy convergence; in eV                                   \n");
			fprintf(fp, "  EDIFFG    = -2E-2        # Ionic convergence; eV/AA                                        \n");
		}
		if (!strcmp(argv[i], "dos")) {
			fprintf(fp, "! Density of States\n");
			fprintf(fp, "  ISPIN     =  1           # 1: non-spin-polarized calculations are performed.               \n");
			fprintf(fp, "# ISPIN     =  2           # 2: spin-polarized calculations (collinear) are performed.       \n");
			fprintf(fp, "  IBRION    = -1                                                                             \n");
			fprintf(fp, "  NSW       =  0                                                                             \n");
			fprintf(fp, "  ISMEAR    = -5                                                                             \n");
			fprintf(fp, "  LORBIT    =  11                                                                            \n");
			fprintf(fp, "  NEDOS     =  1000                                                                          \n");
			fprintf(fp, "# ICHARG    =  11          # Non-selfconsistent calculations; superposition charge densities \n");
		}
		if (!strcmp(argv[i], "bcd")) {
			fprintf(fp, "! Basic Charge Density\n");
			fprintf(fp, "  LAECHG    = .TRUE.       # Write core charge into CHGCAR file                              \n");
			fprintf(fp, "  LCHARG    = .TRUE.       # Write CHGCAR file                                               \n");
		}
		if (!strcmp(argv[i], "scd")) {
			fprintf(fp, "! Spin Charge Density\n");
			fprintf(fp, "  ISPIN     =  2           # Spin polarised DFT                                              \n");
			fprintf(fp, "  LAECHG    = .TRUE.       # Write core charge into CHGCAR file                              \n");
			fprintf(fp, "  LCHARG    = .TRUE.       # Write CHGCAR file                                               \n");
		}
		if (!strcmp(argv[i], "pcd")) {
			fprintf(fp, "! Partial Charge Density\n");
			fprintf(fp, "  ISPIN     =  2           # Spin polarised DFT                                              \n");
			fprintf(fp, "  LAECHG    = .TRUE.       # Write core charge into CHGCAR file                              \n");
			fprintf(fp, "  LCHARG    = .TRUE.       # Write CHGCAR file                                               \n");
		}
		if (!strcmp(argv[i], "pcd_ik")) {
			fprintf(fp, "! Decomposed Charge Density\n");
			fprintf(fp, "  ISTART    =  1           # Job: 0 - new  1 - cont  2 - samecut                             \n");
			fprintf(fp, "  ICHARG    =  1           # Read charge : 1 - file 2 - atom 10 - const                      \n");
			fprintf(fp, "  LPARD     = .TRUE.       # Activate decomposed charge density                              \n");
			fprintf(fp, "  LSEPB     = .TRUE.       # Separately write PARCHG.nb by every band or not                 \n");
			fprintf(fp, "  LSEPK     = .TRUE.       # Separately write PARCHG.nk by every kpoint or not               \n");
			fprintf(fp, "\n");
			fprintf(fp, "! Method I : Partial Charge for the specified BANDSand KPOINTS\n");
			fprintf(fp, "  IBAND     =  20 21 22 23 # Set this parameters manually                                    \n");
			fprintf(fp, "  KPUSE     =  1 2 3 4     # Set this parameters manually                                    \n");
			fprintf(fp, "\n");
			fprintf(fp, "! ********************** Notes **********************\n");
			fprintf(fp, "! (1) Copy IBZKPT as KPOINTS for static calculation,\n");
			fprintf(fp, "! (2) Band structure calculation.\n");
		}
		if (!strcmp(argv[i], "pcd_en")) {
			fprintf(fp, "! Decomposed Charge Density\n");
			fprintf(fp, "  ISTART    =  1           # Job: 0 - new  1 - cont  2 - samecut                             \n");
			fprintf(fp, "  ICHARG    =  1           # Read charge : 1 - file 2 - atom 10 - const                      \n");
			fprintf(fp, "  LPARD     = .TRUE.       # Activate decomposed charge density                              \n");
			fprintf(fp, "  LSEPB     = .TRUE.       # Separately write PARCHG.nb by every band or not                 \n");
			fprintf(fp, "  LSEPK     = .TRUE.       # Separately write PARCHG.nk by every kpoint or not               \n");
			fprintf(fp, "\n");
			fprintf(fp, "! Method II : Partial Charge in the energy rang of[-10.3 - 5.1]\n");
			fprintf(fp, "  EINT      = -10.3 - 5.1  # Set this parameters manually                                    \n");
			fprintf(fp, "\n");
			fprintf(fp, "! ********************** Notes **********************\n");
			fprintf(fp, "! (1) Copy IBZKPT as KPOINTS for static calculation,\n");
			fprintf(fp, "! (2) Band structure calculation.\n");
		}
		if (!strcmp(argv[i], "pcd_ef")) {
			fprintf(fp, "! Decomposed Charge Density\n");
			fprintf(fp, "  ISTART    =  1           # Job: 0 - new  1 - cont  2 - samecut                             \n");
			fprintf(fp, "  ICHARG    =  1           # Read charge : 1 - file 2 - atom 10 - const                      \n");
			fprintf(fp, "  LPARD     = .TRUE.       # Activate decomposed charge density                              \n");
			fprintf(fp, "  LSEPB     = .TRUE.       # Separately write PARCHG.nb by every band or not                 \n");
			fprintf(fp, "  LSEPK     = .TRUE.       # Separately write PARCHG.nk by every kpoint or not               \n");
			fprintf(fp, "\n");
			fprintf(fp, "! Method III : Partial Charge in the energy rang of[EF - 1 - EF]\n");
			fprintf(fp, "  NBMOD     = -3                                \n");
			fprintf(fp, "  EINT      = -1           # Set this parameters manually                                    \n");
			fprintf(fp, "\n");
			fprintf(fp, "! ********************** Notes **********************\n");
			fprintf(fp, "! (1) Copy IBZKPT as KPOINTS for static calculation,\n");
			fprintf(fp, "! (2) Band structure calculation.\n");
		}
		if (!strcmp(argv[i], "wfn")) {
			fprintf(fp, "! Wave Function Calculation\n");
			fprintf(fp, "  LAECHG    = .TRUE.       # Write core charge into CHGCAR file                              \n");
			fprintf(fp, "  LCHARG    = .TRUE.       # Write CHGCAR file                                               \n");
		}
		if (!strcmp(argv[i], "tlp")) {
			fprintf(fp, "! Total Local Potential\n");
			fprintf(fp, "  NPAR      =  4           # Number of compute cores that work on an individual orbital      \n");
			fprintf(fp, "  ISTART    =  1           # Job: 0 - new  1 - cont  2 - samecut                             \n");
			fprintf(fp, "  ICHARG    =  1           # How VASP constructs the initial charge density                  \n");
			fprintf(fp, "  LWAVE     = .TRUE.       # Whether the wavefunctions are written to the WAVECAR file       \n");
			fprintf(fp, "  LCHARG    = .TRUE.       # LCHARG determines whether the charge densities(files CHGCAR and CHG) are written\n");
			fprintf(fp, "  LVTOT     = .TRUE.       # Whether the total local potential is written to the LOCPOT file \n");
			fprintf(fp, "  LVHAR     = .FALSE.      # Whether the electrostatic potential is written to the LOCPOT file\n");
			fprintf(fp, "  LELF      = .FALSE.      # Whether to create an ELFCAR file or not                         \n");
			fprintf(fp, "  ENCUT     =  400         # Cutoff energy for the planewave basis set in eV                 \n");
			fprintf(fp, "  ISMEAR    =  1           # How the partial occupancies are set for each orbital.           \n");
			fprintf(fp, "  SIGMA     =  0.2         # The width of the smearing in eV                                 \n");
			fprintf(fp, "  EDIFF     =  1E-6        # Global break condition for the electronic SC - loop             \n");
			fprintf(fp, "  NELMIN    =  5           # Specifies the minimum number of electronic SCF steps            \n");
			fprintf(fp, "  NELM      =  300         # Maximum number of electronic SC(selfconsistency) steps          \n");
			fprintf(fp, "  GGA       =  PE          # Type of generalized - gradient - approximation one wishes to use\n");
			fprintf(fp, "  LREAL     =  Auto        # Whether the projection operators are evaluated in real - space or in reciprocal space\n");
			fprintf(fp, "  ICORELEVEL=  2                                                                             \n");
			fprintf(fp, "  CLNT      =  1           # Species                                                         \n");
			fprintf(fp, "  CLN       =  1           # Main quantum number of excited core electron                    \n");
			fprintf(fp, "  CLL       =  0                                                                             \n");
			fprintf(fp, "  CLZ       =  0.5                                                                           \n");
		}
		if (!strcmp(argv[i], "tdm"))
		{
			fprintf(fp, "! Transition Dipole Moment\n");
			fprintf(fp, "  ISTART    =  1			# Job: 0 - new  1 - cont  2 - samecut                             \n");
			fprintf(fp, "  LREAL     = .FALSE.		# Whether the projection operators are evaluated in real - space or in reciprocal space\n");
			fprintf(fp, "  PREC      =  Normal		# Precision level                                                 \n");
			fprintf(fp, "  LWAVE     = .TRUE.		# Whether the wavefunctions are written to the WAVECAR file       \n");
			fprintf(fp, "  LCHARG    = .TRUE.		# LCHARG determines whether the charge densities(files CHGCAR and CHG) are written\n");
			fprintf(fp, "  ADDGRID   = .TRUE.	    # Increase grid; helps GGA convergence                            \n");
			fprintf(fp, "  ENCUT     =  400			# Cutoff energy for the planewave basis set in eV                 \n");
			fprintf(fp, "  ICHARG    =  11			# How VASP constructs the initial charge density                  \n");
			fprintf(fp, "\n");
			fprintf(fp, "  ISMEAR    =  0			# How the partial occupancies are set for each orbital.           \n");
			fprintf(fp, "  SIGMA     =  0.05		# The width of the smearing in eV                                 \n");
			fprintf(fp, "  NELM      =  60			# Maximum number of electronic SC(selfconsistency) steps          \n");
			fprintf(fp, "  NELMIN    =  6			# Specifies the minimum number of electronic SCF steps            \n");
			fprintf(fp, "  EDIFF     =  1E-07		# Global break condition for the electronic SC - loop             \n");
		}
		if (!strcmp(argv[i], "los"))
		{
			fprintf(fp, "! Linear Optical Spectrums\n");
			fprintf(fp, "  PREC      =  Normal		# Precision level                                                 \n");
			fprintf(fp, "  ENCUT     =  250.0		# Cutoff energy for the planewave basis set in eV                 \n");
			fprintf(fp, "  ALGO      =  BSE		    # Electronic Minimisation Algorithm                               \n");
			fprintf(fp, "  ISMEAR    =  0			# How the partial occupancies are set for each orbital            \n");
			fprintf(fp, "  SIGMA     =  0.05		# The width of the smearing in eV                                 \n");
			fprintf(fp, "  NBANDS    =  	        # Set this parameters manually, number of bands to include        \n");
			fprintf(fp, "  NOMEGA    =  72		    # Specifies the number of (imaginary) frequency and imaginary time grid points\n");
			fprintf(fp, "  LSPECTRAL = .TRUE.       # Specifies to use the spectral method                            \n");
			fprintf(fp, "  LOPTICS   = .TRUE.       # Switch on optic calculations                                    \n");
			fprintf(fp, "  OMEGAMAX  =  10	        # Specifies the maximum frequency for the dense part of the frequency grid\n");
			fprintf(fp, "  NBANDSO   =  16	        # Determines how many occupied orbitals are included              \n");
			fprintf(fp, "  NBANDSV   =  16	        # Determines how many unoccupied orbitals are included            \n");
			fprintf(fp, "  NEDOS     =  3000	    # Number of grid points in DOS                                    \n");
		}
		if (!strcmp(argv[i], "esp"))
		{
			fprintf(fp, "! Electrostatic Potential\n");
			fprintf(fp, "  LWAVE     = .TRUE.		# Whether the wavefunctions are written to the WAVECAR file       \n");
			fprintf(fp, "  LCHARG    = .TRUE.		# LCHARG determines whether the charge densities(files CHGCAR and CHG) are written\n");
			fprintf(fp, "  LVHAR     = .TRUE.       # Whether the electrostatic potential is written to the LOCPOT file\n");
			fprintf(fp, "  LELF      = .FALSE.      # Whether to create an ELFCAR file or not                         \n");
			fprintf(fp, "  LORBIT    =  11          # PAW radii for projected DOS                                     \n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("Written %s file!\n", filename.c_str());
	copy_newinp(filename.c_str(), "NEWINP");
	return filename;
}

void INCAR_fix(char name[], vector<string> value)
{
	string fix = "";
	for (int i = 0; i < value.size(); i++)
		fix += i == value.size() - 1 ? value[i] : value[i] + " ";
	FILE* fp1 = fopen("INCAR", "r");
	FILE* fp2 = fopen("_remove_Buf", "w");
	char newline[1024];
	char buf1[1024];
	char buf2[1024];
	int exi = 0;
	if (fp1 == NULL)
	{
		printf("INCAR IS NOT EXIST!\n");
		return;
	}
	while (fgets(buf1, 1024, fp1) != NULL)
	{
		char tmp[100];
		for (int i = 0; i < strlen(buf1); i++)
			tmp[i] = buf1[i];
		char* key = strtok(tmp, "=");
		key = strtok(key, " ");
		if (strcmp(key, name))
			fputs(buf1, fp2);
		else
		{
			sprintf(newline, "%s = %s\n", name, fix.c_str());
			fputs(newline, fp2);
			exi = 1;
		}
	}
	if (!exi)
	{
		sprintf(newline, "%s = %s\n", name, fix.c_str());
		fputs(newline, fp2);
	}
	fclose(fp1);
	fclose(fp2);
	FILE* fp3 = fopen("INCAR", "w");
	FILE* fp4 = fopen("_remove_Buf", "r");
	while (fgets(buf2, 1024, fp4) != NULL)
	{
		fputs(buf2, fp3);
	}
	fclose(fp3);
	fclose(fp4);
	remove("_remove_Buf");
}

void INCAR_replace(int argc, char* argv[])
{
	if (argc < 4)
		return;
	int pos = 2;
	int flag = 0;
	vector<vector<string> >value;
	vector<string> tmp;
	vector<char*> keyword;
	//VASPMATE -i_rep A 1 2 3 | B c d e |...
	while (pos <= argc)
	{
		if (pos == argc || (flag && !strcmp(argv[pos], "-")) || !strcmp(argv[pos], "/"))
		{
			flag = 0;
			value.push_back(tmp);
			tmp.clear();
			pos++;
			continue;
		}
		if (!flag)
		{
			keyword.push_back(argv[pos]);
			flag = 1;
		}
		else
			tmp.push_back(argv[pos]);
		pos++;
	}
	if (keyword.size() != value.size())
	{
		printf("Input parameter may be error!\n");
		return;
	}
	for (int i = 0; i < keyword.size(); i++)
		INCAR_fix(keyword[i], value[i]);
}

void INCAR_delete(int argc, char* argv[])
{
	vector<string> keyword;
	for (int i = 2; i < argc; i++)
		keyword.push_back(argv[i]);
	FILE* fp1 = fopen("INCAR", "r");
	if (fp1 == NULL)
	{
		printf("No INCAR file is found!\n");
		return;
	}
	FILE* fp2 = fopen("_remove_Buf", "w");
	char buf1[1024];
	char buf2[1024];
	while (fgets(buf1, 1024, fp1) != NULL)
	{
		int flag = 0;
		for (int i = 0; i < keyword.size(); i++)
			if (strstr(buf1, keyword[i].c_str()) != NULL)
				flag = 1;
		if (!flag)
			fputs(buf1, fp2);
	}
	fclose(fp1);
	fclose(fp2);
	FILE* fp3 = fopen("INCAR", "w");
	FILE* fp4 = fopen("_remove_Buf", "r");
	while (fgets(buf2, 1024, fp4) != NULL)
	{
		fputs(buf2, fp3);
	}
	fclose(fp3);
	fclose(fp4);
	remove("_remove_Buf");
}

void INCAR_remove(int argc, char* argv[])
{
	vector<string> keyword;
	for (int i = 2; i < argc; i++)
		keyword.push_back(argv[i]);
	FILE* fp1 = fopen("INCAR", "r");
	if (fp1 == NULL)
	{
		printf("No INCAR file is found!\n");
		return;
	}
	FILE* fp2 = fopen("_remove_Buf", "w");
	char newline[50];
	char buf1[1024];
	char buf2[1024];
	while (fgets(buf1, 1024, fp1) != NULL)
	{
		int flag = 0;
		for (int i = 0; i < keyword.size(); i++)
			if (strstr(buf1, keyword[i].c_str()) != NULL)
				flag = 1;
		if (flag)
		{
			snprintf(newline, sizeof(newline),"#%s", buf1);
			fputs(newline, fp2);
			fputs("\n", fp2);
		}
		else
			fputs(buf1, fp2);
	}
	fclose(fp1);
	fclose(fp2);
	FILE* fp3 = fopen("INCAR", "w");
	FILE* fp4 = fopen("_remove_Buf", "r");
	while (fgets(buf2, 1024, fp4) != NULL)
	{
		fputs(buf2, fp3);
	}
	fclose(fp3);
	fclose(fp4);
	remove("_remove_Buf");
}

void INCAR_append(int argc, char* argv[])
{
	int pos = 2;
	int flag = 0;
	vector<vector<string> >value;
	vector<string> tmp;
	vector<char*> keyword;
	//VASPMATE -i_rep A 1 2 3 | B c d e |...
	while (pos <= argc)
	{
		if (pos == argc || (flag && !strcmp(argv[pos], "-")) || !strcmp(argv[pos], "/"))
		{
			flag = 0;
			value.push_back(tmp);
			tmp.clear();
			pos++;
			continue;
		}
		if (!flag)
		{
			keyword.push_back(argv[pos]);
			flag = 1;
		}
		else
			tmp.push_back(argv[pos]);
		pos++;
	}
	if (keyword.size() != value.size())
	{
		printf("Input parameter may be error!\n");
		return;
	}
	FILE* fp = fopen("INCAR", "at+");
	for (int i = 0; i < keyword.size(); i++)
	{
		fprintf(fp, "  %s = ", keyword[i]);
		for (int j = 0; j < value[i].size(); j++)
		{
			fprintf(fp, "  %s ", value[i][j].c_str());
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void pcd_model(int argc, char* argv[])
{
	if (!strcmp("-ik", argv[2]))
	{
		if (access("INCAR", 0) == -1 || GetInfoINCAR("KPUSE").size() == 0)
		{
			char* list[] = { NULL ,NULL,"pcd_ik" };
			string filename = write_INCAR(3, list);
			copy_file("INCAR", filename.c_str());
		}
		vector<string> ik;
		int position = 0;
		for (int i = 3; i < argc; i++)
		{
			if (!strcmp("-ib", argv[i]))
			{
				position = i;
				break;
			}
			ik.push_back(argv[i]);
		}
		INCAR_fix("KPUSE", ik);
		if (position)
		{
			vector<string> ib;
			for (int i = position + 1; i < argc; i++)
				ib.push_back(argv[i]);
			INCAR_fix("IBAND", ib);
		}
	}
	if (!strcmp("-ib", argv[2]))
	{
		if (access("INCAR", 0) == -1 || GetInfoINCAR("IBAND").size() == 0)
		{
			char* list[] = { NULL ,NULL,"pcd_ik" };
			string filename = write_INCAR(3, list);
			copy_file("INCAR", filename.c_str());
		}
		vector<string> ib;
		int position = 0;
		for (int i = 3; i < argc; i++)
		{
			if (!strcmp("-ik", argv[i]))
			{
				position = i;
				break;
			}
			ib.push_back(argv[i]);
		}
		INCAR_fix("IBAND", ib);
		if (position)
		{
			vector<string> ik;
			for (int i = position + 1; i < argc; i++)
				ik.push_back(argv[i]);
			INCAR_fix("KPUSE", ik);
		}
	}
	if (!strcmp("-en", argv[2]))
	{
		char* list[] = { NULL ,NULL,"pcd_en" };
		string filename = write_INCAR(3, list);
		copy_file("INCAR", filename.c_str());
		if (GetInfoINCAR("NBMOD").size() != 0)
			printf("The current mode used may not be method II,please check!\n");
		else if (argc < 5)
			printf("You should set a Energy range(two values)\n");
		vector<string> en;
		for (int i = 3; i < argc; i++)
			en.push_back(argv[i]);
		INCAR_fix("EINT", en);
	}
	if (!strcmp("-ef", argv[2]))
	{
		char* list[] = { NULL ,NULL,"pcd_ef" };
		string filename = write_INCAR(3, list);
		copy_file("INCAR", filename.c_str());
		if (GetInfoINCAR("NBMOD") != "-3")
			printf("The current mode used may not be method III,please check!\n");
		else if (argc > 4)
			printf("You should set only one energy value\n");
		vector<string> ef;
		for (int i = 3; i < argc; i++)
			ef.push_back(argv[i]);
		INCAR_fix("EINT", ef);
	}
}

int addIVDWdefault(const char file[], string val)
{
	FILE* fp = fopen(file, "r");
	if (fp == nullptr)
	{
		printf("%s IS NOT EXIST!\n", file);
		return -1;
	}
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	bool p[3] = { 0 };
	has_vacuum_slab_or_not(pos, p);
	if (p[1] == 1 || p[1] == 1 || p[1] == 1)
	{
		copy_newinp("INCAR","temINP");
		INCAR_fix("IVDW", { val });
		copy_newinp("INCAR","NEWINP");
		copy_newinp("temINP","INCAR");
		remove("temINP");		
		printf("The structure may be 2D system! Add IVDW = 12 automatically!\n");
		return 1;
	}
	else
	{
		printf("The structure may be not 2D system!\n");
		return 0;
	}
}

void addIVDW_table(int ivdw_numb)
{
	copy_newinp("INCAR","temINP");
	if (ivdw_numb == 0)
		printf("Please enter the parameter after '-t' !\n");
	else if (ivdw_numb == 1) //DFT-D3(BJ)
	{
		INCAR_fix("IVDW", { "12" });
		INCAR_fix("VDW_S6", { "1.0" });
		INCAR_fix("VDW_S8", { "0.7875" });
		INCAR_fix("VDW_A1", { "0.4289" });
		INCAR_fix("VDW_A2", { "4.4407" });
		printf("ADD DFT-D3(BJ) automatically!\n");
	}
	else if (ivdw_numb == 2) //DFT-D2(G)
	{
		INCAR_fix("IVDW", { "1" });
		INCAR_fix("LVDW", { ".TRUE." });
		INCAR_fix("VDW_S6", { "0.75" });
		INCAR_fix("VDW_SCALING", { "0.75" });
		printf("ADD DFT-D2(G) automatically!\n");
	}
	else if (ivdw_numb == 3) //DFT-D3(zero)
	{
		INCAR_fix("IVDW", { "11" });
		INCAR_fix("VDW_S6", { "1.0" });
		INCAR_fix("VDW_S8", { "0.722" });
		INCAR_fix("VDW_SR", { "1.217" });
		printf("ADD DFT-D3(zero) automatically!\n");
	}
	else if (ivdw_numb == 4) //optB86b-vdw
	{
		INCAR_fix("LUSE_VDW", { ".TRUE." });
		INCAR_fix("LASPH", { ".TRUE." });
		INCAR_fix("AGGAC", { "0.0000" });
		INCAR_fix("GGA", { "MK" });
		INCAR_fix("PARAM1", { "0.1234" });
		INCAR_fix("PARAM2", { "1.0000" });
		printf("ADD optB86b-vdw automatically!\n");
	}
	else if (ivdw_numb == 5) //optB88-vdw
	{
		INCAR_fix("LUSE_VDW", { ".TRUE." });
		INCAR_fix("LASPH", { ".TRUE." });
		INCAR_fix("AGGAC", { "0.0000" });
		INCAR_fix("GGA", { "BO" });
		INCAR_fix("PARAM1", { "0.1833333333" });
		INCAR_fix("PARAM2", { "0.2200000000" });
		printf("ADD optB88-vdw automatically!\n");
	}
	else if (ivdw_numb == 6) //optPBE-vdw
	{
		INCAR_fix("LUSE_VDW", { ".TRUE." });
		INCAR_fix("LASPH", { ".TRUE." });
		INCAR_fix("AGGAC", { "0.0000" });
		INCAR_fix("GGA", { "OR" });
		printf("ADD optPBE-vdw automatically!\n");
	}
	else if (ivdw_numb == 7) //revPBE-vdw
	{
		INCAR_fix("LUSE_VDW", { ".TRUE." });
		INCAR_fix("LASPH", { ".TRUE." });
		INCAR_fix("AGGAC", { "0.0000" });
		INCAR_fix("GGA", { "RE" });
		printf("ADD revPBE-vdw automatically!\n");
	}
	copy_newinp("INCAR","NEWINP");
	copy_newinp("temINP","INCAR");
	remove("temINP");
}

map<string, LDAU_DATA> aflow_ldau
{
{"H",{"0","0","0"}},
{"He",{"0","0","0"}},
{"Li",{"0","0","0"}},
{"Be",{"0","0","0"}},
{"B",{"0","0","0"}},
{"C",{"0","0","0"}},
{"N",{"0","0","0"}},
{"O",{"0","0","0"}},
{"F",{"0","0","0"}},
{"Ne",{"0","0","0"}},
{"Na",{"0","0","0"}},
{"Mg",{"0","0","0"}},
{"Al",{"0","0","0"}},
{"Si",{"0","0","0"}},
{"P",{"0","0","0"}},
{"S",{"0","0","0"}},
{"Cl",{"0","0","0"}},
{"Ar",{"0","0","0"}},
{"K",{"0","0","0"}},
{"Ca",{"0","0","0"}},
{"Sc",{"2","2.9","0"}},
{"Ti",{"2","4.4","0"}},
{"V", {"2","2.7","0"}},
{"Cr",{"2","3.5","0"}},
{"Mn",{"2","4","0"}},
{"Fe",{"2","4.6","0"}},
{"Co",{"2","5","0"}},
{"Ni",{"2","5.1","0"}},
{"Cu",{"2","4","0"}},
{"Zn",{"2","7.5","0"}},
{"Ga",{"2","3.9","0"}},
{"Ge",{"0","0","0"}},
{"As",{"0","0","0"}},
{"Se",{"0","0","0"}},
{"Br",{"0","0","0"}},
{"Kr",{"0","0","0"}},
{"Rb",{"0","0","0"}},
{"Sr",{"0","0","0"}},
{"Y", {"0","0","0"}},
{"Zr",{"0","0","0"}},
{"Nb",{"2","2.1","0"}},
{"Mo",{"2","2.4","0"}},
{"Tc",{"2","2.7","0"}},
{"Ru",{"2","3","0"}},
{"Rh",{"2","3.3","0"}},
{"Pd",{"2","3.6","0"}},
{"Ag",{"2","5.8","0"}},
{"Cd",{"2","2.1","0"}},
{"In",{"2","1.9","0"}},
{"Sn",{"2","3.5","0"}},
{"Sb",{"0","0","0"}},
{"Te",{"0","0","0"}},
{"I", {"0","0","0"}},
{"Xe",{"0","0","0"}},
{"Cs",{"0","0","0"}},
{"Ba",{"0","0","0"}},
{"La",{"3","8.1","0.6"}},
{"Ce",{"3","7","0.7"}},
{"Pr",{"3","6.5","1"}},
{"Nd",{"3","7.2","1"}},
{"Pm",{"0","0","0"}},
{"Sm",{"3","7.4","1"}},
{"Eu",{"3","6.4","1"}},
{"Gd",{"3","6.7","0.7"}},
{"Tb",{"0","0","0"}},
{"Dy",{"3","5.6","0"}},
{"Ho",{"0","0","0"}},
{"Er",{"0","0","0"}},
{"Tm",{"3","7","1"}},
{"Yb",{"3","7","0.7"}},
{"Lu",{"3","4.8","0.9"}},
{"Hf",{"0","0","0"}},
{"Ta",{"2","2","0"}},
{"W", {"2","2.2","0"}},
{"Re",{"2","2.4","0"}},
{"Os",{"2","2.6","0"}},
{"Ir",{"2","2.8","0"}},
{"Pt",{"2","3","0"}},
{"Au",{"2","4","0"}},
{"Hg",{"0","0","0"}},
{"Tl",{"0","0","0"}},
{"Pb",{"1","0","0"}},
{"Bi",{"1","0","0"}},
{"Po",{        }},
{"At",{        }},
{"Rn",{        }},
{"Fr",{        }},
{"Ra",{        }},
{"Ac",{        }},
{"Th",{"3","5","0"}},
{"Pa",{"0","0","0"}},
{"U",{"3","4","0"}},
{"Np",{        }},
{"Pu",{"0","0","0"}},
{"Am",{        }},
{"Cm",{        }},
{"Bk",{        }},
{"Cf",{        }},
{"Es",{        }},
{"Fm",{        }},
{"Md",{        }},
{"No",{        }},
{"Lr",{        }},
{"Rf",{        }},
{"Db",{        }},
{"Sg",{        }},
{"Bh",{        }},
{"Hs",{        }},
{"Mt",{        }}
};

LDAU::LDAU(const char table[], const char file[])
{
	LDAUTYPE = 2; // default value
	LMAXMIX = 2; // default value
	FILE* fp = fopen(file, "r");
	if (fp == nullptr)
	{
		printf("No input file is found!\n");
		return;
	}
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	//fix ldau value according to table
	if (table != nullptr)
	{
		FILE* fpt = fopen(table, "r");
		if (fpt == nullptr)
		{
			printf("%s IS NOT EXIST!\n", table);
			printf("Use default LDAU value!\n");
		}
		else
		{
			char buf[1024];
			while (fgets(buf, 1024, fpt))
			{
				char elem[10], L[10], U[10], J[10];
				if (sscanf(buf, "%s%s%s%s", elem, L, U, J) != 4)
				{
					remove(std::begin(buf), std::end(buf), '\n');
					printf("[%s] in [%s] file may be wrong format!\n", buf, table);
					continue;
				}
				if (!aflow_ldau.count(elem))
				{
					printf("elem [%s] in [%s] file may be wrong elemnt!\n", elem, table);
					continue;
				}
				else
				{
					aflow_ldau[elem].J = J;
					aflow_ldau[elem].L = L;
					aflow_ldau[elem].U = U;
				}
			}
		}
		fclose(fpt);
	}
	data.resize(pos.nant[1]);
	for (int i = 0; i < pos.nant[1]; i++)
	{
		if (aflow_ldau.count(pos.elemsym[i]))
		{
			data[i] = aflow_ldau[pos.elemsym[i]];
		}
		else
		{
			printf("Error! %s is not an element type!\n", pos.elemsym[i]);
			return;
		}
	}
	for (int i = 0; i < data.size(); i++)
		LMAXMIX = max(LMAXMIX, 2 * atoi(data[i].L.c_str()));
}

void LDAU::AddLDAU_default(const char file[])
{
	FILE* fp1 = fopen("INCAR", "r");
	if (fp1 == NULL)
	{
		printf("INCAR IS NOT EXIST!\n");
		return;
	}
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		return;
	}
	fclose(fp);
	copy_newinp("INCAR","temINP");
	INCAR_fix("LDAU", { ".TRUE." });
	INCAR_fix("LDAUTYPE", { to_string(LDAUTYPE) });
	vector<string> l(data.size()), u(data.size()), j(data.size());
	for (int i = 0; i < data.size(); i++)
	{
		l[i] = data[i].L;
		u[i] = data[i].U;
		j[i] = data[i].J;
	}
	INCAR_fix("LDAUL", l);
	INCAR_fix("LDAUU", u);
	INCAR_fix("LDAUJ", j);
	INCAR_fix("LMAXMIX", { to_string(LMAXMIX) });
	copy_newinp("INCAR","NEWINP");
	copy_newinp("temINP","INCAR");
	remove("temINP");
	fclose(fp1);
	printf("Add LDAU value to NEWINP automatically!\n");
}

map<string, double> aflow_mag
{
{"H",0.001},
{"He",0},
{"Li",0.001},
{"Be",0},
{"B",0.001},
{"C",0.004},
{"N",-0.012},
{"O",0.328},
{"F",0.107},
{"Ne",0},
{"Na",0},
{"Mg",0.005},
{"Al",0.002},
{"Si",-0.013},
{"P",0.007},
{"S",-0.002},
{"Cl",0.03},
{"Ar",0},
{"K",0},
{"Ca",-0.031},
{"Sc",-0.147},
{"Ti",1.345},
{"V",1.047},
{"Cr",3.199},
{"Mn",4.297},
{"Fe",1.15},
{"Co",1.816},
{"Ni",1.435},
{"Cu",0.732},
{"Zn",0.002},
{"Ga",0.013},
{"Ge",-0.013},
{"As",0.024},
{"Se",-0.205},
{"Br",-0.004},
{"Kr",0},
{"Rb",-0.001},
{"Sr",0.009},
{"Y",-0.069},
{"Zr",0.075},
{"Nb",-0.425},
{"Mo",0.252},
{"Tc",2.491},
{"Ru",0.422},
{"Rh",-0.048},
{"Pd",1.667},
{"Ag",-0.005},
{"Cd",0.004},
{"In",-0.017},
{"Sn",0.019},
{"Sb",0.002},
{"Te",-0.002},
{"I",0.002},
{"Xe",0.107},
{"Cs",-0.001},
{"Ba",0.002},
{"La",0.15},
{"Ce",0.981},
{"Pr",2.002},
{"Nd",3.023},
{"Pm",4.102},
{"Sm",5.581},
{"Eu",6.927},
{"Gd",7.39},
{"Tb",-0.099},
{"Dy",-0.113},
{"Ho",-0.123},
{"Er",-0.133},
{"Tm",1.021},
{"Yb",0.796},
{"Lu",-0.177},
{"Hf",0.021},
{"Ta",-0.431},
{"W",-0.537},
{"Re",0.894},
{"Os",-0.112},
{"Ir",-0.136},
{"Pt",0.025},
{"Au",0.002},
{"Hg",0.001},
{"Tl",0.006},
{"Pb",-0.015},
{"Bi",-0.03},
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

Mag::Mag(bool noncol, const char table[]): noncolliner(noncol)
{
	FILE* fp = fopen("INPOS", "r");
	if (fp == nullptr)
	{
		printf("INPOS IS NOT EXIST!\n");
		return;
	}
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	//fix magmom value according to table
	if (table != nullptr)
	{
		FILE* fpt = fopen(table, "r");
		if (fpt == nullptr)
		{
			printf("%s IS NOT EXIST!\n", table);
			printf("Use default MAG value!\n");
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
				if (!aflow_mag.count(elem))
				{
					printf("elem [%s] in [%s] file may be wrong elemnt!\n", elem, table);
					continue;
				}
				else
				{
					aflow_mag[elem] = mag;
				}
			}
		}
		fclose(fpt);
	}
	data.resize(pos.nant[0]);
	int cnt = 0;
	for (int i = 0; i < pos.nant[1]; i++)
	{
		if (aflow_mag.count(pos.elemsym[i]))
		{
			for (int j = 0; j < pos.typenum[i]; j++)
			{
				if(!noncolliner)
					data[cnt++] = to_string(aflow_mag[pos.elemsym[i]]);
				else
					data[cnt++] = "3*" + to_string(aflow_mag[pos.elemsym[i]]);
			}
		}
		else
		{
			printf("Error! %s is not an element type!\n", pos.elemsym[i]);
			return;
		}
	}
}

void Mag::AddMag_default(bool noncol)
{
	FILE* fp1 = fopen("INCAR", "r");
	if (fp1 == NULL)
	{
		printf("INCAR IS NOT EXIST!\n");
		return;
	}
	else
	{
		copy_newinp("INCAR","temINP");
		INCAR_fix("ISPIN", { "2" });
		if (noncol == 1)
			INCAR_fix("LNONCOLLINE", { ".TURE." });
		else if (noncol == 0)
			INCAR_fix("LNONCOLLINE", { ".FALSE." });
		INCAR_fix("MAGMON", data);
		copy_newinp("INCAR","NEWINP");
		copy_newinp("temINP","INCAR");
		remove("temINP");
		printf("Add magnetic moment to INCAR automatically!\n");
		fclose(fp1);
	}
}
