#include"incar.h"
#include"tools.h"
using namespace std;

void print_global(FILE* fp)
{
	fprintf(fp, "Global Parameters\n");
	fprintf(fp, " ISTART =  0            (Read existing wavefunction; if there)\n");
	fprintf(fp, " # ISPIN =  1            (Spin polarised DFT)\n");
	fprintf(fp, " ICHARG =  2         (Non-self-consistent: GGA/LDA band structures)\n");
	fprintf(fp, " LREAL  = .FALSE.       (Projection operators: automatic)\n");
	fprintf(fp, " ENCUT  =  520        (Cut-off energy for plane wave basis set, in eV)\n");
	fprintf(fp, " PREC   =  Accurate       (Precision level)\n");
	fprintf(fp, " #LWAVE  = .TRUE.        (Write WAVECAR or not)\n");
	fprintf(fp, " #LCHARG = .TRUE.        (Write CHGCAR or not)\n");
	fprintf(fp, " ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)\n");
	fprintf(fp, " # LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)\n");
	fprintf(fp, " # LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)\n");
	fprintf(fp, " # NELECT =             (No. of electrons: charged cells; be careful)\n");
	fprintf(fp, " # LPLANE = .TRUE.      (Real space distribution; supercells)\n");
	fprintf(fp, " #NPAR   = 4           (Max is no. nodes; don't set for hybrids)\n");
	fprintf(fp, " # NWRITE = 2           (Medium-level output)\n");
	fprintf(fp, " # KPAR   = 2           (Divides k-grid into separate groups)\n");
	fprintf(fp, " # NGX    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
	fprintf(fp, " # NGY    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
	fprintf(fp, " # NGZ    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
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
		if (!strcmp(argv[i], "stc") || !strcmp(argv[i], "sst")) {
			print_global(fp);
			fprintf(fp, "#Static State Calculation							     			\n");
			fprintf(fp, "NSW    =  0          (number of ionic steps)                             \n");
			fprintf(fp, "ISMEAR =  0            (gaussian smearing method)              \n");
			fprintf(fp, "SIGMA  =  0.1         (please check the width of the smearing)\n");
			fprintf(fp, "IBRION =  -1            (Algorithm: 0-MD; 1-Quasi-New; 2-CG)                \n");
			fprintf(fp, "#LORBIT =  11           (PAW radii for projected DOS)           \n");
			fprintf(fp, "#NEDOS  =  2001         (DOSCAR points)                         \n");
			fprintf(fp, "NELM   =  60           (Max electronic SCF steps)              \n");
			fprintf(fp, "        (SCF energy convergence; in eV)         \n");
			fprintf(fp, "EDIFF  =  1E-06                                                           \n");
			fprintf(fp, "EDIFFG = -2E-2      (Ionic convergence; eV/AA)                          \n");
		}
		if (!strcmp(argv[i], "rlx") || !strcmp(argv[i], "lar")) {
			print_global(fp);
			fprintf(fp, "#Lattice Atomic Relaxtion													\n");
			fprintf(fp, "ISMEAR =  0            (gaussian smearing method)              \n");
			fprintf(fp, "SIGMA  =  0.1          (please check the width of the smearing)\n");
			fprintf(fp, "NSW    =  300          (number of ionic steps)                             \n");
			fprintf(fp, "ISMEAR =  0            (gaussian smearing method )                         \n");
			fprintf(fp, "SIGMA  =  0.05         (please check the width of the smearing)            \n");
			fprintf(fp, "IBRION =  2            (Algorithm: 0-MD; 1-Quasi-New; 2-CG)                \n");
			fprintf(fp, "ISIF   =  3            (optimize atomic coordinates and lattice parameters)\n");
			fprintf(fp, "        (SCF energy convergence; in eV)         \n");
			fprintf(fp, "EDIFFG = -2E-2         (Ionic convergence; eV/AA)                          \n");
			fprintf(fp, "EDIFF  =  1E-06                                                           \n");
			fprintf(fp, "PREC   =  Accurate     (Precision level)                                   \n");
		}
		if (!strcmp(argv[i], "mds")) {
			print_global(fp);
			fprintf(fp, "#Electronic Relaxation																 \n");
			fprintf(fp, " ISMEAR =  0                                                                            \n");
			fprintf(fp, " SIGMA  =  0.05                                                                         \n");
			fprintf(fp, " EDIFF  =  1E-06                                                                        \n");
			fprintf(fp, "                                                                                        \n");
			fprintf(fp, "#Molecular Dynamics                                                                     \n");
			fprintf(fp, " IBRION =  0            (Activate MD)                                                   \n");
			fprintf(fp, " NSW    =  100          (Max ionic steps)                                               \n");
			fprintf(fp, " EDIFFG = -1E-02        (Ionic convergence; eV/A)                                       \n");
			fprintf(fp, " POTIM  =  1            (Timestep in fs)                                                \n");
			fprintf(fp, " SMASS  =  0            (MD Algorithm: -3-microcanonical ensemble; 0-canonical ensemble)\n");
			fprintf(fp, " ! TEBEG  =     100     (Start temperature K)                                           \n");
			fprintf(fp, " ! TEEND  =     100     (Final temperature K)                                           \n");
			fprintf(fp, " ! MDALGO =  1          (Andersen Thermostat)                                           \n");
			fprintf(fp, " ! ISYM   =  0          (Symmetry: 0=none; 2=GGA; 3=hybrids)                            \n");
		}
		if (!strcmp(argv[i], "sar")) {
			print_global(fp);
			fprintf(fp, "#Electronic Relaxation															  \n");
			fprintf(fp, " NELM   =  60           (Max electronic SCF steps)                               \n");
			fprintf(fp, " NELMIN =  6            (Min electronic SCF steps)                               \n");
			fprintf(fp, " EDIFF  =  1E-06        (SCF energy convergence; in eV)                          \n");
			fprintf(fp, " # GGA  =  PS           (PBEsol exchange-correlation)                            \n");
			fprintf(fp, "                                                                                 \n");
			fprintf(fp, "#Ionic Relaxation                                                                \n");
			fprintf(fp, " NSW    =  100          (Max ionic steps)                                        \n");
			fprintf(fp, " IBRION =  2            (Algorithm: 0-MD; 1-Quasi-New; 2-CG)                     \n");
			fprintf(fp, " ISIF   =  2            (Stress/relaxation: 2-Ions, 3-Shape/Ions/V, 4-Shape/Ions)\n");
			fprintf(fp, " EDIFFG = -2E-02        (Ionic convergence; eV/AA)                               \n");
			fprintf(fp, " # ISYM =  2            (Symmetry: 0=none; 2=GGA; 3=hybrids)                     \n");
		}
		if (!strcmp(argv[i], "mpc")) {
			fprintf(fp, "#Magnetic Proproty Calculation												 \n");
			fprintf(fp, "ISPIN      =  2        (Spin polarised DFT)                                     \n");
			fprintf(fp, "! MAGMOM   =           (Set this parameters manually)                           \n");
			fprintf(fp, "LASPH      = .TRUE.    (Non-spherical elements; d/f convergence)                \n");
			fprintf(fp, "GGA_COMPAT = .FALSE.   (Apply spherical cutoff on gradient field)               \n");
			fprintf(fp, "VOSKOWN    =  1        (Enhances the magnetic moments and the magnetic energies)\n");
			fprintf(fp, "LMAXMIX    =  4        (For d elements increase LMAXMIX to 4, f: LMAXMIX = 6)   \n");
			fprintf(fp, "! AMIX       =  0.2    (Mixing parameter to control SCF convergence)            \n");
			fprintf(fp, "# BMIX       =  0.0001 (Mixing parameter to control SCF convergence)            \n");
			fprintf(fp, "# AMIX_MAG   =  0.4    (Mixing parameter to control SCF convergence)            \n");
			fprintf(fp, "# BMIX_MAG   =  0.0001 (Mixing parameter to control SCF convergence)            \n");
		}
		if (!strcmp(argv[i], "soc")) {
			fprintf(fp, "#Spin-Orbit Coupling Calculation\n");
			fprintf(fp, "LSORBIT    = .TRUE.    (Activate SOC)\n");
			fprintf(fp, "GGA_COMPAT = .FALSE.   (Apply spherical cutoff on gradient field)\n");
			fprintf(fp, "VOSKOWN    =  1        (Enhances the magnetic moments and the magnetic energies)\n");
			fprintf(fp, "LMAXMIX    =  4        (For d elements increase LMAXMIX to 4, f: LMAXMIX = 6)\n");
			fprintf(fp, "! SAXIS    =  0 0 1    (Direction of the magnetic field)\n");
			fprintf(fp, "! MAGMOM   =  0 0 3    (Set this parameters manually, Local magnetic moment parallel to SAXIS, 3*NIONS*1.0 for non-collinear magnetic systems)\n");
			fprintf(fp, "! NBANDS   =           (Set this parameters manually, 2 * number of bands of collinear-run) \n");
		}
		if (!strcmp(argv[i], "hse")) {
			fprintf(fp, "#Hybrid function HSE06					                                                              \n");
			fprintf(fp, "LHFCALC= .TRUE.       (Activate HF)                                                                  \n");
			fprintf(fp, "AEXX   =  0.25        (25% HF exact exchange, adjusted this value to reproduce experimental band gap)\n");
			fprintf(fp, "HFSCREEN= 0.2         (Switch to screened exchange; e.g. HSE06)                                      \n");
			fprintf(fp, "ALGO   =  ALL         (Electronic Minimisation Algorithm; ALGO=58)                                   \n");
			fprintf(fp, "TIME   =  0.4         (Timestep for IALGO5X)                                                         \n");
			fprintf(fp, "PRECFOCK= N           (HF FFT grid)                                                                  \n");
			fprintf(fp, "! NKRED    = 2        (Reduce k-grid-even only, see also NKREDX, NKREDY and NKREDZ)                  \n");
			fprintf(fp, "# HFLMAX   = 4        (HF cut-off: 4d, 6f)                                                           \n");
			fprintf(fp, "# LDIAG    = .TRUE.   (Diagnolise Eigenvalues)                                                       \n");
		}
		if (!strcmp(argv[i], "vdw")) {
			fprintf(fp, "#DFT-D3 Correction			                                   \n");
			fprintf(fp, " IVDW   =  11           (DFT-D3 method of method with no damping)\n");
		}
		if (!strcmp(argv[i], "sic")) {
			fprintf(fp, "#DFT+U Calculation				                    \n");
			fprintf(fp, " LDAU   = .TRUE.        (Activate DFT+U)           \n");
			fprintf(fp, " LDATYPE=  2            (Dudarev; only U-J matters)\n");
			fprintf(fp, " LDAUL  =  2 -1         (Orbitals for each species)\n");
			fprintf(fp, " LDAUU  =  2  0         (U for each species)		\n");
			fprintf(fp, " LDAUJ  =  0  0         (J for each species)		\n");
			fprintf(fp, " LMAXMIX=  4            (Mixing cut-off; 4-d, 6-f) \n");
		}
		if (!strcmp(argv[i], "ecc")) {
			print_global(fp);
			fprintf(fp, "#Elastic constants Calculation					                                 \n");
			fprintf(fp, "# NGZ    = 500         (FFT grid mesh density for nice charge/potential plots)  \n");
			fprintf(fp, "IBRION    =  6         (Determine the Hessian matrix)                           \n");
			fprintf(fp, "NFREE     =  4         (How many displacements are used for each direction; 2-4)\n");
			fprintf(fp, "ISIF      =  3         (Stress/relaxation: 3-Shape/Ions/V)                      \n");
			fprintf(fp, "NSW       =  1         (Max ionic steps)                                        \n");
			fprintf(fp, "EDIFFG = -2E-2         (Ionic convergence; eV/AA)                          \n");
			fprintf(fp, "EDIFF  =  1E-06                                                           \n");
			//fprintf(fp, "PREC      =  Accurate      (High level)                                             \n");
		}
		if (!strcmp(argv[i], "bca")) {
			fprintf(fp, "#Bader Charge Analysis										\n");
			fprintf(fp, "LAECHG     = .TRUE.    (Write core charge into CHGCAR file)\n");
			fprintf(fp, "LCHARG     = .TRUE.    (Write CHGCAR file)                 \n");
		}
		if (!strcmp(argv[i], "elf")) {
			fprintf(fp, "#Electron Localization Function				              \n");
			fprintf(fp, "ISTART =  1            (Read existing wavefunction; if there)\n");
			fprintf(fp, "LELF   = .TRUE.        (Activate ELF)                        \n");
		}
		if (!strcmp(argv[i], "fpm")) {
			fprintf(fp, "#Frozen Phonon Method											  \n");
			fprintf(fp, "ISMEAR =  0            (Gaussian smearing)                       \n");
			fprintf(fp, "SIGMA  =  0.01         (Smearing value in eV)                          \n");
			fprintf(fp, "IBRION =  -1           (Ions are not moved)                            \n");
			fprintf(fp, "EDIFF  =  1E-06        (SCF energy convergence; in eV)                 \n");
			fprintf(fp, "PREC   =  Accurate     (Precision level)                               \n");
			fprintf(fp, "! ENCUT  =  520        (Cut-off energy for plane wave basis set, in eV)\n");
			fprintf(fp, "IALGO  =  38           (Davidson block iteration scheme)               \n");
			fprintf(fp, "LREAL  = .FALSE.       (Projection operators: false)                   \n");
			fprintf(fp, "LWAVE  = .FLASE.       (Write WAVECAR or not)                          \n");
			fprintf(fp, "LCHARG = .FLASE.       (Write CHGCAR or not)                           \n");
			fprintf(fp, "ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)          \n");
		}
		if (!strcmp(argv[i], "dfp")) {
			fprintf(fp, "#Density functional Perturbation										\n");
			fprintf(fp, "ISMEAR =  0            (Gaussian smearing)					            \n");
			fprintf(fp, "SIGMA  =  0.05         (Smearing value in eV)                          \n");
			fprintf(fp, "IBRION =  8            (determines the Hessian matrix using DFPT)      \n");
			fprintf(fp, "EDIFF  =  1E-06        (SCF energy convergence; in eV)                 \n");
			fprintf(fp, "PREC   =  Accurate     (Precision level)                               \n");
			fprintf(fp, "# ENCUT  =  520        (Cut-off energy for plane wave basis set, in eV)\n");
			fprintf(fp, "IALGO  =  38           (Davidson block iteration scheme)               \n");
			fprintf(fp, "LREAL  = .FALSE.       (Projection operators: false)                   \n");
			fprintf(fp, "LWAVE  = .FLASE.       (Write WAVECAR or not)                          \n");
			fprintf(fp, "LCHARG = .FLASE.       (Write CHGCAR or not)                           \n");
			fprintf(fp, "ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)          \n");
		}
		if (!strcmp(argv[i], "neb")) {
			print_global(fp);
			fprintf(fp, "#Budged Elastic Band (NEB)\n");
			fprintf(fp, "IMAGES =  5            \n");
			fprintf(fp, "NSW    =  500          (number of ionic steps)  \n");
			fprintf(fp, "ISMEAR =  0            (gaussian smearing method)\n");
			fprintf(fp, "SIGMA  =  0.05         (please check the width of the smearing)\n");
			fprintf(fp, "IBRION =  3            (do MD with a zero time step)\n");
			fprintf(fp, "POTIM  =  0            (Zero time step so that VASP does not move the ions)\n");
			fprintf(fp, "SPRING =  -5.0         (spring force (eV/A2) between images\n");
			fprintf(fp, "LCLIMB =  .TRUE.       (turn on the climbing image algorithm)\n");
			fprintf(fp, "ICHAIN =  0            (Indicates which method to run. NEB (ICHAIN=0) is the default)\n");
			fprintf(fp, "IOPT   =  1            (LBFGS = Limited-memory Broyden-Fletcher-Goldfarb-Shanno)\n");
		}
		if (!strcmp(argv[i], "tsd")) {
			print_global(fp);
			fprintf(fp, "#The dimer																						\n");
			fprintf(fp, "NSW    =  500          (number of ionic steps)                                                 \n");
			fprintf(fp, "IBRION =  3            (do MD with a zero time step)                                           \n");
			fprintf(fp, "POTIM  =  0            (Zero time step so that VASP does not move the ions)                    \n");
			fprintf(fp, "ICHAIN =  2            (Use the dimer method required for the latest code)                     \n");
			fprintf(fp, "DdR    =  0.005        (The dimer separation (twice the distance between images)               \n");
			fprintf(fp, "DRotMax  =  1          (Maximum number of rotation steps per translation step)                 \n");
			fprintf(fp, "DFNMin =  0.01         (Magnitude of the rotational force below which the dimer is not rotated)\n");
			fprintf(fp, "DFNMax =  1.0          (Magnitude of the rotational force below which dimer rotation stops)    \n");
			fprintf(fp, "IOPT   =  2            (CG = Conjugate Gradient)                                               \n");
		}
		if (!strcmp(argv[i], "pbs")) {
			fprintf(fp, "Global Parameters\n");
			fprintf(fp, " ISTART =  1            (Read existing wavefunction; if there)\n");
			fprintf(fp, " # ISPIN =  1            (Spin polarised DFT)\n");
			fprintf(fp, " ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)\n");
			fprintf(fp, " LREAL  = .FALSE.       (Projection operators: automatic)\n");
			fprintf(fp, " ENCUT  =  520        (Cut-off energy for plane wave basis set, in eV)\n");
			fprintf(fp, " PREC   =  Accurate       (Precision level)\n");
			fprintf(fp, " #LWAVE  = .TRUE.        (Write WAVECAR or not)\n");
			fprintf(fp, " #LCHARG = .TRUE.        (Write CHGCAR or not)\n");
			fprintf(fp, " ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)\n");
			fprintf(fp, " # LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)\n");
			fprintf(fp, " # LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)\n");
			fprintf(fp, " # NELECT =             (No. of electrons: charged cells; be careful)\n");
			fprintf(fp, " # LPLANE = .TRUE.      (Real space distribution; supercells)\n");
			fprintf(fp, " #NPAR   = 4           (Max is no. nodes; don't set for hybrids)\n");
			fprintf(fp, " # NWRITE = 2           (Medium-level output)\n");
			fprintf(fp, " # KPAR   = 2           (Divides k-grid into separate groups)\n");
			fprintf(fp, " # NGX    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
			fprintf(fp, " # NGY    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
			fprintf(fp, " # NGZ    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
			fprintf(fp, "#Static State Calculation							     			\n");
			fprintf(fp, "NSW    =  0          (number of ionic steps)                             \n");
			fprintf(fp, "ISMEAR =  0            (gaussian smearing method)              \n");
			fprintf(fp, "SIGMA  =  0.1         (please check the width of the smearing)\n");
			fprintf(fp, "IBRION =  -1            (Algorithm: 0-MD; 1-Quasi-New; 2-CG)                \n");
			fprintf(fp, "LORBIT =  11           (PAW radii for projected DOS)           \n");
			fprintf(fp, "#NEDOS  =  2001         (DOSCAR points)                         \n");
			fprintf(fp, "NELM   =  60           (Max electronic SCF steps)              \n");
			fprintf(fp, "        (SCF energy convergence; in eV)         \n");
			fprintf(fp, "EDIFF  =  1E-06                                                           \n");
			fprintf(fp, "EDIFFG = -2E-2      (Ionic convergence; eV/AA)                          \n");
		}
		if (!strcmp(argv[i], "bbs")) {
			fprintf(fp, "Global Parameters\n");
			fprintf(fp, " ISTART =  1            (Read existing wavefunction; if there)\n");
			fprintf(fp, " # ISPIN =  1            (Spin polarised DFT)\n");
			fprintf(fp, " ICHARG =  11         (Non-self-consistent: GGA/LDA band structures)\n");
			fprintf(fp, " LREAL  = .FALSE.       (Projection operators: automatic)\n");
			fprintf(fp, " ENCUT  =  520        (Cut-off energy for plane wave basis set, in eV)\n");
			fprintf(fp, " PREC   =  Accurate       (Precision level)\n");
			fprintf(fp, " #LWAVE  = .TRUE.        (Write WAVECAR or not)\n");
			fprintf(fp, " #LCHARG = .TRUE.        (Write CHGCAR or not)\n");
			fprintf(fp, " ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)\n");
			fprintf(fp, " # LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)\n");
			fprintf(fp, " # LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)\n");
			fprintf(fp, " # NELECT =             (No. of electrons: charged cells; be careful)\n");
			fprintf(fp, " # LPLANE = .TRUE.      (Real space distribution; supercells)\n");
			fprintf(fp, " #NPAR   = 4           (Max is no. nodes; don't set for hybrids)\n");
			fprintf(fp, " # NWRITE = 2           (Medium-level output)\n");
			fprintf(fp, " # KPAR   = 2           (Divides k-grid into separate groups)\n");
			fprintf(fp, " # NGX    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
			fprintf(fp, " # NGY    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
			fprintf(fp, " # NGZ    = 500         (FFT grid mesh density for nice charge/potential plots)\n");
			fprintf(fp, "#Static State Calculation							     			\n");
			fprintf(fp, "NSW    =  0          (number of ionic steps)                             \n");
			fprintf(fp, "ISMEAR =  0            (gaussian smearing method)              \n");
			fprintf(fp, "SIGMA  =  0.1         (please check the width of the smearing)\n");
			fprintf(fp, "IBRION =  -1            (Algorithm: 0-MD; 1-Quasi-New; 2-CG)                \n");
			fprintf(fp, "LORBIT =  10           (PAW radii for projected DOS)           \n");
			fprintf(fp, "#NEDOS  =  2001         (DOSCAR points)                         \n");
			fprintf(fp, "NELM   =  60           (Max electronic SCF steps)              \n");
			fprintf(fp, "        (SCF energy convergence; in eV)         \n");
			fprintf(fp, "EDIFF  =  1E-06                                                           \n");
			fprintf(fp, "EDIFFG = -2E-2      (Ionic convergence; eV/AA)                          \n");
		}
		if (!strcmp(argv[i], "dos")) {
			fprintf(fp, "#Density of States\n");
			fprintf(fp, "IBRION = -1\n");
			fprintf(fp, "NSW = 0  \n");
			fprintf(fp, "ISMEAR = -5\n");
			fprintf(fp, "LORBIT = 11\n");
			fprintf(fp, "NEDOS = 1000\n");
		}
		if (!strcmp(argv[i], "bcd")) {
			fprintf(fp, "#Bader Charge Density                   \n");
			fprintf(fp, "LAECHG = .TRUE.    (Write core charge into CHGCAR file)\n");
			fprintf(fp, "LCHARG = .TRUE.    (Write CHGCAR file)           \n");
		}
		if (!strcmp(argv[i], "scd")) {
			fprintf(fp, "#Bader Charge Density                   \n");
			fprintf(fp, "ISPIN =  2            (Spin polarised DFT)\n");
			fprintf(fp, "LAECHG = .TRUE.    (Write core charge into CHGCAR file)\n");
			fprintf(fp, "LCHARG = .TRUE.    (Write CHGCAR file)           \n");
		}
		if (!strcmp(argv[i], "pcd")) {
			fprintf(fp, "#Bader Charge Density                    \n");
			fprintf(fp, "ISPIN = 2      (Spin polarised DFT)        \n");
			fprintf(fp, "LAECHG = .TRUE.    (Write core charge into CHGCAR file)\n");
			fprintf(fp, "LCHARG = .TRUE.    (Write CHGCAR file)            \n");
		}
		if (!strcmp(argv[i], "pcd_ik")) {
			fprintf(fp, "#Decomposed Charge Density                                               \n");
			fprintf(fp, "ISTART = 1        (Job: 0 - new  1 - cont  2 - samecut)                  \n");
			fprintf(fp, "ICHARG = 1        (Read charge : 1 - file 2 - atom 10 - const)           \n");
			fprintf(fp, "LPARD = .TRUE.        (Activate decomposed charge density)               \n");
			fprintf(fp, "LSEPB = .TRUE.        (Separately write PARCHG.nb by every band or not)  \n");
			fprintf(fp, "LSEPK = .TRUE.        (Separately write PARCHG.nk by every kpoint or not)\n");
			fprintf(fp, "                                      \n");
			fprintf(fp, "# Method I : Partial Charge for the specified BANDSand KPOINTS        \n");
			fprintf(fp, "IBAND = 20 21 22 23   (Set this parameters manually)            \n");
			fprintf(fp, "KPUSE = 1 2 3 4       (Set this parameters manually)            \n");
			fprintf(fp, "                                      \n");
			fprintf(fp, "# * **********Notes * ************                      \n");
			fprintf(fp, "# (1) Copy IBZKPT as KPOINTS for static calculation,            \n");
			fprintf(fp, "# (2) Band structure calculation.                      \n");
		}
		if (!strcmp(argv[i], "pcd_en")) {
			fprintf(fp, "#Decomposed Charge Density                                               \n");
			fprintf(fp, "ISTART = 1        (Job: 0 - new  1 - cont  2 - samecut)                  \n");
			fprintf(fp, "ICHARG = 1        (Read charge : 1 - file 2 - atom 10 - const)           \n");
			fprintf(fp, "LPARD = .TRUE.        (Activate decomposed charge density)               \n");
			fprintf(fp, "LSEPB = .TRUE.        (Separately write PARCHG.nb by every band or not)  \n");
			fprintf(fp, "LSEPK = .TRUE.        (Separately write PARCHG.nk by every kpoint or not)\n");
			fprintf(fp, "                                      \n");
			fprintf(fp, "# Method II : Partial Charge in the energy rang of[-10.3 - 5.1]      \n");
			fprintf(fp, "EINT = -10.3 - 5.1    (Set this parameters manually)            \n");
			fprintf(fp, "                                      \n");
			fprintf(fp, "# * **********Notes * ************                      \n");
			fprintf(fp, "# (1) Copy IBZKPT as KPOINTS for static calculation,            \n");
			fprintf(fp, "# (2) Band structure calculation.                      \n");
		}
		if (!strcmp(argv[i], "pcd_ef")) {
			fprintf(fp, "#Decomposed Charge Density                                               \n");
			fprintf(fp, "ISTART = 1        (Job: 0 - new  1 - cont  2 - samecut)                  \n");
			fprintf(fp, "ICHARG = 1        (Read charge : 1 - file 2 - atom 10 - const)           \n");
			fprintf(fp, "LPARD = .TRUE.        (Activate decomposed charge density)               \n");
			fprintf(fp, "LSEPB = .TRUE.        (Separately write PARCHG.nb by every band or not)  \n");
			fprintf(fp, "LSEPK = .TRUE.        (Separately write PARCHG.nk by every kpoint or not)\n");
			fprintf(fp, "                                      \n");
			fprintf(fp, "# Method III : Partial Charge in the energy rang of[EF - 1 - EF]      \n");
			fprintf(fp, "NBMOD = -3                                \n");
			fprintf(fp, "EINT = -1            (Set this parameters manually)            \n");
			fprintf(fp, "                                      \n");
			fprintf(fp, "# * **********Notes * ************                      \n");
			fprintf(fp, "# (1) Copy IBZKPT as KPOINTS for static calculation,            \n");
			fprintf(fp, "# (2) Band structure calculation.                      \n");
		}
		if (!strcmp(argv[i], "wfn")) {
			fprintf(fp, "#Bader Charge Density                   \n");
			fprintf(fp, "LAECHG = .TRUE.    (Write core charge into CHGCAR file)\n");
			fprintf(fp, "LCHARG = .TRUE.    (Write CHGCAR file)           \n");
		}
		if (!strcmp(argv[i], "tlp")) {
			fprintf(fp, "#Total Local Potential                    \n");
			fprintf(fp, "NPAR = 4 # number of compute cores that work on an individual orbital								 \n");
			fprintf(fp, "ISTART = 1 # whether or not to read the WAVECAR file.												 \n");
			fprintf(fp, "ICHARG = 1 # how VASP constructs the initial charge density.										 \n");
			fprintf(fp, "LWAVE = .TRUE. # whether the wavefunctions are written to the WAVECAR file							 \n");
			fprintf(fp, "LCHARG = .TRUE. # LCHARG determines whether the charge densities(files CHGCAR and CHG) are written.	 \n");
			fprintf(fp, "LVTOT = .TRUE. # whether the total local potential is written to the LOCPOT file					 \n");
			fprintf(fp, "LVHAR = .FALSE. # whether the electrostatic potential is written to the LOCPOT file					 \n");
			fprintf(fp, "LELF = .FALSE. # whether to create an ELFCAR file or not.											 \n");
			fprintf(fp, "ENCUT = 400 # cutoff energy for the planewave basis set in eV										 \n");
			fprintf(fp, "ISMEAR = 1 # how the partial occupancies are set for each orbital.									 \n");
			fprintf(fp, "SIGMA = 0.2 # the width of the smearing in eV.														 \n");
			fprintf(fp, "EDIFF = 1E-6 # global break condition for the electronic SC - loop									 \n");
			fprintf(fp, "NELMIN = 5 # specifies the minimum number of electronic SCF steps.									 \n");
			fprintf(fp, "NELM = 300 # maximum number of electronic SC(selfconsistency) steps									 \n");
			fprintf(fp, "GGA = PE # type of generalized - gradient - approximation one wishes to use.						 \n");
			fprintf(fp, "LREAL = Auto # whether the projection operators are evaluated in real - space or in reciprocal space.\n");
			fprintf(fp, "ICORELEVEL = 2																						 \n");
			fprintf(fp, "CLNT = 1   # species																				 \n");
			fprintf(fp, "CLN = 1  # main quantum number of excited core electron												 \n");
			fprintf(fp, "CLL = 0																								 \n");
			fprintf(fp, "CLZ = 0.5																							 \n");
		}
		if (!strcmp(argv[i], "tdm"))
		{
			fprintf(fp, "#Transition Dipole Moment                   \n");
			fprintf(fp, "ISTART = 1			 \n");
			fprintf(fp, "LREAL = F			 \n");
			fprintf(fp, "PREC = Normal		 \n");
			fprintf(fp, "LWAVE = .TRUE.		 \n");
			fprintf(fp, "LCHARG = .TRUE.		 \n");
			fprintf(fp, "ADDGRID = .TRUE.	 \n");
			fprintf(fp, "ENCUT = 400			 \n");
			fprintf(fp, "ICHARG = 11			 \n");
			fprintf(fp, "					 \n");
			fprintf(fp, "ISMEAR = 0			 \n");
			fprintf(fp, "SIGMA = 0.05		 \n");
			fprintf(fp, "NELM = 60			 \n");
			fprintf(fp, "NELMIN = 6			 \n");
			fprintf(fp, "EDIFF = 1E-07		 \n");
		}
		if (!strcmp(argv[i], "los"))
		{
			fprintf(fp, "#Linear Optical Spectrums                 \n");
			fprintf(fp, "PREC = Normal	  \n");
			fprintf(fp, "ENCUT = 250.0	  \n");
			fprintf(fp, "ALGO = BSE		  \n");
			fprintf(fp, "ISMEAR = 0		  \n");
			fprintf(fp, "SIGMA = 0.05	  \n");
			fprintf(fp, "NBANDS = 216	  \n");
			fprintf(fp, "NOMEGA = 72		  \n");
			fprintf(fp, "LSPECTRAL = .TRUE.\n");
			fprintf(fp, "LOPTICS = .TRUE.  \n");
			fprintf(fp, "OMEGAMAX = 10	  \n");
			fprintf(fp, "NBANDSO = 16	  \n");
			fprintf(fp, "NBANDSV = 16	  \n");
			fprintf(fp, "NEDOS = 3000	  \n");
		}
		if (!strcmp(argv[i], "esp"))
		{
			fprintf(fp, "#Electrostatic Potential                 \n");
			fprintf(fp, "LWAVE = T  \n");
			fprintf(fp, "LCHARG = T \n");
			fprintf(fp, "LVTOT = F  \n");
			fprintf(fp, "LVHAR = T  \n");
			fprintf(fp, "LELF = F   \n");
			fprintf(fp, "LORBIT = 11\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return filename;
}

void INCAR_fix(char name[], vector<string> value)
{
	string fix = "";
	for (int i = 0; i < value.size(); i++)
		fix += i == value.size() - 1 ? value[i] : value[i] + " ";
	FILE* fp1 = fopen("INCAR", "r");
	FILE* fp2 = fopen("Buf", "w");
	char newline[50];
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
		if (strstr(buf1, name) == NULL)
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
	FILE* fp4 = fopen("Buf", "r");
	while (fgets(buf2, 1024, fp4) != NULL)
	{
		fputs(buf2, fp3);
	}
	fclose(fp3);
	fclose(fp4);
	remove("Buf");
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
		if (pos == argc || (flag && !strcmp(argv[pos], "-")))
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
	if (argc < 2)
		return;
	vector<string> keyword;
	for (int i = 2; i < argc; i++)
		keyword.push_back(argv[i]);
	FILE* fp1 = fopen("INCAR", "r");
	if (fp1 == NULL)
		return;
	FILE* fp2 = fopen("Buf", "w");
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
	FILE* fp4 = fopen("Buf", "r");
	while (fgets(buf2, 1024, fp4) != NULL)
	{
		fputs(buf2, fp3);
	}
	fclose(fp3);
	fclose(fp4);
	remove("Buf");
}

void INCAR_remove(int argc, char* argv[])
{
	if (argc < 2)
		return;
	vector<string> keyword;
	for (int i = 2; i < argc; i++)
		keyword.push_back(argv[i]);
	FILE* fp1 = fopen("INCAR", "r");
	if (fp1 == NULL)
		return;
	FILE* fp2 = fopen("Buf", "w");
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
			sprintf(newline, "#%s", buf1);
			fputs(newline, fp2);
		}
		else
			fputs(buf1, fp2);
	}
	fclose(fp1);
	fclose(fp2);
	FILE* fp3 = fopen("INCAR", "w");
	FILE* fp4 = fopen("Buf", "r");
	while (fgets(buf2, 1024, fp4) != NULL)
	{
		fputs(buf2, fp3);
	}
	fclose(fp3);
	fclose(fp4);
	remove("Buf");
}

void INCAR_append(int argc, char* argv[])
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
		if (pos == argc || (flag && !strcmp(argv[pos], "^")))
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
	if (fp == NULL)
		return;
	for (int i = 0; i < keyword.size(); i++)
	{
		fprintf(fp, "%s = ", keyword[i]);
		for (int j = 0; j < value[i].size(); j++)
		{
			fprintf(fp, "%s ", value[i][j].c_str());
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void pcd_model(int argc, char* argv[])
{
	if (strcmp("--i", argv[1]) || strcmp("-pcd", argv[2]))
		return;
	if (!strcmp("-ik", argv[3]))
	{
		if (access("INCAR", 0) == -1 || GetInfoINCAR("KPUSE").size() == 0)
		{
			const char* list[] = { NULL ,NULL,"pcd_ik" };
			string filename = write_INCAR(3, const_cast<char**> (list));
			copy_file("INCAR", filename.c_str());
		}
		vector<string> ik;
		int position = 0;
		for (int i = 4; i < argc; i++)
		{
			if (!strcmp("-ib", argv[i]))
			{
				position = i;
				break;
			}
			ik.push_back(argv[i]);
		}
		INCAR_fix(const_cast<char*>("KPUSE"), ik);
		if (position)
		{
			vector<string> ib;
			for (int i = position + 1; i < argc; i++)
				ib.push_back(argv[i]);
			INCAR_fix(const_cast<char*>("IBAND"), ib);
		}
	}
	if (!strcmp("-ib", argv[3]))
	{
		if (access("INCAR", 0) == -1 || GetInfoINCAR("IBAND").size() == 0)
		{
			const char* list[] = { NULL ,NULL,"pcd_ik" };
			string filename = write_INCAR(3, const_cast<char**> (list));
			copy_file("INCAR", filename.c_str());
		}
		vector<string> ib;
		int position = 0;
		for (int i = 4; i < argc; i++)
		{
			if (!strcmp("-ik", argv[i]))
			{
				position = i;
				break;
			}
			ib.push_back(argv[i]);
		}
		INCAR_fix(const_cast < char*>("IBAND"), ib);
		if (position)
		{
			vector<string> ik;
			for (int i = position + 1; i < argc; i++)
				ik.push_back(argv[i]);
			INCAR_fix(const_cast < char*>("KPUSE"), ik);
		}
	}
	if (!strcmp("-en", argv[3]))
	{
		const char* list[] = { NULL ,NULL,"pcd_en" };
		string filename = write_INCAR(3, const_cast<char**> (list));
		copy_file("INCAR", filename.c_str());
		if (GetInfoINCAR("NBMOD").size() != 0)
			printf("The current mode used may not be method II,please check!\n");
		else if (argc < 6)
			printf("You should set a Energy range(two values)\n");
		vector<string> en;
		for (int i = 4; i < argc; i++)
			en.push_back(argv[i]);
		INCAR_fix(const_cast < char*>("EINT"), en);
	}
	if (!strcmp("-ef", argv[3]))
	{
		const char* list[] = { NULL ,NULL,"pcd_ef" };
		string filename = write_INCAR(3, const_cast<char**> (list));
		copy_file("INCAR", filename.c_str());
		if (GetInfoINCAR("NBMOD") != "-3")
			printf("The current mode used may not be method III,please check!\n");
		else if (argc > 5)
			printf("You should set only one energy value\n");
		vector<string> ef;
		for (int i = 4; i < argc; i++)
			ef.push_back(argv[i]);
		INCAR_fix(const_cast<char*>("EINT"), ef);
	}
}