#include"../../include/VASPMATE_include/VASPMATE_main.h"
#include"../../include/VASPMATE_include/plotWave.h"
#include"../../include/VASPMATE_include/structure_operator.h"
#include"../../include/VASPMATE_include/potcar.h"
#include"../../include/VASPMATE_include/bskpt.h"
#include"../../include/VASPMATE_include/dos.h"
#include"../../include/VASPMATE_include/incar.h"
#include"../../include/VASPMATE_include/band.h"
#include"../../include/VASPMATE_include/chgcar.h"
#include"../../include/VASPMATE_include/acneb.h"
#include"../../include/VASPMATE_include/cif.h"
#include"../../include/VASPMATE_include/Eigenval.h"
#include"../../include/VASPMATE_include/shermo.h"
#include"../../include/VASPMATE_include/outcar.h"
#include"../../include/VASPMATE_include/potential.h"
#include"../../include/VASPMATE_include/spa_plot.h"
#include"../../include/VASPMATE_include/emc.h"
#include"../../include/VASPMATE_include/elastic.h"
#include"../../include/VASPMATE_include/elastic_DL.h"
#include"../../include/VASPMATE_include/magn.h"
#include"../../include/VASPMATE_include/tools.h"
#include"../../include/VASPMATE_include/Magcouple.h"
#include"../../include/VASPMATE_include/incarstd.h"
#include"../../include/VASPMATE_include/clean.h"
#include"../../include/VASPMATE_include/aimd.h"
#include"../../include/VASPMATE_include/enth.h"
#include"../../include/VASPMATE_include/sqlitedb.h"
#include"../../include/VASPMATE_include/RaW4db.h"
#include"../../include/VASPMATE_include/vel.h"

int VASPMATE_main(int argc, char* argv[])
{
	//start
	if (argc == 1)
	{
		printf("***---------------------------------------VASPMATE Version 2.0.0 (2024.10.28)----------------------------------------***\n");
		printf("***     An integrated user-interface program for high-throughput first principles computations through VASP code.    ***\n");
		printf("***             Copyright[c] 2022-2024, Beihang University by Zhaocheng Pan, Zhuoye Hu and Ruifeng Zhang.            ***\n");
		printf("***                              Please send bugs and suggestions to zrfcms@buaa.edu.cn                              ***\n");
		printf("***------------------------------------------------------------------------------------------------------------------***\n");
		return 0;
	}
	//version VASPMATE -v/-version
	if (!strcmp("--v", argv[1]) || !strcmp("--version", argv[1]))
	{
		printf("***----------------------------------------VASPMATE Version 2.0.0 (2024.10.28)---------------------------------------***\n");
		printf("***     An integrated user-interface program for high-throughput first principles computations through VASP code.    ***\n");
		printf("***            Copyright[c] 2022-2024, Beihang University by Zhaocheng Pan, Zhuoye Hu and Ruifeng Zhang.             ***\n");
		printf("***------------------------------------------------------------------------------------------------------------------***\n");
		return 0;
	}
	//license VASPMATE -l/-license
	if (!strcmp("--l", argv[1]) || !strcmp("--license", argv[1]))
	{
		printf("***------------------------------------------------------------------------------------------------------------------***\n");
		printf("***     An integrated user-interface program for high-throughput first principles computations through VASP code.    ***\n");
		printf("***            Copyright[c] 2022-2024, Beihang University by Zhaocheng Pan, Zhuoye Hu and Ruifeng Zhang.             ***\n");
		printf("***----------------------------------------VASPMATE Version 2.0.0 (2024.10.28)---------------------------------------***\n");
		printf("***               This program is currently copyrighted and distributed free of charge for academic,                 ***\n");
		printf("***                    scientific and educational and non-commercial users with our permission.                      ***\n");
		printf("***                 You are welcome to redistribute it under certain conditions with our permission.                 ***\n");
		printf("***                         Part of these terms may be changed without prior announcement.                           ***\n");
		printf("***                    This program is provided as is without any expressed or implied warranty.                     ***\n");
		printf("***------------------------------------------------------------------------------------------------------------------***\n");
		return 0;
	}
	//license VASPMATE -h/-help
	if (!strcmp("--h", argv[1]) || !strcmp("--help", argv[1]))
	{
		printf("************************************************************************************************************************\n");
		printf("***                                    The syntax format and rules for VASPMATE:                                     ***\n");
		printf("***               VASPMATE <--mode/--module> (inputfile) (outputfile) <-option/-parameter> [list/value]              ***\n");
		printf("************************************************************************************************************************\n");
		printf("***------------------------POSCAR-------------------------|*|--------------------------INCAR-------------------------***\n");
		printf("*** VASPMATE --prim         # Generate Primitive cell     |*| VASPMATE --i/inp -option  # Load Templates             ***\n");
		printf("*** VASPMATE --unit         # Generate Unit cell          |*|          option = rlx/stc/nsc/mds/sar/mag/opt/soc/hse/ ***\n");
		printf("*** VASPMATE --super        # Generate Super cell         |*|          vdw/sic/ecc/bca/elf/fpm/dfp/neb/tsd/dos/pbs/  ***\n");
		printf("*** VASPMATE --ieee         # Generate IEEE cell          |*|          bbs/bcd/scd/pcd/tlp/esp/wfn/tdm/los/...       ***\n");
		printf("*** VASPMATE --vol          # Get cell volume             |*| VASPMATE --ia/id/irm/irp  # keyword/parameter revision ***\n");
		printf("*** VASPMATE --cell         # Get cell info               |*| VASPMATE --is/istd # keyword/parameter standardization ***\n");
		printf("*** VASPMATE --num          # Get atom number             |*| VASPMATE --iu/ldau # LDA+U setting                     ***\n");
		printf("*** VASPMATE --atom         # Get atom info               |*| VASPMATE --im/imag # Magnetic moment setting           ***\n");
		printf("*** VASPMATE --symm         # Get symmetry                |*| VASPMATE --iv/ivdw # VDW+D setting                     ***\n");
		printf("*** VASPMATE --pos2cif/cif2pos  # Transform CIF&POS       |*|-------------------------KPOINTS------------------------***\n");
		printf("*** VASPMATE --fixc/fixa/fixe   # Fix atomic positions    |*| VASPMATE --k/km/kpt/kmesh  # Set KPOINTS by Input      ***\n");
		printf("*** VASPMATE --cartes/direct    # Format Transformation   |*| VASPMATE --ka              # Auto set KPOINTS by KPPRA ***\n");
		printf("*** VASPMATE --sortc/sorte      # Sort atom order         |*| VASPMATE --kv              # Auto set KPOINTS by KSPAC ***\n");
		printf("*** VASPMATE --movc/movd        # Move atom position      |*|--------------------------POTCAR------------------------***\n");
		printf("*** VASPMATE --proj/redef   # Cell Projection/Orientation |*| VASPMATE --pot          # Generate POTCAR via POSCAR   ***\n");
		printf("*** VASPMATE --affine/alias # Cell Deformation            |*| VASPMATE --pote         # Generate POTCAR via Element  ***\n");
		printf("************************************************************************************************************************\n");
		printf("***--------------------------------------------------Module-related--------------------------------------------------***\n");
		printf("*** VASPMATE --std3d/std2d  # Standardize POSCAR for band |*| VASPMATE --3dka/3dkv/3dkm  # 3D band kpoints           ***\n");
		printf("*** VASPMATE --ka3d/ka2d    # Create KPATH for band       |*| VASPMATE --3dbs            # 3D band derivation        ***\n");
		printf("*** VASPMATE --band         # Band-Structure related      |*| VASPMATE --db              # Database Creation         ***\n");
		printf("*** VASPMATE --dos          # DOS-related                 |*| VASPMATE --db2js           # Database Transformation   ***\n");
		printf("*** VASPMATE --bader        # Bader-charge-related        |*| VASPMATE --thermo          # Thermo correction         ***\n");
		printf("*** VASPMATE --neb          # NEB-related                 |*| VASPMATE --enth            # Formation enthalpy        ***\n");
		printf("*** VASPMATE --vcd          # Charge-density-related      |*| VASPMATE --elas/elae       # Elastic Property          ***\n");
		printf("*** VASPMATE --pcd          # Partial-charge-density      |*| VASPMATE --opti            # Optical Property          ***\n");
		printf("*** VASPMATE --wfn          # Wave-function-related       |*| VASPMATE --magn            # Magnetic Property         ***\n");
		printf("*** VASPMATE --fska/fskv/fskm  # Fermi surface kpoints    |*| VASPMATE --sto             # Stochastics Model         ***\n");
		printf("*** VASPMATE --fsxd/fs         # Fermi surface derivation |*| VASPMATE --evo             # Evolutionary Model        ***\n");
		printf("************************************************************************************************************************\n");
		return 0;
	}
	// VASPMATE --cif2pos file1(INCIF) file2(NEWPOS)
	if (!strcmp("--cif2pos", argv[1]))
	{
		if (argc < 2)
			return -1;
		else if (argc == 2)
			TranCIFToPOSCAR("INCIF", "NEWPOS");
		else if (argc == 3)
			TranCIFToPOSCAR(argv[2], "NEWPOS");
		else if (argc == 4)
			TranCIFToPOSCAR(argv[2], argv[3]);
		else
			return -1;
		return 0;
	}
	//VASPMATE --pos2cif file1(INPOS) file2(NEWCIF)
	if (!strcmp("--pos2cif", argv[1]))
	{
		if (argc < 2)
			return -1;
		else if (argc == 2)
			TranPOSCARToCIF("INPOS", "NEWCIF");
		else if (argc == 3)
			TranPOSCARToCIF(argv[2], "NEWCIF");
		else if (argc == 4)
			TranPOSCARToCIF(argv[2], argv[3]);
		else
			return -1;
		return 0;
	}
	// VASPMATE --chg2cub file1(chgcar) file2(cube)
	if (!strcmp("--chg2cub", argv[1]))
	{
		if (argc < 4)
			return -1;
		CHGCAR chg;
		chg.ChgcarToCube(argv[2], argv[3]);
		return 0;
	}
	if (!strcmp("--cub2chg", argv[1]))
	{
		if (argc < 4)
			return -1;
		Cube cube(argv[2]);
		if (argc == 4)
			cube.CubeToChgcar(argv[3], 0);
		else
			cube.CubeToChgcar(argv[3], atoi(argv[4]));
		return 0;
	}
	// VASPMATE --prim file1(INPOS) file2(PRIMPOS)
	// VASPMATE --prim file1 file2
	if (!strcmp("--prim", argv[1]))
	{
		if (argc == 2)
			EV_primitive("INPOS", "PRIMPOS");
		if (argc == 3)
			EV_primitive(argv[2], "PRIMPOS");
		if (argc == 4)
			EV_primitive(argv[2], argv[3]);
		if (argc > 4)
			return -1;
		return 0;
	}
	//VASPMATE --unit file1(INPOS) file2(UNITPOS)
	// VASPMATE --unit file1 file2
	if (!strcmp("--unit", argv[1]))
	{
		if (argc == 2)
			EV_unitcell("INPOS", "UNITPOS");
		if (argc == 3)
			EV_unitcell(argv[2], "UNITPOS");
		if (argc == 4)
			EV_unitcell(argv[2], argv[3]);
		if (argc > 4)
			return -1;
		return 0;
	}

	// VASPMATE --symm file(INPOS)
	if (!strcmp("--symm", argv[1]))
	{
		if (argc == 2)
			EV_get_symmetry("INPOS");
		if (argc == 3)
			EV_get_symmetry(argv[2]);
		if (argc > 3)
			return -1;
		return 0;
	}

	//VASPMATE --super file1(INPOS) file2(SUPEPOS) (-np) [sn1 sn2 sn3]
	if (!strcmp("--super", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-np", 3 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-np" ) && p == 2)
			{
				int super[3] = { atoi(argv[p + 1]),atoi(argv[p + 2]) ,atoi(argv[p + 3]) };
				EV_supercell("INPOS", "SUPERPOS", super);
			}
			else if (!strcmp(argv[p], "-np" ) && p == 3)
			{
				int super[3] = { atoi(argv[p + 1]),atoi(argv[p + 2]) ,atoi(argv[p + 3]) };
				EV_supercell(argv[2], "SUPERPOS", super);
			}
			else if (!strcmp(argv[p], "-np" ) && p == 4)
			{
				int super[3] = { atoi(argv[p + 1]),atoi(argv[p + 2]) ,atoi(argv[p + 3]) };
				EV_supercell(argv[2], argv[3], super);
			}
			else 
				return -1;
		}
		else if (check_par == 2)
		{
			if (argc == 7)
			{
				int super[3] = { atoi(argv[4]),atoi(argv[5]) ,atoi(argv[6]) };
				EV_supercell(argv[2], argv[3], super);
			}
			else if (argc == 6)
			{
				int super[3] = { atoi(argv[3]),atoi(argv[4]) ,atoi(argv[5]) };
				EV_supercell(argv[2], "SUPERPOS", super);
			}
			else if (argc == 5)
			{
				int super[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
				EV_supercell("INPOS", "SUPERPOS", super);
			}
			else
				return -1;
		}
		return 0;
	}

    // VASPMATE --affine file1(INPOS) file2(AFFPOS) -txx/-tyy/-tzz [strain]
    // VASPMATE --affine file1(INPOS) file2(AFFPOS) -sxy/-syz/-szx [strain]
    // VASPMATE --affine file1(INPOS) file2(AFFPOS) -pxy/-pyz/-pzx [strain]
    // VASPMATE --affine file1(INPOS) file2(AFFPOS) -tension [xx/yy/zz] [init_strain step_length step_num
    // VASPMATE --affine file1(INPOS) file2(AFFPOS) -simshear/-purshear [xy/yz/zx] [init_strain step_length step_num
	if (!strcmp("--affine", argv[1]))
	{
		int p; int p_single;
		vector<const char*> par_single = {"-txx", "-tyy", "-tzz", "-sxy", "-syz", "-szx", "-pxy", "-pyz", "-pzx"};
		vector<const char*> par_mulple = {"-tension", "-simshear", "-purshear"};
		int check_single = check_mulpara(argc, argv, par_single, 1 , &p_single);
		int check_mulple = check_mulpara(argc, argv, par_mulple, 4 , &p);
		if (check_single == 0 || check_mulple == 0)
			return -1;
		else if (check_mulple == 1 && check_single != 1)
		{
			if (p == 2)
				EV_affine("INPOS", "AFFPOS" , argv[p], argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atoi(argv[p + 4]));
			else if (p == 3)
				EV_affine(argv[2], "AFFPOS" , argv[p], argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atoi(argv[p + 4]));
			else if (p == 4)
				EV_affine(argv[2], argv[3] , argv[p], argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atoi(argv[p + 4]));
			else
				return -1;
		}
		else if (check_single == 1 && check_mulple != 1)
		{
            char model[10], mode[10];
            int matched = sscanf(argv[p_single], "-%1s%2s", model, mode);
            if (strcmp(model, "t") == 0) 
                strcpy(model, "-tension");
			else if (strcmp(model, "s") == 0) 
                strcpy(model, "-simshear");
			else if (strcmp(model, "p") == 0) 
                strcpy(model, "-purshear");
			else
                return -1;
			if (p_single == 2)
				EV_affine("INPOS", "AFFPOS" , model, mode, atof(argv[p_single + 1]), atof(argv[p_single + 1]), 0);
			else if (p_single == 3)
				EV_affine(argv[2], "AFFPOS" , model, mode, atof(argv[p_single + 1]), atof(argv[p_single + 1]), 0);
			else if (p_single == 4)
				EV_affine(argv[2], argv[3] , model, mode, atof(argv[p_single + 1]), atof(argv[p_single + 1]), 0);
			else
				return -1;
		}
		else
			return -1;
		return 0;
	}

	// VASPMATE --alias file1(INPOS) file2(ALIPOS) -tensi [xx/yy/zz] [istart iend ispac] [value]
	// VASPMATE --alias file1(INPOS) file2(ALIPOS) -tensi [xx/yy/zz] [strain position]
	// VASPMATE --alias file1(INPOS) file2(ALIPOS) -shear [xy/yz/zx] [istart1 iend1 ispac1] [istart2 iend2 ispac2] [value]
	// VASPMATE --alias file1(INPOS) file2(ALIPOS) -shear [xy/yz/zx] [strain position]
	if (!strcmp("--alias", argv[1]))
	{
		int p_ten;int p_she; int check_vec = 2; int p_vec;
		vector<const char*> par_ten = {"-txx", "-tyy", "-tzz"};
		vector<const char*> par_she = {"-sxy", "-syz", "-szx"};
		int check_vten = check_mulpara(argc, argv, par_ten, 2 , &p_ten);
		int check_vshe = check_mulpara(argc, argv, par_she, 3 , &p_she);
		int ten; int she; int check_par = 2; int p;
		int check_ten = check_para(argc, argv, "-tensi", 5 , &ten);
		int check_she = check_para(argc, argv, "-shear", 8 , &she);
		if ( check_ten == 1 || check_she == 1 )
		{
			check_par = 1;
			if (check_ten == 1)
				p = ten;
			if (check_she == 1)
				p = she;
		}
		if ( check_vten == 1 || check_vshe == 1 )
		{
			check_vec = 1;
			if (check_vten == 1)
				p_vec = p_ten;
			if (check_vshe == 1)
				p_vec = p_she;
		}
		if ((check_par == 2 && check_vec == 2) || (check_vec == 0 || check_vec == 0))
			return -1;
		else if (check_par == 1 && check_vec != 1)
		{
			if (p == 2)
			{
				if (check_ten == 1)
					EV_alias_tensile("INPOS", "ALIPOS", argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atof(argv[p + 4]), atof(argv[p + 5]));
				if (check_she == 1)
					EV_alias_shear("INPOS", "ALIPOS", argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atof(argv[p + 4]), atof(argv[p + 5]),
						atof(argv[p + 6]), atof(argv[p + 7]), atof(argv[p + 8]));
			}
			else if (p == 3)
			{
				if (check_ten == 1)
					EV_alias_tensile(argv[2], "ALIPOS", argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atof(argv[p + 4]), atof(argv[p + 5]));
				if (check_she == 1)
					EV_alias_shear(argv[2], "ALIPOS", argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atof(argv[p + 4]), atof(argv[p + 5]),
						atof(argv[p + 6]), atof(argv[p + 7]), atof(argv[p + 8]));
			}
			else if (p == 4)
			{
				if (check_ten == 1)
					EV_alias_tensile(argv[2], argv[3], argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atof(argv[p + 4]), atof(argv[p + 5]));
				if (check_she == 1)
					EV_alias_shear(argv[2], argv[3], argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), atof(argv[p + 4]), atof(argv[p + 5]),
						atof(argv[p + 6]), atof(argv[p + 7]), atof(argv[p + 8]));
			}
			else
				return -1;
		}
		else if (check_par != 1 && check_vec == 1)
		{
			char model[10], mode[10];
			int matched = sscanf(argv[p_vec], "-%1s%2s", model, mode);			
			if (p_vec == 2)
			{
				if (check_vten == 1)
					EV_alias_tensile("INPOS", "ALIPOS", mode, atof(argv[p_vec + 1]), atof(argv[p_vec + 1]), 1, atof(argv[p_vec + 2]));
				if (check_vshe == 1)
					EV_alias_shear("INPOS", "ALIPOS", mode, atof(argv[p_vec + 1]), atof(argv[p_vec + 1]), 1, atof(argv[p_vec + 2]),
						atof(argv[p_vec + 2]), 0, atof(argv[p_vec + 3]));
			}
			else if (p_vec == 3)
			{
				if (check_vten == 1)
					EV_alias_tensile(argv[2], "ALIPOS", mode, atof(argv[p_vec + 1]), atof(argv[p_vec + 1]), 1, atof(argv[p_vec + 2]));
				if (check_vshe == 1)
					EV_alias_shear(argv[2], "ALIPOS", mode, atof(argv[p_vec + 1]), atof(argv[p_vec + 1]), 1, atof(argv[p_vec + 2]),
						atof(argv[p_vec + 2]), 0, atof(argv[p_vec + 3]));
			}
			else if (p_vec == 4)
			{
				if (check_vten == 1)
					EV_alias_tensile(argv[2], argv[3], mode, atof(argv[p_vec + 1]), atof(argv[p_vec + 1]), 1, atof(argv[p_vec + 2]));
				if (check_vshe == 1)
					EV_alias_shear(argv[2], argv[3], mode, atof(argv[p_vec + 1]), atof(argv[p_vec + 1]), 1, atof(argv[p_vec + 2]),
						atof(argv[p_vec + 2]), 0, atof(argv[p_vec + 3]));
			}
			else
				return -1;
		}
		return 0;
	}

	// VASPMATE --proj file1(INPOS) file2(PROJPOS) -rot [rotx roty rotz]
	// VASPMATE --proj file1(INPOS) file2(PROJPOS) -ind [pvh pvk pvl] [uvu uvv uvw]
	// VASPMATE --proj file1(INPOS) file2(PROJPOS) -mat [mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33]
	if (!strcmp("--proj", argv[1]))
	{
		int rot; int ind; int mat; int check_par = 0; int p;
		int check_rot = check_para(argc, argv, "-rot", 3 , &rot);
		int check_ind = check_para(argc, argv, "-ind", 6 , &ind);
		int check_mat = check_para(argc, argv, "-mat", 9 , &mat);
		if ( check_rot == 1 || check_ind == 1 || check_mat == 1 )
		{
			check_par = 1;
			if (check_rot == 1)
				p = rot;
			if (check_ind == 1)
				p = ind;
			if (check_mat == 1)
				p = mat;
		}
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (check_rot == 1)
			{
				double rot[3] = { atof(argv[p + 1]),atof(argv[p + 2]) ,atof(argv[p + 3]) };
				if (p == 2)
					EV_rotproj("INPOS", "PROJPOS", rot);
				else if (p == 3)
					EV_rotproj(argv[2], "PROJPOS", rot);
				else if (p == 4)
					EV_rotproj(argv[2], argv[3], rot);
				else
					return -1;
			}
			if (check_ind == 1)
			{
				if (p == 2)
					EV_indproj("INPOS", "PROJPOS", atoi(argv[p + 1]), atoi(argv[p + 2]), atoi(argv[p + 3]), atoi(argv[p + 4]), atoi(argv[p + 5]), atoi(argv[p + 6]));
				else if (p == 3)
					EV_indproj(argv[2], "PROJPOS", atoi(argv[p + 1]), atoi(argv[p + 2]), atoi(argv[p + 3]), atoi(argv[p + 4]), atoi(argv[p + 5]), atoi(argv[p + 6]));
				else if (p == 4)
					EV_indproj(argv[2], argv[3], atoi(argv[p + 1]), atoi(argv[p + 2]), atoi(argv[p + 3]), atoi(argv[p + 4]), atoi(argv[p + 5]), atoi(argv[p + 6]));
				else
					return -1;
			}
			if (check_mat == 1)
			{
				double projmat[3][3] = { atof(argv[p + 1]),atof(argv[p + 2]) ,atof(argv[p + 3]),
					atof(argv[p + 4]),atof(argv[p + 5]) ,atof(argv[p + 6]), atof(argv[p + 7]),atof(argv[p + 8]) ,atof(argv[p + 9]) };
				if (p == 2)
					EV_matproj("INPOS", "PROJPOS", projmat);
				else if (p == 3)
					EV_matproj(argv[2], "PROJPOS", projmat);
				else if (p == 4)
					EV_matproj(argv[2], argv[3], projmat);
				else
					return -1;
			}
		}
		return 0;
	}

	//??
	if (!strcmp(argv[1], "--ads"))
	{
		if (argc == 6)
			adsorbent(argv[2], argv[3], argv[4], atof(argv[5]));
		return 0;
	}
	// VASPMATE --redef file1(INPOS) file2(REDEPOS) (-par) [vect11 vect12 vect13 vect21 vect22 vect23 vect31 vect32 vect33] [a1 x1 a2 x2]
	if (!strcmp("--redef", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 13 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
			{
				int rot[3][3] = { atoi(argv[p + 1]), atoi(argv[p + 2]), atoi(argv[p + 3]),atoi(argv[p + 4]), atoi(argv[p + 5]), atoi(argv[p + 6]),
					atoi(argv[p + 7]), atoi(argv[p + 8]), atoi(argv[p + 9]) };
				int a1, x1, a2, x2;
				if (!strcmp("a", argv[p + 10])) a1 = 0;
				else if (!strcmp("b", argv[p + 10])) a1 = 1;
				else if (!strcmp("c", argv[p + 10])) a1 = 2;

				if (!strcmp("x", argv[p + 11])) x1 = 0;
				else if (!strcmp("y", argv[p + 11])) x1 = 1;
				else if (!strcmp("z", argv[p + 11])) x1 = 2;

				if (!strcmp("a", argv[p + 12])) a2 = 0;
				else if (!strcmp("b", argv[p + 12])) a2 = 1;
				else if (!strcmp("c", argv[p + 12])) a2 = 2;

				if (!strcmp("xy", argv[p + 13])) x2 = 1 - x1;
				else if (!strcmp("xz", argv[p + 13])) x2 = 2 - x1;
				else if (!strcmp("yz", argv[p + 13])) x2 = 3 - x1;
				EV_redefine("INPOS", "REDEPOS", rot, 0.001, a1, x1, a2, x2, 1);
			}
			else if (!strcmp(argv[p], "-par" ) && p == 3)
			{
				int rot[3][3] = { atoi(argv[p + 1]), atoi(argv[p + 2]), atoi(argv[p + 3]),atoi(argv[p + 4]), atoi(argv[p + 5]), atoi(argv[p + 6]),
					atoi(argv[p + 7]), atoi(argv[p + 8]), atoi(argv[p + 9]) };
				int a1, x1, a2, x2;
				if (!strcmp("a", argv[p + 10])) a1 = 0;
				else if (!strcmp("b", argv[p + 10])) a1 = 1;
				else if (!strcmp("c", argv[p + 10])) a1 = 2;

				if (!strcmp("x", argv[p + 11])) x1 = 0;
				else if (!strcmp("y", argv[p + 11])) x1 = 1;
				else if (!strcmp("z", argv[p + 11])) x1 = 2;

				if (!strcmp("a", argv[p + 12])) a2 = 0;
				else if (!strcmp("b", argv[p + 12])) a2 = 1;
				else if (!strcmp("c", argv[p + 12])) a2 = 2;

				if (!strcmp("xy", argv[p + 13])) x2 = 1 - x1;
				else if (!strcmp("xz", argv[p + 13])) x2 = 2 - x1;
				else if (!strcmp("yz", argv[p + 13])) x2 = 3 - x1;
				EV_redefine(argv[2], "REDEPOS", rot, 0.001, a1, x1, a2, x2, 1);
			}
			else if (!strcmp(argv[p], "-par" ) && p == 4)
			{
				int rot[3][3] = { atoi(argv[p + 1]), atoi(argv[p + 2]), atoi(argv[p + 3]),atoi(argv[p + 4]), atoi(argv[p + 5]), atoi(argv[p + 6]),
					atoi(argv[p + 7]), atoi(argv[p + 8]), atoi(argv[p + 9]) };
				int a1, x1, a2, x2;
				if (!strcmp("a", argv[p + 10])) a1 = 0;
				else if (!strcmp("b", argv[p + 10])) a1 = 1;
				else if (!strcmp("c", argv[p + 10])) a1 = 2;

				if (!strcmp("x", argv[p + 11])) x1 = 0;
				else if (!strcmp("y", argv[p + 11])) x1 = 1;
				else if (!strcmp("z", argv[p + 11])) x1 = 2;

				if (!strcmp("a", argv[p + 12])) a2 = 0;
				else if (!strcmp("b", argv[p + 12])) a2 = 1;
				else if (!strcmp("c", argv[p + 12])) a2 = 2;

				if (!strcmp("xy", argv[p + 13])) x2 = 1 - x1;
				else if (!strcmp("xz", argv[p + 13])) x2 = 2 - x1;
				else if (!strcmp("yz", argv[p + 13])) x2 = 3 - x1;
				EV_redefine(argv[2], argv[3], rot, 0.001, a1, x1, a2, x2, 1);
			}
			else
				return -1;
		}
		else if (check_par == 2)
		{
			if (argc == 17)
			{
				int rot[3][3] = { atoi(argv[4]),atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
					atoi(argv[8]), atoi(argv[9]), atoi(argv[10]),atoi(argv[11]), atoi(argv[12]) };
				int a1, x1, a2, x2;
				if (!strcmp("a", argv[13])) a1 = 0;
				else if (!strcmp("b", argv[13])) a1 = 1;
				else if (!strcmp("c", argv[13])) a1 = 2;

				if (!strcmp("x", argv[14])) x1 = 0;
				else if (!strcmp("y", argv[14])) x1 = 1;
				else if (!strcmp("z", argv[14])) x1 = 2;

				if (!strcmp("a", argv[15])) a2 = 0;
				else if (!strcmp("b", argv[15])) a2 = 1;
				else if (!strcmp("c", argv[15])) a2 = 2;

				if (!strcmp("xy", argv[16])) x2 = 1 - x1;
				else if (!strcmp("xz", argv[16])) x2 = 2 - x1;
				else if (!strcmp("yz", argv[16])) x2 = 3 - x1;
				EV_redefine(argv[2], argv[3], rot, 0.001, a1, x1, a2, x2, 1);
			}
			else if (argc == 13)
			{
				int rot[3][3] = { atoi(argv[4]),atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
					atoi(argv[8]), atoi(argv[9]), atoi(argv[10]),atoi(argv[11]), atoi(argv[12]) };
				EV_redefine(argv[2], argv[3], rot, 0.001, 0, 0, 0, 0, 0);
			}
			else if (argc == 11)
			{
				int rot[3][3] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
					atoi(argv[8]), atoi(argv[9]), atoi(argv[10]) };
				EV_redefine("INPOS", "REDEPOS", rot, 0.001, 0, 0, 0, 0, 0);
			}
			else if (argc == 15)
			{
				int rot[3][3] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
					atoi(argv[8]), atoi(argv[9]), atoi(argv[10]) };
				int a1, x1, a2, x2;
				if (!strcmp("a", argv[11])) a1 = 0;
				else if (!strcmp("b", argv[11])) a1 = 1;
				else if (!strcmp("c", argv[11])) a1 = 2;

				if (!strcmp("x", argv[12])) x1 = 0;
				else if (!strcmp("y", argv[12])) x1 = 1;
				else if (!strcmp("z", argv[12])) x1 = 2;

				if (!strcmp("a", argv[13])) a2 = 0;
				else if (!strcmp("b", argv[13])) a2 = 1;
				else if (!strcmp("c", argv[13])) a2 = 2;

				if (!strcmp("xy", argv[14])) x2 = 1 - x1;
				else if (!strcmp("xz", argv[14])) x2 = 2 - x1;
				else if (!strcmp("yz", argv[14])) x2 = 3 - x1;
				EV_redefine("INPOS", "REDEPOS", rot, 0.001, a1, x1, a2, x2, 1);
			}
			else
				return -1;
		}
		return 0;
	}
	// VASPMATE --ieee file1(INPOS) file2(IEEEPOS)
	if (!strcmp("--ieee", argv[1]))
	{
		if (argc == 4)
			EV_recell(argv[2], argv[3]);
		else if (argc == 3)
			EV_recell(argv[2], "IEEEPOS");
		else if (argc == 2)
			EV_recell("INPOS", "IEEEPOS");
		else
			return -1;
		return 0;
	}
	// VASPMATE --fixc file2(INPOS) file2(FIXPOS) (-par) axis m n F F F (fix m-n)
	if (!strcmp("--fixc", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 6 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
			{
				char fix[3] = { argv[p + 4][0],argv[p + 5][0],argv[p + 6][0] };
				EV_fixatomcoor("INPOS", "FIXPOS", argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), fix);
			}
			else if (!strcmp(argv[p], "-par" ) && p == 3)
			{
				char fix[3] = { argv[p + 4][0],argv[p + 5][0],argv[p + 6][0] };
				EV_fixatomcoor(argv[2], "FIXPOS", argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), fix);
			}
			else if (!strcmp(argv[p], "-par" ) && p == 4)
			{
				char fix[3] = { argv[p + 4][0],argv[p + 5][0],argv[p + 6][0] };
				EV_fixatomcoor(argv[2], argv[3], argv[p + 1], atof(argv[p + 2]), atof(argv[p + 3]), fix);
			}
			else
				return -1;
		}
		else if (check_par == 2)
		{
			// VASPMATE --fixc axis m n F F F
			if (argc == 8)
			{
				char fix[3] = { argv[5][0],argv[6][0],argv[7][0] };
				EV_fixatomcoor("INPOS", "FIXPOS", argv[2], atof(argv[3]), atof(argv[4]), fix);
			}
			//VASPMATE --fixc file1 file2 axis m n F F F
			else if (argc == 10)
			{
				char fix[3] = { argv[7][0],argv[8][0],argv[9][0] };
				EV_fixatomcoor(argv[2], argv[3], argv[4], atof(argv[5]), atof(argv[6]), fix);
			}
			else
				return -1;
		}
		return 0;
	}
	// VASPMATE --fixa file2 file2 1 2 3 ... F F F 
	if (!strcmp("--fixa", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 4 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			char** argv_ = (char**)malloc((argc-1)*sizeof(char*));
			for (int i = 0, j = 0; i < argc; i++) 
				if (i != p)
					argv_[j++] = argv[i];
			if (!strcmp(argv[p], "-par" ) && p == 2)
			{
				char fix[3] = { argv[p + 2][0],argv[p + 3][0],argv[p + 4][0] };
				EV_fixatomindex("INPOS", "FIXPOS", argc - 1, argv_, fix);
			}
			else if (!strcmp(argv[p], "-par" ) && p == 3)
			{
				char fix[3] = { argv[p + 2][0],argv[p + 3][0],argv[p + 4][0] };
				EV_fixatomindex(argv[2], "FIXPOS", argc - 1, argv_, fix);
			}
			else if (!strcmp(argv[p], "-par" ) && p == 4)
			{
				char fix[3] = { argv[p + 2][0],argv[p + 3][0],argv[p + 4][0] };
				EV_fixatomindex(argv[2], argv[3], argc - 1, argv_, fix);
			}
			else
				return -1;			
		}
		else if (check_par == 2)
		{
			char fix[3] = { argv[argc - 3][0],argv[argc - 2][0],argv[argc - 1][0] };
			EV_fixatomindex(argv[2], argv[3], argc, argv, fix);
		}
		return 0;
	}
	// VASPMATE --fixe file2 file2 A B C... F F F 
	if (!strcmp("--fixe", argv[1]))
	{
		int p;
		vector<string> ele;
		int check_par = check_para(argc, argv, "-par", 4 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			int find_slash = p; int fix_num = 0;
			const int MAX_FIX_NUM = 100;
			char fix[MAX_FIX_NUM][3] = {{0}};
			while (find_slash < argc)
			{
				if (!strcmp(argv[find_slash], "/"))
				{
					ele.push_back(argv[find_slash - 4]);
					fix[fix_num][0] = argv[find_slash - 3][0];
					fix[fix_num][1] = argv[find_slash - 2][0];
					fix[fix_num][2] = argv[find_slash - 1][0];
					fix_num++;
				}
				else if (find_slash == argc - 1)
				{
					ele.push_back(argv[find_slash - 3]);
					fix[fix_num][0] = argv[find_slash - 2][0];
					fix[fix_num][1] = argv[find_slash - 1][0];
					fix[fix_num][2] = argv[find_slash - 0][0];
					fix_num++;
					break;
				}
				find_slash++;
			}
			if (!strcmp(argv[p], "-par" ) && p == 2)
				/*for (int i = p + 1; i < argc - 3; i++)
					ele.push_back(argv[i]);
				char fix[3] = { argv[argc - 3][0],argv[argc - 2][0],argv[argc - 1][0] };*/
				EV_fixatomele("INPOS", "FIXPOS", ele, fix);
			if (!strcmp(argv[p], "-par" ) && p == 3)
				EV_fixatomele(argv[2], "FIXPOS", ele, fix);
			if (!strcmp(argv[p], "-par" ) && p == 4)
				EV_fixatomele(argv[2], argv[3], ele, fix);
		}
		else if (check_par == 2)
		{
			int find_slash = 4; int fix_num = 0;
			const int MAX_FIX_NUM = 100;
			char fix[MAX_FIX_NUM][3] = {{0}};
			while (find_slash < argc)
			{
				if (!strcmp(argv[find_slash], "/"))
				{
					ele.push_back(argv[find_slash - 4]);
					fix[fix_num][0] = argv[find_slash - 3][0];
					fix[fix_num][1] = argv[find_slash - 2][0];
					fix[fix_num][2] = argv[find_slash - 1][0];
					fix_num++;
				}
				else if (find_slash == argc - 1)
				{
					ele.push_back(argv[find_slash - 3]);
					fix[fix_num][0] = argv[find_slash - 2][0];
					fix[fix_num][1] = argv[find_slash - 1][0];
					fix[fix_num][2] = argv[find_slash - 0][0];
					fix_num++;
					break;
				}
				find_slash++;
			}
			EV_fixatomele(argv[2], argv[3], ele, fix);
		}
		return 0;
	}
	// VASPMATE --ufix file
	if (!strcmp("--ufix", argv[1]))
	{
		if (argc == 2)
			EV_cleanfix("INPOS");
		else if (argc == 3)
			EV_cleanfix(argv[2]);
		else
			return -1;
		return 0;
	}
	// VASPMATE --swap file(INPOS) file2(SWAPOS) (-elem) [elem1 elem2...]
	if (!strcmp("--swap", argv[1]))
	{
		int e;
		vector<string> ele;
		int check_ele = check_para(argc, argv, "-elem", "-e", 1 , &e); //The fourth parameter is the minimum number of input parameters required.
		if (check_ele == 0)
			return -1;
		else if (check_ele == 1)
		{
			for (int i = e + 1; i < argc; i++)
				ele.push_back(argv[i]);
			if (e == 2)
				EV_swapele("INPOS", "SWAPOS", ele);
			else if (e == 3)
				EV_swapele(argv[2], "SWAPOS", ele);
			else if (e == 4)
				EV_swapele(argv[2], argv[3], ele);
			else
				return -1;
		}
		else if (check_ele == 2)
		{
			for (int i = 4; i < argc; i++)
				ele.push_back(argv[i]);
			EV_swapele(argv[2], argv[3], ele);
		}
		return 0;
	}
	// VASPMATE --cartes file1(INPOS) file2(NEWPOS)
	if (!strcmp("--cartes", argv[1]))
	{
		if (argc == 2)
			EV_direct_to_carts("INPOS", "NEWPOS");
		else if (argc == 3)
			EV_direct_to_carts(argv[2], "NEWPOS");
		else if (argc == 4)
			EV_direct_to_carts(argv[2], argv[3]);
		else
			return -1;
		return 0;
	}
	// VASPMATE --direct file1(INPOS) file2(NEWPOS)
	if (!strcmp("--direct", argv[1]))
	{
		if (argc == 2)
			EV_carts_to_direct("INPOS", "NEWPOS");
		else if (argc == 3)
			EV_carts_to_direct(argv[2], "NEWPOS");
		else if (argc == 4)
			EV_carts_to_direct(argv[2], argv[3]);
		else
			return -1;
		return 0;
	}
	// VASPMATE --movc file1(INPOS) file2(MOVEPOS) (-par) [axis] [min max] [dx dy dz]
	if (!strcmp("--movc", argv[1]))
	{
		char label[3] = { '\0' }; strcpy(label, "-c");
		/*for (int i = 0; i < argc; i++) 
		{	
			if (!strcmp("-c", argv[i]) || !strcmp("c", argv[i]))
				{strcpy(label, "-c"); break;}
			else if (!strcmp("-d", argv[i]) || !strcmp("d", argv[i]))
				{strcpy(label, "-d"); break;}
		}*/
		int p;
		int check_par = check_para(argc, argv, "-par", 6 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
			{
				/*if (label[0] == '\0')
				{
					double move[3] = { atof(argv[p + 4]),atof(argv[p + 5]),atof(argv[p + 6]) };
					EV_move("INPOS", "MOVEPOS", label, argv[p + 1], move, atof(argv[p + 2]), atof(argv[p + 3]));
				}*/
				double move[3] = { atof(argv[p + 4]),atof(argv[p + 5]),atof(argv[p + 6]) };
				EV_move("INPOS", "MOVEPOS", label, argv[p + 1], move, atof(argv[p + 2]), atof(argv[p + 3]));
			}
			else if (!strcmp(argv[p], "-par" ) && p == 3)
			{
				double move[3] = { atof(argv[p + 4]),atof(argv[p + 5]),atof(argv[p + 6]) };
				EV_move(argv[2], "MOVEPOS", label, argv[p + 1], move, atof(argv[p + 2]), atof(argv[p + 3]));
			}
			else if (!strcmp(argv[p], "-par" ) && p == 4)
			{
				double move[3] = { atof(argv[p + 4]),atof(argv[p + 5]),atof(argv[p + 6]) };
				EV_move(argv[2], argv[3], label, argv[p + 1], move, atof(argv[p + 2]), atof(argv[p + 3]));
			}
			else
				return -1;
		}
		else if (check_par == 2)
		{
			//VASPMATE --movc axis mmin mmax dx dy dz 
			if (argc == 8) {
				double move[3] = { atof(argv[5]),atof(argv[6]),atof(argv[7]) };
				EV_move("INPOS", "MOVEPOS", label, argv[2], move, atof(argv[3]), atof(argv[4]));
			}
			//VASPMATE --movc file1 file2 axis mmin mmax dx dy dz 
			else if (argc == 10) {
				double move[3] = { atof(argv[7]),atof(argv[8]),atof(argv[9]) };
				EV_move(argv[2], argv[3], label, argv[4], move, atof(argv[5]), atof(argv[6]));
			}
			else
				return -1;
		}
		return 0;
	}

	// VASPMATE --movd file1(INPOS) file2(MOVEPOS) (-par) [axis] [min max] [da db dc]
	if (!strcmp("--movd", argv[1]))
	{
		char label[3] = { '\0' }; strcpy(label, "-d");
		int p;
		int check_par = check_para(argc, argv, "-par", 6 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
			{
				double move[3] = { atof(argv[p + 4]),atof(argv[p + 5]),atof(argv[p + 6]) };
				EV_move("INPOS", "MOVEPOS", label, argv[p + 1], move, atof(argv[p + 2]), atof(argv[p + 3]));
			}
			else if (!strcmp(argv[p], "-par" ) && p == 3)
			{
				double move[3] = { atof(argv[p + 4]),atof(argv[p + 5]),atof(argv[p + 6]) };
				EV_move(argv[2], "MOVEPOS", label, argv[p + 1], move, atof(argv[p + 2]), atof(argv[p + 3]));
			}
			else if (!strcmp(argv[p], "-par" ) && p == 4)
			{
				double move[3] = { atof(argv[p + 4]),atof(argv[p + 5]),atof(argv[p + 6]) };
				EV_move(argv[2], argv[3], label, argv[p + 1], move, atof(argv[p + 2]), atof(argv[p + 3]));
			}
			else 
				return -1;
		}
		else if (check_par == 2)
		{
			//VASPMATE --movd axis mmin mmax dx dy dz 
			if (argc == 8) {
				double move[3] = { atof(argv[5]),atof(argv[6]),atof(argv[7]) };
				EV_move("INPOS", "MOVEPOS", label, argv[2], move, atof(argv[3]), atof(argv[4]));
			}
			//VASPMATE --movd file1 file2 axis mmin mmax dx dy dz 
			else if (argc == 10) {
				double move[3] = { atof(argv[7]),atof(argv[8]),atof(argv[9]) };
				EV_move(argv[2], argv[3], label, argv[4], move, atof(argv[5]), atof(argv[6]));
			}
			else 
				return -1;
		}
		return 0;
	}
	// VASPMATE --sortc file1(INPOS) file2(SORTPOS) mode
	if (!strcmp("--sortc", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 1 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
				EV_atomsort_coord("INPOS", "SORTPOS", argv[p + 1]);
			else if (!strcmp(argv[p], "-par" ) && p == 3)
				EV_atomsort_coord(argv[2], "SORTPOS", argv[p + 1]);
			else if (!strcmp(argv[p], "-par" ) && p == 4)
				EV_atomsort_coord(argv[2], argv[3], argv[p + 1]);
			else
				return -1;
		}
		else if (check_par == 2)
		{
			if (argc == 3)
				EV_atomsort_coord("INPOS", "SORTPOS", argv[2]);
			else if (argc == 5)
				EV_atomsort_coord(argv[2], argv[3], argv[4]);
			else
				return -1;
		}
		return 0;
	}
	// VASPMATE --sorte file1(INPOS) file2(SORTPOS) ele1 ele2
	if (!strcmp("--sorte", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 1 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
				EV_atomsort_element("INPOS", "SORTPOS", argv[p + 1], argv[p + 2]);
			else if (!strcmp(argv[p], "-par" ) && p == 3)
				EV_atomsort_element(argv[2], "SORTPOS", argv[p + 1], argv[p + 2]);
			else if (!strcmp(argv[p], "-par" ) && p == 4)
				EV_atomsort_element(argv[2], argv[3], argv[p + 1], argv[p + 2]);
			else
				return -1;
		}
		else if (check_par == 2)
		{
			if (argc == 4)
				EV_atomsort_element("INPOS", "SORTPOS", argv[2], argv[3]);
			else if (argc == 6)
				EV_atomsort_element(argv[2], argv[3], argv[4], argv[5]);
			else
				return -1;
		}
		return 0;
	}
	//INCAR
	if (!strcmp("--i", argv[1]) || !strcmp("--inp", argv[1]))
	{
		write_INCAR(argc, argv);
		return 0;
	}
	if (!strcmp("--pcd", argv[1]))
	{
		pcd_model(argc, argv);
		return 0;
	}
	if (!strcmp("--i_append", argv[1]) || !strcmp("--i_app", argv[1]) || !strcmp("--ia", argv[1]))
	{
		if (argc < 4)
			return -1;
		INCAR_append(argc, argv);
		return 0;
	}
	if (!strcmp("--i_delete", argv[1]) || !strcmp("--i_del", argv[1]) || !strcmp("--id", argv[1]))
	{
		if (argc < 3)
			return -1;
		INCAR_delete(argc, argv);
		return 0;
	}
	if (!strcmp("--i_remove", argv[1]) || !strcmp("--i_rem", argv[1]) || !strcmp("--irm", argv[1]))
	{
		if (argc < 3)
			return -1;
		INCAR_remove(argc, argv);
		return 0;
	}
	if (!strcmp("--i_replace", argv[1]) || !strcmp("--i_rep", argv[1]) || !strcmp("--irp", argv[1]))
	{
		if (argc < 4)
			return -1;
		INCAR_replace(argc, argv);
		return 0;
	}
	//VASPMATE --is/--istd (--i_std or --i_standard)
	if (!strcmp("--is", argv[1]) || !strcmp("--istd", argv[1]) || !strcmp("--i_std", argv[1]) || !strcmp("--i_standard", argv[1]))
	{
		inpstd::incarstd();
		return 0;
	}
	//Add default ldau to incar according to element
	if (!strcmp("--ldau", argv[1]) || !strcmp("--iu", argv[1]))
	{
		int a; int t;
		int check_a = check_para(argc, argv, "-a", 0 , &a); //1 2
		int check_t = check_para(argc, argv, "-t", 1 , &t); //0 1 2
		char table[100];
		if (check_t == 2 && check_a == 2)
		{
			if (argc == 2)
			{
				LDAU ldau(nullptr, "INPOS");
				ldau.AddLDAU_default("INPOS");
			}
			else if (argc == 3)
			{
				LDAU ldau(nullptr, argv[2]);
				ldau.AddLDAU_default(argv[2]);
			}
			else
				return -1;
		}
		else if (check_t == 1 && check_a == 2 )
		{
			strcpy(table, argv[t + 1]);
			if (a == 2 || t == 2)
			{
				LDAU ldau(table, "INPOS");
				ldau.AddLDAU_default("INPOS");
			}
			else if (a == 3 || t == 3)
			{
				LDAU ldau(table, argv[2]);
				ldau.AddLDAU_default(argv[2]);
			}
			else
				return -1;
		}
		else if (check_t == 0 && check_a == 2 )
		{
			printf("TABLE NAME BEHIND \"-t\" IS NOT EXIST!\n");
			printf("Use default LDAU value!\n");
			if (a == 2 || t == 2)
			{
				LDAU ldau(nullptr, "INPOS");
				ldau.AddLDAU_default("INPOS");
			}
			else if (a == 3 || t == 3)
			{
				LDAU ldau(nullptr, argv[2]);
				ldau.AddLDAU_default(argv[2]);
			}
			else
				return -1;
		}
		else if (check_a == 1)
		{
			if (a == 2 || t == 2)
			{
				LDAU ldau(nullptr, "INPOS");
				ldau.AddLDAU_default("INPOS");
			}
			else if (a == 3 || t == 3)
			{
				LDAU ldau(nullptr, argv[2]);
				ldau.AddLDAU_default(argv[2]);
			}
			else
				return -1;
		}
		else
			return -1;
		return 0;
	}
	//Add default magnetic moment to incar according to element
	if (!strcmp("--imag", argv[1]) || !strcmp("--im", argv[1]))
	{
		int p_vec;
		vector<const char*> par_vec = {"-a", "-clm", "-sfm", "-sfmw", "-afm", "-nfm"};
		int check_vec = check_mulpara(argc, argv, par_vec, 0 , &p_vec);
		char mode[10];
		if (argc == 2 || check_vec != 1)
			strcpy(mode, "a");
		else if (check_vec == 1)
			int matched = sscanf(argv[p_vec], "-%4s", mode);
		int t;
		int check_table = check_para(argc, argv, "-t", 1 , &t);
		cout << mode <<endl;
		cout << p_vec <<endl;
		if (check_table == 1)
		{
			if (p_vec == 3)
			{
				magorder magorder(argv[2], "INCAR", argv[t + 1]);
				magorder.incar_generate(mode);
			}
			else
			{
				magorder magorder("INPOS", "INCAR", argv[t + 1]);
				magorder.incar_generate(mode);
			}		
		}
		else if (check_table == 0)
		{
			printf("TABLE NAME BEHIND \"-t\" IS NOT EXIST!\n");
			printf("Use default MAG value!\n");
			if (p_vec == 3)
			{
				magorder magorder(argv[2], "INCAR");
				magorder.incar_generate(mode);
			}
			else
			{
				magorder magorder("INPOS", "INCAR");
				magorder.incar_generate(mode);
			}	
		}
		else
		{
			if (p_vec == 3)
			{
				magorder magorder(argv[2], "INCAR");
				magorder.incar_generate(mode);
			}
			else
			{
				magorder magorder("INPOS", "INCAR");
				magorder.incar_generate(mode);
			}	
		}
		return 0;
	}
	//Add IVDW moment to incar according to element
	if (!strcmp("--ivdw", argv[1])|| !strcmp("--iv", argv[1]))
	{
		int ivdw_numb = 0;
		int a = 0;int t = 0;
		int check_a = check_para(argc, argv, "-a", 0 , &a); //The fourth parameter is the minimum number of input parameters required.
		int check_t = check_para(argc, argv, "-t", 1 , &t); //The fourth parameter is the minimum number of input parameters required.
		if (check_t == 0)
			return -1;
		else if (check_t == 1)
		{
			if (!strcmp(argv[t + 1], "d3b"))
				ivdw_numb = 1; //DFT-D3(BJ)
			else if (!strcmp(argv[t + 1], "d2"))
				ivdw_numb = 2; //DFT-D2(G)
			else if (!strcmp(argv[t + 1], "d3z"))
				ivdw_numb = 3; //DFT-D3(zero)
			else if (!strcmp(argv[t + 1], "b86"))
				ivdw_numb = 4; //optB86b-vdw
			else if (!strcmp(argv[t + 1], "b88"))
				ivdw_numb = 5; //optB88-vdw
			else if (!strcmp(argv[t + 1], "pbe"))
				ivdw_numb = 6; //optPBE-vdw
			else if (!strcmp(argv[t + 1], "rpbe"))
				ivdw_numb = 7; //revPBE-vdw
			else
				ivdw_numb = 0;
			addIVDW_table(ivdw_numb);	
		}
		else
		{
			if (a == 3)
				addIVDWdefault(argv[2]);
			else
				addIVDWdefault();
		}
		return 0;
	}
	//KPOINTS
	// VASPMATE --ka 8000 G
	if (!strcmp("--ka", argv[1]) || !strcmp("--kppra", argv[1]) || !strcmp("--kpta", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 2 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
				E_ikppra("INPOS", argv[p + 2][0], atoi(argv[p + 1]));
			else if (!strcmp(argv[p], "-par" ) && p == 3)
				E_ikppra(argv[2], argv[p + 2][0], atoi(argv[p + 1]));
			else
				return -1;
		}
		else if (check_par == 2)
		{
			if (argc == 2)
				E_ikppra("INPOS", 'G', 1000);
			else if (argc == 3)
				E_ikppra("INPOS", 'G', atoi(argv[2]));
			else if (argc == 4)
				E_ikppra("INPOS", argv[3][0], atoi(argv[2]));
			else if (argc == 5) //VASPMATE --ka INPOS 8000 G
				E_ikppra(argv[2], argv[4][0], atoi(argv[3]));
			else
				return -1;
		}
		return 0;
	}
	// VASPMATE --kv 0.5 G 
	if (!strcmp("--kv", argv[1]) || !strcmp("--kspac", argv[1]) || !strcmp("--kptv", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 2 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
				E_ikspac("INPOS", argv[p + 2][0], atof(argv[p + 1]));
			else if (!strcmp(argv[p], "-par" ) && p == 3)
				E_ikspac(argv[2], argv[p + 2][0], atof(argv[p + 1]));
			else
				return -1;
		}
		else if (check_par == 2)
		{
			if (argc == 2)
				E_ikspac("INPOS", 'G', 0.5);
			else if (argc == 3)
				E_ikspac("INPOS", 'G', atof(argv[2]));
			else if (argc == 4)
				E_ikspac("INPOS", argv[3][0], atof(argv[2]));
			else if (argc == 5) //VASPMATE --kv INPOS 0.5 G
				E_ikspac(argv[2], argv[4][0], atof(argv[3]));
			else
				return -1;
		}
		return 0;
	}
	//VASPMATE --k 1 1 1 G 
	if (!strcmp("--k", argv[1]) || !strcmp("--kmesh", argv[1]) || !strcmp("--km", argv[1]) || !strcmp("--kpt", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (argc == 7)
			{
				int imesh[3] = { atoi(argv[p + 1]),atoi(argv[p + 2]) ,atoi(argv[p + 3]) };
				w_kpt(argv[p + 4][0], imesh);
			}			
			else 
			{
				return -1;
			}
		}
		else if (check_par == 2)
		{
			if (argc == 6)
			{
				int imesh[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
				w_kpt(argv[5][0], imesh);
			}
			else
			{
				return -1;
			}
		}
		return 0;
	}
	// VASPMATE --pot 
	// VASPMATE --pot -PBE 
	// VASPMATE --pot -PBE s
	if (!strcmp("--pot", argv[1]))
	{
		vector<string> label;
		if (argc == 2)
			pot_merge("INPOS", "-PBE", label);
		else if (argc == 3)
			pot_merge("INPOS", argv[2], label);
		else if (argc > 3)
		{
			int t; int suff ;
			int check_t = check_para(argc, argv, "-t", 0 , &t); //The fourth parameter is the minimum number of input parameters required.
			int check_suff = check_para(argc, argv, "-suff", 0 , &suff); //The fourth parameter is the minimum number of input parameters required.
			if (check_suff == 0 || check_t == 0)
				printf("Error: Too few arguments, please enter the correct number of parameters!\n");
			else if (check_t == 1 && check_suff == 1)
			{
				if (t == (argc - 3)) //VASPMATE --pot -t PBE -suff (this is wrong!!)
					pot_merge("INPOS", argv[t + 1], label);
				else //VASPMATE --pot -t PBE -suff sv
				{
					if (t == 2) //VASPMATE --pot -t PBE -suff sv GW
					{
						for (int i = suff + 1; i < argc; i++)
							label.push_back(argv[i]);
						pot_merge("INPOS", argv[t + 1], label);	
					}
					else if (t == 3) //VASPMATE --pot INPOS -t PBE -suff s
					{
						for (int i = suff + 1; i < argc; i++)
							label.push_back(argv[i]);
						pot_merge(argv[2], argv[t + 1], label);	
					}
				}
			}			
			else if (check_t == 1 && check_suff == 2)
			{
				if (t == 2) //VASPMATE --pot -t PBE
					pot_merge("INPOS", argv[t + 1], label);
				else //VASPMATE --pot INPOS -t PBE
					pot_merge(argv[2], argv[t + 1], label);	
			}
			else if (check_t == 2 && check_suff == 1) 
			{
				if (suff == 4) //VASPMATE --pot INPOS -PBE -suff sv GW
				{
					for (int i = suff + 1; i < argc; i++)
						label.push_back(argv[i]);
					pot_merge(argv[2], argv[3], label);	
				}
				else if (suff == 3) //VASPMATE --pot -PBE -suff s
				{
					for (int i = suff + 1; i < argc; i++)
						label.push_back(argv[i]);
					pot_merge("INPOS", argv[2], label);	
				}
			}
			else if (check_t == 2 && check_suff == 2) //VASPMATE --pot INPOS PBE
				pot_merge(argv[2], argv[3], label);
		}
		else
			return -1;
		return 0;
	}
	// VASPMATE --pote -type [elem1 elem2...] (-suff postfix1 postfix2…)
	// VASPMATE --pote -t [type] -e [elem1 elem2...] (-suff postfix1 postfix2…)
	if (!strcmp("--pote", argv[1]))
	{
		vector<string> label;
		vector<string> element;
		int t; int suff; int e;
		int check_t = check_para(argc, argv, "-t", 0 , &t); 
		int check_suff = check_para(argc, argv, "-suff", 0 , &suff); 
		int check_ele = check_para(argc, argv, "-e", 0 , &e);
		if (check_suff == 0 || check_t == 0 || check_ele == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_t == 1 && check_suff == 1 && check_ele == 1)
		{
			if (strcmp(argv[t + 1], "-PBE") && strcmp(argv[t + 1], "-LDA") && strcmp(argv[t + 1], "-GGA") && strcmp(argv[t + 1], "PBE") && strcmp(argv[t + 1], "LDA")  && strcmp(argv[t + 1], "GGA"))
			{ //VASPMATE --pote -t -e B N C -suff sv GW
				for (int i = e + 1; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else
			{ //VASPMATE --pote -t GGA -e B N C -suff sv GW
				for (int i = e + 1; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, argv[t + 1], label);
			}
		}
		else if (check_t == 1 && check_suff == 2 && check_ele == 1)
		{
			if (strcmp(argv[3], "-PBE") && strcmp(argv[3], "-LDA") && strcmp(argv[3], "-GGA") && strcmp(argv[3], "PBE") && strcmp(argv[3], "LDA")  && strcmp(argv[3], "GGA"))
			{ //VASPMATE --pote -t -e B N C
				for (int i = e + 1; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else
			{ //VASPMATE --pote -t PBE -e B N C
				for (int i = e + 1; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, argv[t + 1], label);
			}				
		}
		else if (check_t == 1 && check_suff == 1 && check_ele == 2) //VASPMATE --pote -t PBE B N C -suff sv GW
		{
			if (strcmp(argv[2], "-PBE") && strcmp(argv[2], "-LDA") && strcmp(argv[2], "-GGA") && strcmp(argv[2], "PBE") && strcmp(argv[2], "LDA")  && strcmp(argv[2], "GGA"))
			{ //VASPMATE --pote -t B N C -suff sv GW
				for (int i = t + 2; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else
			{ //VASPMATE --pote -t PBE B N C -suff sv GW
				for (int i = t + 2; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, argv[2], label);
			}			
		}
		else if (check_t == 1 && check_suff == 2 && check_ele == 2) //VASPMATE --pote -t PBE B N C -suff sv GW
		{
			if (strcmp(argv[2], "-PBE") && strcmp(argv[2], "-LDA") && strcmp(argv[2], "-GGA") && strcmp(argv[2], "PBE") && strcmp(argv[2], "LDA")  && strcmp(argv[2], "GGA"))
			{ //VASPMATE --pote -t B N C -suff sv GW
				for (int i = t + 2; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else
			{ //VASPMATE --pote -t PBE B N C -suff sv GW
				for (int i = t + 2; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, argv[2], label);
			}
		}
		else if (check_t == 2 && check_suff == 1 && check_ele == 1)
		{
			if (strcmp(argv[2], "-PBE") && strcmp(argv[2], "-LDA") && strcmp(argv[2], "-GGA") && strcmp(argv[2], "PBE") && strcmp(argv[2], "LDA")  && strcmp(argv[2], "GGA"))
			{
				for (int i = e + 1; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else //VASPMATE --pote -PBE -e B N C -suff sv GW
			{
				for (int i = e + 1; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, argv[2], label);
			}
		}
		else if (check_t == 2 && check_suff == 2 && check_ele == 1)
		{
			if (strcmp(argv[2], "-PBE") && strcmp(argv[2], "-LDA") && strcmp(argv[2], "-GGA") && strcmp(argv[2], "PBE") && strcmp(argv[2], "LDA")  && strcmp(argv[2], "GGA"))
			{
				for (int i = e + 1; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else //VASPMATE --pote -PBE -e B N C
			{
				for (int i = e + 1; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, argv[2], label);
			}
		}
		else if (check_t == 2 && check_suff == 2 && check_ele == 2)
		{
			if (strcmp(argv[2], "-PBE") && strcmp(argv[2], "-LDA") && strcmp(argv[2], "-GGA") && strcmp(argv[2], "PBE") && strcmp(argv[2], "LDA")  && strcmp(argv[2], "GGA"))
			{ //VASPMATE --pote B N C
				for (int i = 2; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else //VASPMATE --pote -PBE B N C
			{
				for (int i = 3; i < argc; i++)
					element.push_back(argv[i]);
				pot_merge_element(element, argv[2], label);
			}
		}
		else if (check_t == 2 && check_suff == 1 && check_ele == 2)
		{
			if (strcmp(argv[2], "-PBE") && strcmp(argv[2], "-LDA") && strcmp(argv[2], "-GGA") && strcmp(argv[2], "PBE") && strcmp(argv[2], "LDA")  && strcmp(argv[2], "GGA"))
			{
				for (int i = 3; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, "-PBE", label);
			}
			else //VASPMATE --pote -PBE B N C -suff sv GW
			{
				for (int i = 3; i < suff; i++)
					element.push_back(argv[i]);
				for (int i = suff + 1; i < argc; i++)
					label.push_back(argv[i]);
				pot_merge_element(element, argv[2], label);
			}
		}
		return 0;
	}
	// VASPMATE --check
	if (!strcmp("--check", argv[1]))
	{
		if (!strcmp("--check", argv[1]))
    	{
			if(argc == 2)
				check();
			else if(!strcmp("-in", argv[2]))
			{
				int filenumb = 0;
				if(argc == 3)
					check();
				else if(!strcmp("imp", argv[3])){
					filenumb = 1;
					check_one(filenumb);
				}
				else if(!strcmp("pos", argv[3])){
					filenumb = 2;
					check_one(filenumb);
				}
				else if(!strcmp("pot", argv[3])){
					filenumb = 3;
					check_one(filenumb);
				}
				else if(!strcmp("kpt", argv[3])){
					filenumb = 4;
					check_one(filenumb);
				}
				else
					return -1;
			}
			else if(!strcmp("-imp", argv[2]))
				check_one(1);
			else if(!strcmp("pos", argv[2]))
				check_one(2);
			else if(!strcmp("pot", argv[2]))
				check_one(3);
			else if(!strcmp("kpt", argv[2]))
				check_one(4);
			else if(!strcmp("-out", argv[2]))
			{
				if(argc == 3)
					RelaxJudgeConvergence();
				else if(!strcmp("rlx", argv[3]))
					RelaxJudgeConvergence();
				else if(!strcmp("stc", argv[3]))
					StaticJudgeConvergence();
				else if(!strcmp("md", argv[3]))
					mdJudgeConvergence();
				else
					return -1;			
			}
			else if(!strcmp("rlx", argv[2]))
				RelaxJudgeConvergence();
			else if(!strcmp("stc", argv[2]))
				StaticJudgeConvergence();
			else if(!strcmp("md", argv[2]))
				printf("This feature is still in the testing phase, please stay tuned for subsequent updates!\n");
			else
				return -1;
    	}
		return 0;
	}
	// VASPMATE --std3d file1(INPOS) file2(STD3POS)
	if (!strcmp("--std3d", argv[1]))
	{
		if (argc == 2)
			bandskpt_3d("INPOS", "STD3POS");
		else if (argc == 3)
			bandskpt_3d(argv[2], "STD3POS");
		else if (argc == 4)
			bandskpt_3d(argv[2], argv[3]);
		else
			return -1;
		return 0;
	}
	// VASPMATE --std2d file1(INPOS) file2(STD2POS)
	if (!strcmp("--std2d", argv[1]))
	{
		if (argc == 2)
			bandskpt_2d("INPOS", "STD2POS");
		else if (argc == 3)
			bandskpt_2d(argv[2], "STD2POS");
		else if (argc == 4)
			bandskpt_2d(argv[2], argv[3]);
		else
			return -1;
		return 0;
	}
	// VASPMATE --ka3d file1(INPOS) file2(STD3POS) kppra(20)
	if (!strcmp("--ka3d", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 1 , &p); 
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
				bandskpt_3d("INPOS", "STD3POS", atoi(argv[p + 1]));
			else if (!strcmp(argv[p], "-par" ) && p == 3)
				bandskpt_3d(argv[2], "STD3POS", atoi(argv[p + 1]));
			else if (!strcmp(argv[p], "-par" ) && p == 4)
				bandskpt_3d(argv[2], argv[3], atoi(argv[p + 1]));
			else
				return -1;
			return 0;
		}
		else if (check_par == 2)
		{
			if (argc == 2)
				bandskpt_3d("INPOS", "STD3POS", 20);
			else if (argc == 3)
				bandskpt_3d("INPOS", "STD3POS", atoi(argv[2]));
			else if (argc == 4)
				bandskpt_3d(argv[2], argv[3], 20);
			else if (argc == 5)
				bandskpt_3d(argv[2], argv[3], atoi(argv[4]));
			else
				return -1;
			return 0;
		}
	}
	if (!strcmp("--ka2d", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 1 , &p); 
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2)
				bandskpt_2d("INPOS", "STD2POS", atoi(argv[p + 1]));
			else if (!strcmp(argv[p], "-par" ) && p == 3)
				bandskpt_2d(argv[2], "STD2POS", atoi(argv[p + 1]));
			else if (!strcmp(argv[p], "-par" ) && p == 4)
				bandskpt_2d(argv[2], argv[3], atoi(argv[p + 1]));
			else
				return -1;
			return 0;
		}
		else if (check_par == 2)
		{		
			if (argc == 2)
				bandskpt_2d("INPOS", "STD2POS", 20);
			if (argc == 3)
				bandskpt_2d("INPOS", "STD2POS", atoi(argv[2]));
			if (argc == 4)
				bandskpt_2d(argv[2], argv[3], 20);
			if (argc == 5)
				bandskpt_2d(argv[2], argv[3], atoi(argv[4]));
			return 0;
		}
	}
	//HSE KPOINTS
	//VASPMATE --kahse 8000 0.05 G
	if (!strcmp(argv[1], "--kahse"))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); 
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2) //VASPMATE --kahse -par 8000 0.05 G
				HSE_mesh("INPOS", atoi(argv[p + 1]), 0, NULL, "ka", atof(argv[p + 2]), argv[p + 3][0]);
			else if (!strcmp(argv[p], "-par" ) && p == 2) //VASPMATE --kahse INPOS -par 8000 0.05 G
				HSE_mesh(argv[2], atoi(argv[p + 1]), 0, NULL, "ka", atof(argv[p + 2]), argv[p + 3][0]);
		}
		else if (check_par == 2)
		{
			if (argc == 2)
				HSE_mesh("INPOS", 8000, 0, NULL, "ka", 0.05, 'G');
			else if (argc == 4)
				HSE_mesh("INPOS", atoi(argv[2]), 0, NULL, "ka", atof(argv[3]), 'G');
			else if (argc == 5)
				HSE_mesh("INPOS", atoi(argv[2]), 0, NULL, "ka", atof(argv[3]), argv[4][0]);
			else if (argc == 6) //VASPMATE --kahse INPOS 8000 0.05 G
				HSE_mesh(argv[2], atoi(argv[3]), 0, NULL, "ka", atof(argv[4]), argv[5][0]);
		}
		return 0;
	}
	//VASPMATE --kvhse 0.5 0.05 G
	if (!strcmp(argv[1], "--kvhse"))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); 
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			if (!strcmp(argv[p], "-par" ) && p == 2) //VASPMATE --kvhse -par 8000 0.05 G
				HSE_mesh("INPOS", atoi(argv[p + 1]), 0, NULL, "kv", atof(argv[p + 2]), argv[p + 3][0]);
			else if (!strcmp(argv[p], "-par" ) && p == 2) //VASPMATE --kvhse INPOS -par 8000 0.05 G
				HSE_mesh(argv[2], atoi(argv[p + 1]), 0, NULL, "kv", atof(argv[p + 2]), argv[p + 3][0]);
		}
		else if (check_par == 2)
		{
			if (argc == 2)
				HSE_mesh("INPOS", 0, 0.5, NULL, "kv", 0.05, 'G');
			else if (argc == 4)
				HSE_mesh("INPOS", 0, atof(argv[2]), NULL, "kv", atof(argv[3]), 'G');
			else if (argc == 5)
				HSE_mesh("INPOS", 0, atof(argv[2]), NULL, "kv", atof(argv[3]), argv[4][0]);
			else if (argc == 6) //VASPMATE --kvhse INPOS -par 0.5 0.05 G
				HSE_mesh(argv[2], atoi(argv[3]), 0, NULL, "kv", atof(argv[4]), argv[5][0]);
		}
		return 0;
	}
	//VASPMATE --khse 1 1 1 0.05 G
	if (!strcmp(argv[1], "--kmhse"))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); 
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			int kpt[] = { atoi(argv[p + 1]), atoi(argv[p + 2]), atoi(argv[p + 3]) };
			if (!strcmp(argv[p], "-par" ) && p == 2) //VASPMATE --kmhse -par 7 7 1 0.05 G
				HSE_mesh("INPOS", 0, 0, kpt, "kpt", atof(argv[p + 4]), argv[p + 5][0]);

		}
		else if (check_par == 2)
		{
			if (argc == 2)
			{
				int kpt[] = { 1,1,1 };
				HSE_mesh("INPOS", 0, 0, kpt, "kpt", 0.05, 'G');
			}
			else if (argc == 6)
			{
				int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
				HSE_mesh("INPOS", 0, 0, kpt, "kpt", atof(argv[5]), 'G');
			}
			else if (argc == 7)
			{
				int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
				HSE_mesh("INPOS", 0, 0, kpt, "kpt", atof(argv[5]), argv[6][0]);
			}
		}
		return 0;
	}
	//	DOS
	if (!strcmp("--dos", argv[1]))
	{
		DOS dos;
		dos.GetDos(argc, argv);
		return 0;
	}
	// band
	if (!strcmp("--band", argv[1]))
	{
		BAND band;
		band.getband(argc, argv);
		return 0;
	}
	//effective mass
	if (!strcmp("--em", argv[1]))
	{
		EMC emc("EMC.in");
		if (argc == 5 && !strcmp("-b", argv[3]))
			emc.get_emc_option(atoi(argv[4]));
		else
			emc.get_emc_option();
	}
	//bader
	if (!strcmp("--bader", argv[1]))
	{
		CHGCAR chg;
		chg.operator_bader(argc, argv);
		return 0;
	}
	//neb
	if (!strcmp("--neb", argv[1]))
	{
		neb(argc, argv);
		return 0;
	}
	//wavefun
	if (!strcmp("--wfn", argv[1]))
	{
		wavefun(argc, argv);
		return 0;
	}
	//cdd
	if (!strcmp("--vcd", argv[1]))
	{
		CHGCAR chg;
		if (!strcmp("-split", argv[2]))
		{
			if (argc == 3)
				chg.readchgcar("CHGCAR");
			else
				chg.readchgcar(argv[3]);
			chg.writechgcar("CHGTOT.vasp", 0);
			chg.writespincar("CHGSPIN.vasp");
			chg.writespinUp_Dwcar("CHGSPIN_UP.vasp", "CHGSPIN_DW.vasp");
		}
		else if (!strcmp("-sum", argv[2]))
			chg.series_oper(argc, argv, 1);
		else if (!strcmp("-diff", argv[2]))
			chg.series_oper(argc, argv, 0);
		else
			return -1;
		return 0;
	}
	// Fermi surface 
	if (!strcmp("--fskv", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			return -1;
		else if (check_par == 1)
		{
			if (argc == 4)
				GetFermiMesh("INPOS", 0, atof(argv[p + 1]), NULL, "kv", 'G');
			if (argc == 5)
				GetFermiMesh("INPOS", 0, atof(argv[p + 1]), NULL, "kv", argv[p + 2][0]);
			if (argc == 6)
				GetFermiMesh(argv[2], 0, atof(argv[p + 1]), NULL, "kv", argv[p + 2][0]);
		}
		else if (check_par == 2)
		{
			if (argc == 3)
				GetFermiMesh("INPOS", 0, atof(argv[2]), NULL, "kv", 'G');
			if (argc == 4)
				GetFermiMesh("INPOS", 0, atof(argv[2]), NULL, "kv", argv[3][0]);
			if (argc == 5)
				GetFermiMesh(argv[2], 0, atof(argv[3]), NULL, "kv", argv[4][0]);
		}
		return 0;
	}
	if (!strcmp("--fska", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			if (argc == 4)
				GetFermiMesh("INPOS", atoi(argv[p + 1]), 0, NULL, "ka", 'G');
			if (argc == 5)
				GetFermiMesh("INPOS", atoi(argv[p + 1]), 0, NULL, "ka", argv[p + 2][0]);
			if (argc == 6)
				GetFermiMesh(argv[2], atoi(argv[p + 1]), 0, NULL, "ka", argv[p + 2][0]);
		}
		else if (check_par == 2)
		{
			if (argc == 3)
				GetFermiMesh("INPOS", atoi(argv[2]), 0, NULL, "ka", 'G');
			if (argc == 4)
				GetFermiMesh("INPOS", atoi(argv[2]), 0, NULL, "ka", argv[3][0]);
			if (argc == 5)
				GetFermiMesh(argv[2], atoi(argv[3]), 0, NULL, "ka", argv[4][0]);
		}
		return 0;
	}
	if (!strcmp("--fskm", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			int mesh[3] = { atoi(argv[p + 1]),atoi(argv[p + 2]) ,atoi(argv[p + 3]) };
			if (argc == 6)
				GetFermiMesh("INPOS", 0, 0, mesh, "kpt", 'G');
			if (argc == 7)
				GetFermiMesh("INPOS", 0, 0, mesh, "kpt", argv[p + 4][0]);
		}
		else if (check_par == 2)
		{
			int mesh[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
			if (argc == 5)
				GetFermiMesh("INPOS", 0, 0, mesh, "kpt", 'G');
			if (argc == 6)
				GetFermiMesh("INPOS", 0, 0, mesh, "kpt", argv[5][0]);
		}
		return 0;
	}
	if (!strcmp("--fsxd", argv[1]))
	{
		string LORBIT = GetInfoINCAR("LORBIT");
		if (LORBIT.length() == 0)
			LORBIT = "0";
		string ISPIN = GetInfoINCAR("ISPIN");
		if (ISPIN.length() == 0)
			ISPIN = "1";
		DOS dos;
		dos.readDoscar(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
		double efermi = dos.fermi_energy();
		EIGENVAL eigen("EIGENVAL", efermi);
		set<int> select_band_index;
		if (argc > 2 && !strcmp(argv[2], "-ib"))
		{
			for (int i = 3; i < argc; i++)
				select_band_index.insert(atoi(argv[i]));
		}
		eigen.TranEigenToXcrysden(select_band_index, efermi);
		return 0;
	}
	if (!strcmp("--fs", argv[1]))
	{
		string LORBIT = GetInfoINCAR("LORBIT");
		if (LORBIT.length() == 0)
			LORBIT = "0";
		string ISPIN = GetInfoINCAR("ISPIN");
		if (ISPIN.length() == 0)
			ISPIN = "1";
		DOS dos;
		dos.readDoscar(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
		double efermi = dos.fermi_energy();
		BAND band;
		band.readPROCAR(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
		vector<vector<vector<vector<double> > >	> ion_dos = band.getion_dos();
		vector<vector<vector<vector<double> > >	> ion_dos_up = band.getion_dos_up();
		vector<vector<vector<vector<double> > >	> ion_dos_dw = band.getion_dos_dw();
		EIGENVAL eigen("EIGENVAL", efermi);
		eigen.TranEigenToFermiSurface(efermi, LORBIT, ion_dos, ion_dos_up, ion_dos_dw, argc, argv);
		return 0;
	}
	// 3D bandstructure 
	if (!strcmp("--3dkv", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			if (argc == 4)
				Get3DbandMesh("INPOS", 0, atof(argv[p + 1]), NULL, "kv", 'G');
			if (argc == 5)
				Get3DbandMesh("INPOS", 0, atof(argv[p + 1]), NULL, "kv", argv[p + 2][0]);
			if (argc == 6)
				Get3DbandMesh(argv[2], 0, atof(argv[p + 1]), NULL, "kv", argv[p + 2][0]);
		}
		else if (check_par == 2)		
		{
			if (argc == 3)
				Get3DbandMesh("INPOS", 0, atof(argv[2]), NULL, "kv", 'G');
			if (argc == 4)
				Get3DbandMesh("INPOS", 0, atof(argv[2]), NULL, "kv", argv[3][0]);
			if (argc == 5)
				Get3DbandMesh(argv[2], 0, atof(argv[3]), NULL, "kv", argv[4][0]);
		}
		return 0;
	}
	if (!strcmp("--3dka", argv[1]))
	{
		int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			if (argc == 4)
				Get3DbandMesh("INPOS", atoi(argv[p + 1]), 0, NULL, "ka", 'G');
			if (argc == 5)
				Get3DbandMesh("INPOS", atoi(argv[p + 1]), 0, NULL, "ka", argv[p + 2][0]);
			if (argc == 6)
				Get3DbandMesh(argv[2], atoi(argv[p + 1]), 0, NULL, "ka", argv[p + 2][0]);
		}
		else if (check_par == 2)
		{
			if (argc == 3)
				Get3DbandMesh("INPOS", atoi(argv[2]), 0, NULL, "ka", 'G');
			if (argc == 4)
				Get3DbandMesh("INPOS", atoi(argv[2]), 0, NULL, "ka", argv[3][0]);
			if (argc == 5)
				Get3DbandMesh(argv[2], atoi(argv[3]), 0, NULL, "ka", argv[4][0]);
		}
		return 0;
	}
	if (!strcmp("--3dkm", argv[1]))
	{
			int p;
		int check_par = check_para(argc, argv, "-par", 0 , &p); //The fourth parameter is the minimum number of input parameters required.
		if (check_par == 0)
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
		else if (check_par == 1)
		{
			int mesh[3] = { atoi(argv[p + 1]),atoi(argv[p + 2]) ,atoi(argv[p + 3]) };
			if (argc == 6)
				Get3DbandMesh("INPOS", 0, 0, mesh, "kpt", 'G');
			if (argc == 7)
				Get3DbandMesh("INPOS", 0, 0, mesh, "kpt", argv[5][0]);
		}
		else if (check_par == 2)
		{
			int mesh[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
			if (argc == 5)
				Get3DbandMesh("INPOS", 0, 0, mesh, "kpt", 'G');
			if (argc == 6)
				Get3DbandMesh("INPOS", 0, 0, mesh, "kpt", argv[5][0]);
		}
		return 0;
	}
	if (!strcmp("--3dbs", argv[1]))
	{
		vector<int> band_index;
		if (argc > 2 && !strcmp(argv[2], "-ib"))
		{
			for (int i = 3; i < argc; i++)
				band_index.push_back(atoi(argv[i]));
		}
		string LORBIT = GetInfoINCAR("LORBIT");
		if (LORBIT.length() == 0)
			LORBIT = "0";
		string ISPIN = GetInfoINCAR("ISPIN");
		if (ISPIN.length() == 0)
			ISPIN = "1";
		DOS dos;
		dos.readDoscar(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
		double efermi = dos.fermi_energy();
		EIGENVAL eigen("EIGENVAL", efermi);
		eigen.get3dband(band_index, efermi);
		return 0;
	}
	//unfold band 
	//VASPMATE --unka 8000 0.05 G redef.in
	if (!strcmp(argv[1], "--unka"))
	{
		if (argc == 2)
			getunfoldkpoints("INPOS", 8000, 0, NULL, "ka", 0.05, 'G', "redef.in");
		else if (argc == 4)
			getunfoldkpoints("INPOS", atoi(argv[2]), 0, NULL, "ka", atof(argv[3]), 'G', "redef.in");
		else if (argc == 5)
			getunfoldkpoints("INPOS", atoi(argv[2]), 0, NULL, "ka", atof(argv[3]), argv[4][0], "redef.in");
		return 0;
	}
	//VASPMATE --unkv 0.5 0.05 G redef.in
	if (!strcmp(argv[1], "--unkv"))
	{
		if (argc == 2)
			getunfoldkpoints("INPOS", 0, 0.5, NULL, "kv", 0.05, 'G', "redef.in");
		else if (argc == 4)
			getunfoldkpoints("INPOS", 0, atof(argv[2]), NULL, "kv", atof(argv[3]), 'G', "redef.in");
		else if (argc == 5)
			getunfoldkpoints("INPOS", 0, atof(argv[2]), NULL, "kv", atof(argv[3]), argv[4][0], "redef.in");
		return 0;
	}
	//VASPMATE --unkm 1 1 1 0.05 G redef.in
	if (!strcmp(argv[1], "--unkm"))
	{
		if (argc == 2)
		{
			int kpt[] = { 1,1,1 };
			getunfoldkpoints("INPOS", 0, 0, kpt, "kpt", 0.05, 'G', "redef.in");
		}
		else if (argc == 6)
		{
			int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
			getunfoldkpoints("INPOS", 0, 0, kpt, "kpt", atof(argv[5]), 'G', "redef.in");
		}
		else if (argc == 7)
		{
			int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
			getunfoldkpoints("INPOS", 0, 0, kpt, "kpt", atof(argv[5]), argv[6][0], "redef.in");
		}
		return 0;
	}
	if (!strcmp(argv[1], "--unfold"))
	{
		string LORBIT = GetInfoINCAR("LORBIT");
		if (LORBIT.length() == 0)
			LORBIT = "0";
		string ISPIN = GetInfoINCAR("ISPIN");
		if (ISPIN.length() == 0)
			ISPIN = "1";
		DOS dos;
		dos.readDoscar(atoi(ISPIN.c_str()), atoi(LORBIT.c_str()));
		double efermi = dos.fermi_energy();
		EIGENVAL eigen("EIGENVAL", efermi);
		eigen.getunfold();
	}
	if (!strcmp("--thermo", argv[1]))
	{
		shermo(argc, argv);
		return 0;
	}
	//VASPMATE --del POSCAR POSCAR1 z 10 20
	if (!strcmp("--del", argv[1]))
	{
		if (argc != 7)
			return -1;
		FILE* fp1 = fopen(argv[2], "r");
		if (fp1 == NULL)
		{
			printf("%s is not exist!\n");
			return -1;
		}
		POSCAR pos;
		readposcar(fp1, pos);
		direct_to_carts(pos.iflg, pos.vec, pos.xyz, pos.nant);
		fclose(fp1);
		int axis = -1;
		if (!strcmp(argv[4], "x"))
			axis = 0;
		else if (!strcmp(argv[4], "y"))
			axis = 1;
		else if (!strcmp(argv[4], "y"))
			axis = 2;
		else
		{
			printf("Axis is wrong! Please check!\n");
			return -1;
		}
		double h1 = atof(argv[5]);
		double h2 = atof(argv[6]);
		int flag = 1;
		while (flag && pos.nant[0] != 0)
		{
			flag = 0;
			for (int i = 0; i < pos.nant[0]; i++)
			{
				if (pos.xyz[i][axis] > h1 && pos.xyz[i][axis] < h2)
				{
					flag = 1;
					pos = delete_atom(pos, i);
					break;
				}
			}
		}
		FILE* fp2 = fopen(argv[3], "w");
		savposcar(fp2, pos);
		fclose(fp2);
		return 0;
	}
	if (!strcmp("--en", argv[1]))
	{
		double energy = get_energy();
		FILE* fp = fopen("OUTCAR_ENERGY", "w");
		fprintf(fp,"%lf\n",energy);
		fclose(fp);
		return 0;
	}
	if (!strcmp("--enth", argv[1]))
	{
		return enth::operator_enth(argc, argv);
	}
	if (!strcmp("--vol", argv[1]))
	{
		char file[1024] = "NULL";
		if (argc == 2)
			strcpy(file, "INPOS");
		else
			strcpy(file, argv[2]);
		FILE* fp = fopen(file, "r");
		if (fp == NULL)
		{
			printf("No input file is found!\n");
			return 0;
		}
		POSCAR pos;
		readposcar(fp, pos);
		fclose(fp);
		double vol = volume(pos.vec);
		printf("%lf\n", vol);
		return 0;
	}
	if (!strcmp("--num", argv[1]))
	{
		char file[1024] = "NULL";
		if (argc == 2)
			strcpy(file, "INPOS");
		else
			strcpy(file, argv[2]);
		FILE* fp = fopen(file, "r");
		if (fp == NULL)
		{
			printf("No input file is found!\n");
			return 0;
		}
		POSCAR pos;
		readposcar(fp, pos);
		fclose(fp);
		FILE* fp1 = fopen("ATOM_NUM", "w");
		fprintf(fp1,"%d (total)", pos.nant[0]);
		printf("%d (total)", pos.nant[0]);
		for (int i = 0; i < pos.nant[1]; i++)
		{
			fprintf(fp1, " %d %s", pos.typenum[i], pos.elemsym[i]);
			printf(" %d %s", pos.typenum[i], pos.elemsym[i]);
		}
		fclose(fp1);
		fprintf(fp1, "\n");
		printf("\n");
		printf("Written ATOM_NUM file!\n");
		return 0;
	}
	if (!strcmp("--cell", argv[1]))
	{
		char file[1024] = "NULL";
		if (argc == 2)
			strcpy(file, "INPOS");
		else
			strcpy(file, argv[2]);
		FILE* fp = fopen(file, "r");
		if (fp == NULL)
		{
			printf("No input file is found!\n");
			return 0;
		}
		POSCAR pos;
		readposcar(fp, pos);
		fclose(fp);
		double* cell_angle = getangle(pos.vec);
		double* cell_vectors = getcevec(pos.vec);
		FILE* fp_ = fopen("CELLCAR", "w");
		fprintf(fp_, "The lattice constant is:\n          %f\n" , pos.latt);
		fprintf(fp_, "The angles between the axes are:\n          α= %0.2f° , β= %0.2f° , γ= %0.2f°\n" , cell_angle[0] , cell_angle[1], cell_angle[2]);
		fprintf(fp_, "The lattice vectors are:\n          a= %f Å , b= %f Å , c= %f Å\n" , cell_vectors[0] , cell_vectors[1], cell_vectors[2] );
		fclose(fp_);
		delete[] cell_angle;
		delete[] cell_vectors;
		Output_to_terminal("CELLCAR");
		printf("Written CELLCAR file!\n");	
		return 0;
	}
	if (!strcmp("--atom", argv[1]))
	{
		char file[1024] = "NULL";
		if (argc == 2)
			strcpy(file, "INPOS");
		else
			strcpy(file, argv[2]);
		FILE* fp = fopen(file, "r");
		if (fp == NULL)
		{
			printf("No input file is found!\n");
			return 0;
		}
		POSCAR pos;
		readposcar(fp, pos);
		fclose(fp);
		FILE* fp_ = fopen("ATOMCAR", "w");
		fprintf(fp_, "In the %s file, there are %d types of elements with their quantities as follows:\n" , file, pos.nant[1]);
		for(int i = 0; i < pos.nant[1]; i++)
			fprintf(fp_, "   %-2s" , pos.elemsym[i]);
		fprintf(fp_, "\n");
		for(int i = 0; i < pos.nant[1]; i++)
			fprintf(fp_, "   %-2d" , pos.typenum[i]);
		fprintf(fp_, "\n\n");
		if(pos.iflg == 0)
			fprintf(fp_, "The type of coordinates is: Direct.\n\n");
		else if(pos.iflg == 1)
			fprintf(fp_, "The type of coordinates is: Cartes.\n\n");
		fprintf(fp_, "The coordinates of each atom are as follows:\n");
		int current_atom = 0;
		for (int i = 0; i < pos.nant[2]; ++i)
			for (int j = 0; j < pos.typenum[i]; ++j)
				if (current_atom < pos.nant[0])
				{
					fprintf(fp, "   %-2s  %16.9lf %16.9lf %16.9lf\n", 
						pos.xyz[current_atom][0], pos.xyz[current_atom][1], pos.xyz[current_atom][2], pos.elemsym[i]);
					current_atom++;
				}
		fclose(fp_);
		Output_to_terminal("ATOMCAR");
		printf("Written ATOMCAR file!\n");
		return 0;
	}
	if (!strcmp("--hull", argv[1]))
	{
		hull_point hp;
		hp.spa_convexhull();
		return 0;
	}
	if (!strcmp("--elas", argv[1]))
	{
		//stress
		if (!strcmp("-g", argv[2]) || !strcmp("-generate", argv[2]))
		{
			if (argc == 3)
			{
				elas::elastic elas;
				elas.generate();
			}
			else if (argc > 3)
			{
				vector<double> strain;
				for (int i = 2; i < argc; i++)
					strain.push_back(atof(argv[i]));
				elas::elastic elas(strain);
				elas.generate();
			}
		}
		else if (!strcmp("-d", argv[2]) || !strcmp("-derive", argv[2]))
		{
			if (argc == 3)
			{
				elas::elastic elas;
				elas.getstress();
				elas.calculate();
			}
			else if (argc > 3)
			{
				vector<double> strain;
				for (int i = 2; i < argc; i++)
					strain.push_back(atof(argv[i]));
				elas::elastic elas(strain);
				elas.getstress();
				elas.calculate();
			}
		}
	}
	if (!strcmp("--elae", argv[1]) && strcmp("ss", argv[2]) != 0 )
	{
		//energy
		if (!strcmp("-g", argv[2]) || !strcmp("-generate", argv[2]))
		{
			if (argc == 3)
			{
				elas::elastic_en elas_energy;
				elas_energy.generate();
			}
			else if (argc > 3)
			{
				vector<double> strain_energy;
				for (int i = 2; i < argc; i++)
					strain_energy.push_back(atof(argv[i]));
				elas::elastic_en elas_energy(strain_energy);
				elas_energy.generate();
			}
		}
		else if (!strcmp("-d", argv[2]) || !strcmp("-derive", argv[2]))
		{
			if (argc == 3)
			{
				elas::elastic_en elas_energy;
				elas_energy.getenergy();
				elas_energy.calculate();
			}
			else if (argc > 3)
			{
				vector<double> strain_energy;
				for (int i = 2; i < argc; i++)
					strain_energy.push_back(atof(argv[i]));
				elas::elastic_en elas_energy(strain_energy);
				elas_energy.getenergy();
				elas_energy.calculate();
			}
		}
		return 0;
	}
	else if (!strcmp("--elae", argv[1]) && strcmp("ss", argv[2]) == 0)
	{
		int numDL = atoi(argv[4]);
		if (!strcmp("-g", argv[3]) || !strcmp("-generate", argv[3]))
		{
			if (argc == 5)
			{
				elas_DL::elastic_DL elasDL(numDL);
				elasDL.generate();
			}
			else if (argc > 5)
			{
				vector<double> strain_energyDL;
				for (int i = 4; i < argc; i++)
					strain_energyDL.push_back(atof(argv[i]));
				elas_DL::elastic_DL elasDL(numDL, strain_energyDL);
				elasDL.generate();
			}
		}
		else if (!strcmp("-d", argv[3]) || !strcmp("-derive", argv[3]))
		{
			if (argc == 5)
			{
				elas_DL::elastic_DL elasDL(numDL);
				elasDL.getenergy();
				elasDL.calculate();	
			}
			else if (argc > 5)
			{
				vector<double> strain_energyDL;
				for (int i = 4; i < argc; i++)
					strain_energyDL.push_back(atof(argv[i]));
				elas_DL::elastic_DL elasDL(numDL, strain_energyDL);
				elasDL.getenergy();
				elasDL.calculate();	
			}
		}
		return 0;
	}
	if (!strcmp("--db", argv[1]))
	{
		int a; int b; int inc; int json;
		int check_b = check_para(argc, argv, "-b", 0 , &b);
		int check_a = check_para(argc, argv, "-a", 0 , &a);
		int check_json = check_para(argc, argv, "-js", 0 , &json);
		int check_inc = check_para(argc, argv, "-inc", "-include", 1 , &inc);
		vector<string> file_name;
		if(check_inc == 1)
		{
			for (int i = inc + 1; i < argc; i++)
				file_name.push_back(std::string(argv[i]));
		}
		_outcar::OUTCAR out;
		if (check_json == 1)
		{
			if (argc == 3)
			{
				out.output("log.sdata", "-a", file_name);
				TranLogdataToJson("log.sdata","log.json");
			}
			else if (argc == 4)
			{
				if (check_b == 1 && check_a != 1)
				{
					out.output("log.sdata", "-b", file_name);
					TranLogdataToJson("log.sdata","log.json");
				}
				else if (check_a == 1 && check_b != 1)
				{
					out.output("log.sdata", "-a", file_name);
					TranLogdataToJson("log.sdata","log.json");	
				}
				else
				{
					out.output("log.sdata", "-a", file_name);
					TranLogdataToJson("log.sdata", argv[2]);
				}
			}
			else if (argc > 4)
			{
				if(check_b == 1)
					out.output("log.sdata", "-b", file_name);
				else
					out.output("log.sdata", "-a", file_name);
				TranLogdataToJson("log.sdata", argv[2]);
			}
		}
		else
		{
			if (argc == 2)
				out.output("log.sdata", "-a", file_name);
			else if (argc == 3)
			{
				if (check_b == 1 && check_a != 1)
					out.output("log.sdata", "-b", file_name);
				else if (check_a == 1 && check_b != 1)
					out.output("log.sdata", "-a", file_name);
				else
					out.output(argv[2], "-a", file_name);
			}
			else if (argc > 3)
			{
				if(check_b == 1)
					out.output(argv[2], "-b", file_name);
				else
					out.output(argv[2], "-a", file_name);
			}
		}
		return 0;
	}
	if (!strcmp("--db2js", argv[1]))
	{
		if(argc == 3)
			TranLogdataToJson(argv[2]);
		else if(argc == 4)
			TranLogdataToJson(argv[2],argv[3]);
		else if(argc == 2)
			TranLogdataToJson();
		return 0;
	}
	if (!strcmp("--js2db", argv[1]))
	{
		if(argc == 3)
			TranJsonToLogdata(argv[2]);
		else if(argc == 4)
			TranJsonToLogdata(argv[2],argv[3]);
		else if(argc == 2)
			TranJsonToLogdata();
		return 0;
	}
	if (!strcmp("--dat2csv", argv[1]))
	{
		if(argc == 4)
			TranFileToCsv(argv[2],argv[3]);
		else
			return -1;
		return 0;
	}
	if (!strcmp("--csv2dat", argv[1]))
	{
		if(argc == 4)
			TranCsvToFile(argv[2],argv[3]);
		else
			return -1;
		return 0;
	}
	if (!strcmp("--dbplus", argv[1]))
	{
		int json; int check_json = check_para(argc, argv, "-js", 0 , &json);
		if( check_json == 1)
		{
			if(argc < 5)
				return -1;
			else
			{
				vector<string> filename;
				for(int i = 2; i < argc - 1; i++)
					filename.push_back(string(argv[i]));
				PlusJson(filename);
			}
		}
		else
		{
			if(argc < 4)
				return -1;
			else
			{
				vector<string> filename;
				for(int i = 3; i < argc; i++)
					filename.push_back(string(argv[i]));
				PlusLog(filename, argv[2]);
			}
		}
		return 0;
	}
	if (!strcmp("--opti", argv[1]))
	{
		auto diel = get_dielectric();
		if (!strcmp("-lop", argv[2]))
		{
			for (int i = 0; i < argc; i++)
			{
				if (!strcmp("-2d", argv[i]))
					get_linear_optical_spectrums_2d(diel);
				else if (!strcmp("-3d", argv[i]))
					get_linear_optical_spectrums_3d(diel);
			}
		}
		return 0;
	}
	if (!strcmp("--mole", argv[1]))
	{
		moledynamic mole;
		if (!strcmp("-pcf", argv[2]))
			mole.poltPCDAT();
		else if (!strcmp("-en", argv[2]))
			mole.poltEnergy();
		return 0;
	}
	//VASPMATE --amag file1(INPOS) file2(INCAR) -t [Table] -mode a/fm/afm/fimwyck/fimelem/afmwyck
	/*if (!strcmp("--amag", argv[1]))
	{
		int mo;
		int check_mode = check_para(argc, argv, "-mode", 1 , &mo); //The fourth parameter is the minimum number of input parameters required.
		char** mode=(char**)malloc(6*sizeof(char*));
		if (check_mode == 0)
		{
			printf("Error: Too few arguments, please enter the correct number of parameters!\n");
			return 0;
		}
		else if (check_mode == 2)
		{
			printf("Note: It has been detected that the mode for generating magnetic moments has not been set; by default, VASPMATE will be generated for all (using the “-a” method).\n");
			mode[0] = "a";
			mode[1] = NULL;
		}	
		else if (check_mode == 1)
		{
			if (!strcmp("a", argv[mo + 1]))
			{
				mode[0] = "a";
				mode[1] = NULL;
			}
			else
			{
				for (int i = mo + 1; i < argc; i++)
					mode[i - mo - 1] = argv[i];
				mode [argc - mo - 1] = NULL;
			}
		}
		int t;
		int check_table = check_para(argc, argv, "-t", 1 , &t); //The fourth parameter is the minimum number of input parameters required.
		if (check_table == 1)
			amagorder amagorder("INPOS", "INCAR", argv[t + 1], mode);
		else
			amagorder amagorder("INPOS", "INCAR", nullptr, mode);
		return 0;
	}
	*/
	if (!strcmp("--magn", argv[1]))
	{
		int p; int p_vec; int number_p;
		vector<const char*> par_vec = {"-a", "-sfm", "-fm", "-afm", "-nfm"};
		int check_vec = check_mulpara(argc, argv, par_vec, 0 , &p_vec);
		int check_num = check_number(argc, argv, &number_p);
		char mode[10];
		int matched = sscanf(argv[p_vec], "-%4s", mode);		
		if (argc == 2 || check_vec != 1)
			strcpy(mode, "a");
		int t;
		int check_table = check_para(argc, argv, "-t", 1 , &t);
		if (!strcmp("-g", argv[2]) || !strcmp("-generate", argv[2]))
		{
			int t;
			int check_table = check_para(argc, argv, "-t", 1 , &t);
			if (check_num == 2) // no number
			{
				if (check_table == 1)
				{
					magorder magorder("INPOS", "INCAR", argv[t + 1]);
					magorder.generate(mode, 1);
				}
				else
				{
					magorder magorder;
					magorder.generate(mode, 1);
				}
			}
			else if(check_num == 1)
			{
				if (check_table == 1)
				{
					magorder magorder("INPOS", "INCAR", argv[t + 1]);
					magorder.generate(mode, atoi(argv[number_p]));
				}
				else
				{
					magorder magorder("INPOS", "INCAR", nullptr);
					magorder.generate(mode, atoi(argv[number_p]));
				}
			}
		}
		else if (!strcmp("-d", argv[2]) || !strcmp("-derive", argv[2]))
		{
			magorder magorder;
			magorder.derive();
		}
		return 0;
	}
	if (!strcmp("--mcoup", argv[1]))
	{
		if (!strcmp("-g", argv[2]) || !strcmp("-generate", argv[2]))
		{
			Magcouple Magcouple;
			Magcouple.generate();
		}
		else if (!strcmp("-d", argv[2]) || !strcmp("-derive", argv[2]))
		{
			Magcouple Magcouple;
			Magcouple.derive();
		}
		return 0;
	}
	if (!strcmp("--clean", argv[1]))
	{
		return clean_operat(argc, argv);
	}
	if (!strcmp("--mds", argv[1]))
	{
		return aimd::operator_aimd(argc, argv);
	}
	if (!strcmp("--dbs", argv[1]))
	{
		return VMdb::operator_db(argc, argv);
	}
	if (!strcmp("--vel", argv[1]))
	{
		return vel::operator_vel(argc, argv);
	}
	if (!strcmp("--db2s", argv[1]))
	{
		try{
			_RaW4db::RaW4db VASPMATE_collect;
			VASPMATE_collect.operator_db_collect(argc, argv);
		}catch(const logic_error& e)
		{
			cerr << e.what() << endl;
		}
		return 0;
	}
	return 1;
}
