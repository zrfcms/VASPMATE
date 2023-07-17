/*===============================================================
*   Copyright[c] 2022-2023, Z. C. Pan and R. F. Zhang
*
*   This file is part of the
*   VASPMATE - An efficient program for high - throughput first principles 
	computations as partner of VASP code
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*================================================================*/
#include"../include/plotWave.h"
#include"../include/structure_operator.h"
#include"../include/potcar.h"
#include"../include/bskpt.h"
#include"../include/dos.h"
#include"../include/incar.h"
#include"../include/band.h"
#include"../include/chgcar.h"
#include"../include/acneb.h"
#include"../include/cif.h"
#include"../include/Eigenval.h"
#include"../include/shermo.h"
#include"../include/outcar.h"
#include"../include/potential.h"
#include"../include/spa_plot.h"
#include"../include/emc.h"

int main(int argc, char* argv[])
{
	//start
	if (argc == 1)
	{
		printf("VASPMATE Version 1.1.0 (2023.4.1)\n");
		printf("An efficient program for high - throughput first principles computations as partner of VASP code.\n");
		printf("Copyright[c] 2022 - 2023, Beihang University, by Zhaocheng Pan and Ruifeng Zhang\n");
		printf("Please send bugsand suggestions to zrfcms@buaa.edu.cn\n");
		return 0;
	}
	//version
	if (!strcmp("-v", argv[1]) || !strcmp("-version", argv[1]))
	{
		printf(" VASPMATE version 1.1.0 (2023.4.1)\n\n");
		printf(" VASPMATE - an efficient program for simple and efficient\n");
		printf(" high-throughput computations as partner of VASP code.\n");
		printf(" Copyright[c] 2022 - 2023, Beihang University, by Z. C. Pan and R. F. Zhang\n");
		return 0;
	}
	//license
	if (!strcmp("-l", argv[1]) || !strcmp("-license", argv[1]))
	{
		printf(" VASPMATE - an efficient program for simple and efficient\n");
		printf(" high-throughput computations as partner of VASP code.\n");
		printf(" Copyright[c] 2022 - 2023, Beihang University, by Z. C. Pan and R. F. Zhang\n");
		printf(" ---------------------version 1.0.0--------------------- \n");
		printf(" This program is currently copyrighted and distributed free of charge for academic,\n");
		printf(" scientific and educational and non-commercial users with our permission. You are\n");
		printf(" welcome to redistribute it under certain conditions with our permission. Part of\n");
		printf(" these terms may be changed without prior announcement. This program is provided \n");
		printf(" \"as is\" without any expressed or implied warranty.\n");
	}
	// VASPMATE --cif2pos file1(cif) file2(poscar)
	if (!strcmp("--cif2pos", argv[1]))
	{
		if (argc < 4)
			return -1;
		TranCIFToPOSCAR(argv[2], argv[3]);
		return 0;
	}
	if (!strcmp("--pos2cif", argv[1]))
	{
		if (argc < 4)
			return -1;
		TranPOSCARToCIF(argv[2], argv[3]);
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
	// VASPMATE --prim (INPOS) (PRIMPOS)
	// VASPMATE --prim file1 file2
	if (!strcmp("--prim", argv[1]))
	{
		if (argc == 2)
			EV_primitive("INPOS", "PRIMPOS");
		if (argc == 4)
			EV_primitive(argv[2], argv[3]);
	}
	// VASPMATE --unit (INPOS) (UNITPOS)
	// VASPMATE --unit file1 file2
	if (!strcmp("--unit", argv[1]))
	{
		if (argc == 2)
			EV_unitcell("INPOS", "UNITPOS");
		if (argc == 4)
			EV_unitcell(argv[2], argv[3]);
		return 0;
	}

	// VASPMATE --symm file
	if (!strcmp("--symm", argv[1]) && argc == 3)
	{
		EV_get_symmetry(argv[2]);
		return 0;
	}
	if (!strcmp("--symm", argv[1]) && argc == 2)
	{
		EV_get_symmetry("INPOS");
		return 0;
	}

	// VASPMATE --super file1 file2 a b c
	// VASPMATE --super (INPOS) (SUPERPOS) a b c
	if (!strcmp("--super", argv[1]))
	{
		if (argc == 7)
		{
			int super[3] = { atoi(argv[4]),atoi(argv[5]) ,atoi(argv[6]) };
			EV_supercell(argv[2], argv[3], super);
		}
		if (argc == 5)
		{
			int super[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
			EV_supercell("INPOS", "SUPERPOS", super);
		}
		return 0;
	}

	// VASPMATE --affine (AFFPOS0) -simshear(-purshear,-tension)  xx(yy zz xy xz yz) init_strain step_length step_num
	if (!strcmp("--affine", argv[1]) && argc == 7)
	{
		EV_affine("AFFPOS0", argv[2], argv[3], atof(argv[4]), atof(argv[5]), atoi(argv[6]));
		return 0;
	}

	// VASPMATE --alias (ALIPOS0) (ALIPOS_dz) -tensi mode(xx,yy,zz) istartz iendz ispacingz [zvalue]
	// VASPMATE --alias (ALIPOS0) (ALIPOS_dz) -shear mode(xy,xz,yz) istartx iendx ispacingx istarty iendy ispacingy [zvalue]
	if (!strcmp("--alias", argv[1]))
	{
		if (!strcmp("-tension", argv[2]) && argc == 8)
			EV_alias_tensile("ALIPOS0", argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
		if (!strcmp("-shear", argv[2]) && argc == 11)
			EV_alias_shear("ALIPOS0", argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]),
				atof(argv[8]), atof(argv[9]), atof(argv[10]));
		return 0;
	}

	// VASPMATE --proj -rot file1(INPOS) file2(PROJPOS) rotx rotx rotz
	// VASPMATE --proj -ind file1(INPOS) file2(PROJPOS) pvh pvk pvl uvu uvv uvw
	// VASPMATE --proj -mat file1(INPOS) file2(PROJPOS) mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33
	if (!strcmp("--proj", argv[1]))
	{
		if (!strcmp("-rot", argv[2]) && argc == 8)
		{
			double rot[3] = { atof(argv[5]),atof(argv[6]) ,atof(argv[7]) };
			EV_rotproj(argv[3], argv[4], rot);
		}
		if (!strcmp("-rot", argv[2]) && argc == 6)
		{
			double rot[3] = { atof(argv[3]),atof(argv[4]) ,atof(argv[5]) };
			EV_rotproj("INPOS", "PROJPOS", rot);
		}
		if (!strcmp("-ind", argv[2]) && argc == 11)
		{
			EV_indproj(argv[3], argv[4], atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]));
		}
		if (!strcmp("-ind", argv[2]) && argc == 9)
		{
			EV_indproj("INPOS", "PROJPOS", atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]));
		}
		if (!strcmp("-mat", argv[2]) && argc == 14)
		{
			double projmat[3][3] = { atof(argv[5]),atof(argv[6]) ,atof(argv[7]),
			 atof(argv[8]),atof(argv[9]) ,atof(argv[10]), atof(argv[11]),atof(argv[12]) ,atof(argv[13]) };
			EV_matproj(argv[3], argv[4], projmat);
		}
		if (!strcmp("-mat", argv[2]) && argc == 12)
		{
			double projmat[3][3] = { atof(argv[3]) ,atof(argv[4]),atof(argv[5]),atof(argv[6]) ,atof(argv[7]),
			 atof(argv[8]),atof(argv[9]) ,atof(argv[10]), atof(argv[11]) };
			EV_matproj("INPOS", "PROJPOS", projmat);
		}
		return 0;
	}
	if (!strcmp(argv[1], "--ads"))
	{
		if (argc == 6)
			adsorbent(argv[2], argv[3], argv[4], atof(argv[5]));
		return 0;
	}
	// VASPMATE --redef (INPOS) (REDPOS) vect11 vect12 vect13 vect21 vect22 vect23 vect31 vect32 vect33 a1 x1 a2 x2
	if (!strcmp("--redef", argv[1]))
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
		if (argc == 13)
		{
			int rot[3][3] = { atoi(argv[4]),atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
				atoi(argv[8]), atoi(argv[9]), atoi(argv[10]),atoi(argv[11]), atoi(argv[12]) };
			EV_redefine(argv[2], argv[3], rot, 0.001, 0, 0, 0, 0, 0);
		}
		if (argc == 11)
		{
			int rot[3][3] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]),atoi(argv[5]), atoi(argv[6]), atoi(argv[7]),
				atoi(argv[8]), atoi(argv[9]), atoi(argv[10]) };
			EV_redefine("INPOS", "REDPOS", rot, 0.001, 0, 0, 0, 0, 0);
		}
		if (argc == 15)
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
			EV_redefine("INPOS", "REDPOS", rot, 0.001, a1, x1, a2, x2, 1);
		}
		return 0;
	}
	// VASPMATE --ieee file1(INPOS) file2(IEEEPOS)
	if (!strcmp("--ieee", argv[1]))
	{
		if (argc == 4)
			EV_recell(argv[2], argv[3]);
		if (argc == 2)
			EV_recell("INPOS", "IEEEPOS");
		return 0;
	}
	// VASPMATE --fixc file2(INPOS) file2(FIXPOS) axis m n F F F (fix m-n)
	if (!strcmp("--fixc", argv[1]))
	{
		// VASPMATE --fixc axis m n F F F
		if (argc == 8)
		{
			char fix[3] = { argv[5][0],argv[6][0],argv[7][0] };
			EV_fixatomcoor("INPOS", "FIXPOS", argv[2], atof(argv[3]), atof(argv[4]), fix);
		}
		//VASPMATE --fixc file1 file2 axis m n F F F
		if (argc == 10)
		{
			char fix[3] = { argv[7][0],argv[8][0],argv[9][0] };
			EV_fixatomcoor(argv[2], argv[3], argv[4], atof(argv[5]), atof(argv[6]), fix);
		}
		return 0;
	}
	// VASPMATE --fixa file2 file2 1 2 3 ... F F F 
	if (!strcmp("--fixa", argv[1]))
	{
		char fix[3] = { argv[argc - 3][0],argv[argc - 2][0],argv[argc - 1][0] };
		EV_fixatomindex(argv[2], argv[3], argc, argv, fix);
	}
	// VASPMATE --fixe file2 file2 A B C... F F F 
	if (!strcmp("--fixe", argv[1]))
	{
		vector<string> ele;
		for (int i = 4; i < argc - 3; i++)
			ele.push_back(argv[i]);
		char fix[3] = { argv[argc - 3][0],argv[argc - 2][0],argv[argc - 1][0] };
		EV_fixatomele(argv[2], argv[3], ele, fix);
		return 0;
	}
	// VASPMATE --ufix file
	if (!strcmp("--ufix", argv[1]))
	{
		if (argc == 2)
			EV_cleanfix("INPOS");
		if (argc == 3)
			EV_cleanfix(argv[2]);
		return 0;
	}
	// VASPMATE --cartes file1(INPOS) file2(NEWPOS)
	if (!strcmp("--cartes", argv[1]))
	{
		if (argc == 2)
			EV_direct_to_carts("INPOS", "NEWPOS");
		if (argc == 4)
			EV_direct_to_carts(argv[2], argv[3]);
		return 0;
	}
	// VASPMATE --direct file1(INPOS) file2(NEWPOS)
	if (!strcmp("--direct", argv[1]))
	{
		if (argc == 2)
			EV_carts_to_direct("INPOS", "NEWPOS");
		if (argc == 4)
			EV_carts_to_direct(argv[2], argv[3]);
		return 0;
	}
	// VASPMATE --move (-c/-d) file1(INPOS) file2(MOVEPOS) axis mmin mmax dx dy dz 
	if (!strcmp("--move", argv[1]))
	{
		char label[3];
		if (!strcmp("-c", argv[2]))
			strcpy(label, "-c");
		else if (!strcmp("-d", argv[2]))
			strcpy(label, "-d");
		//VASPMATE --move axis mmin mmax dx dy dz 
		if (argc == 8) {
			double move[3] = { atof(argv[5]),atof(argv[6]),atof(argv[7]) };
			EV_move("INPOS", "MOVEPOS", label, argv[2], move, atof(argv[3]), atof(argv[4]));
		}
		//VASPMATE --move -c/-d axis mmin mmax dx dy dz 
		if (argc == 9) {
			double move[3] = { atof(argv[6]),atof(argv[7]),atof(argv[8]) };
			EV_move("INPOS", "MOVEPOS", label, argv[3], move, atof(argv[4]), atof(argv[5]));
		}
		//VASPMATE --move file1 file2 axis mmin mmax dx dy dz 
		if (argc == 10) {
			double move[3] = { atof(argv[7]),atof(argv[8]),atof(argv[9]) };
			EV_move(argv[2], argv[3], label, argv[4], move, atof(argv[5]), atof(argv[6]));
		}
		//VASPMATE --move -c/-d file1 file2 axis mmin mmax dx dy dz 
		if (argc == 11) {
			double move[3] = { atof(argv[8]),atof(argv[9]),atof(argv[10]) };
			EV_move(argv[3], argv[4], label, argv[5], move, atof(argv[6]), atof(argv[7]));
		}
		return 0;
	}
	// VASPMATE --sortc file1(INPOS) file2(SORTPOS) mode
	if (!strcmp("--sortc", argv[1]))
	{
		if (argc == 3)
			EV_atomsort_coord("INPOS", "SORTPOS", argv[2]);
		if (argc == 5)
			EV_atomsort_coord(argv[2], argv[3], argv[4]);
		return 0;
	}
	// VASPMATE --sorte file1(INPOS) file2(SORTPOS) ele1 ele2
	if (!strcmp("--sorte", argv[1]))
	{
		if (argc == 4)
			EV_atomsort_element("INPOS", "SORTPOS", argv[2], argv[3]);
		if (argc == 6)
			EV_atomsort_element(argv[2], argv[3], argv[4], argv[5]);
		return 0;
	}
	//INCAR
	if (!strcmp("--i", argv[1]))
	{
		if (!strcmp(argv[2], "-pcd"))
			pcd_model(argc, argv);
		else
			write_INCAR(argc, argv);
		return 0;
	}
	if (!strcmp("--i_append", argv[1]) || !strcmp("--i_app", argv[1]))
	{
		INCAR_append(argc, argv);
		return 0;
	}
	if (!strcmp("--i_delete", argv[1]) || !strcmp("--i_del", argv[1]))
	{
		INCAR_delete(argc, argv);
		return 0;
	}
	if (!strcmp("--i_remove", argv[1]) || !strcmp("--i_rem", argv[1]))
	{
		INCAR_remove(argc, argv);
		return 0;
	}
	if (!strcmp("--i_replace", argv[1]) || !strcmp("--i_rep", argv[1]))
	{
		INCAR_replace(argc, argv);
		return 0;
	}
	//KPOINTS
	// VASPMATE --ka 8000 G
	if (!strcmp("--ka", argv[1]) || !strcmp("--kppra", argv[1]) || !strcmp("--kpta", argv[1]))
	{
		if (argc == 2)
			E_ikppra("INPOS", 'G', 1000);
		else if (argc == 3)
			E_ikppra("INPOS", 'G', atoi(argv[2]));
		else if (argc == 4)
			E_ikppra("INPOS", argv[3][0], atoi(argv[2]));
		return 0;
	}
	// VASPMATE --kv 0.5 G 
	if (!strcmp("--kv", argv[1]) || !strcmp("--kspac", argv[1]) || !strcmp("--kptv", argv[1]))
	{
		if (argc == 2)
			E_ikspac("INPOS", 'G', 0.5);
		else if (argc == 3)
			E_ikspac("INPOS", 'G', atof(argv[2]));
		else if (argc == 4)
			E_ikspac("INPOS", argv[3][0], atof(argv[2]));
		return 0;
	}
	//VASPMATE --k 1 1 1 G 
	if (!strcmp("--k", argv[1]) || !strcmp("--kmesh", argv[1]) || !strcmp("--km", argv[1]))
	{
		if (argc == 6)
		{
			int imesh[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
			w_kpt(argv[5][0], imesh);
		}
		else
		{
			int mesh[3] = { 1,1,1 };
			w_kpt('G', mesh);
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
			pot_merge("POSCAR", const_cast<char*>("-PBE"), label);
		if (argc == 3)
			pot_merge("POSCAR", argv[2], label);
		if (argc > 3)
		{
			for (int i = 3; i < argc; i++)
				label.push_back(argv[i]);
			pot_merge("POSCAR", argv[2], label);
		}
		return 0;
	}
	// VASPMATE --pote B N
	// VASPMATE --pote -PBE B N
	// VASPMATE --pote -PBE B N sv 
	if (!strcmp("--pote", argv[1]))
	{
		vector<string> label;
		vector<string> element;
		vector<string> postfix = { "s","d","h","sv","pv","GW" };
		int tail = argc;
		if (strcmp(argv[2], "-PBE") && strcmp(argv[2], "-LDA") && strcmp(argv[2], "-GGA"))
		{
			for (int i = 2; i < argc; i++)
			{
				if (find(postfix.begin(), postfix.end(), argv[i]) != postfix.end())
					tail = i;
			}
			for (int i = 2; i < tail; i++)
				element.push_back(argv[i]);
			for (int i = tail; i < argc; i++)
				label.push_back(argv[i]);
			pot_merge_element(element, const_cast<char*>("-PBE"), label);
		}
		else
		{
			for (int i = 3; i < argc; i++)
			{
				if (find(postfix.begin(), postfix.end(), argv[i]) != postfix.end())
					tail = i;
			}
			for (int i = 3; i < tail; i++)
				element.push_back(argv[i]);
			for (int i = tail; i < argc; i++)
				label.push_back(argv[i]);
			pot_merge_element(element, argv[2], label);
		}
		return 0;
	}
	// VASPMATE --check
	if (!strcmp("--check", argv[1]))
		check();
	// VASPMATE --kpt3d file1(INPOS) file2(PRIMPOS) kppra(20)
	if (!strcmp("--kpt3d", argv[1]))
	{
		if (argc == 2)
			bandskpt_3d("INPOS", "PRIMPOS", 20);
		if (argc == 3)
			bandskpt_3d("INPOS", "PRIMPOS", atoi(argv[2]));
		if (argc == 4)
			bandskpt_3d(argv[2], argv[3], 20);
		if (argc == 5)
			bandskpt_3d(argv[2], argv[3], atoi(argv[4]));
		return 0;
	}
	if (!strcmp("--kpt2d", argv[1]))
	{
		if (argc == 2)
			bandskpt_2d("INPOS", "PRIMPOS", 20);
		if (argc == 3)
			bandskpt_2d("INPOS", "PRIMPOS", atoi(argv[2]));
		if (argc == 4)
			bandskpt_2d(argv[2], argv[3], 20);
		if (argc == 5)
			bandskpt_2d(argv[2], argv[3], atoi(argv[4]));
		return 0;
	}
	//HSE KPOINTS
	//VASPMATE --kahse 8000 0.05 G
	if (!strcmp(argv[1], "--kahse"))
	{
		if (argc == 2)
			HSE_mesh("INPOS", 8000, 0, NULL, const_cast<char*>("ka"), 0.05, 'G');
		else if (argc == 4)
			HSE_mesh("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), atof(argv[3]), 'G');
		else if (argc == 5)
			HSE_mesh("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), atof(argv[3]), argv[4][0]);
		return 0;
	}
	//VASPMATE --kvhse 0.5 0.05 G
	if (!strcmp(argv[1], "--kvhse"))
	{
		if (argc == 2)
			HSE_mesh("INPOS", 0, 0.5, NULL, const_cast<char*>("kv"), 0.05, 'G');
		else if (argc == 4)
			HSE_mesh("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), atof(argv[3]), 'G');
		else if (argc == 5)
			HSE_mesh("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), atof(argv[3]), argv[4][0]);
		return 0;
	}
	//VASPMATE --khse 1 1 1 0.05 G
	if (!strcmp(argv[1], "--kmhse"))
	{
		if (argc == 2)
		{
			int kpt[] = { 1,1,1 };
			HSE_mesh("INPOS", 0, 0, kpt, const_cast<char*>("kpt"), 0.05, 'G');
		}
		else if (argc == 6)
		{
			int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
			HSE_mesh("INPOS", 0, 0, kpt, const_cast<char*>("kpt"), atof(argv[5]), 'G');
		}
		else if (argc == 7)
		{
			int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
			HSE_mesh("INPOS", 0, 0, kpt, const_cast<char*>("kpt"), atof(argv[5]), argv[6][0]);
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
		char file[1024] = "EMC.in";
		int idx = -1;
		for (int i = 2; i < argc; i++)
		{
			if (!strcmp("-f", argv[i]) && i + 1 < argc)
				strcpy(file, argv[i + 1]);
			if (!strcmp("-b", argv[i]) && i + 1 < argc)
				idx = atoi(argv[i + 1]);
		}
		EMC emc(file);
		idx == -1 ? emc.get_emc_option() : emc.get_emc_option(idx);
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
	if (!strcmp("--wfun", argv[1]))
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
		if (!strcmp("-sum", argv[2]))
			chg.series_oper(argc, argv, 1);
		else if (!strcmp("-diff", argv[2]))
			chg.series_oper(argc, argv, 0);
		return 0;
	}
	// Fermi surface 
	if (!strcmp("--fskv", argv[1]))
	{
		if (argc == 3)
			GetFermiMesh("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), 'G');
		if (argc == 4)
			GetFermiMesh("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), argv[3][0]);
		return 0;
	}
	if (!strcmp("--fska", argv[1]))
	{
		if (argc == 3)
			GetFermiMesh("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), 'G');
		if (argc == 4)
			GetFermiMesh("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), argv[3][0]);
		return 0;
	}
	if (!strcmp("--fskm", argv[1]))
	{
		int mesh[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
		if (argc == 5)
			GetFermiMesh("INPOS", 0, 0, mesh, const_cast<char*>("kpt"), 'G');
		if (argc == 6)
			GetFermiMesh("INPOS", 0, 0, mesh, const_cast<char*>("kpt"), argv[5][0]);
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
		if (argc > 2 && !strcmp(argv[2], "-b"))
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
		if (argc == 3)
			Get3DbandMesh("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), 'G');
		if (argc == 4)
			Get3DbandMesh("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), argv[3][0]);
		return 0;
	}
	if (!strcmp("--3dka", argv[1]))
	{
		if (argc == 3)
			Get3DbandMesh("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), 'G');
		if (argc == 4)
			Get3DbandMesh("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), argv[3][0]);
		return 0;
	}
	if (!strcmp("--3dkm", argv[1]))
	{
		int mesh[3] = { atoi(argv[2]),atoi(argv[3]) ,atoi(argv[4]) };
		if (argc == 5)
			Get3DbandMesh("INPOS", 0, 0, mesh, const_cast<char*>("kpt"), 'G');
		if (argc == 6)
			Get3DbandMesh("INPOS", 0, 0, mesh, const_cast<char*>("kpt"), argv[5][0]);
		return 0;
	}
	if (!strcmp("--3dbs", argv[1]))
	{
		vector<int> band_index;
		if (argc > 2 && !strcmp(argv[2], "-b"))
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
			getunfoldkpoints("INPOS", 8000, 0, NULL, const_cast<char*>("ka"), 0.05, 'G', "redef.in");
		else if (argc == 4)
			getunfoldkpoints("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), atof(argv[3]), 'G', "redef.in");
		else if (argc == 5)
			getunfoldkpoints("INPOS", atoi(argv[2]), 0, NULL, const_cast<char*>("ka"), atof(argv[3]), argv[4][0], "redef.in");
		return 0;
	}
	//VASPMATE --unkv 0.5 0.05 G redef.in
	if (!strcmp(argv[1], "--unkv"))
	{
		if (argc == 2)
			getunfoldkpoints("INPOS", 0, 0.5, NULL, const_cast<char*>("kv"), 0.05, 'G', "redef.in");
		else if (argc == 4)
			getunfoldkpoints("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), atof(argv[3]), 'G', "redef.in");
		else if (argc == 5)
			getunfoldkpoints("INPOS", 0, atof(argv[2]), NULL, const_cast<char*>("kv"), atof(argv[3]), argv[4][0], "redef.in");
		return 0;
	}
	//VASPMATE --unkm 1 1 1 0.05 G redef.in
	if (!strcmp(argv[1], "--unkm"))
	{
		if (argc == 2)
		{
			int kpt[] = { 1,1,1 };
			getunfoldkpoints("INPOS", 0, 0, kpt, const_cast<char*>("kpt"), 0.05, 'G', "redef.in");
		}
		else if (argc == 6)
		{
			int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
			getunfoldkpoints("INPOS", 0, 0, kpt, const_cast<char*>("kpt"), atof(argv[5]), 'G', "redef.in");
		}
		else if (argc == 7)
		{
			int kpt[] = { atoi(argv[2]), atoi(argv[3]), atoi(argv[4]) };
			getunfoldkpoints("INPOS", 0, 0, kpt, const_cast<char*>("kpt"), atof(argv[5]), argv[6][0], "redef.in");
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

	if (!strcmp("--del", argv[1]))
	{
		FILE* fp1 = fopen(argv[2], "r");
		POSCAR pos;
		readposcar(fp1, pos);
		direct_to_carts(pos.iflg, pos.vec, pos.xyz, pos.nant);
		fclose(fp1);
		double h = atof(argv[4]);
		int flag = 1;
		while (flag && pos.nant[0] != 0)
		{
			flag = 0;
			for (int i = 0; i < pos.nant[0]; i++)
			{
				if (pos.xyz[i][2] < h)
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
	}
	if (!strcmp("--en", argv[1]))
	{
		double energy = get_energy();
		FILE* fp = fopen("OUTCAR_ENERGY", "w");
		fprintf(fp,"%lf\n",energy);
		fclose(fp);
	}
	if (!strcmp("--enth", argv[1]))
	{
		vector<double> en;
		for (int i = 2; i < argc; i++)
			en.push_back(atof(argv[i]));
		double enth = get_Enthalpy_formation(en);
		printf("%lf\n", enth);
	}
	if (!strcmp("--vol", argv[1]))
	{
		char file[1024] = "NULL";
		if (argc == 2)
			strcpy(file, "INPOS");
		else
			strcpy(file, argv[2]);
		FILE* fp = fopen(file, "r");
		POSCAR pos;
		readposcar(fp, pos);
		fclose(fp);
		double vol = volume(pos.vec);
		printf("%lf\n", vol);
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
			printf("%s is not exist!\n", file);
			return -1;
		}
		POSCAR pos;
		readposcar(fp, pos);
		fclose(fp);
		FILE* fp1 = fopen("ATOM_NUM", "w");
		fprintf(fp1,"%d\n", pos.nant[0]);
		fclose(fp1);
		return 0;
	}
	if (!strcmp("--hull", argv[1]))
	{
		hull_point hp;
		if (argc == 2)
			hp.spa_convexhull();
		else
			hp.spa_convexhull(argv[2]);
		return 0;
	}
}