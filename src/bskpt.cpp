#include"../include/bskpt.h"
#include"../include/tools.h"
static SpglibError spglib_error_code = SPGLIB_SUCCESS;

static map<int, bool> time_reversal_Data = {
   {1, 0},
   {2, 1},
   {3, 0},
   {4, 0},
   {5, 0},
   {6, 0},
   {7, 0},
   {8, 0},
   {9, 0},
   {10, 1},
   {11, 1},
   {12, 1},
   {13, 1},
   {14, 1},
   {15, 1},
   {16, 0},
   {17, 0},
   {18, 0},
   {19, 0},
   {20, 0},
   {21, 0},
   {22, 0},
   {23, 0},
   {24, 0},
   {25, 0},
   {26, 0},
   {27, 0},
   {28, 0},
   {29, 0},
   {30, 0},
   {31, 0},
   {32, 0},
   {33, 0},
   {34, 0},
   {35, 0},
   {36, 0},
   {37, 0},
   {38, 0},
   {39, 0},
   {40, 0},
   {41, 0},
   {42, 0},
   {43, 0},
   {44, 0},
   {45, 0},
   {46, 0},
   {47, 1},
   {48, 1},
   {49, 1},
   {50, 1},
   {51, 1},
   {52, 1},
   {53, 1},
   {54, 1},
   {55, 1},
   {56, 1},
   {57, 1},
   {58, 1},
   {59, 1},
   {60, 1},
   {61, 1},
   {62, 1},
   {63, 1},
   {64, 1},
   {65, 1},
   {66, 1},
   {67, 1},
   {68, 1},
   {69, 1},
   {70, 1},
   {71, 1},
   {72, 1},
   {73, 1},
   {74, 1},
   {75, 0},
   {76, 0},
   {77, 0},
   {78, 0},
   {79, 0},
   {80, 0},
   {81, 0},
   {82, 0},
   {83, 1},
   {84, 1},
   {85, 1},
   {86, 1},
   {87, 1},
   {88, 1},
   {89, 0},
   {90, 0},
   {91, 0},
   {92, 0},
   {93, 0},
   {94, 0},
   {95, 0},
   {96, 0},
   {97, 0},
   {98, 0},
   {99, 0},
   {100, 0},
   {101, 0},
   {102, 0},
   {103, 0},
   {104, 0},
   {105, 0},
   {106, 0},
   {107, 0},
   {108, 0},
   {109, 0},
   {110, 0},
   {111, 0},
   {112, 0},
   {113, 0},
   {114, 0},
   {115, 0},
   {116, 0},
   {117, 0},
   {118, 0},
   {119, 0},
   {120, 0},
   {121, 0},
   {122, 0},
   {123, 1},
   {124, 1},
   {125, 1},
   {126, 1},
   {127, 1},
   {128, 1},
   {129, 1},
   {130, 1},
   {131, 1},
   {132, 1},
   {133, 1},
   {134, 1},
   {135, 1},
   {136, 1},
   {137, 1},
   {138, 1},
   {139, 1},
   {140, 1},
   {141, 1},
   {142, 1},
   {143, 0},
   {144, 0},
   {145, 0},
   {146, 0},
   {147, 1},
   {148, 1},
   {149, 0},
   {150, 0},
   {151, 0},
   {152, 0},
   {153, 0},
   {154, 0},
   {155, 0},
   {156, 0},
   {157, 0},
   {158, 0},
   {159, 0},
   {160, 0},
   {161, 0},
   {162, 1},
   {163, 1},
   {164, 1},
   {165, 1},
   {166, 1},
   {167, 1},
   {168, 0},
   {169, 0},
   {170, 0},
   {171, 0},
   {172, 0},
   {173, 0},
   {174, 0},
   {175, 1},
   {176, 1},
   {177, 0},
   {178, 0},
   {179, 0},
   {180, 0},
   {181, 0},
   {182, 0},
   {183, 0},
   {184, 0},
   {185, 0},
   {186, 0},
   {187, 0},
   {188, 0},
   {189, 0},
   {190, 0},
   {191, 1},
   {192, 1},
   {193, 1},
   {194, 1},
   {195, 0},
   {196, 0},
   {197, 0},
   {198, 0},
   {199, 0},
   {200, 1},
   {201, 1},
   {202, 1},
   {203, 1},
   {204, 1},
   {205, 1},
   {206, 1},
   {207, 0},
   {208, 0},
   {209, 0},
   {210, 0},
   {211, 0},
   {212, 0},
   {213, 0},
   {214, 0},
   {215, 0},
   {216, 0},
   {217, 0},
   {218, 0},
   {219, 0},
   {220, 0},
   {221, 1},
   {222, 1},
   {223, 1},
   {224, 1},
   {225, 1},
   {226, 1},
   {227, 1},
   {228, 1},
   {229, 1},
   {230, 1},
};

void get_cell_params(double vec[3][3], double& a, double& b, double& c, double& alpha, double& beta, double& gamma)
{
	a = sqrt(vec[0][0] * vec[0][0] + vec[0][1] * vec[0][1] + vec[0][2] * vec[0][2]);
	b = sqrt(vec[1][0] * vec[1][0] + vec[1][1] * vec[1][1] + vec[1][2] * vec[1][2]);
	c = sqrt(vec[2][0] * vec[2][0] + vec[2][1] * vec[2][1] + vec[2][2] * vec[2][2]);
	alpha = acos((vec[2][0] * vec[1][0] + vec[2][1] * vec[1][1] + vec[2][2] * vec[1][2]) / b / c);
	beta = acos((vec[0][0] * vec[2][0] + vec[0][1] * vec[2][1] + vec[0][2] * vec[2][2]) / a / c);
	gamma = acos((vec[0][0] * vec[1][0] + vec[0][1] * vec[1][1] + vec[0][2] * vec[1][2]) / a / b);
}

Bravais find_bravais(Centering centering, int spacegroup_number)
{
	Bravais bravais;
	if (spacegroup_number == 1 || spacegroup_number == 2)
	{
		bravais = TRL;
	}
	if ((spacegroup_number >= 3 && spacegroup_number < 16) && ((centering == A_FACE) || (centering == B_FACE) || (centering == C_FACE)))
	{
		bravais = MCLC;
	}
	if ((spacegroup_number >= 3 && spacegroup_number < 16) && (centering == PRIMITIVE))
	{
		bravais = MCL;
	}
	if (spacegroup_number >= 143 && spacegroup_number < 195)
	{
		vector<int> hR = { 146,148,150,161,166,167 };
		if (find(hR.begin(), hR.end(), spacegroup_number) != hR.end())
			bravais = RHL;
		else
			bravais = HEX;
	}
	if (spacegroup_number >= 16 && spacegroup_number < 75 && ((centering == A_FACE) || (centering == B_FACE) || (centering == C_FACE)))
	{
		bravais = ORCC;
	}
	if (spacegroup_number >= 16 && spacegroup_number < 75 && (centering == BODY))
	{
		bravais = ORCI;
	}
	if (spacegroup_number >= 16 && spacegroup_number < 75 && (centering == FACE))
	{
		bravais = ORCF;
	}
	if (spacegroup_number >= 16 && spacegroup_number < 75 && (centering == PRIMITIVE))
	{
		bravais = ORC;
	}
	if (spacegroup_number >= 75 && spacegroup_number < 143 && (centering == BODY))
	{
		bravais = BCT;
	}
	if (spacegroup_number >= 75 && spacegroup_number < 143 && (centering == PRIMITIVE))
	{
		bravais = TET;
	}
	if (spacegroup_number >= 195 && spacegroup_number < 230 && (centering == BODY))
	{
		bravais = BCC;
	}
	if (spacegroup_number >= 195 && spacegroup_number < 230 && (centering == FACE))
	{
		bravais = FCC;
	}
	if (spacegroup_number >= 195 && spacegroup_number < 230 && (centering == PRIMITIVE))
	{
		bravais = CUB;
	}
	return bravais;
}

vector<Symmpoint> get_path(string path, map<char, Symmpoint> path_point)
{
	vector<Symmpoint>highsymmpoint;
	for (int i = 0; i < 2 * path.length() - 2; i++)
	{
		if (path[i / 2] == '|' || path[i / 2 + 1] == '|')
			continue;
		if (i % 2 == 0)
		{
			highsymmpoint.push_back(path_point[path[i / 2]]);
		}
		if (i % 2 == 1)
		{
			highsymmpoint.push_back(path_point[path[i / 2 + 1]]);
		}
	}
	return highsymmpoint;
}

vector<Symmpoint> find_highsymmpoint(Bravais bravais, double vec[3][3], int spacegroup_num, bool time_reversal)
{
	vector<Symmpoint> point;
	map<char, Symmpoint> path_point;
	string path;
	string s;
	// a,b,c,alpha,beta,gamma : length and angles of convention lattice vectors
	double a, b, c, alpha, beta, gamma;
	get_cell_params(vec, a, b, c, alpha, beta, gamma);
	if (bravais == CUB)//Q->X1
	{
		s = "GRMXQ";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, 0.5, 0.5, "R");
		Symmpoint p3(0.5, 0.5, 0, "M");
		Symmpoint p4(0, 0.5, 0, "X");
		Symmpoint p5(0.5, 0, 0, "X_1");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		vector<int> space = { 195,198,200,201,205 };
		vector<int>::iterator it = find(space.begin(), space.end(), spacegroup_num);
		if (it != space.end())
			path = "GXMGRX|RMQ";
		else
			path = "GXMGRX|RM";
	}
	if (bravais == FCC)// Q->W2
	{
		s = "GXLWQKU";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, 0, 0.5, "X");
		Symmpoint p3(0.5, 0.5, 0.5, "L");
		Symmpoint p4(0.5, 0.25, 0.75, "W");
		Symmpoint p5(0.75, 0.25, 0.5, "W_2");
		Symmpoint p6(0.375, 0.375, 0.75, "K");
		Symmpoint p7(0.625, 0.25, 0.625, "U");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		point.push_back(p6);
		point.push_back(p7);
		vector<int> space = { 196,202,204 };
		vector<int>::iterator it = find(space.begin(), space.end(), spacegroup_num);
		if (it != space.end())
			path = "GXU|KGLWXQ";
		else
			path = "GXU|KGLWX";
	}
	if (bravais == BCC)
	{
		s = "GHPN";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, -0.5, 0.5, "H");
		Symmpoint p3(0.25, 0.25, 0.25, "P");
		Symmpoint p4(0, 0, 0.5, "N");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		path = "GHNGPH|PN";
	}
	if (bravais == TET)
	{
		s = "GZMARX";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0, 0, 0.5, "Z");
		Symmpoint p3(0.5, 0.5, 0, "M");
		Symmpoint p4(0.5, 0.5, 0.5, "A");
		Symmpoint p5(0, 0.5, 0.5, "R");
		Symmpoint p6(0, 0.5, 0, "X");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		point.push_back(p6);
		path = "GXMGZRAZ|XR|MA";
	}
	if (bravais == BCT)
	{
		if (c < a)//BCT1 Q->Z0
		{
			double x1 = (1 + c * c / a / a) / 4;
			s = "GMXPZQN";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(-0.5, 0.5, 0.5, "M");
			Symmpoint p3(0, 0, 0.5, "X");
			Symmpoint p4(0.25, 0.25, 0.25, "P");
			Symmpoint p5(x1, x1, -x1, "Z");
			Symmpoint p6(-x1, 1 - x1, x1, "Z_0");
			Symmpoint p7(0, 0.5, 0, "N");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			path = "GXMGZ|QM|XPNG";
		}
		if (c >= a)//BCT2 Q->S0 W->G
		{
			double x1 = (1 + a * a / c / c) / 4;
			double x2 = a * a / (2 * c * c);
			s = "GMXPNQSRW";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, 0.5, -0.5, "M");
			Symmpoint p3(0, 0, 0.5, "X");
			Symmpoint p4(0.25, 0.25, 0.25, "P");
			Symmpoint p5(0, 0.5, 0, "N");
			Symmpoint p6(-x1, x1, x1, "S_0");
			Symmpoint p7(x1, 1 - x1, -x1, "S");
			Symmpoint p8(-x2, x2, 0.5, "R");
			Symmpoint p9(0.5, 0.5, -x2, "G");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			path = "GXPNGMS|QG|XR|WM";
		}
	}
	if (bravais == ORC)
	{
		s = "GXZUYSTR";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, 0, 0, "X");
		Symmpoint p3(0, 0, 0.5, "Z");
		Symmpoint p4(0.5, 0, 0.5, "U");
		Symmpoint p5(0, 0.5, 0, "Y");
		Symmpoint p6(0.5, 0.5, 0, "S");
		Symmpoint p7(0, 0.5, 0.5, "T");
		Symmpoint p8(0.5, 0.5, 0.5, "R");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		point.push_back(p6);
		point.push_back(p7);
		point.push_back(p8);
		path = "GXSYGZURTZ|XU|YT|SR";
	}
	if (bravais == ORCF)
	{
		if ((1 / a / a) > (1 / b / b + 1 / c / c))//ORCF1 S->∑0 U->U0 A->A0 C->C0
		{
			double x1 = (1 + a * a / b / b - a * a / c / c) / 4;
			double x2 = (1 + a * a / b / b + a * a / c / c) / 4;
			s = "GTZYSUACL";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(1, 0.5, 0.5, "T");
			Symmpoint p3(0.5, 0.5, 0, "Z");
			Symmpoint p4(0.5, 0, 0.5, "Y");
			Symmpoint p5(0, x2, x2, "SIGMA_0");
			Symmpoint p6(1, 1 - x2, 1 - x2, "U_0");
			Symmpoint p7(0.5, 0.5 + x1, x1, "A_0");
			Symmpoint p8(0.5, 0.5 - x1, 1 - x1, "C_0");
			Symmpoint p9(0.5, 0.5, 0.5, "L");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			path = "GYTZGS|UT|YC|AZ|GL";
		}
		else if ((1 / c / c) > (1 / b / b + 1 / a / a))//ORCF2 D->Ʌ0 Q->Q0 W->G0 H->H0
		{
			double x1 = (1 + c * c / a / a - c * c / b / b) / 4;
			double x2 = (1 + c * c / a / a + c * c / b / b) / 4;
			s = "GTZYDQWHL";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0, 0.5, 0.5, "T");
			Symmpoint p3(0.5, 0.5, 1, "Z");
			Symmpoint p4(0.5, 0, 0.5, "Y");
			Symmpoint p5(x2, x2, 0, "LAMBDA_0");
			Symmpoint p6(1 - x2, 1 - x2, 1, "Q_0");
			Symmpoint p7(0.5 - x1, 1 - x1, 0.5, "G_0");
			Symmpoint p8(0.5 + x1, x1, 0.5, "H_0");
			Symmpoint p9(0.5, 0.5, 0.5, "L");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			path = "GTZYGD|QZ|TW|HY|GL";
		}
		else//ORCF3 A->A0 C->C0 B->B0 D->D0 W->G0 H->H0
		{
			double x1 = (1 + a * a / b / b - a * a / c / c) / 4;
			double x2 = (1 + b * b / a / a - b * b / c / c) / 4;
			double x3 = (1 + c * c / b / b - c * c / a / a) / 4;
			s = "GTZYACBDQHL";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0, 0.5, 0.5, "T");
			Symmpoint p3(0.5, 0.5, 0, "Z");
			Symmpoint p4(0.5, 0, 0.5, "Y");
			Symmpoint p5(0.5, 0.5 + x1, x1, "A_0");
			Symmpoint p6(0.5, 0.5 - x1, 1 - x1, "C_0");
			Symmpoint p7(0.5 + x2, 0.5, x2, "B_0");
			Symmpoint p8(0.5 - x2, 0.5, 1 - x2, "D_0");
			Symmpoint p9(x3, 0.5 + x3, 0.5, "G_0");
			Symmpoint p10(1 - x3, 0.5 - x3, 0.5, "H_0");
			Symmpoint p11(0.5, 0.5, 0.5, "L");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			path = "GYC|AZB|DTQ|HY|TGZ|GL";
		}
	}
	if (bravais == ORCI)
	{
		if (c >= a && c >= b)//ORCI1 Q->∑0 F->F2 Y->Y0 U->U0 L->L0 M->M0 J->J0
		{
			double x1 = (1 + a * a / c / c) / 4;
			double x2 = (1 + b * b / c / c) / 4;
			double x3 = (b * b - a * a) / c / c / 4;
			double x4 = (a * a + b * b) / c / c / 4;
			s = "GXSRTWQFYULMJ";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, 0.5, -0.5, "X");
			Symmpoint p3(0.5, 0, 0, "S");
			Symmpoint p4(0, 0.5, 0, "R");
			Symmpoint p5(0, 0, 0.5, "T");
			Symmpoint p6(0.25, 0.25, 0.25, "W");
			Symmpoint p7(-x1, x1, x1, "SIGMA_0");
			Symmpoint p8(x1, 1 - x1, -x1, "F_2");
			Symmpoint p9(x2, -x2, x2, "Y_0");
			Symmpoint p10(1 - x2, x2, -x2, "U_0");
			Symmpoint p11(-x4, x4, 0.5 - x3, "L_0");
			Symmpoint p12(x4, -x4, 0.5 + x3, "M_0");
			Symmpoint p13(0.5 - x3, 0.5 + x3, -x4, "J_0");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			path = "GXF|QGY|UX|GRWSGTW";
		}
		if (a >= c && a >= b)//ORCI2 Y->Y0 U->U2 L->Ʌ0 Q->G2 E->K2 Y->K4
		{
			double x1 = (1 + b * b / a / a) / 4;
			double x2 = (1 + c * c / a / a) / 4;
			double x3 = (c * c - b * b) / a / a / 4;
			double x4 = (c * c + b * b) / a / a / 4;
			s = "GXSRTWYULQKEY";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(-0.5, 0.5, 0.5, "X");
			Symmpoint p3(0.5, 0, 0, "S");
			Symmpoint p4(0, 0.5, 0, "R");
			Symmpoint p5(0, 0, 0.5, "T");
			Symmpoint p6(0.25, 0.25, 0.25, "W");
			Symmpoint p7(x1, -x1, x1, "Y_0");
			Symmpoint p8(-x1, x1, 1 - x1, "U_2");
			Symmpoint p9(x2, x2, -x2, "LAMBDA_0");
			Symmpoint p10(-x2, 1 - x2, x2, "G_2");
			Symmpoint p11(0.5 - x3, -x4, x4, "K");
			Symmpoint p12(0.5 + x3, x4, -x4, "K_2");
			Symmpoint p13(-x4, 0.5 - x3, 0.5 + x3, "K_4");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			path = "GXU|YGL|QX|GRWSGTW";
		}
		if (b >= c && b >= a)//ORCI3 M->∑0 F->F0 L->Ʌ0 Q->G0 V->V0 H->H0 E->H2
		{
			double x1 = (1 + c * c / b / b) / 4;
			double x2 = (1 + a * a / b / b) / 4;
			double x3 = (a * a - c * c) / b / b / 4;
			double x4 = (c * c + a * a) / b / b / 4;
			s = "GXSRTWMFLQVHE";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, -0.5, 0.5, "X");
			Symmpoint p3(0.5, 0, 0, "S");
			Symmpoint p4(0, 0.5, 0, "R");
			Symmpoint p5(0, 0, 0.5, "T");
			Symmpoint p6(0.25, 0.25, 0.25, "W");
			Symmpoint p7(-x2, x2, x2, "SIGMA_0");
			Symmpoint p8(x2, -x2, 1 - x2, "F_0");
			Symmpoint p9(x1, x1, -x1, "LAMBDA_0");
			Symmpoint p10(1 - x1, -x1, x1, "G_0");
			Symmpoint p11(x4, 0.5 - x3, -x4, "V_0");
			Symmpoint p12(-x4, 0.5 + x3, x4, "H_0");
			Symmpoint p13(0.5 + x3, -x4, 0.5 - x3, "H_2");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			path = "GXF|MGL|QX|GRWSGTW";
		}
	}
	if (bravais == ORCC)
	{
		vector<int> oA = { 38,39,40,41 };
		//oC
		if ((find(oA.begin(), oA.end(), spacegroup_num) == oA.end() && a < b))//M->∑0 C->C0 A->A0 E->E0
		{
			double x1 = (1 + a * a / b / b) / 4;
			s = "GYTZSRMCAE";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(-0.5, 0.5, 0, "Y");
			Symmpoint p3(-0.5, 0.5, 0.5, "T");
			Symmpoint p4(0, 0, 0.5, "Z");
			Symmpoint p5(0, 0.5, 0, "S");
			Symmpoint p6(0, 0.5, 0.5, "R");
			Symmpoint p7(x1, x1, 0, "SIGMA_0");
			Symmpoint p8(-x1, 1 - x1, 0, "C_0");
			Symmpoint p9(x1, x1, 0.5, "A_0");
			Symmpoint p10(-x1, 1 - x1, 0.5, "E_0");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			path = "GYC|MGZA|ETY|GSRZT";
		}
		if ((find(oA.begin(), oA.end(), spacegroup_num) == oA.end() && a >= b))//Q->T2 W->Z2 E->R2 D->Δ0 U->G0 I->G2
		{
			double x1 = (1 + b * b / a / a) / 4;
			s = "GYTQZWSREDFBUI";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, 0.5, 0, "Y");
			Symmpoint p3(0.5, 0.5, 0.5, "T");
			Symmpoint p4(0.5, 0.5, -0.5, "T_2");
			Symmpoint p5(0, 0, 0.5, "Z");
			Symmpoint p6(0, 0, -0.5, "Z_2");
			Symmpoint p7(0, 0.5, 0, "S");
			Symmpoint p8(0, 0.5, 0.5, "R");
			Symmpoint p9(0, 0.5, -0.5, "R_2");
			Symmpoint p10(-x1, x1, 0, "Delta_0");
			Symmpoint p11(x1, 1 - x1, 0, "F_0");
			Symmpoint p12(-x1, x1, 0.5, "B_0");
			Symmpoint p13(-x1, x1, -0.5, "B_2");
			Symmpoint p14(x1, 1 - x1, 0.5, "G_0");
			Symmpoint p15(x1, 1 - x1, -0.5, "G_2");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			point.push_back(p14);
			point.push_back(p15);
			path = "GYF|DGZB|UTY|GSRZT";
		}
		//oA
		if (find(oA.begin(), oA.end(), spacegroup_num) != oA.end() && b < c)
		{
			double x1 = (1 + b * b / c / c) / 4;
			s = "GYTZSRMCAE";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(-0.5, 0.5, 0, "Y");
			Symmpoint p3(-0.5, 0.5, 0.5, "T");
			Symmpoint p4(0, 0, 0.5, "Z");
			Symmpoint p5(0, 0.5, 0, "S");
			Symmpoint p6(0, 0.5, 0.5, "R");
			Symmpoint p7(x1, x1, 0, "SIGMA_0");
			Symmpoint p8(-x1, 1 - x1, 0, "C_0");
			Symmpoint p9(x1, x1, 0.5, "A_0");
			Symmpoint p10(-x1, 1 - x1, 0.5, "E_0");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			path = "GYC|MGZA|ETY|GSRZT";
		}
		if (find(oA.begin(), oA.end(), spacegroup_num) != oA.end() && b >= c)
		{
			double x1 = (1 + c * c / b / b) / 4;
			s = "GYTQZWSREDFBUI";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, 0.5, 0, "Y");
			Symmpoint p3(0.5, 0.5, 0.5, "T");
			Symmpoint p4(0.5, 0.5, -0.5, "T_2");
			Symmpoint p5(0, 0, 0.5, "Z");
			Symmpoint p6(0, 0, -0.5, "Z_2");
			Symmpoint p7(0, 0.5, 0, "S");
			Symmpoint p8(0, 0.5, 0.5, "R");
			Symmpoint p9(0, 0.5, -0.5, "R_2");
			Symmpoint p10(-x1, x1, 0, "Delta_0");
			Symmpoint p11(x1, 1 - x1, 0, "F_0");
			Symmpoint p12(-x1, x1, 0.5, "B_0");
			Symmpoint p13(-x1, x1, -0.5, "B_2");
			Symmpoint p14(x1, 1 - x1, 0.5, "G_0");
			Symmpoint p15(x1, 1 - x1, -0.5, "G_2");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			point.push_back(p14);
			point.push_back(p15);
			path = "GYF|DGZB|UTY|GSRZT";
		}
	}
	if (bravais == HEX)//Q->H2
	{
		s = "GAKHQML";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0, 0, 0.5, "A");
		Symmpoint p3(1 / 3.0, 1 / 3.0, 0, "K");
		Symmpoint p4(1 / 3.0, 1 / 3.0, 0.5, "H");
		Symmpoint p5(1 / 3.0, 1 / 3.0, -0.5, "H_2");
		Symmpoint p6(0.5, 0, 0, "M");
		Symmpoint p7(0.5, 0, 0.5, "L");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		point.push_back(p6);
		point.push_back(p7);
		vector<int> space = { 143,144,145,147,149,151,153,157,159,162,163 };
		vector<int>::iterator it = find(space.begin(), space.end(), spacegroup_num);
		if (it != space.end())
			path = "GMKGALHA|LM|HKQ";
		else
			path = "GMKGALHA|LM|HK";
	}
	if (bravais == RHL)
	{
		if (sqrt(3) * a < sqrt(2) * c)//Q->L2 W->L4 E->F2 S->S0 R->S2 Y->S4 U->S6 H->H0 I->H2 O->H4 P->H6
			//M->M0 D->M2 J->M4 K->M6 Z->M8
		{
			double x1 = a * a / c / c / 4;
			double x2 = 5.0 / 6 - 2 * x1;
			double x3 = 1.0 / 3 + x1;
			s = "GTLQWFESRYUHIOPMDJKZ";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, 0.5, 0.5, "T");
			Symmpoint p3(0.5, 0, 0, "L");
			Symmpoint p4(0, -0.5, 0, "L_2");
			Symmpoint p5(0, 0, -0.5, "L_4");
			Symmpoint p6(0.5, 0, 0.5, "F");
			Symmpoint p7(0.5, 0.5, 0, "F_2");
			Symmpoint p8(x3, -x3, 0, "S_0");
			Symmpoint p9(1 - x3, 0, x3, "S_2");
			Symmpoint p10(x3, 0, -x3, "S_4");
			Symmpoint p11(1 - x3, x3, 0, "S_6");
			Symmpoint p12(0.5, -1 + x2, 1 - x2, "H_0");
			Symmpoint p13(x2, 1 - x2, 0.5, "H_2");
			Symmpoint p14(x2, 0.5, 1 - x2, "H_4");
			Symmpoint p15(0.5, 1 - x2, -1 + x2, "H_6");
			Symmpoint p16(x3, -1 + x2, x3, "M_0");
			Symmpoint p17(1 - x3, 1 - x2, 1 - x3, "M_2");
			Symmpoint p18(x2, x3, x3, "M_4");
			Symmpoint p19(1 - x3, 1 - x3, 1 - x2, "M_6");
			Symmpoint p20(x3, x3, -1 + x2, "M_8");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			point.push_back(p14);
			point.push_back(p15);
			point.push_back(p16);
			point.push_back(p17);
			point.push_back(p18);
			point.push_back(p19);
			point.push_back(p20);
			path = "GTI|HLGS|RFG";
		}
		if (sqrt(3) * a >= sqrt(2) * c)//P->P0 Q->P2 R->R0 W->M2
		{
			double x1 = 1.0 / 6 - c * c / a / a / 9;
			double x2 = 1.0 / 2 - 2 * x1;
			double x3 = 1.0 / 2 + x1;
			s = "GTPQRMWLF";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, -0.5, 0.5, "T");
			Symmpoint p3(x2, -1 + x2, x2, "P_0");
			Symmpoint p4(x2, x2, x2, "P_2");
			Symmpoint p5(1 - x2, -x2, -x2, "R_0");
			Symmpoint p6(1 - x3, -x3, 1 - x3, "M");
			Symmpoint p7(x3, -1 + x3, -1 + x3, "M_2");
			Symmpoint p8(0.5, 0, 0, "L");
			Symmpoint p9(0.5, -0.5, 0, "F");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			path = "GLTP|QGF";
		}
	}
	if (bravais == MCL)//Q->B2 W->Y2 R->C2 T->D2 U->H2 I->H4 O->M2 P->M4
	{
		double x1 = (1 + (a / c) * cos(beta)) / sin(beta) / sin(beta) / 2;
		double x2 = 1.0 / 2 + x1 * c * cos(beta) / a;
		s = "GZBQYWCRDTAEHUIMOP";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0, 0.5, 0, "Z");
		Symmpoint p3(0, 0, 0.5, "B");
		Symmpoint p4(0, 0, -0.5, "B_2");
		Symmpoint p5(0.5, 0, 0, "Y");
		Symmpoint p6(-0.5, 0, 0, "Y_2");
		Symmpoint p7(0.5, 0.5, 0, "C");
		Symmpoint p8(-0.5, 0.5, 0, "C_2");
		Symmpoint p9(0, 0.5, 0.5, "D");
		Symmpoint p10(0, 0.5, -0.5, "D_2");
		Symmpoint p11(-0.5, 0, 0.5, "A");
		Symmpoint p12(-0.5, 0.5, 0.5, "E");
		Symmpoint p13(-x1, 0, 1 - x2, "H");
		Symmpoint p14(-1 + x1, 0, x2, "H_2");
		Symmpoint p15(-x1, 0, -x2, "H_4");
		Symmpoint p16(-x1, 0.5, 1 - x2, "M");
		Symmpoint p17(-1 + x1, 0.5, x2, "M_2");
		Symmpoint p18(-x1, 0.5, -x2, "M_4");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		point.push_back(p6);
		point.push_back(p7);
		point.push_back(p8);
		point.push_back(p9);
		point.push_back(p10);
		point.push_back(p11);
		point.push_back(p12);
		point.push_back(p13);
		point.push_back(p14);
		point.push_back(p15);
		point.push_back(p16);
		point.push_back(p17);
		point.push_back(p18);
		path = "GZDBGAEZRWG";
	}
	if (bravais == MCLC)
	{
		double para = -a * cos(beta) / c + a * a * sin(beta) * sin(beta) / b / b;
		if (b < a * sin(beta))//Y->Y2 Q->Y4 M->M2 W->V2 L->L2 E->C2 R->C4 T->D2 U->E2 I->E4
		{
			double x1 = (2 + a / c * cos(beta)) / sin(beta) / sin(beta) / 4;
			double x2 = 1.0 / 2 - 2 * x1 * c * cos(beta) / a;
			double x3 = 3.0 / 4 - b * b / a / a / sin(beta) / sin(beta) / 4;
			double x4 = x3 - (3.0 / 4 - x3) * a * cos(beta) / c;
			s = "GYQAMVWLCERDTEUI";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(-0.5, 0.5, 0, "Y_2");
			Symmpoint p3(0.5, -0.5, 0, "Y_4");
			Symmpoint p4(0, 0, 0.5, "A");
			Symmpoint p5(-0.5, 0.5, 0.5, "M_2");
			Symmpoint p6(0.5, 0, 0, "V");
			Symmpoint p7(0, 0.5, 0, "V_2");
			Symmpoint p8(0, 0.5, 0.5, "L_2");
			Symmpoint p9(1 - x3, 1 - x3, 0, "C");
			Symmpoint p10(-1 + x3, x3, 0, "C_2");
			Symmpoint p11(x3, -1 + x3, 0, "C_4");
			Symmpoint p12(-1 + x4, x4, 0.5, "D");
			Symmpoint p13(1 - x4, 1 - x4, 0.5, "D_2");
			Symmpoint p14(-1 + x1, 1 - x1, 1 - x2, "E");
			Symmpoint p15(-x1, x1, x2, "E_2");
			Symmpoint p16(x1, -x1, 1 - x2, "E_4");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			point.push_back(p14);
			point.push_back(p15);
			point.push_back(p16);
			path = "GC|EYGMD|TAG|LGW";
		}
		if (b >= a * sin(beta) && para <= 1)//V->V2 L->L2 Q->F2 W->F4 E->H2 R->H4 P->G T->G2 U->G4 I->G6
		{
			double x1 = (1 + a * a / b / b) / 4;
			double x2 = -a * c * cos(beta) / b / b / 2;
			double x3 = (a * a / b / b + (1 + a / c * cos(beta)) / sin(beta) / sin(beta)) / 4;
			double x4 = 1.0 / 2 - 2 * x3 * c * cos(beta) / a;
			double x5 = 1 + x3 - 2 * x1;
			double x6 = x4 - 2 * x2;
			s = "GYAMVLFQWHERPTUI";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, 0.5, 0, "Y");
			Symmpoint p3(0, 0, 0.5, "A");
			Symmpoint p4(0.5, 0.5, 0.5, "M");
			Symmpoint p5(0, 0.5, 0, "V_2");
			Symmpoint p6(0, 0.5, 0.5, "L_2");
			Symmpoint p7(-1 + x5, 1 - x5, 1 - x6, "F");
			Symmpoint p8(1 - x5, x5, x6, "F_2");
			Symmpoint p9(x5, 1 - x5, 1 - x6, "F_4");
			Symmpoint p10(-x3, x3, x4, "H");
			Symmpoint p11(x3, 1 - x3, 1 - x4, "H_2");
			Symmpoint p12(x3, -x3, 1 - x4, "H_4");
			Symmpoint p13(-x1, x1, x2, "G");
			Symmpoint p14(x1, 1 - x1, -x2, "G_2");
			Symmpoint p15(x1, -x1, -x2, "G_4");
			Symmpoint p16(1 - x1, x1, x2, "G_6");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			point.push_back(p14);
			point.push_back(p15);
			point.push_back(p16);
			path = "GYMAG|LGV";
		}
		if (b >= a * sin(beta) && para > 1)//M->M2 Q->V2 L->L2 W->I2 E->K2 R->K4 T->H2 U->H4 A->N2 O->N4 P->N6
		{
			double x1 = (a * a / b / b + (1 + a / c * cos(beta)) / sin(beta) / sin(beta)) / 4;
			double x2 = 1 - x1 * b * b / a / a;
			double x3 = 1.0 / 2 - 2 * x1 * c * cos(beta) / a;
			double x4 = x3 / 2 + a * a / b / b / 4 + a * c * cos(beta) / b / b / 2;
			double x5 = 2 * x4 - x1;
			double x6 = c / 2 / a / cos(beta) * (1 - 4 * x5 + a * a * sin(beta) * sin(beta) / b / b);
			double x7 = -1.0 / 4 + x6 / 2 - x1 * c * cos(beta) / a;
			s = "GYAMVQLIWKERHTUNAOP";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0.5, 0.5, 0, "Y");
			Symmpoint p3(0, 0, 0.5, "A");
			Symmpoint p4(-0.5, 0.5, 0.5, "M_2");
			Symmpoint p5(0.5, 0, 0, "V");
			Symmpoint p6(0, 0.5, 0, "V_2");
			Symmpoint p7(0, 0.5, 0.5, "L_2");
			Symmpoint p8(-1 + x2, x2, 0.5, "I");
			Symmpoint p9(1 - x2, 1 - x2, 0.5, "I_2");
			Symmpoint p10(-x5, x5, x6, "K");
			Symmpoint p11(-1 + x5, 1 - x5, 1 - x6, "K_2");
			Symmpoint p12(1 - x5, x5, x6, "K_4");
			Symmpoint p13(-x1, x1, x3, "H");
			Symmpoint p14(x1, 1 - x1, 1 - x3, "H_2");
			Symmpoint p15(x1, -x1, 1 - x3, "H_4");
			Symmpoint p16(-x4, x4, x7, "N");
			Symmpoint p17(x4, 1 - x4, -x7, "N_2");
			Symmpoint p18(x4, -x4, -x7, "N_4");
			Symmpoint p19(1 - x4, x4, x7, "N_6");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			point.push_back(p10);
			point.push_back(p11);
			point.push_back(p12);
			point.push_back(p13);
			point.push_back(p14);
			point.push_back(p15);
			point.push_back(p16);
			point.push_back(p17);
			point.push_back(p18);
			point.push_back(p19);
			path = "GAW|IMGY|LGQ";
		}
	}
	if (bravais == TRL)
	{
		if (alpha >= Pi / 2 && beta >= Pi / 2 && gamma >= Pi / 2)//reduce cell is all-obtuse
		{
			s = "GZYXVUTR";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0, 0, 0.5, "Z");
			Symmpoint p3(0, 0.5, 0, "Y");
			Symmpoint p4(0.5, 0, 0, "X");
			Symmpoint p5(0.5, 0.5, 0, "V");
			Symmpoint p6(0.5, 0, 0.5, "U");
			Symmpoint p7(0, 0.5, 0.5, "T");
			Symmpoint p8(0.5, 0.5, 0.5, "R");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			path = "GX|YGZ|RGT|UGV";
		}
		if (alpha < Pi / 2 && beta < Pi / 2 && gamma < Pi / 2)//reduce cell is all-acute Q->Y2 V->V2 U->U2 T->T2 R->R2
		{
			s = "GZYQXVUTR";
			Symmpoint p1(0, 0, 0, "GAMMA");
			Symmpoint p2(0, 0, 0.5, "Z");
			Symmpoint p3(0, 0.5, 0, "Y");
			Symmpoint p4(0, -0.5, 0, "Y_2");
			Symmpoint p5(0.5, 0, 0, "X");
			Symmpoint p6(0.5, -0.5, 0, "V_2");
			Symmpoint p7(-0.5, 0, 0.5, "U_2");
			Symmpoint p8(0, -0.5, 0.5, "T_2");
			Symmpoint p9(-0.5, -0.5, 0.5, "R_2");
			point.push_back(p1);
			point.push_back(p2);
			point.push_back(p3);
			point.push_back(p4);
			point.push_back(p5);
			point.push_back(p6);
			point.push_back(p7);
			point.push_back(p8);
			point.push_back(p9);
			path = "GX|YGZ|RGT|UGV";
		}
	}
	for (int i = 0; i < s.size(); i++)
		path_point.insert(pair <char, Symmpoint>(s[i], point[i]));
	write_highsymmpoint(point, path, path_point, time_reversal);
	return 	get_path(path, path_point);
}

void write_kpt(vector<Symmpoint> highsymmpoint, int kppra, bool time_reversal)
{
	FILE* fp = fopen("NEWKPATH", "w");
	fprintf(fp, "K-Path Generated by VASPMATE\n");
	fprintf(fp, "	%d\n", kppra);
	fprintf(fp, "Line-Mode\n");
	fprintf(fp, "Reciprocal\n");
	for (int i = 0; i < highsymmpoint.size(); i++)
	{
		fprintf(fp, "	%10f   %10f   %10f	%s\n", highsymmpoint[i].b1, highsymmpoint[i].b2, highsymmpoint[i].b3, highsymmpoint[i].symbol.c_str());
		if (i % 2 == 1)
			fprintf(fp, "\n");
	}
	string gamma = "GAMMA";
	string time_reversal_symbol;
	if (!time_reversal)
	{
		for (int i = 0; i < highsymmpoint.size(); i++)
		{
			if (highsymmpoint[i].symbol == gamma)
			{
				time_reversal_symbol = highsymmpoint[i].symbol;
				fprintf(fp, "	%10f   %10f   %10f	%s\n", highsymmpoint[i].b1, highsymmpoint[i].b2, highsymmpoint[i].b3, time_reversal_symbol.c_str());
			}
			else
			{
				time_reversal_symbol = highsymmpoint[i].symbol + "'";
				fprintf(fp, "	%10f   %10f   %10f	%s\n", -highsymmpoint[i].b1, -highsymmpoint[i].b2, -highsymmpoint[i].b3, time_reversal_symbol.c_str());
			}
			if (i % 2 == 1)
				fprintf(fp, "\n");
		}

	}
	fclose(fp);
	printf("Written NEWKPATH file!\n");
}

void write_highsymmpoint(vector<Symmpoint> point, string path, map<char, Symmpoint> path_point, bool time_reversal)
{
	FILE* fp = fopen("HIGH_SYMMETRY_POINT", "w");
	fprintf(fp, "High - symmetry points(in fractional coordinates).\n");
	if (!time_reversal)
		fprintf(fp, "Note on time reversal : The second half of the path is required only if the system does not have time - reversal symmetry\n");
	for (int i = 0; i < point.size(); i++)
	{
		fprintf(fp, "	%10f   %10f   %10f  %s\n", point[i].b1, point[i].b2, point[i].b3, point[i].symbol.c_str());
	}
	string gamma = "GAMMA";
	if (!time_reversal)
	{
		for (int i = 0; i < point.size(); i++)
		{
			if (point[i].symbol == gamma)
				continue;
			else
			{
				string time_reversal_symbol = point[i].symbol + "'";
				fprintf(fp, "	%10f   %10f   %10f  %s\n", -point[i].b1, -point[i].b2, -point[i].b3, time_reversal_symbol.c_str());
			}
		}
	}
	string suggest_path;
	for (int i = 0; i < path.size(); i++)
	{
		if (path[i] == '|')
			suggest_path += path[i];
		else
		{
			if (i < path.size() - 1 && path[i + 1] != '|')
				suggest_path += path_point[path[i]].symbol + '-';
			else
				suggest_path += path_point[path[i]].symbol;
		}
	}
	if (!time_reversal)
	{
		suggest_path += '|';
		for (int i = 0; i < path.size(); i++)
		{
			if (path[i] == '|')
				suggest_path += path[i];
			else
			{
				string time_reversal_symbol = path_point[path[i]].symbol + "'";
				if (i < path.size() - 1 && path[i + 1] != '|')
					suggest_path += time_reversal_symbol + '-';
				else
					suggest_path += time_reversal_symbol;
			}
		}
	}
	fprintf(fp, "suggest path : %s\n", suggest_path.c_str());
	fclose(fp);
	printf("Written HIGH_SYMMETRY_POINT file!\n");
}

static void set_cell(double lattice[3][3],
	double position[][3],
	int types[],
	Cell* cell)
{
	int i;

	mat_copy_matrix_d3(lattice, cell->lattice);
	for (i = 0; i < cell->size; i++) {
		types[i] = cell->types[i];
		mat_copy_vector_d3(position[i], cell->position[i]);
	}
}

Cell* Transform_to_Primitive(const Cell* cell,
	SPGCONST double trans_mat[3][3],
	const Centering centering,
	int spacegroup_number,
	const double symprec)
{
	int* mapping_table;
	double tmat[3][3], tmat_inv[3][3], prim_lat[3][3];
	Cell* primitive;

	mapping_table = NULL;
	primitive = NULL;

	mat_inverse_matrix_d3(tmat_inv, trans_mat, 0);

	static double a_mat[3][3] = { {    1,    0,    0},
							 {    0, 1. / 2,-1. / 2},
							 {    0, 1. / 2, 1. / 2} };
	static double c_mat[3][3] = { { 1. / 2, 1. / 2,    0},
								 {-1. / 2, 1. / 2,    0},
								 {    0,    0,    1} };
	static double r_mat[3][3] = { { 2. / 3,-1. / 3,-1. / 3 },
								 { 1. / 3, 1. / 3,-2. / 3 },
								 { 1. / 3, 1. / 3, 1. / 3 } };
	static double i_mat[3][3] = { {-1. / 2, 1. / 2, 1. / 2 },
								 { 1. / 2,-1. / 2, 1. / 2 },
								 { 1. / 2, 1. / 2,-1. / 2 } };
	static double f_mat[3][3] = { {    0, 1. / 2, 1. / 2 },
								 { 1. / 2,    0, 1. / 2 },
								 { 1. / 2, 1. / 2,    0 } };
	static double mC_mat[3][3] = { 0.5, -0.5, 0, 0.5, 0.5, 0, 0, 0, 1 };
	static double oA_mat[3][3] = { 0, 0, 1, 0.5, 0.5, 0, -0.5, 0.5, 0 };

	if (spacegroup_number == 38 || spacegroup_number == 39 || spacegroup_number == 40 || spacegroup_number == 41)//oA
		mat_multiply_matrix_d3(tmat, tmat_inv, oA_mat);
	else {
		switch (centering) {
		case PRIMITIVE:
			mat_copy_matrix_d3(tmat, tmat_inv);
			break;
		case A_FACE:
			if ((spacegroup_number >= 3 && spacegroup_number < 16))//mC
				mat_multiply_matrix_d3(tmat, tmat_inv, mC_mat);
			else
				mat_multiply_matrix_d3(tmat, tmat_inv, a_mat);
			break;
		case C_FACE:
			if ((spacegroup_number >= 3 && spacegroup_number < 16))//mC
				mat_multiply_matrix_d3(tmat, tmat_inv, mC_mat);
			else
				mat_multiply_matrix_d3(tmat, tmat_inv, c_mat);
			break;
		case FACE:
			mat_multiply_matrix_d3(tmat, tmat_inv, f_mat);
			break;
		case BODY:
			mat_multiply_matrix_d3(tmat, tmat_inv, i_mat);
			break;
		case R_CENTER:
			mat_multiply_matrix_d3(tmat, tmat_inv, r_mat);
			break;
		default:
			goto err;
		}
	}

	if ((mapping_table = (int*)malloc(sizeof(int) * cell->size)) == NULL) {
		warning_print("spglib: Memory could not be allocated ");
		goto err;
	}

	mat_multiply_matrix_d3(prim_lat, cell->lattice, tmat);
	primitive = cel_trim_cell(mapping_table, prim_lat, cell, symprec);

	free(mapping_table);
	mapping_table = NULL;

	return primitive;

err:
	return NULL;
}

static int Find_Primitive_Cell(double lattice[3][3],
	double position[][3],
	int types[],
	const int num_atom,
	const double symprec,
	double angle_tolerance)
{
	int num_prim_atom;
	Centering centering;
	SpglibDataset* dataset;
	Cell* primitive, * bravais;

	double identity[3][3] = { { 1, 0, 0 },
							 { 0, 1, 0 },
							 { 0, 0, 1 } };

	num_prim_atom = 0;
	dataset = NULL;
	primitive = NULL;
	bravais = NULL;

	if ((dataset = spgat_get_dataset_with_hall_number(lattice,
		position,
		types,
		num_atom,
		0,
		symprec,
		angle_tolerance)) == NULL) {
		return 0;
	}
	int spacegroup_number = dataset->spacegroup_number;

	SpacegroupType spgtype;
	spgtype = spgdb_get_spacegroup_type(dataset->hall_number);

	if ((centering = spgtype.centering) == CENTERING_ERROR) {
		goto err;
	}

	if ((bravais = cel_alloc_cell(dataset->n_std_atoms)) == NULL) {
		spg_free_dataset(dataset);
		return 0;
	}

	cel_set_cell(bravais,
		lattice,
		position,
		types);

	spg_free_dataset(dataset);

	primitive = Transform_to_Primitive(bravais,
		identity,
		centering,
		spacegroup_number,
		symprec);
	cel_free_cell(bravais);
	bravais = NULL;

	if (primitive == NULL) {
		goto err;
	}

	set_cell(lattice, position, types, primitive);
	num_prim_atom = primitive->size;

	cel_free_cell(primitive);
	primitive = NULL;

	return num_prim_atom;

err:
	spglib_error_code = SPGERR_CELL_STANDARDIZATION_FAILED;
	return 0;
}

int* ir_reciprocal_mesh(const char* file, int ikppra, double kspac, char ikmesh,
	int kpt[3], char label[], vector<int>& weight, vector<vector<double> >& irr_coor)
{
	int* mesh = (int*)malloc(sizeof(int) * 3);
	POSCAR pos;
	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return mesh;
	}
	readposcar(fp, pos);//vec : lattice of convention cell
	fclose(fp);
	_matrix kvec;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			kvec.mat[i][j] = pos.vec[i][j];
	_matrix rec = RecMat(kvec);
	bool has_vacuum_slab[3];
	has_vacuum_slab_or_not(pos, has_vacuum_slab);
	if (!strcmp(label, "ka"))
		mesh = kpta(rec, ikppra, pos.nant[0], ikmesh, has_vacuum_slab);
	if (!strcmp(label, "kv"))
		mesh = kpvf(rec, kspac, pos.nant[0], ikmesh, has_vacuum_slab);
	if (!strcmp(label, "kpt"))
	{
		mesh[0] = kpt[0];
		mesh[1] = kpt[1];
		mesh[2] = kpt[2];
	}
	int grid_address[mesh[0] * mesh[1] * mesh[2]][3];
	int grid_mapping_table[mesh[0] * mesh[1] * mesh[2]];
	int is_shift[3];
	if (ikmesh == 'G')
	{
		is_shift[0] = 0;
		is_shift[1] = 0;
		is_shift[2] = 0;
	}
	else
	{
		is_shift[0] = 1;
		is_shift[1] = 1;
		is_shift[2] = 1;
	}
	int* types = (int*)malloc(sizeof(int) * pos.nant[0]);
	translate_typenum_type(types, pos.nant[1], pos.typenum);
	double n_vec[3][3];
	transpose_matrix(pos.vec, n_vec);
	SpglibDataset* dataset;
	dataset = spg_get_dataset(n_vec, pos.xyz, types, pos.nant[0], SYMPREAC);
	int spacegroup_number = dataset->spacegroup_number;
	int num_ir = spg_get_ir_reciprocal_mesh(grid_address,
		grid_mapping_table,
		mesh,
		is_shift,
		1,
		n_vec,
		pos.xyz,
		types,
		pos.nant[0],
		SYMPREAC);
	weight.resize(num_ir);
	vector<int> nonrepet;
	for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
	{
		if (grid_mapping_table[i] == i)
			nonrepet.push_back(i);
	}
	for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
	{
		for (int j = 0; j < num_ir; j++)
		{
			if (grid_mapping_table[i] == nonrepet[j])
				weight[j]++;
		}
	}
	for (int i = 0; i < mesh[0] * mesh[1] * mesh[2]; i++)
	{
		if (grid_mapping_table[i] == i)
		{
			vector<double> temp(3);
			for (int j = 0; j < 3; j++)
			{
				temp[j] = (double)grid_address[i][j] / mesh[j];
				if (temp[j] > 0.5)
					temp[j] -= 1;
			}
			irr_coor.push_back(temp);
		}
	}
	delete[] types;
	return mesh;
}

void HSE_mesh(const char file[], int ikppra, double kspac, int kpt[], char label[], double kpra, char ikmesh)
{
	vector<int> weight;
	vector<vector<double> > irr_coor;
	int* mesh = ir_reciprocal_mesh(file, ikppra, kspac, ikmesh, kpt, label, weight, irr_coor);
	char   title[MAX_NCHAR];        // Title of system
	double latt;                    // Scaling factor
	int    ifix;                    // Select 0 or normal 1
	int    iflg;                    // Direct 0 or Cartes 1
	int    nant[2];                 // number of atoms and types
	int    typenum[MAX_NELEM];      // number of atoms of each type
	int    elemnum[MAX_NELEM];      // atomic number
	char   elemsym[MAX_NELEM][3];   // atomic symbol
	double vec[3][3];               // lattice vector
	double xyz[MAX_NATOM][3];       // atomic coordinate
	char   fix[MAX_NATOM][3];       // Selective fix on each atom
	FILE* fp1 = fopen(file, "r");
	if (fp1 == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);//vec : lattice of convention cell
	fclose(fp1);
	FILE* fp2 = fopen("NEWKPATH", "r");
	if (fp2 == NULL)
	{
		printf("NEWKPATH IS NOT EXIST!\n");
		return;
	}
	char buf[1024];
	int line = 0;
	vector<vector<double> > high_symmetry;
	while (fgets(buf, 1024, fp2) != NULL)
	{
		line++;
		if (strspn(buf, " \t\n") == strlen(buf))
			continue;
		if (line > 4)
		{
			vector<double> hs(3, 0);
			sscanf(buf, "%lf%lf%lf", &hs[0], &hs[1], &hs[2]);
			high_symmetry.push_back(hs);
		}
	}
	fclose(fp2);
	double rec_vec[3][3];
	RecMat(vec, rec_vec);
	vector<vector<double> > carts_kpt(high_symmetry.size(), vector<double>(3));
	for (int i = 0; i < high_symmetry.size(); i++)
	{
		double p[3] = { high_symmetry[i][0],high_symmetry[i][1],high_symmetry[i][2] };
		for (int j = 0; j < 3; j++)
		{
			carts_kpt[i][j] = p[0] * rec_vec[0][j] + p[1] * rec_vec[1][j] + p[2] * rec_vec[2][j];
		}
	}
	vector<int> point(high_symmetry.size() / 2);
	for (int i = 0; i < point.size(); i++)
	{
		double length = 0;
		for (int j = 0; j < 3; j++)
			length += (carts_kpt[2 * i + 1][j] - carts_kpt[2 * i][j]) * (carts_kpt[2 * i + 1][j] - carts_kpt[2 * i][j]);
		point[i] = max(5, (int)(sqrt(length) / kpra));
	}
	vector<vector<double> > band_kpt;
	for (int i = 0; i < point.size(); i++)
	{
		for (int j = 0; j < point[i]; j++)
		{
			double x = high_symmetry[2 * i][0] + (high_symmetry[2 * i + 1][0] - high_symmetry[2 * i][0]) / (point[i] - 1) * j;
			double y = high_symmetry[2 * i][1] + (high_symmetry[2 * i + 1][1] - high_symmetry[2 * i][1]) / (point[i] - 1) * j;
			double z = high_symmetry[2 * i][2] + (high_symmetry[2 * i + 1][2] - high_symmetry[2 * i][2]) / (point[i] - 1) * j;
			vector<double> temp = { x,y,z };
			band_kpt.push_back(temp);
		}
	}
	FILE* fp = fopen("NEWKPT", "w");
	//mesh0 mesh1 mesh2 num_scf num_hse num_path num_path0 ...
	fprintf(fp, "%d	%d	%d	%d", mesh[0], mesh[1], mesh[2], irr_coor.size());
	fprintf(fp, "	%d	%d", band_kpt.size(), point.size());
	for (int i = 0; i < point.size(); i++)
		fprintf(fp, "	%d", point[i]);
	fprintf(fp, "          # Parameters to Generate KPOINTS (Do NOT Edit This Line)\n");
	fprintf(fp, "	%d\n", irr_coor.size() + band_kpt.size());
	fprintf(fp, "Reciprocal lattice\n");
	for (int i = 0; i < irr_coor.size(); i++)
		fprintf(fp, "    %lf	%lf	%lf	%d\n", irr_coor[i][0], irr_coor[i][1], irr_coor[i][2], weight[i]);
	for (int i = 0; i < band_kpt.size(); i++)
		fprintf(fp, "    %lf	%lf	%lf	%d\n", band_kpt[i][0], band_kpt[i][1], band_kpt[i][2], 0);
	fclose(fp);
	printf("Written NEWKPT file!\n");
}

void bandskpt_3d(const char* file1, const char* file2, int kppra)
{
	const char unitfile[] = "UNITCELL";
	EV_unitcell(file1, unitfile);//translate file1 to convention cell
	char   title[MAX_NCHAR];        // Title of system
	double latt;                    // Scaling factor
	int    ifix;                    // Select 0 or normal 1
	int    iflg;                    // Direct 0 or Cartes 1
	int    nant[2];                 // number of atoms and types
	int    typenum[MAX_NELEM];      // number of atoms of each type
	int    elemnum[MAX_NELEM];      // atomic number
	char   elemsym[MAX_NELEM][3];   // atomic symbol
	double vec[3][3];               // lattice vector
	double xyz[MAX_NATOM][3];       // atomic coordinate
	char   fix[MAX_NATOM][3];       // Selective fix on each atom
	FILE* fp = fopen(unitfile, "r");
	if (fp == NULL)
		return;
	readposcar(fp, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);//vec : lattice of convention cell
	fclose(fp);
	double prim_vec[3][3];
	double getdata_vec[3][3];
	transpose_matrix(vec, getdata_vec);
	double getdata_xyz[MAX_NATOM][3];
	for (int i = 0; i < nant[0]; i++)
		for (int j = 0; j < 3; j++)
			getdata_xyz[i][j] = xyz[i][j];
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	SpglibDataset* dataset;
	dataset = spg_get_dataset(getdata_vec, getdata_xyz, types, nant[0], SYMPREAC);
	SpacegroupType spacegrouptype = spgdb_get_spacegroup_type(dataset->hall_number);
	Centering centering = spacegrouptype.centering;
	int spacegroup_number = dataset->spacegroup_number;
	bool time_reversal = time_reversal_Data[spacegroup_number];
	Bravais bravais = find_bravais(centering, spacegroup_number);
	vector<Symmpoint> highsymmpoint;
	if (bravais == TRL)
	{
		//First step : cell that is Niggli reduced in reciprocal space
		double reciprocal_cell_orig[3][3];
		RecMat(vec, reciprocal_cell_orig);
		//This is Niggli-reduced
		double niggli_vec[3][3];
		transpose_matrix(reciprocal_cell_orig, niggli_vec);
		spg_niggli_reduce(niggli_vec, SYMPREAC);
		double reciprocal_cell2[3][3];
		transpose_matrix(niggli_vec, reciprocal_cell2);
		double real_cell2[3][3];
		RecMat(reciprocal_cell2, real_cell2);
		double ka2, kb2, kc2, kalpha2, kbeta2, kgamma2;
		get_cell_params(reciprocal_cell2, ka2, kb2, kc2, kalpha2, kbeta2, kgamma2);
		double para1 = fabs(kb2 * kc2 * cos(kalpha2));
		double para2 = fabs(kc2 * ka2 * cos(kbeta2));
		double para3 = fabs(ka2 * kb2 * cos(kgamma2));
		double M2[3][3] = { 0 };
		if (para1 <= para2 && para1 <= para3)
		{
			M2[0][2] = 1;
			M2[1][0] = 1;
			M2[2][1] = 1;
		}
		if (para2 <= para1 && para2 <= para3)
		{
			M2[0][1] = 1;
			M2[1][2] = 1;
			M2[2][0] = 1;
		}
		if (para3 <= para1 && para3 <= para2)
		{
			M2[0][0] = 1;
			M2[1][1] = 1;
			M2[2][2] = 1;
		}
		// First change of vectors to have | ka3 kb3 cosgamma3 | smallest
		double real_cell2_T[3][3];
		transpose_matrix(real_cell2, real_cell2_T);
		double real_cell3_temp[3][3];
		multi_mat(real_cell2_T, M2, real_cell3_temp);
		double real_cell3[3][3];
		transpose_matrix(real_cell3_temp, real_cell3);
		double reciprocal_cell3[3][3];
		RecMat(real_cell3, reciprocal_cell3);
		double ka3, kb3, kc3, kalpha3, kbeta3, kgamma3;
		get_cell_params(reciprocal_cell3, ka3, kb3, kc3, kalpha3, kbeta3, kgamma3);
		double coskalpha3 = cos(kalpha3);
		double coskbeta3 = cos(kbeta3);
		double coskgamma3 = cos(kgamma3);
		double M3[3][3] = { 0 };
		if (coskalpha3 > 0.0 && coskbeta3 > 0.0 && coskgamma3 > 0.0)  // 1a
		{
			M3[0][0] = 1;
			M3[1][1] = 1;
			M3[2][2] = 1;
		}
		else if (coskalpha3 <= 0.0 && coskbeta3 <= 0.0 && coskgamma3 <= 0.0)  // 1b
		{
			M3[0][0] = 1;
			M3[1][1] = 1;
			M3[2][2] = 1;
		}
		else if (coskalpha3 > 0.0 && coskbeta3 <= 0.0 && coskgamma3 <= 0.0)  // 2a
		{
			M3[0][0] = 1;
			M3[1][1] = -1;
			M3[2][2] = -1;
		}
		else if (coskalpha3 <= 0.0 && coskbeta3 > 0.0 && coskgamma3 > 0.0)  // 2b
		{
			M3[0][0] = 1;
			M3[1][1] = -1;
			M3[2][2] = -1;
		}
		else if (coskalpha3 <= 0.0 && coskbeta3 > 0.0 && coskgamma3 <= 0.0)  // 3a
		{
			M3[0][0] = -1;
			M3[1][1] = 1;
			M3[2][2] = -1;
		}
		else if (coskalpha3 > 0.0 && coskbeta3 <= 0.0 && coskgamma3 > 0.0)  // 3b
		{
			M3[0][0] = -1;
			M3[1][1] = 1;
			M3[2][2] = -1;
		}
		else if (coskalpha3 <= 0.0 && coskbeta3 <= 0.0 && coskgamma3 > 0.0)  // 4a
		{
			M3[0][0] = -1;
			M3[1][1] = -1;
			M3[2][2] = 1;
		}
		else if (coskalpha3 > 0.0 && coskbeta3 > 0.0 && coskgamma3 <= 0.0)  // 4b
		{
			M3[0][0] = -1;
			M3[1][1] = -1;
			M3[2][2] = 1;
		}
		double real_cell3_T[3][3];
		transpose_matrix(real_cell3, real_cell3_T);
		double real_vec_final_temp[3][3];
		multi_mat(real_cell3_T, M3, real_vec_final_temp);
		double real_vec_final[3][3];
		transpose_matrix(real_vec_final_temp, real_vec_final);
		double reciprocal_cell_final[3][3];
		RecMat(real_vec_final, reciprocal_cell_final);
		highsymmpoint = find_highsymmpoint(bravais, reciprocal_cell_final, spacegroup_number, time_reversal);
		//direct -> carts in convention vec
		direct_to_carts(iflg, vec, xyz, nant);
		//carts ->direct in real_vec_final
		carts_to_direct(iflg, real_vec_final, xyz, nant);
		transpose_matrix(real_vec_final, prim_vec);
	}
	else
	{
		transpose_matrix(vec, prim_vec);
		highsymmpoint = find_highsymmpoint(bravais, vec, spacegroup_number, time_reversal);
	}
	int num_primitive_atom = Find_Primitive_Cell(prim_vec, xyz, types, nant[0], SYMPREAC, -1.0);
	if (num_primitive_atom)
	{
		double vec_prim[3][3];
		transpose_matrix(prim_vec, vec_prim);
		nant[0] = num_primitive_atom;
		int* primitive_typenum = (int*)malloc(sizeof(int) * nant[1]);
		translate_type_typenum(types, nant, primitive_typenum);
		FILE* fp2 = fopen(file2, "w");
		savposcar(fp2, title, latt, ifix, iflg, nant, primitive_typenum, elemnum, elemsym, vec_prim, xyz, fix);
		fclose(fp2);
		write_kpt(highsymmpoint, kppra, time_reversal);
		remove("structure_operator");
		free(types);
		return;
	}
	else
	{
		free(types);
		return;
	}
}

vector<Symmpoint> get_2Dpath(double a1, double a2, double cosgamma)
{
	vector<Symmpoint> point;
	map<char, Symmpoint> path_point;
	string path;
	string s;
	if (Isequal(a1, a2, SYMPREAC) && Isequal(cosgamma, -0.5, SYMPREAC)) //Hexagonal
	{
		s = "GKM";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(1.0 / 3, 1.0 / 3, 0, "K");
		Symmpoint p3(0.5, 0, 0, "M");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		path = "GMKG";
	}
	if (Isequal(a1, a2, SYMPREAC) && Isequal(cosgamma, 0, SYMPREAC)) //Square
	{
		s = "GMX";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, 0.5, 0, "M");
		Symmpoint p3(0.5, 0, 0, "X");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		path = "GXMG";
	}
	if (!Isequal(a1, a2, SYMPREAC) && Isequal(cosgamma, 0, SYMPREAC)) //Rectangular
	{
		s = "GSXY";
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, 0.5, 0, "S");
		Symmpoint p3(0.5, 0, 0, "X");
		Symmpoint p4(0, 0.5, 0, "Y");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		path = "GXSYGS";
	}
	if (Isequal(a1 / 2, a2 * cosgamma, SYMPREAC) && !Isequal(cosgamma, 0, SYMPREAC)) //Centered Rectangular
	{
		s = "GXQCH"; //Q->H1
		double x1 = (1.0 - a1 * cosgamma / a2) / 2 / (1 - cosgamma * cosgamma);
		double x2 = 1.0 / 2 - x1 * a2 * cosgamma / a1;
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, 0, 0, "X");
		Symmpoint p3(1 - x1, x2, 0, "H_1");
		Symmpoint p4(0.5, 0.5, 0, "C");
		Symmpoint p5(x1, 1 - x2, 0, "H");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		path = "GXQCHG";
	}
	if (!Isequal(a1, a2, SYMPREAC) && !Isequal(cosgamma, 0, SYMPREAC)) //Oblique
	{
		s = "GXQCHY";//Q->H1
		double x1 = (1.0 - a1 * cosgamma / a2) / 2 / (1 - cosgamma * cosgamma);
		double x2 = 1.0 / 2 - x1 * a2 * cosgamma / a1;
		Symmpoint p1(0, 0, 0, "GAMMA");
		Symmpoint p2(0.5, 0, 0, "X");
		Symmpoint p3(1 - x1, x2, 0, "H_1");
		Symmpoint p4(0.5, 0.5, 0, "C");
		Symmpoint p5(x1, 1 - x2, 0, "H");
		Symmpoint p6(0, 0.5, 0, "Y");
		point.push_back(p1);
		point.push_back(p2);
		point.push_back(p3);
		point.push_back(p4);
		point.push_back(p5);
		point.push_back(p6);
		path = "GXQCHYG";
	}
	for (int i = 0; i < s.size(); i++)
		path_point.insert(pair <char, Symmpoint>(s[i], point[i]));
	write_highsymmpoint(point, path, path_point, 1);
	return 	get_path(path, path_point);
}

void bandskpt_2d(const char* file1, const char* file2, int kppra)
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
	double xyz[MAX_NATOM][3];       // atomic coordinate
	char   fix[MAX_NATOM][3];       // Selective fix on each atom
	FILE* fp = fopen(file1, "r");
	if (fp == NULL)
	{
		printf("%s IS NOT EXIST!\n", file1);
		return;
	}
	readposcar(fp, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);
	fclose(fp);
	double ntemp_vec[3][3];
	transpose_matrix(vec, ntemp_vec);
	int* types = (int*)malloc(sizeof(int) * nant[0]);
	translate_typenum_type(types, nant[1], typenum);
	int num_primitive_atom = spg_find_primitive(ntemp_vec, xyz, types, nant[0], SYMPREAC);
	double prim_vec[3][3];
	transpose_matrix(ntemp_vec, prim_vec);
	double a1 = sqrt(prim_vec[0][0] * prim_vec[0][0] + prim_vec[0][1] * prim_vec[0][1] + prim_vec[0][2] * prim_vec[0][2]);
	double a2 = sqrt(prim_vec[1][0] * prim_vec[1][0] + prim_vec[1][1] * prim_vec[1][1] + prim_vec[1][2] * prim_vec[1][2]);
	double cosgamma = (prim_vec[0][0] * prim_vec[1][0] + prim_vec[0][1] * prim_vec[1][1] + prim_vec[0][2] * prim_vec[1][2]) / a1 / a2;
	vector<Symmpoint> highsymmpoint = get_2Dpath(a1, a2, cosgamma);
	write_kpt(highsymmpoint, kppra, 1);
	free(types);
}

void GetFermiMesh(const char file[], int kppra, double kspac, int mesh[3], char label[], char ikmesh)
{
	vector<int> weight;
	vector<vector<double> > irr_coor;
	int* imesh = ir_reciprocal_mesh(file, kppra, kspac, ikmesh, mesh, label, weight, irr_coor);
	FILE* fp1 = fopen("FERMIKPT", "w");
	fprintf(fp1, "%c	%d	%d	%d	# Parameters to Generate KPOINTS (Don't Edit This Line)\n", ikmesh, imesh[0], imesh[1], imesh[2]);
	fprintf(fp1, "	%d\n", irr_coor.size());
	fprintf(fp1, "Reciprocal lattice\n");
	for (int i = 0; i < irr_coor.size(); i++)
		fprintf(fp1, "    %lf	%lf	%lf	%d\n", irr_coor[i][0], irr_coor[i][1], irr_coor[i][2], weight[i]);
	fclose(fp1);
	printf("Written FERMIKPT file!\n");
}

void Get3DbandMesh(const char file[], int kppra, double kspac, int kpt[3], char label[], char ikmesh)
{
	FILE* fp = fopen(file, "r");
	POSCAR pos;
	readposcar(fp, pos);
	fclose(fp);
	int* mesh = (int*)malloc(sizeof(int) * 3);
	_matrix kvec;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			kvec.mat[i][j] = pos.vec[i][j];
	_matrix rec = RecMat(kvec);
	bool has_vacuum_slab[3];
	has_vacuum_slab_or_not(pos, has_vacuum_slab);
	if (!strcmp(label, "ka"))
		mesh = kpta(rec, kppra, pos.nant[0], ikmesh, has_vacuum_slab);
	if (!strcmp(label, "kv"))
		mesh = kpvf(rec, kspac, pos.nant[0], ikmesh, has_vacuum_slab);
	if (!strcmp(label, "kpt"))
	{
		mesh[0] = kpt[0];
		mesh[1] = kpt[1];
		mesh[2] = kpt[2];
	}
	int num_kpt = mesh[0] * mesh[1] * 1; //2D
	vector<vector<double> > position(num_kpt, vector<double>(3));
	vector<int> weight(num_kpt, 1);
	int cnt = 0;
	for (int j = 0; j < mesh[1]; j++)
		for (int i = 0; i < mesh[0]; i++)
		{
			position[cnt][0] = -0.5 + 1.0 / (mesh[0] - 1) * i;
			position[cnt][1] = -0.5 + 1.0 / (mesh[1] - 1) * j;
			position[cnt][2] = 0;
			if (i == mesh[0] - 1 || j == mesh[1] - 1)
				weight[cnt] = 0;
			cnt++;
		}
	FILE* fp1 = fopen("3DbandKPT", "w");
	fprintf(fp1, "KMesh	%d	%d	%d \n", mesh[0], mesh[1], 1);
	fprintf(fp1, "	%d\n", num_kpt);
	fprintf(fp1, "Reciprocal lattice\n");
	for (int i = 0; i < num_kpt; i++)
		fprintf(fp1, "	%lf	%lf	%lf	%d\n", position[i][0], position[i][1], position[i][2], weight[i]);
	fclose(fp1);
	printf("Wriiten 3DbandKPT file!\n");
}

void get_tran_matrix(const char tranfile[],double tran[3][3])
{
	FILE* fp = fopen(tranfile, "r");
	char buf[1024];
	if (fp == NULL)
	{
		printf("%s is not exist\n", tranfile);
		return;
	}
	fgets(buf, 1024, fp);
	for (int i = 0; i < 3; i++)
	{
		if (fgets(buf, 1024, fp) == NULL)
		{
			printf("Can't find tran_matrix from %s!\n", tranfile);
			return;
		}
		else
			sscanf(buf, "%lf%lf%lf", &tran[i][0], &tran[i][1], &tran[i][2]);
	}
	fclose(fp);
}

void getunfoldkpoints(const char file[], int ikppra, double kspac, int kpt[], char label[], double kpra, char ikmesh, const char tranfile[])
{
	vector<int> weight;
	vector<vector<double> > irr_coor;
	int* mesh = ir_reciprocal_mesh(file, ikppra, kspac, ikmesh, kpt, label, weight, irr_coor);
	char   title[MAX_NCHAR];        // Title of system
	double latt;                    // Scaling factor
	int    ifix;                    // Select 0 or normal 1
	int    iflg;                    // Direct 0 or Cartes 1
	int    nant[2];                 // number of atoms and types
	int    typenum[MAX_NELEM];      // number of atoms of each type
	int    elemnum[MAX_NELEM];      // atomic number
	char   elemsym[MAX_NELEM][3];   // atomic symbol
	double vec[3][3];               // lattice vector
	double xyz[MAX_NATOM][3];       // atomic coordinate
	char   fix[MAX_NATOM][3];       // Selective fix on each atom
	FILE* fp1 = fopen(file, "r");
	if (fp1 == NULL)
	{
		printf("%s IS NOT EXIST!\n", file);
		return;
	}
	readposcar(fp1, title, latt, ifix, iflg, nant, typenum, elemnum, elemsym, vec, xyz, fix);//vec : lattice of supercell
	fclose(fp1);
	FILE* fp2 = fopen("NEWKPATH", "r");
	if (fp2 == NULL)
	{
		printf("NEWKPATH IS NOT EXIST!\n");
		return;
	}
	char buf[1024];
	int line = 0;
	vector<vector<double> > high_symmetry;
	while (fgets(buf, 1024, fp2) != NULL)
	{
		line++;
		if (strspn(buf, " \t\n") == strlen(buf))
			continue;
		if (line > 4)
		{
			vector<double> hs(3, 0);
			sscanf(buf, "%lf%lf%lf", &hs[0], &hs[1], &hs[2]);
			high_symmetry.push_back(hs);
		}
	}
	fclose(fp2);
	//Transfer supercell lattice
	double tran[3][3];
	get_tran_matrix(tranfile, tran);
	double inv_tran[3][3];
	brinv(tran, inv_tran);
	multi_mat(inv_tran, vec, vec);
	double rec_vec[3][3];
	RecMat(vec, rec_vec);
	vector<vector<double> > carts_kpt(high_symmetry.size(), vector<double>(3));
	for (int i = 0; i < high_symmetry.size(); i++)
	{
		double p[3] = { high_symmetry[i][0],high_symmetry[i][1],high_symmetry[i][2] };
		for (int j = 0; j < 3; j++)
		{
			carts_kpt[i][j] = p[0] * rec_vec[0][j] + p[1] * rec_vec[1][j] + p[2] * rec_vec[2][j];
		}
	}
	vector<int> point(high_symmetry.size() / 2);
	for (int i = 0; i < point.size(); i++)
	{
		double length = 0;
		for (int j = 0; j < 3; j++)
			length += (carts_kpt[2 * i + 1][j] - carts_kpt[2 * i][j]) * (carts_kpt[2 * i + 1][j] - carts_kpt[2 * i][j]);
		point[i] = max(5, (int)(sqrt(length) / kpra));
	}
	vector<vector<double> > band_kpt;//kmesh_primcell
	for (int i = 0; i < point.size(); i++)
	{
		for (int j = 0; j < point[i]; j++)
		{
			double x = high_symmetry[2 * i][0] + (high_symmetry[2 * i + 1][0] - high_symmetry[2 * i][0]) / (point[i] - 1) * j;
			double y = high_symmetry[2 * i][1] + (high_symmetry[2 * i + 1][1] - high_symmetry[2 * i][1]) / (point[i] - 1) * j;
			double z = high_symmetry[2 * i][2] + (high_symmetry[2 * i + 1][2] - high_symmetry[2 * i][2]) / (point[i] - 1) * j;
			vector<double> temp = { x,y,z };
			band_kpt.push_back(temp);
		}
	}
	vector<vector<double> > supercell_band_kpt(band_kpt.size(), vector<double>(3)); //kmesh_supercell
	vector<vector<int> > gvector_supercell(band_kpt.size(), vector<int>(3));
	for (int i = 0; i < band_kpt.size(); i++)
	{
		double p[3] = { band_kpt[i][0],band_kpt[i][1],band_kpt[i][2] };
		for (int j = 0; j < 3; j++)
		{
			supercell_band_kpt[i][j] = p[0] * tran[0][j] + p[1] * tran[1][j] + p[2] * tran[2][j];
			gvector_supercell[i][j] = (int)(supercell_band_kpt[i][j] + 0.5);
			supercell_band_kpt[i][j] -= gvector_supercell[i][j];
		}
	}
	FILE* fp = fopen("NEWKPT", "w");
	//mesh0 mesh1 mesh2 num_scf num_hse num_path num_path0 ...
	fprintf(fp, "%d	%d	%d	%d", mesh[0], mesh[1], mesh[2], irr_coor.size());
	fprintf(fp, "	%d	%d", supercell_band_kpt.size(), point.size());
	for (int i = 0; i < point.size(); i++)
		fprintf(fp, "	%d", point[i]);
	fprintf(fp, "          # Parameters to Generate KPOINTS (Do NOT Edit This Line)\n");
	fprintf(fp, "	%d\n", irr_coor.size() + supercell_band_kpt.size());
	fprintf(fp, "Reciprocal lattice\n");
	for (int i = 0; i < irr_coor.size(); i++)
		fprintf(fp, "    %lf    %lf    %lf     %d\n", irr_coor[i][0], irr_coor[i][1], irr_coor[i][2], weight[i]);
	for (int i = 0; i < band_kpt.size(); i++)
		fprintf(fp, "    %lf    %lf    %lf     %d\n", supercell_band_kpt[i][0], supercell_band_kpt[i][1], supercell_band_kpt[i][2], 0);
	fclose(fp);
	printf("Written NEWKPT file!\n");
	FILE* fp3 = fopen("KPOINTS_MAPPING_TABLE", "w");
	fprintf(fp3, "Unfolded K-Points in Primitive-Cell                  Folded K-Points in Super-Cell                  G-vectors for Super-Cell\n");
	for (int i = 0; i < band_kpt.size(); i++)
		fprintf(fp3, "    %lf    %lf    %lf                  %lf    %lf    %lf                  %d     %d     %d\n", 
			band_kpt[i][0], band_kpt[i][1], band_kpt[i][2],
			supercell_band_kpt[i][0], supercell_band_kpt[i][1], supercell_band_kpt[i][2], 
			gvector_supercell[i][0], gvector_supercell[i][1], gvector_supercell[i][2]);
	fclose(fp3);
	printf("Written KPOINTS_MAPPING_TABLE file!\n");
}
