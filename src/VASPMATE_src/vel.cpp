#include"../../include/VASPMATE_include/vel.h"
#include"../../include/VASPMATE_include/mapdata.h"

using namespace std;
using namespace vel;

int vel::operator_vel(int argc, char* argv[])
{
    if (argc == 2)
        write_vel(300,"INPOS");
    else if (!strcmp(argv[2], "INPOS")) //VASPMATE --vel INPOS -T [temperature]
    {
        int par_T;
        int check_T = check_para(argc, argv, "-T", 1, &par_T);
        if(check_T == 0)
            return -1;
        else if (check_T == 1)
            write_vel(atof(argv[par_T + 1]),"INPOS");
        else
            return -1;
    }
	else if(!strcmp(argv[2], "-T")) //VASPMATE --vel -T [temperature]
    {
        int par_T;
        int check_T = check_para(argc, argv, "-T", 1, &par_T);
        if(check_T == 0)
            return -1;
        else if (check_T == 1)
            write_vel(atof(argv[par_T + 1]),"INPOS");
        else
            return -1;
    }
	return 0;
}
double cal_vel(double temperature, double mass)
 {  
	std::random_device rd;
    std::default_random_engine gen(rd());
	const double kB = 1.38064852e-23;
    std::normal_distribution<double> dist(0.0, 1.0);  
    return dist(gen)*sqrt(kB*temperature/mass);  
}  //generate a velocity component that follows the Maxwell Boltzmann distribution and return it
void copyfileContent(const char* inpos,const char* velopos)
{
	
	
};
void vel::write_vel(double temperature, const char* inpos)
{
	FILE* fp1 = fopen(inpos,"r");
	POSCAR pos;
	readposcar(fp1,pos);                 // Selective fix on each atom
	fclose(fp1);
    FILE* fp = fopen("VELOPOS", "w");
	savposcar(fp,pos);
	fprintf(fp,"\n");
    for (int j = 0; j < pos.nant[1]; j++)
	{
		for (int l = 0; l < pos.typenum[j]*3; l++)
		{
			double velocity = cal_vel(temperature, atomicMasses[pos.elemsym[j]]*1.66e-27);
			fprintf(fp,"%12.6f",velocity);
			fprintf(fp,"\t");
			if ((l + 1)%3 == 0)
				fprintf(fp,"\n");
		}
	}
    fclose(fp);
    printf("Write VELOPOS file!\n");
};
 