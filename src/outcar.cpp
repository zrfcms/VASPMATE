#include<outcar.h>
double get_energy()
{
	double energy = 0;
	FILE* fp_ = fopen("OUTCAR", "r");
	if (fp_ == NULL)
	{
		printf("OUTCAR IS NOT EXIST!\n");
		return energy;
	}
	char buf[1024];
	while (fgets(buf, 1024, fp_) != NULL)
	{
		if (strstr(buf, "energy(sigma->0)") != NULL)
		{
			sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%lf", &energy);
		}
	}
	fclose(fp_);
	return energy;
}