/*
 * Plot KS pseudo wavefunction at specific kpoint, specific band.
 * An example that uses wavehigh class
 */

#include"../include/plotWave.h"

using ionizing::WAVEHIGH;

/*
	VASPMATE --wfun [options]
	-b -k -h -l
*/
void wavefun(int argc, char* argv[])
{
	bool
		is_list_info = false,
		is_print_example = false,
		is_real = false;
	int
		ikpoint = 1,
		iband = 1,
		ispin = 0;

	std::string prefix{ "wfn" };
	std::string wavecar_fname{ "WAVECAR" };
	std::string poscar_fname{ "POSCAR" };

	for (int i = 2; i < argc; i++)
	{
		if (!strcmp(argv[i], "-k") || !strcmp(argv[i], "-kpoint"))
			ikpoint = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "-b") || !strcmp(argv[i], "-band"))
			iband = atoi(argv[i + 1]);
		if (!strcmp(argv[i], "-e") || !strcmp(argv[i], "-example"))
			is_print_example = 1;
		if (!strcmp(argv[i], "-l") || !strcmp(argv[i], "-list"))
			is_list_info = 1;
	}

	if (is_print_example) {
		std::cout << "Example:" << std::endl
			<< "    For wave plot: " << argv[0]
			<< " -wfun -kpoint=1 -band=1" << std::endl
			<< "       Or shortly: " << argv[0]
			<< " -wfun -k 1 -b 1" << std::endl
			<< "    For wave info: " << argv[0] << " -list or " << argv[0] << " -l" << std::endl;
			return ;
	}

	WAVEHIGH wave{ wavecar_fname.c_str() };

	if (is_list_info) {
		wave.printInfo(std::cout);
	}

	for (; ispin < wave.getInfo()._nSpin; ispin++)
		wave.plotWave(ispin, ikpoint - 1, iband - 1, poscar_fname.c_str(), prefix.c_str(), is_real);
}


