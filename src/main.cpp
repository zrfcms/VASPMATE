#include <cstring>
#include <iostream>
#include"../include/VASPMATE_include/VASPMATE_main.h"
#include "../include/Stochastics_include/Stochastics_main.h"
#include "../include/Evolutionary_include/Evolutionary_main.h"
#include "../include/help_include/help.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    if (argc > 1 && !strcmp("--sto", argv[1]))
    {
        //sto
        Stochastics_main(argc, argv);
    }
    else if (argc > 1 && !strcmp("--evo", argv[1]))
    {
        //evo
        Evolutionary_main(argc, argv);
    }
    else
    {
        int flag_help = VASPMATE_main(argc, argv);
        //cout << flag_help << endl;
        if (flag_help != 0)
        {
            vaspmate_help(flag_help, argc, argv);
        }
    }
    return 0;
}
