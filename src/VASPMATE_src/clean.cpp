#include"../../include/VASPMATE_include/clean.h"

using namespace std;

void rm_file(const vector<string> &file_vector)
{
    for (const auto &file_name : file_vector) 
    {
        if (remove(file_name.c_str()) == 0)
            cout << "Successfully deleted file: " << file_name << endl;
    }
}

void rm_other_file()
{
    const char* dir_path = ".";
    DIR* dir = opendir(dir_path);
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) 
    {
        string file_name(entry->d_name);
        if (file_name == "." || file_name == "..")
            continue;
        struct stat statbuf;
        if (stat(file_name.c_str(), &statbuf) == 0 && S_ISREG(statbuf.st_mode)) 
        {
            // If the file is not INCAR, POSCAR, KPOINTS, or POTCAR, remove it
            if (file_name != "INCAR" && file_name != "POSCAR" && file_name != "KPOINTS" && file_name != "POTCAR" &&
                file_name.find(".sh") == string::npos && file_name.find(".ksh") == string::npos) 
            {
                if (remove(file_name.c_str()) == 0) 
                   cout << "Successfully deleted file: " << file_name << endl;
            }
        }
    }
    closedir(dir);
}

int clean_operat(int argc, char* argv[])
{
	if (argc == 2)
    {
        rm_other_file();
    }
    else if (!strcmp("-ssc", argv[2]))
	{
        vector<string> file_vector = {""};
        rm_file(file_vector);
	}
    else if (!strcmp("-nsc", argv[2]))
	{
        vector<string> file_vector = {""};
        rm_file(file_vector);
	}
    else if (!strcmp("-rlx", argv[2]))
	{
        vector<string> file_vector = {"PROCAR","DOSCAR","IBZKPT","EIGENVAL","CHGCAR","CHG","PCDAT","XDATCAR"};
        rm_file(file_vector);
	}
    else if (!strcmp("-mds", argv[2]))
	{
        vector<string> file_vector = {"PROCAR","DOSCAR","IBZKPT","EIGENVAL","CHGCAR","CHG","PCDAT"};
        rm_file(file_vector);
	}
    else if (!strcmp("-dos", argv[2]))
	{
        vector<string> file_vector = {"PROCAR","IBZKPT","EIGENVAL","CHGCAR","CHG","PCDAT","XDATCAR"};
        rm_file(file_vector);
	}
    else if (!strcmp("-band", argv[2]))
	{
        vector<string> file_vector = {"DOSCAR","IBZKPT","EIGENVAL","CHGCAR","CHG","PCDAT","XDATCAR"};
        rm_file(file_vector);
	}
    else if (!strcmp("-chg", argv[2]))
	{
        vector<string> file_vector = {"PROCAR","DOSCAR","IBZKPT","EIGENVAL","PCDAT","XDATCAR"};
        rm_file(file_vector);
	}
    else if (!strcmp("-elf", argv[2]))
	{
        vector<string> file_vector = {"PROCAR","DOSCAR","IBZKPT","EIGENVAL","CHGCAR","CHG","PCDAT","XDATCAR"};
        rm_file(file_vector);
	}
    else
        return -1;
    return 0;
}