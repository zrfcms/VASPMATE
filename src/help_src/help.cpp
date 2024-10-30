#include"../../include/help_include/help.h"
using std::cout;
using std::multimap;
using std::endl;
using std::string;

multimap<string, string> help_command = {
    {"--v","VASPMATE --v/--version"},
    {"--version","VASPMATE --v/--version"},
    {"--l","VASPMATE --l/--license"},
    {"--license","VASPMATE --l/--license"},
    {"--h","VASPMATE --h/--help"},
    {"--help","VASPMATE --h/--help"},
    {"--cif2pos", "VASPMATE --cif2pos file1(INCIF) file2(NEWPOS)"},
    {"--pos2cif", "VASPMATE --pos2cif file1(INPOS) file2(NEWCIF)"},
    {"--del","VASPMATE --del"},
    {"--chg2cub", "VASPMATE --chg2cub file1(chgcar) file2(cube)"},
    {"--prim", "VASPMATE --prim file1(INPOS) file2(PRIMPOS)"},
    {"--unit", "VASPMATE --unit file1(INPOS) file2(UNITPOS)"},
    {"--symm", "VASPMATE --symm file(INPOS)"},
    {"--num", "VASPMATE --num file(INPOS)"},
    {"--atom", "VASPMATE --atom file(INPOS)"},
    {"--cell", "VASPMATE --cell file(INPOS)"},
    {"--vol", "VASPMATE --vol file(INPOS)"},
    {"--super", "VASPMATE --super file1(INPOS) file2(SUPEPOS) (-np) [sn1 sn2 sn3]"},
    {"--affine", "VASPMATE --affine file1(INPOS) file2(AFFPOS) -txx/-tyy/-tzz [strain]"},
    {"--affine", "VASPMATE --affine file1(INPOS) file2(AFFPOS) -sxy/-syz/-szx [strain]"},
    {"--affine", "VASPMATE --affine file1(INPOS) file2(AFFPOS) -pxy/-pyz/-pzx [strain]"},
    {"--affine", "VASPMATE --affine file1(INPOS) file2(AFFPOS) -tension [xx/yy/zz] [init_strain step_length step_num]"},
    {"--affine", "VASPMATE --affine file1(INPOS) file2(AFFPOS) -simshear/-purshear [xy/yz/zx] [init_strain step_length step_num]"},
    {"--alias", "VASPMATE --alias file1(INPOS) file2(ALIPOS) -txx/-tyy/-tzz [strain] [position]"},
    {"--alias", "VASPMATE --alias file1(INPOS) file2(ALIPOS) -sxy/-syz/-szx [strain1] [strain2] [position]"},
    {"--alias", "VASPMATE --alias file1(INPOS) file2(ALIPOS) -tensi [xx/yy/zz] [istart iend ispac] [position]"},
    {"--alias", "VASPMATE --alias file1(INPOS) file2(ALIPOS) -shear [xy/yz/zx] [istart1 iend1 ispac1] [istart2 iend2 ispac2] [position]"},
    {"--proj", "VASPMATE --proj file1(INPOS) file2(PROJPOS) -rot [rotx roty rotz]"},
    {"--proj", "VASPMATE --proj file1(INPOS) file2(PROJPOS) -ind [pvh pvk pvl] [uvu uvv uvw]"},
    {"--proj", "VASPMATE --proj file1(INPOS) file2(PROJPOS) -mat [mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33]"},
    {"--redef", "VASPMATE --redef file1(INPOS) file2(REDEPOS) (-par) [vect11 vect12 vect13 vect21 vect22 vect23 vect31 vect32 vect33] [a1 x1 a2 x2]"},
    {"--ieee", "VASPMATE --ieee file1(INPOS) file2(IEEEPOS)"},
    {"--fixc", "VASPMATE --fixc file1(INPOS) file2(FIXPOS) (-par) [axis] [m n] [fix-tag]"},
    {"--fixa", "VASPMATE --fixa file1(INPOS) file2(FIXPOS) (-par) [number] [fix-tag]"},
    {"--fixe", "VASPMATE --fixe file1(INPOS) file2(FIXPOS) (-par) [element1] [fix-tag] / [element2] [fix-tag] / "},
    {"--swap", "VASPMATE --swap file(INPOS) file2(SWAPOS) (-elem) [elem1 elem2...]"},
    {"--ufix", "VASPMATE --ufix file(INPOS)"},
    {"--cartes", "VASPMATE --cartes file1(INPOS) file2(NEWPOS)"},
    {"--direct", "VASPMATE --direct file1(INPOS) file2(NEWPOS)"},
    {"--sortc", "VASPMATE --sortc file1(INPOS) file2(SORTPOS) (-par) [axis]"},
    {"--sorte", "VASPMATE --sorte file1(INPOS) file2(SORTPOS) (-par) [element-list]"},
    {"--movc", "VASPMATE --movc file1(INPOS) file2(MOVEPOS) (-par) [axis] [min max] [dx dy dz]"},
    {"--movd", "VASPMATE --movd file1(INPOS) file2(MOVEPOS) (-par) [axis] [min max] [dx dy dz]"},
    {"--i", "VASPMATE --i [option-list]"},
    {"--inp", "VASPMATE --inp [option-list]"},
    {"--ia", "VASPMATE --ia [keyword] [value] / [keyword] [value] / … ] (--i_app or --i_append)"},
    {"--id", "VASPMATE --id [keyword-list] (--i_del or --i_delete)"},
    {"--irm", "VASPMATE --irm [keyword-list] (--i_rem or --i_remove)"},
    {"--irp", "VASPMATE --irp [keyword] [value] / [keyword] [value] / … ] (--i_rep or --i_replace)"},
    {"--iu", "VASPMATE --iu (INPOS) (-a)"},
    {"--iu", "VASPMATE --iu (INPOS) (-t) [Table]"},
    {"--ldau", "VASPMATE --ldau (INPOS) (-a)"},
    {"--ldau", "VASPMATE --ldau (INPOS) (-t) [Table]"},
    {"--iv", "VASPMATE --iv (INPOS) (-a)"},
    {"--iv", "VASPMATE --iv (INPOS) (-t) (d3b/d2/d3z/b86/b88/pbe/rpbe)"},
    {"--ivdw", "VASPMATE --ivdw (INPOS) (-a)"},
    {"--ivdw", "VASPMATE --ivdw (INPOS) (-t) (d3b/d2/d3z/b86/b88/pbe/rpbe)"},
    {"--k", "VASPMATE --k (-par) [k1 k2 k3] [kmesh] [kscheme]"},
    {"--km", "VASPMATE --km (-par) [k1 k2 k3] [kmesh] [kscheme]"},
    {"--kpt", "VASPMATE --kpt (-par) [k1 k2 k3] [kmesh] [kscheme]"},
    {"--kmseh", "VASPMATE --kmesh (-par) [k1 k2 k3] [kmesh] [kscheme]"},
    {"--ka", "VASPMATE --ka file(INPOS) (-par) [kppra] [kscheme]"},
    {"--kv", "VASPMATE --kv file(INPOS) (-par) [kspac] [kscheme]"},
    {"--pot", "VASPMATE --pot file(INPOS)"},
    {"--pot", "VASPMATE --pot file(INPOS) -type (-suff postfix1 postfix2…)"},
    {"--pot", "VASPMATE --pot file(INPOS) -t type (-suff postfix1 postfix2…)"},
    {"--pote", "VASPMATE --pote -type [elem1 elem2...] (-suff [postfix1 postfix2…])"},
    {"--pote", "VASPMATE --pote -t [type] -e [elem1 elem2...] (-suff [postfix1 postfix2…])"},
    {"--check", "VASPMATE --check (-in) (inp/pos/kpt/pot)"},
    {"--check", "VASPMATE --check (-out) (rlx/stc/md)"},
    {"--clean", "VASPMATE --clean (-ssc/-nsc/-rlx/-mds)"},
    {"--std3d", "VASPMATE --std3d file1(INPOS) file2(STD3POS)"},
    {"--std2d", "VASPMATE --std2d file1(INPOS) file2(STD2POS)"},
    {"--ka3d", "VASPMATE --ka3d file(INPOS) (-par) [value]"},
    {"--ka2d", "VASPMATE --ka2d file(INPOS) (-par) [value]"},
    {"--band", "VASPMATE --band -b"},
    {"--band", "VASPMATE --band -a"},
    {"--band", "VASPMATE --band -e"},
    {"--band", "VASPMATE --band -s [atom-index] [element-index]"},
    {"--band", "VASPMATE --band -sa [atom-index]"},
    {"--band", "VASPMATE --band -m/-ma/-me [atom-index]/[element-index]"},
    {"--band", "VASPMATE --band -o/-oa/-oe [atom-index]/[element-index] [orbit-index]"},
    {"--band", "VASPMATE --band -id/-index"},
    {"--band", "VASPMATE --band -bg"},
    {"--band", "VASPMATE --band -em [point] [direct]…(-nb [band-index]) (-np [fitting-point])"},
    {"--kahse", "VASPMATE --kahse (INPOS) (-par) [kppra] [resolution] [kscheme]"},
    {"--kvhse", "VASPMATE --kvhse (INPOS) (-par) [kspac] [resolution] [kscheme]"},
    {"--kmhse", "VASPMATE --kmhse (-par) k_1 k_2 k_3 [resolution] [kscheme]"},
    {"--band", "VASPMATE --band -hb"},
    {"--band", "VASPMATE --band -ha"},
    {"--band", "VASPMATE --band -he"},
    {"--band", "VASPMATE --band -hs [atom-index] [element-index]"},
    {"--band", "VASPMATE --band -hsa [atom-index]"},
    {"--band", "VASPMATE --band -hse [element-index]"},
    {"--band", "VASPMATE --band -hm/-hma/-hme [atom-index]/[element-index]"},
    {"--band", "VASPMATE --band -ho/-hoa/-hoe [atom-index]/[element-index]"},
    {"--band", "VASPMATE --band -hbg"},
    {"--dos", "VASPMATE --dos -t [none]"},
    {"--dos", "VASPMATE --dos -a [none]"},
    {"--dos", "VASPMATE --dos -e [none]"},
    {"--dos", "VASPMATE --dos -s [atom-index]"},
    {"--dos", "VASPMATE --dos -sa [atom-index]"},
    {"--dos", "VASPMATE --dos -se [element-index]"},
    {"--dos", "VASPMATE --dos -m/-me/-ma [select-list]"},
    {"--dos", "VASPMATE --dos -o/-oe/-oa [atom&orbit-list]"},
    {"--dos", "VASPMATE --dos -bc"},
    {"--bader", "VASPMATE --bader -comb (AECCAR0) (AECCAR2) (factor1) (factor2)"},
    {"--bader", "VASPMATE --bader -d/-derive/-calc"},
    {"--neb", "VASPMATE --neb -sim file1(INIPOS) file2(FINPOS)"},
    {"--neb", "VASPMATE --neb -ins (number) (-line/-idpp)"},
    {"--neb", "VASPMATE --neb -out"},
    {"--vcd", "VASPMATE --vcd -split file(CHGCAR)"},
    {"--vcd", "VASPMATE --vcd -sum [File_list]"},
    {"--vcd", "VASPMATE --vcd -diff [File_list]"},
    {"--pcd", "VASPMATE --pcd -ib/-ik/-en/-ef [value]"},
    {"--wfn", "VASPMATE --wfn -ik [kpoint_index] -ib [band_index]"},
    {"--fska", "VASPMATE --fska (INPOS) (-par) [kppra] [kscheme]"},
    {"--fskv", "VASPMATE --fskv (INPOS) (-par) [kspac] [kscheme]"},
    {"--dat2csv", "VASPMATE --file2csv inputfile outputfile"},
    {"--db2js", "VASPMATE --db2js inputfile outputfile"},
    {"--js2db", "VASPMATE --js2db inputfile outputfile"},
    {"--csv2dat", "VASPMATE --csv2file inputfile outputfile"},
    {"--mds", "VASPMATE --mds -nve -ts [steps]"},
    {"--mds", "VASPMATE --mds -nvt -T [temperature] -ts [steps]"},
    {"--mds", "VASPMATE --mds -npt -T [temperature] -p [pressure] -ts [steps]"},
    {"--mds", "VASPMATE --mds -et"},
    {"--mds", "VASPMATE --mds -mt"},
    {"--mds", "VASPMATE --mds -ns [time step numbers]"},
    {"--dbs", "VASPMATE --dbs database(vasp.db)"},
    {"--dbs", "VASPMATE --dbs vasp.db -tables"},
    {"--dbs", "VASPMATE --dbs vasp.db -schema [table]"},
    {"--dbs", "VASPMATE --dbs vasp.db -table [table-name] -create [keyword] [type]"},
    {"--dbs", "VASPMATE --dbs vasp.db -table [table-name] -alter [table-name-new]"},
    {"--dbs", "VASPMATE --dbs vasp.db -table [table-name] -drop"},
    {"--dbs", "VASPMATE --dbs vasp.db -table [table-name] -set [keyword] [type]"},
    {"--dbs", "VASPMATE --dbs vasp.db -table [table-name1] -insert [keyword] [values]"},
    {"--dbs", "VASPMATE --dbs vasp.db -table [table-name1] -delete [condition]"},
    {"--dbs", "VASPMATE --dbs vasp.db -plus vasp1.db"},
    {"--dbs", "VASPMATE --dbs vasp.db -table [table-name] -merge [table-name1 table-name2]"},
    {"--db2s", "VASPMATE --db2s name.db -table [table_name] (-a/-b) (-inc/include file_name)"},
};

//https://blog.csdn.net/tissar/article/details/86475879
class levenshtein
{
public:
    static int compare( const std::string& s1, const std::string& s2 )
    {
        // create two work vectors of integer distances
        const int m = s1.size();
        const int n = s2.size();
        int* v1 = new int[n+1];
        int* v2 = new int[n+1];
        
        // initialize v1 (the previous row of distances)
        // this row is A[0][i]: edit distance for an empty s1
        // the distance is just the number of characters to delete from s2
        for( int i = 0; i <= n; ++i ) {
            v1[i] = i;
        }
        
        // calculate v2 (current row distances) from the previous row v1
        for( int i = 0; i < m; ++i ) {
            // first element of v2 is A[i+1][0]
            // edit distance is delete (i+1) chars from s to match empty s2
            v2[0] = i+1;
        
            // use formula to fill in the rest of the row
            for( int j = 0; j < n; ++j ) {
                // calculating costs for A[i+1][j+1]
                int deletionCost = v1[j+1] + 1;
                int insertionCost = v2[j] + 1;
                int substitutionCost = v1[j];
                if( s1[i] != s2[j] ) {
                    ++ substitutionCost;
                }
                v2[j+1] = min3( deletionCost, insertionCost, substitutionCost );
            }
            // copy v2 (current row) to v1 (previous row) for next iteration
            swap( v1, v2 );
        }
        
        // after the last swap, the results of v2 are now in v1
        int retval = v1[n];
        delete []v1;
        delete []v2;
        return retval;
    }
private:
    static int min3( int a, int b, int c )
    {
        if ( a < b ) {
            return a < c ? a : c ;
        }
        else {
            return b < c ? b : c ;
        }
    }
    static void swap( int*& a, int*& b )
    {
        int* tmp = a;
        a = b;
        b = tmp;
    }
};

void vaspmate_help(int flag_help, int argc, char* argv[])
{
    string check_command = string(argv[1]);
    if (flag_help == 1) //VASPMATE did not recognize your command
    {
        bool find_similar = false;
        cout << "VASPMATE did not recognize your command, you may have meant to input: " << endl << endl;
        for (const auto& pair : help_command) 
        {
            int d = levenshtein::compare( pair.first, check_command );
            if (d <= 1) 
            {
                cout << pair.second << endl;
                find_similar = true;
            }
        }
        if (find_similar == false)
            cout << "Sorry, no similar command was found in VASPMATE!" << endl;
    }
    else if (flag_help == -1)
    {
        cout << "VASPMATE suspects that the format of the command line you entered may be incorrect, please refer to the proper format and re-enter!" << endl << endl;
        for (const auto& pair : help_command) 
        {
            if (!pair.first.empty() && pair.first == check_command) 
            {
                cout << pair.second << endl;
            }
        }
    }
}