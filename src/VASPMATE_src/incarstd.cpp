#include"../../include/VASPMATE_include/incarstd.h"
using std::cout;
using std::endl;
using std::vector;
using std::string;
using namespace inpstd;

void incarstd::std_incar_Content(const vector<string> &incar_con)
{
    std::unordered_set<string> keywords;
    for (const auto &line : incar_con) 
    {
        string processedLine = line;
        // Keywords to be ignored at first , delete the line
        if (processedLine.find_first_not_of(" \t") == processedLine.find("#") || processedLine.find_first_not_of(" \t") == processedLine.find("!") 
         || processedLine.find_first_not_of(" \t") == processedLine.find("(")|| processedLine.find("=") == string::npos)
        {
            continue;
        }
        // Keywords to be ignored at mid , delete the comment behind the Keywords
        size_t pos = processedLine.find("#");
        if (pos == string::npos)
            pos = processedLine.find("!");
        if (pos == string::npos)
            pos = processedLine.find("(");
        if (pos != string::npos)
            processedLine = processedLine.substr(0, pos);
        processedLine = processedLine.substr(0, processedLine.find_last_not_of(" \t\n\r") + 1);
        processedLine = processedLine.substr(processedLine.find_first_not_of(" \t"));
        string keyword;
        std::istringstream lineStream(processedLine);
        lineStream >> keyword;
        if (keywords.find(keyword) == keywords.end()) 
        {
            keywords.insert(keyword);
            incar_Content_std.push_back(processedLine);
        }
    }
}

void incarstd::toget_chunk(const vector<string> &incar_con)
{
    for (const auto &line : incar_con)
    {
        std::istringstream lineStream(line);
        string keyword;
        lineStream >> keyword;
        string category;
        if (       keyword == "ISTART"  || keyword == "ICHARG"  || keyword == "ALGO"    || keyword == "TIME"    || keyword == "EDIFF"   || 
                   keyword == "ENCUT"   || keyword == "KSPACING"|| keyword == "KGAMMA"  || keyword == "NCORE"   || keyword == "LREAL"   ||
                   keyword == "PREC"    || keyword == "ADDGRID" || keyword == "SYSTEM"  || keyword == "PRECFOCK"|| keyword == "ISMEAR"  ||
                   keyword == "SIGMA"   ) {
            category = "Initialization";

        } else if (keyword == "LCHARG"  || keyword == "LWAVE"   || keyword == "LVTOT"   || keyword == "LVHAR") {
            category = "File input/output";

        } else if (keyword == "NSW"     || keyword == "IBRION"  || keyword == "ISIF"    || keyword == "EDIFFG"  || keyword == "SYMPREC" ||
                   keyword == "ISYM"    || keyword == "SMASS"   || keyword == "PSTRESS" ) {
            category = "Structure Relaxation";

        } else if (keyword == "NELM"    || keyword == "LORBIT"  || keyword == "NELMIN") {
            category = "Electronic Self-consistency";

        } else if (keyword == "MAGMOM"  || keyword == "ISPIN"   || keyword == "LNONCOLLINEAR"||keyword == "LSORBIT"|| keyword == "SAXIS" ||
                   keyword == "VOSKOWN") {
            category = "Magnetism";

        } else if (keyword == "LOPTICS" || keyword == "CSHIFT"  || keyword == "OMEGAMAX"|| keyword == "OMEGAMIN") {
            category = "Optics";

        } else if (keyword == "NBANDS"  || keyword == "HFALGO"  || keyword == "LPARD"   || keyword == "LSCALAPACK"|| keyword == "LMAXPAW"||
                   keyword == "LASYNC"  || keyword == "NPCG"    || keyword == "LMAXFOCK") {
            category = "Band Structure";

        } else if (keyword == "EMIN"    || keyword == "EMAX"    || keyword == "NEDOS"   || keyword == "DEGAS"   || keyword == "NFREE") {
            category = "DOS";

        } else if (keyword == "HFSCREEN"|| keyword == "LHFCALC" || keyword == "AEXX"    || keyword == "NKREDZ"  || keyword == "NKREDY"  ||
                   keyword == "ENCUTFOCK"|| keyword == "NKRED"  || keyword == "NKREDX"   ) {
            category = "Hybrid Functionals";

        } else if (keyword == "LUSE_VDW"|| keyword == "VDW_ALPHA"|| keyword == "VDW_C6"|| keyword == "LVDWSCS" || keyword == "VDW_RADIUS"||
                   keyword == "VDW_S6"  || keyword == "VDW_SR"   || keyword == "VDW_D" || keyword == "VDW_C6AU"|| keyword == "VDW_ALPHA" ||
                   keyword == "VDW_C6"  || keyword == "VDW_R0AU" || keyword == "VDW_R0"|| keyword == "LVDW_EWALD") {
            category = "VDW";

        } else if (keyword == "LDAUU"   || keyword == "LDAUJ"    || keyword == "LDAUL" || keyword == "LMAXFOCKAE"|| keyword == "LPARDIFF"||
                   keyword == "LMAXMIX"){
            category = "LDA+U";

        } else if (keyword == "LELF"    || keyword == "NPAR"     || keyword == "KPAR") {
            category = "Parallelization";

        } else if (keyword == "POTIM"   || keyword == "IORDER"   || keyword == "NSIM"  || keyword == "LPLANE") {
            category = "Molecular Dynamics";

        } else {
            category = "Other";
        }
        chunk[category].push_back(line);
    }
}
void incarstd::write_std_incar(const char file[])
{
    FILE* fp = fopen(file, "w");
    fprintf(fp, "! Standardize the INCAR using VASPMATE!\n");
    for (const auto &category : chunk_order)
    {
        auto it = chunk.find(category);
        if (it != chunk.end())
        {
            fprintf(fp, "# %s\n", it->first.c_str());
            for (const auto &line : it->second)
            {
                fprintf(fp, "  %s\n", line.c_str());
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
    printf("Written %s file!\n", file);
}
incarstd::incarstd(const char file[])
{
    std::ifstream ifs(file);
    if (!ifs.is_open()) 
    {
        cout << "NO "<< file << " IS NOT EXIST!" << endl;
        return;
    }
    string line;
    while (std::getline(ifs, line)) 
    {
        incar_Content.push_back(line);
    }
    std_incar_Content(incar_Content); //incar_Content_std will be used below
    toget_chunk(incar_Content_std);
    write_std_incar();
}