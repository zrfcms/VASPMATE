#pragma once
#ifndef _INPSTD_
#define _INPSTD_

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <map>
#include <algorithm>
using std::cout;
using std::endl;
using std::map;
using std::vector;
using std::string;

namespace inpstd
{
	class incarstd
	{
	public:
        incarstd(const char file[] = "INCAR");
        void std_incar_Content(const vector<string> &incar_con);
        void toget_chunk(const vector<string> &incar_con);
        void write_std_incar(const char file[] = "STDINP");
	private:
        vector<string> incar_Content;
        vector<string> incar_Content_std;
        map<string, vector<string>> chunk;
        vector<string> chunk_order = {
        "Initialization", "File input/output", "Structure Relaxation", "Electronic Self-consistency", "Magnetism", "Optics", 
        "Band Structure", "DOS", "Hybrid Functionals", "VDW", "LDA+U", "Parallelization", "Molecular Dynamics", "Other"};
	};
};

#endif