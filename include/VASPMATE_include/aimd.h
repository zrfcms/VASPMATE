#pragma once
#ifndef _AIMD_
#define _AIMD_

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<map>
#include<string>
#include<string.h>
#include"tools.h"
#include"read_write.h"
using namespace std;
namespace aimd
{
    int operator_aimd(int argc, char* argv[]);

    void write_nve(int steps);
    void write_nvt(double temp, int steps);
    void write_npt(double temp, double press, int steps);

    void get_et(FILE* fp);
    void get_mt(FILE* fp);

    void get_postep(FILE* fp, vector<int> time_step);
}

#endif