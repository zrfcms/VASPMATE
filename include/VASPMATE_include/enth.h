#pragma once
#ifndef _ENTH_
#define _ENTH_

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<map>
#include<string>
#include<string.h>
#include"tools.h"
#include"read_write.h"

using namespace std;
namespace enth
{
    int operator_enth(int argc, char* argv[]);
    void micid_convexhull(string path_, vector<double> ref_en);
    void Formation_Enthalpy(vector<double> ref_en);
}

#endif