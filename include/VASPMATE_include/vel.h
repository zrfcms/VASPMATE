#pragma once
#ifndef _VEL_
#define _VEL_

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<cmath>
#include<vector>
#include<map>
#include<random>
#include"tools.h"
#include<string>
#include<string.h>
#include"read_write.h"
using namespace std;
namespace vel
{
    int operator_vel(int argc, char* argv[]);

    void write_vel(double temperature,const char* inpos); 
};
void copyfileContent(const char* inpos,const char* velopos);

#endif
