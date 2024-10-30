#pragma once
#define _POSIX_C_SOURCE 2
#ifndef _VCLEAN_
#define _VCLEAN_

#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<string>
#include<cstring>
#include<dirent.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>

void rm_file(const std::vector<std::string> &file_vector);

void rm_other_file();

int clean_operat(int argc, char* argv[]);

#endif