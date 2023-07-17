#pragma once
#ifndef _NEB_
#define _NEB_

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<sys/stat.h>
#include"read_write.h"
double neb_similarity(const char name1[50], const char name2[50]);
void neb_Points(const char name1[50], const char name2[50], int N);
int neb_data(int N);
int neb_movie(int N);
void neb(int argc, char* argv[]);
#endif