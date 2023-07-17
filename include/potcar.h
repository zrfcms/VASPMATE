#pragma once
#ifndef _POTCAR_
#define _POTCAR_

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<string>
#include<unistd.h>
#include<vector>
#include<algorithm>
#include<sys/uio.h>
#include<read_write.h>
#include<sys/types.h>
#include<pwd.h>
using namespace std;
void pot_merge(const char file[], char mode[], vector<string> label);

void pot_merge_element(vector<string> element, char mode[], vector<string> label);

void check();

#endif