#pragma once
#ifndef _INCAR_
#define _INCAR_

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<vector>
#include<string>

std::string write_INCAR(int argc, char* argv[]);

void INCAR_fix(char name[], std::vector<std::string> value);

void INCAR_replace(int argc, char* argv[]);

void INCAR_delete(int argc, char* argv[]);

void INCAR_remove(int argc, char* argv[]);

void INCAR_append(int argc, char* argv[]);

void pcd_model(int argc, char* argv[]);

#endif