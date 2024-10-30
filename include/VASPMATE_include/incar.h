#pragma once
#ifndef _INCAR_
#define _INCAR_

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<map>
#include"read_write.h"
std::string write_INCAR(int argc, char* argv[]);

void INCAR_fix(char name[], std::vector<std::string> value);

void INCAR_replace(int argc, char* argv[]);

void INCAR_delete(int argc, char* argv[]);

void INCAR_remove(int argc, char* argv[]);

void INCAR_append(int argc, char* argv[]);

void pcd_model(int argc, char* argv[]);

struct LDAU_DATA
{
	std::string L;
	std::string U;
	std::string J;
};

class LDAU
{
public:
	LDAU(const char table[] = nullptr, const char file[] = "INPOS");
	void AddLDAU_default(const char file[] = "INPOS");
private:
	std::vector<LDAU_DATA> data;
	int LMAXMIX;
	int LDAUTYPE;
};

class Mag
{
public:
	Mag(bool noncol = 0, const char table[] = nullptr);
	void AddMag_default(bool noncol = 0);
private:
	std::vector<string> data;
	bool noncolliner;
};

int addIVDWdefault(const char file[] = "INPOS", string val = "12");
void addIVDW_table(int ivdw_numb = 0);
#endif