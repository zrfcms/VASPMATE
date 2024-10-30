#pragma once
#ifndef _VM_dbs_
#define _VM_dbs_

#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<sqlite3.h>
#include<sstream>
#include"tools.h"
#include"read_write.h"
using namespace std;

namespace VMdb
{
    extern const char db_table_head[];
    extern const char db_table_name[];

    extern const char db_last_id[];
    
    int operator_db(int argc, char* argv[]);

    int open_database(sqlite3 **db, const char *database_name);
    int create_table(sqlite3 **db, const char *table_name, const char *table_attribute);
    int insert_data(sqlite3 **db, const char *table_name, const char *attr, vector<string> msg);
    int insert_data(sqlite3 **db, const char *table_name, string attr, string msg);

    int query_table(sqlite3 **db);
    int query_table(sqlite3 **db, const char *table_name);

    int alter_table(sqlite3 **db, const char *old_table, const char *new_table);
    int alter_column(sqlite3 **db, const char *table_name, const char *column_name, const char *column_definition);

    int update_table(sqlite3 **db, const char *table_name, const char *set_command, const char *condition);
    int get_last_id(sqlite3 **db);

    int delete_table(sqlite3 **db, const char *table_name);
    int delete_data(sqlite3 **db, const char *table_name, const char *condition);

    int combine_dbtable(sqlite3 **db, const char *table_name, const char *db2);
    int plus_dbtable(sqlite3 **db, const char *database2);

    int find_data(sqlite3 **db, const char *table_name, const char *condition);

    int write_callback(void* data, int argc, char** argv, char** azColName);

    string standard_str(const char* input); //"apple,banana,cherry" -> 'apple','banana','cherry'
    vector<string> standard_vec(const char* input); //"apple,banana,cherry" -> <'apple','banana','cherry'>
    string getno_keyColumns(sqlite3 **db, const string& tableName);

    bool issafe_input(const char *str);

    int alter_column_no_print(sqlite3 **db, const char *table_name, const char *column_name, const char *column_definition);   //for RaW4db
    int insert_data_no_print(sqlite3 **db, const char *table_name, const char *attr, vector<string> msg);
    int update_table_no_print(sqlite3 **db, const char *table_name, const char *set_command, const char *condition);
}

#endif