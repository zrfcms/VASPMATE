#include"../../include/VASPMATE_include/sqlitedb.h"
#include"../../include/VASPMATE_include/mapdata.h"

const char VMdb::db_table_head[] = "CREATE TABLE ";
const char VMdb::db_table_name[] = "ID INTEGER PRIMARY KEY AUTOINCREMENT," \
    "PATH VARCHAR DEFAULT NULL," \
    "NAME VARCHAR," \
    "STRUCT TEXT DEFAULT NULL," \
    "Classification INTEGER DEFAULT NULL," \
    "Crystal_system INTEGER DEFAULT NULL," \
    "Pearson_symbol VARCHAR DEFAULT NULL," \
    "International_symbol VARCHAR DEFAULT NULL," \
    "Band_gap REAL DEFAULT NULL," \
    "Formation_energy REAL DEFAULT NULL," \
    "Elastic_stability_conditions INTEGER DEFAULT NULL," \
    "Total_magnetization REAL DEFAULT NULL," \
    "Lattice_a REAL DEFAULT NULL," \
    "Lattice_b REAL DEFAULT NULL," \
    "Lattice_c REAL DEFAULT NULL," \
    "Lattice_alpha REAL DEFAULT NULL," \
    "Lattice_beta REAL DEFAULT NULL," \
    "Lattice_gamma REAL DEFAULT NULL," \
    "Lattice_volume REAL DEFAULT NULL," \
    "System_size INTEGER DEFAULT NULL," \
    "Hall_symbol VARCHAR DEFAULT NULL," \
    "International_number INTEGER DEFAULT NULL," \
    "Pointgroup_symbol VARCHAR DEFAULT NULL," \
    "Bandgap_type INTEGER DEFAULT NULL," \
    "Bond_type INTEGER DEFAULT NULL," \
    "Fermi_energy REAL DEFAULT NULL," \
    "VBM_location REAL DEFAULT NULL," \
    "CBM_location REAL DEFAULT NULL," \
    "Band_structure TEXT DEFAULT NULL," \
    "Density_of_state TEXT DEFAULT NULL," \
    "Thermodynamic_Stability INTEGER DEFAULT NULL," \
    "Total_energy REAL DEFAULT NULL," \
    "Cohesive_energy REAL DEFAULT NULL," \
    "Stiffness_tensor TEXT DEFAULT NULL," \
    "Compliance_tensor TEXT DEFAULT NULL," \
    "Voigt_Youngs_modulus REAL DEFAULT NULL," \
    "Voigt_shear_modulus REAL DEFAULT NULL," \
    "Voigt_bulk_modulus REAL DEFAULT NULL," \
    "Voigt_Poisson_ratio REAL DEFAULT NULL," \
    "Reuss_Youngs_modulus REAL DEFAULT NULL," \
    "Reuss_shear_modulus REAL DEFAULT NULL," \
    "Reuss_bulk_modulus REAL DEFAULT NULL," \
    "Reuss_Poisson_ratio REAL DEFAULT NULL," \
    "Hill_Youngs_modulus REAL DEFAULT NULL," \
    "Hill_shear_modulus REAL DEFAULT NULL," \
    "Hill_bulk_modulus REAL DEFAULT NULL," \
    "Hill_Poisson_ratio REAL DEFAULT NULL," \
    "Pugh_ratio REAL DEFAULT NULL," \
    "Cauchy_pressure REAL DEFAULT NULL," \
    "Chung_Buessem_anisotropy_index REAL DEFAULT NULL," \
    "Universal_elastic_anisotropy_index REAL DEFAULT NULL," \
    "Magnetism_ordering INTEGER DEFAULT NULL," \
    "Atomic_magnetic_moment TEXT DEFAULT NULL," \
    "Potential_type VARCHAR DEFAULT NULL," \
    "Functional_type VARCHAR DEFAULT NULL," \
    "Precise VARCHAR DEFAULT NULL," \
    "Energy_cutoff REAL DEFAULT NULL," \
    "Minimization_algorithm VARCHAR DEFAULT NULL," \
    "Electronic_convergence VARCHAR DEFAULT NULL," \
    "Self_consistent VARCHAR DEFAULT NULL," \
    "Intergration_scheme VARCHAR DEFAULT NULL," \
    "Spin_polarization VARCHAR DEFAULT NULL," \
    "Spin_orbit_coupling VARCHAR DEFAULT NULL," \
    "Relaxation VARCHAR DEFAULT NULL," \
    "Pullay_stress REAL DEFAULT NULL," \
    "Ionic_update VARCHAR DEFAULT NULL," \
    "Ionic_convergence VARCHAR DEFAULT NULL," \
    "MetaGGA_type VARCHAR DEFAULT NULL," \
    "Hybrid_type VARCHAR DEFAULT NULL," \
    "DFT_U_type VARCHAR DEFAULT NULL," \
    "VDW_D_type VARCHAR DEFAULT NULL," \
    "Solvation_model VARCHAR DEFAULT NULL," \
    "Dipole_correction VARCHAR DEFAULT NULL," \
    "Energy_in_OUTCAR TEXT DEFAULT NULL," \
    "Stress_in_OUTCAR TEXT DEFAULT NULL," \
    "Force_in_OUTCAR TEXT DEFAULT NULL," \
    "Charge_in_OUTCAR TEXT DEFAULT NULL," \
    "Magnitization_in_OUTCAR TEXT DEFAULT NULL," \
    "Elasticity_in_OUTCAR TEXT DEFAULT NULL," \
    "Author VARCHAR DEFAULT NULL," \
    "Affiliation VARCHAR DEFAULT NULL," \
    "Email VARCHAR DEFAULT NULL," \
    "Date VARCHAR DEFAULT NULL," \
    "Source VARCHAR DEFAULT NULL," \
    "Source_ID VARCHAR DEFAULT NULL," \
    "Reference TEXT DEFAULT NULL";

const char VMdb::db_last_id[] = "SELECT last_insert_rowid()";

int VMdb::operator_db(int argc, char* argv[])
{
    int _table; int p_vec;
	int check_table = check_para(argc, argv, "-table", 1 , &_table);
    vector<const char*> par_vec = {"-tables", "-schema", "-create", "-alter", "-set", "-drop", "-insert", "-delete", "-combine", "-update", "-plus", "-select"};
	int check_vec = check_mulpara(argc, argv, par_vec, 0 , &p_vec);
    sqlite3 *db; int db_return = 0;
    if(argc == 2)
        db_return = open_database(&db, "vasp.db");
    else
        db_return = open_database(&db, argv[2]);
    if(check_vec != 2)
    {
        if(!strcmp("-tables", argv[p_vec]))
            db_return = query_table(&db);
        else if(!strcmp("-schema", argv[p_vec]))
        {
            if(p_vec == argc - 1)
                db_return = query_table(&db, "vaspinfo");
            else
                db_return = query_table(&db, argv[p_vec + 1]);
        }
        else if(!strcmp("-create", argv[p_vec]))
        {
            if(p_vec == argc -1)
                db_return = create_table(&db, argv[_table + 1], db_table_name);
            else
            {
                int totalSize = 1;
                for (int i = p_vec + 1; i < argc; i++)
                    totalSize = totalSize + strlen(argv[i]) + 1;
                char *creat_type = (char *)malloc(totalSize * sizeof(char));
                creat_type[0] = '\0';
                // CREATE TABLE table_name(column1 datatype,column2 datatype,PRIMARY KEY( one or more columns ));
                for (int i = p_vec + 1; i < argc; i++)
                {
                    strcat(creat_type, argv[i]);
                    strcat(creat_type, " ");
                }
                db_return = create_table(&db, argv[_table + 1], creat_type);
                free(creat_type);
            }
        }
        else if(!strcmp("-alter", argv[p_vec]))
            db_return = alter_table(&db, argv[_table + 1], argv[p_vec + 1]);
        else if(!strcmp("-update", argv[p_vec]))
            db_return = update_table(&db, argv[_table + 1], argv[p_vec + 1], argv[p_vec + 2]);
        else if(!strcmp("-set", argv[p_vec]))
        {
            int totalSize = 1;
            for (int i = p_vec + 1; i < argc; i++)
                totalSize = totalSize + strlen(argv[i]) + 1;
            char *creat_type = (char *)malloc(totalSize * sizeof(char));
            creat_type[0] = '\0';
            for (int i = p_vec + 1; i < argc; i++)
            {
                strcat(creat_type, argv[i]);
                strcat(creat_type, " ");
            }
            db_return = alter_column(&db, argv[_table + 1], argv[p_vec + 1], creat_type);
            free(creat_type);
        }
        else if(!strcmp("-insert", argv[p_vec]))
        {
            string stand_attr = standard_str(argv[p_vec + 1]);
            vector<string> stand_msg = standard_vec(argv[p_vec + 2]);
            db_return = insert_data(&db, argv[_table + 1], stand_attr.c_str(), stand_msg);
        }
        else if(!strcmp("-drop", argv[p_vec]))
            db_return = delete_table(&db, argv[_table + 1]);
        else if(!strcmp("-combine", argv[p_vec]))
            db_return = combine_dbtable(&db, argv[_table + 1], argv[p_vec + 1]);
        else if(!strcmp("-plus", argv[p_vec]))
            db_return = plus_dbtable(&db, argv[p_vec + 1]);
        else if(!strcmp("-delete", argv[p_vec]))
        {
            int totalSize = 1;
            for (int i = p_vec + 1; i < argc; i++)
                totalSize = totalSize + strlen(argv[i]) + 1;
            char *creat_type = (char *)malloc(totalSize * sizeof(char));
            creat_type[0] = '\0';
            // CREATE TABLE table_name(column1 datatype,column2 datatype,PRIMARY KEY( one or more columns ));
            for (int i = p_vec + 1; i < argc; i++)
            {
                strcat(creat_type, argv[i]);
                    strcat(creat_type, " ");
            } 
            db_return = delete_data(&db, argv[_table + 1], creat_type);
            free(creat_type);
        }
        else if(!strcmp("-select", argv[p_vec]))
        {
            int totalSize = 1;
            for (int i = p_vec + 1; i < argc; i++)
                totalSize = totalSize + strlen(argv[i]) + 1;
            char *creat_type = (char *)malloc(totalSize * sizeof(char));
            creat_type[0] = '\0';
            // CREATE TABLE table_name(column1 datatype,column2 datatype,PRIMARY KEY( one or more columns ));
            for (int i = p_vec + 1; i < argc; i++)
            {
                strcat(creat_type, argv[i]);
                strcat(creat_type, " ");
            }      
            db_return = find_data(&db, argv[_table + 1], creat_type);
        }
        else
        {
            sqlite3_close(db);
            return -1;
        }
    }
    sqlite3_close(db);
    return 0;
}

int VMdb::open_database(sqlite3 **db, const char *database_name)
{
    int len;
    len = sqlite3_open(database_name, db);
    if(len)
    {
        printf("Open database name %s failure.\n", database_name);
        return -1;
    }
    printf("Open a sqlite3 database name %s successfully!\n", database_name);
    return 0;
}

int VMdb::create_table(sqlite3 **db, const char *table_name, const char *table_attribute)
{
    char *zErrMsg=NULL;
    if (!issafe_input(table_name) || !issafe_input(table_attribute))
    {
        printf("Warning: Dangerous input! %s or %s\n", table_name, table_attribute);
        return -1;
    }
    int totalSize = strlen(db_table_head) + strlen(table_name) + strlen(table_attribute) + strlen("''()") + 1;
    char *table_command = (char *)malloc(totalSize * sizeof(char)); 
    // CREATE TABLE table_name(column1 datatype,column2 datatype,PRIMARY KEY( one or more columns ));
    if (table_command)
        snprintf(table_command, totalSize, "%s'%s'(%s)", db_table_head, table_name, table_attribute);
    //cout << table_command <<endl;
    if(sqlite3_exec(*db, table_command, NULL, NULL, &zErrMsg) != SQLITE_OK)
        printf("SQL warning: %s\n", zErrMsg);
    else
        printf("Create table %s successfully\n", table_name);
    sqlite3_free(zErrMsg);
    free(table_command);
    return 0;
}

int VMdb::insert_data(sqlite3 **db, const char *table_name, const char *attr, vector<string> msg)
{
    sqlite3_stmt *stmt;
    if (!issafe_input(table_name))
    {
        printf("Warning: Dangerous input! %s\n", table_name);
        return -1;
    }
    string question_mark = "?";
    for(int i = 0; i < msg.size() - 1; i++)
        question_mark += ",?";
    string sql = "INSERT INTO " + string(table_name) + " (" + string(attr) + ") VALUES (" + question_mark+ ");";
    if (sqlite3_prepare_v2(*db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK)
    {
        printf("Failed to prepare statement: %s\n", sqlite3_errmsg(*db));
        return -1;
    }
    for (int i = 0; i < msg.size(); ++i)
        sqlite3_bind_text(stmt, i + 1, msg[i].c_str(), -1, SQLITE_TRANSIENT);
    if (sqlite3_step(stmt) != SQLITE_DONE)
        printf("SQL error: %s\n", sqlite3_errmsg(*db));
    else
    {
        /*printf("Insert ");
        for (int i = 0; i < msg.size(); ++i)
            printf("%s ", msg[i].c_str());
        printf("of %s to table %s successfully\n", attr, table_name);*/
        printf("Insert data successfully\n");
    }
    sqlite3_finalize(stmt);
    return 0;
}

int VMdb::insert_data(sqlite3 **db, const char *table_name, string attr, string msg)
{
    sqlite3_stmt *stmt;
    if (!issafe_input(table_name))
    {
        printf("Warning: Dangerous input! %s\n", table_name);
        return -1;
    }
    string sql = "INSERT INTO " + string(table_name) + " (" + attr + ") VALUES (?);";
    if (sqlite3_prepare_v2(*db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK)
    {
        printf("Failed to prepare statement: %s\n", sqlite3_errmsg(*db));
        return -1;
    }
    sqlite3_bind_text(stmt, 1, msg.c_str(), -1, SQLITE_TRANSIENT);
    sqlite3_finalize(stmt);
    return 0;
}

int VMdb::query_table(sqlite3 **db)
{ 
    char *zErrMsg=NULL;
    int totalSize = strlen("SELECT tbl_name FROM sqlite_master WHERE type = 'table';") + 1;
    char *query_command = (char *)malloc(totalSize * sizeof(char));
    if (query_command)
        snprintf(query_command, totalSize, "SELECT tbl_name FROM sqlite_master WHERE type = 'table';");
    if(sqlite3_exec(*db, query_command, VMdb::write_callback, NULL, &zErrMsg) != SQLITE_OK)
    {
        printf("SQL error: %s\n", zErrMsg);
    }
    free(query_command);
    sqlite3_free(zErrMsg);
    return 0;
}

int VMdb::query_table(sqlite3 **db, const char *table_name)
{
    char *zErrMsg=NULL;
    if (!issafe_input(table_name))
    {
        printf("Warning: Dangerous input! %s\n", table_name);
        return -1;
    }
    int totalSize = strlen("SELECT sql FROM sqlite_master WHERE type='table' AND name='';") + strlen(table_name) + 1;
    char *query_command = (char *)malloc(totalSize * sizeof(char));
    if (query_command)
        snprintf(query_command, totalSize, "SELECT sql FROM sqlite_master WHERE type='table' AND name='%s';", table_name);
    if(sqlite3_exec(*db, query_command, VMdb::write_callback, NULL, &zErrMsg) != SQLITE_OK)
    {
        printf("SQL error: %s\n", zErrMsg);
    }
    free(query_command);
    sqlite3_free(zErrMsg);
    return 0;
}

int VMdb::delete_table(sqlite3 **db, const char *table_name)
{
    char *zErrMsg=NULL;
    if (!issafe_input(table_name))
    {
        printf("Warning: Dangerous input! %s\n", table_name);
        return -1;
    }
    int totalSize = strlen("DROP TABLE '';") + strlen(table_name) + 1;
    char *delete_command = (char *)malloc(totalSize * sizeof(char));
    if (delete_command)
        snprintf(delete_command, totalSize, "DROP TABLE '%s';", table_name);
    //cout << delete_command << endl;
    if(sqlite3_exec(*db, delete_command, NULL, NULL, &zErrMsg) != SQLITE_OK)
    {
        printf("SQL error: %s\n", zErrMsg);
    }
    else
        printf("Delete table %s successfully!\n", table_name);
    sqlite3_free(zErrMsg);
    return 0;
}

int VMdb::delete_data(sqlite3 **db, const char *table_name, const char *condition)
{
    sqlite3_stmt *stmt;
    if (!issafe_input(table_name))
    {
        printf("Warning: Dangerous input! %s\n", table_name);
        return -1;
    }
    string sql = "DELETE FROM " + string(table_name) + " WHERE " + string(condition) + ";";
    if (sqlite3_prepare_v2(*db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK) {
        printf("Failed to prepare statement: %s\n", sqlite3_errmsg(*db));
        return -1;
    }
    int step_result = sqlite3_step(stmt);
    if (step_result != SQLITE_DONE)
    {
        printf("SQL error: %s\n", sqlite3_errmsg(*db));
        sqlite3_finalize(stmt);
        return -1;
    } 
    else 
    {
        int affected_rows = sqlite3_changes(*db);
        if (affected_rows == 0)
            printf("No records match the specified condition [%s] in table %s!\n", condition, table_name);
        else
            printf("Deleted %d records from table %s where condition [%s]!\n", affected_rows, table_name, condition);
    }
    return 0;
}

int VMdb::find_data(sqlite3 **db, const char *table_name, const char *condition)
{
    sqlite3_stmt *stmt;
    if (!issafe_input(table_name))
    {
        printf("Warning: Dangerous input! %s\n", table_name);
        return -1;
    }
    string sql = "SELECT * FROM " + string(table_name) + " WHERE " + string(condition) + ";";
    if (sqlite3_prepare_v2(*db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK) {
        printf("Failed to prepare statement: %s\n", sqlite3_errmsg(*db));
        return -1;
    }
    int step_result = sqlite3_step(stmt);
    if (step_result == SQLITE_ROW) 
    {
        do {
            int num_cols = sqlite3_column_count(stmt);
            for (int i = 0; i < num_cols; ++i) {
                const char* col_name = sqlite3_column_name(stmt, i);
                const char* col_data = (const char*)sqlite3_column_text(stmt, i);
                if (col_data) {
                    printf("%s: %s\n", col_name, col_data);
                } else {
                    printf("%s: NULL\n", col_name);
                }
            }
            printf("\n");
        } while ((step_result = sqlite3_step(stmt)) == SQLITE_ROW);
    } else if (step_result == SQLITE_DONE) {
        printf("Warning: No rows matched the condition!\n");
    } else {
        printf("SQL error: %s\n", sqlite3_errmsg(*db));
    }
    sqlite3_finalize(stmt);
    return 0;
}

int VMdb::alter_table(sqlite3 **db, const char *old_table, const char *new_table)
{
    char *zErrMsg = NULL;
    if (!issafe_input(old_table) || !issafe_input(new_table))
    {
        printf("Warning: Dangerous input! %s or %s\n", old_table, new_table);
        return -1;
    }
    string sql = "ALTER TABLE " + string(old_table) + " RENAME TO '" + string(new_table) + "';";
    if (sqlite3_exec(*db, sql.c_str(), NULL, NULL, &zErrMsg) != SQLITE_OK) 
    {
        printf("SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return -1;
    }
    printf("Change table name from %s to %s successfully!\n", old_table, new_table);
    return 0;
}

int VMdb::alter_column(sqlite3 **db, const char *table_name, const char *column_name, const char *column_definition)
{
    char *zErrMsg = NULL;
    if (!issafe_input(column_name) || !issafe_input(column_definition))
    {
        printf("Warning: Dangerous input! %s or %s\n", column_name, column_definition);
        return -1;
    }
    std::string sql = "ALTER TABLE " + std::string(table_name) + " ADD COLUMN " + std::string(column_name) + " " + std::string(column_definition) + ";";
    if (sqlite3_exec(*db, sql.c_str(), NULL, NULL, &zErrMsg) != SQLITE_OK) 
    {
        printf("SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return -1;
    }
    printf("Add column %s(%s) in table %s successfully!\n", column_name, column_definition, table_name);
    return 0;
}

int VMdb::update_table(sqlite3 **db, const char *table_name, const char *set_command, const char *condition)
{
    sqlite3_stmt *stmt;
    string sql = "UPDATE " + string(table_name) + " SET " + string(set_command) + " WHERE " + string(condition) + ";";
    if (sqlite3_prepare_v2(*db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK) {
        printf("Failed to prepare statement: %s\n", sqlite3_errmsg(*db));
        return -1;
    }
    if (sqlite3_step(stmt) != SQLITE_DONE)
        printf("SQL error: %s\n", sqlite3_errmsg(*db));
    else
        printf("Update table successfully!\n");
    sqlite3_finalize(stmt);
    return 0; 
}

int VMdb::get_last_id(sqlite3 **db)
{
    sqlite3_stmt *stmt;
    int lastId = -1;

    if (sqlite3_prepare_v2(*db, db_last_id, -1, &stmt, 0) == SQLITE_OK)
        if (sqlite3_step(stmt) == SQLITE_ROW)
            lastId = sqlite3_column_int(stmt, 0);
    sqlite3_finalize(stmt);
    return lastId;
}

int VMdb::combine_dbtable(sqlite3 **db, const char *table_name, const char *db2_path)
{
    char *zErrMsg = NULL;
    int rc;
    if (!issafe_input(table_name))
    {
        printf("Warning: Dangerous input! %s\n", table_name);
        return -1;
    }
    string selectQuery = getno_keyColumns(db, string(table_name));
    string attachCmd = "ATTACH DATABASE ? AS db2;";
    sqlite3_stmt* attachStmt;
    if (sqlite3_prepare_v2(*db, attachCmd.c_str(), -1, &attachStmt, NULL) == SQLITE_OK) 
    {
        sqlite3_bind_text(attachStmt, 1, db2_path, -1, SQLITE_STATIC);
        rc = sqlite3_step(attachStmt);
        sqlite3_finalize(attachStmt);
        if (rc != SQLITE_DONE) 
        {
            printf("Failed to attach database %s\n", db2_path);
            return -1;
        }
    } 
    else 
    {
        printf("Failed to prepare attach database statement\n");
        return -1;
    }
    string mergeCmd = "INSERT INTO " + string(table_name) + "(" + selectQuery + ") " +
                        " SELECT " + selectQuery + " FROM db2." + string(table_name) + ";";
    rc = sqlite3_exec(*db, mergeCmd.c_str(), NULL, NULL, &zErrMsg);
    if (rc != SQLITE_OK) 
    {
        printf("SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return -1;
    }
    string detachCmd = "DETACH DATABASE db2;";
    sqlite3_exec(*db, detachCmd.c_str(), NULL, NULL, &zErrMsg);
    printf("Successfully merged table %s from %s.\n", table_name, db2_path);
    return 0;
}

int VMdb::plus_dbtable(sqlite3 **db, const char *db2_path)
{
    sqlite3* db2;
    char *zErrMsg = NULL;
    int rc = sqlite3_open(db2_path, &db2);
    if (rc != SQLITE_OK) 
    {
        printf("Cannot open database: %s\n", sqlite3_errmsg(db2));
        return -1;
    }
    string attachCmd = "ATTACH DATABASE ? AS db2;";
    sqlite3_stmt* attachStmt;
    if (sqlite3_prepare_v2(*db, attachCmd.c_str(), -1, &attachStmt, NULL) == SQLITE_OK) 
    {
        sqlite3_bind_text(attachStmt, 1, db2_path, -1, SQLITE_STATIC);
        rc = sqlite3_step(attachStmt);
        sqlite3_finalize(attachStmt);
        if (rc != SQLITE_DONE) 
        {
            printf("Failed to attach database %s\n", db2_path);
            return -1;
        }
    }
    else
    {
        printf("Failed to prepare attach database statement\n");
        return -1;
    }
    sqlite3_stmt* stmt;
    string getTableNamesSQL = "SELECT name FROM sqlite_master WHERE type='table';";
    if (sqlite3_prepare_v2(db2, getTableNamesSQL.c_str(), -1, &stmt, nullptr) == SQLITE_OK)
    {
        while (sqlite3_step(stmt) == SQLITE_ROW)
        {
            string tableName = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
            if (tableName == "sqlite_sequence")
                continue;
            string checkTableSQL = "SELECT name FROM sqlite_master WHERE type='table' AND name='" + tableName + "';";
            sqlite3_stmt* checkStmt;
            if (sqlite3_prepare_v2(*db, checkTableSQL.c_str(), -1, &checkStmt, nullptr) == SQLITE_OK)
            {
                if (sqlite3_step(checkStmt) == SQLITE_ROW)
                {
                    string selectQuery = getno_keyColumns(&db2, string(tableName));
                    string mergeDataSQL = "INSERT INTO " + string(tableName) + "(" + selectQuery + ") " +
                        " SELECT " + selectQuery + " FROM db2." + string(tableName) + ";";
                    rc = sqlite3_exec(*db, mergeDataSQL.c_str(), NULL, NULL, &zErrMsg);
                    if (rc != SQLITE_OK)
                    {
                        printf("Failed to merge table: %s, error: %s\n", tableName.c_str(), zErrMsg);
                        sqlite3_free(zErrMsg);
                    }
                }
                else
                {
                    string mergeDataSQL = "CREATE TABLE IF NOT EXISTS " + tableName +
                                        " AS SELECT * FROM db2." + tableName + ";";
                    rc = sqlite3_exec(*db, mergeDataSQL.c_str(), NULL, NULL, &zErrMsg);
                    if (rc != SQLITE_OK)
                    {
                        printf("Failed to merge table: %s, error: %s\n", tableName.c_str(), zErrMsg);
                        sqlite3_free(zErrMsg);
                    }
                }
                sqlite3_finalize(checkStmt);
            }
        }
        sqlite3_finalize(stmt);
    } 
    else 
    {
        printf("Failed to prepare table listing statement.\n");
        sqlite3_close(db2);
        return -1;
    }
    string detachCmd = "DETACH DATABASE db2;";
    sqlite3_exec(*db, detachCmd.c_str(), NULL, NULL, &zErrMsg);
    sqlite3_close(db2);
    printf("Successfully merged %s.\n", db2_path);
    return 0;
}

int VMdb::insert_data_no_print(sqlite3 **db, const char *table_name, const char *attr, vector<string> msg)
{
    sqlite3_stmt *stmt;
    if (!issafe_input(table_name))
    {
        return -1;
    }
    string question_mark = "?";
    for(int i = 0; i < msg.size() - 1; i++)
        question_mark += ",?";
    string sql = "INSERT INTO " + string(table_name) + " (" + string(attr) + ") VALUES (" + question_mark+ ");";
    if (sqlite3_prepare_v2(*db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK)
    {
        return -1;
    }
    for (int i = 0; i < msg.size(); ++i)
        sqlite3_bind_text(stmt, i + 1, msg[i].c_str(), -1, SQLITE_TRANSIENT);
    if (sqlite3_step(stmt) != SQLITE_DONE)
        return -1;
    sqlite3_finalize(stmt);
    return 0;
}

int VMdb::alter_column_no_print(sqlite3 **db, const char *table_name, const char *column_name, const char *column_definition)
{
    char *zErrMsg = NULL;
    if (!issafe_input(column_name) || !issafe_input(column_definition))
    {
        return -1;
    }
    std::string sql = "ALTER TABLE " + std::string(table_name) + " ADD COLUMN " + std::string(column_name) + " " + std::string(column_definition) + ";";
    if (sqlite3_exec(*db, sql.c_str(), NULL, NULL, &zErrMsg) != SQLITE_OK) 
    {
        sqlite3_free(zErrMsg);
        return -1;
    }
    return 0;
}

int VMdb::update_table_no_print(sqlite3 **db, const char *table_name, const char *set_command, const char *condition)
{
    sqlite3_stmt *stmt;
    string sql = "UPDATE " + string(table_name) + " SET " + string(set_command) + " WHERE " + string(condition) + ";";
    if (sqlite3_prepare_v2(*db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK)
        return -1;
    if (sqlite3_step(stmt) != SQLITE_DONE)
        return -1;
    sqlite3_finalize(stmt);
    return 0; 
}

int VMdb::write_callback(void* data, int argc, char** argv, char** azColName)
{
    for (int i = 0; i < argc; i++) 
    {
        cout << azColName[i] << ": " << argv[i] << endl;
    }
    return 0;
}

string VMdb::standard_str(const char* input)
{
    string inputStr(input);  
    stringstream ss(inputStr);  
    string item;
    vector<string> parts;   
    while (std::getline(ss, item, ','))
        parts.push_back("'" + item + "'");  
    std::ostringstream result;  
    for (size_t i = 0; i < parts.size(); ++i) 
    {  
        if (i != 0)
            result << ",";
        result << parts[i];
    }
    return result.str();
}

vector<string> VMdb::standard_vec(const char* input)
{
    string inputStr(input);  
    stringstream ss(inputStr);  
    string item;
    vector<string> parts;   
    while (std::getline(ss, item, ','))
        parts.push_back(item);
    return parts;
}

string VMdb::getno_keyColumns(sqlite3 **db, const string& tableName)
{
    vector<string> columns;
    string query = "PRAGMA table_info(" + tableName + ");";
    sqlite3_stmt* stmt;
    int rc = sqlite3_prepare_v2(*db, query.c_str(), -1, &stmt, nullptr);
    while ((rc = sqlite3_step(stmt)) == SQLITE_ROW)
    {
        string columnName(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1)));
        int pk = sqlite3_column_int(stmt, 5); // Primary Key flag
        if (pk == 0)
            columns.push_back(columnName);
    }
    sqlite3_finalize(stmt);
    string selectQuery;
    for (size_t i = 0; i < columns.size(); ++i)
    {
        if(i != columns.size() - 1)
            selectQuery +=columns[i]  + ",";
        else
            selectQuery +=columns[i];
    }
    return selectQuery;
}

bool VMdb::issafe_input(const char *str)
{
    while (*str)
    {
        if (!isalnum((unsigned char)*str) && *str != '_' && *str != ' '&& *str != ',')
            return false;
        str++;
    }
    return true;
}
