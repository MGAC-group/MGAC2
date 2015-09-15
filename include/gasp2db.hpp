#pragma once
#include "gasp2common.hpp"
#include "gasp2struct.hpp"
#include "gasp2param.hpp"
#include "gasp2pop.hpp"
#include "sqlite3.h"

using namespace std;

class GASP2db {
public:
	GASP2db(string name);

	int init(); //sets up initial tables
	int init(string name);
	int shutdown();

	//takes a population and updates the database table
	//if a structure does not exist it is added
	//if a structure does exist, it is updated
	//this necessarily overwrites structures with the same UID
	//since UID collision is improbable, this will not contribute
	//to a meaningful loss of structures
	int update(GASP2pop pop);

	GASP2pop getAll();

	GASP2pop getGroup(int index);
	GASP2pop getCluster(int index);
	GASP2pop getIndv(UUID u);

private:
	string path;
	sqlite3 *dbconn;



};

