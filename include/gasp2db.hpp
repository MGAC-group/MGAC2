#pragma once
#include "gasp2common.hpp"
#include "gasp2struct.hpp"
#include "gasp2param.hpp"
#include "gasp2pop.hpp"
#include "sqlite3.h"

using namespace std;


/*
 * A POTENTIAL PROBLEM:
 * sqlite3 doesn't like NFS; there are corruption issues.
 * so we write to a local file and then when we KNOW that
 * the lock mode is in a safe state, we copy to the nfs
 * destination so that data is not lost. this also
 * supports the idea of being able to access output
 * files in the middle of a run, so that two concurrent
 * accesses doesn't ruin the day
 *
 * a reasonable copy mode is the following:
 *  {
 *    std::ifstream  src("from.ogv", std::ios::binary);
 *    std::ofstream  dst("to.ogv",   std::ios::binary);
 *
 *    dst << src.rdbuf();
 *  }
 *
 *
 *  whether or not this is important depends pretty strongly
 *  on how large the DB is, probably. also depends on commit
 *  frequency. we can probably expect problems if the FS goes
 *  south, but in that situation any r/w will fail.
 *
 *
 *
 * Table 1: Structures
 * contains the main structure list in stored form w/blobs
 * organized by UUID
 *
 * NOPE NOPE NOPE
 * TABLE 2 NOT INCLUDED IN ACTUAL SOFTWARE, XML IS JUST STORED INSTEAD
 * Table 2: Storage types
 * contains an xml meta format for the coordinate blobs
 * organized by storage form
 *
 * Table 3: Input files
 * contains the input files; corresponds to input
 * organized by input file number
 *
 *
 *
 */

class GASP2db {
public:
	GASP2db(); //name of DB

	void init(); //sets up initial tables
	void initTable(string name);
	int connect(); //opens the DB
	int disconnect(); //closes the DB
	int load(string name); //loads the table path, then checks to see if it can be opened

	//takes a population and updates the database table
	//if a structure does not exist it is added
	//if a structure does exist, it is updated
	//in an update, the second entry of a structure DB is overwritten
	//the first entry is reserved for initial commit after fitcell
	//is completed, to track differences between fitcell and opt
	//this necessarily overwrites structures with the same UID
	//since UID collision is improbable, this will not contribute
	//to a meaningful loss of structures
	bool create(GASP2pop pop, string table); //creates main record
	bool update(GASP2pop pop, string table); //updates second record

	GASP2pop getAll();

	GASP2pop getSpcGroup(int best, int index, string name);
	//GASP2pop getCluster(int index);
	//GASP2pop getIndv(UUID u);
	GASP2pop getGen(int best, int gen, string name);
	GASP2pop getBest(int best, string name);

	GASP2pop getIncomplete(string name);

	//input tables
	int addInput(string infile);

private:
	bool openState;
	string path;
	sqlite3 *dbconn;

	GASP2pop getxml(string sql);

	//GASP2struct getItem(UUID id);
	//void updateItem(UUID id, GASP2struct item); //can create or update

};

