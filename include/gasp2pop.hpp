#pragma once
#include "gasp2common.hpp"
#include "gasp2param.hpp"
#include "gasp2struct.hpp"

using namespace std;

typedef enum GAselection {
	Roulette=0, Pattern,
}GAselection;

typedef enum GAscaling {
	Con, Lin, Exp,
}GAscaling;

typedef enum GAmode {
	Steady, Strict
}GAmode;

typedef enum GAtype {
	Elite, Classic,
}GAtype;


class GASP2pop {
private:
	vector<GASP2struct> structures;
	vector<double> scaling;

public:

	int size() { return structures.size(); }
	GASP2struct *indv(int n) { if(n < size()) return &structures[n]; else return &structures[0];};
	void clear() { structures.clear(); scaling.clear(); }

	//sorts the population on energy or volume
	//for volume, structure with volume closest
	//to the expected volume is best structure
	void energysort();
	void volumesort();
	//sort based on symmetry operations
	void symmsort();

	GASP2pop subpop(int start, int subsize);

	void init(GASP2struct s, int size, Spacemode mode, Index spcg = 1);

	//generate a new population by crossing the
	//existing population (for classic)
	//using whatever selection scheme is available.
	GASP2pop newPop(int size, Spacemode smode, GAselection mode=Roulette);

	//performs a full crossing of all members of
	//the population and generates a new one
	GASP2pop fullCross(Spacemode mode);
	GASP2pop fullCross(Spacemode mode, GASP2pop alt);

	//performs an inplace cross,
	//where individuals in the population are mixed
	//with each other
	//this is the actual crossing algorithm for a classic GA
	GASP2pop inplaceCross(Spacemode mode);


	//add N members to the population (for elitism)
	void addIndv(int add);
	//adds individuals from another population
	//to the existing population
	void addIndv(GASP2pop add);
	void addIndv(GASP2struct t){ structures.push_back(t); };
	void mergeIndv(GASP2pop add, int ind);

	//checks to see if a structure is compelte
	//if not, then it is added into the output pop
	//and opt is to False for that structure

	GASP2pop completeCheck();

	//remove N worst members from the population (for elitism)
	//returns the removed members
	//removes based on sort order,
	//so, last n structures get removed regardless of whether
	//the pop is sorted or not
	GASP2pop remIndv(int n);

	//removes structures from the pop that exceed the
	//volume limits (based on scaled scores);
	GASP2pop volLimit() {GASP2pop t; return volLimit(t);};
	GASP2pop symmLimit(int limit) {GASP2pop t; return symmLimit(t, limit);};
	GASP2pop volLimit(GASP2pop &bad);
	GASP2pop symmLimit(GASP2pop &bad, int limit);
	GASP2pop spacebin(int binsize, int binsave = 230);
	void spacebinV(vector<GASP2pop> &bins, int binsize, int binsave = 230);


	void scale(double con, double lin, double exp);
	//produces a new structure list which contains only
	//the unique best structures of the list.
	//structures where the energy is known (completed)
	//take precedence over fitcelled structures
	//the lowest energy structure of a set takes precedence
//	GASP2pop unique();

	//mutates the population
	void mutate(double rate, Spacemode mode);

	//fitcell/eval functions
	void runFitcell(int threads); //threaded
	void runSymmetrize(int threads);
	void runEval(string hosts, GASP2param params, bool (*eval)(vector<GASP2molecule>&, GASP2cell&, double&, double&, double&, time_t&, string, GASP2param)); //serial parallel
	//hosts is a string in the form of a machinefile in the /tmp dir
	//using a UUID; so, /tmp/UUID.hosts
	//should prevent hiccups with existing file handles


	//parsing handlers
	string saveXML();
	bool parseXML(string name, string &errorstring);
	bool loadXMLrestart(tinyxml2::XMLElement *elem, string& errorstring);

	bool writeCIF(string name);

	//remove duplicate UUIDs
	//since the chance of a collision is so low, is not problematic
	void dedup(int max);


};
