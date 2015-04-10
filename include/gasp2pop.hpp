#pragma once
#include "gasp2common.hpp"
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
	void clear() { structures.clear(); scaling.clear(); }

	//sorts the population on energy or volume
	//for volume, structure with volume closest
	//to the expected volume is best structure
	void energysort();
	void volumesort();


	void init(GASP2struct s, int size);

	//generate a new population by crossing the
	//existing population (for classic)
	//using whatever selection scheme is available.
	GASP2pop newPop(int size, GAselection mode=Roulette);

	//performs a full crossing of all members of
	//the population and generates a new one
	GASP2pop fullCross();

	//add N members to the population (for elitism)
	void addIndv(int add);
	//adds individuals from another population
	//to the existing population
	void addIndv(GASP2pop add);

	//remove N worst members from the population (for elitism)
	//return the removed members
	//removes based on sort order,
	//so, last n structures get removed regardless of whether
	//the pop is sorted or not
	GASP2pop remIndv(int n);

	//removes structures from the pop that exceed the
	//volume limits (based on scaled scores);
	GASP2pop volLimit(GASP2pop &bad=nullptr);

	//produces a new structure list which contains only
	//the unique best structures of the list.
	//structures where the energy is known (completed)
	//take precedence over fitcelled structures
	//the lowest energy structure of a set takes precedence
//	GASP2pop unique();

	//mutates the population
	void mutate(double rate);

	//fitcell/eval functions
	GASP2pop runFitcell(int threads); //threaded
	GASP2pop runEval(string hosts); //serial parallel
	//hosts is a string in the form of a machinefile in the /tmp dir
	//using a UUID; so, /tmp/UUID.hosts
	//should prevent hiccups with existing file handles


	//parsing handlers
	string saveXML();
	bool parseXML(string name, string &errorstring);
	bool loadXMLrestart(tinyxml2::XMLElement *elem, string& errorstring);

	bool writeCIF(string name);

private:
	vector<double> scale(double con, double lin, double exp);


};
