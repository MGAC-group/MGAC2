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

public:

	int size() { return structures.size(); }

	//sorts the population on energy or volume
	//for volume, structure with volume closest
	//to the expected volume is best structure
	void energysort();
	void volumesort();

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
	//can remove
	GASP2pop remIndv(int sub);

	//removes structures from the pop that exceed the
	//volume limits
	GASP2pop volLimit(double min, double max);

	//produces a new structure list which contains only
	//the unique best structures of the list.
	//structures where the energy is known (completed)
	//take precedence over fitcelled structures
	//the lowest energy structure of a set takes precedence
	GASP2pop unique();

	//mutates the population
	void mutate(float rate);

	//fitcell
	GASP2pop runFitcell();
	GASP2pop runEval();

	//parsing handlers
	bool saveXML(string &name);
	bool parseXML(string &name);
	bool loadXMLrestart(tinyxml2::XMLElement *elem, string& errorstring);


};
