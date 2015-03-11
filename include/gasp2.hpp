#pragma once
#include "gasp2common.hpp"
#include "gasp2param.hpp"
#include "gasp2struct.hpp"
#include "gasp2pop.hpp"
//#include "gasp2qe.hpp"


using namespace std;

class GASP2control {
public:
	GASP2control(int ID, string infile);
	GASP2control(time_t start, int size, string input, string restart="");
	void server_prog();
	void client_prog();

	void getHostInfo();

private:
	//procedural variables
	time_t starttime;
	string infile;
	string restart;
	int worldSize;
	int ID;

	string hostname;
	int nodethreads;

	GASP2param params;
	GASP2struct root; //base structure which all other structures are derived from
	vector<GASP2pop> populations;


	bool parseInput(tinyxml2::XMLDocument *doc, string & errors);

};
