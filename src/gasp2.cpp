#include "gasp2.hpp"

using namespace std;

//this constructor sets up the client
//side control scheme, and the placeholders
//for threads, hostfiles, and copy structures.

GASP2control::GASP2control(int ID, string infile) {

	//parse the infile
	tinyxml2::XMLDocument doc;
	if(!doc.LoadFile(infile.c_str()))
		parseInput(&doc);
	else
		exit(1);


}

//this constructor is invoked for a server
//controller. the server acts as both client
//and server
GASP2control::GASP2control(time_t start, int size, string input, string restart) {
	starttime = start;
	worldSize = size;
	this->restart = restart;
	infile = input;
	int ID = 0;

	//parse the infile with output handling!
	tinyxml2::XMLDocument doc;
	doc.LoadFile(infile.c_str());
	if(doc.ErrorID() == 0)
		parseInput(&doc);
	else {
		cout << "!!! There was a problem with opening the input file! Aborting... " << endl;
		exit(1);
	}
	//parse the restart if it exists



}

void GASP2control::server_prog() {

}


void GASP2control::client_prog() {

}

void GASP2control::parseInput(tinyxml2::XMLDocument *doc) {
	params.parseXML(doc);


}
