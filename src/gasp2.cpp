#include "gasp2.hpp"

using namespace std;

//this constructor sets up the client
//side control scheme, and the placeholders
//for threads, hostfiles, and copy structures.

GASP2control::GASP2control(int ID, string infile) {

	this->ID = ID;

	//parse the infile
	tinyxml2::XMLDocument doc;
	string errorstring;
	if(!doc.LoadFile(infile.c_str()))
		parseInput(&doc, errorstring);
	else
		exit(1);

	//rgen is in gasp2common
	//kinda wonky, but is okay
	rgen.seed(params.seed+ID);

}

//this constructor is invoked for a server
//controller ONLY. the server acts as both client
//and server by dispatching pthreads in addition
//to managing the other clients
GASP2control::GASP2control(time_t start, int size, string input, string restart) {
	starttime = start;
	worldSize = size;
	this->restart = restart;
	infile = input;
	ID = 0;

	//parse the infile with output handling!
	tinyxml2::XMLDocument doc;
	string errorstring;
	doc.LoadFile(infile.c_str());
	if(doc.ErrorID() == 0) {
		if(parseInput(&doc, errorstring)==false) {
			cout << "There was an error in the input file: " << errorstring << endl;
			exit(1);
		}

	}
	else {
		cout << "!!! There was a problem with opening the input file!" << endl;
		cout << "Check to see if the file exists or if the XML file" << endl;
		cout << "is properly formed, with tags formatted correctly." << endl;
		cout << "Aborting... " << endl;
		exit(1);
	}


	//parse the restart if it exists

	rgen.seed(params.seed);

}

void GASP2control::getHostInfo() {

	FILE * p =  popen("cat /proc/cpuinfo | grep \"physical id\" | sort | uniq | wc -l","r");
	FILE * c =  popen("cat /proc/cpuinfo | grep \"core id\" | sort | uniq | wc -l","r");
	char p_str[10], c_str[10];

	fgets(p_str, 10, p);
	fgets(c_str, 10, c);

	stringstream ss;
	int phys, core;

	ss << p_str; ss >> phys; ss.clear();
	ss << c_str; ss >> core; ss.clear();

	nodethreads = phys*core;
	char _host[1024];
	gethostname(_host,1024);
	hostname = _host;

}

void GASP2control::server_prog() {
	getHostInfo();
}


void GASP2control::client_prog() {
	getHostInfo();
}

bool GASP2control::parseInput(tinyxml2::XMLDocument *doc, string& errors) {

	//pass the doc to gasp2param
	if(params.parseXML(doc, errors)==false)
		return false;
	if(root.parseXMLDoc(doc, errors)==false)
		return false;


	for(int i = 0; i < 3; i++) {
		root.init();

	}

	if(ID == 0) {
		params.logParams();
		root.logStruct();
	}
	return true;




}
