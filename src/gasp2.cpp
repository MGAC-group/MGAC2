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

//	for(int i = 0; i < 2; i++) {
//		root.unfitcell();
//		root.init();
//		root.fitcell();
//		if(!root.cifOut("fitcelldebug.cif"))
//			cout << "Bad file for fitcelldebug!\n";
//		cout << "i" << i << endl;
//	}


//	root.init();
//	root.fitcell();
//	if(!root.cifOut("fitcelldebug.cif"))
//		cout << "Bad file for fitcelldebug!\n";
//
//	for(int i = 0; i < 10; i++) {
//		root.unfitcell();
//		root.mutateStruct(0.10);
//		root.fitcell();
//		if(!root.cifOut("fitcelldebug.cif"))
//			cout << "Bad file for fitcelldebug!\n";
//	}


//	root.init();
//	for(int i = 2; i < 3; i++) {
//		root.unfitcell();
//		root.overrideSpacegroup(i);
//		root.fitcell();
//		if(!root.cifOut("fitcelldebug.cif"))
//			cout << "Bad file for fitcelldebug!\n";
//	}
//
//	root.unfitcell();

	//cout << root.serializeXML() << endl;

	getHostInfo();

	rootpop.init(root, params.popsize);




	GASP2pop temp = rootpop.newPop(1);

	temp.runFitcell(8);
	temp.volumesort();

	string hosts("blahblahblahblah");
	//temp.remIndv(temp.size() - 1);
	temp.runEval(hosts, params, QE::empty);
	temp.scale(0.1,1.0,5.0);

	cout << "energysort\n";
	temp.energysort();
	cout << "volumesort\n";
	temp.volumesort();




	//
//	cout << "out size: " << out.size() << endl;
//	out.addIndv(10);
//	cout << "out size: " << out.size() << endl;
//	out.remIndv(5);
//	cout << "out size: " << out.size() << endl;


//	//out.runFitcell(8);
//	//out.volumesort();
//
//	//GASP2pop out2 = out.volLimit();
//
//	//cout << "good" << endl;
//	//out2.volumesort();
//
//
//	//out2.remIndv(out2.size()-2);
//	//cout << "out2 size: " << out2.size() << endl;
//	//out2.volumesort();
//
//	//GASP2pop out3 = out2.fullCross();
//
//	out3.runFitcell(8);
//	out3.writeCIF("fitcelldebug.cif");
//	out3.volumesort();
//	cout << "out2 size: " << out2.size() << endl;
//	cout << "out3 size: " << out3.size() << endl;
//	out3.addIndv(out2);
//	cout << "out3 size: " << out3.size() << endl;

	//string s = out4.saveXML();
	//cout << s << endl;

	//out.writeCIF("fitcelldebug.cif");


	//string out = rootpop.saveXML();


	//cout << out << endl;



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

	if(ID == 0) {
		params.logParams();
		root.logStruct();
	}
	return true;




}
