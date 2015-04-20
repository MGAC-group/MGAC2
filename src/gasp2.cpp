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

	//collect the info from all hosts
	Host h; h.hostname = hostname;
	h.threads = nodethreads;
	hostlist.push_back(h);
	for(int i = 1; i < worldSize; i++) {
		recvHost(h.hostname, h.threads, i);
		hostlist.push_back(h);
	}

	//generate the local machinefile info
	UUID u; u.generate();
	string localmachinefile = "/tmp/machinefile-";
	localmachinefile.append(u.toStr());

	ofstream outf;
	outf.open(localmachinefile.c_str(), ofstream::out);
	if(outf.fail()) {
		cout << "Could not write machinefile!\n";
		exit(1);
	}
	string name = makeMachinefile({0,1});
	outf << name <<endl;
	outf.close();



	rootpop.init(root, params.popsize);
	GASP2pop temp = rootpop.newPop(1);

	temp.runFitcell(8);
	temp.volumesort();

	string hosts("blahblahblahblah");
	temp.remIndv(temp.size() - 1);
	temp.writeCIF("pre.cif");
	temp.runEval(localmachinefile, params, QE::runQE);
	temp.indv(0)->forceOK();
	temp.writeCIF("post.cif");

	//clean up files
	remove(localmachinefile.c_str());
}


void GASP2control::client_prog() {
	getHostInfo();
	//cout << "client prog: " << ID << endl;
	sendHost(hostname, nodethreads, 0);

//	UUID u; u.generate();
//	string localmachinefile = "/tmp/machinefile-";
//	localmachinefile.append(u.toStr());
//
//	remove(localmachinefile.c_str());
	while(true) {
        chrono::milliseconds t(200);
        this_thread::sleep_for(t);
	};

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

bool GASP2control::sendHost(string host, int procs, int target) {
	int t = host.length();
	//send size info
	MPI_Send(&t, 1, MPI_INT, target,0,MPI_COMM_WORLD);

	//send info
	MPI_Send(host.c_str(),t,MPI_CHAR,target,0,MPI_COMM_WORLD);
	//send proc info
	MPI_Send(&procs, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
	return true;
}

bool GASP2control::recvHost(string &host, int &procs, int target) {
	int ierr;
	int v = 0;
	//recv size info
	MPI_Recv(&v, 1, MPI_INT, target,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	char * buff = new char[v];
	//recv hostname
	MPI_Recv((void *) buff, v, MPI_CHAR, target,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	host = buff;
	//recv proc count

	MPI_Recv(&procs, 1, MPI_INT, target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//cout << "procs:" << procs << endl;
	delete [] buff;

	return true;
}

string GASP2control::makeMachinefile(vector<int> slots) {
	stringstream out;
	out.str("");

	for(int i = 0; i < slots.size(); i++) {
		if(slots[i] < hostlist.size()) {
			for(int j = 0; j < hostlist[slots[i]].threads; j++)
				out << hostlist[slots[i]].hostname << endl;
		}
	}
	return out.str();
}
