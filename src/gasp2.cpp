#include "gasp2.hpp"

using namespace std;

//this constructor sets up the client
//side control scheme, and the placeholders
//for threads, hostfiles, and copy structures.

std::mutex eval_mut;

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
//	checkSeed(params.seed);
	MPI_Status s;
	MPI_Recv(&params.seed, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &s);
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
	if(restart.length() > 0) {
		doc.LoadFile(restart.c_str());
		if(doc.ErrorID() == 0) {
			tinyxml2::XMLElement * pop = doc.FirstChildElement("mgac")->FirstChildElement("pop");
			if(!rootpop.loadXMLrestart(pop, errorstring)) {
				cout << "There was an error in the restart file: " << errorstring << endl;
				MPI_Abort(1,MPI_COMM_WORLD);
			}

		}
		else {
			cout << "!!! There was a problem with opening the RESTART file!" << endl;
			cout << "Check to see if the file exists or if the XML file" << endl;
			cout << "is properly formed, with tags formatted correctly." << endl;
			cout << "Aborting... " << endl;
			MPI_Abort(1,MPI_COMM_WORLD);
		}
	}



	//random gen stuff
//	checkSeed(params.seed);

	for(int i = 1; i < worldSize; i++)
		MPI_Send(&params.seed, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

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


	Instruction recv, send;
	//this vector indicates the owner of all nodes
	//so, if ownerlist[i] = 0, node i is owned by
	//the server thread
	//0-worldsize means the node is running
	//IDLE (-1) means the node is idle
	//DOWN (-2) means the node is in a bad state
	//and is not responding to pings
	vector<int> ownerlist(worldSize);
	//this vector contains a ping test for each node

	vector<bool> pinged(worldSize);
	//the polling time is 200 ms
	chrono::milliseconds t(200);
	chrono::milliseconds timeout(0);

	//initialize ownerlist
	ownerlist[0] = IDLE;

	//establish who is paying attention initially
	for(int i = 1; i < worldSize; i++) {
		//try to get a signal for ~5 seconds
		sendIns(Ping, i);
		cout << mark() << "Sending initial ping to client " << i << endl;
	}
	this_thread::sleep_for(chrono::seconds(1));
	for(int i = 1; i < worldSize; i++) {
		recvIns(recv, i);
		if(recv & Ackn) {
			cout << mark() << "Client " << i << " ack'd" << endl;
			ownerlist[i] = IDLE;
		}
		else
			ownerlist[i] = DOWN;

	}


	//do the initialization
	if(restart.length() == 0) {
		rootpop.init(root, params.popsize, params.spacemode);
	}
	else
		cout << mark() << "Starting from restart population \"" << restart << "\"\n";
	writePop(rootpop, "root", 0);

	int replace = static_cast<int>(static_cast<double>(params.popsize)*params.replacement);

	//calc the total number of processing elements
	int totalPE = 0;
	for(int i = 0; i < worldSize; i++) {
		totalPE += hostlist[i].threads;

	}


	//pick the calculation mode
	Instruction evalmode = None;
	if(params.calcmethod == "qe")
		evalmode = (Instruction) (evalmode | DoQE);
//	else if(params.calcmethod == "fitcell")
//		evalmode = (Instruction) (evalmode | DoFitcell);

	GASP2pop evalpop, lastpop;
	lastpop = rootpop;
	vector<GASP2pop> split;

	vector<future<bool>> popsend(worldSize);
	future<bool> eval;
	future<bool> writer;

	MPI_Status m;

	if(params.mode == "stepwise") {
		for(int step = 0; step < params.generations; step++) {
			//scale, cross and mutate
			if(step > 0)
				lastpop.scale(params.const_scale, params.lin_scale, params.exp_scale);
			if(params.type == "elitism") {
				evalpop = lastpop.newPop(replace, params.spacemode);
				evalpop.addIndv(lastpop);
			}
			else if (params.type == "classic") {
				evalpop = lastpop.fullCross(params.spacemode);
			}
			evalpop.mutate(params.mutation_prob, params.spacemode);

			cout << mark() << "Crossover size: " << evalpop.size() << endl;


			//f
//			int count = evalpop.size() / totalPE;
//			int remainder = evalpop.size() % totalPE;
//			int remleft = remainder;
//			int index = 0;
			split.clear();
			split.resize(worldSize);
			evalpop.symmsort();

			//distribute the structures evenly so
			//that one host doesn't get too overloaded.
			for(int i = 0; i < evalpop.size(); i++) {
				split[i % worldSize].addIndv(*evalpop.indv(i));
			}
			cout << mark() << "Pop sizes: ";
			for(int i = 0; i < worldSize; i++)
				cout << split[i].size() << " ";
			cout << endl;


			cout << mark() << "Sending populations" << endl;
			//--addendum-- this may not be true because of MPI_thread_init
			//so it's still a safe assumption, probably
			//because of how MPI handles order of messages
			//blocking sends from the master thread must be issued
			//in the thread and not-asynchronously. multiple sends
			//issued in multiple threads creates potential for a
			//race condition, because the messages dispatched from
			//the master thread are out of order
			for(int i = 1; i < worldSize; i++) {
				sendIns(PopAvail, i);
				sendPop(split[i], i);
			}

			cout << mark() << "Performing fitcell" << endl;

			//send order to perform fitcell
			for(int i = 1; i < worldSize; i++) {
				sendIns(DoFitcell, i);
				ownerlist[i] = i;
			}
			//do the local fitcell
			split[0].runFitcell(hostlist[0].threads);


			//wait for completion of threads
			bool complete = false;
			while(!complete) {
				int recvd;
				//collect messages
				complete = true;
				for(int i = 1; i < worldSize; i++) {
					MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					while(recvd){
						recvIns(recv, i);
						if(recv & Busy)
							complete = false;
						//handle each message appropriately
						MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					}
				}
				for(int i = 1; i < worldSize; i++)
					sendIns(Ping, i);

				this_thread::sleep_for(t);
			}

			//receive final populations
			cout << mark() << "Retrieving populations" << endl;
			for(int i = 1; i < worldSize; i++) {
				sendIns(SendPop, i);
				recvPop(&split[i], i);
				cout << "i:" << i << endl;
			}



			cout << mark() << "Fitcell complete" << endl;


			//combine, then save fitcell
			lastpop = evalpop;
			evalpop.clear();
			for(int i = 0; i < worldSize; i++)
				evalpop.addIndv(split[i]);
			evalpop.volumesort();
			writePop(evalpop, "fitcell", step);




			//evaluate




			//sort and reduce
			evalpop.energysort();
			if(params.type == "elitism") {
				evalpop.remIndv(replace);
			}
			else if (params.type == "classic") {
				; //?
			}
			lastpop = evalpop;
			//save pops
			writePop(evalpop, "final", step);

		}

	}
	else if(params.mode == "steadystate") {
		cout << "Steadystate not yet implemented" << endl;
	}


//	ofstream outf;
//	outf.open(localmachinefile.c_str(), ofstream::out);
//	if(outf.fail()) {
//		cout << "Could not write machinefile!\n";
//		exit(1);
//	}
//	string name = makeMachinefile({0,1});
//	outf << name <<endl;
//	outf.close();
//	rootpop.init(root, params.popsize);
//	GASP2pop temp = rootpop.newPop(1);

	//clean up files
	remove(localmachinefile.c_str());
	for(int i = 1; i < worldSize; i++)
		sendIns(Shutdown, i);

	cout << mark() << "Shutting down" << endl;
	return;

}


bool GASP2control::runEvals(Instruction i, GASP2pop p, string machinefilename) {

	//get the hostfile stuff
//	string machinefile; int junk;
//	recvHost(machinefile, junk, 0);
//	ofstream outf;
//	outf.open(machinefilename.c_str(), ofstream::out);
//	if(outf.fail()) {
//		cout << "Could not write machinefile!\n";
//		exit(1);
//	}
//	outf << machinefile << endl;
//	outf.close();

	//only one of these will execute
	if(i & DoFitcell) {
		eval_mut.lock();
		p.runFitcell(nodethreads);
		eval_mut.unlock();
	}
	if(i & DoCharmm)
		cout << "Charmm not implemented!" << endl;
	if(i & DoQE)
		p.runEval(machinefilename, params, QE::runQE);
	if(i & DoCustom)
		cout << "Custom eval not implemented!" << endl;

	//cout << "eval finished on client " << ID << endl;
	return true;

}

void GASP2control::client_prog() {

	getHostInfo();
	//cout << "client prog: " << ID << endl;
	sendHost(hostname, nodethreads, 0);

	UUID u; u.generate();
	string localmachinefile = "/tmp/machinefile-";
	localmachinefile.append(u.toStr());


	MPI_Status m;
	bool busy = false;
	bool popsent = true;
	int recvd;
	GASP2pop evalpop;
	Instruction i, ack;
	//future<bool> popsend, eval;
	chrono::milliseconds timeout(0);
	chrono::milliseconds t(200);

	bool sendQueued = false;
	bool evalQueued = false;
	bool recvQueued = false;
	int target = 0;

	while(true) {
		i = None; ack = None;

		MPI_Iprobe(0, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		while(recvd) {
			recvIns(i, 0);

			//cout << "client received: " << i << endl;
			if(i & Shutdown) {
				remove(localmachinefile.c_str());
				return;
			}

			if(i & Ping){
				ack = (Instruction) (ack | Ackn);
			}

			if(i & SendPop || sendQueued) {
				if(eval_mut.try_lock()) {
					eval_mut.unlock();
					sendPop(evalpop, 0);
					sendQueued = false;
				}
				else {
					sendQueued = true;
					ack = (Instruction) (ack | Busy);
				}
			}

			if( ((i & PopAvail) || recvQueued) && !sendQueued ) {
				if(eval_mut.try_lock()) {
					eval_mut.unlock();
					//sendPop(evalpop, 0);
					recvPop(&evalpop, 0);
					recvQueued = false;
				}
				else {
					recvQueued = true;
					ack = (Instruction) (ack | Busy);
				}
			}

			if(i & (DoFitcell | DoCharmm | DoQE | DoCustom) || evalQueued){
				//if(!eval.valid() && eval_mut.try_lock()) {
				if(eval_mut.try_lock()) {
					eval_mut.unlock();
					//cout << "Evalpop size: " << evalpop.size() << endl;
					//eval = async(launch::async, &GASP2control::runEvals, this, i, evalpop, localmachinefile);
					async(launch::async, &GASP2control::runEvals, this, i, evalpop, localmachinefile);
					evalQueued = false;
				}
				else {
					//cout << mark() << "ERROR ON CLIENT " << ID << ": COULD NOT LAUNCH EVAL THREAD!" << endl;
					//MPI_Abort(MPI_COMM_WORLD,1);
					cout << mark() << "Eval thread is locked, retrying in ~200 ms..." << endl;
					evalQueued = true;
				}

			}

			if(eval_mut.try_lock()) {
				eval_mut.unlock();
			}
			else {
				ack = (Instruction) (ack | Busy);
			}
//
//					if(popsend.valid() && popsend.wait_for(timeout)==future_status::ready) {
//						popsend.get();
//					}
//					if(eval_mut.try_lock() && eval.valid()) {// && eval.valid() && eval.wait_for(timeout)==future_status::ready) {
//						eval.get();
//						eval_mut.unlock();
//						break;
//					}




			if(ack != None)
				sendIns(ack,0);

			MPI_Iprobe(0, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		}

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
	MPI_Send(&t, 1, MPI_INT, target,HOSTS,MPI_COMM_WORLD);

	//send info
	MPI_Send(host.c_str(),t,MPI_CHAR,target,HOSTS,MPI_COMM_WORLD);
	//send proc info
	MPI_Send(&procs, 1, MPI_INT, target, HOSTS, MPI_COMM_WORLD);
	return true;
}

bool GASP2control::recvHost(string &host, int &procs, int target) {
	int ierr;
	int v = 0;
	//recv size info
	MPI_Recv(&v, 1, MPI_INT, target,HOSTS,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	char * buff = new char[v];
	//recv hostname
	MPI_Recv((void *) buff, v, MPI_CHAR, target,HOSTS,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	host = buff;
	//recv proc count

	MPI_Recv(&procs, 1, MPI_INT, target, HOSTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//cout << "procs:" << procs << endl;
	delete [] buff;

	return true;
}


bool GASP2control::sendPop(GASP2pop p, int target) {
	string pop = p.saveXML();
	//cout << "Pop" << endl << pop << endl << endl;

	pop += "\0";
	int t = pop.length();

	cout << mark() << " Sendpop to " << target << " launched, ID " << ID << ", size: " << t <<  endl;
	//send size info
	MPI_Send(&t, 1, MPI_INT, target,POP,MPI_COMM_WORLD);

	//send info
	MPI_Send(pop.c_str(),t,MPI_CHAR,target,POP,MPI_COMM_WORLD);

	cout << mark() << " Sendpop " << target << " finished, ID " << ID << endl;

	return true;
}

bool GASP2control::recvPop(GASP2pop *p, int target) {
	string pop;
	int v = 0;

	//recv size info
	MPI_Recv(&v, 1, MPI_INT, target,POP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	cout << mark() << " Recvpop from " << target << " launched, ID " << ID << ", size: " << v << endl;
	//char * buff = new char[v+1];
	//recv hostname
	pop.resize(v+1);

	MPI_Recv((void *) pop.data(), v, MPI_CHAR, target,POP,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	//buff[v] = '\0';
	//pop = buff;

	//cout << "Pop" << endl << pop << endl << endl;

	//delete [] buff;

	string error;
	p->parseXML(pop,error);

	//cout << "Error" <<  error << endl;

	cout << mark() << " Recvpop " << target << " finished, ID " << ID << ", structsize: " << p->size() << endl;

	return true;
}


bool GASP2control::sendIns(Instruction i, int target) {
	//cout << "sending ins to target " << target << endl;

	MPI_Request r;
	MPI_Isend(&i, 1, MPI_INT, target, CONTROL, MPI_COMM_WORLD, &r);
	return true;
}

bool GASP2control::recvIns(Instruction &i, int target) {
	//cout << "receiving ins from target " << target << endl;
	MPI_Request r;
	MPI_Irecv(&i, 1, MPI_INT, target, CONTROL, MPI_COMM_WORLD, &r);
	//cout << "procs:" << procs << endl;

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


bool GASP2control::writePop(GASP2pop pop, string tag, int step) {

	stringstream outname;
	ofstream outf;
	outname.str("");
	outname << params.outputfile << "_" << tag << "_" << setw(3) << setfill('0') << step << ".pop";
	cout << mark() << "Saving " << tag << " \"" << outname.str() << "\"..." << endl;

	outf.open(outname.str().c_str(), ofstream::out);
	if(outf.fail()) {
		cout << mark() << "ERROR: COULD NOT OPEN FILE FOR SAVING! continuing sadly..." << endl;
		return false;
	}
	else {
		outf << "<mgac>\n" << pop.saveXML() << endl << "</mgac>\n";
		outf.close();
		cout << mark() << "Save complete" << endl;
	}

	return true;
}
