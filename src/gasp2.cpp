#include "gasp2.hpp"

using namespace std;

//this constructor sets up the client
//side control scheme, and the placeholders
//for threads, hostfiles, and copy structures.

std::mutex eval_mut;
std::mutex longeval_mut;
std::atomic<bool> save_state;
std::atomic<bool> completed;

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


//this function is awful and really needs to be cleaned up. badly
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
		else {
			cout << mark() << "Client " << i << " is down!" << endl;
			ownerlist[i] = DOWN;
		}

	}


	//do the initialization
	if(restart.length() == 0) {
		rootpop.init(root, params.popsize, params.spacemode, params.group);
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
			cout << mark() << "Starting new generation " << step << endl;
			//scale, cross and mutate
			cout << mark() << "Scaling..." << endl;
			if(step > 0)
				lastpop.scale(params.const_scale, params.lin_scale, params.exp_scale);
			cout << mark() << "Crossing..." << endl;

			if(params.type == "elitism") {
				evalpop = lastpop.newPop(replace, params.spacemode);
				evalpop.addIndv(lastpop);
			}

			else if (params.type == "classic") {
				lastpop.addIndv(rootpop);
				evalpop = lastpop.fullCross(params.spacemode);
			}
			evalpop.mutate(params.mutation_prob, params.spacemode);

			cout << mark() << "Crossover size: " << evalpop.size() << endl;


			evalpop = evalpop.symmLimit(params.symmlimit);


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
			cout << mark() << "Server fitcell finished" << endl;


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
				//cout << "i:" << i << endl;
			}

			cout << mark() << "Fitcell complete" << endl;


			//combine, then save fitcell
			lastpop = evalpop;
			evalpop.clear();
			for(int i = 0; i < worldSize; i++)
				evalpop.addIndv(split[i]);
			evalpop.volumesort();
			//writePop(evalpop, "fitcell", step);

			//evaluate
			evalpop.symmsort();
			GASP2pop bad, good, restart;
			good = evalpop.volLimit(bad);
			cout << mark() << "Volume limited population size:" << good.size() << endl;

			good.volumesort();
			writePop(good, "vollimit", step);
			future<bool> serverthread;

			if(good.size() > 0) {

				//establish who is paying attention initially
				for(int i = 1; i < worldSize; i++) {

					//clear the message queue for safety
					int recvd;
					recv = None;
					MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					while(recvd){
						recvIns(recv, i);
						MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					}

					//try to get a signal for ~5 seconds
					sendIns(Ping, i);
					cout << mark() << "Sending initial ping to client " << i << endl;
				}
				this_thread::sleep_for(chrono::seconds(2));
				for(int i = 1; i < worldSize; i++) {
					recvIns(recv, i);
					if(recv & Ackn) {
						cout << mark() << "Client " << i << " ack'd" << endl;
						ownerlist[i] = IDLE;
					}
					else {
						cout << mark() << "Client " << i << " is down!" << endl;
						ownerlist[i] = DOWN;
					}

				}

				if(params.calcmethod == "qe") {

					cout << mark() << "Beginning QE evaluation" << endl;
					good.symmsort();
					int avail = 0;
					//var to store reference index of sent structure
					vector<GASP2pop> localpops(worldSize);
					vector<int> evaldind(worldSize);
					GASP2pop evald = good;

					for(int launched=0,completesteps=0; ; ) {
						cout << mark() << "start of launch block" << endl;

						//test node 0
						if(ownerlist[0] == 0) {
							if(eval_mut.try_lock()) {
								eval_mut.unlock();
								serverthread.get();
								ownerlist[0]=IDLE;
								evald.mergeIndv(localpops[0],evaldind[0]);
								completesteps++;
							}
						}

						//establish active/idle nodes, receive pops
						cout << mark() << "pinging" << endl;
						for(int i = 1; i < worldSize; i++) {
							if(ownerlist[i] == 1)
								sendIns(Ping, i);
						}
						this_thread::sleep_for(chrono::seconds(1));
						for(int i = 1; i < worldSize; i++) {
							if(ownerlist[i] == i) {
								int recvd;
								recv = None;
								MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
								while(recvd){
									recvIns(recv, i);

									cout << mark() << "server QE received from " << i << ": " << recv << endl;

										if(recv & PopAvail) {
											sendIns(SendPop, i);
											recvPop(&localpops[i], i);
		//										cout << mark() << "gasp server pop energies" << endl;
		//										for(int j = 0; j < evalpop.size(); j++)
		//											cout << "energy: " << localpops[i].indv(j)->getEnergy() << endl;
											evald.mergeIndv(localpops[i],evaldind[i]);
											restart = evald;
											restart.addIndv(bad);
											writePop(restart, "restart", 0);
										}
										if(recv & Busy) {
											cout << "marking busy " << i << endl;
											ownerlist[i] = i;
										}
										else if(recv & Complete) {
											sendIns(SendPop,i);
											recvPop(&localpops[i],i);
											evald.mergeIndv(localpops[i],evaldind[i]);
											restart = evald;
											restart.addIndv(bad);
											writePop(restart, "restart", 0);

											cout << mark() << "Client " << i << " finished eval" << endl;
											ownerlist[i] = IDLE;
											completesteps++;
										}

									MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
								}
							}
						}


						//cleanup node control lists
						for(int i = 0; i < worldSize; i++) {
							if(ownerlist[ownerlist[i]] == IDLE)
								ownerlist[i] = IDLE;
						}

						avail = 0;
						for(int i = 0; i < worldSize; i++)
							if(ownerlist[i] == IDLE)
								avail++;

						cout << mark() << "structures avail block: " << avail << endl;
						cout << "Ownerlist: ";
						for(int i = 0; i < worldSize; i++)
							cout << ownerlist[i] << " ";
						cout << endl;


						//structures available
						while(avail > 0 && launched < good.size()) {
							//see how many nodes are needed
							int given;
							int needed=1*good.indv(launched)->getSymmcount();
							if(needed > avail)
								given = avail;
							//see if there will be leftovers and how much is next
							//if there are leftovers just give the extras to this
							//evaluation (lazy scheduling, but it works)
							else {
								int next = 0;
								if(launched + 1 < good.size())
									next = 1*good.indv(launched+1)->getSymmcount();
								if(next > (avail-needed))
									given = avail;
								else if (next == 0)
									given = avail;
								else
									given = needed;

							}

							cout << mark() << "given is " << given << endl;

							vector<int> slots;
							slots.clear();
							int first = IDLE;
							for(int i = 0; i < worldSize; i++) {
								if(ownerlist[i] == IDLE) {
									if(first == IDLE)
										first = i;
									slots.push_back(i);
									ownerlist[i] = first;
								}
								if(slots.size() >= given)
									break;
							}


							cout << mark() << "first is " << first << endl;
							if(first == IDLE)
								break;

							avail -= given;
							cout << mark() << "avail is " << avail << endl;

							string hostfile = makeMachinefile(slots);

							if(first == 0) {
								writeHost(localmachinefile, hostfile);
							}
							else {
								sendIns(GetHost, first);
								sendHost(hostfile, 0, first);
							}


							evaldind[first] = launched;
							localpops[first] = good.subpop(launched,1);

							if(first == 0) {
								if(eval_mut.try_lock()) {
									cout << mark() << "Launching QE on server with "<< given << "nodes..."<< endl;
									eval_mut.unlock();
									//cout << "Evalpop size: " << evalpop.size() << endl;
									//eval = async(launch::async, &GASP2control::runEvals, this, i, evalpop, localmachinefile);
									serverthread = std::async(launch::async, &GASP2control::runEvals, this, DoQE, &localpops[first], localmachinefile);
								}
								launched++;
								break;
							}
							else {
								cout << mark() << "Launching QE on client " << first << "with "<< given << "nodes..." << endl;
								sendIns(PopAvail, first);
								sendPop(localpops[first], first);
								sendIns(DoQE, first);
								launched++;
								break;
							}


						}

						cout << mark() << "launched is " << launched << endl;
						//test for completion (all nodes idle)
						bool flag = true;
						for(int i = 0; i < worldSize; i++) {
							if(ownerlist[i] != IDLE && ownerlist[i] != DOWN )
								flag = false;
						}
						if(flag && (launched >= good.size()) && (completesteps >= good.size()) )
							break;

						this_thread::sleep_for(chrono::seconds(2));

					}

					evald.energysort();
					writePop(evald, "evald", step);
					good = evald;
					good.addIndv(bad);
					evalpop = good;

					cout << mark() << "End of QE evalution" << endl;

				} //end QE eval block

			} //if good size > 0
			else {
				cout << mark() << "No structures to be energy evaluated, continuing..." << endl;

			}



			//sort and reduce
			evalpop.energysort();
			if(params.type == "elitism") {
				evalpop.remIndv(replace);
			}
			else if (params.type == "classic") {
				evalpop.remIndv(evalpop.size()-params.popsize);
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


bool GASP2control::runEvals(Instruction i, GASP2pop* p, string machinefilename) {

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
	eval_mut.lock();
	//only one of these will execute
	if(i & DoFitcell) {
		p->runFitcell(nodethreads);
	}
	if(i & DoCharmm)
		cout << "Charmm not implemented!" << endl;
	if(i & DoQE) {
		p->runEval(machinefilename, params, QE::runQE);
	}
	if(i & DoCustom)
		cout << "Custom eval not implemented!" << endl;

	eval_mut.unlock();

	cout << "Eval finished on client " << ID << endl;
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
	GASP2pop evalpop, temppop;
	Instruction i, ack, evalmode = (Instruction)(0u);
	future<bool> popsend, eval;
	chrono::milliseconds timeout(0);
	chrono::milliseconds t(200);

	bool queueSend = false;
	bool queueEval = false;
	bool queueRecv = false;
	int target = 0;

	auto start = chrono::steady_clock::now();

	string hostfile;
	int itemp;

	save_state=false;
	completed = false;

	while(true) {
		i = None; ack = None;

		//cout << "client cycle " << this->ID << endl;

		MPI_Iprobe(0, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		while(recvd) {
			recvIns(i, 0);

			cout << "client " << this->ID << "received: " << i << endl;
			if(i & Shutdown) {
				remove(localmachinefile.c_str());
				return;
			}

			if(i & Ping){
				ack = (Instruction) (ack | Ackn);
			}

			if(i & GetHost) {
				recvHost(hostfile, itemp, 0);
				writeHost(localmachinefile, hostfile);
			}

			if(i & PopAvail)
				queueRecv = true;

			if(i & SendPop)
				queueSend = true;

			//do not allow the evalmode to change
			if(i & (DoFitcell | DoCharmm | DoQE | DoCustom) && !queueEval) {
				evalmode = i;
				queueEval = true;
			}

			if(eval_mut.try_lock()) {
				eval_mut.unlock();
				//cout << mark() << "sending Complete" << endl;
				//ack = (Instruction) (ack | Complete);
			}
			else {
				ack = (Instruction) (ack | Busy);
			}

			if(ack != None) {
				cout << "client " << ID << " sent ack " << ack << endl;
				sendIns(ack,0);
			}

			MPI_Iprobe(0, CONTROL, MPI_COMM_WORLD, &recvd, &m);

		}

		//execute work
		//A SEND ALWAYS TAKES PRECEDENCE OVER A RECEIVE

		if(queueSend) {
			//see what the current evalmode
			cout << "client " << ID << ": send queued" << endl;
			if(evalmode & DoFitcell) {
				if(eval_mut.try_lock()) {
					sendPop(evalpop, 0);
					queueSend = false;
					eval_mut.unlock();
				}
				else {
					cout << "client " << ID << ": cannot send, thread still working (retrying in 200 ms)" << endl;
				}
			}
			else if(evalmode & DoQE) {
				//evaluation not happening
				if(eval_mut.try_lock()) {
					sendPop(evalpop, 0);
					queueSend = false;
					eval_mut.unlock();
				}
				else {
					//in this instance the client sends info
					sendPop(temppop, 0);
					queueSend = false;
				}

			}
		}

		if(queueRecv) {
			cout << "client " << ID << ": recv queued" << endl;
			if(!queueSend) {
				if(eval_mut.try_lock()) {
					eval_mut.unlock();
					recvPop(&evalpop, 0);
					queueRecv = false;
				}
				else {
					cout << "client " << ID << ":Eval thread is locked, cannot send yet..." << endl;
				}

			}
			else {
				cout << "client " << ID << ": send queued, cannot receive yet" << endl;
			}
		}

		if(queueEval) {
			cout << "client " << ID << ": eval queued" << endl;
			if(eval_mut.try_lock()) {
				eval_mut.unlock();
				temppop = evalpop;
				async(launch::async, &GASP2control::runEvals, this, evalmode, &evalpop, localmachinefile);
				queueEval = false;
			}
			else {
				cout << mark() << "Eval thread is locked, retrying in ~200 ms..." << ID << endl;
			}
		}

		//store the temppop
		if( save_state && longeval_mut.try_lock()) {
			cout << "queueing intermediate pop info on client " << ID << endl;
//					cout << mark() << "gasp client pop energies" << endl;
//					for(int i = 0; i < evalpop.size(); i++)
//						cout << "energy: " << evalpop.indv(i)->getEnergy() << endl;
			temppop = evalpop;
			sendIns(PopAvail, 0);
			save_state=false;
			longeval_mut.unlock();
		}

		if(completed) {
			sendIns(Complete, 0);
			completed = false;
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

void GASP2control::writeHost(string name, string data) {
		ofstream outf;
		outf.open(name.c_str(), ofstream::out);
		if(outf.fail()) {
			cout << "Could not write machinefile!\n";
			exit(1);
		}
		outf << data <<endl;
		outf.close();
}
