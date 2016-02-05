#include "gasp2.hpp"

using namespace std;

//this constructor sets up the client
//side control scheme, and the placeholders
//for threads, hostfiles, and copy structures.

std::mutex eval_mut;
std::mutex longeval_mut;
std::atomic<bool> save_state;
std::atomic<bool> completed;

//client version constructor
GASP2control::GASP2control(int ID) {

	this->ID = ID;


	string infile;
	int temp;

	//get the infile via text receive
	recvHost(infile, temp, 0);

	//parse the infile
	tinyxml2::XMLDocument doc;
	string errorstring;
	if(!doc.Parse(infile.c_str()))
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


	ID = 0;
	starttime = start;
	worldSize = size;

	startstep=0;

	this->restart = restart;
	infile = "";


	//parse the infile with output handling!
	if(input.length() < 1 && restart.length() < 1) {
		cout << "GASP2control was passed an empty restart AND an empty infile!" << endl;
		exit(1);
	}

	tinyxml2::XMLDocument doc;
	string errorstring;



	//parse the restart if it exists
	if(restart.length() > 0) {
		params.outputfile = this->restart;
		if(!db.load(restart)) {
			cout << "There was an error loading the restart file!" << endl;
			MPI_Abort(1,MPI_COMM_WORLD);
			exit(1);
		}
		//get the last infile
		infile = db.getLastInput(startstep);
		startstep += 1;

	}


	//read the input file if it exists
	//this version of the infile supercedes the restart input
	if(input.length() > 0) {
		infile = get_file_contents(input.c_str());
	}


	//process the input files
	if( doc.Parse(infile.c_str()) ) {
		cout << "!!! There was a problem with parsing the input file!" << endl;
		cout << "Check to see if the file exists or if the XML file" << endl;
		cout << "is properly formed, with tags formatted correctly." << endl;
		cout << "Aborting... " << endl;
		exit(1);
	}
	if(parseInput(&doc, errorstring)==false) {
		cout << "There was an error in the input file: " << errorstring << endl;
		exit(1);
	}

	//if there wasn't a restart, do the initial setups
	if(restart.length() < 1) {
		string suffix = ".sq3";
		db.load(params.outputfile + suffix);
		db.init();
	}

	//add a new input file to the startup list
	db.addInput(infile, start);
	db.updateTime(start, startstep);

	//random gen stuff
//	checkSeed(params.seed);

	for(int i = 1; i < worldSize; i++)
		sendHost(infile, 0, i);

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
	hostname[1023] = '\0';
	gethostname(_host,1023);
	hostname = _host;
	cout << "Hostname ID " << ID << ":" << hostname << endl;

}


//handles the restart logic

//FIXME FIXME FIXME
void GASP2control::setup_restart() {

	if(restart.length() > 0) {

		cout << mark() << "Starting from restart population \"" << restart << "\"\n";


		GASP2pop incomplete = db.getIncomplete("structs");
		if(incomplete.size() > 0) {
			server_qe(incomplete, startstep);
		}
		incomplete.clear();
		GASP2pop complete = db.getAll("structs");

		complete.spacebinV(bins,params.popsize);

	}
}

//		//the rootpop might have zero energy structures
//		//those need to be evaluated first
//
//		rootpop = rootpop.symmLimit(params.symmlimit);
//		rootpop.runSymmetrize(hostlist[0].threads);
//		GASP2pop zeros;
//		GASP2pop evald = rootpop.energysplit(zeros);
//		zeros = zeros.volLimit();
//
//		if(params.outputmode == "xml") {
//			writePopXML(zeros, "zeros", 0);
//			writePopXML(evald, "res-evald", 0);
//		}
//
//		if(zeros.size() > 0) {
//			server_qe(zeros,startstep);
//			evald.addIndv(zeros);
//		}
//		rootpop = evald;
//
//		//implicit: rootpop is the restart pop if a restart is specified
//
////		if (params.type == "clustered")
////			rootpop.spacebinCluster(hostlist[0].threads, bins, clusters, params);
////		else
//		rootpop.spacebinV(bins, params.popsize);
//		if(params.outputmode == "xml") {
//			writePopXML(rootpop, "root", 0);
//		}
//	}
//	else {
//		//rootpop is not initialized
//		server_randbuild();
//		rootpop = randpop;
//	}
//}


void GASP2control::server_fitcell(GASP2pop &pop) {
	//distribute the structures evenly so
	//that one host doesn't get too overloaded.
	//now takes into account IDLE/DOWN state
	chrono::milliseconds t(200);
	chrono::milliseconds timeout(0);
	MPI_Status m;
	Instruction recv, send;
	vector<GASP2pop> split(worldSize);
	cout << mark() << "Fitcell starting: " << pop.size() << endl;

	int j = 0;
	for(int i = 0; i < pop.size(); i++) {
		for( ; ; j=((j+1) % worldSize) ) {
			if(ownerlist[j] == IDLE) {
				split[j].addIndv(*pop.indv(i));
				j=((j+1) % worldSize);
				break;
			}
		}
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
		if(sendPop(split[i], i) == false)
			cout << mark() << "Could not send population "<< i << " during fitcell!";
	}

	cout << mark() << "Performing fitcell" << endl;

	//send order to perform fitcell
	for(int i = 1; i < worldSize; i++) {
		sendIns(DoFitcell, i);
		//ownerlist[i] = i;
	}
	//do the local fitcell
	split[0].runFitcell(hostlist[0].threads);
	cout << mark() << "Server fitcell finished" << endl;

	//this region coudl results in a deadlock.....must be careful
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
		if(recvPop(&split[i], i) == false)
			cout << mark() << "Could not receive population "<< i << " during fitcell!";;
		//cout << "i:" << i << endl;
	}

	pop.clear();
	for(int i = 0; i < worldSize; i++)
		pop.addIndv(split[i]);
	split.clear();

	cout << mark() << "Fitcell complete" << endl;

}


void GASP2control::server_qe(GASP2pop &pop, int step) {
		int completed_structs;
		int launched_structs;
		MPI_Status m;
		Instruction recv, send;
		auto restart_timer = chrono::steady_clock::now();

		cout << mark() << "Beginning QE evaluation" << endl;
		pop.symmsort();
		int avail = 0;
		//var to store reference index of sent structure

		future<bool> serverthread;

		vector<GASP2pop> localpops(worldSize);
		vector<int> evaldind(worldSize);
		GASP2pop evald = pop;
		bool idle = false;
		save_state=false;

		ownerlist[0] = IDLE;
		eval_mut.unlock();


		vector<int> evallist(pop.size());
		cout << "evallist size " << evallist.size() << endl;


		for(int i = 0; i < evallist.size(); i++)
			evallist[i] = 0;

		cout << mark() << " evallist: ";
		for(int i = 0; i < evallist.size(); i++)
			cout << evallist[i] << " ";
		cout << endl;

		//for(int launched=0,completesteps=0; ; ) {
		while(true) {
			//cout << mark() << "start of launch block" << endl;
			bool queue_restart_write = false;
			bool complete_flag = false;


//AML: I am eliminating serverside threads for QE for the time being

			//test node 0
//						if(ownerlist[0] == 0) {
//							if(serverthread.valid() && serverthread.wait_for(timeout)==future_status::ready) {
//								//eval_mut.unlock();
//								serverthread.get();
//								ownerlist[0]=IDLE;
//								queue_restart_write=true;
//								completesteps++;
//							}
//						}
//						//save node 0
//						if( save_state && longeval_mut.try_lock()) {
//							if(chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - restart_timer).count() >= 180) {
//								cout << "queueing intermediate pop info on client " << ID << endl;
//								evald.mergeIndv(localpops[0],evaldind[0]);
//								complete_flag = true;
//								queue_restart_write = true;
//								save_state=false;
//							}
//							longeval_mut.unlock();
//						}



			//establish completion state
			for(int i = 1; i < worldSize; i++) {
				if(ownerlist[i] == i)
					sendIns(Ping, i);
			}
			this_thread::sleep_for(chrono::seconds(1));
			for(int i = 1; i < worldSize; i++) {
				if(ownerlist[i] == i) {
					int recvd;
					MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					while(recvd){
						recv = None;
						recvIns(recv, i);

						//cout << mark() << "server QE received from " << i << ": " << recv << endl;
							if(recv & Busy) {
								//cout << "marking busy " << i << endl;
								ownerlist[i] = i;
							}
							else if(recv & Complete) {
								cout << mark() << "Client " << i << " finished eval" << endl;
								ownerlist[i] = GETPOP;
								complete_flag = true;
							}
						MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					}
				}
			}


			//receive block
			if( chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - restart_timer).count() >= 180 || complete_flag) {
				for(int i = 1; i < worldSize; i++) {
					if(ownerlist[i] == GETPOP || ownerlist[i] == i) {
						sendIns(SendPop, i);
						this_thread::sleep_for(chrono::seconds(1)); //HACKHACK this is totally evil and I am an evil person for putting it here
						if(recvPop(&localpops[i], i)) {
							evald.mergeIndv(localpops[i],evaldind[i]);
							queue_restart_write = true;
							if(ownerlist[i] == GETPOP) {
								ownerlist[i] = IDLE;
								evallist[evaldind[i]] = 2;
							}
						}
						else {
							cout << mark() << "Server receive failed on client " << i << endl;
							if(ownerlist[i] == GETPOP) {
								cout << mark() << "A final structure was lost on " << ID << "..." << endl;
							}
							else {
								cout << mark() << "An intermediate send failed on " << ID << "..." << endl;
							}
							evallist[evaldind[i]] = 0;
							ownerlist[i] = DOWN;
						}
					}
				}
			}

			//check the down nodes and re-idle them if they are okay
			for(int i = 1; i < worldSize; i++) {
				if(ownerlist[i] == DOWN) {
					//cout << mark() << "Testing " << i << "..." << endl;
					sendIns(Ping, i);
					int recvd;
					MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					while(recvd){
						recv = None;
						recvIns(recv, i);
						if(recv & Ackn) {
							ownerlist[i] = IDLE;
							//break;
						}
						MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
					}
				}
			}

			//output temporary save file
			if(queue_restart_write) {
				GASP2pop restart = evald;
				//restart.energysort();
//				if(params.outputmode == "xml") {
//					writePopXML(restart, "restart", step);
//				}
				db.update(restart, "structs");
				restart_timer = chrono::steady_clock::now();
				queue_restart_write = false;
			}


			//cleanup node control lists
			for(int i = 0; i < worldSize; i++) {
				if(ownerlist[ownerlist[i]] == IDLE || ownerlist[ownerlist[i]] == DOWN)
					ownerlist[i] = IDLE;
			}
			avail = 0;
			for(int i = 0; i < worldSize; i++)
				if(ownerlist[i] == IDLE)
					avail++;




			//if there are clients available, try three times to launch a job
			for(int n = 0; n < 3; n++) {
				if(avail > 0) {

					cout << mark() << " evallist: ";
					for(int i = 0; i < evallist.size(); i++)
						cout << evallist[i] << " ";
					cout << endl;

					int ind = -1;
					for(int i = 0; i < evallist.size();i++) {
						if(evallist[i] == 0)  {
							ind = i;
							break;
						}
					}
					if(ind == -1) break;

					cout << mark() << "Evaluating structure " << ind << ", ID " << pop.indv(ind)->getID().toStr() << endl;
					cout << mark() << "structures avail: " << avail << endl;
					cout << mark() << " ownerlist: ";
					for(int i = 0; i < worldSize; i++)
						cout << ownerlist[i] << " ";
					cout << endl;

					//see how many nodes are needed
					int given;
					int needed;
					int next_needed;

					//we need as many nodes as there are symmops
					//if P1, we take two nodes to avoid serverthread
					if(pop.indv(ind)->getSymmcount() == 1)
						needed = 2;
					else
						needed = pop.indv(ind)->getSymmcount();

					//we see how many nodes the next structure needs
					int next_ind = -1;
					for(int i = ind+1; i < evallist.size();i++) {
						if(evallist[i] == 0)  {
							next_ind = i;
							break;
						}
					}

					if(next_ind >= 0) {
						if(pop.indv(next_ind)->getSymmcount() == 1)
							next_needed = 2;
						else
							next_needed = pop.indv(next_ind)->getSymmcount();
					}
					else
						next_needed = 0;

					//in the event we do not have enough nodes
					//to fulfill the symmetry needs
					if(needed > worldSize)
						given = avail;

					if(next_needed > (avail-needed))
						given = avail;
					else if (next_needed == 0)
						given = avail;
					else
						given = needed;

					//limit the maximum number of nodes
					//to eight, because of efficiency concerns
					if(given > 8)
						given = 8;

					cout << mark() << "needed is " << needed << endl;
					cout << mark() << "next_needed is " << next_needed << endl;
					cout << mark() << "given is " << given << endl;


					//get the slots
					vector<int> slots;
					slots.clear();
					int first = IDLE;
					for(int i = worldSize-1; i >= 0; i--) {
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
//								if(first == IDLE)
//									break;

					avail -= given;
					cout << mark() << "avail is " << avail << endl;

					string hostfile;
					hostfile.clear();
					hostfile = makeMachinefile(slots);

					if(first == 0) {
						//writeHost(localmachinefile, hostfile);
					}
					else {
						sendIns(GetHost, first);
						if(!sendHost(hostfile, 0, first)) {
							cout << mark() << "Sever host send failed" << endl;
							ownerlist[first] = DOWN;
						}
					}

					if(ownerlist[first] == first) { //prevents launch if first node goes DOWN
						evaldind[first] = ind;
						localpops[first] = pop.subpop(ind,1);

						if(first == 0) {
//										if(eval_mut.try_lock()) {
//											eval_mut.unlock();
//											if(!serverthread.valid()) {
//												cout << mark() << "Launching QE on server with "<< given << "nodes..."<< endl;
//												serverthread = std::async(launch::async, &GASP2control::runEvals, this, DoQE, &localpops[first], localmachinefile);
//												launched++;
//											}
//											else {
//												cout << mark() << "serverthread stuck" << endl;
//											}
//
//										}
//										else {
//											cout << mark() << "server thread is locked! marking DOWN" << endl;
//											ownerlist[0] = DOWN;
//										}
							cout << mark() << "first is 0! something went wrong..." << endl;
							ownerlist[0] = IDLE;
						}
						else {
							cout << mark() << "Launching QE on client " << first << "with "<< given << "nodes..." << endl;
							sendIns(PopAvail, first);
							if(!sendPop(localpops[first], first)) {
								cout << mark() << "Server send failed" << endl;
								ownerlist[first] = DOWN;
								evallist[ind] = 0;
							}
							else {
								sendIns(DoQE, first);
								//launched++;
								evallist[ind] = 1;
							}
						}
					}

					launched_structs = 0;
					completed_structs = 0;
					for(int i = 0; i < evallist.size(); i++) {
						if( evallist[i] == 1 ) launched_structs++;
						if( evallist[i] == 2 ) completed_structs++;
					}
					launched_structs += completed_structs;

					cout << mark() << launched_structs  << " launched" << endl;
					cout << mark() << completed_structs  << " completed" << endl;

					cout << mark() << " ownerlist: ";
					for(int i = 0; i < worldSize; i++)
						cout << ownerlist[i] << " ";
					cout << endl;

				} //if avail structures and not all structures launched



			} //for three steps



			launched_structs = 0;
			completed_structs = 0;
			for(int i = 0; i < evallist.size(); i++) {
				if( evallist[i] == 1 ) launched_structs++;
				if( evallist[i] == 2 ) completed_structs++;
			}
			launched_structs += completed_structs;


			idle = true;
			for(int i = 0; i < worldSize; i++) {
				if(ownerlist[i] != IDLE && ownerlist[i] != DOWN)
					idle = false;
			}


			cout << flush;
			if( ( completed_structs >= pop.size() ) && idle)
				break;
			this_thread::sleep_for(chrono::seconds(3));

		}

		cout << mark() << "Sorting after QE" << endl;
		evald.energysort();
//		if(params.outputmode == "xml") {
//			writePopXML(evald, "evald", step);
//		}
		db.update(evald, "structs");
		pop = evald;

		cout << mark() << "End of QE evalution" << endl;

	//end QE eval block

}


void GASP2control::ownerlist_update() {

	MPI_Status m;
	Instruction recv, send;

	for(int i = 1; i < worldSize; i++) {
		//ownerlist[i] = IDLE;

		//clear the message queue for safety
		int recvd;
		recv = None;
		int queuelength = 0;
		MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		while(recvd){
			recvIns(recv, i);
			queuelength++;
			MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		}
		cout << mark() << "ID " << ID << ": messages were in the queue: " << queuelength << endl;

		//try to get a signal for ~5 seconds
		sendIns(Ping, i);
		cout << mark() << "Sent initial ping to client " << i << endl;
	}
	this_thread::sleep_for(chrono::seconds(5));
	for(int i = 1; i < worldSize; i++) {
		int recvd;
		recv = None;
		ownerlist[i] = DOWN;
		MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		while(recvd){
			recvIns(recv, i);
			if(recv & Ackn) {
				cout << mark() << "Client " << i << " ack'd" << endl;
				ownerlist[i] = IDLE;
				//break;
			}
			MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		}


		if(ownerlist[i] == DOWN) {
			cout << mark() << "Client " << i << " is down!" << endl;

		}

	}

}


void GASP2control::server_randbuild(int gen) {
	cout << mark() << "Generating new randpop" << endl;
	randpop.clear();
	int factor = 1;
	if(params.spacemode != Spacemode::Single) {
		factor = 4;
	}
	randpop.init(root, params.popsize*factor, params.spacemode, params.group);
	randpop.setGen(gen);
}

// STATS CODE
//				//record statistics about spacegroup generation
//				vector<int> stats;
//				stats.resize(231);
//				for(int i = 0; i < stats.size(); i++)
//					stats[i] = 0;
//				for(int i = 0; i < evalpop.size(); i++)
//					stats[evalpop.indv(i)->getSpace()]+=1;
//				ofstream statfile;
//				stringstream statname;
//				statname << params.outputfile << "_stats_" << setw(3) << setfill('0') << step;
//				statfile.open(statname.str().c_str(), ofstream::out);
//				for(int i = 1; i < stats.size(); i++)
//					statfile << i << " " << stats[i] << endl;
//				statfile.close();


GASP2pop GASP2control::server_popbuild(int gen) {

	int replace = static_cast<int>(static_cast<double>(params.popsize)*params.replacement);

	GASP2pop outpop;

	//used to remove duplicates in clustering
	GASP2param specialp = params;
	specialp.clusterdiff = 0.9;


	//elitism, equivalent to MGAC1
	if(params.type == "elitism") {

		outpop = randpop;

		for(int i = 0; i < params.binlimit; i++) {

			bins[i].scale(params.const_scale, params.lin_scale, params.exp_scale);
			outpop.addIndv(bins[i].newPop(replace, params.spacemode));
			outpop = outpop.symmLimit(params.symmlimit);
		}
	}
	//full cross with random pop insertion
	else if (params.type == "fullcross" ) {

		outpop = randpop;
		outpop.addIndv( randpop.fullCross(params.spacemode) );

		for(int i = 0; i < params.binlimit; i++) {
			outpop.addIndv(bins[i].fullCross(params.spacemode));
			outpop.addIndv(bins[i].fullCross(params.spacemode, rootpop));
			outpop = outpop.symmLimit(params.symmlimit);
		}

	}
	//this algorithm is used in the precluster
	//generally speaking, the idea is to obtain a set of volume limited structures
	//that are SIGNIFICANTLY DIFFERENT from each other, such that every structure
	//in the cluster set has less than some percentage difference with all other structures
	else if (params.type == "precluster") {

		outpop = randpop;
		randpop.scale(1.0,0.0,0.0);
		//if(params.spacemode == Spacemode::Single) {
		//	outpop.addIndv( randpop.newPop(replace,params.spacemode) );
		//}
		//else{
		outpop.addIndv( randpop.fullCross(params.spacemode) );
		//}

		//in single spacegroup mode, a crossing size of the population size squared is acceptable
		//because the dimensionality is lower. no limits are set for the same reason
//		if(params.spacemode == Spacemode::Single) {
//			if(clusters[0].size() < (params.popsize*2)) {
//			  outpop.addIndv(clusters[0].newPop(params.popsize*2,params.spacemode));
//			  outpop.addIndv(clusters[0].newPop(randpop,params.popsize*2,params.spacemode));
//			}
//		}
		//in multispacegroup,the popsize for each cluster needs to be limited because the population
		//expansion can be very large. once a cluster reached the limit, only new structures resulting
		//randpop crossings will be added.
//		else {
			//GASP2pop indivs[460];
			GASP2pop indivs[230];
			//int selfself[230],
			int selfrand[230], total = 0;
			for(int i = 0; i < 230; i++) {
				if(clusters[i].size() < (params.popsize)) {
					if(clusters[i].size() < std::floor(std::sqrt( params.popsize))) {
						//selfself[i]=clusters[i].size()*clusters[i].size();
						selfrand[i]=clusters[i].size()*randpop.size();
					}
					else {
						//selfself[i]=(replace);
						selfrand[i]=(replace);
					}
					if(selfrand[i] > (replace))
						selfrand[i] = (replace);
					//partials.addIndv(bins[i]);
				}
				else {
					//selfself[i] = 0;
					selfrand[i] = 0;
				}
				//total += (selfself[i] + selfrand[i]);
				total += selfrand[i];

			}

			//cout << mark() << " a" << endl;
			for(int i = 0; i < 230; i++) {
				//cout << mark() << "popbuild i: " << i << endl;
				//if(clusters[i].size() > 1)
					//indivs[i*2] = clusters[i].newPop(selfself[i],params.spacemode);
				indivs[i] = clusters[i].newPop(randpop,selfrand[i],params.spacemode);
			}
			//cout << mark() << " b" << endl;
			outpop.reserve(total);
			//cout << mark() << " c" << endl;
			for(int i = 0; i < 230; i++) {
				outpop.addIndv(indivs[i]);
			//	outpop.addIndv(indivs[i*2+1]);
			}
			//cout << mark() << " d" << endl;

		//}

		outpop = outpop.symmLimit(params.symmlimit);
		outpop = outpop.spacebinUniques(hostlist[0].threads, clusters, params);

	}
	//this mode is used for clustered structure generation.
	else if (params.type == "clustered") {

		outpop = randpop;
		//outpop.addIndv( randpop.fullCross(params.spacemode) );
		if(params.spacemode == Spacemode::Single) {
			outpop.addIndv( randpop.newPop(replace,params.spacemode) );
		}
		else{
			outpop.addIndv( randpop.fullCross(params.spacemode) );
		}
		randpop.scale(1.0,0.0,0.0);

		for(int i = 0; i < params.binlimit; i++) {
			bins[i].scale(params.const_scale, params.lin_scale, params.exp_scale);
			outpop.addIndv(bins[i].newPop(replace,params.spacemode));
			outpop.addIndv(bins[i].newPop(randpop,replace,params.spacemode));
		}

		outpop = outpop.symmLimit(params.symmlimit);
		//outpop = outpop.spacebinUniques(hostlist[0].threads, bins, specialp);

	}
	else if(params.type == "finaleval") {
		outpop = rootpop.completeCheck();
		outpop = outpop.symmLimit(params.symmlimit);
		outpop.runSymmetrize(hostlist[0].threads);
	}

	outpop = outpop.symmLimit(params.symmlimit);
	outpop.setGen(gen);
	return outpop;
}


void GASP2control::down_check() {

	MPI_Status m;
	Instruction recv, send;

	//re - check the down nodes and re-idle them if they are okay
	for(int i = 1; i < worldSize; i++) {
		if(ownerlist[i] == DOWN) {
			sendIns(Ping, i);
			int recvd;
			MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
			while(recvd){
				recv = None;
				recvIns(recv, i);
				if(recv != None) {
					ownerlist[i] = IDLE;
					//break;
				}
				MPI_Iprobe(i, CONTROL, MPI_COMM_WORLD, &recvd, &m);
			}
		}
	}

	cout << mark() << " ownerlist end: ";
	for(int i = 0; i < worldSize; i++)
		cout << ownerlist[i] << " ";
	cout << endl;

	int downcount = 0;
	// if there are more than 3 down nodes at the end of the process, then we kill
	for(int i = 0; i < worldSize; i++) {
		if(ownerlist[i] == DOWN) {
			downcount++;
		}
	}
	if(downcount >= params.downlimit) {
		cout << mark() << "At the end of the QE evaluation, there were down nodes. Killing this job..." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		exit(1);
	}

}

void GASP2control::server_popcombine(GASP2pop pop) {

	if(params.type == "elitism") {
		//pop.remIndv(replace);
		pop.spacebinV(bins, params.popsize);
		//writePop(evalpop, "final", step);
	}
	else if (params.type == "fullcross" || params.type == "clustered") {
		pop.spacebinV(bins, params.popsize);
	}
	else if (params.type == "precluster") {
		pop.spacebinCluster(hostlist[0].threads, bins, clusters, params);
	}

	if(params.spacemode != Spacemode::Single) {
		cout << mark() << "Best : ";
		for(int n = 0; n < params.binlimit; n++) {
			if(bins[n].size() > 0) {
				cout << spacegroupNames[bins[n].indv(0)->getSpace()] <<",";
			}
			else {
				cout <<"null,";
			}
		}
		cout << endl;
	}




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

	MPI_Barrier(MPI_COMM_WORLD);

	//generate the local machinefile info
	UUID u; u.generate();
	string localmachinefile = "/tmp/machinefile-";
	localmachinefile.append(u.toStr());

	//this vector indicates the owner of all nodes
	//so, if ownerlist[i] = 0, node i is owned by
	//the server thread
	//0-worldsize means the node is running
	//IDLE (-1) means the node is idle
	//DOWN (-2) means the node is in a bad state
	//and is not responding to pings
	//vector<int> ownerlist(worldSize);
	//this vector contains a ping test for each node

	//vector<bool> pinged(worldSize);
	//the polling time is 200 ms

	//initialize ownerlist
	ownerlist.resize(worldSize);
	for(int i = 0; i < worldSize; i++)
		ownerlist[i] = IDLE;


	GASP2pop evalpop, lastpop, bestpop;
	//vector<GASP2pop> bins(230), clusterbins(230);
	bins.resize(230); clusters.resize(230);
//	if(params.spacemode != Spacemode::Single) {
//		bins.resize(230); clusters.resize(230);
//	}
//	else {
//		bins.resize(1);	clusters.resize(1);
//		//VERY IMPORTANT! if binlimit is not set right
//		//all of the binning algorithms will segfault
//		params.binlimit=1;
//	}
	for(int i = 0; i < clusters.size(); i++)
		clusters[i].clear();
	for(int i = 0; i < bins.size(); i++)
		bins[i].clear();

	db.initTable("structs");
	db.initTable("badfitcell");

	setup_restart();


	//pick the calculation mode
	evalmode = None;
	if(params.calcmethod == "qe")
		evalmode = (Instruction) (evalmode | DoQE);
//	else if(params.calcmethod == "fitcell")
//		evalmode = (Instruction) (evalmode | DoFitcell);



	vector<future<bool>> popsend(worldSize);
	future<bool> eval;
	future<bool> writer;




////////////////////START PRECLUSTER////////////////////////
	if(params.precompute > 0 && restart.length() <= 0) {
		int replace = static_cast<int>(static_cast<double>(params.popsize)*params.replacement);
		cout << mark() << "Starting precluster" << endl;
		string temp = params.type;
		db.initTable("precluster");

		bestpop.clear();
		GASP2pop prerand;
		for(int pcstep = 0; pcstep < params.precompute; pcstep++) {
			cout << mark() << "Precluster step " << pcstep << endl;
			params.type = "precluster";
			GASP2pop precompute;
			server_randbuild(0);
			precompute = server_popbuild(0);
			server_fitcell(precompute);
			precompute = precompute.volLimit();
			server_popcombine(precompute);

			params.type = temp;

			//TODO: add bin stats collection
			for(int n = 0; n < clusters.size(); n++) {
				bestpop.addIndv(clusters[n]);
			}
			cout << mark() << "Cluster size " << bestpop.size() << endl;
//			if(params.outputmode == "xml") {
//				writePopXML(bestpop, "precluster", 0);
//			}
			db.create(bestpop, "precluster");
			bestpop.clear();


		}


		db.updateTime(time(0), startstep);

		//bins = clusters;
		GASP2pop pre;
		for(int i = 0; i < bins.size(); i++) {
			//we trim to keep things under control
			//clusters[i].addIndv(bins[i]);
			if(clusters[i].size() > replace)
				clusters[i].remIndv(clusters[i].size() - replace);
			pre.addIndv(clusters[i]);
			clusters[i].clear();
		}
		cout << mark() << "Precluster QE evaluation starting" << endl;
		db.create(pre,"structs");
		pre.runSymmetrize(hostlist[0].threads);
		server_qe(pre, 0);

		server_popcombine(pre);
		bestpop.clear();
//		for(int n = 0; n < bins.size(); n++) {
//			bestpop.addIndv(bins[n]);
//		}
////		if(params.outputmode == "xml") {
////			writePopXML(bestpop, "precluster-eval", 0);
////		}
//		bestpop.clear();

		cout << mark() << "Precluster finished" << endl;

	}

////////////////////END PRECLUSTER////////////////////////


////////////////////START GENERATION EVALS/////////////////////////
	if(params.mode == "stepwise") {
		//cout << mark() << "STARTING STEP: " << startstep << endl;
		for(int step = startstep; step < params.generations; step++) {

			db.updateTime(time(0), step);

			GASP2pop bad, restart, good, tempstore;

			cout << mark() << "Starting new generation " << step << endl;

			//we always reinitialize the rootpop at each generation
			server_randbuild(step);
//			if(params.outputmode == "xml") {
//				writePopXML(randpop, "rand", step);
//			}

			//scale, cross and mutate
			cout << mark() << "Scaling and Crossing..." << endl;
			evalpop.clear();


			evalpop.addIndv(server_popbuild(step));

			if(params.type == "finaleval") {
				cout << mark() << "Performing final evaluation" << endl;
				good = evalpop;
			}
			else {
				cout << mark() << "Mutating..." << endl;
				evalpop.mutate(params.mutation_prob, params.spacemode);

				good.clear(); bad.clear();
				good = evalpop.symmLimit(bad, params.symmlimit);
				//good.addIndv(partials);
				//good.symmsort();

				server_fitcell(good);

				cout << mark() << "Sorting" << endl;
				evalpop = good;

				good.clear(); bad.clear();
				good = evalpop.volLimit(bad);
				db.create(good, "structs");
				//db.create(bad, "badfitcell"); //this was a bad idea, too much data
			}

			cout << mark() << "Candidate popsize:" << good.size() << endl;

			good.volumesort();
			if(params.outputmode == "xml") {
				writePopXML(good, "vollimit", step);
			}

			if(good.size() > 0) {

				if(params.type != "clustered" && params.type != "fitcell") {
					ownerlist_update();
				}
				if(params.calcmethod == "qe") {
					server_qe(good, step);
					cout << mark() << "Processing outliers" << endl;
					bad.clear();
					GASP2pop ok;
					GASP2pop bad = good.outliers(ok);
					if(bad.size() > 0)
						server_qe(bad,step);
					good=ok;
					good.addIndv(bad);
				}
				evalpop = good;


			} //if good size > 0
			else {
				cout << mark() << "No structures to be energy evaluated, continuing..." << endl;

			}

			//sort and reduce
			cout << mark() << "Energy sort" << endl;
			evalpop.energysort();

			if(params.type == "finaleval") {
//				cout << mark() << "Writing final output" << endl;
//				if(params.outputmode == "xml") {
//					writePopXML(evalpop, "fulleval", step);
//				}
				cout << mark() << "Final evaluation completed and written, exiting..." << endl;
				break;
			}
			else {
				cout << mark() << "Recombining populations" << endl;
				server_popcombine(evalpop);
//				bestpop.clear();
//				for(int n = 0; n < bins.size(); n++) {
//					bestpop.addIndv(bins[n]);
//				}
//				if(params.outputmode == "xml") {
//					writePopXML(bestpop, "final", step);
//				}
//				bestpop.clear();
			}


			//do a final check to see if the nodes are down
			down_check();

			db.updateTime(time(0), step);

		}

	}
	else if(params.mode == "steadystate") {
		cout << "Steadystate not yet implemented" << endl;
	}
////////////////////END GENERATION EVALS/////////////////////////

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
//	cout << "runeval evalmode " << (int) i << endl;
	//only one of these will execute
	if(i & DoFitcell) {
		p->runFitcell(nodethreads);
		completed = true;
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

	MPI_Barrier(MPI_COMM_WORLD);

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
	chrono::milliseconds t(100);

	bool queueSend = false;
	bool queueEval = false;
	bool queueRecv = false;
	bool queueHost = false;
	int target = 0;

	auto start = chrono::steady_clock::now();

	string hostfile;
	int itemp;

	save_state=false;
	completed = false;

	while(true) {
		 ack = None;

		//cout << "client cycle " << this->ID << endl;

		MPI_Iprobe(0, CONTROL, MPI_COMM_WORLD, &recvd, &m);
		while(recvd) {
			i = None;
			recvIns(i, 0);

			//cout << mark() << "client " << this->ID << " received: " << i << endl;
			if(i & Shutdown) {
				remove(localmachinefile.c_str());
				return;
			}

			if(i & Ping){
				ack = (Instruction) (ack | Ackn);
			}

			if(i & GetHost) {
				if(recvHost(hostfile, itemp, 0))
					writeHost(localmachinefile, hostfile);
				else
					remove(localmachinefile.c_str());
			}

			if(i & PopAvail)
				queueRecv = true;

			if(i & SendPop)
				queueSend = true;

			//do not allow the evalmode to change
			if( (i & (DoFitcell | DoCharmm | DoQE | DoCustom)) && !queueEval) {
				evalmode = i;
				//cout << "evalmode " << (int) evalmode << endl;
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

			if(ack > None) {
				//cout << "client " << ID << " sent ack " << ack << endl;
				sendIns(ack,0);
			}

			MPI_Iprobe(0, CONTROL, MPI_COMM_WORLD, &recvd, &m);

		}

		//execute work
		//A SEND ALWAYS TAKES PRECEDENCE OVER A RECEIVE

		if(queueSend) {
			//see what the current evalmode
			//cout << "client " << ID << ": send queued" << endl;
			if(evalmode & DoFitcell) {
				if(eval_mut.try_lock()) {
					if(!sendPop(evalpop, 0))
						cout << mark() << "Sendpop DoFitcell " << ID << " failed!" << endl;;
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
					if(!sendPop(evalpop, 0))
						cout << mark() << "Sendpop DoQEev " << ID << " failed!" << endl;
					queueSend = false;
					eval_mut.unlock();
				}
				else {
					//in this instance the client sends info
					if(!sendPop(temppop, 0))
						cout << mark() << "Sendpop DoQE " << ID << " failed!" << endl;;
					queueSend = false;
				}

			}
		}

		if(queueRecv) {
			//cout << "client " << ID << ": recv queued" << endl;
			if(!queueSend) {
				if(eval_mut.try_lock()) {
					eval_mut.unlock();
					if(!recvPop(&evalpop, 0))
						cout << mark() << "Recvpop " << ID << " failed!" << endl;;
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
			//cout << "client " << ID << ": eval queued " << evalmode << endl;
			if(eval_mut.try_lock()) {
				eval_mut.unlock();
				temppop = evalpop;
				if(!eval.valid()) {
					eval = std::async(launch::async, &GASP2control::runEvals, this, evalmode, &evalpop, localmachinefile);
					this_thread::sleep_for(t);
					queueEval = false;
				}
				else {
					cout << mark() << "thread " << ID << " is not valid?" << endl;
				}
			}
			else {
				cout << mark() << "Eval thread is locked, retrying in ~200 ms..." << ID << endl;
			}
		}

		//store the temppop
		if( save_state && longeval_mut.try_lock()) {
			//cout << "queueing intermediate pop info on client " << ID << endl;
//					cout << mark() << "gasp client pop energies" << endl;
//					for(int i = 0; i < evalpop.size(); i++)
//						cout << "energy: " << evalpop.indv(i)->getEnergy() << endl;
			temppop = evalpop;
			//sendIns(PopAvail, 0);
			save_state=false;
			longeval_mut.unlock();
		}


		if(eval.valid() && eval.wait_for(timeout)==future_status::ready){
			eval.get();
			//if(completed == false) {
				//cout << mark() << "Client " << ID << ": weird exit on thread occurred..." << endl;
			completed = true;
			//}
		}
		if(completed) {

			sendIns(Complete, 0);
			completed = false;
		}
		cout << flush;
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
	MPI_Request m;

	cout << mark() << " Sendhost to " << target << " launched, ID " << ID << ", size: " << t <<  endl;
	//cout << mark() << "hostinfo sender size " << t << endl << endl << host << endl << endl;

	//send size info
	MPI_Issend(&t, 1, MPI_INT, target,HOSTS1,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;

	//send info
	MPI_Issend(host.c_str(),t,MPI_CHAR,target,HOSTS2,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;
	//send proc info



	MPI_Issend(&procs, 1, MPI_INT, target, HOSTS3, MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;

	cout << mark() << " Sendhost " << target << " finished, ID " << ID << endl;

	return true;
}

bool GASP2control::recvHost(string &host, int &procs, int target) {
	int ierr;
	int v = 0;
	MPI_Request m;
	//recv size info
	MPI_Irecv(&v, 1, MPI_INT, target,HOSTS1,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;

	char * buff = new char[v+1];
	buff[v] = '\0';
	//cout << mark() << "hostsize:" << v << endl;
	cout << mark() << " Recvhost from " << target << " launched, ID " << ID << ", size: " << v << endl;
	//recv hostname
	MPI_Irecv((void *) buff, v, MPI_CHAR, target,HOSTS2,MPI_COMM_WORLD,&m);
	if(!testReq(m, 120))
		return false;
	host = buff;
	//recv proc count
	//cout << mark() << "hostinfo recveiver" << endl << endl << host << endl << endl;

	MPI_Irecv(&procs, 1, MPI_INT, target, HOSTS3, MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;
	//cout << "procs:" << procs << endl;
	delete [] buff;
	cout << mark() << " Recvhost " << target << " finished, ID " << ID << endl;
	return true;
}

//this tests a request for completion every
//five seconds for up to a minute. if the request
//fails to complete in that time, it is cancelled.
bool GASP2control::testReq(MPI_Request &m, int t) {
	int ok;

	MPI_Test(&m, &ok, MPI_STATUS_IGNORE);
	for(int i = 0; i < t*10; i++) {
		if(ok==true) break;
		this_thread::sleep_for(chrono::milliseconds(100));
		MPI_Test(&m, &ok, MPI_STATUS_IGNORE);
	}
	if(ok==false) {
		cout << mark() << "A testReq failed on ID " << ID << endl;
		MPI_Cancel(&m);
		MPI_Request_free(&m);
		return false;
	}
	return true;

}

bool GASP2control::sendPop(GASP2pop p, int target) {
	string pop = p.saveXML();
	//cout << "Pop" << endl << pop << endl << endl
	MPI_Request m;

	pop += "\0";
	int t = pop.length();

	//cout << mark() << " Sendpop to " << target << " launched, ID " << ID << ", size: " << t <<  endl;
	//send size info
	MPI_Issend(&t, 1, MPI_INT, target,POP1,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;


	//send info
	MPI_Issend(pop.c_str(),t,MPI_CHAR,target,POP2,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;

	//cout << mark() << " Sendpop " << target << " finished, ID " << ID << endl;

	return true;
}

bool GASP2control::recvPop(GASP2pop *p, int target) {
	string pop;
	int v = 0;
	MPI_Request m;
	bool ok;

	//recv size info
	MPI_Irecv(&v, 1, MPI_INT, target,POP1,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;
	//cout << mark() << " Recvpop from " << target << " launched, ID " << ID << ", size: " << v << endl;
	//char * buff = new char[v+1];
	//recv hostname
	pop.resize(v+1);

	MPI_Irecv((void *) pop.data(), v, MPI_CHAR, target,POP2,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;

	//buff[v] = '\0';
	//pop = buff;

	//cout << "Pop" << endl << pop << endl << endl;

	//delete [] buff;

	string error;
	p->parseXML(pop,error);

	//cout << "Error" <<  error << endl;

	//cout << mark() << " Recvpop " << target << " finished, ID " << ID << ", structsize: " << p->size() << endl;

	return true;
}


bool GASP2control::sendIns(Instruction i, int target) {
	//cout << "sending ins to target " << target << endl;

	MPI_Request r;
	MPI_Issend(&i, 1, MPI_INT, target, CONTROL, MPI_COMM_WORLD, &r);
	//if(!testReq(r, 120))
	//	return false;
	//MPI_Request_free(&r);
	return true;
}

bool GASP2control::recvIns(Instruction &i, int target) {
	//cout << "receiving ins from target " << target << endl;
	MPI_Request r;
	MPI_Irecv(&i, 1, MPI_INT, target, CONTROL, MPI_COMM_WORLD, &r);
	if(!testReq(r, 120))
		return false;
	//cout << "procs:" << procs << endl;
	//MPI_Request_free(&r);
	return true;
}


string GASP2control::makeMachinefile(vector<int> slots) {
	stringstream out;
	out.str("");

	for(int i = 0; i < slots.size(); i++) {
		if(slots[i] < hostlist.size()) {
			for(int j = 0; j < (hostlist[slots[i]].threads - 1); j++)
				out << hostlist[slots[i]].hostname << endl;
		}
	}
	return out.str();
}


bool GASP2control::writePopXML(GASP2pop pop, string tag, int step) {

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
		outf << data.c_str(); // <<endl;
		outf.close();
}
