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
GASP2control::GASP2control(time_t start, int size, string input, string restart, int _startstep) {
	starttime = start;
	startstep = _startstep;
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
			rootpop.energysort();
			rootpop.remIndv(rootpop.size() - params.popsize);

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
	hostname[1023] = '\0';
	gethostname(_host,1023);
	hostname = _host;
	cout << "Hostname ID " << ID << ":" << hostname << endl;

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

	//vector<bool> pinged(worldSize);
	//the polling time is 200 ms
	chrono::milliseconds t(200);
	chrono::milliseconds timeout(0);

	//initialize ownerlist
	for(int i = 0; i < worldSize; i++)
		ownerlist[i] = IDLE;

//	//establish who is paying attention initially
//	for(int i = 1; i < worldSize; i++) {
//		//try to get a signal for ~5 seconds
//		sendIns(Ping, i);
//		cout << mark() << "Sending initial ping to client " << i << endl;
//	}
//	this_thread::sleep_for(chrono::seconds(1));
//	for(int i = 1; i < worldSize; i++) {
//		recvIns(recv, i);
//		if(recv & Ackn) {
//			cout << mark() << "Client " << i << " ack'd" << endl;
//			ownerlist[i] = IDLE;
//		}
//		else {
//			cout << mark() << "Client " << i << " is down!" << endl;
//			ownerlist[i] = DOWN;
//		}
//
//	}


	GASP2pop evalpop, lastpop, bestpop, partials;
	vector<GASP2pop> bins(230);

	vector<GASP2pop> split;

	if(restart.length() > 0) {
		lastpop = rootpop;
		rootpop.clear();
		cout << mark() << "Starting from restart population \"" << restart << "\"\n";
		if(params.spacemode != Spacemode::Single) {
			//bestpop = lastpop;
			lastpop.spacebinV(bins, params.popsize);
		}
	}
	else {
		lastpop.init(root, params.popsize, params.spacemode, params.group);
	}

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



	vector<future<bool>> popsend(worldSize);
	future<bool> eval;
	future<bool> writer;
	vector<int> stats;

	//auto restart_timer;


	MPI_Status m;


	if(params.mode == "stepwise") {
		for(int step = startstep; step < params.generations; step++) {


////////////////////START GEN BUILD/////////////////////////
			GASP2pop bad, restart, good;


			cout << mark() << "Starting new generation " << step << endl;

			//we always reinitialize the rootpop at each generation
			//this is designed to precent
			cout << mark() << "Generating new rootpop..." << endl;
			rootpop.clear();
			rootpop.init(root, params.popsize, params.spacemode, params.group);
			writePop(rootpop, "root", step);

			//scale, cross and mutate
			cout << mark() << "Scaling..." << endl;
			if(step > 0)
				lastpop.scale(params.const_scale, params.lin_scale, params.exp_scale);

			cout << mark() << "Crossing..." << endl;



			evalpop.clear();
			partials.clear();
			if(params.type == "elitism") {
				evalpop = lastpop.newPop(replace, params.spacemode);
				evalpop.addIndv(lastpop);
			}

			else if (params.type == "classic") {

				evalpop = rootpop;
				evalpop.addIndv( rootpop.fullCross(params.spacemode) );

				if(params.spacemode == Spacemode::Single || step == 0) {
					evalpop.addIndv(lastpop.fullCross(params.spacemode));
					evalpop.addIndv(lastpop.fullCross(params.spacemode, rootpop));
					//partials.addIndv(lastpop);
				}
				else {
					for(int i = 0; i < params.binlimit; i++) {
						evalpop.addIndv(bins[i].fullCross(params.spacemode));
						evalpop.addIndv(bins[i].fullCross(params.spacemode, rootpop));
						//partials.addIndv(bins[i]);
					}
				}


			}

			//partials = partials.completeCheck();

			cout << mark() << "Mutating..." << endl;
			evalpop.mutate(params.mutation_prob, params.spacemode);

			cout << mark() << "Crossover size: " << evalpop.size() << endl;

			stats.clear();
			stats.resize(231);
			for(int i = 0; i < stats.size(); i++)
				stats[i] = 0;
			for(int i = 0; i < evalpop.size(); i++)
				stats[evalpop.indv(i)->getSpace()]+=1;
			ofstream statfile;
			stringstream statname;
			statname << params.outputfile << "_stats_" << setw(3) << setfill('0') << step;
			statfile.open(statname.str().c_str(), ofstream::out);
			for(int i = 1; i < stats.size(); i++)
				statfile << i << " " << stats[i] << endl;
			statfile.close();

////////////////////END POP BUILD/////////////////////////


			good.clear();
			bad.clear();
			good = evalpop.symmLimit(bad, params.symmlimit);
			//good.addIndv(partials);

			split.clear();
			split.resize(worldSize);
			good.symmsort();

////////////////////START FITCELL/////////////////////////

			//distribute the structures evenly so
			//that one host doesn't get too overloaded.
			//now takes into account IDLE/DOWN state
			int j = 0;
			for(int i = 0; i < evalpop.size(); i++) {
				for( ; ; j=((j+1) % worldSize) ) {
					if(ownerlist[j] == IDLE) {
						split[j].addIndv(*good.indv(i));
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

			cout << mark() << "Fitcell complete" << endl;


////////////////////END FITCELL/////////////////////////

			//combine, then save fitcell
			cout << mark() << "Sorting" << endl;
			lastpop = evalpop;
			good.clear();
			for(int i = 0; i < worldSize; i++)
				good.addIndv(split[i]);
//			if (params.type == "classic")
//				good.addIndv(bad);

			evalpop = good;
			//evalpop.volumesort();
			//writePop(evalpop, "fitcell", step);

			//evaluate
			evalpop.symmsort();
			//GASP2pop bad, good, restart;

			good.clear();
			bad.clear();
			good = evalpop.volLimit(bad);
			//add the incomplete structures to the partials to
			//bypass the fitcell step


			cout << mark() << "Candidate popsize:" << good.size() << endl;

			good.volumesort();
			writePop(good, "vollimit", step);


			if(good.size() > 0) {
				future<bool> serverthread;



				//establish who is paying attention initially
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


////////////////////START QE DISPATCH/////////////////////////
				if(params.calcmethod == "qe") {

					int completed_structs;
					int launched_structs;


					auto restart_timer = chrono::steady_clock::now();

					cout << mark() << "Beginning QE evaluation" << endl;
					good.symmsort();
					int avail = 0;
					//var to store reference index of sent structure
					vector<GASP2pop> localpops(worldSize);
					vector<int> evaldind(worldSize);
					GASP2pop evald = good;
					bool idle = false;
					save_state=false;

					ownerlist[0] = IDLE;
					eval_mut.unlock();


					vector<int> evallist(good.size());
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

									cout << mark() << "server QE received from " << i << ": " << recv << endl;
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
								cout << mark() << "Testing " << i << "..." << endl;
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
							restart = evald;
							restart.energysort();
							writePop(restart, "restart", step);
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

								cout << mark() << "Evaluating structure " << ind << endl;
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
								if(good.indv(ind)->getSymmcount() == 1)
									needed = 2;
								else
									needed = good.indv(ind)->getSymmcount();

								//we see how many nodes the next structure needs
								int next_ind = -1;
								for(int i = ind+1; i < evallist.size();i++) {
									if(evallist[i] == 0)  {
										next_ind = i;
										break;
									}
								}

								if(next_ind >= 0) {
									if(good.indv(next_ind)->getSymmcount() == 1)
										next_needed = 2;
									else
										next_needed = good.indv(next_ind)->getSymmcount();
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
								//to 2x the symmlimit
								if(given > (params.symmlimit * 2))
									given = params.symmlimit * 2;

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
									localpops[first] = good.subpop(ind,1);

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
						if( ( completed_structs >= good.size() ) && idle)
							break;
						this_thread::sleep_for(chrono::seconds(3));

					}

					cout << mark() << "Sorting after QE" << endl;
					evald.energysort();
					writePop(evald, "evald", step);
					good = evald;



					cout << mark() << "End of QE evalution" << endl;

				} //end QE eval block

				good.addIndv(bad);
				evalpop = good;

			} //if good size > 0
			else {
				cout << mark() << "No structures to be energy evaluated, continuing..." << endl;

			}



			//sort and reduce
			evalpop.energysort();
			if(params.type == "elitism") {
				evalpop.remIndv(replace);
				lastpop = evalpop;
				writePop(evalpop, "final", step);
			}
			else if (params.type == "classic") {

//				writePop(bestpop, "best", step);
				if(params.spacemode == Spacemode::Single) {
					lastpop.addIndv(evalpop);
					lastpop.energysort();
					lastpop.remIndv(lastpop.size()-params.popsize);
					writePop(lastpop, "final", step);
				}
				else {
					cout << mark() << "Sorting bins..." << endl;
					bestpop.clear();
					evalpop.spacebinV(bins, params.popsize);
					//TODO: add bin stats collection
					for(int n = 0; n < bins.size(); n++) {
						bestpop.addIndv(bins[n]);
					}

					cout << mark() << "Best spacegroups: ";
					for(int n = 0; n < params.binlimit; n++) {
						if(bins[n].size() > 0) {
							cout << spacegroupNames[bins[n].indv(0)->getSpace()] <<",";

						}
						else {
							cout <<"null,";
						}
					}
					cout << endl;
					//bestpop.addIndv(evalpop);
					//bestpop = bestpop.spacebin(params.popsize, 50);
					//lastpop = bestpop;
					//bestpop.energysort();
					writePop(bestpop, "final", step);
					bestpop.clear();
				}
//				evalpop = bestpop;
//				evalpop.energysort();
//				evalpop.remIndv(evalpop.size()-params.popsize);

			}


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
//

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
	chrono::milliseconds t(200);

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

			cout << mark() << "client " << this->ID << " received: " << i << endl;
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
			cout << "client " << ID << ": send queued" << endl;
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
					if(!sendPop(temppop, 0))
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
			cout << "client " << ID << ": recv queued" << endl;
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
			cout << "client " << ID << ": eval queued " << evalmode << endl;
			if(eval_mut.try_lock()) {
				eval_mut.unlock();
				temppop = evalpop;
				if(!eval.valid()) {
					eval = std::async(launch::async, &GASP2control::runEvals, this, evalmode, &evalpop, localmachinefile);
					this_thread::sleep_for(t);
					queueEval = false;
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
	for(int i = 0; i < t; i++) {
		if(ok==true) break;
		this_thread::sleep_for(chrono::seconds(1));
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

	cout << mark() << " Sendpop to " << target << " launched, ID " << ID << ", size: " << t <<  endl;
	//send size info
	MPI_Issend(&t, 1, MPI_INT, target,POP1,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;


	//send info
	MPI_Issend(pop.c_str(),t,MPI_CHAR,target,POP2,MPI_COMM_WORLD, &m);
	if(!testReq(m, 120))
		return false;

	cout << mark() << " Sendpop " << target << " finished, ID " << ID << endl;

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
	cout << mark() << " Recvpop from " << target << " launched, ID " << ID << ", size: " << v << endl;
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

	cout << mark() << " Recvpop " << target << " finished, ID " << ID << ", structsize: " << p->size() << endl;

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
		outf << data.c_str(); // <<endl;
		outf.close();
}
