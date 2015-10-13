#include "gasp2.hpp"
#include "optionparser.h"

using namespace std;

const string version = "2.0-alpha";

//front material for option parsing
struct Arg: public option::Arg
{
	  static option::ArgStatus Required(const option::Option& option, bool msg)
	  {
	    if (option.arg != 0)
	      return option::ARG_OK;

	    if (msg) cout << "Option '" << option.name << "' requires argument" << endl;
	    return option::ARG_ILLEGAL;
	  }
	  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
	  {
	    if (option.arg != 0 && option.arg[0] != 0)
	      return option::ARG_OK;

	    if (msg) cout << "Option '" << option.name << "' requires argument" << endl;
	    return option::ARG_ILLEGAL;
	  }
};

enum optIndex {INPUT,HELP,RESTART,STEP,SPACEGROUPS,CONVERT,COMBINE,SIZE,TEMPLATE,PLANE,RECLUSTER,THREADS };

const option::Descriptor usage[] =
{
		{INPUT,0,"i","input",Arg::Required,"-i,--input  The XML input file for running"},
		{HELP,0,"h","help",option::Arg::Optional,"-h,--help  Displays this help message"},
		{RESTART,0,"r","restart",Arg::NonEmpty,"-r,--restart  Optional argument for a restart file"},
		{STEP,0,"S","step",Arg::NonEmpty,"-S,--step Optional argument used to specify starting step"},
		{SPACEGROUPS,0,"l","spacegroups",Arg::Optional,"-l,--spacegroups  List valid spacegroups"},
		{CONVERT, 0, "c", "cif", Arg::NonEmpty, "-c, --cif  Name of output file for cif; takes the input (-i) and turns it into a cif"},
		{COMBINE, 0, "m","merge",Arg::NonEmpty, "-m, --merge Combine multiple files to form a single population file (comma delimited)"},
		{SIZE, 0, "s","size", Arg::NonEmpty, "-s, --size Used in conjunction with merge to denote the size of merged population"},
		{TEMPLATE, 0, "t","template",Arg::NonEmpty, "-t, --template An XML molecule template from a cif; if a plane is given then rotation and other values will be checked"},
		{PLANE, 0, "p", "plane",Arg::NonEmpty, "-p, --plane Specifies the three atom plane to be used for a template (comma delimited)"},
		{RECLUSTER, 0, "g", "recluster",Arg::NonEmpty, "-g, --recluster Reclusters a population using a standard input for input and the population to be reclustered as the g argument"},
		{THREADS, 0, "j", "threads",Arg::NonEmpty, "-j, --threads Number of threads to use (ie, clustering)"},
		{ 0, 0, 0, 0, 0, 0 }
};


int main( int argc, char* argv[] ) {

	int mpithreading;
	//MPI_Init(&argc,&argv);
	MPI_Init_thread(&argc,&argv, MPI_THREAD_FUNNELED, &mpithreading);
	if(mpithreading == MPI_THREAD_FUNNELED)
		cout << "Threading requested was given!" << endl;

	int size, ID;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);



    //option parsing
    argc-=(argc>0); argv+=(argc>0);
    option::Stats  stats(usage, argc, argv);
    option::Option options[stats.options_max], buffer[stats.buffer_max];
    option::Parser parse(usage, argc, argv, options, buffer);

    if(parse.error())
    	return 1;

    if(options[HELP] || argc == 0) {
    	option::printUsage(cout, usage);
    	return 0;
    }

    //spacegroup stuff
	if(!loadSpaceGroups()) {
		cout << "There was a problem loading the spacegroup file!\n";
		exit(1);
	}
	int threads=1;
	if(options[THREADS]) {
		threads = std::stoi(string(options[THREADS].arg));
		if(threads < 1) {
			cout << "Less than 1 thread specified, exiting..." << endl;
			exit(1);
		}
	}

    if(options[SPACEGROUPS]) {
    	cout << "Valid groups are the following (can use index or string names):" << endl << endl;
    	cout << " Index  Name" << endl;
    	for(int i = 1; i < spacegroupNames.size(); i++)
    		cout << " " << setw(4) << i << ":  " << spacegroupNames[i] << endl;

    	return 0;
    }

    if(options[COMBINE]) {
    	int size=300;
    	if(options[SIZE]) {
    		//TODO: Need to put something here to set size
    		size = std::stoi(string(options[SIZE].arg));
    	}

    	if(options[INPUT] && options[COMBINE].arg != NULL) {
    		GASP2pop temppop, finalpop;
    		ofstream outf;
    		string files = options[INPUT].arg;
    		vector<string> filelist=split(files,',');
    		string errorstring;
    		for(int i = 0; i < filelist.size(); i++) {
    			temppop.clear();
				tinyxml2::XMLDocument doc;
				doc.LoadFile(filelist[i].c_str());
				if(doc.ErrorID() == 0) {
					tinyxml2::XMLElement * pop = doc.FirstChildElement("mgac")->FirstChildElement("pop");
					if(!temppop.loadXMLrestart(pop, errorstring)) {
						cout << "There was an error in " << filelist[i] <<": " << errorstring << endl;
						exit(1);//MPI_Abort(1,MPI_COMM_WORLD);
					}
					finalpop.addIndv(temppop);
				}
				else {
					cout << "!!! There was a problem with opening input file \"" << filelist[i] << "\"!" << endl;

					cout << "Check to see if the file exists or if the XML file" << endl;
					cout << "is properly formed, with tags formatted correctly." << endl;
					cout << "Aborting... " << endl;
					exit(1);//MPI_Abort(1,MPI_COMM_WORLD);
				}
    		}
    		//remove(options[CONVERT].arg);
    		finalpop.energysort();
    		finalpop.dedup(size);
    		finalpop.runSymmetrize(2);

    		outf.open(options[COMBINE].arg, ofstream::out);
    		if(outf.fail()) {
    			cout << mark() << "ERROR: COULD NOT OPEN FILE FOR SAVING! exiting sadly..." << endl;
    			exit(1);
    		}
    		else {
    			outf << "<mgac>\n" << finalpop.saveXML() << endl << "</mgac>\n";
    			outf.close();
    		}

    		//finalpop.writeCIF(string(options[COMBINE].arg));
    		cout << mark() << "Files successfully merged and written" << endl;

    	}

    	return 0;
    }

    if(options[CONVERT]) {
    	if(options[CONVERT].arg == NULL) {
    		return 0;
    	}
    	if(options[INPUT] && options[CONVERT].arg != NULL) {

    		GASP2pop temppop;
    		string errorstring;
    		tinyxml2::XMLDocument doc;
    		doc.LoadFile(options[INPUT].arg);
    		cout << mark() << "xml loaded" << endl;
    		if(doc.ErrorID() == 0) {
    			tinyxml2::XMLElement * pop = doc.FirstChildElement("mgac")->FirstChildElement("pop");
    			if(!temppop.loadXMLrestart(pop, errorstring)) {
    				cout << "There was an error in the restart file: " << errorstring << endl;
    				exit(1);//MPI_Abort(1,MPI_COMM_WORLD);
    			}
    		}
    		else {
    			cout << "!!! There was a problem with opening the input file!" << endl;
    			cout << "Check to see if the file exists or if the XML file" << endl;
    			cout << "is properly formed, with tags formatted correctly." << endl;
    			cout << "Aborting... " << endl;
    			exit(1);MPI_Abort(1,MPI_COMM_WORLD);
    		}
    		remove(options[CONVERT].arg);
    		doc.Clear();
    		cout << mark() << "pop loaded" << endl;
    		//temppop.energysort();
    		//cout << mark() << "sorted" << endl;
    		temppop.runSymmetrize(threads);
    		cout << mark() << "symmed" << endl;
    		temppop.writeCIF(string(options[CONVERT].arg));
    		cout << mark() << "File successfully converted" << endl;
    	}
    	else {
    		cout << "Requires -i for the input file!" << endl;

    	}


    	return 0;
    }

    if(options[TEMPLATE]) {
    	if(options[INPUT] && options[TEMPLATE].arg != NULL) {
    		GASP2struct st;
    		cout << mark() << "converting" << endl;
    		if(options[PLANE]) {
    			st.readCifMol(options[INPUT].arg, options[TEMPLATE].arg, options[PLANE].arg);
    		}
    		else{
    			st.readCifMol(options[INPUT].arg, options[TEMPLATE].arg);
    		}


    	}
    	else {
    		cout << "Requires -i for the input file!" << endl;
    	}
    	return 0;
    }

    if(options[RECLUSTER]) {
    	if(options[RECLUSTER].arg == NULL) {
    		return 0;
    	}
    	if(options[INPUT] && options[RESTART].arg != NULL && options[RECLUSTER].arg != NULL) {

    		GASP2control client(string(options[INPUT].arg));
    		GASP2param p = client.getParams();
    		cout << "got one" << endl;
    		GASP2pop temppop;
    		string errorstring;
    		tinyxml2::XMLDocument doc;
    		doc.LoadFile(options[RESTART].arg);
    		if(doc.ErrorID() == 0) {
    			tinyxml2::XMLElement * pop = doc.FirstChildElement("mgac")->FirstChildElement("pop");
    			if(!temppop.loadXMLrestart(pop, errorstring)) {
    				cout << "There was an error in the restart file: " << errorstring << endl;
    				exit(1);//MPI_Abort(1,MPI_COMM_WORLD);
    			}
    		}
    		else {
    			cout << "!!! There was a problem with opening the input file!" << endl;
    			cout << "Check to see if the file exists or if the XML file" << endl;
    			cout << "is properly formed, with tags formatted correctly." << endl;
    			cout << "Aborting... " << endl;
    			exit(1);MPI_Abort(1,MPI_COMM_WORLD);
    		}

    		cout << mark() << "precluster" << endl;
    		GASP2pop clusters;
    		temppop.clusterReset();
    		temppop.cluster(clusters, p, threads);
    		cout << mark() << "cluster done, size: " << clusters.size() << endl;
    		//cout << mark() <<"group count: " << temppop.assignClusterGroups(p, threads) << endl;
    		cout << mark() << "writing dists" << endl;
    		temppop.allDistances(p, threads);
    		cout << mark() <<"done"<< endl;


    		ofstream outf;
    		outf.open(options[RECLUSTER].arg, ofstream::out);
    		if(outf.fail()) {
    			cout << mark() << "ERROR: COULD NOT OPEN FILE FOR SAVING! exiting sadly..." << endl;
    			exit(1);
    		}
    		else {
    			outf << "<mgac>\n" << temppop.saveXML() << endl << "</mgac>\n";
    			outf.close();
    		}

    		//temppop.writeCIF(string(options[RECLUSTER].arg));
    		cout << mark() << "File successfully converted" << endl;
    	}
    	else {
    		cout << "Requires -i for the input file!" << endl;
    	}


    	return 0;
    }



    if(!options[INPUT]) {
    	cout << "An input file is required! Use -i to specify." << endl;
    	option::printUsage(cout, usage);
    	return 0;
    }
    string infile;
    infile.assign(options[INPUT].arg);



    time_t prog_start;
    //startup control for server/client
    if(ID==0) {
    	prog_start = time(0);
    	char * dt = ctime(&prog_start);
    	cout << "------------------------------------------------------------" << endl;
    	cout << "MGAC v"<< version << " startup at " << dt;
    	cout << endl;
    	cout << "Users of this program should cite: " << endl;
    	cout << endl << endl;
    	cout << "------------------------------------------------------------" << endl << endl;
    	cout << "Using input file " << infile << endl << endl;

    	if( options[RESTART] && options[RESTART].arg != NULL) {
    		if(options[STEP] && options[STEP].arg != NULL) {
    			int steps;
    			stringstream stepconvert;
    			stepconvert << options[STEP].arg;
    			stepconvert >> steps;

    			GASP2control server(prog_start,size, infile, options[RESTART].arg, steps);
    			server.server_prog();
    		}
    		//ignore if steps isn't correct
    		else {
    			GASP2control server(prog_start,size, infile, options[RESTART].arg);
    			server.server_prog();
    		}
    	}
    	else {
    		GASP2control server(prog_start,size, infile);
    		server.server_prog();
    	}
    }
    else {
    	GASP2control client(ID, infile);
    	client.client_prog();
    }

    //program termination
    if(ID==0) {
    	time_t final = time(0);
    	char *dt = ctime(&final);
    	cout << endl;
    	cout << "------------------------------------------------------------" << endl;
    	cout << "MGAC completed at " << dt;
    	cout << "------------------------------------------------------------" << endl;
    }

    MPI_Finalize();
	return 0;
}
