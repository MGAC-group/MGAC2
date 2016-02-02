#include "gasp2.hpp"
#include "optionparser.h"

using namespace std;

const string version = "2.0-pre-beta";

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

enum optIndex {INPUT,HELP,RESTART,SPACEGROUPS,CONVERT,SIZE,TEMPLATE,PLANE,THREADS };

const option::Descriptor usage[] =
{
		{INPUT,0,"i","input",Arg::Required,"-i,--input  The XML input file for running"},
		{HELP,0,"h","help",option::Arg::Optional,"-h,--help  Displays this help message"},
		{RESTART,0,"r","restart",Arg::NonEmpty,"-r,--restart  Optional argument for a restart file"},
		//{STEP,0,"S","step",Arg::NonEmpty,"-S,--step Optional argument used to specify starting step"},
		{SPACEGROUPS,0,"l","spacegroups",Arg::Optional,"-l,--spacegroups  List valid spacegroups"},
		{CONVERT, 0, "c", "cif", Arg::NonEmpty, "-c, --cif  Name of output file for cif; takes the input (-i) and turns it into a cif"},
		//{COMBINE, 0, "m","merge",Arg::NonEmpty, "-m, --merge Combine multiple files to form a single population file (comma delimited)"},
		{SIZE, 0, "s","size", Arg::NonEmpty, "-s, --size Used in conjunction with cif conversion to determine pop size"},
		{TEMPLATE, 0, "t","template",Arg::NonEmpty, "-t, --template Creates an XML molecule template from a cif; if a plane is given then rotation and other values will be checked"},
		{PLANE, 0, "p", "plane",Arg::NonEmpty, "-p, --plane Specifies the three atom plane to be used for a template (comma delimited)"},
		//{RECLUSTER, 0, "g", "recluster",Arg::NonEmpty, "-g, --recluster Reclusters a population using a standard input for input and the population to be reclustered as the g argument"},
		{THREADS, 0, "j", "threads",Arg::NonEmpty, "-j, --threads Number of threads to use (ie, clustering)"},
		{ 0, 0, 0, 0, 0, 0 }
};


int main( int argc, char* argv[] ) {

	int mpithreading;
	//MPI_Init(&argc,&argv);
	MPI_Init_thread(&argc,&argv, MPI_THREAD_FUNNELED, &mpithreading);
//	if(mpithreading == MPI_THREAD_FUNNELED)
//		cout << "Threading requested was given!" << endl;

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


    if(options[CONVERT]) {
    	if(options[CONVERT].arg == NULL) {
    		return 0;
    	}
    	if(options[INPUT] && options[CONVERT].arg != NULL) {

    		size = 0;
    		if(options[SIZE]) {
    			size = std::stoi(string(options[SIZE].arg));
    		}

    		//TODO: need a way to set table name via arg

    		GASP2pop temppop;

    		GASP2db db;
    		db.load(options[INPUT].arg);
    		if(size == 0)
    			temppop = db.getAll("structs");
    		else
    			temppop = db.getBest(size,"structs");


    		remove(options[CONVERT].arg);
    		//doc.Clear();
    		cout << mark() << "pop loaded" << endl;
    		temppop.energysort();
    		cout << mark() << "sorted" << endl;
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


    time_t prog_start;
    //startup control for server/client
    if(ID==0) {

        if(!options[INPUT] && !options[RESTART] ) {
        	cout << "An input or restart file is required! Use -i or -r to specify." << endl;
        	option::printUsage(cout, usage);
        	MPI_Abort(MPI_COMM_WORLD, 1);
        	return 0;
        }


        string infile = "", restart = "";

        if(options[INPUT] && options[INPUT].arg != NULL)
        	infile.assign(options[INPUT].arg);
        if( options[RESTART] && options[RESTART].arg != NULL)
        	restart.assign(options[RESTART].arg);


    	prog_start = time(0);
    	char * dt = ctime(&prog_start);
    	cout << "------------------------------------------------------------" << endl;
    	cout << "MGAC v"<< version << " startup at " << dt;
    	cout << endl;
    	cout << "Users of this program should cite: " << endl;
    	cout << endl << endl;
    	cout << "------------------------------------------------------------" << endl << endl;
    	cout << "Using input file " << infile << endl << endl;


		GASP2control server(prog_start,size, infile, restart);
		server.server_prog();

    }
    else {
    	GASP2control client(ID);
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
