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

enum optIndex {INPUT,HELP,RESTART,SPACEGROUPS,CONVERT };

const option::Descriptor usage[] =
{
		{INPUT,0,"i","input",Arg::Required,"-i,--input  The XML input file for running"},
		{HELP,0,"h","help",option::Arg::Optional,"-h,--help  Displays this help message"},
		{RESTART,0,"r","restart",Arg::NonEmpty,"-r,--restart  Optional argument for a restart file"},
		{SPACEGROUPS,0,"ls","spacegroups",Arg::Optional,"-ls,--spacegroups  List valid spacegroups"},
		{CONVERT, 0, "c", "cif", Arg::NonEmpty, "-c, --cif  Name of output file for cif; takes the input (-i) and turns it into a cif"},
		{ 0, 0, 0, 0, 0, 0 }
};


int main( int argc, char* argv[] ) {

	int mpithreading;
	MPI_Init(&argc,&argv);
	//MPI_Init_thread(&argc,&argv, MPI_THREAD_MULTIPLE, &mpithreading);
	if(mpithreading == MPI_THREAD_MULTIPLE)
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

    		GASP2pop temppop;
    		string errorstring;
    		tinyxml2::XMLDocument doc;
    		doc.LoadFile(options[INPUT].arg);
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
    		temppop.writeCIF(string(options[CONVERT].arg));
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
    		GASP2control server(prog_start,size, infile, options[RESTART].arg);
    		server.server_prog();
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
