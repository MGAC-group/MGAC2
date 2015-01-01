#include "gasp2.hpp"
#include "optionparser.h"

using namespace std;

const string version = "0.1";

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

enum optIndex {INPUT,HELP,RESTART };

const option::Descriptor usage[] =
{
		{INPUT,0,"i","input",Arg::Required,"--input  The XML input file for running"},
		{HELP,0,"h","help",option::Arg::Optional,"--help  Displays this help message"},
		{RESTART,0,"r","restart",Arg::NonEmpty,"--restart  Optional argument for a restart file"},
		{ 0, 0, 0, 0, 0, 0 }
};


int main( int argc, char* argv[] ) {

	MPI_Init(&argc,&argv);

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

    if(!options[INPUT]) {
    	cout << "An input file is required! Use -i to specify." << endl;
    	option::printUsage(cout, usage);
    	return 0;
    }
    string infile;
    infile.assign(options[INPUT].arg);


    time_t prog_start;
    //startup control
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

    		GASP2control server(prog_start,size-1, infile, options[RESTART].arg);
    		server.server_prog();
    	}
    	else {
    		GASP2control server(prog_start,size-1, infile);
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


	return 0;
}
