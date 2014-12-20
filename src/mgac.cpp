#include "gasp2.hpp"
#include "optionparser.h"
#include <mpi.h>

using namespace std;

const string version = "0.1";


enum optIndex {RESTART, INPUT, UNKNOWN, HELP };

const option::Descriptor usage[] =
{
		{INPUT,0,"i","input",option::Arg::Optional,"--input  The XML input file for running"},


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

    string infile;
    if(options[INPUT] && options[INPUT].arg != NULL) {
    	infile(options[INPUT].arg);
    }
    else {

    	exit(1);
    }

    //startup control
    if(ID==0) {
    	time_t prog_start = time(0);
    	char * dt = ctime(&prog_start);
    	cout << "------------------------------------------------------------" << endl;
    	cout << "MGAC v"<< version << " startup at " << dt;
    	cout << endl << endl;
    	cout << "Users of this program should cite: " << endl;
    	cout << endl << endl;
    	cout << "------------------------------------------------------------" << endl;

    	if( options[RESTART] && )
    		GASP2control server(prog_start,size-1, infile, )
    	else
    		GASP2control server(prog_start, infile);

    	server.server_prog();

    }
    else {

    	GASP2control client();
    	client.client_prog();

    }

    //program termination
    if(ID==0) {
    	time_t final = time(0) - prog_start;
    	char *dt = ctime(&final);
    	cout << "------------------------------------------------------------" << endl;
    	cout << "MGAC completed at " << final;
    	cout << "------------------------------------------------------------" << endl;
    }


	return 0;
}
