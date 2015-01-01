#include "gasp2common.hpp"
#include "gasp2param.hpp"
#include "gasp2pop.hpp"
#include "gasp2struct.hpp"
#include "gasp2qe.hpp"


using namespace std;

class GASP2control {
public:
	GASP2control(int ID, string infile);
	GASP2control(time_t start, int size, string input, string restart="");
	void server_prog();
	void client_prog();



private:
	//procedural variables
	time_t starttime;
	string infile;
	string restart;
	int worldSize;
	int ID;

	GASP2param params;

	bool parseInput(tinyxml2::XMLDocument *doc, string & errors);

};
