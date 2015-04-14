#include "gasp2qe.hpp"

using namespace std;

namespace QE {

	bool empty(vector<GASP2molecule>& mols, GASP2cell& unit, double& energy, double& force, double&pressure, time_t&time, string hostfile, GASP2param params)
	{
		uniform_real_distribution<> d(0,1);
		chrono::milliseconds thread_wait(100);
		cout << "Running the energy eval.." << endl;
		this_thread::sleep_for(thread_wait);
		energy = -255.0 + d(rgen);
		return true;
	};

	bool runQE(vector<GASP2molecule>& mols, GASP2cell& unit, double& energy, double& force, double&pressure, time_t&time, string hostfile, GASP2param params) {

		//write the inputfile

		//build the launcher line scripts

		//launch via popen2

		//scrape the output and parse/check

		//cleanup inputs

	}

}
