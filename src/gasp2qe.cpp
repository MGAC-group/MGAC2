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

		//setup the temp files and scratch paths
		UUID runID;
		runID.generate();
		string name = runID.toStr();
		string prefix = params.QEoutdir + "/" + name;
		string infile = prefix + ".in";
		//string outfile = prefix + ".out";


		//write the inputfile
		stringstream input;
		input << std::fixed << std::setprecision(6);
		input << "&control" << endl;
		input << "  calculation = " << params.QEcalculation << endl;
		input << "  restart_mode = " << params.QErestart_mode << endl;
		input << "  prefix = " << (params.QEprefix + "-" + name) << endl;
		input << "  tstress = " << params.QEtstress << endl;
		input << "  tprnfor = " << params.QEtprnfor << endl;
		input << "  nstep = " << params.QEnstep << endl;
		input << "  pseudo_dir = '" << params.QEpseudo_dir << "'"<< endl;
		input << "  outdir = '" << params.QEoutdir << "'" << endl;
		input << "  wf_collect = " << params.QEwf_collect << endl;
		input << "  verbosity = " << params.QEverbosity << endl;
		input << "  etot_conv_thr = " << params.QEetot_conv_thr << endl;
		input << "  forc_conv_thr = " << params.QEforc_conv_thr << endl;
		input << "/" << endl;
		input << "&system" << endl;
		input << "  ibrav = 14" << endl;
		input << "    A = " << unit.a << endl;
		input << "    B = " << unit.b << endl;
		input << "    C = " << unit.c << endl;
		input << "    cosBC = " << cos(unit.alpha) << endl;
		input << "    cosAC = " << cos(unit.beta) << endl;
		input << "    cosAB = " << cos(unit.gamma) << endl;

		int nat = 0;
		for(int i = 0; i < mols.size(); i++) {
			nat += mols[i].atoms.size();
		}
		int ntyp = params.QEpseudos.size();

		input << "  nat = " << nat << endl;
		input << "  ntyp = " << ntyp << endl;
		input << "  ecutwfc = " << params.QEecutwfc << endl;
		input << "  ecutrho = " << params.QEecutrho << endl;
		input << "  spline_ps = " << params.QEspline_ps << endl;
		//have to make a decision about this...what version to support?
		//this only appears to work in 5.1, not 5.0.2 or less
		//input << "  vdw_corr = " << params.QEvdw_corr << endl;
		input << "  london = .true." << endl;
		input << "/" << endl;
		input << "&electrons" << endl;
		input << "  conv_thr = " << params.QEconv_thr << endl;
		input << "/" << endl;
		input << "&ions" << endl;
		input << "/" << endl;
		input << "&cell" << endl;
		input << "  cell_dynamics = " << params.QEcell_dynamics << endl;
		input << "  press_conv_thr = " << params.QEpress_conv_thr << endl;
		input << "/" << endl << endl;

		input << "ATOMIC_SPECIES" << endl;
		for(int i = 0; i < ntyp; i++) {
			input << params.QEpseudos[i] << endl;
		}
		input << endl;

		Vec3 temp;
		Mat3 toFrac = cartToFrac(unit);
		input << "ATOMIC_POSITIONS (crystal)" << endl;
		for (int i = 0; i < mols.size(); i++) {
			for (int j = 0; j < mols[i].atoms.size(); j++) {
				input << getElemName(mols[i].atoms[j].type) << " ";
				temp = toFrac*mols[i].atoms[j].pos + mols[i].pos;
				input << temp[0] << " ";
				input << temp[1] << " ";
				input << temp[2] << " (1 1 1)" << endl;
			}
		}
		input << endl;

		input << "K_POINTS " << params.QEk_points << endl;
		input << params.QEk_point_spec << endl;


		ofstream outf;
		outf.open(infile.c_str(), ofstream::out);
		if(outf.fail()) {
			cout << "A problem was detected when trying to write a QE input file!" << endl;
			return false;
		}

		outf << input.str();
		outf.close();




		//build the launcher line scripts
		stringstream launcher;

		launcher << params.QEpreamble <<";";
		launcher << params.QEmpirunpath << " ";

		launcher << params.QEpath << " -inp " << infile;

		//launch via popen2
		FILE * out;
		pid_t t;
		out = popen2(launcher.str().c_str(), t);

		//search strings for solving issues
        string energy_scan = "!    total energy              =";
        string force_scan = "Total force =";
        string press_scan = "total   stress  (Ry/bohr**3)                   (kbar)     P=";
        string any_error = "Error in routine";
        string fft_error = "Not enough space allocated for radial FFT: try restarting with a larger cell_factor";
        string complete_scan = "This run was terminated on";
        string early_term = "The maximum number of steps has been reached.";



		//scrape the output and parse/check
		//starting point
		char buff[1024];
		string output = "", junk;
		output.reserve(2000000);
		int length, oldlength = 0;
		size_t pos;
		stringstream ss;
		int status;

		auto start = chrono::steady_clock::now();

	    while( chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start).count()
	    		< 600 ) {//params.QEscftimeout) {
	        while(fgets(buff, sizeof(buff), out)) {
	          //cout << buff;
	        	output.append(buff);
	        }
	        length = output.size();
	        if(length > oldlength) {
	        	pos = output.rfind(energy_scan);
	        	if(pos!=string::npos) {
	        		ss.clear();
	        		 ss.str(output.substr(pos, 160));

	        		 cout << ss.str() << endl;
	        		 ss >> junk >> junk >> junk >> junk;
	        		 cout << "junk: " << junk << endl;
	        		 ss >> energy;
	        		 cout << "New energy: " << energy << endl;
	        	}




	        	oldlength = length;
	        }

	        pid_t result = waitpid(t, &status, WNOHANG);
	        if(result == 0) cout << "alive" << endl;
	        else if(result == -1) cout << "problems" << endl;
	        else { cout << "done!" << endl; break; }

	        //cout << buff << endl;

	        //sbuff.append(buff);
	        //if(sbuff.find("Error")!=string::npos) {
	        //  flag = true; pclose(in);
	        //}
	        chrono::milliseconds t(200);
	        this_thread::sleep_for(t);
	    }

	    pclose2(t,out);
		//cleanup inputs

	}

}
