#include "gasp2qe.hpp"

using namespace std;

namespace QE {

	bool empty(vector<GASP2molecule>& mols, GASP2cell& unit, double& energy, double& force, double&pressure, time_t&time, string hostfile, GASP2param params)
	{
		uniform_real_distribution<> d(0,1);
		//cout << "Running the energy eval.." << endl;
		this_thread::sleep_for(chrono::seconds(5));
		energy = -255.0 + d(rgen);
		return true;
	};

	bool runQE(vector<GASP2molecule>& mols, GASP2cell& unit, double &energy, double &force, double&pressure, time_t&time, string hostfile, GASP2param params) {

		bool outstat = true;

		//setup the temp files and scratch paths
		UUID runID;
		runID.generate();
		string name = runID.toStr();
		string prefix = params.QEoutdir + "/" + name + "/";
		string infile = prefix + name + ".in";
		//string outfile = prefix + ".out";
		int result = mkdir(prefix.c_str(), 0700);
		if(result == -1) {
			cout << "Cannot make QE temp directory! Aborting..." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
			exit(1);
		}


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
		input << "  outdir = '" << prefix << "'" << endl;
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
		launcher << "-machinefile " << hostfile << " ";
		launcher << params.QEpath << " -inp " << infile;

		//launch via popen2
		FILE * out;
		pid_t t;
		out = popen2(launcher.str().c_str(), t);

		//search strings for QE output
        const string energy_scan = "!    total energy              =";
        const string force_scan = "Total force =";
        const string press_scan = "total   stress  (Ry/bohr**3)                   (kbar)     P=";
        const string any_error = "Error in routine";
        const string fft_error = "Not enough space allocated for radial FFT: try restarting with a larger cell_factor";
        const string complete_scan = "This run was terminated on";
        const string early_term = "The maximum number of steps has been reached.";
        const string cell_param = "CELL_PARAMETERS";
        const string atom_pos = "ATOMIC_POSITIONS";



		//scrape the output and parse/check
		//starting point
		char buff[1024];
		string output = "", junk;
		output.reserve(2000000);
		int length, oldlength = 0;
		size_t pos, epos = 0, fpos = 0, ppos = 0, cp_pos = 0, ap_pos = 0;
		stringstream ss;
		int status;
		Mat3 toCart;
		Vec3 va,vb,vc;
		double alat;
		int numprocesses;

		save_state=false;

		auto start = chrono::steady_clock::now();

// debug QE output
//		outf.open("qedebug.out", ofstream::out);
//		if(outf.fail()) {
//			cout << "A problem was detected when trying to write a QE debug file!" << endl;
//			return false;
//		}

        chrono::seconds start_wait(15);
        this_thread::sleep_for(start_wait);

	    while(true) {
	    	if(chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start).count()
		    		>= params.QEscftimeout) {
	    		cout << mark() << "Timeout limit reached for a QE run!\n";
	    		outstat = false;
	    		break;
	    	}

	        while(fgets(buff, sizeof(buff), out)) {
	        	//outf << buff;
	        	output.append(buff);
	        }
	        length = output.size();
	        if(length > oldlength) {

	        	while(!longeval_mut.try_lock()) {
	        		this_thread::sleep_for(chrono::milliseconds(200));
	        	}

	        	//energy - four junk fields
	        	pos = output.rfind(energy_scan);
	        	if(pos!=string::npos && pos > epos) {
	        		ss.clear();
	        		ss.str(output.substr(pos, 160));
	        		ss >> junk >> junk >> junk >> junk;
	        		ss >> energy;
	        		//the energy used here does not take into account basic stoich
	        		//only the number of spacegroup operations
	        		energy /= static_cast<double>(spacegroups[unit.spacegroup].R.size());
	        		energy *= RydToKCalMol;
	        		//cout << "New energy: " << energy << endl;
	        		epos = pos;
	        	}

	        	//force - three junk fields
	        	pos = output.rfind(force_scan);
	        	if(pos!=string::npos && pos > fpos) {
	        		ss.clear();
	        		ss.str(output.substr(pos, 160));
	        		ss >> junk >> junk >> junk;
	        		ss >> force;
	        		//cout << "New force: " << force << endl;
	        		fpos = pos;
	        	}

	        	//press - five junk fields
	        	pos = output.rfind(press_scan);
	        	if(pos!=string::npos && pos > ppos) {
	        		ss.clear();
	        		ss.str(output.substr(pos, 160));
	        		ss >> junk >> junk >> junk >> junk >> junk;
	        		ss >> pressure;
	        		//cout << "New pressure: " << pressure << endl;
	        		ppos = pos;
	        	}

	        	//get the cell params
	        	pos = output.rfind(cell_param);
	        	if(pos!=string::npos && pos > cp_pos) {
	        		ss.clear();
	        		ss.str(output.substr(pos,80*4));
	        		//get alat
	        		ss >> junk >> junk >> alat >> junk;

	        		//get vectors
	        		ss >> va[0] >> va[1] >> va[2];
	        		ss >> vb[0] >> vb[1] >> vb[2];
	        		ss >> vc[0] >> vc[1] >> vc[2];

//	        		cout << alat << endl << va << vb << vc << endl;;
	        		unit.a = alat * AUtoAng * len(va);
	        		unit.b = alat * AUtoAng * len(vb);
	        		unit.c = alat * AUtoAng * len(vc);
	        		unit.gamma = angle(va, vl_0, vb);
	        		unit.beta = angle(va, vl_0, vc);
	        		unit.alpha = angle(vb, vl_0, vc);
	        		//cout << unit.a << " " << unit.b << " " << unit.c << endl;
	        		//cout << unit.alpha <<" "<<unit.beta <<" "<<unit.gamma << endl;

	        		toCart = fracToCart(unit);
	        		cp_pos = pos;
	        	}

	        	//get the coordinates
	        	pos = output.rfind(atom_pos);
	        	if(pos!=string::npos && pos > ap_pos) {
	        		ss.clear();
	        		ss.str(output.substr(pos,80+80*nat));
	        		//cout << ss.str() << endl;
	        		//cout
	        		//clear the first line
	        		ss >> junk >> junk;
	        		for (int i = 0; i < mols.size(); i++) {
	        			for (int j = 0; j < mols[i].atoms.size(); j++) {
	        				ss >> junk;
	        				ss >> temp[0] >> temp[1] >> temp[2];
	        				//cout << temp << endl;
	        				mols[i].atoms[j].pos = toCart*temp;
	        				//temp = toFrac*mols[i].atoms[j].pos + mols[i].pos;
	        			}

	        		}
	        		ap_pos = pos;
	        		save_state=true;
	        	}




	        	//fft error - not an issue, just needs
	        	//a restart if appropriate
	        	pos = output.rfind(fft_error);
	        	if(pos!=string::npos) {
	        		cout << mark() << "FFT error detected in QE run!\n";
	        		outstat = false;
	        		longeval_mut.unlock();
	        		break;
	        	}

	        	//unknown error - this is bad, we want to know
	        	//when it happens and we want to exit as fast as possible
	        	//this could leave stray processes, but it's not clear
	        	pos = output.rfind(any_error);
	        	if(pos!=string::npos) {
	        		cout << mark() << "Unknown error detected in QE run! Output:\n";
	        		cout << output.substr(pos, 2048) << endl << endl;
	        		pclose2(t,out);
	        		//send emergency shutdown signal instead of exit?
	        		exit(1);
	        	}

	        	//early term
	        	pos = output.rfind(early_term);
	        	if(pos!=string::npos) {
	        		cout<< mark() << "QE finished with maximum number of SCF cycles.\n";
	        		outstat = false;
	        		longeval_mut.unlock();
	        		break;
	        	}

	        	//normal finish
	        	pos = output.rfind(complete_scan);
	        	if(pos!=string::npos) {
	        		cout << mark()  << "QE completed the run\n";
	        		longeval_mut.unlock();
	        		break;
	        	}

	        	oldlength = length;
	        	longeval_mut.unlock();
	        }



	        pid_t result = waitpid(t, &status, WNOHANG);
	        if(result == 0) ;//cout << "alive" << endl;
	        else if(result == -1) cout << "problems" << endl;
	        else { cout << "done!" << endl; break; }

	        //cout << buff << endl;

	        //sbuff.append(buff);
	        //if(sbuff.find("Error")!=string::npos) {
	        //  flag = true; pclose(in);
	        //}
	        chrono::milliseconds t(500);
	        this_thread::sleep_for(t);
	    }

		//outf.close();


	    pclose2(t,out);
		//cleanup inputs
	    string rem = "rm -rf " + prefix;
	    system(rem.c_str());

        //add the time:
	    int duration = chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start).count();
	    cout << mark() << "QE completed in: " << duration << endl;
        time += duration;


        completed = true;
	    return outstat;
	}

}
