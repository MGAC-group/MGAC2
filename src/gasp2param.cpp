#include "gasp2param.hpp"

using namespace std;

GASP2param::GASP2param() {

	//this is where the bulk of defaults are set

	//crystal params
	//distcell = 0.05;// ang
	//intradist = 0.5;// ang
	//maxvol = 100.0; //100.0
	//minvol = -50.0; //-50.0

	//run params
	calcmethod = "fitcell"; //fitcell, qe, custom
	mode = "generational"; //steadystate or generational
	type = "elitism"; //elitism or classic
	selector = "roulette"; //roulette or pattern
	scaling = "linear"; //constant, linear, exponential
	cross_prob = 1.0; //1.0
	replacement = 1.0; //percentage, 1.0
	mutation_prob = 0.01; //0.01
	generations = 50; //50
	popsize = 30; //30
	seed = 0;
	test_seed = false; //if true; sets seed to 0
	const_scale = 1.0;
	lin_scale = 1.0;
	exp_scale = 0.5;
	outputmode = "xml";
	outputfile = "population_data";

	//qe params
	QEcalculation = "vc-relax"; //vc-relax
	QErestart_mode = "from_scratch"; //from_scratch
	QEtstress = ".true."; //.true.;
	QEtprnfor = ".true."; //.true.;
	QEnstep = 20; //70
	QEwf_collect = ".true."; //.true.
	QEverbosity = ".true."; //high?
	QEetot_conv_thr = "1.0D-3"; //1.0D-3
	QEforc_conv_thr = "1.0D-2"; //1.0D-2
	QEpress_conv_thr = "0.5D0"; //0.5D0
	QEecutwfc = "55"; //55
	QEecutrho = "550"; //550
	QEspline_ps = ".true."; //.true.
	QEvdw_corr = "DFT-D"; //.true.
	QEconv_thr = "1.D-7"; //1.D-7
	QEcell_dynamics = "bfgs"; //bfgs
	QEk_points = "automatic"; //automatic
	QEk_point_spec = "2 2 2  1 1 1"; //2 2 2   1 1 1
	QErestart_limit = 3; //3
	QEscftimeout = 6000; //time in seconds\

	QEpreamble = "";

}


bool GASP2param::parseXML(tinyxml2::XMLDocument *doc, string& errorstring) {

	double dtemp;
	const char * stemp;
	int itemp;

	//setdefs();

//*************************
//
// 		CRYSTAL TAG
//
//*************************
	tinyxml2::XMLElement *crystal = doc->FirstChildElement("mgac")->FirstChildElement("crystal");
	if(crystal) {
//AML: Wrong place
//		//		//expectvol
//		if(!crystal->QueryDoubleAttribute("expectvol", &dtemp)) {
//			expectvol = dtemp;
//			if(expectvol <= 0.0) {
//				errorstring = "expectvol must be greater than 0.0!\n"; return false;
//			}
//		}
//		else {
//			errorstring = "The expected volume (expectvol) was not specified but is required!\n"; return false;
//		}
//
//		//distcell
//		if(!crystal->QueryDoubleAttribute("distcell", &dtemp)) {
//			distcell = dtemp;
//			if(distcell <= 0.0) {
//				errorstring = "distcell must be greater than 0.0!\n";
//				return false;
//			}
//		}
//
//		//intradist
//		if(!crystal->QueryDoubleAttribute("intradist", &dtemp)) {
//			intradist = dtemp;
//			if(intradist <= 0.0) {
//				errorstring = "intradist must be greater than 0.0!\n";
//				return false;
//			}
//		}
//
//		//maxvol
//		if(!crystal->QueryDoubleAttribute("maxvol", &dtemp))
//			maxvol = dtemp;
//
//		//minvol
//		if(!crystal->QueryDoubleAttribute("minvol", &dtemp))
//			minvol = dtemp;
//
//		if(minvol > maxvol) {
//			errorstring = "minvol cannot be greater than maxvol! (defaults are max=+100.0, min=-50.0)\n";
//			return false;
//		}

		//name
		stemp = crystal->Attribute("name");
		if(stemp)
			name.append(stemp);
		else {
			errorstring = "A name was not specified but is required!\n"; return false;
		}

//AML: wrong place
//		//spacegroup
//		stemp = crystal->Attribute("spacegroup");
//		if(stemp)
//			spacegroup.append(stemp);
//		else {
//			errorstring = "A spacegroup was not specified but is required!\n"; return false;
//		}

	}
	else {
		errorstring = "The 'crystal' tag is required but was not found!\n";
		return false;
	}

//*************************
//
// 		RUN TAG
//
//*************************
	tinyxml2::XMLElement *run = doc->FirstChildElement("mgac")->FirstChildElement("run");
	if(run) {
		//calcmethod
		stemp = run->Attribute("calcmethod");
		if(stemp) {
			calcmethod = "";
			calcmethod.append(stemp);
			if (calcmethod=="fitcell" ||
				calcmethod=="qe" ||
				calcmethod=="custom") { ; }
			else {
				errorstring = "The calcmethod specified does not match a valid method!\n";
				return false;
			}
		}
		else {
			errorstring = "A calcmethod was not specified but is required!\n"; return false;
		}


		//mode
		stemp = run->Attribute("mode");
		if(stemp) {
			mode = "";
			mode.append(stemp);
			if (mode=="steadystate" || mode =="stepwise") { ; }
			else {
				errorstring = "The run mode specified does not match a valid type!\n";
				return false;
			}
		}

		//type
		stemp = run->Attribute("type");
		if(stemp) {
			type = "";
			type.append(stemp);
			if (type=="elitism" || type=="classic") { ; }
			else {
				errorstring = "The run type specified does not match a valid type!\n";
				return false;
			}
		}

		//selector
		stemp = run->Attribute("selector");
		if(stemp) {
			selector = "";
			selector.append(stemp);

			if (selector=="roulette") { ; }
			else if (selector=="pattern") {
				errorstring = "The pattern selector is not implemented yet!\n";
				return false;
			}
			else {
				errorstring = "The selector specified does not match a valid type!\n";
				return false;
			}
		}

		//scaling
		stemp = run->Attribute("scaling");
		if(stemp) {
			scaling = "";
			scaling.append(stemp);
			if (scaling=="linear" ||
				scaling=="constant" ||
				scaling=="exponential") { ; }
			else {
				errorstring = "The scaling mode specified does not match a valid type!\n";
				return false;
			}
		}

		//outputfile
		stemp = run->Attribute("outputfile");
		if(stemp) {
			outputfile = "";
			outputfile.append(stemp);
		}

		//outputmode
		stemp = run->Attribute("outputmode");
		if(stemp) {
			outputmode = "";
			outputmode.append(stemp);
			if (outputmode=="xml" ||
				outputmode=="hdf5" ||
				outputmode=="XML" ||
				outputmode=="HDF5" ||
				outputmode=="h5" ) { ; }
			else {
				errorstring = "The population output mode specified does not match a valid type!\n";
				return false;
			}
		}


		//cross_prob
		if(!run->QueryDoubleAttribute("cross_prob", &dtemp)) {
			cross_prob = dtemp;
			if(cross_prob <= 0.0) {
				errorstring = "cross_prob must be greater than 0.0!\n";
				return false;
			}
		}

		//replacement
		if(!run->QueryDoubleAttribute("replacement", &dtemp)) {
			replacement = dtemp;
			if(replacement <= 0.0) {
				errorstring = "replacement must be greater than 0.0!\n";
				return false;
			}
		}

		//mutation
		if(!run->QueryDoubleAttribute("mutation", &dtemp)) {
			mutation_prob = dtemp;
			if(mutation_prob < 0.0) {
				errorstring = "mutation must be non-negative!\n";
				return false;
			}
		}

		//generations
		if(!run->QueryIntAttribute("generations", &itemp)) {
			generations = itemp;
			if(generations <= 0) {
				errorstring = "Generations must be greater than 0!\n";
				return false;
			}
		}

		//popsize
		if(!run->QueryIntAttribute("popsize", &itemp)) {
			popsize = itemp;
			if(popsize <= 0) {
				errorstring = "Population size must be greater than 0!\n";
				return false;
			}
		}

		//seed
		if(!run->QueryIntAttribute("seed", &itemp)) {
			seed = itemp;
			if(seed < -1) {
				errorstring = "The seed must be an unsigned integer value, or -1!\n";
				return false;
			}
		}

		//test_seed
		stemp = run->Attribute("test_seed");
		if(stemp)
			test_seed = (stemp == "true");
		if(test_seed)
			seed = 0;

		//const_scale
		if(!run->QueryDoubleAttribute("const_scale", &dtemp)) {
			const_scale = dtemp;
			if(const_scale < 0.0) {
				errorstring = "const_scale must be non-negative!\n";
				return false;
			}
		}

		//lin_scale
		if(!run->QueryDoubleAttribute("lin_scale", &dtemp)) {
			lin_scale = dtemp;
			if(lin_scale < 0.0) {
				errorstring = "lin_scale must be non-negative!\n";
				return false;
			}
		}

		//exp_scale
		if(!run->QueryDoubleAttribute("exp_scale", &dtemp)) {
			exp_scale = dtemp;
			if(exp_scale < 0.0) {
				errorstring = "exp_scale must be non-negative!\n";
				return false;
			}
		}

	}
	else {
		errorstring = "The 'run' tag is required but was not found!\n";
		return false;
	}


//*************************
//
// 		QE TAG
//
//*************************
	if(calcmethod == "qe") {

		tinyxml2::XMLElement * pseudo = doc->FirstChildElement("mgac")->FirstChildElement("qe")->FirstChildElement("pseudo");
		string ps, elem, mass, name;
		while(pseudo) {


			//elem
			stemp = pseudo->Attribute("elem");
			if(stemp)
				elem = stemp;
			else {
				errorstring = "An element must be specified for the pseudopotential!"; return false;
			}

			//mass
			stemp = pseudo->Attribute("mass");
			if(stemp)
				mass = stemp;
			else {
				errorstring = "A mass must be specified for the pseudopotential!"; return false;
			}

			//name
			stemp = pseudo->Attribute("name");
			if(stemp)
				name = stemp;
			else {
				errorstring = "A filename must be specified for the pseudopotential!"; return false;
			}

			ps = elem + " " + mass + " " + name;
			QEpseudos.push_back(ps);
			pseudo = pseudo->NextSiblingElement("pseudo");
		}


		tinyxml2::XMLElement *qe = doc->FirstChildElement("mgac")->FirstChildElement("qe")->FirstChildElement("param");
		if(qe) {

			//cout << "QE tags are NOT being validated yet!" << endl;

			//QEpath
			stemp = qe->Attribute("qepath");
			if(stemp)
				QEpath = stemp;
			else {
				errorstring = "The path to QE must be specified!"; return false;
			}


			//QEmpirunpath
			stemp = qe->Attribute("mpirunpath");
			if(stemp)
				QEmpirunpath = stemp;
			else {
				errorstring = "The path to mpirun must be specified!"; return false;
			}


			//QEprefix
			stemp = qe->Attribute("prefix");
			if(stemp)
				QEprefix = stemp;
			else {
				errorstring = "A prefix for the QE runs must be specified!"; return false;
			}


			//QEpseudo_dir
			stemp = qe->Attribute("pseudo_dir");
			if(stemp)
				QEpseudo_dir = stemp;
			else {
				errorstring = "The path to the QE pseudopotentials must be specified!"; return false;
			}


			//QEoutdir
			stemp = qe->Attribute("outdir");
			if(stemp)
				QEoutdir = stemp;
			else {
				errorstring = "A path to a scratch directory is required!"; return false;
			}


//			//
//			stemp = qe->Attribute("");
//			if(stemp)
//				 = stemp;
//			else {
//				errorstring = ""; return false;
//			}

			//QEpreamble
			stemp = qe->Attribute("preamble");
			if(stemp)
				QEpreamble = stemp;

			//QEcalculation
			stemp = qe->Attribute("calculation");
			if(stemp)
				QEcalculation = stemp;

			//QErestart_mode
			stemp = qe->Attribute("restart_mode");
			if(stemp)
				QErestart_mode = stemp;

			//QEtstress
			stemp = qe->Attribute("tstress");
			if(stemp)
				QEtstress = stemp;

			//QEtprnfor
			stemp = qe->Attribute("tprnfor");
			if(stemp)
				QEtprnfor = stemp;

			//QEnstep
			if(!qe->QueryIntAttribute("nstep", &itemp)) {
				QEnstep = itemp;
				if(QEnstep < 1) {
					errorstring = "The number of QE steps must be greater than zero!\n";
					return false;
				}
			}


			//QEwf_collect
			stemp = qe->Attribute("wf_collect");
			if(stemp)
				QEwf_collect = stemp;

			//QEverbosity
			stemp = qe->Attribute("verbosity");
			if(stemp)
				QEverbosity = stemp;

			//QEetot_conv_thr
			stemp = qe->Attribute("etot_conv_thr");
			if(stemp)
				QEetot_conv_thr = stemp;

			//QEforc_conv_thr
			stemp = qe->Attribute("forc_conv_thr");
			if(stemp)
				QEforc_conv_thr = stemp;

			//QEpress_conv_thr
			stemp = qe->Attribute("press_conv_thr");
			if(stemp)
				QEpress_conv_thr = stemp;

			//QEecutwfc
			stemp = qe->Attribute("ecutwfc");
			if(stemp)
				QEecutwfc = stemp;

			//QEecutrho
			stemp = qe->Attribute("ecutrho");
			if(stemp)
				QEecutrho = stemp;

			//QEspline_ps
			stemp = qe->Attribute("spline_ps");
			if(stemp)
				QEspline_ps = stemp;

			//QEvdw_corr
			stemp = qe->Attribute("vdw_corr");
			if(stemp)
				QEvdw_corr = stemp;

			//QEconv_thr
			stemp = qe->Attribute("conv_thr");
			if(stemp)
				QEconv_thr = stemp;

			//QEcell_dynamics
			stemp = qe->Attribute("cell_dynamics");
			if(stemp)
				QEcell_dynamics = stemp;

			//QEk_points
			stemp = qe->Attribute("k_points");
			if(stemp)
				QEk_points = stemp;

			//QEk_points_spec
			stemp = qe->Attribute("k_points_spec");
			if(stemp)
				QEk_point_spec = stemp;

			//QErestart_limit
			if(!qe->QueryIntAttribute("restart_limit", &itemp)) {
				QErestart_limit = itemp;
				if(QEnstep < 0) {
					errorstring = "The number of QE restarts must non-negative!\n";
					return false;
				}
			}

			//QEscftimeout
			if(!qe->QueryIntAttribute("scftimeout", &itemp)) {
				QEscftimeout = itemp;
				if(QEnstep < 20) {
					errorstring = "The scftimeout must be at least 20 seconds!\n";
					return false;
				}
			}

		}
		else {
			errorstring = "QE running mode was selected but no associated tag was found!\n";
		}
	}


	return true;

}

void GASP2param::logParams() {
	cout << "The following running parameters are being used for this prediction:" << endl << endl;

	if(test_seed)
		cout << "THIS RUN IS BEING PERFORMED AS A VALIDATION TEST USING SEED 0!" << endl;

	cout << "   Crystal name: " << name << endl;
	cout << "   Calculation method: " << calcmethod << endl;
	cout << endl;
	//cout << "   Crystal Parameters:" << endl;
	//cout << "      spacegroup   = " << setw(12) << spacegroup << endl;
	cout << setprecision(2) << fixed;
	//cout << "      distcell     = " << setw(12) << distcell << endl;
	//cout << "      intradist    = " << setw(12) << intradist << endl;
	//cout << "      expectvol    = " << setw(12) << expectvol << endl;
	//cout << "      minvol       = " << setw(12) << minvol << endl;
	//cout << "      maxvol       = " << setw(12) << maxvol << endl;
	cout << endl;
	cout << "   Running Parameters:" << endl;
	cout << setprecision(0) << fixed;
	cout << "      popsize      = " << setw(12) << popsize << endl;
	cout << "      generations  = " << setw(12) << generations << endl;
	cout << "      seed         = " << setw(12) << seed << endl;

	cout << setprecision(3) << fixed << endl;
	cout << "      replacement  = " << setw(12) << replacement << endl;
	cout << "      mutation     = " << setw(12) << mutation_prob << endl;
	cout << "      crossover    = " << setw(12) << cross_prob << endl;

	cout << "      mode         = " << setw(12) << mode << endl;
	cout << "      type         = " << setw(12) << type << endl;
	cout << "      selector     = " << setw(12) << selector << endl;
	cout << "      scaling      = " << setw(12) << scaling << endl;
	cout << setprecision(2) << fixed;
	if(scaling == "linear") {
		cout << "      lin_scale    = " << setw(12) << lin_scale << endl;
		cout << "      const_scale  = " << setw(12) << const_scale << endl;
	}
	else if (scaling == "exponential") {
		cout << "      exp_scale    = " << setw(12) << exp_scale << endl;
		cout << "      const_scale  = " << setw(12) << const_scale << endl;
	}

	if(calcmethod == "qe") {
		cout << endl;
		cout << "   QE Parameters:" << endl;
		cout << setprecision(0) << fixed;
		cout << "      calculation    =  " << setw(12) << QEcalculation << endl;
		cout << "      restart_mode   =  " << setw(12) << QErestart_mode << endl;
		cout << "      nstep          =  " << setw(12) << QEnstep << endl;
		cout << "      restart_limit  =  " << setw(12) << QErestart_limit << endl;
		cout << "      scftimeout     =  " << setw(12) << QEscftimeout << endl;
		cout << "      conv_thr       =  " << setw(12) << QEconv_thr << endl;
		cout << "      etot_conv_thr  =  " << setw(12) << QEetot_conv_thr << endl;
		cout << "      forc_conv_thr  =  " << setw(12) << QEforc_conv_thr << endl;
		cout << "      press_conv_thr =  " << setw(12) << QEpress_conv_thr << endl;
		cout << "      ecutwfc        =  " << setw(12) << QEecutwfc << endl;
		cout << "      ecutrho        =  " << setw(12) << QEecutrho << endl;
		cout << "      k_points       =  " << setw(12) << QEk_points << endl;
		cout << "      k_point_spec   =  " << setw(12) << QEk_point_spec << endl;
		cout << "      prefix         =  " << setw(12) << QEprefix << endl;
		cout << "      tstress        =  " << setw(12) << QEtstress << endl;
		cout << "      tprnfor        =  " << setw(12) << QEtprnfor << endl;
		cout << "      wf_collect     =  " << setw(12) << QEwf_collect << endl;
		cout << "      verbosity      =  " << setw(12) << QEverbosity << endl;
		cout << "      spline_ps      =  " << setw(12) << QEspline_ps << endl;
		cout << "      vdw_corr       =  " << setw(12) << QEvdw_corr << endl;
		cout << "      cell_dynamics  =  " << setw(12) << QEcell_dynamics << endl;


		cout << endl;
		cout << "   QE Paths: " << endl;
		cout << "      qepath:     " << QEpath << endl;
		cout << "      pseudo_dir: " << QEpseudo_dir << endl;
		cout << "      outdir:     " << QEoutdir << endl;
		cout << "      mpirunpath: " << QEmpirunpath << endl;

		cout << endl;
		cout << "   QE Preamble: " << endl << QEpreamble << endl;


		cout << endl;
		cout << "   QE Pseudopotentials:" << endl;
		for(int i = 0; i < QEpseudos.size(); i++) {
			cout << "      " << QEpseudos[i] << endl;
		}
		cout << endl;
	}






}
