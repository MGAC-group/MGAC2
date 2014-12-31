#include "gasp2param.hpp"

using namespace std;

GASP2param::GASP2param() {

	//this is where the bulk of defaults are set

	//crystal params
	double distcell = 0.05;// ang
	double intradist = 0.5;// ang
	double maxvol = 100.0; //100.0
	double minvol = -50.0; //-50.0

	//run params
	string calcmethod = "fitcell"; //fitcell, qe, custom
	string mode = "generational"; //steadystate or generational
	string type = "elitism"; //elitism or classic
	string selector = "roulette"; //roulette or pattern
	string scaling = "linear"; //constant, linear, exponential
	double cross_prob = 1.0; //1.0
	double replacement = 1.0; //percentage, 1.0
	double mutation_prob = 0.01; //0.01
	int generations = 50; //50
	int popsize = 30; //30
	int seed = 0;
	bool test_seed = false; //if true; sets seed to 0
	int clientseed = 0; //internal
	double const_scale = 1.0;
	double lin_scale = 1.0;
	double exp_scale = 0.5;


	//qe params
	string QEcalculation = "vc-relax"; //vc-relax
	string QErestart_mode = "from_scratch"; //from_scratch
	string QEtstress = ".true."; //.true.;
	string QEtprnfor = ".true."; //.true.;
	int QEnstep = 70; //70
	string QEwf_collect = ".true."; //.true.
	string QEverbosity = ".true."; //high?
	string QEetot_conv_thr = "1.0D-3"; //1.0D-3
	string QEforc_conv_thr = "1.0D-2"; //1.0D-2
	string QEpress_conv_thr = "0.5D0"; //0.5D0
	string QEecutwfc = "55"; //55
	string QEecutrho = "550"; //550
	string QEspline_ps = ".true."; //.true.
	string QElondon = ".true."; //.true.
	string QEconv_thr = "1.D-7"; //1.D-7
	string QEcell_dynamics = "bfgs"; //bfgs
	string QEk_points = "automatic"; //automatic
	string QEk_point_spec = "2 2 2  1 1 1"; //2 2 2   1 1 1
	int QErestart_limit = 3; //3
	int QEscftimeout = 6000; //time in seconds

}

bool GASP2param::parseXML(tinyxml2::XMLDocument *doc) {

	tinyxml2::XMLElement *crystal = doc->FirstChildElement("mgac")->FirstChildElement("crystal");

	double temp;

	if(!crystal->QueryDoubleAttribute("expectvol", &temp))
		cout << "Expectvol is " << temp << endl;
	else
		cout << "There was an error parsing!" << endl;


	return true;

}

void GASP2param::logParams() {


}
