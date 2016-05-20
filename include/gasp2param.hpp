#pragma once

/*
 * The intention of gasp2param is to provide a place for
 * parameters to be collected in a sensible way. GASP2param
 * is declared as public to avoid having to write
 * a new handler for every new parameter, and because GASP2param
 * is more correctly placed at the same peer level as GASP2control.
 *
 * GASP2param should only really be used for parameters which are
 * not nicely transfered over MPI or things which are only interpreted
 * once during program setup, although it can contain other things fine.
 *
 * To add a parameter it should be added to the header, a default
 * added to the constructor, a mechanism for adding from
 * both XML and HDF5.
 *
 * The convention is a non-prefixed name is a general parameter,
 * a QE parameter prefixed with QE, a custom with C, charmm with
 * CHARMM, etc etc.
 */
#include "gasp2common.hpp"
using namespace std;



class GASP2param {

public:
	//front matter
	GASP2param();
	~GASP2param() {};

	//crystal params
	//double distcell; //0.05 ang
	//double intradist; //0.5 ang
	string name; //REQUIRED name of structure
	//double expectvol; //REQUIRED
	//double maxvol; //100.0
	//double minvol; //-50.0
	string spacegroup; //REQUIRED ?

	//run params
	string calcmethod; //REQUIRED: fitcell, qe, custom
	string mode; //steadystate or stepwise or final
	string type; //elitism or fullcross or finaleval
	string selector; //roulette or pattern
	string scaling; //constant, linear, exponential
	double cross_prob; //1.0
	double replacement; //percentage, 1.0
	double mutation_prob; //0.01
	int generations; //50
	int popsize; //30
	int seed; //0
	bool test_seed; //if true; sets seed to 0
	double const_scale; //scaling parameters
	double lin_scale; // ..
	double exp_scale;// ..
	string outputmode; //SQLITE ONLY
	string outputfile; //name of outputfile
	Spacemode spacemode; //single, limited, or full spacegroup mode
	Index group; //single spacegroup specification
	int symmlimit; //limit to number of symm groups
	int binlimit; //number of bins to limit search on
	int downlimit; //number of downed nodes
	int precompute;
	double clusterdiff; //percent difference between two structs
	double clustersize;
	double chebyshevlimit;
	//the following parameters are for clustering
	double dihstep; //minimum difference between dihedrals
	double rotstep; // rotations difference
	double ratiostep; //ratio difference
	double posstep; //fractional position difference
	double angstep; //cell angle differences



	//qe params
	string QEpath; //REQUIRED
	string QEmpirunpath; //REQUIRED
	string QEpreamble; //in case something else needs to be executed before qe is invoked
	string QEcalculation; //vc-relax
	string QErestart_mode; //from_scratch
	string QEprefix; //REQUIRED
	string QEtstress; //.true.;
	string QEtprnfor; //.true.;
	int QEnstep; //20
	string QEpseudo_dir; //REQUIRED
	string QEoutdir; //REQUIRED
	string QEwf_collect; //.true.
	string QEverbosity; //high?
	string QEetot_conv_thr; //1.0D-3
	string QEforc_conv_thr; //1.0D-2
	string QEpress_conv_thr; //0.5D0
	string QEecutwfc; //55
	string QEecutrho; //550
	string QEspline_ps; //.true.
	string QEvdw_corr; //DFT-D
	string QEconv_thr; //1.D-7
	string QEcell_dynamics; //bfgs
	string QEk_points; //automatic
	string QEk_point_spec; //2 2 2   1 1 1
	int QErestart_limit; //2
	int QEscftimeout; //time in seconds
	vector<string> QEpseudos; //REQUIRED
//	int QEscalemin;
// 	int QEscalemax;

	//CHARMM params
	//AML: at some point I may add in CHARMM support

	//custom params
	string Cscratch; //REQUIRED location of scratch handling
	string Cscript; //REQUIRED name of conversion/running script

	//parsers and writers
	bool parseXML( tinyxml2::XMLDocument *doc, string& errorstring );
	void logParams();


};
