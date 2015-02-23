#include "gasp2common.hpp"

using namespace std;



typedef enum StructError {
	OKStruct=0,
	OptBadBond,
	OptBadAng,
	OptBadDih,
	FitcellBadDih,
	NoFitcell,


}StructError;

struct GASP2atom {
	Vec3 pos; //ALWAYS fractional
	Elem type;
	Index label; //index to a string in "names";
};

struct GASP2bond {
	Index a,b; //indices to atoms
	double minLen, maxLen; //ANGSTROMS
	double len; //ANGSTROMS
};

struct GASP2angle {
	Index a,b,c;
	double maxAng,minAng; //RADIANS
	double ang; //RADIANS
};

struct GASP2dihedral {
	Index a,b,c,d; //indices to atoms;
	vector<Index> update; //list of atoms to update
	double maxAng, minAng; //RADIANS
	double ang; //RADIANS - XMLOUT
};

struct GASP2molecule {
	Vec3 pos; //ALWAYS fractional, refers to centroid XMLOUT
	Mat3 rot; //the rotation matrix of the molecule XMLOUT
	//Why rotation matrix? Space is cheap (realtively speaking)
	//
	vector<GASP2atom> atoms;

	Index label; //index to a string in "names";


	Index symm; //index of symmetry operation, 0 is always x,y,z
	Mat3 symmR; //given symmop rotation
	Vec3 symmT; //given symmop translation

	//bonds and angles are constraints;
	//at this time they only play a role in post-opt checking
	//dihedrals are modified during the fitcell operation
	vector<GASP2bond> bonds;
	vector<GASP2angle> angles;
	vector<GASP2dihedral> dihedrals;

	//Placement: when there is more than one molecule in the
	//asymmetric unit cell, then placement matters.
	//stoichiometric ratios must be preserved, otherwise energy comparisons
	//are somewhat useless. the presumption is that this is known, or can be
	//be guessed somewhow (ie, not our problem yet).
	//if plindex = 0, then pos is a fractional coordinate in the cell
	//if plindex > 0, then pos is a unit vector, with the origin at the centroid
	//of the previous index molecule. so, mol1 is 0,0,0->p1,p2,p3;
	//mol2 is p1,p2,p3->p1+pos1,p2+pos2,p3+po3, etc.
	//this results in much better clustering
	//using fitcell without QE optimization will also lead a set of structures
	//a poster I read suggests that the relative energies of cocrystal with different
	//ratios may be elucidated easily, provided the energy calculations are accurate
	//enough, which should be the case with MGAC.
	//the operators for multiple molecules shoudl be easy enough
	//by simply randomly picking molecules from the set for mixing, and then assigning
	//based on stoichiometry.

	double expectvol;
	int mincount;
	int maxcount;
	Index plindex; //placement index for multiple molecules
};

struct GASP2cell {
	double a,b,c; //ANGSTROMS
	double alpha,beta,gamma; //RADIANS
	double ratA,ratB,ratC; //UNIT-LESS
	Index label; //index to a string in "names";
	Index spacegroup; //spacegroup number (same as in IUCR tables)
	Vec3 maxFrac; //fractional limit of molecules;
	Vec3 minFrac;

};



//AML: NO POINTERS. We are not writing copy constructors.

/*
 * Workflows:
 *
 * Read from input:
 * parseXML finds the crystal root and passes it to the
 *
 *
 *
 */

class GASP2struct {
public:
	GASP2struct();
private:
	GASP2cell unit,builtunit,finalunit;
	vector<GASP2molecule> molecules; //contains native unsymmed molecules ONLY
	vector<GASP2molecule> builtcell; //contains symmed molecules from fitcell
	vector<GASP2molecule> finalcell; //contains the post optimization molecules;

	//AML:this is the only pointer allowed (not really a pointer)

	//It performs the optimization (say, a QE vc-relax);
	//what ever the function does, the only thing it is allowed to change
	//is the unit cell parameters and the position of atoms
	//changing other values (like the dihedral values) is performed by check()
	bool (*eval)(vector<GASP2molecule>&, GASP2cell&, double&);

public:


	bool fitcell(double distcell, double intradist); //puts the molecules in the built unit cell
	bool check(); //checks for violation of constraints (usually after opt)
	bool evaluate();
	void setEval(bool (*e)(vector<GASP2molecule>&, GASP2cell&, double&) ) {eval = e;};

	//also resets other parametewrs

	double getVolume(); //returns the unit cell volume

	tinyxml2::XMLElement* serializeXML(); //commits structure to XML element
	bool parseXML(tinyxml2::XMLDocument *doc, string& errorstring ); //reads structure info from an XML file to molecules
	bool parseXMLStruct(tinyxml2::XMLElement *elem, string& errorstring); //reads a structure that was transmitted
	bool readH5() {return true;};
	bool writeH5() {return true;};

	bool MPICopy(int target) {return true;}; //sends a struct via MPI to a client
	bool MPIRecv() {return true;}; //corresponding receive via MPI

	void initialize();

	//void mutate(double rate);
	//void cross(GASP2struct partner, GASP2struct &childA, GASP2struct &childB);



	bool completed() {return complete;};
	bool opted() {return didOpt;};
	bool fitcelled() {return didFitcell;};
	bool rejected() {return (finalstate != OKStruct);};
	void setOpt() { doOpt = true; };
	double getEnergy() {return finEnergy;};

private:
	//control flags
	bool doOpt;
	bool didFitcell;
	bool didOpt;
	bool complete; //true if opt completed and no further opt is required;
	StructError finalstate;
	Index crylabel;




	//final energy of system
	double energy;

	//phylogeny
	uuid_t ID;
	uuid_t parentA;
	uuid_t parentB;


private:
	//a vector of names shared by all GASP2structs
	//it doesn't make much sense to copy strings all the time
	//but only when necessary; these names are not shared with the
	//and therefore names should only accessed when running on the
	//server_program
	static vector<string> names;
	Index newName(const char* name);
	Index newName(string name);

	//spacegroup things
	//void setSymmOp(GASP2molecule &mol);
	//void enforceCrystalType();




};
