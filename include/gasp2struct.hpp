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
	Vec3 pos; //always saved as absolute angstrom coord
	Elem type;
	NIndex label; //index to a string in "names";
};

struct GASP2bond {
	NIndex label;
	Index a,b; //indices to atoms in atom VECTOR, not the name index
	double minLen, maxLen; //ANGSTROMS
	double len; //ANGSTROMS
};

struct GASP2angle {
	NIndex label;
	Index a,b,c;
	double maxAng,minAng; //DEGREES
	double ang; //DEGREES
};

struct GASP2dihedral {
	NIndex label;
	Index a,b,c,d; //indices to atoms;
	vector<Index> update; //list of atoms to update
	double maxAng, minAng; //DEGREES
	double ang; //DEGREES - XMLOUT
	void clear() {
		update.clear();
	}
};

struct GASP2molecule {
	Vec3 pos; //ALWAYS fractional, refers to centroid XMLOUT
	Mat3 rot; //the rotation matrix of the molecule XMLOUT
	//Why rotation matrix? Space is cheap (realtively speaking)
	//
	vector<GASP2atom> atoms;
	Index p1,p2,p3; //atom plane indices

	NIndex label; //index to a string in "names";


	Index symm; //index of symmetry operation, 0 is always x,y,z
	Mat3 symmR; //given symmop rotation; INTERNAL
	Vec3 symmT; //given symmop translation; INTERNAL

	//bonds and angles are constraints;
	//at this time they only play a role in post-opt checking
	//dihedrals are modified during the fitcell operation
	vector<GASP2bond> bonds;
	vector<GASP2angle> angles;
	vector<GASP2dihedral> dihedrals;

	void clear() {
		atoms.clear();
		dihedrals.clear();
		angles.clear();
		bonds.clear();
	}
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
	Index plindex; //placement index for multiple molecules
};

//stoichiometric composition
//handles a single molecule
//
struct GASPcomp {
	NIndex mol; //refers to the index of the molecule
	int min, max;
	int count;
};

struct GASP2cell {
	double a,b,c; //ANGSTROMS
	double alpha,beta,gamma; //DEGREES
	double ratA,ratB,ratC; //UNIT-LESS
	Index spacegroup; //spacegroup number (same as in IUCR tables)
	Vec3 maxFrac; //fractional limit of molecules; INTERNAL
	Vec3 minFrac; //INTERNAL
	vector<GASPcomp> stoich;
	void clear() {
		stoich.clear();
	}

};



//AML: NO POINTERS. We are not writing copy constructors.

class GASP2struct {
public:
	GASP2struct();
private:

	//main data constructs
	GASP2cell unit;
	vector<GASP2molecule> molecules;

	//evaluator function pointer (only pointer allowed!)
	//It performs the optimization (say, a QE vc-relax);
	//what ever the function does, the only thing it is allowed to change
	//is the unit cell parameters and the position of atoms
	//changing other values (like the dihedral values) is performed internally by check()
	bool (*eval)(vector<GASP2molecule>&, GASP2cell&, double&);

public:

	//structure manipulators
	bool fitcell(); //puts the molecules in the built unit cell
	bool unfitcell(); //reduces the structure to only fundamental units
	bool check(); //checks for violation of constraints (usually after opt)
	bool evaluate();
	void setEval(bool (*e)(vector<GASP2molecule>&, GASP2cell&, double&) ) {eval = e;};

	//I/O handlers
	tinyxml2::XMLElement* serializeXML(); //commits structure to XML element
	bool parseXMLDoc(tinyxml2::XMLDocument *doc, string& errorstring ); //reads structure info from an XML file to molecules
	bool parseXMLStruct(tinyxml2::XMLElement *elem, string& errorstring); //reads a structure that was transmitted
	//bool readH5() {return true;};
	//bool writeH5() {return true;};
	void logStruct();

	//MPI handlers
	bool MPICopy(int target) {return true;}; //sends a struct via MPI to a client
	bool MPIRecv() {return true;}; //corresponding receive via MPI

	//genetic operators
	//void initialize();
	//void mutate(double rate);
	//void cross(GASP2struct partner, GASP2struct &childA, GASP2struct &childB);

	//getters
	bool completed() {return complete;};
	bool opted() {return didOpt;};
	bool fitcelled() {return isFitcell;};
	bool rejected() {return (finalstate != OKStruct);};
	double getEnergy() {return energy;};
	double getVolume(); //returns the unit cell volume

private:
	//values
	double interdist;
	double intradist;
	double maxvol;
	double minvol;


	//control flags
	bool isFitcell;
	bool didOpt;
	bool complete; //true if opt completed and no further opt is required;
	StructError finalstate;
	NIndex crylabel;
	time_t time;
	int steps;

	//final energy of system
	double energy;
	double force;
	double pressure;

	//identifiers
	UUID ID;
	UUID parentA;
	UUID parentB;


private:
	//a vector of names shared by all GASP2structs
	//it doesn't make much sense to copy strings all the time
	//but only when necessary; these names are not shared with the
	//and therefore names should only accessed when running on the
	//server_program
	static vector<string> names;
	NIndex newName(const char* name);
	NIndex newName(string name);
	//Looks up an atom index
	Index atomLookup(NIndex nameInd, GASP2molecule mol);


	//spacegroup things
	void setSymmOp(GASP2molecule &mol);
	void enforceCrystalType();

	//parsing auxiliaries
	bool readAtom(tinyxml2::XMLElement *elem, string& errorstring, GASP2atom &at);
	bool readDih(tinyxml2::XMLElement *elem, string& errorstring, GASP2dihedral &dih, GASP2molecule mol);
	bool readBond(tinyxml2::XMLElement *elem, string& errorstring, GASP2bond &bond, GASP2molecule mol);
	bool readAngle(tinyxml2::XMLElement *elem, string& errorstring, GASP2angle &ang, GASP2molecule mol);
	bool readMol(tinyxml2::XMLElement *elem, string& errorstring, GASP2molecule &mol);
	bool readCell(tinyxml2::XMLElement *elem, string& errorstring, GASP2cell &cell);
	bool readInfo(tinyxml2::XMLElement *elem, string& errorstring);




};
