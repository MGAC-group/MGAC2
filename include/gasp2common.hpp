#define _GLIBCXX_USE_NANOSLEEP
#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <ctime>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include "tinyxml2.h"
#include "svl/SVL.h"
#include "mpi.h"
#include <uuid/uuid.h>
#include <random>
#include <algorithm>
#include <limits>
#include <future>
#include <thread>
#include <fcntl.h>
#include <signal.h>
#include <sys/wait.h>

#define READ 0
#define WRITE 1

extern "C" {
	extern char _binary_spacegroups_xml_start;
	extern char _binary_spacegroups_xml_end;
}




#define PI 3.1415926535897932384626433

using namespace std;




typedef unsigned int Index;
typedef unsigned int NIndex;

typedef enum struct Elem {
	UNK=0,C,H,N,O,P,Cl,F,S,Br,
}Elem;

typedef enum struct Centering {
	UNK=0, P, C, I, F, R, H,
}Centering;

typedef enum struct Schoenflies {
	 UNK=0,Cn, Cnv, Cnh, Sn, Dn, Dnh, Dnd,
	T, Th, O, Td, Oh,
}Schoenflies;

typedef enum struct Axisnum {
	UNK=0,One, Two, Three, Four, Six,
}Axisnum;

typedef struct cryGroup {
	Schoenflies s;
	Axisnum a;
	Centering c; //centering
	vector<int> indices;

}cryGroup;

typedef enum Spacemode {
	Limited=0, //Excludes Tetrahedral/Octahedral schoenflies and F centerings
	Single, //Single spacegroup only
	Full, //All spacegroups (including T,O and F)
}Spacemode;


typedef enum struct Lattice {
	UNK=0,
	Cubic,
	Tetragonal,
	Orthorhombic,
	Hexagonal,
	Rhombohedral,
	Monoclinic,
	Triclinic,
}Lattice;

typedef struct Spgroup {
	vector<Mat3> R;
	vector<Vec3> T;
	Lattice L;

}Symops;

typedef double Centerfl; //encodes P C/A/B/N/D/E I, F, R, H depending on schoenflies/axis
typedef double Subtype; //encodes glide and screw axes based on group number

//EXTERNALS
extern vector<Spgroup> spacegroups;
extern vector<string> spacegroupNames;
extern mt19937_64 rgen;
extern vector<cryGroup> groupgenes;

double rad(double deg);
double deg(double rad);

double angle ( const Vec3 & A, const Vec3 & B, const Vec3 & C );

double dihedral ( const Vec3 & A, const Vec3 & B, const Vec3 & C, const Vec3 & D );

Mat3 stabilize(Mat3 m);

double rcov( const Elem type );

double vdw( Elem type );

int indexSelect(double n, int size);

Vec3 modVec3 (Vec3 in);

Elem getElemType(string in);

string getElemName(Elem in);

Axisnum getAxis(string in);
string getAxis(Axisnum in);

Lattice getLattice(string in);

Schoenflies getSchoenflies(string in);
string getSchoenflies(Schoenflies in);

vector<string> split(string in, char delim=' ');


bool loadSpaceGroups();
void parseSymop ( string symm, Mat3 &symmR, Vec3 &symmT );

string tfconv(bool var);

//special stuff
FILE * popen2(const char *command, pid_t &pid);
int pclose2(pid_t pid, FILE *outfp);

//wrapper for libuuid
class UUID {
private:
        uuid_t uuid;
public:
        UUID(void);
        UUID(const UUID &u);
        void generate();
        void clear();
        string toStr();
        void frStr(string in);
        UUID& operator=(UUID u);
        friend bool operator==(const UUID &u, const UUID &v);
};


