#include "gasp2common.hpp"

//using namespace std;

//static data
vector<Spgroup> spacegroups;
vector<string> spacegroupNames;

mt19937_64 rgen;

//this table is where the different spacegroup genes are encoded.
//there are only two restrictions, which are that
//D4d and D6d are not allowed. Those are caught in setSpacegroup
//as part of GASP2struct
vector<cryGroup> groupgenes = {
//Axis order 1
	{ Schoenflies::Cn, Axisnum::One, Centering::P, {1,} },
	{ Schoenflies::Cnv, Axisnum::One, Centering::P, {6,7} },
	{ Schoenflies::Cnv, Axisnum::One, Centering::C, {8,9} },
	{ Schoenflies::Cnh, Axisnum::One, Centering::P, {6,7} },
	{ Schoenflies::Cnh, Axisnum::One, Centering::C, {8,9} },
	{ Schoenflies::Sn, Axisnum::One, Centering::P, {6,7} },
	{ Schoenflies::Sn, Axisnum::One, Centering::C, {8,9} },
	{ Schoenflies::Dn, Axisnum::One, Centering::P, {3,4} },
	{ Schoenflies::Dn, Axisnum::One, Centering::C, {5,} },
	{ Schoenflies::Dnh, Axisnum::One, Centering::P, {25,26,27,28,29,30,31,32,33,34} },
	{ Schoenflies::Dnh, Axisnum::One, Centering::C, {35,36,37} },
	{ Schoenflies::Dnh, Axisnum::One, Centering::F, {42,43} },
	{ Schoenflies::Dnh, Axisnum::One, Centering::I, {44,45,46} },
	{ Schoenflies::Dnd, Axisnum::One, Centering::P, {10,11,13,14} },
	{ Schoenflies::Dnd, Axisnum::One, Centering::C, {12,15} },
//Axis order 2
	{ Schoenflies::Cn, Axisnum::Two, Centering::P, {3,4} },
	{ Schoenflies::Cn, Axisnum::Two, Centering::C, {5,} },
	{ Schoenflies::Cnv, Axisnum::Two, Centering::P, {25,26,27,28,29,30,31,32,33,34} },
	{ Schoenflies::Cnv, Axisnum::Two, Centering::C, {35,36,37} },
	{ Schoenflies::Cnv, Axisnum::Two, Centering::F, {42,43} },
	{ Schoenflies::Cnv, Axisnum::Two, Centering::I, {44,45,46} },
	{ Schoenflies::Cnh, Axisnum::Two, Centering::P, {10,11,13,14} },
	{ Schoenflies::Cnh, Axisnum::Two, Centering::C, {12,15} },
	{ Schoenflies::Sn, Axisnum::Two, Centering::P, {2,} },
	{ Schoenflies::Dn, Axisnum::Two, Centering::P, {16,17,18,19} },
	{ Schoenflies::Dn, Axisnum::Two, Centering::C, {20,21} },
	{ Schoenflies::Dn, Axisnum::Two, Centering::F, {22,} },
	{ Schoenflies::Dn, Axisnum::Two, Centering::I, {23,24} },
	{ Schoenflies::Dnh, Axisnum::Two, Centering::P, {47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62} },
	{ Schoenflies::Dnh, Axisnum::Two, Centering::C, {63,64,65,66,67,68} },
	{ Schoenflies::Dnh, Axisnum::Two, Centering::F, {69,70} },
	{ Schoenflies::Dnh, Axisnum::Two, Centering::I, {71,72,73,74} },
	{ Schoenflies::Dnd, Axisnum::Two, Centering::P, {111,112,113,114,115,116,117,118} },
	{ Schoenflies::Dnd, Axisnum::Two, Centering::I, {119,120,121,122} },
//Axis order 3
	{ Schoenflies::Cn, Axisnum::Three, Centering::P, {143,144,145} },
	{ Schoenflies::Cn, Axisnum::Three, Centering::R, {146,} },
	{ Schoenflies::Cnv, Axisnum::Three, Centering::P, {156,157,158,159} },
	{ Schoenflies::Cnv, Axisnum::Three, Centering::R, {160,161} },
	{ Schoenflies::Cnh, Axisnum::Three, Centering::P, {174,} },
	{ Schoenflies::Sn, Axisnum::Three, Centering::P, {174,} },
	{ Schoenflies::Dn, Axisnum::Three, Centering::P, {149,150,151,152,153,154} },
	{ Schoenflies::Dn, Axisnum::Three, Centering::R, {155,} },
	{ Schoenflies::Dnh, Axisnum::Three, Centering::P, {187,188,189,190} },
	{ Schoenflies::Dnd, Axisnum::Three, Centering::P, {162,163,164,165} },
	{ Schoenflies::Dnd, Axisnum::Three, Centering::R, {166,167} },
//Axis order 4
	{ Schoenflies::Cn, Axisnum::Four, Centering::P, {75,76,77,78} },
	{ Schoenflies::Cn, Axisnum::Four, Centering::I, {79,80} },
	{ Schoenflies::Cnv, Axisnum::Four, Centering::P, {99,100,101,102,103,104,105,106} },
	{ Schoenflies::Cnv, Axisnum::Four, Centering::I, {107,108,109,110} },
	{ Schoenflies::Cnh, Axisnum::Four, Centering::P, {83,84,85,86} },
	{ Schoenflies::Cnh, Axisnum::Four, Centering::I, {87,88} },
	{ Schoenflies::Sn, Axisnum::Four, Centering::P, {81,} },
	{ Schoenflies::Sn, Axisnum::Four, Centering::I, {82,} },
	{ Schoenflies::Dn, Axisnum::Four, Centering::P, {89,90,91,92,93,94,95,96} },
	{ Schoenflies::Dn, Axisnum::Four, Centering::I, {97,98} },
	{ Schoenflies::Dnh, Axisnum::Four, Centering::P, {123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138} },
	{ Schoenflies::Dnh, Axisnum::Four, Centering::I, {139,140,141,142} },
	//{ Schoenflies::Dnd, Axisnum::Four, Centering::P, {} }, not valid
//Axis order 6
	{ Schoenflies::Cn, Axisnum::Six, Centering::P, {168,169,170,171,172,173} },
	{ Schoenflies::Cnv, Axisnum::Six, Centering::P, {183,184,185,186} },
	{ Schoenflies::Cnh, Axisnum::Six, Centering::P, {175,176} },
	{ Schoenflies::Sn, Axisnum::Six, Centering::P, {147,} }, //C3i
	{ Schoenflies::Sn, Axisnum::Six, Centering::R, {148,} }, //C3i
	{ Schoenflies::Dn, Axisnum::Six, Centering::P, {177,178,179,180,181,182} },
	{ Schoenflies::Dnh, Axisnum::Six, Centering::P, {191,192,193,194} },
	//{ Schoenflies::Dnd, Axisnum::Six, Centering::P, {} }, not valid
//Extended groups (T and O groups)
	{ Schoenflies::T, Axisnum::UNK, Centering::P, {195,198} },
	{ Schoenflies::T, Axisnum::UNK, Centering::F, {196,} },
	{ Schoenflies::T, Axisnum::UNK, Centering::I, {197,199} },
	{ Schoenflies::Td, Axisnum::UNK, Centering::P, {215,218} },
	{ Schoenflies::Td, Axisnum::UNK, Centering::F, {216,219} },
	{ Schoenflies::Td, Axisnum::UNK, Centering::I, {217,220} },
	{ Schoenflies::Th, Axisnum::UNK, Centering::P, {200,201,205} },
	{ Schoenflies::Th, Axisnum::UNK, Centering::F, {202,203} },
	{ Schoenflies::Th, Axisnum::UNK, Centering::I, {204,206} },
	{ Schoenflies::O, Axisnum::UNK, Centering::P, {207,208,212,213} },
	{ Schoenflies::O, Axisnum::UNK, Centering::F, {209,210} },
	{ Schoenflies::O, Axisnum::UNK, Centering::I, {211,214} },
	{ Schoenflies::Oh, Axisnum::UNK, Centering::P, {221,222,223,224} },
	{ Schoenflies::Oh, Axisnum::UNK, Centering::F, {225,226,227,228} },
	{ Schoenflies::Oh, Axisnum::UNK, Centering::I, {229,230} },
};


double rad(double deg) {
	return deg/180.0*PI;
}

double deg(double rad) {
	return rad/PI*180.0;
}



double angle ( const Vec3 & A, const Vec3 & B, const Vec3 & C ) {
    double cos_angle = dot(A-B,C-B) / ( len(A-B) * len(C-B) );

    if(cos_angle > 1.0) cos_angle = 1.0;

    return ( acos(cos_angle) );
}



double dihedral ( const Vec3 & A, const Vec3 & B, const Vec3 & C, const Vec3 & D ) {
	Vec3 a,b,c;
	double x,y;
	a = B-A; b = C-B; c = D-C;
	x = dot(norm(b)*a,cross(b,c));
	y = dot(cross(a,b),cross(b,c));
	return atan2(x,y);
}


//AML: Used to re-orthogonalize a rotation matrix
//resulting from the multiplication of two rotation
//matrices. Taken from Direction Cosine Matrix IMU:
//Theory by William Premerlani and Paul Bizard
Mat3 stabilize(Mat3 m) {
	Vec3 x = m[0];
	Vec3 y = m[1];

	double error = dot(x,y);
	Vec3 xort = x-(error/2.0)*y;
	Vec3 yort = y-(error/2.0)*x;
	Vec3 zort = cross(xort,yort);

	Mat3 out;
	out[0] = 0.5*(3.0-dot(xort,xort))*xort;
	out[1] = 0.5*(3.0-dot(yort,yort))*yort;
	out[2] = 0.5*(3.0-dot(zort,zort))*zort;

	return out;
}


//AML: Adapted from Gabriel's code
//return a matrix that the columns 1th and 2th are 2 vectors that
//defines the molecular plane. if(det(matrix3Result==1) then OK
//else anything is wrong
Mat3 molecularPlane(const Vec3 at1, const Vec3 at2, const Vec3 at3)
		{

	Vec3 v1 = at2 - at1;
	Vec3 v2 = at3 - at1;
	Vec3 v3 = cross(v1,v2);

	if (len(v1) > 0.1 && len(v3) > 0.1) {
		v1 = norm(v1);
		v3 = norm(v3);
		v2 = cross(v3,v1);
	} else {
		if (len(v1) > 0.1)
			v1 = norm(v1);
		if (len(v2) > 0.1)
			v2 = norm(v2);
		if (len(v3) > 0.1)
			v3 = norm(v3);
	}

	Mat3 result;
	result[0] = v1;
	result[1] = v2;
	result[2] = v3;
	return result;
}



/// Reference: Cambridge Crystallographic Data Centre
double rcov( const Elem type ) {
    if ( type == Elem::H )
        return 0.23;
    else if ( type == Elem::C )
        return 0.68;
    else if ( type == Elem::N  )
        return 0.68;
    else if ( type == Elem::O  )
        return 0.68;
    else if ( type == Elem::F  )
        return 0.64;
    else if ( type == Elem::S  )
        return 1.02;
    else if ( type == Elem::Cl  )
        return 0.99;
    else if ( type == Elem::P  )
        return 1.05;
    else if ( type == Elem::Br  )
        return 1.21;
}


//Van der Waals radii are from
//Rowland, R.S., J. Phys. Chem. (2006) 100, pp 7384-7391
//excepting P, which is:
//Bondi, J. Phys. Chem. (1964) 68, pp 441
double vdw( Elem type ) {
    if ( type == Elem::H )
        return 1.09;
    else if ( type == Elem::C )
        return 1.75;
    else if ( type == Elem::N  )
        return 1.61;
    else if ( type == Elem::O  )
        return 1.56;
    else if ( type == Elem::F  )
        return 1.44;
    else if ( type == Elem::S  )
        return 1.79;
    else if ( type == Elem::Cl  )
        return 1.74;
    else if ( type == Elem::Br  )
        return 1.85;
    //else if ( type == I  )
    //    return 2.00;
    else if ( type == Elem::P  )
        return 1.80;

}

int indexSelect(double n, int size) {
	double i = 1.0 / static_cast<double>(size);
	int index = 0;
	double s = 0.0;
	while(true) {
		if(s > 1.0)
			break;
		if( (n > s) && (n < (s+i)) ) {
			return index;
		}
		else {
			s += i;
			index++;
		}
	}
	return 0;
}

Elem getElemType(string in) {
	Elem type;
	//UNK=0,C,H,N,O,P,Cl,F,S,Br,
	if(in == "C")
		type = Elem::C;
	else if(in == "H")
		type = Elem::H;
	else if(in == "N")
		type = Elem::N;
	else if(in == "O")
		type = Elem::O;
	else if(in == "P")
		type = Elem::P;
	else if(in == "Cl")
		type = Elem::Cl;
	else if(in == "F")
		type = Elem::F;
	else if(in == "S")
		type = Elem::S;
	else if(in == "Br")
		type = Elem::Br;
	else
		type = Elem::UNK;

	return type;
}

string getElemName(Elem in) {
	string type;
	//UNK=0,C,H,N,O,P,Cl,F,S,Br,
	if(in == Elem::C)
		type = "C";
	else if(in == Elem::H)
		type = "H";
	else if(in == Elem::N)
		type = "N";
	else if(in == Elem::O)
		type = "O";
	else if(in == Elem::P)
		type = "P";
	else if(in == Elem::Cl)
		type = "Cl";
	else if(in == Elem::F)
		type = "F";
	else if(in == Elem::S)
		type = "S";
	else if(in == Elem::Br)
		type = "Br";
	else
		type = "UNK";

	return type;


}

Axisnum getAxis(string in) {
	Axisnum type;
	if(in == string("One"))
		type = Axisnum::One;
	else if(in == "Two")
		type = Axisnum::Two;
	else if(in == "Three")
		type = Axisnum::Three;
	else if(in == "Four")
		type = Axisnum::Four;
	else if(in == "Six")
		type = Axisnum::Six;
	else
		type = Axisnum::UNK;
	return type;
}

string getAxis(Axisnum in) {
	string type;
	if(in == Axisnum::One)
		type = "One";
	else if(in == Axisnum::Two)
		type = "Two";
	else if(in == Axisnum::Three)
		type = "Three";
	else if(in == Axisnum::Four)
		type = "Four";
	else if(in == Axisnum::Six)
		type = "Six";
	else
		type = "UNK";

	return type;
}

Schoenflies getSchoenflies(string in) {
	Schoenflies type;
	if(in == "Cn")
		type = Schoenflies::Cn;
	else if(in == "Cnv")
		type = Schoenflies::Cnv;
	else if(in == "Cnh")
		type = Schoenflies::Cnh;
	else if(in == "Sn")
		type = Schoenflies::Sn;
	else if(in == "Dn")
		type = Schoenflies::Dn;
	else if(in == "Dnd")
		type = Schoenflies::Dnd;
	else if(in == "Dnh")
		type = Schoenflies::Dnh;
	else if(in == "T")
		type = Schoenflies::T;
	else if(in == "Th")
		type = Schoenflies::Th;
	else if(in == "O")
		type = Schoenflies::O;
	else if(in == "Td")
		type = Schoenflies::Td;
	else if(in == "Oh")
		type = Schoenflies::Oh;
	else
		type = Schoenflies::UNK;

	return type;
}

string getSchoenflies(Schoenflies in) {
	string type;
	if(in == Schoenflies::Cn)
		type = "Cn";
	else if(in == Schoenflies::Cnv)
		type = "Cnv";
	else if(in == Schoenflies::Cnh)
		type = "Cnh";
	else if(in == Schoenflies::Sn)
		type = "Sn";
	else if(in == Schoenflies::Dn)
		type = "Dn";
	else if(in == Schoenflies::Dnd)
		type = "Dnd";
	else if(in == Schoenflies::Dnh)
		type = "Dnh";
	else if(in == Schoenflies::T)
		type = "T";
	else if(in == Schoenflies::Th)
		type = "Th";
	else if(in == Schoenflies::O)
		type = "O";
	else if(in == Schoenflies::Td)
		type = "Td";
	else if(in == Schoenflies::Oh)
		type = "Oh";
	else
		type = "UNK";
	return type;
}

Lattice getLattice(string in) {
	Lattice type;
	if(in == "cubic")
		type = Lattice::Cubic;
	else if(in == "tetragonal")
		type = Lattice::Tetragonal;
	else if(in == "orthorhombic")
		type = Lattice::Orthorhombic;
	else if(in == "hexagonal")
		type = Lattice::Hexagonal;
	else if(in == "monoclinic")
		type = Lattice::Monoclinic;
	else if(in == "triclinic")
		type = Lattice::Triclinic;
	else if(in == "rhombohedral")
		type = Lattice::Rhombohedral;
	else
		type = Lattice::UNK;

	return type;
}

vector<string> split(string in, char delim) {
	vector<string> out;
	istringstream ss(in);
	while(getline(ss,in, delim))
		out.push_back(in);
	return out;
}

//libUUID wrapper class stuff
UUID::UUID(void) {
		generate();
};
UUID::UUID(const UUID &u) {
		uuid_copy(this->uuid, u.uuid);
}
void UUID::generate() {
		uuid_generate(uuid);
};
void UUID::clear() {
		uuid_clear(uuid);
}
string UUID::toStr() {
		char buff[36];
		uuid_unparse(uuid, buff);
		string s = buff;
		return s;
}
void UUID::frStr(string in) {
		char buff[37];
		in.copy(buff, 36);
		buff[36] = '\0';
		uuid_parse(buff, uuid);
}
UUID& UUID::operator=(UUID u) {
		uuid_copy(this->uuid, u.uuid);
		return *this;
}

bool operator==(const UUID &u, const UUID &v) {
	return (uuid_compare(u.uuid,v.uuid) == 0);
}


bool loadSpaceGroups() {
	string strtemp;
	const char* stemp;
	Mat3 rot; Vec3 trans;

	tinyxml2::XMLDocument spacedoc;
	tinyxml2::XMLElement *group;

	Spgroup spg;
	//make empty groups in a couple of places;
	spacegroups.push_back(spg);
	spacegroupNames.push_back("NULL");

    const char * s = &_binary_spacegroups_xml_start;
    const char * e = &_binary_spacegroups_xml_end;
    size_t size = e-s;
    if(spacedoc.Parse(s, size)) {
    	cout << "A problem was encountered when loading the spacegroup file!\n";
    	return false;
    }
    tinyxml2::XMLElement *elem = spacedoc.FirstChildElement("spglist")->FirstChildElement("spg");
    if(!elem) {
    	cout << "A problem was encountered when loading the spacegroup file!\n";
    	return false;
    }

    while(elem) {
    	stemp = elem->Attribute("sym");
		if(stemp) {
			strtemp = stemp;
			spacegroupNames.push_back(strtemp);
		}
		spg.R.clear();
		spg.T.clear();
    	stemp = elem->Attribute("sys");
		if(stemp) {
			strtemp = stemp;
			spg.L = getLattice(strtemp);
			if(spg.L == Lattice::UNK) {
				cout << "An unknown lattice type was specified in the spacegroup file!\n";
				return false;
			}
		}
		group = elem->FirstChildElement("op");
		while(group) {
			strtemp = group->GetText();
			parseSymop(strtemp, rot, trans);
			spg.R.push_back(rot);
			spg.T.push_back(trans);
			group = group->NextSiblingElement("op");
		}
		spacegroups.push_back(spg);
		elem = elem->NextSiblingElement("spg");

    }

    cout << endl;
    return true;

}

//adapted from the original MGAC code (written by
void parseSymop ( string symm, Mat3 &symmR, Vec3 &symmT ) {

    double    aux_coord;
    Vec3    aux_vect;
    vector<string> symms;

	symms = split(symm, ',');


	for (int j = 0 ; j < 3 ; j++ ) {

		aux_vect = vl_0;
		aux_coord = 0.0;

		if ( symms[j].find("-x") != string::npos )
			aux_vect[0] = (-1.0);
		else if ( ( symms[j].find("+x") != string::npos ) || ( symms[j].find("x") != string::npos ) )
			aux_vect[0] = (1.0);
		else
			aux_vect[0] = (0.0);

		if ( symms[j].find("-y") != string::npos )
			aux_vect[1] = (-1.0);
		else if ( ( symms[j].find("+y") != string::npos ) || ( symms[j].find("y") != string::npos ) )
			aux_vect[1] = (1.0);
		else
			aux_vect[1] = (0.0);

		if ( symms[j].find("-z") != string::npos )
			aux_vect[2] = (-1.0);
		else if ( ( symms[j].find("z") != string::npos ) || ( symms[j].find("z") != string::npos ) )
			aux_vect[2] = (1.0);
		else
			aux_vect[2] = (0.0);

		if ( symms[j].find("+1/2") != string::npos )
			aux_coord = 1.0/2.0;
		else if ( symms[j].find("-1/2") != string::npos )
			aux_coord = -1.0/2.0;
		else if ( symms[j].find("+1/3") != string::npos )
			aux_coord = 1.0/3.0;
		else if ( symms[j].find("+2/3") != string::npos )
			aux_coord = 2.0/3.0;
		else if ( symms[j].find("+1/4") != string::npos )
			aux_coord = 1.0/4.0;
		else if ( symms[j].find("-1/4") != string::npos )
			aux_coord = -1.0/4.0;
		else if ( symms[j].find("+3/4") != string::npos )
			aux_coord = 3.0/4.0;
		else if ( symms[j].find("+5/4") != string::npos )
			aux_coord = 5.0/4.0;
		else if ( symms[j].find("+1/6") != string::npos )
			aux_coord = 1.0/6.0;
		else if ( symms[j].find("-1/6") != string::npos )
			aux_coord = -1.0/6.0;
		else if ( symms[j].find("+5/6") != string::npos )
			aux_coord = 5.0/6.0;
		else if ( symms[j].find("+7/6") != string::npos )
			aux_coord = 7.0/6.0;
		else
			aux_coord = 0.0;

		/* Set matrix R and vector T parameters */
		symmR[j] = aux_vect;
		symmT[j] = aux_coord;

	}

}

string tfconv(bool var) {
	if(var)
		return "true";
	else
		return "false";
}


