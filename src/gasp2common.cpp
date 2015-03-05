#include "gasp2common.hpp"

using namespace std;

//static datas
tinyxml2::XMLDocument spacegroups;
vector<string> spacegroupNames;


//double dist ( const Vec3 & A, const Vec3 & B ) {
//    return sqrt( (A-B).dot(A-B) );
//}


///*** Plane angle ***///

double angle ( const Vec3 & A, const Vec3 & B, const Vec3 & C ) {
    double cos_angle = dot(A-B,C-B) / ( len(A-B) * len(C-B) );

    if(cos_angle > 1.0) cos_angle = 1.0;

    return ( acos(cos_angle) );
}

///*** Dihedral angle ***///

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
    if ( type == H )
        return 0.23;
    else if ( type == C )
        return 0.68;
    else if ( type == N  )
        return 0.68;
    else if ( type == O  )
        return 0.68;
    else if ( type == F  )
        return 0.64;
    else if ( type == S  )
        return 1.02;
    else if ( type == Cl  )
        return 0.99;
    else if ( type == P  )
        return 1.05;
    else if ( type == Br  )
        return 1.21;
}


//Van der Waals radii are from
//Rowland, R.S., J. Phys. Chem. (2006) 100, pp 7384-7391
//excepting P, which is:
//Bondi, J. Phys. Chem. (1964) 68, pp 441
double vdw( Elem type ) {
    if ( type == H )
        return 1.09;
    else if ( type == C )
        return 1.75;
    else if ( type == N  )
        return 1.61;
    else if ( type == O  )
        return 1.56;
    else if ( type == F  )
        return 1.44;
    else if ( type == S  )
        return 1.79;
    else if ( type == Cl  )
        return 1.74;
    else if ( type == Br  )
        return 1.85;
    //else if ( type == I  )
    //    return 2.00;
    else if ( type == P  )
        return 1.80;

}

Elem getElemType(string in) {
	Elem type;
	//UNK=0,C,H,N,O,P,Cl,F,S,Br,
	if(in == "C")
		type = C;
	else if(in == "H")
		type = H;
	else if(in == "N")
		type = N;
	else if(in == "O")
		type = O;
	else if(in == "P")
		type = P;
	else if(in == "Cl")
		type = Cl;
	else if(in == "F")
		type = F;
	else if(in == "S")
		type = S;
	else if(in == "Br")
		type = Br;
	else
		type = UNK;

	return type;
}

string getElemName(Elem in) {
	string type;
	//UNK=0,C,H,N,O,P,Cl,F,S,Br,
	if(in == C)
		type = "C";
	else if(in == H)
		type = "H";
	else if(in == N)
		type = "N";
	else if(in == O)
		type = "O";
	else if(in == P)
		type = "P";
	else if(in == Cl)
		type = "Cl";
	else if(in == F)
		type = "F";
	else if(in == S)
		type = "S";
	else if(in == Br)
		type = "Br";
	else
		type = "UNK";

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

    const char * s = &_binary_spacegroups_xml_start;
    const char * e = &_binary_spacegroups_xml_end;
    size_t size = e-s;
    if(spacegroups.Parse(s, size)) {
    	cout << "A problem was encountered when loading the spacegroup file!\n";
    	return false;
    }
    tinyxml2::XMLElement *elem = spacegroups.FirstChildElement("spglist")->FirstChildElement("spg");
    if(!elem) {
    	cout << "A problem was encountered when loading the spacegroup file!\n";
    	return false;
    }
    while(elem) {
    	stemp = elem->Attribute("sym");
		if(stemp) {
			strtemp = stemp;
			//cout << "Sym is " << strtemp << " ";
			spacegroupNames.push_back(strtemp);
		}
		elem = elem->NextSiblingElement("spg");
    }

    cout << endl;
    return true;

}




