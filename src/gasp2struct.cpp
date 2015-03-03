#include "gasp2struct.hpp"


vector<string> GASP2struct::names;

GASP2struct::GASP2struct() {
	molecules.clear();
	isFitcell = false;
	didOpt = false;
	complete = false;
	finalstate = OKStruct;
	crylabel = 0;
	time = 0;
	steps = 0;

	energy = 0.0;
	force = 0.0;
	pressure = 0.0;
	interdist = 0.0;
	intradist = 0.0;
	maxvol = 0.0;
	minvol = 0.0;

	ID.generate();
	parentA.clear();
	parentB.clear();



}

double GASP2struct::getVolume() {
	double phi = 2 * cos(unit.alpha) * cos(unit.beta) * cos(unit.gamma);
    return ( 2 * unit.a * unit.b * unit.c *
             sqrt(1 - cos(unit.alpha)*cos(unit.alpha)
            		- cos(unit.beta)*cos(unit.beta)
            		- cos(unit.gamma)*cos(unit.gamma)
            		+ phi) );

}

bool GASP2struct::fitcell() {


	return true;
}

bool GASP2struct::check() {
	Index a,b,c,d;
	double val;

	for(int i = 0; i < molecules.size(); i++) {
		//bonds
		for(int j = 0; j < molecules[i].bonds.size(); j++) {
			a = molecules[i].bonds[j].a;
			b = molecules[i].bonds[j].b;
			val = len(molecules[i].atoms[a].pos - molecules[i].atoms[b].pos);
			if(val > molecules[i].bonds[j].maxLen || val < molecules[i].bonds[j].minLen)
				finalstate = OptBadBond;
			molecules[i].bonds[j].len = val;
		}
		//angles
		for(int j = 0; j < molecules[i].angles.size(); j++) {
			a = molecules[i].angles[j].a;
			b = molecules[i].angles[j].b;
			c = molecules[i].angles[j].c;
			val = angle(molecules[i].atoms[a].pos,molecules[i].atoms[b].pos,molecules[i].atoms[c].pos);
			if(val > molecules[i].angles[j].maxAng || val < molecules[i].angles[j].minAng)
				finalstate = OptBadAng;
			molecules[i].angles[j].ang = val;
		}
		//dihedrals
		for(int j = 0; j < molecules[i].dihedrals.size(); j++) {
			a = molecules[i].dihedrals[j].a;
			b = molecules[i].dihedrals[j].b;
			c = molecules[i].dihedrals[j].c;
			d = molecules[i].dihedrals[j].d;
			val = dihedral(molecules[i].atoms[a].pos,molecules[i].atoms[b].pos,
					molecules[i].atoms[c].pos,molecules[i].atoms[d].pos);
			if(val > molecules[i].dihedrals[j].maxAng || val < molecules[i].dihedrals[j].minAng)
				finalstate = OptBadDih;
			molecules[i].dihedrals[j].ang = val;
		}
		//plane rotations


		//position

	}
	return true;
}

bool GASP2struct::evaluate() {
	if(!isFitcell) {
		finalstate = NoFitcell;
		return false;
	}

	//never evaluate something that has already been evaluated
	if(!didOpt) {
		complete = eval(molecules, unit, energy);
		didOpt = true;
	}

	return complete;
};

//searches for a name in the namelist
//if the name is not found, the name is added.
NIndex GASP2struct::newName(string name) {
	NIndex n = names.size();
	for(int i = 0; i < n; i++) {
		if(name == names[i])
			return i;
	}
	names.push_back(name);
	return n;
}

NIndex GASP2struct::newName(const char * name) {
	string temp = name;
	return newName(temp);
}

//searches for a name index that matches the given index
//in the atoms of a molecule. returns -1 if not found.
Index GASP2struct::atomLookup(NIndex nameInd, GASP2molecule mol) {
	for(int i = 0; i < mol.atoms.size(); i++) {
		if(nameInd == mol.atoms[i].label)
			return i;
	}
	return -1;
}


bool GASP2struct::parseXMLDoc(tinyxml2::XMLDocument *doc, string& errorstring ) {

	tinyxml2::XMLElement *crystal = doc->FirstChildElement("mgac")->FirstChildElement("crystal");
	return parseXMLStruct(crystal, errorstring);

}

bool GASP2struct::parseXMLStruct(tinyxml2::XMLElement * crystal, string & errorstring) {

	double dtemp;
	const char * stemp;
	string strtemp;
	int itemp;

	//parse the crystal
	//at this point we assume only one crystal in the specification
	//multiple spacegroups will likely be handled by a delimited list
	//shoudl that mode of operation ever be implemented

	if(!crystal) {
		errorstring = "No crystal structure specified!\n";
		return false;
	}
	else {
		//label
		if(! crystal->QueryIntAttribute("name",&itemp) ) {
			if(itemp > 0)
				crylabel = itemp;
		} else {
			stemp = crystal->Attribute("name");
			if(stemp)
				crylabel =  newName(stemp);
			else {
				errorstring = "A crystal name was not specified!\n";
				return false;
			}
		}


		//interdist
		if(!crystal->QueryDoubleAttribute("interdist", &dtemp)) {
			interdist = dtemp;
		}
		else {
			errorstring = "interdist is required but not specified!.\n";
			return false;
		}
		if(interdist <= 0.0) {
			errorstring = "interdist is out of bounds (interdist > 0.0).\n";
			return false;
		}

		//intradist
		if(!crystal->QueryDoubleAttribute("intradist", &dtemp)) {
			intradist = dtemp;
		}
		else {
			errorstring = "intradist is required but not specified!.\n";
			return false;
		}
		if(intradist <= 0.0) {
			errorstring = "intradist is out of bounds (intradist > 0.0).\n";
			return false;
		}

		//maxvol
		if(!crystal->QueryDoubleAttribute("maxvol", &dtemp)) {
			maxvol = dtemp;
		}
		else {
			errorstring = "maxvol is required but not specified!.\n";
			return false;
		}

		//minvol
		if(!crystal->QueryDoubleAttribute("minvol", &dtemp)) {
			minvol = dtemp;
		}
		else {
			errorstring = "minvol is required but not specified!.\n";
			return false;
		}
		if(minvol > maxvol) {
			errorstring = "minvol is out of bounds (minvol > maxvol).\n";
			return false;
		}

		//parse the molecules
		GASP2molecule tempmol;
		tinyxml2::XMLElement *molecule, *item;
		molecule = crystal->FirstChildElement("molecule");
		if(!molecule) {
			errorstring = "No molecules specified in the structure!\n";
			return false;
		}
		while(molecule) {
			if(!readMol(molecule, errorstring, tempmol))
				return false;
			//push the molecule(s) onto the list
			this->molecules.push_back(tempmol);
			tempmol.clear();
			molecule = molecule->NextSiblingElement("molecule");
		}
		//parse the unitcell info (if present)
		GASP2cell tempcell;

		item = crystal->FirstChildElement("cell");
		if(item) {
			if(!readCell(item, errorstring, tempcell))
				return false;
			else
				unit = tempcell;
		}

		//parse the structure info (if present)

		item = crystal->FirstChildElement("info");
		if(item) {
			if(!readInfo(item, errorstring))
				return false;
		}
	}

	return true;
}

bool GASP2struct::readAtom(tinyxml2::XMLElement *elem, string& errorstring, GASP2atom &at) {

	const char * stemp;
	int itemp;
	string strtemp;
	double dtemp;


	if(elem) {
		//atom elem
		stemp = elem->Attribute("elem");
		if(stemp) {
			strtemp = stemp;
			at.type = getElemType(strtemp);
		}
		else {
			errorstring = "An atom element was not specified!\n";
			return false;
		}
		if(at.type == UNK) {
			errorstring = "An unknown or unsupported element name was used!\n";
			return false;
		}

		//atom label
		//first search for an index, then parse as string
		if(! elem->QueryIntAttribute("title",&itemp) ) {
			if(itemp > 0)
				at.label = itemp;
		} else {
			stemp = elem->Attribute("title");
			if(stemp)
				at.label =  newName(stemp);
			else {
				errorstring = "An atom title name was not specified!\n";
				return false;
			}
		}

		//x coord (fractional)
		if(!elem->QueryDoubleAttribute("x", &dtemp)) {
			at.pos[0] = dtemp;
		}
		else {
			errorstring = "No value was given for an X coordinate!\n";
			return false;
		}
		//y coord (fractional)
		if(!elem->QueryDoubleAttribute("y", &dtemp)) {
			at.pos[1] = dtemp;
		}
		else {
			errorstring = "No value was given for a Y coordinate!\n";
			return false;
		}
		//z coord (fractional)
		if(!elem->QueryDoubleAttribute("z", &dtemp)) {
			at.pos[2] = dtemp;
		}
		else {
			errorstring = "No value was given for a Z coordinate!\n";
			return false;
		}

	return true;
	}
	else {
		errorstring = "Bad atom! Something is wrong...";

		return false;
	}

}

bool GASP2struct::readDih(tinyxml2::XMLElement *elem, string& errorstring, GASP2dihedral &dih, GASP2molecule mol) {
	int itemp;
	const char * stemp;
	string strtemp;
	vector<string> vstemp;
	double dtemp;


	//title
	//first search for an index, then parse as string
	if(! elem->QueryIntAttribute("title",&itemp) ) {
		if(itemp > 0)
			dih.label = itemp;
	} else {
		stemp = elem->Attribute("title");
		if(stemp)
			dih.label =  newName(stemp);
		else {
			errorstring = "A dihedral title name was not specified!\n";
			return false;
		}
	}

	//angle
	stemp = elem->Attribute("angle");
	if(stemp) {
		strtemp = stemp;
		vstemp = split(strtemp);
		int size = vstemp.size();
		if(size != 4) {
			errorstring = "The dihedral angle specification for "+names[dih.label]+" does not have 4 elements!\n";
			return false;
		}
		dih.a = atomLookup(newName(vstemp[0]), mol);
		dih.b = atomLookup(newName(vstemp[1]), mol);
		dih.c = atomLookup(newName(vstemp[2]), mol);
		dih.d = atomLookup(newName(vstemp[3]), mol);
		if(dih.a < 0 || dih.b < 0 || dih.c < 0 || dih.d < 0) {
			errorstring = "An atom designation in dihedral "+names[dih.label]+" did not match any atom in the molecule!\n";
			return false;
		}
	}

	//update
	dih.update.clear();
	stemp = elem->Attribute("angle");
	if(stemp) {
		strtemp = stemp;
		vstemp = split(strtemp);
		int size = vstemp.size();
		if(size < 1) {
			errorstring = "The dihedral angle update specification for "+names[dih.label]+" does not have any elements!\n";
			return false;
		}
		for(int i = 0; i < size; i++) {
			Index t = atomLookup(newName(vstemp[i]), mol);
			if(t < 0) {
				errorstring = "An atom designation in dihedral update "+names[dih.label]+" did not match any atom in the molecule!\n";
				return false;
			}
			dih.update.push_back(t);
		}
	}

	//min
	dih.minAng = -180.0;
	if(!elem->QueryDoubleAttribute("min", &dtemp)) {
		dih.minAng = dtemp;
	}
	if(dih.minAng < -180.0 || dih.minAng > 180.0) {
		errorstring = "The minimum angle in dihedral "+names[dih.label]+" is out of bounds (-180.0 < ang < 180.0).\n";
		return false;
	}

	//max
	dih.maxAng = 180.0;
	if(!elem->QueryDoubleAttribute("max", &dtemp)) {
		dih.maxAng = dtemp;
	}
	if(dih.maxAng < -180.0 || dih.maxAng > 180.0) {
		errorstring = "The maximum angle in dihedral "+names[dih.label]+" is out of bounds (-180.0 < ang < 180.0).\n";
		return false;
	}


	//val
	dih.ang = 0.0;
	if(!elem->QueryDoubleAttribute("val", &dtemp)) {
		dih.ang = dtemp;
	}
	if(dih.ang < -180.0 || dih.ang > 180.0) {
		errorstring = "The angle value in dihedral "+names[dih.label]+" is out of bounds (-180.0 < ang < 180.0).\n";
		return false;
	}

	return true;

}

bool GASP2struct::readBond(tinyxml2::XMLElement *elem, string& errorstring, GASP2bond &bond, GASP2molecule mol) {
	int itemp;
	const char * stemp;
	string strtemp;
	vector<string> vstemp;
	double dtemp;


	//title
	//first search for an index, then parse as string
	if(! elem->QueryIntAttribute("title",&itemp) ) {
		if(itemp > 0)
			bond.label = itemp;
	} else {
		stemp = elem->Attribute("title");
		if(stemp)
			bond.label =  newName(stemp);
		else {
			errorstring = "A bond title name was not specified!\n";
			return false;
		}
	}

	//atoms
	stemp = elem->Attribute("atoms");
	if(stemp) {
		strtemp = stemp;
		vstemp = split(strtemp);
		int size = vstemp.size();
		if(size != 2) {
			errorstring = "The bond specification for "+names[bond.label]+" does not have 2 elements!\n";
			return false;
		}
		bond.a = atomLookup(newName(vstemp[0]), mol);
		bond.b = atomLookup(newName(vstemp[1]), mol);
		if(bond.a < 0 || bond.b < 0) {
			errorstring = "An atom designation in bond "+names[bond.label]+" did not match any atom in the molecule!\n";
			return false;
		}
	}

	//min
	bond.minLen = -180.0;
	if(!elem->QueryDoubleAttribute("min", &dtemp)) {
		bond.minLen = dtemp;
	}
	if(bond.minLen <= 0.0) {
		errorstring = "The minimum length in bond "+names[bond.label]+" is out of bounds (len > 0.0).\n";
		return false;
	}

	//max
	bond.maxLen = 180.0;
	if(!elem->QueryDoubleAttribute("max", &dtemp)) {
		bond.maxLen = dtemp;
	}
	if(bond.maxLen <= 0.0 || bond.maxLen < bond.minLen) {
		errorstring = "The maximum length in bond "+names[bond.label]+" is out of bounds (len > 0.0, max > min).\n";
		return false;
	}

	//val
	bond.len = 0.0;
	if(!elem->QueryDoubleAttribute("val", &dtemp)) {
		bond.len = dtemp;
	}

	return true;
}

bool GASP2struct::readAngle(tinyxml2::XMLElement *elem, string& errorstring, GASP2angle &ang, GASP2molecule mol) {
	int itemp;
	const char * stemp;
	string strtemp;
	vector<string> vstemp;
	double dtemp;


	//title
	//first search for an index, then parse as string
	if(! elem->QueryIntAttribute("title",&itemp) ) {
		if(itemp > 0)
			ang.label = itemp;
	} else {
		stemp = elem->Attribute("title");
		if(stemp)
			ang.label =  newName(stemp);
		else {
			errorstring = "An angle title name was not specified!\n";
			return false;
		}
	}

	//angle
	stemp = elem->Attribute("angle");
	if(stemp) {
		strtemp = stemp;
		vstemp = split(strtemp);
		int size = vstemp.size();
		if(size != 3) {
			errorstring = "The angle specification for "+names[ang.label]+" does not have 3 elements!\n";
			return false;
		}
		ang.a = atomLookup(newName(vstemp[0]), mol);
		ang.b = atomLookup(newName(vstemp[1]), mol);
		ang.c = atomLookup(newName(vstemp[2]), mol);
		if(ang.a < 0 || ang.b < 0 || ang.c < 0) {
			errorstring = "An atom designation in angle "+names[ang.label]+" did not match any atom in the molecule!\n";
			return false;
		}
	}

	//min
	ang.minAng = 0.0;
	if(!elem->QueryDoubleAttribute("min", &dtemp)) {
		ang.minAng = dtemp;
	}
	if(ang.minAng < 0.0 || ang.minAng > 180.0) {
		errorstring = "The minimum angle in angle "+names[ang.label]+" is out of bounds (0.0 < ang < 180.0).\n";
		return false;
	}

	//max
	ang.maxAng = 180.0;
	if(!elem->QueryDoubleAttribute("max", &dtemp)) {
		ang.maxAng = dtemp;
	}
	if(ang.maxAng < 0.0 || ang.maxAng > 180.0) {
		errorstring = "The maximum angle in angle "+names[ang.label]+" is out of bounds (0.0 < ang < 180.0).\n";
		return false;
	}


	//val
	ang.ang = 0.0;
	if(!elem->QueryDoubleAttribute("val", &dtemp)) {
		ang.ang = dtemp;
	}

	return true;
}

bool GASP2struct::readMol(tinyxml2::XMLElement *elem, string& errorstring, GASP2molecule &mol) {

	double dtemp;
	const char * stemp;
	string strtemp;
	vector<string> vstemp;
	int itemp;

	GASP2atom tempatom;
	GASP2bond tempbond;
	GASP2angle tempangle;
	GASP2dihedral tempdih;
	tinyxml2::XMLElement *item;

	//name
	if(! elem->QueryIntAttribute("name",&itemp) ) {
		if(itemp > 0)
			mol.label = itemp;
	} else {
		stemp = elem->Attribute("name");
		if(stemp)
			mol.label =  newName(stemp);
		else {
			errorstring = "A molecule name was not specified!\n";
			return false;
		}
	}

	//plindex
	if(! elem->QueryIntAttribute("plind",&itemp) ) {
		if(itemp > 0)
			mol.plindex = itemp;
	}

	//symm
	if(! elem->QueryIntAttribute("symm",&itemp) ) {
		if(itemp > 0)
			mol.symm = itemp;
	}

	//expectvol
	if(!elem->QueryDoubleAttribute("expectvol", &dtemp)) {
		mol.expectvol = dtemp;
	}
	else {
		errorstring = "No expected volume was specified for "+names[mol.label]+"!\n";
		return false;
	}
	if(mol.expectvol <= 0.0) {
		errorstring = "The expected volume for "+names[mol.label]+" is out of bounds (vol > 0.0).\n";
		return false;
	}

	//parse the atoms
	//MUST COME BEFORE DIHEDRALS, etc
	item = elem->FirstChildElement("atom");
	if(!item) {
		errorstring = "No atoms specified in molecule!\n";
		return false;
	}
	while(item) {
		if(! readAtom(item, errorstring, tempatom) )
			return false;
		//commit the atom to the list
		mol.atoms.push_back(tempatom);
		item = item->NextSiblingElement("atom");
	}

	//plane
	stemp = elem->Attribute("plane");
	if(stemp) {
		strtemp = stemp;
		vstemp = split(strtemp);
		int size = vstemp.size();
		if(size != 3) {
			errorstring = "The plane specification for "+names[mol.label]+" does not have 3 elements!\n";
			return false;
		}
		mol.p1 = atomLookup(newName(vstemp[0]), mol);
		mol.p2 = atomLookup(newName(vstemp[1]), mol);
		mol.p3 = atomLookup(newName(vstemp[2]), mol);
		if(mol.p1 < 0 || mol.p2 < 0 || mol.p3 < 0) {
			errorstring = "An atom designation in plane for "+names[mol.label]+" did not match any atom in the molecule!\n";
			return false;
		}
	}

	//parse the dihedrals
	item = elem->FirstChildElement("dihedral");
	while(item) {
		if(! readDih(item, errorstring, tempdih, mol) )
			return false;
		mol.dihedrals.push_back(tempdih);
		item = item->NextSiblingElement("dihedral");
		tempdih.clear();
	}

	//parse the bonds
	item = elem->FirstChildElement("bond");
	while(item) {
		if(! readBond(item, errorstring, tempbond, mol) )
			return false;
		mol.bonds.push_back(tempbond);
		item = item->NextSiblingElement("bond");
	}

	//parse the angles
	item = elem->FirstChildElement("angle");
	while(item) {
		if(! readAngle(item, errorstring, tempangle, mol) )
			return false;
		mol.angles.push_back(tempangle);
		item = item->NextSiblingElement("angle");
	}

	//rotation matrix
	item = elem->FirstChildElement("rot");
	if(item) {
		stemp = item->GetText();
		if(stemp) {
			strtemp = stemp;
			vstemp = split(strtemp);
			int size = vstemp.size();
			if(size != 9) {
				errorstring = "The rotation matrix for "+names[mol.label]+" does not have 9 elements!\n";
				return false;
			}
			for(int i = 0; i < size; i++) {
				istringstream ( vstemp[i] ) >> dtemp;
				mol.rot[i/3][i%3] = dtemp;
			}
		}
	}

	//position vector
	item = elem->FirstChildElement("pos");
	if(item) {
		stemp = item->GetText();
		if(stemp) {
			strtemp = stemp;
			vstemp = split(strtemp);
			int size = vstemp.size();
			if(size != 3) {
				errorstring = "The position vector for "+names[mol.label]+" does not have 3 elements!\n";
				return false;
			}
			for(int i = 0; i < size; i++) {
				istringstream ( vstemp[i] ) >> dtemp;
				mol.pos[i] = dtemp;
			}
		}
	}
	return true;
}

bool GASP2struct::readCell(tinyxml2::XMLElement *elem, string& errorstring, GASP2cell &cell) {

	double dtemp;
	int itemp;
	const char * stemp;
	string strtemp;

	//a
	cell.a = 0.01;
	if(!elem->QueryDoubleAttribute("a", &dtemp)) {
		cell.a = dtemp;
	}
	if(cell.a <= 0.0) {
		errorstring = "Cell length A must be greater than 0!\n";
		return false;
	}

	//b
	cell.b = 0.01;
	if(!elem->QueryDoubleAttribute("b", &dtemp)) {
		cell.b = dtemp;
	}
	if(cell.b <= 0.0) {
		errorstring = "Cell length B must be greater than 0!\n";
		return false;
	}


	//c
	cell.c = 0.01;
	if(!elem->QueryDoubleAttribute("c", &dtemp)) {
		cell.c = dtemp;
	}
	if(cell.c <= 0.0) {
		errorstring = "Cell length C must be greater than 0!\n";
		return false;
	}


	//alpha
	cell.alpha = 90.0;
	if(!elem->QueryDoubleAttribute("al", &dtemp)) {
		cell.alpha = dtemp;
	}
	if(cell.alpha < -180.0 || cell.alpha > 180.0) {
		errorstring = "Cell alpha is out of bound (-180.0 < x < 180.0)!\n";
		return false;
	}


	//beta
	cell.beta = 90.0;
	if(!elem->QueryDoubleAttribute("bt", &dtemp)) {
		cell.beta = dtemp;
	}
	if(cell.beta < -180.0 || cell.beta > 180.0) {
		errorstring = "Cell beta is out of bound (-180.0 < x < 180.0)!\n";
		return false;
	}


	//gamma
	cell.gamma = 90.0;
	if(!elem->QueryDoubleAttribute("gm", &dtemp)) {
		cell.gamma = dtemp;
	}
	if(cell.gamma < -180.0 || cell.gamma > 180.0) {
		errorstring = "Cell gamma is out of bound (-180.0 < x < 180.0)!\n";
		return false;
	}


	//rA
	cell.ratA = 0.01;
	if(!elem->QueryDoubleAttribute("rA", &dtemp)) {
		cell.ratA = dtemp;
	}
	if(cell.ratA <= 0.0) {
		errorstring = "Cell ratio A must be greater than 0!\n";
		return false;
	}


	//rB
	cell.ratB = 0.01;
	if(!elem->QueryDoubleAttribute("rB", &dtemp)) {
		cell.ratB = dtemp;
	}
	if(cell.ratB <= 0.0) {
		errorstring = "Cell ratio B must be greater than 0!\n";
		return false;
	}


	//rC
	cell.ratC = 0.01;
	if(!elem->QueryDoubleAttribute("rC", &dtemp)) {
		cell.ratC = dtemp;
	}
	if(cell.ratC <= 0.0) {
		errorstring = "Cell ratio C must be greater than 0!\n";
		return false;
	}


	//spacegroup
	cell.spacegroup = 1;
	if(!elem->QueryIntAttribute("spacegroup", &itemp)) {
		if(itemp > 0 && itemp <= 230)
			cell.spacegroup = itemp;
		else {
			errorstring = "A non-valid spacegroup identifier was given (out of bounds).\n";
			return false;
		}
	}
	else {
		stemp = elem->Attribute("spacegroup");
		strtemp = stemp;
		cell.spacegroup = 1;
		for(int i = 0; i < 230; i++) {
//			if(strtemp == ) {
//				cell.spacegroup = i;
//				break;
//			}
		}
		if(cell.spacegroup == 0) {
			errorstring = "A non-valid string identifier was given for a spacegroup value.\n";
			return false;
		}

	}


	//stoichiometry


	return true;
}

bool GASP2struct::readInfo(tinyxml2::XMLElement *elem, string& errorstring) {


	const char * stemp;
	string strtemp;
	UUID uuidtemp;
	double dtemp;
	int itemp;


	//ID
	stemp = elem->Attribute("id");
	if(stemp) {
		strtemp = stemp;
		if(strtemp.size() != 36) {
			errorstring = "A given UUID does not have the right length!\n";
			return false;
		}
		uuidtemp.frStr(strtemp);
		ID = uuidtemp;
	}

	//ParentA
	stemp = elem->Attribute("parentA");
	if(stemp) {
		strtemp = stemp;
		if(strtemp.size() != 36) {
			errorstring = "A given UUID does not have the right length!\n";
			return false;
		}
		uuidtemp.frStr(strtemp);
		parentA = uuidtemp;
	}

	//ParentB
	stemp = elem->Attribute("parentB");
	if(stemp) {
		strtemp = stemp;
		if(strtemp.size() != 36) {
			errorstring = "A given UUID does not have the right length!\n";
			return false;
		}
		uuidtemp.frStr(strtemp);
		parentB = uuidtemp;
	}

	//opt flag
	stemp = elem->Attribute("opt");
	if(stemp) {
		strtemp = stemp;
		if(strtemp != "true" && strtemp != "false" && strtemp != "t" && strtemp != "f")	{
			errorstring = "Opt flag is not true or false!\n";
			return false;
		}
		didOpt = false;
		if(strtemp == "true")
			didOpt = true;
	}

	//complete flag
	stemp = elem->Attribute("complete");
	if(stemp) {
		strtemp = stemp;
		if(strtemp != "true" && strtemp != "false" && strtemp != "t" && strtemp != "f")	{
			errorstring = "Complete flag is not true or false!\n";
			return false;
		}
		complete = false;
		if(strtemp == "true")
			complete = true;
	}

	//fitcell flag
	stemp = elem->Attribute("fitcell");
	if(stemp) {
		strtemp = stemp;
		if(strtemp != "true" && strtemp != "false" && strtemp != "t" && strtemp != "f")	{
			errorstring = "Fitcell flag is not true or false!\n";
			return false;
		}
		isFitcell = false;
		if(strtemp == "true")
			isFitcell = true;
	}

	//energy
	energy = 0.0;
	if(!elem->QueryDoubleAttribute("energy", &dtemp)) {
		energy = dtemp;
	}

	//force
	force = 0.0;
	if(!elem->QueryDoubleAttribute("force", &dtemp)) {
		force = dtemp;
	}

	//pressure
	pressure = 0.0;
	if(!elem->QueryDoubleAttribute("pressure", &dtemp)) {
		pressure = dtemp;
	}

	//error
	stemp = elem->Attribute("error");
	if(stemp) {
		strtemp = stemp;
		if(strtemp == "OKStruct")
			finalstate = OKStruct;
		else if(strtemp == "OptBadBond")
			finalstate = OptBadBond;
		else if(strtemp == "OptBadAng")
			finalstate = OptBadAng;
		else if(strtemp == "OptBadDih")
			finalstate = OptBadDih;
		else if(strtemp == "FitcellBadDih")
			finalstate = FitcellBadDih;
		else if(strtemp == "NoFitcell")
			finalstate = NoFitcell;
		else {
			errorstring = "An structure error was given, but it does not match any known error types!\n";
			return false;
		}
	}

	//time
	time = 0;
	if(!elem->QueryIntAttribute("time", &itemp)) {
		time = itemp;
	}

	//steps
	steps = 0;
	if(!elem->QueryIntAttribute("steps", &itemp)) {
		steps = itemp;
	}

	return true;

}

void GASP2struct::logStruct() {
	cout << endl << endl;
	cout << setprecision(8) << fixed;
	cout << "Structure id: " << ID.toStr() << endl;
	cout << "Parent A id : " << parentA.toStr() << endl;
	cout << "Parent B id : " << parentB.toStr() << endl << endl;

//	cout << "Name index:" << endl;
//	for(int i = 0; i< names.size(); i++) {
//		cout << "  " << i << ": " << names[i] << endl;
//	}
//	cout << endl;

	cout << "Number of molecules: " << molecules.size() << endl << endl;;

	for(int i = 0; i < molecules.size(); i++) {
		cout << "  Molecule " << names[molecules[i].label] << endl;
		for(int j = 0; j < molecules[i].atoms.size(); j++) {
			GASP2atom at = molecules[i].atoms[j];
			cout << "  At:" <<j<<" "<< names[at.label]<<" "<<getElemName(at.type)<<" "<<at.pos<<endl;
		}
		for(int j = 0; j < molecules[i].dihedrals.size(); j++) {
			GASP2dihedral dh = molecules[i].dihedrals[j];
			cout << "  Dih " << names[dh.label]<<" "<<dh.a<<" "<<dh.b<<" "<<dh.c<<" "<<dh.d<<" "<<endl;
		}
		for(int j = 0; j < molecules[i].angles.size(); j++) {
			GASP2angle an = molecules[i].angles[j];
			cout << "  Ang " << names[an.label]<<" "<<an.a<<" "<<an.b<<" "<<an.c<<endl;
		}
		for(int j = 0; j < molecules[i].bonds.size(); j++) {
			GASP2bond bn = molecules[i].bonds[j];
			cout << "  Bon " << names[bn.label]<<" "<<bn.a<<" "<<bn.b<<endl;
		}
		cout << endl;
	}
	cout << endl;

}



