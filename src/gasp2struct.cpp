#include "gasp2struct.hpp"


vector<string> GASP2struct::names;

GASP2struct::GASP2struct() {
	molecules.clear();
	builtcell.clear();
	finalcell.clear();
	doOpt = false;
	didFitcell = false;
	didOpt = false;
	complete = false;
	finalstate = OKStruct;

}

double GASP2struct::getVolume() {
	double phi = 2 * cos(unit.alpha) * cos(unit.beta) * cos(unit.gamma);
    return ( 2 * unit.a * unit.b * unit.c *
             sqrt(1 - cos(unit.alpha)*cos(unit.alpha)
            		- cos(unit.beta)*cos(unit.beta)
            		- cos(unit.gamma)*cos(unit.gamma)
            		+ phi) );

}

bool GASP2struct::fitcell(double distcell, double intradist) {


	return true;
}

bool GASP2struct::check() {
	Index a,b,c,d;
	double val;

	for(int i = 0; i < finalcell.size(); i++) {
		//bonds
		for(int j = 0; j < finalcell[i].bonds.size(); j++) {
			a = finalcell[i].bonds[j].a;
			b = finalcell[i].bonds[j].b;
			val = len(finalcell[i].atoms[a].pos - finalcell[i].atoms[b].pos);
			if(val > finalcell[i].bonds[j].maxLen || val < finalcell[i].bonds[j].minLen)
				finalstate = OptBadBond;
			finalcell[i].bonds[j].len = val;
		}
		//angles
		for(int j = 0; j < finalcell[i].angles.size(); j++) {
			a = finalcell[i].angles[j].a;
			b = finalcell[i].angles[j].b;
			c = finalcell[i].angles[j].c;
			val = angle(finalcell[i].atoms[a].pos,finalcell[i].atoms[b].pos,finalcell[i].atoms[c].pos);
			if(val > finalcell[i].angles[j].maxAng || val < finalcell[i].angles[j].minAng)
				finalstate = OptBadAng;
			finalcell[i].angles[j].ang = val;
		}
		//dihedrals
		for(int j = 0; j < finalcell[i].dihedrals.size(); j++) {
			a = finalcell[i].dihedrals[j].a;
			b = finalcell[i].dihedrals[j].b;
			c = finalcell[i].dihedrals[j].c;
			d = finalcell[i].dihedrals[j].d;
			val = dihedral(finalcell[i].atoms[a].pos,finalcell[i].atoms[b].pos,
					finalcell[i].atoms[c].pos,finalcell[i].atoms[d].pos);
			if(val > finalcell[i].dihedrals[j].maxAng || val < finalcell[i].dihedrals[j].minAng)
				finalstate = OptBadDih;
			finalcell[i].dihedrals[j].ang = val;
		}
		//plane rotations


		//position

	}
	return true;
}

bool GASP2struct::evaluate() {
	if(didFitcell) {
		finalcell = builtcell;
		finalunit = unit;
	}
	else {
		finalstate = NoFitcell;
		return false;
	}


	//never evaluate something that has already been evaluated
	if(doOpt & !didOpt) {
		complete = eval(finalcell, finalunit, energy);
		didOpt = true;
	}

	return complete;
};

Index GASP2struct::newName(string name) {
	Index n = names.size();
	names.push_back(name);
	return n;
}

Index GASP2struct::newName(const char * name) {
	Index n = names.size();
	string temp(name);
	names.push_back(temp);
	return n;
}


bool GASP2struct::parseXML(tinyxml2::XMLDocument *doc, string& errorstring ) {

	double dtemp;
	const * char stemp;
	string strtemp;
	int itemp;

	//parse the crystal
	//at this point we assume only one crystal in the specification
	//multiple spacegroups will likely be handled by a delimited list
	//shoudl that mode of operation ever be implemented
	tinyxml2::XMLElement *crystal = doc->FirstChildElement("mgac")->FirstChildElement("crystal");
	if(crystal) {
		//label
		stemp = crystal->Attribute("name");
		if(stemp)
			crylabel =  newName(stemp);
		else {
			errorstring = "A crystal name was not specified!\n";
			return false;
		}
		//spacegroup, minFrac, maxFrac





	//parse the molecules
		GASP2molecule tempmol;
		GASP2atom tempatom;
		GASP2bond tempbond;
		GASP2angle tempangle;
		GASP2dihedral tempdih;


		tinyxml2::XMLElement *molecule, *item;
		molecule = crystal->FirstChildElement("molecule");
		while(molecule) {
			//name
			stemp = molecule->Attribute("name");
			if(stemp)
				tempmol.label =  newName(stemp);
			else {
				errorstring = "A molecule name was not specified!\n";
				return false;
			}

			//thought: how to deal with variable number of molecules?

			//count
			if(!molecule->QueryIntAttribute("count", &itemp)) {
				tempmol.count = itemp;
				if(tempmol.count <= 0) {
					errorstring = "The number of molecules must be greater than 0!\n";
					return false;
				}
			}


			//parse the atoms
			//MUST COME BEFORE DIHEDRALS, etc

			item = molecule->FirstChildElement("atom");
			if(!item) {
				errorstring = "No atoms specified in molecule!\n";
				return false;
			}
			while(item) {

				//atom elem
				stemp = item->Attribute("elem");
				if(stemp)
					tempmol.type = getElemType(stemp);
				else {
					errorstring = "An atom element was not specified!\n";
					return false;
				}
				if(tempmol.type == UNK) {
					errorstring = "An unknown or unsupported element name was used!\n";
					return false;
				}

				//atom label
				stemp = item->Attribute("title");
				if(stemp)
					tempatom.label =  newName(stemp);
				else {
					errorstring = "An atom title name was not specified!\n";
					return false;
				}

				//x coord (fractional)
				if(!item->QueryDoubleAttribute("x", &dtemp)) {
					tempatom.pos[0] = dtemp;
				}
				else {
					errorstring = "No value was given for an X coordinate!\n";
					return false;
				}
				//y coord (fractional)
				if(!item->QueryDoubleAttribute("y", &dtemp)) {
					tempatom.pos[1] = dtemp;
				}
				else {
					errorstring = "No value was given for a Y coordinate!\n";
					return false;
				}
				//z coord (fractional)
				if(!item->QueryDoubleAttribute("z", &dtemp)) {
					tempatom.pos[2] = dtemp;
				}
				else {
					errorstring = "No value was given for a Z coordinate!\n";
					return false;
				}


				//commit the atom to the list
				tempmol.atoms.push_back(tempatom);
				item = molecule->NextSiblingElement("atom");
			}

			//parse the dihedrals
			item = molecule->FirstChildElement("dihedral");
			while(item) {


				tempmol.atoms.push_back(tempdih);
				item = molecule->NextSiblingElement("dihedral");
			}

			//parse the bonds
			item = molecule->FirstChildElement("bond");
			while(item) {


				tempmol.atoms.push_back(tempbond);
				item = molecule->NextSiblingElement("bond");
			}

			//parse the angles
			item = molecule->FirstChildElement("angle");
			while(item) {


				tempmol.atoms.push_back(tempangle);
				item = molecule->NextSiblingElement("angle");
			}

			//push the molecule(s) onto the list
			for(int i = 0; i < tempmol.count; i++)
				this->molecules.push_back(tempmol);
			molecule = molecule->NextSiblingElement("molecule")
		}



	}

	return true;
}

