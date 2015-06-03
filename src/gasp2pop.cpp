#include "gasp2pop.hpp"

using namespace std;

struct {
	bool operator() (GASP2struct a, GASP2struct b) {
		return a.getEnergy() < b.getEnergy();
	}
} ecomp;

struct {
	bool operator() (GASP2struct a, GASP2struct b) {
		double ascore = a.getVolScore();
		double bscore = b.getVolScore();
		if(ascore < 0.0)
			ascore = numeric_limits<double>::max();
		if(bscore < 0.0)
			bscore = numeric_limits<double>::max();
		return ascore < bscore;
	}
} vcomp;

struct {
	bool operator() (GASP2struct a, GASP2struct b) {
		return a.getSymmcount() > b.getSymmcount();
	}
} symmcomp;

//sorts by energy first
//if the structure does not have an energy (ie,
// energy = 0.0) then it is sorted by volume with
//all other no-energy structures
void GASP2pop::energysort() {
	std::sort(structures.begin(), structures.end(), ecomp);
	for(int i = 0; i < structures.size(); i++) {
		if(structures[i].getEnergy() >= 0.0) {
			std::sort(structures.begin() + i, structures.end(), vcomp);
			break;
		}
	}

//	for(int i = 0; i < size(); i++)
//		cout << structures[i].getEnergy() << endl;
}

void GASP2pop::volumesort() {
	std::sort(structures.begin(), structures.end(), vcomp);

//	for(int i = 0; i < size(); i++)
//		cout << structures[i].getVolScore() << endl;
}

//sorts by number of symmetry operations, with
//larger numbers coming first
void GASP2pop::symmsort() {
	std::sort(structures.begin(), structures.end(), symmcomp);

}

GASP2pop GASP2pop::subpop(int start, int subsize) {
	GASP2pop out;
	int end = start + subsize;
	if( end > size() )
		end = size();
	for(int i = start; i < end; i++) {
		out.structures.push_back(structures[i]);
	}
	return out;
}

void GASP2pop::init(GASP2struct s, int size, Spacemode mode, Index spcg) {
	for(int i = 0; i < size; i++) {
		structures.push_back(s);
		//try ten times to get a good init
		for(int j = 0; j < 10; j++)
			if(structures.back().init(mode, spcg)) break;
	}
	scale(1.0,0.0,0.0);

}


//newPop assumes that the scaling is valid if the size of the
//current population size matches the scaling size
//if the sizes are different, it applies a constant scaling
//size the ranking of the population is unknown
GASP2pop GASP2pop::newPop(int size, Spacemode smode, GAselection mode) {
	GASP2pop out;
	GASP2struct a,b;
	out.structures.reserve(size);
	int selA, selB; //selection indices
	if(mode == Roulette) {
		if(scaling.size() != structures.size()) {
			cout << mark() << "WARNING: Scaling automatically applied on population..." << endl;
			scale(1.0,0.0,0.0);
		}
		discrete_distribution<int> d(scaling.begin(), scaling.end());
		for(int i = 0; i < size; i++) {
			selA = d(rgen);
			selB = d(rgen);
			while(selA == selB)
				selB = d(rgen);
			structures[selA].crossStruct(structures[selB], a,b, 0.5, smode);
			out.structures.push_back(a);
			//cout << "selA/selB: " << selA << "/" << selB << endl;
		}
	}
	else { //pattern
		cout << "Pattern mode not implemented!" << endl;
	}
	return out;
}

GASP2pop GASP2pop::fullCross(Spacemode mode) {
	GASP2pop out;
	GASP2struct a,b;
	int size = structures.size();
	out.structures.reserve(size*size - size);

	for(int i = 0; i < size; i++) {
		for(int j = i; j < size; j++) {
			if(i==j) continue;
			structures[i].crossStruct(structures[j],a,b, 1.0, mode);
			out.structures.push_back(a);
			out.structures.push_back(b);
		}
	}
	return out;
}

void GASP2pop::addIndv(int add) {

	int init = structures.size();
	structures.reserve(init + add);
	//this kind of implies the pop has at least one struct
	for(int i = 0; i < add; i++)
		structures.push_back(structures[0]);
	for(int i = init; i < (init+add); i++)
		structures[i].init();

}

void GASP2pop::addIndv(GASP2pop add) {
	structures.reserve(structures.size() + add.structures.size());
	structures.insert(structures.end(), add.structures.begin(), add.structures.end());

}

void GASP2pop::mergeIndv(GASP2pop add, int ind) {
	if(add.size() > 1)
		return;
	//structures.reserve(structures.size() + add.structures.size());
	//structures.insert(structures.begin()+ind, add.structures.begin(), add.structures.end());
	structures[ind] = add.structures[0];

}

GASP2pop GASP2pop::remIndv(int n) {
	GASP2pop bad;
	for(int i = 0; i < n; i++) {
		bad.structures.push_back(structures.back());
		structures.pop_back();
	}
	return bad;
}


GASP2pop GASP2pop::volLimit(GASP2pop &bad) {
	GASP2pop ok;

	for(int i = 0; i < size(); i++) {
		if(structures[i].minmaxVol() && structures[i].fitcelled() && !structures[i].rejected() )
			ok.structures.push_back(structures[i]);
		else
			bad.structures.push_back(structures[i]);
	}

	return ok;
}

GASP2pop GASP2pop::symmLimit(GASP2pop &bad, int limit) {
	GASP2pop ok;

	for(int i = 0; i < size(); i++) {
		if(structures[i].getSymmcount() < limit)
			ok.structures.push_back(structures[i]);
		else
			bad.structures.push_back(structures[i]);
	}

	return ok;
}

//scaling must be manually after ranking

//scaling obeys a three part weighted equation:
// scale(score) = c + l*lin_score + e*exp_score;
//c is a constant base value
//l and e are weights
//lin_score is the normalized difference with the minimum score:
// 1 - [(score - min)/(max/min)]
//exp_score multiplies lin_score by 10 and passes it through
//and exponential filter: e^(-0.5*10.0*lin_score)
//the volume and exponential components are added together
//since they are both normalized on the range of values
//present for both volume and energy. structures with
//both energy and volume will always be better than
void GASP2pop::scale(double con, double lin, double exp) {
	scaling.clear();

	if(structures.size() == 1) {
		scaling.push_back(1.0);
		return;
	}

	if(lin == 0.0 && exp == 0.0) {
		for(int i = 0; i < size(); i++)
			scaling.push_back(1.0);
		return;
	}


	double val, diffE, diffV;
	vector<double> vol, ener;
	double minE = numeric_limits<double>::max(), maxE = -1.0 * numeric_limits<double>::max();
	double minV = numeric_limits<double>::max(), maxV = -1.0 * numeric_limits<double>::max();

	for(int i = 0; i < structures.size(); i++) {
		vol.push_back(structures[i].getVolScore());
		//cout << "vol " << vol.back() << endl;
		ener.push_back(structures[i].getEnergy());

		if(vol.back() >= 0.0) {
			if(vol.back() > maxV) maxV = vol.back();
			if(vol.back() < minV) minV = vol.back();
		}
		//we assume energies will usually be negative, since E = 0.0 is more or less impossible
		//since some structures to no have energies calculated, we do not want to influence
		//the energy
		if(ener.back() < 0.0) {
			if(ener.back() > maxE) maxE = ener.back();
			if(ener.back() < minE) minE = ener.back();
		}
	}

	//cout << "minV/maxV " << minV << "/" << maxV << endl;
	//cout << "minE/maxE " << minE << "/" << maxE << endl;

	double volscore, enerscore;
	diffE = maxE-minE;
	diffV = maxV-minV;
	for(int i = 0; i < structures.size(); i++) {

		val = con; //assign a constant minimum to all values
		if(vol[i] >= 0.0) {
			volscore = (vol[i] - minV)/diffV;
			val += lin*(1.0-volscore); //add the linear component
			val += exp*(std::exp(-0.5*10.0*volscore));//add the exponential component
			//cout << "vol " << volscore;
		}
		if(ener[i] < 0.0) { //add energy components if the energy is valid
			enerscore = (ener[i] - minE)/diffE;
			//cout << " enerscore " << enerscore;
			val += lin*(1.0-enerscore) + exp*std::exp(-0.5*10.0*enerscore);
		}

		scaling.push_back(val);
		//cout << " scale " << scaling[i] << endl;
	}
}

//this determines whether a structure is mutated
//if a mutation happens, then 0.15 of all genes
//within a structure will be mutated, on average
//whether or not this is a slient mutation is irrrelevant
//the point is that the structure is modified
//enough that it still partially resembles the initial structure
void GASP2pop::mutate(double rate, Spacemode mode) {
	//place boundaries on rate
	if(rate >= 1.0)
		rate = 1.0;
	if(rate < 0.0)
		rate = 0.0;

	std::bernoulli_distribution p(rate);
	for(int i = 0; i < structures.size(); i++) {
		if(p(rgen))
			structures[i].mutateStruct(0.15, mode);
	}

}

string GASP2pop::saveXML() {
	string open("<pop>\n");
	string close("</pop>\n");
	string out("");
	for(int i = 0; i < structures.size(); i++) {
		out += structures[i].serializeXML();
	}
	return open + out + close;
}

bool GASP2pop::parseXML(string name, string &errorstring) {

	tinyxml2::XMLDocument doc;
	doc.Parse(name.c_str());




	if(doc.ErrorID() == 0) {
		tinyxml2::XMLElement *pop = doc.FirstChildElement("pop");
		if(this->loadXMLrestart(pop,errorstring))
			return false;
	}

	return true;
}

bool GASP2pop::loadXMLrestart(tinyxml2::XMLElement *pop, string & errorstring) {
	GASP2struct temp;
	if(pop) {
		tinyxml2::XMLElement * crystal = pop -> FirstChildElement("crystal");
		if(!crystal) {
			errorstring = "Something went wrong when parsing the population! There may be no structures!\n";
			return false;
		}
		//throw away anything in the structures.
		structures.clear();
		while(crystal) {
			if(!temp.parseXMLStruct(crystal, errorstring))
				return false;
			//push the molecule(s) onto the list
			this->structures.push_back(temp);
			temp.clear();
			crystal = crystal->NextSiblingElement("crystal");
		}
	}
	else
		return false;

	return true;
}

bool GASP2pop::writeCIF(string name) {

	ofstream outf;
	outf.open(name.c_str(), ofstream::out | ofstream::app);
	if(outf.fail())
		return false;

	string temp;

	for(int i = 0; i < size(); i++) {
		temp = "";
		if(structures[i].cifString(temp))
			outf << temp << endl;
	}

	outf.close();
	return true;
}

void GASP2pop::runFitcell(int threads) {
	//GASP2pop out; out.structures = this->structures;

	if (size() < threads)
		threads = size();
	//setup the futures
	vector<future<bool>> futures(threads);
	chrono::milliseconds timeout(0);
	chrono::milliseconds thread_wait(20);

	//cout << "size: " << size() << endl;

	int thread_run = 0;
	//for all the structures
	for(int i = 0; i < size(); ) {

		//launch
		for(int j = 0; j < threads; j++) {
			if(!futures[j].valid() && (i < size()) && thread_run < threads ) {
				//cout << "instance i: " << i << endl;
				futures[j] = async(launch::async, &GASP2struct::fitcell, &structures[i], 300.0);
				i++;
				thread_run++;
			}
		}
		//cleanup
		for(int j = 0; j < threads; j++) {
			if(futures[j].valid() && futures[j].wait_for(timeout)==future_status::ready){
				futures[j].get();
				thread_run--;
			}
		}
		//wait so we don't burn cycles
		this_thread::sleep_for(thread_wait);

	}

	for(int j = 0; j < threads; j++) {
		if(futures[j].valid()) {
			//cout << "waiting on j: " << j << endl;
			futures[j].wait();
			futures[j].get();
			//cout << "got j: " << j << endl;
		}
	}


	//cout << structures[0].serializeXML() << endl;
	cout << mark() << "where am I?" << endl;
	//return out;
}




void GASP2pop::runEval(string hosts, GASP2param p, bool (*eval)(vector<GASP2molecule>&, GASP2cell&, double&, double&, double&, time_t&, string, GASP2param)) {

	future<bool> thread;
	bool result;

	chrono::milliseconds timeout(0);
	chrono::milliseconds thread_wait(100);

	for(int i = 0; i < size(); ) {

		//launch
		if(!thread.valid()) {
			if(!structures[i].completed()) {
				structures[i].setEval(eval);
				thread = async(launch::async, &GASP2struct::evaluate, &this->structures[i], hosts, p);
				i++;
			}
		}
		//cleanup
		if(thread.wait_for(timeout)==future_status::ready)
			thread.get();

		//wait so we don't burn cycles
		this_thread::sleep_for(thread_wait);

	}
	if(thread.valid()) {
		thread.wait();
		thread.get();
	}


//	cout << mark() << "pop energies" << endl;
//	for(int i = 0; i < size(); i++)
//		cout << "energy: " << structures[i].getEnergy() << endl;


}


