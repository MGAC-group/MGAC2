#include "gasp2pop.hpp"

using namespace std;

//sorting functions for different structure sorting mthods

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
		return a.getContacts() > b.getContacts();
	}
} ccomp;

struct {
	bool operator() (GASP2struct a, GASP2struct b) {
		return a.getPseudoenergy() < b.getPseudoenergy();
	}
} psecomp;


struct {
	bool operator() (GASP2struct a, GASP2struct b) {
		return a.getSymmcount() > b.getSymmcount();
	}
} symmcomp;

struct {
	bool operator() (GASP2pop a, GASP2pop b) {
		double ae = 1.0, be = 1.0;
		if(a.size() >= 1)
			ae = a.indv(0)->getEnergy();
		if(b.size() >= 1)
			be = b.indv(0)->getEnergy();
		return ae < be;
	}
} popcomp;



//sorts by energy first
//if the structure does not have an energy (ie,
// energy = 0.0) then it is sorted by volume with
//all other no-energy structures
void GASP2pop::energysort() {
	std::sort(structures.begin(), structures.end(), ecomp);
	for(int i = 0; i < structures.size(); i++) {
		if(structures[i].getEnergy() >= 0.0) {
			std::sort(structures.begin() + i, structures.end(), psecomp);
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


//extract a subpopulation of the given population
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

//randomly initialize the population
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
	if(structures.size() < 2) {
		return out;
	}

	//cout << "newpop" << endl;
	out.structures.reserve(size*2);
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
			if(selA == selB)
				selB = d(rgen);
			structures[selA].crossStruct(structures[selB], a,b, 0.6, smode);
			out.structures.push_back(a);
			out.structures.push_back(b);
			//cout << "selA/selB: " << selA << "/" << selB << endl;
		}
	}
	else { //pattern
		cout << "Pattern mode not implemented!" << endl;
	}
	return out;
}


//generate a new population from two different populations
GASP2pop GASP2pop::newPop(GASP2pop alt, int size, Spacemode smode, GAselection mode) {
	GASP2pop out;
	GASP2struct a,b;
	//cout << "alt newpop" << endl;
	if(structures.size() < 2 || alt.structures.size() < 2) {
		return out;
	}

	out.structures.reserve(2*size);
	int selA, selB; //selection indices
	if(mode == Roulette) {
		if(scaling.size() != structures.size()) {
			cout << mark() << "WARNING: Scaling automatically applied on population..." << endl;
			scale(1.0,0.0,0.0);
		}
		if(alt.scaling.size() != alt.structures.size()) {
			cout << mark() << "WARNING: Scaling automatically applied on population..." << endl;
			alt.scale(1.0,0.0,0.0);
		}
//		scaling.reserve(structures.size());
//		for(int n = 0; n < structures.size(); n++)
//			scaling.push_back(1.0);
//		alt.scaling.reserve(alt.structures.size());
//		for(int n = 0; n < alt.structures.size(); n++)
//			alt.scaling.push_back(1.0);


		discrete_distribution<int> d(scaling.begin(), scaling.end());
		discrete_distribution<int> c(alt.scaling.begin(), alt.scaling.end());
		for(int i = 0; i < size; i++) {
			selA = d(rgen);
			selB = c(rgen);
//			while(selA == selB)
//				selB = d(rgen);
			structures[selA].crossStruct(alt.structures[selB], a,b, 0.6, smode);
			out.structures.push_back(a);
			out.structures.push_back(b);
			//cout << "selA/selB: " << selA << "/" << selB << endl;
		}

	}
	else { //pattern
		//AML: I had plans for a pattern method but never finished implementing it
		cout << "Pattern mode not implemented!" << endl;
	}
	return out;
}

//cross all structures with every other structure in the population
GASP2pop GASP2pop::fullCross(Spacemode mode) {
	GASP2pop out;
	GASP2struct a,b;
	int size = structures.size();
	out.structures.reserve(size*size - size);

	for(int i = 0; i < size; i++) {
		for(int j = i; j < size; j++) {
			if(i==j) continue;
			structures[i].crossStruct(structures[j],a,b, 0.5, mode);
			out.structures.push_back(a);
			out.structures.push_back(b);
		}
	}
	return out;
}

//alternative fullcross for two pops;
//cross every structure in one population with every structure in the other
GASP2pop GASP2pop::fullCross(Spacemode mode, GASP2pop alt) {
	GASP2pop out;
	GASP2struct a,b;
	int size = structures.size();
	int altsize = alt.size();
	out.structures.reserve(size*altsize);

	for(int i = 0; i < size; i++) {
		for(int j = 0; j < altsize; j++) {
			structures[i].crossStruct(alt.structures[j],a,b, 0.5, mode);
			out.structures.push_back(a);
			out.structures.push_back(b);
		}
	}
	return out;
}

//shuffle structures in place
//used for a more classic GA where genes are preserved in the population
GASP2pop GASP2pop::inplaceCross(Spacemode mode) {
	GASP2pop out;
	GASP2struct a,b;
	int size = structures.size();
	out.structures.reserve(size);

	std::shuffle(structures.begin(), structures.end(), rgen);

	for(int i = 0; (i < size) && ((i+1) < size); i+=2) {

			structures[i].crossStruct(structures[i+1],a,b, 0.2, mode);
			out.structures.push_back(a);
			out.structures.push_back(b);
	}
	return out;
}

//add random structures to the population
void GASP2pop::addIndv(int add) {

	int init = structures.size();
	structures.reserve(init + add);
	//this kind of implies the pop has at least one struct
	for(int i = 0; i < add; i++)
		structures.push_back(structures[0]);
	for(int i = init; i < (init+add); i++)
		structures[i].init();

}

//combine populations
void GASP2pop::addIndv(GASP2pop add) {
	//structures.reserve(structures.size() + add.structures.size());
	structures.insert(structures.end(), add.structures.begin(), add.structures.end());

}

//used to merge structures in GASP2control server program
void GASP2pop::mergeIndv(GASP2pop add, int ind) {
	if(add.size() > 1)
		return;
	//structures.reserve(structures.size() + add.structures.size());
	//structures.insert(structures.begin()+ind, add.structures.begin(), add.structures.end());
	structures[ind] = add.structures[0];

}

//checks to see is a structure is complete and then sets the opt to false
//AML: pretty sure this is a defunct function
GASP2pop GASP2pop::completeCheck() {
	GASP2pop output;
	for(int i = 0; i < size(); i++) {
		//if( structures[i].opted() && !(structures[i].completed()) ) {
			structures[i].setUnopted();
			output.addIndv(structures[i]);
		//}
	}

	return output;

}

//removes the last N individuals from the population
GASP2pop GASP2pop::remIndv(int n) {
	GASP2pop bad;
	for(int i = 0; i < n; i++) {
		bad.structures.push_back(structures.back());
		structures.pop_back();
	}
	return bad;
}

//splits structures depending if they are negative or positive
GASP2pop GASP2pop::energysplit(GASP2pop &zero) {
	GASP2pop ok;

	for(int i = 0; i < size(); i++) {
		if(structures[i].getEnergy() < 0.0 )
			ok.structures.push_back(structures[i]);
		else
			zero.structures.push_back(structures[i]);
	}

	return ok;
}

//filters structures on whether or not they pass the structure volume limits
//volume limits are established from minmaxVol()
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


//filters structures based on symmetry operation in the spacegroup
GASP2pop GASP2pop::symmLimit(GASP2pop &bad, int limit) {
	GASP2pop ok;

	for(int i = 0; i < size(); i++) {
		if(structures[i].getSymmcount() <= limit)
			ok.structures.push_back(structures[i]);
		else
			bad.structures.push_back(structures[i]);
	}

	return ok;
}

//this function sorts structures into spacegroup bins and then
//determines the unique set of structures per the genetic unique determination method
//in AML dissertation, Ch. 6
GASP2pop GASP2pop::spacebinUniques(int threads, vector<GASP2pop> clusterbins, GASP2param p, int binsave) {

	vector<GASP2pop> tempbin(230);
	GASP2pop out;
	//sort into bins

	for(int i = 0; i < size(); i++) {
		int group = structures[i].getSpace();
		tempbin[group-1].addIndv(structures[i]);
	}

	cout << mark() << "spacebinUnique, sorting: " << size() << endl;
	int total = 0;
	for(int i = 0; i < 230; i++) {
		//bins[i].energysort();
		//bins[i].dedup(300); //FIXME: hardcoded value of 300 dedups, may not be safe
		int s = tempbin[i].size();
		if(s <= 0) continue;
		int pre = out.size();
		out.addIndv(tempbin[i].getUniques(clusterbins[i], p, threads));
		total += (out.size() - pre);
		if( s > 0 )
			cout  << mark() << "spcg " << i+1 << ", size " << s  << ", uniq " << (out.size() - pre) << endl;
		//bins[i].stripClusters(clusterbins[i].size(),25);
		//cout << "post-strip" << endl;
		//bins[i].remIndv(bins[i].size() - 50);
		//clusterbins[i].assignClusterGroups(p);
	}

	cout << mark() << "total uniques: " << total << endl;

	return out;
}

//this is the same as spacebinUniques, but uses the clustering algorithm instead
void GASP2pop::spacebinCluster(int threads, vector<GASP2pop> &bins, vector<GASP2pop> &clusterbins, GASP2param p, int binsave) {

	vector<GASP2pop> tempbin(230);
	//reorder the bins in correct order
	for(int i = 0; i < bins.size(); i++) {
		if(bins[i].size() >= 1) {
			int index = bins[i].indv(0)->getSpace() - 1;
			tempbin[index] = bins[i];
		}
	}
	bins = tempbin;
	tempbin.clear();

	cout << "post-sbcsort" << endl;
	for(int i = 0; i < size(); i++) {
		int group = structures[i].getSpace();
		bins[group-1].addIndv(structures[i]);
	}
	cout << "post-add" << endl;

	for(int i = 0; i < 230; i++) {
		//bins[i].energysort();
		//bins[i].dedup(300); //FIXME: hardcoded value of 300 dedups, may not be safe
		bins[i].cluster(clusterbins[i], p, threads);
		if(clusterbins[i].size() > 0)
			cout  << mark() << "post-cluster " << i+1  << ":" << clusterbins[i].size() << endl;
		//bins[i].stripClusters(clusterbins[i].size(),25);
		//cout << "post-strip" << endl;
		//bins[i].remIndv(bins[i].size() - 50);
		//clusterbins[i].assignClusterGroups(p);
	}
	cout << mark() << "end-it" << endl;
	//std::sort(bins.begin(), bins.end(), popcomp);

}


//semi complicated function
//sorts all the structures by spacegroup,
//then removes all except the best binsize structures
//for each spacegroup, then combines them into a single
//pop for reuse
void GASP2pop::spacebinV(vector<GASP2pop> &bins, int binsize, int binsave) {

	vector<GASP2pop> tempbin(230);
	//reorder the bins in correct order
	for(int i = 0; i < bins.size(); i++) {
		if(bins[i].size() >= 1) {
			int index = bins[i].indv(0)->getSpace() - 1;
			tempbin[index] = bins[i];
		}
	}
	bins = tempbin;
	tempbin.clear();

	for(int i = 0; i < size(); i++) {
		int group = structures[i].getSpace();
		bins[group-1].addIndv(structures[i]);
	}

	for(int i = 0; i < 230; i++) {
		bins[i].energysort();
		//bins[i].dedup(binsize);
	}

	std::sort(bins.begin(), bins.end(), popcomp);

	for(int i = 0; i < binsave; i++) {
		if(bins[i].size() > binsize)
			bins[i].remIndv(bins[i].size() - binsize);
	}


}

//the original bin sorting function
GASP2pop GASP2pop::spacebin(int binsize, int binsave) {
	vector<GASP2pop> bins(230);
	GASP2pop out;


	for(int i = 0; i < size(); i++) {
		int group = structures[i].getSpace();
		bins[group-1].addIndv(structures[i]);
	}

	for(int i = 0; i < 230; i++) {
		bins[i].energysort();
		//bins[i].dedup(binsize);
	}

	std::sort(bins.begin(), bins.end(), popcomp);

	for(int i = 0; i < binsave; i++) {
		if(bins[i].size() > binsize)
			bins[i].remIndv(bins[i].size() - binsize);
		out.addIndv(bins[i]);
	}

	return out;

}

//scaling must be manually performedafter ranking

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
//whether or not this is a silent mutation is irrelevant
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


//save this population as XML
string GASP2pop::saveXML() {
	string open("<pop>\n");
	string close("</pop>\n");
	string out("");
	for(int i = 0; i < structures.size(); i++) {
		out += structures[i].serializeXML();
	}
	return open + out + close;
}

//parse an XML definition of the population
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


//loads an xml restart function
//although the xml restart method is defunct, this function is still used
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

//write a CIF to string
bool GASP2pop::writeCIF(string name) {

	ofstream outf;
	outf.open(name.c_str(), ofstream::out | ofstream::app);
	if(outf.fail())
		return false;

	string temp;

	for(int i = 0; i < size(); i++) {
		temp = "";
		if(structures[i].cifString(temp, i))
			outf << temp << endl;
	}

	outf.close();
	return true;
}

//run a multi-threaded fitcell on all structures in the population
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
	//cout << mark() << "where am I?" << endl;
	//return out;
}

//helper function for runSymmetrize
bool multiSymm(int start, int end, GASP2pop* self) {
	for(int index = start; index < end; index ++) {
		self->indv(index)->simplesymm();
	}
	return true;
}

//applies symmetry operation to structures after fitcell has
//already been performed in multi-threaded mode
//used especially for restarts and .cif conversion
void GASP2pop::runSymmetrize(int threads) {
	//GASP2pop out; out.structures = this->structures;

	int thread_run = 0;
	//for all the structures
	if (structures.size() < threads)
		threads = structures.size();
	if(threads == 0) {
		//cout << "zero thread" << endl;
		return;
	}

	if(size() < 2)  {
		structures[0].setClustergroup(0);
		return;
	}



	//setup the futures
	vector<future<bool>> futures(threads);
	vector<int> finalindex(threads + 1);
	chrono::milliseconds timeout(0);
	chrono::milliseconds thread_wait(20);

	//cout << "got there" << endl;
//		for(int i = 0; i < structures.size(); ) {
//			//get the distance to the current centroid
	int s,e;
	//cout << "size " << size() << endl;
	int threadavg = (size() / threads);
	//cout << "threadavg " << threadavg << endl;

	for(int i = 0; i < threads; i++)
		finalindex[i] = i*threadavg;
	finalindex[threads] = size();

	//for(int k = 0; k < finalindex.size(); k++)
	//	cout << "findex " << finalindex[k] << endl;

	for(int j = 0; j < threads; j++) {
		if(!futures[j].valid() ) {//&& (i < size()) && thread_run < threads ) {
			s = finalindex[j];
			e = finalindex[j+1];
			futures[j] = async(launch::async, multiSymm, s,e, this);
		}
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
	//cout << mark() << "where am I?" << endl;
	//return out;
}



//this function runs the external evaluation function on all threads
//this function is kind of broken though
//only QE really uses it, and when QE is run, only single structures are passed remotely.
void GASP2pop::runEval(string hosts, GASP2param p, bool (*eval)(vector<GASP2molecule>&, GASP2cell&, double&, double&, double&, time_t&, string, GASP2param)) {

//	future<bool> thread;
//	bool result;
//
//	chrono::milliseconds timeout(0);
//	chrono::milliseconds thread_wait(100);
//
	cout << "runeval inner size " << size() << endl;
	for(int i = 0; i < size(); i++) {
		if(!structures[i].completed()) {
			structures[i].setEval(eval);
			structures[i].evaluate(hosts, p);
		}
	}
//
//		//launch
//		if(!thread.valid()) {
//			if(!structures[i].completed()) {
//				structures[i].setEval(eval);
//				thread = async(launch::async, &GASP2struct::evaluate, &this->structures[i], hosts, p);
//				i++;
//			}
//		}
//		//cleanup
//		if(thread.wait_for(timeout)==future_status::ready)
//			thread.get();
//
//		//wait so we don't burn cycles
//		this_thread::sleep_for(thread_wait);
//
//	}
//	if(thread.valid()) {
//		thread.wait();
//		thread.get();
//	}




//	cout << mark() << "pop energies" << endl;
//	for(int i = 0; i < size(); i++)
//		cout << "energy: " << structures[i].getEnergy() << endl;


}

//VERY IMPORTANT! Before deduping population must be energy sorted.


//extract all population members of a particular cluster
GASP2pop GASP2pop::getCluster(int c) {
	GASP2pop out;
	for(int i = 0; i < size(); i++) {
		if(structures[i].getCluster() == c)
			out.addIndv(structures[i]);
	}
	return out;
}

//removes all but the last n structures from the cluster
void GASP2pop::stripClusters(int clusters, int n) {
	if(clusters == 0)
		return;
	vector<GASP2pop> groups(clusters);
	for(int i = 0; i < structures.size(); i++) {
		groups[structures[i].getCluster()].addIndv(structures[i]);
	}
	structures.clear();
	for(int i = 0; i < groups.size(); i++) {
		if(groups[i].size() > n)
			groups[i].structures.erase(groups[i].structures.begin(), groups[i].structures.end() - n);
		this->addIndv(groups[i]);
	}

}

//DO NOT USE THIS FUNCTION
//helper for dedup
bool dedupCompare(int start, int end, GASP2pop *self, GASP2param p) {
	double average, chebyshev, euclid, localchebyshev, localavg, num;
	for(int j = start; j < end; j++) {

		for(int i = 0; i < self->size(); i++) {
			//if(self->indv(i)->getCluster() > -1) continue;
			if(self->indv(i)->simpleCompare(*self->indv(j),p,average,chebyshev,euclid,num)) {
				self->indv(j)->setCluster(1);
				//break;
			}
		}
	}

}

//Used to remove duplicate structures
//is very broken and highly experimental
GASP2pop GASP2pop::dedup(GASP2param p, int threads) {
	double average, chebyshev, euclid, localchebyshev, localavg, num;
	bool result;

	GASP2pop out;
	p.clusterdiff = 0.90;
	clusterReset();


	for(int i = 0; i < size(); i++) {
		if(structures[i].getCluster() > -1) continue;
		for(int j = i+1; j < size(); j++) {
			if(structures[i].simpleCompare(structures[j],p,average,chebyshev,euclid,num)) {
				structures[j].setCluster(1);
			}
		}
	}

//	{
//		if (structures.size() < threads)
//			threads = structures.size();
//		if(threads == 0) {
//			return out;
//		}
//
//		vector<future<bool>> futures(threads);
//		vector<int> finalindex(threads + 1);
//		chrono::milliseconds timeout(0);
//		chrono::milliseconds thread_wait(20);
//
//		int thread_run = 0;
//
//		int s,e;
//		int totalload = ( size() * size() - size() ) / 2;
//		//cout << "size " << size() << endl;
//		int threadavg = (totalload / threads);
//		//cout << "threadavg " << threadavg << endl;
//
//		int val = 0;
//		finalindex[0] = 1;
//		int ind = 1;
//		for(int i = 0; i < size(); i++) {
//			val += (i+1);
//			if(val > threadavg) {
//				finalindex[ind] = i;
//				//cout << "findex " << finalindex[ind] << endl;
//				ind++;
//				val = 0;
//			}
//		}
//		finalindex[threads] = size();
//
//		for(int j = 0; j < threads; j++) {
//			if(!futures[j].valid() ) {//&& (i < size()) && thread_run < threads ) {
//				s = finalindex[j];
//				e = finalindex[j+1];
//				futures[j] = async(launch::async, dedupCompare, s,e, this, p);
//
//			}
//		}
//
//
////		}
//		for(int j = 0; j < threads; j++) {
//			if(futures[j].valid()) {
//				//cout << "waiting on j: " << j << endl;
//				futures[j].wait();
//				futures[j].get();
//				//cout << "got j: " << j << endl;
//			}
//		}
//	}



	int count = 0;
	for(int i = 0; i < size(); i++) {
		if(structures[i].getCluster() < 0)
			count++;
	}
	out.reserve(count);

	for(int i = 0; i < size(); i++) {
		if(structures[i].getCluster() < 0)
			out.addIndv(structures[i]);
	}

	return out;

}


//helper for getuniques
bool multiCompare(int start, int end, GASP2pop *self, GASP2pop *clusters, GASP2param p) {
	double average, chebyshev, euclid;
	double num;

	for(int index = start; index < end; index++) {
		//cout << "index " << index << endl;
		if( self->indv(index)->getCluster() > -1) continue;

		//GASP2struct s = *self.indv(index);
		int cluster = -1;


		//this pass is to find all structures that match an existing known cluster
		for(int j = 0; j < clusters->size(); j++) {
			if(self->indv(index)->simpleCompare(*clusters->indv(j), p, average, chebyshev, euclid, num)) {
				cluster = j;
				break;
			}
		}

		//cout << "cluster " << cluster << endl;

		if(cluster > -1) {
			self->indv(index)->setCluster(cluster);
			continue;
		}
		//cout << "first step " << index << endl;

		//this pass is to catch multiple simlar structures that are new
		for(int j = 0; j < index; j++) {
			//cout << "j" << endl;
			if(j == index) continue;
			if(self->indv(index)->simpleCompare(*self->indv(j), p, average, chebyshev, euclid, num)) {
				//cout << "2styles j " << j << " ind " << index << endl;
				self->indv(index)->setCluster(-2);
			}
			break;
		}

		//cout << "second step " << index << endl;


	}

	return true;
}



//get the unique structures
GASP2pop GASP2pop::getUniques(GASP2pop clusters, GASP2param p, int threads) {
	double average, chebyshev, euclid, localchebyshev, localavg;
	int res;

	GASP2pop out;

	//reset clusters
	//clusterReset();

	for(int n = 0; n < 2; n++) {
		if (structures.size() < threads)
			threads = structures.size();
		if(threads == 0) {
			return out;
		}

		vector<future<bool>> futures(threads);
		vector<int> finalindex(threads + 1);
		chrono::milliseconds timeout(0);
		chrono::milliseconds thread_wait(20);

		int thread_run = 0;

		int s,e;
		int totalload = ( size() * size() - size() ) / 2;
		//cout << "size " << size() << endl;
		int threadavg = (totalload / threads);
		//cout << "threadavg " << threadavg << endl;

		int val = 0;
		finalindex[0] = 0;
		int ind = 1;
		for(int i = 0; i < size(); i++) {
			val += (i+1);
			if(val > threadavg) {
				finalindex[ind] = i;
				//cout << "findex " << finalindex[ind] << endl;
				ind++;
				val = 0;
			}
		}
		finalindex[threads] = size();

		for(int j = 0; j < threads; j++) {
			if(!futures[j].valid() ) {//&& (i < size()) && thread_run < threads ) {
				s = finalindex[j];
				e = finalindex[j+1];
				futures[j] = async(launch::async, multiCompare, s,e, this, &clusters, p);

			}
		}


//		}
		for(int j = 0; j < threads; j++) {
			if(futures[j].valid()) {
				//cout << "waiting on j: " << j << endl;
				futures[j].wait();
				futures[j].get();
				//cout << "got j: " << j << endl;
			}
		}

	}


	for(int i = 0; i < structures.size(); i++) {
		if(structures[i].getCluster() == -1) {
			out.addIndv(structures[i]);
		}
	}
	return out;

}


//assigns structures to pop clusters and adds to new clusters
void GASP2pop::cluster(GASP2pop &clusters, GASP2param p, int threads) {
	double average, chebyshev, euclid, localchebyshev, localavg;
	int res;

	//make two passes to make sure everything is updated nicely
	for(int n = 0; n < 2; n++) {
		//cout << "cluster pass n " << n << endl;
		//update clusters
		if (structures.size() < threads)
			threads = structures.size();
		if(threads == 0) {
			//cout << "zero thread" << endl;
			return;
		}

		//setup the futures
		vector<future<bool>> futures(threads);
		vector<int> finalindex(threads + 1);
		chrono::milliseconds timeout(0);
		chrono::milliseconds thread_wait(20);

		int thread_run = 0;
		//cout << "got there" << endl;
//		for(int i = 0; i < structures.size(); ) {
//			//get the distance to the current centroid
		int s,e;
		int totalload = ( size() * size() - size() ) / 2;
		//cout << "size " << size() << endl;
		int threadavg = (totalload / threads);
		//cout << "threadavg " << threadavg << endl;

		int val = 0;
		finalindex[0] = 0;
		int ind = 1;
		for(int i = 0; i < size(); i++) {
			val += (i+1);
			if(val > threadavg) {
				finalindex[ind] = i;
				//cout << "findex " << finalindex[ind] << endl;
				ind++;
				val = 0;
			}
		}
		finalindex[threads] = size();

		for(int j = 0; j < threads; j++) {
			if(!futures[j].valid() ) {//&& (i < size()) && thread_run < threads ) {
				s = finalindex[j];
				e = finalindex[j+1];
				futures[j] = async(launch::async, multiCompare, s,e, this, &clusters, p);

			}
		}


//		}
		for(int j = 0; j < threads; j++) {
			if(futures[j].valid()) {
				//cout << "waiting on j: " << j << endl;
				futures[j].wait();
				futures[j].get();
				//cout << "got j: " << j << endl;
			}
		}


		//the cluster is unassigned, so we make a new cluster
		for(int i = 0; i < structures.size(); i++) {
			if(structures[i].getCluster() == -1) {
				GASP2struct newcenter = structures[i];
				newcenter.resetID();
				structures[i].setCluster(clusters.size());
				clusters.addIndv(newcenter);
			}
		}

	}

}

//AML: this next function is fundamentally broken

//finds the centroid of a cluster
//FIXME: this does not handle multiple types of molecules right!
//fix needs to generalize all of the clustering algorithms
//FIXME: THIS IS ROYALLY BROKEN; I had to cripple it
//by picking the first structure in the list as the centroid
//probably going to do weird stuff, whatever, fix later
//problems is almost certainly an issue with multi-spacegroup
GASP2struct GASP2pop::cluster_center(int c) {

	//cout << "centering" << endl;
	//set up the new centroid
	GASP2struct center;
	for(int i = 0; i < structures.size(); i++) {
		if(structures[i].getCluster() == c) {
			center = structures[i];
			center.resetID();
			return center;
			break;
		}
	}
	center.resetID();

	//get the initial vector
	vector<double> sum, temp;
	int mol, dih, tempmol, tempdih;
	sum = center.getVector(mol, dih);
	for(int i = 0; i < sum.size(); i++)
		sum[i] = 0.0;
	//cout << "post sum init" << endl;

	//cout << "sumsize:"<<sum.size() << endl;

	int csize = 0;
	for(int i = 0; i < structures.size(); i++) {
		if(structures[i].getCluster() == c) {
			csize++;
			temp.clear();
			temp = structures[i].getVector(tempmol, tempdih);
			//cout << i<<":"<<temp.size() << endl;
			//need to fix silent error
			if(! vectoradd(sum, temp, mol, dih)) {
				cout << "WARNING: cluster centroid function had a problem during sum! Centroid not adjusted" << endl;
				return center;
			}
		}
	}

	//cout << "post summation" << endl;

	vectordiv(sum, static_cast<double>(csize));
	center.setVector(sum, mol, dih);
	//center.fitcell(10.0);

	//cout << "post fitcell" << endl;

	return center;
}


bool multiGroupCompare(int start, int end, GASP2pop *self, GASP2param p) {

	double average, chebyshev, euclid;
	double num;

	for(int index = start; index < end; index++) {
		self->indv(index)->setClustergroup(index);
		for(int j = 0; j < index; j++) {
			self->indv(index)->simpleCompare(*self->indv(j), p, average, chebyshev, euclid, num);
			if( average >= 0.0 )
				if(average < 2.0*p.clustersize && chebyshev < 2.0*p.chebyshevlimit) {
					self->indv(index)->setClustergroup(j);
					break;
				}
		}
	}

}

//assigns all clusters to a cluster group
//based on adjacency via inter-cluster distances
int GASP2pop::assignClusterGroups(GASP2param p, int threads) {

	double average, chebyshev, euclid;

	if (structures.size() < threads)
		threads = structures.size();
	if(threads == 0) {
		//cout << "zero thread" << endl;
		return -1;
	}

	if(size() < 2)  {
		structures[0].setClustergroup(0);
		return 1;
	}



	//setup the futures
	vector<future<bool>> futures(threads);
	vector<int> finalindex(threads + 1);
	chrono::milliseconds timeout(0);
	chrono::milliseconds thread_wait(20);

	int thread_run = 0;
	//cout << "got there" << endl;
//		for(int i = 0; i < structures.size(); ) {
//			//get the distance to the current centroid
	int s,e;
	int totalload = ( size() * size() - size() ) / 2;
	//cout << "size " << size() << endl;
	int threadavg = (totalload / threads);
	//cout << "threadavg " << threadavg << endl;

	int val = 0;
	finalindex[0] = 0;
	int ind = 1;
	for(int i = 0; i < size(); i++) {
		val += (i+1);
		if(val > threadavg) {
			finalindex[ind] = i;
			//cout << "findex " << finalindex[ind] << endl;
			ind++;
			val = 0;
		}
	}
	finalindex[threads] = size();

	//for(int k = 0; k < finalindex.size(); k++)
	//	cout << "findex " << finalindex[k] << endl;

	for(int j = 0; j < threads; j++) {
		if(!futures[j].valid() ) {//&& (i < size()) && thread_run < threads ) {
			s = finalindex[j];
			e = finalindex[j+1];
			futures[j] = async(launch::async, multiGroupCompare, s,e, this, p);
		}
	}



	for(int j = 0; j < threads; j++) {
		if(futures[j].valid()) {
			//cout << "waiting on j: " << j << endl;
			futures[j].wait();
			futures[j].get();
			//cout << "got j: " << j << endl;
		}
	}

	//finalize the groupings
	int virt, actual;
	for(int i = 0; i < size(); i++) {
		virt = structures[i].getClustergroup();
		//the clustergroup assignment is actual
		if(structures[virt].getClustergroup() == virt) {
			structures[i].setClustergroup(virt);
			break;
		}
	}

	//count it
	int count = 0;
	for(int i = 0; i < size(); i++) {
		if(structures[i].getClustergroup() == i)
			count++;

	}
	return count;

}



bool multiDistWrite(int start, int end, GASP2pop *self, GASP2param p, int num) {

	double average, chebyshev, euclid;
	double num2;

	ofstream outf;
	string name = "multidistout" + to_string(num);
	outf.open(name, ofstream::out);
	if(outf.fail()) {
		cout << mark() << "ERROR: COULD NOT OPEN FILE FOR SAVING! exiting sadly..." << endl;
		exit(1);
	}
	for(int index = start; index < end; index++) {
		//self->indv(index)->setClustergroup(index);
		for(int j = 0; j < self->size(); j++) {
			self->indv(index)->simpleCompare(*self->indv(j), p, average, chebyshev, euclid, num2);
			if(average >=0.0)
				outf << index << "," << j << "," << average << "," << chebyshev << "," << num2 << endl;
		}
	}
	outf.close();

}

//assigns all clusters to a cluster group
//based on adjacency via inter-cluster distances
void GASP2pop::allDistances(GASP2param p, int threads) {

	double average, chebyshev, euclid;

	if (structures.size() < threads)
		threads = structures.size();
	if(threads == 0) {
		//cout << "zero thread" << endl;
		return;
	}
	if(size() < 2)
		return;

	//setup the futures
	vector<future<bool>> futures(threads);
	vector<int> finalindex(threads + 1);
	chrono::milliseconds timeout(0);
	chrono::milliseconds thread_wait(20);

	int thread_run = 0;
	//cout << "got there" << endl;
//		for(int i = 0; i < structures.size(); ) {
//			//get the distance to the current centroid
	int s,e;
	int totalload = ( size() * size() - size() ) / 2;
	//cout << "size " << size() << endl;
	int threadavg = (totalload / threads);
	//cout << "threadavg " << threadavg << endl;

	int val = 0;
	finalindex[0] = 0;
	int ind = 1;
	for(int i = 0; i < size(); i++) {
		val += (i+1);
		if(val > threadavg) {
			finalindex[ind] = i;
			//cout << "findex " << finalindex[ind] << endl;
			ind++;
			val = 0;
		}
	}
	finalindex[threads] = size();



	//for(int k = 0; k < finalindex.size(); k++)
	//	cout << "findex " << finalindex[k] << endl;

	for(int j = 0; j < threads; j++) {
		if(!futures[j].valid() ) {//&& (i < size()) && thread_run < threads ) {
			s = finalindex[j];
			e = finalindex[j+1];
			futures[j] = async(launch::async, multiDistWrite, s,e, this, p, j);
		}
	}



	for(int j = 0; j < threads; j++) {
		if(futures[j].valid()) {
			//cout << "waiting on j: " << j << endl;
			futures[j].wait();
			futures[j].get();
			//cout << "got j: " << j << endl;
		}
	}


}


void GASP2pop::clusterReset() {
	for(int i = 0 ; i < size(); i++) {
		structures[i].setCluster(-1);
		structures[i].setClustergroup(-1);
	}

}


//sets the generation of a population
void GASP2pop::setGen(int gen) {

	for(int i = 0; i < size(); i++) {
		structures[i].setGen(gen);
	}

}


//FIXME FIXME FIXME
GASP2pop GASP2pop::outliers(GASP2pop &ok) {

	GASP2pop bad;

	//QE sometimes produces bad energies that are way lower than they ought to be
	//it's not clear why this happens; it's an edge case that isn't worth exploring

	//these energies are usually pretty bad, 2-3 magnitudes lower than then expected
	//energy. they need to be either recalculated or excluded in most situations;
	//typically recalculating the energies has been a successful operation.

	//find the median energy of the population
	double median = structures[size()/2].getEnergy();

	//the threshold is 1/100 the energy of the median
	//for histan this is about 1,600 kJ/mol
	//his threshold is basically well above the expected issue energy level,
	//which can be twice the total energy of the structure
	//(or, about 334,000 kJ/mol for histan)
	double medthresh = median / 100.0;

	//we onyl care about structures that are substantially
	//LOWER than the true energy
	//a structure that is substantuially higher is not unexpected
	//and also not relevant
	for(int i = 0; i < size(); i++) {
		if( (median - structures[i].getEnergy() ) > medthresh )
			bad.addIndv(structures[i]);
		else
			ok.addIndv(structures[i]);
	}

	return bad;

}

