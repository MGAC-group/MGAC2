#include "gasp2pop.hpp"

using namespace std;


void GASP2pop::energysort() {


}

void GASP2pop::volumesort() {


}


//newPop assumes that the scaling is valid if the size of the
//current population size matches the scaling size
//if the sizes are different, it applies a constant scaling
//size the ranking of the population is unknown
GASP2pop GASP2pop::newPop(int size, GAselection mode) {
	GASP2pop out;
	GASP2struct a,b;
	out.structures.reserve(size);
	int selA, selB; //selection indices
	if(mode == Roulette) {
		if(scaling.size() != structures.size())
			scale(1.0,0.0,0.0);
		discrete_distribution<int> d(scaling.begin(), scaling.end());
		for(int i = 0; i < size; i++) {
			selA = d(rgen);
			selB = d(rgen);
			structures[selA].crossStruct(structures[selB], a,b);
			out.structures.push_back(a);
		}
	}
	else { //pattern
		cout << "Pattern mode not implemented!" << endl;
	}
	return out;
}

GASP2pop GASP2pop::fullCross() {
	GASP2pop out;
	GASP2struct a,b;
	int size = structures.size();
	out.structures.reserve(size*size - size);

	for(int i = 0; i < size; i++) {
		for(int j = i; j < size; j++) {
			if(i==j) continue;
			structures[i].crossStruct(structures[j],a,b);
			out.structures.push_back(a);
			out.structures.push_back(b);
		}
	}
	return out;
}

void GASP2pop::addIndv(int add) {

	int init = structures.size();
	structures.resize(init + add);
	for(int i = init; i < (init+add); i++)
		structures[i].init();

}

void GASP2pop::addIndv(GASP2pop add) {
	structures.reserve(structures.size() + add.structures.size());
	structures.insert(structures.end(), add.structures.begin(), add.structures.end());

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
vector<double> scale(double con, double lin, double exp) {
	scaling.clear();
	double val, diffE, diffV;
	vector<double> vol, ener;
	double minE = numeric_limits<double>::max(), maxE = numeric_limits<double>::min();
	double minV = numeric_limits<double>::max(), maxV = numeric_limits<double>::min();

	for(int i = 0; i < structures.size(); i++) {
		vol.push_back(stuctures[i].getVolScore());
		ener.push_back(structures[i].getEnergy());
		if(vol.back() > maxV) maxV = vol.back();
		if(vol.back() < minV) minV = vol.back();
		//we assume energies will usually be negative, since E = 0.0 is more or less impossible
		//since some structures to no have energies calculated, we do not want to influence
		//the energy
		if(ener.back() < 0.0) {
			if(ener.back() > maxE) maxE = ener.back();
			if(ener.back() < minE) minE = ener.back();
		}
	}

	double volscore, enerscore;
	diffE = maxE-minE;
	diffV = maxV-minV;
	for(int i = 0; i < structure.size(); i++) {

		volscore = (vol[i] - minV)/diffV;
		val = con; //assign a constant minimum to all values
		val += lin*(1.0-volscore); //add the linear component
		val += exp*(std::exp(-0.5*10.0*volscore));//add the exponential component
		if(ener[i] < 0.0) { //add energy components if the energy is valid
			enerscore = (ener[i] - minE)/diffE;
			val += lin*(1.0-enerscore) + exp*std::exp(-0.5*10.0*enerscore);
		}

		scaling.push_back(val);
	}


}
