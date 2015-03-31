#include "gasp2pop.hpp"

using namespace std;


void GASP2pop::energysort() {


}

void GASP2pop::volumesort() {


}

GASP2pop GASP2pop::newPop(int size, GAselection mode) {
	GASP2pop out;
	GASP2struct a,b;
	out.structures.reserve(size);
	int selA, selB; //selection indices
	if(mode == Roulette) {
		vector<double> values = this->scale();
		discrete_distribution<int> d(values.begin(), values.end());
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


vector<double> GASP2pop::scale() {
	vector<double> out;




}
