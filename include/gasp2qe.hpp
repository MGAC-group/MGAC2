#pragma once
#include "gasp2common.hpp"
#include "gasp2param.hpp"
#include "gasp2struct.hpp"


using namespace std;

namespace QE {
	bool empty(vector<GASP2molecule>&, GASP2cell&, double&, double&, double&, time_t&, string, GASP2param);

	bool runQE(vector<GASP2molecule>&, GASP2cell&, double&, double&, double&, time_t&, string, GASP2param);
}
