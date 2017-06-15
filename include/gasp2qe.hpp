//
//Copyright (c) 2017 Albert M. Lund (a.m.lund@utah.edu)
//
//This source code is subject to the license found at https://github.com/MGAC-group/MGAC2
//

#pragma once
#include "gasp2common.hpp"
#include "gasp2param.hpp"
#include "gasp2struct.hpp"
#include "gasp2.hpp"


using namespace std;

namespace QE {
	bool empty(vector<GASP2molecule>&, GASP2cell&, double&, double&, double&, time_t&, string, GASP2param);

	bool runQE(vector<GASP2molecule>&, GASP2cell&, double&, double&, double&, time_t&, string, GASP2param);

	//bool completeQE(vector<GASP2molecule>& mols, GASP2cell& unit, double &energy, double &force, double&pressure, time_t&time, string hostfile, GASP2param params) {

}
