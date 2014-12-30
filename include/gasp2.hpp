#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <ctime>
#include <stdlib.h>

#include "gasp2param.hpp"
#include "gasp2pop.hpp"
#include "gasp2struct.hpp"
#include "gasp2qe.hpp"
#include "tinyxml2.h"
#include "svl/SVL.h"
#include "mpi.h"

using namespace std;

class GASP2control {
public:
	GASP2control();
	GASP2control(time_t start, int size, string input, string restart="");
	void server_prog() { cout << "Server says hi!" << endl;};
	void client_prog() { int i = 0; for(;;) i++; };

};
