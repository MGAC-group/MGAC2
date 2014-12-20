#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <stdlib.h>

#include "gasp2param.hpp"
#include "gasp2pop.hpp"
#include "gasp2struct.hpp"
#include "gasp2qe.hpp"
#include "tinyxml2.h"
#include "SVL.h"



class GASP2control {
	GASP2Control();
	GASP2Control(time_t start, int size, string input, string restart="");



};
