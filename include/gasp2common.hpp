#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <ctime>
#include <stdlib.h>
#include "tinyxml2.h"
#include "svl/SVL.h"
#include "mpi.h"
#include <uuid/uuid.h>

#define PI 3.1415926535897932384626433

typedef unsigned int Index;

typedef enum Elem {
	UNK=0,C,H,N,O,P,Cl,F,S,Br,
}Elem;

double angle ( const Vec3 & A, const Vec3 & B, const Vec3 & C );

double dihedral ( const Vec3 & A, const Vec3 & B, const Vec3 & C, const Vec3 & D );

double rcov( const Elem type );

double vdw( Elem type );

Elem getElemType(string in);

