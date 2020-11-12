#pragma once

#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>
#include <cctype>
#include "inmost.h"

#define USE_MPI
#define USE_PARTITIONER
#define USE_PARTITIONER_PARMETIS
#define USE_SOLVER


#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

#define ERR(message) {std::cerr << "Error: " << message  << std::endl; MPI_Finalize(); exit(1);}


#define PMF_PATH "/data90t/geosci/spetrov/INMOST_TESTS/GRID_TO_PMF/test.pmf"


