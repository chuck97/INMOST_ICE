#pragma once

#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <map>
#include <iomanip>
#include <functional>
#include <limits>
#include "inmost.h"
#include <nlohmann/json.hpp>

#define USE_MPI
#define USE_PARTITIONER
#define USE_PARTITIONER_PARMETIS
#define USE_SOLVER


#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

#define INMOST_ICE_ERR(message) {std::cerr << "Error: " << message  << std::endl; MPI_Finalize(); exit(1);}


#define PMF_PATH "/data90t/geosci/spetrov/TRANSPORT_TESTS/PMF_3D_GENERATION/build/grid.pmf"
#define MASS_INT_PATH "../results/MassInt.txt"

// Error code
#define ERRCODE 2

// Log error message 
#define INMOST_ICE_ERR(message) {std::cerr << "Error: " << message  << std::endl; MPI_Finalize(); exit(1);}

////////////////////////////////////////////////////////////////////////////

// Virtual Earth radius
#define EARTH_RADIUS 6400000.0

////////////////////////////////////////////////////////////////////////////

// First scalar field parameters
#define COSINE_BELL_AMPLITUDE 1.0
#define COSINE_BELL_SCALE_FACTOR 0.9
#define COSINE_BELL_RADIUS 0.5
#define COSINE_BELL_BACKGROUND 0.1
#define FIRST_BELL_CENTER_LON M_PI
#define FIRST_BELL_CENTER_LAT M_PI/3.0
#define SECOND_BELL_CENTER_LON M_PI
#define SECOND_BELL_CENTER_LAT -M_PI/3.0

////////////////////////////////////////////////////////////////////////////

// Second scalar field parameters
#define GAUSSIAN_HAT_AMPLITUDE 1.0
#define GAUSSIAN_HAT_WIDTH 5.0
#define FIRST_HAT_CENTER_LON 3.0*M_PI/4.0
#define FIRST_HAT_CENTER_LAT 0.0
#define SECOND_HAT_CENTER_LON 5.0*M_PI/4.0 
#define SECOND_HAT_CENTER_LAT 0.0

////////////////////////////////////////////////////////////////////////////

// Third scalar field parameters
#define CYLINDER_SCALE_FACTOR  1.0
#define CYLINDER_BACKGROUND  0.1
#define CYLINDER_RADIUS 0.5
#define FIRST_CYLINDER_CENTER_LON 3.0*M_PI/4.0
#define FIRST_CYLINDER_CENTER_LAT 0.0
#define SECOND_CYLINDER_CENTER_LON 5.0*M_PI/4.0 
#define SECOND_CYLINDER_CENTER_LAT 0.0



// vareps param 
#define VAREPS 1e-15