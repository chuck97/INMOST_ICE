#pragma once

#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>
#include <cctype>
#include <netcdf.h>
#include <map>
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

#define INMOST_ICE_ERR(message) {std::cerr << "Error: " << message  << std::endl; MPI_Finalize(); exit(1);}


#define PMF_PATH "/data90t/geosci/spetrov/INMOST_TESTS/GRID_TO_PMF/test.pmf"
#define TOPAZ_NC_FILENAME (char*)"/data90t/geosci/spetrov/nc_topaz/dataset-2020.nc"

// Error code
#define ERRCODE 2

// Log error message 
#define ERR(e) {std::cerr << "Error:" <<  nc_strerror(e) << std::endl; MPI_Finalize(); exit(ERRCODE);}

// Time step, time name
#define TIME_STEP_TOPAZ 0
#define TIME_NAME "time"

// Dim names
#define X_NAME_TOPAZ "x"
#define Y_NAME_TOPAZ "y"

// Lat lon names
#define LAT_NAME_TOPAZ "latitude"
#define LON_NAME_TOPAZ "longitude"

// Scale factor name TOPAZ
#define SCALE_FACTOR_NAME_TOPAZ "scale_factor"

// Invalid value name TOPAZ
#define INVALID_VALUE_NAME_TOPAZ "missing_value"

// Offset name TOPAZ
#define OFFSET_NAME_TOPAZ "add_offset"

// Variable names
#define A_NAME_TOPAZ "fice"
#define H_NAME_TOPAZ "hice"
#define U_NAME_TOPAZ "uice"
#define V_NAME_TOPAZ "vice"
#define UW_NAME_TOPAZ "u"
#define VW_NAME_TOPAZ "v"
#define HW_NAME_TOPAZ "ssh"
