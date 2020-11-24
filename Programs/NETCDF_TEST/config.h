#include <netcdf.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>
#include <any>

// This is the name of the data file we will read
#define FILE_NAME (char*)"/data90t/geosci/spetrov/nc_topaz/dataset-ice.nc"
//#define FILE_NAME (char*)"/data90t/geosci/spetrov/nc_topaz/dataset-single.nc"
//#define FILE_NAME (char*)"/data90t/geosci/spetrov/INMOST_ICE/Triangulation/Triangulation_thickening/concentration_data/seaice_conc_monthly_nh_f17_201003_v03r01.nc"

// Error code
#define ERRCODE 2

// Log error message 
#define ERR(e) {std::cerr << "Error:" <<  nc_strerror(e) << std::endl; exit(ERRCODE);}

// Time step, time name
#define TIME_STEP 0
#define TIME_NAME "time"

// Dim names
#define X_NAME "x"
#define Y_NAME "y"
#define Z_NAME "z"
