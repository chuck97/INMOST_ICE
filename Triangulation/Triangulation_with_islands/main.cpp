#include <iostream>
#include <exception>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>

extern "C" 
{
	#include "ani2d_routines.h"
}

static const int    nvmax    = 10'000'000;
static const int    ntmax    = 2*nvmax;
static const int    nbmax    = 1'000'000;
static const double tr_size  = 0.2;

template<typename ValueType>
ValueType* allocate1d(size_t size)
{
	ValueType* array = new ValueType[size];
	return array;
}

template<typename ValueType>
void deallocate1d(ValueType* array)
{
	delete[] array;
}

template<typename ValueType>
ValueType** allocate2d(size_t nrows, size_t ncols)
{
	ValueType** array = new ValueType*[nrows];
	for (size_t i = 0; i < nrows; ++i)
	{
		array[i] = new ValueType[ncols];
	}
	return array;
}

template<typename ValueType>
void deallocate2d(ValueType** array, size_t nrows)
{
	for (size_t i = 0; i < nrows; ++i)
	{
		delete[] array[i];
	} 
}

template<typename ValueType>
ValueType* reshape_2d_2_1d(ValueType** array_2d, size_t nrows, size_t ncols)
{
	ValueType* array_1d = new ValueType[nrows*ncols];
	for(size_t i = 0; i < nrows; ++i)
	{
		for (size_t j = 0; j < ncols; ++j)
		{
			array_1d[i*ncols + j] = array_2d[i][j];
		}
	}
	return array_1d;
}

int main()
{
	std::string input_external_path = "../../Domain/data/arctic_crude_coarsened_external.txt";
	std::string input_islands_path = "../../Domain/data/arctic_crude_coarsened_islands.txt";
	std::ifstream input_external(input_external_path);
	std::ifstream input_islands(input_islands_path);
	
	std::string line;
	std::vector<std::pair<double, double>> points_external;
	std::vector<std::vector<std::pair<double, double>>> points_islands;
	
	
	//READ EXTERNAL BOUNDARY
	
	//skip first line
	getline(input_external, line);
	
	while(getline(input_external, line))
	{
		std::stringstream ss(line);
		double first, second;
		ss >> first >> second;
		points_external.push_back({first, second});
	}
	
	//remove last point (it is the same as first)
	points_external.pop_back();
	
	//READ ISLANDS
	
	int total_island_points = 0;
	std::vector<std::pair<double, double>> current_island;
	
	//skip first line
	getline(input_islands, line);
	
	while(getline(input_islands, line))
	{
		if (line[0] == '#')
		{
			current_island.pop_back();
			total_island_points += current_island.size();
			points_islands.push_back(current_island);
			current_island.clear();
		}
		else
		{
			std::stringstream ss(line);
			double first, second;
			ss >> first >> second;
			current_island.push_back({first, second});
		}
	}
	
	current_island.pop_back();
	total_island_points += current_island.size();
	points_islands.push_back(current_island);
	current_island.clear();
	
	// Make input for AFT routine
	
	int Nbv = points_external.size() + total_island_points;
	int Nbl = Nbv;
	
	std::cout << "boundary verticies: " << Nbv << std::endl;
	std::cout << "boundary edges: " << Nbl << std::endl;
	
	
	double**     bv2d = allocate2d<double>(Nbv, 2);
	int**        bl2d = allocate2d<int>(Nbl, 7);
	double** bltail2d = allocate2d<double>(Nbl, 2);
	
	// fill bv and bl
	size_t k = 0;  // global index
	
	// external points
	
	for (size_t i = 0; i < points_external.size()-1; ++i)
	{
		std::pair<double, double> p = points_external[i];
		bv2d[k][0]  = p.first;
		bv2d[k][1]  = p.second;

 		bl2d[k][0]  = k+1;
 		bl2d[k][1]  = k+2;
 		bl2d[k][2]  = 0;
 		bl2d[k][3]  = -1;
 		bl2d[k][4]  =  k+1;
 		bl2d[k][5]  =  1;
 		bl2d[k][6]  =  0;
 
		bltail2d[k][0] = 0.0;
		bltail2d[k][1] = 0.0;
		++k;
	}
	
	//last point is special
	bv2d[k][0] = points_external.back().first;
	bv2d[k][1] = points_external.back().second;
	
	
	bl2d[k][0] = k+1;
	bl2d[k][1] = 1;
	bl2d[k][2] = 0;
	bl2d[k][3] = -1;
	bl2d[k][4] = k+1;
	bl2d[k][5] = 1;
	bl2d[k][6] = 0;
	
	bltail2d[k][0] = 0.0;
	bltail2d[k][1] = 0.0;
	
	++k;
	
	// islands
	
	for (size_t i = 0; i < points_islands.size(); ++i)
	{
		size_t first_point_global_num = k;
		std::vector<std::pair<double, double>> cuurent_island = points_islands[i];
		
		for (size_t j = 0; j < points_islands[i].size()-1; ++j)
		{
			std::pair<double, double> p = cuurent_island[j];
			bv2d[k][0]  = p.first;
			bv2d[k][1]  = p.second;
	
			bl2d[k][0]  = k+1;
			bl2d[k][1]  = k+2;
			bl2d[k][2]  = 0;
			bl2d[k][3]  = -1;
			bl2d[k][4]  = k+1;
			bl2d[k][5]  = 1;
			bl2d[k][6]  = 0;
	
			bltail2d[k][0] = 0.0;
			bltail2d[k][1] = 0.0;
			++k;
		}
		//last point is special
		bv2d[k][0] = cuurent_island.back().first;
		bv2d[k][1] = cuurent_island.back().second;
		
		bl2d[k][0] = k+1;
		bl2d[k][1] = first_point_global_num+1;
		bl2d[k][2] = 0;
		bl2d[k][3] = -1;
		bl2d[k][4] = k+1;
		bl2d[k][5] = 1;
		bl2d[k][6] = 0;
		
		bltail2d[k][0] = 0.0;
		bltail2d[k][1] = 0.0;
		
		++k;
	}
	
	// reshape input
	double*     bv = reshape_2d_2_1d<double>(bv2d, Nbv, 2);
	int*        bl = reshape_2d_2_1d<int>(bl2d, Nbl, 7);
	double* bltail = reshape_2d_2_1d<double>(bltail2d, Nbl, 2);
		
	 
	//output data
	int nv, nt, nb, nc;
	
	double** vrt2d = allocate2d<double>(nvmax, 2);
	double*  vrt   = reshape_2d_2_1d<double>(vrt2d, nvmax, 2);
	
	int**    tri2d = allocate2d<int>(ntmax, 3);
	int*     tri   = reshape_2d_2_1d<int>(tri2d, ntmax, 3);
	
	int*     labelT = allocate1d<int>(ntmax);
	
	int**    bnd2d = allocate2d<int>(nbmax, 2);
	int*     bnd   = reshape_2d_2_1d<int>(bnd2d, nbmax, 2);
	  
	int*     labelB = allocate1d<int>(nbmax);
	
	double** crv2d  = allocate2d<double>(nbmax, 2);
	double*  crv    = reshape_2d_2_1d<double>(crv2d, nbmax, 2);
	
	int*     iFNC   = allocate1d<int>(nbmax);
	
	double h = tr_size;
	int ierr;
	
	//triangulation
	ierr = aft2dboundary_(&Nbv, bv, &Nbl, bl, bltail, &h,   // geometry
															// mesh data on output
                          &nv, vrt,        					// coordinates of nodes
                          &nt, tri, labelT,   				// triangles and their material
                          &nb, bnd, labelB,    				// boundary edges and their labels
                          &nc, crv, iFNC);        			// curved edges and parametrization (dummy)
    //possible errors
    if (ierr != 0)
    {
		std::cout << "error in function aft2dboundary, error code: " << ierr << std::endl;
		std::terminate();
	}
	
	std::cout << "mesh: number of triangles/vertices: " << nt << " " << nv << std::endl;
	
	if (nv > nvmax)
	{
		std::cout << "to many nodes" << std::endl;
		std::terminate();
	}
	
	if (nt > ntmax)
	{
		std::cout << "to many triangles" << std::endl;
		std::terminate();
	}
	
	if (nb > nbmax)
	{
		std::cout << "to many boundary edges" << std::endl;
		std::terminate();
	}
	
	char output_ps[100] = "full.ps";
	
	//visualization
	graph_(&nv, vrt, &nt, tri, output_ps);
	
	deallocate2d<double>(bv2d, Nbv);
    deallocate1d<double>(bv);
    deallocate2d<int>(bl2d, Nbl);
    deallocate1d<int>(bl);
    deallocate2d<double>(bltail2d, Nbl);
    deallocate1d<double>(bltail);
    deallocate2d<double>(vrt2d, nvmax);
    deallocate1d<double>(vrt);
	deallocate2d<int>(tri2d, ntmax);
	deallocate1d<int>(tri);
	deallocate1d<int>(labelT);
	deallocate2d<int>(bnd2d, nbmax);
	deallocate1d<int>(bnd);
	deallocate1d<int>(labelB);
	deallocate2d<double>(crv2d, nbmax);
	deallocate1d<double>(crv);
	deallocate1d<int>(iFNC);
                            
    return 0;
}
