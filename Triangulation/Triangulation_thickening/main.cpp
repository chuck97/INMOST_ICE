#include <iostream>
#include <exception>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <cmath>
#include "coords_rotation.h"

extern "C" 
{
	#include "ani2d_routines.h"
	#include "mapx_routines.h"
}


static const int    nvmax        = 10'000'000;
static const int    ntmax        = 2*nvmax;
static const int    nbmax        = 1'000'000;
static const double min_tr_size  = 0.07;
static const double max_tr_size  = 0.35;
static const double min_conc     = 0.15;
static const int    min_year     = 2010;
static const int    max_year     = 2019;

static const std::string vrt_output_path = "./triangulation_output/vrt.txt";
static const std::string tri_output_path = "./triangulation_output/tri.txt";
static const std::string bnd_output_path = "./triangulation_output/bnd.txt";

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

template<typename ValueType>
ValueType** reshape_1d_2_2d(ValueType* array_1d, size_t nrows, size_t ncols)
{
	ValueType** array_2d = new ValueType*[nrows];
	for(size_t i = 0; i < nrows; ++i)
	{
		array_2d[i] = new ValueType[ncols];
	}
	
	for (size_t i = 0; i < nrows; ++i)
	{
		for (size_t j = 0; j < ncols; ++j)
		{
			array_2d[i][j] = array_1d[i*ncols + j];
		}
	}
	return array_2d;
}


std::pair<double, double> geo2xy(double lat, double lon)
{
	mapx::mapx_class* the_map;
	the_map = mapx::init_mapx((char*)"./Nps.mpp");
	double x, y;
	int status = mapx::forward_xy_mapx(the_map, lat, lon, &x, &y);
	++status;
	mapx::close_mapx(the_map);
	return(std::make_pair(x*100000.0, y*100000.0));
}

double bilinear_interpolation(double x1, double x2, 
							  double y1, double y2, 
							  double x , double y,
							  double f_x1_y1,
							  double f_x1_y2,
							  double f_x2_y1,
							  double f_x2_y2)
{
	if ((x < x1) or (x > x2) or (y < y1) or (y > y2))
	{
		throw std::invalid_argument("bilinear interpolation failed");
	}
	double f_x_y1 = ((x2 - x)/(x2 - x1))*f_x1_y1 + ((x - x1)/(x2 - x1))*f_x2_y1; //f_x1_y1 + (x - x1)*((f_x2_y1 - f_x1_y1)/(x2 - x1));
	double f_x_y2 = ((x2 - x)/(x2 - x1))*f_x1_y2 + ((x - x1)/(x2 - x1))*f_x2_y2; //f_x1_y2 + (x - x1)*((f_x2_y2 - f_x1_y2)/(x2 - x1));
	double f_x_y  = ((y2 - y)/(y2 - y1))*f_x_y1 + ((y - y1)/(y2 - y1))*f_x_y2;   //f_x_y1  + (y - y1)*((f_x_y2  -  f_x_y1)/(y2 - y1));
	return f_x_y;
}

struct grid_point
{
	double x, y, conc;
};

std::ostream& operator << (std::ostream& out, grid_point p)
{
	out << p.x << " " << p.y << " " << p.conc;
	return out;
}

class XY_conc_grid
{
public:
	void ParseHeader(const std::string& header,
					 size_t& nrows,
					 size_t& ncols,
					 double& step_m,
					 double& xmin,
					 double& xmax,
					 double& ymin,
					 double& ymax)
	{
		std::stringstream ss(header);
		std::string foo;
		std::getline(ss, foo, '=');
		ss >> ncols;
		std::getline(ss, foo, '=');
		ss >> nrows;
		std::getline(ss, foo, '=');
		ss >> xmin;
		std::getline(ss, foo, '=');
		ss >> xmax;
		std::getline(ss, foo, '=');
		ss >> ymin;
		std::getline(ss, foo, '=');
		ss >> ymax;
		std::getline(ss, foo, '=');
		ss >> step_m;
	}
	
	void PrintHeader()
	{
		std::cout << "nrows: " << nrows_ << std::endl;
		std::cout << "ncols: " << ncols_ << std::endl;
		std::cout << "grid step size (meters): " << step_m_ << std::endl;
		std::cout << "x_min (meters): " << xmin_ << std::endl;
		std::cout << "x_max (meters): " << xmax_ << std::endl;
		std::cout << "y_min (meters): " << ymin_ << std::endl;
		std::cout << "y_max (meters): " << ymax_ << std::endl;
	}
	
	void Read_conc_data(std::ifstream& input_concs)
	{
		double current_y_m = ymax_;
		double current_x_m;
		std::string current_line;
		 
		for (size_t j = 0; j < nrows_ + 1; ++j)
		{
			current_x_m = xmin_;
			std::vector<grid_point> current_row;
			getline(input_concs, current_line);
			std::stringstream ss(current_line);
			double current_conc;
			
			for (size_t i = 0; i < ncols_; ++i)
			{
				ss >> current_conc;
				current_conc = std::min(1.0, current_conc);
				current_row.push_back({current_x_m,
									   current_y_m,
									   current_conc});
				current_x_m += step_m_;
			}
			
			ss >> current_conc;
			current_row.push_back({current_x_m,
								   current_y_m,
								   current_conc});
			
			grid_concs_.push_back(current_row);
			
			current_y_m -= step_m_;
		}
	}
	
	XY_conc_grid(const std::string& path)
	{
		std::ifstream input_concs(path);
		std::string line;
		getline(input_concs, line); 
		
		ParseHeader(line, 
					nrows_,
					ncols_,
					step_m_,
					xmin_,
					xmax_,
					ymin_,
					ymax_);
		//PrintHeader();
		
		Read_conc_data(input_concs);
	}
	
	std::vector<std::vector<grid_point>>& concs()
	{
		return grid_concs_;
	}
	
	size_t nrows()
	{
		return nrows_;
	}
	
	size_t ncols()
	{
		return ncols_;
	}
	
	double step_m()
	{
		return step_m_;
	}
	
	double xmin()
	{
		return xmin_;
	}
	
	double xmax()
	{
		return xmax_;
	}

	double ymin()
	{
		return ymin_;
	}
	
	double ymax()
	{
		return ymax_;
	}
	
private:
	
	std::vector<std::vector<grid_point>> grid_concs_;
	size_t nrows_;
	size_t ncols_;
	double step_m_;
	double xmin_;
	double xmax_;
	double ymin_;
	double ymax_;
};

static XY_conc_grid concs("./concentration_data/concentrations_" + std::to_string(min_year) + "_3.txt");


double interpolate_conc(double lat, double lon)
{
	auto coords = geo2xy(lat, lon);
	double x_grid = coords.first;
	double y_grid = coords.second;
	
	if ((x_grid < concs.xmin()) or (x_grid > concs.xmax()) or 
		(y_grid < concs.ymin()) or (y_grid > concs.ymax()))
	{
		return 0.0;
	}
	
	size_t previous_row_num = (size_t)(std::floor((concs.ymax() - y_grid)/concs.step_m()));
	size_t previous_col_num = (size_t)(std::floor((x_grid - concs.xmin())/concs.step_m()));
	
	grid_point left_up, right_up, left_down, right_down;
	
	left_up = concs.concs()[previous_row_num][previous_col_num];
	right_up = concs.concs()[previous_row_num][previous_col_num + 1];
	left_down = concs.concs()[previous_row_num + 1][previous_col_num];
	right_down = concs.concs()[previous_row_num + 1][previous_col_num + 1];
	
	return bilinear_interpolation(left_down.x, right_down.x, 
								   left_down.y, left_up.y, 
								   x_grid , y_grid,
								   left_down.conc,
								   left_up.conc,
								   right_down.conc,
								   right_up.conc);
}

double meshSize(double* XY)
{
	double x = XY[0];
	double y = XY[1];
	
	Euler_rotation_info<double> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF);
	Spherical_Coords<double> model_coords(x, y);
	Spherical_Coords<double> geo_coords = 
		Rotate_Spherical(model_coords, rotation.Get_FORWARD());
		
	
	double lat = geo_coords.Get_y(); 
	double lon = geo_coords.Get_x();
	
	
	double conc = interpolate_conc(lat, lon);
	
	if (lat > 89.0)
	{
		return min_tr_size;
	}
	
	if (conc < min_conc)
	{
		return max_tr_size;
	}
	
	return ( ((min_tr_size - max_tr_size)/(1.0 - min_conc))*(conc - min_conc) + max_tr_size );
}

void merge_dataset(XY_conc_grid& new_concs)
{
	for (size_t i = 0; i < concs.nrows(); ++i)
	{
		for (size_t j = 0; j < concs.ncols(); ++j)
		{
			concs.concs()[i][j].conc = std::max(concs.concs()[i][j].conc,
											  new_concs.concs()[i][j].conc);
		}
	}
}

void add_all_dataset()
{
	{
		std::string new_dataset_path = "./concentration_data/concentrations_" + std::to_string(min_year) + "_4.txt";
		XY_conc_grid new_concs(new_dataset_path);
		merge_dataset(new_concs);
	}
	
	
	for (size_t year = (min_year + 1); year <= max_year; ++year)
	{
		{
			std::string new_dataset_path = "./concentration_data/concentrations_" + std::to_string(year) + "_3.txt";
			XY_conc_grid new_concs(new_dataset_path);
			merge_dataset(new_concs);
		}
		
		{
			std::string new_dataset_path = "./concentration_data/concentrations_" + std::to_string(year) + "_4.txt";
			XY_conc_grid new_concs(new_dataset_path);
			merge_dataset(new_concs);
		}
	} 
} 


void WriteTrInfoToFile(int nv, double* vrt, int nt, int* tri, int nb, int* bnd)
{
	double** vrt2d = reshape_1d_2_2d<double>(vrt, nv, 2);
	int** tri2d    = reshape_1d_2_2d<int>   (tri, nt, 3);
	int** bnd2d    = reshape_1d_2_2d<int>	(bnd, nb, 2);
	
	std::ofstream output_vrt(vrt_output_path);
	std::ofstream output_tri(tri_output_path);
	std::ofstream output_bnd(bnd_output_path);
	
	// vrt coords
	output_vrt << "number of verts: " << nv << std::endl;
	output_vrt.precision(10);
	
	for (size_t i = 0; i < (size_t)nv; ++i)
	{
		output_vrt << std::fixed;
		output_vrt << vrt2d[i][0] << ' ' << vrt2d[i][1] << std::endl;
	}
	
	// tri connenctivity list
	output_tri << "number of triangles: " << nt << std::endl;
	
	for (size_t i = 0; i < (size_t)nt; ++i)
	{
		output_tri << tri2d[i][0] << ' ' << tri2d[i][1] << ' ' << tri2d[i][2] << std::endl;
	}
	
	// bnd edges connectivity list
	output_bnd << "number of boundary edges: " << nb << std::endl;
	
	for (size_t i = 0; i < (size_t)nb; ++i)
	{
		output_bnd << bnd2d[i][0] << ' ' << bnd2d[i][1] << std::endl;
	}
    
    deallocate2d<double>(vrt2d, nv);
    deallocate2d<int>(tri2d, nt);
    deallocate2d<int>(bnd2d, nb);
}


int main()
{
	add_all_dataset();
	
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
	
	double h = min_tr_size;
	int ierr;
	
	ani::registersizefn_((void*)meshSize);
	
	//triangulation
	ierr = ani::aft2dboundary_(&Nbv, bv, &Nbl, bl, bltail, &h,   // geometry
																 // mesh data on output
                          &nv, vrt,        						 // coordinates of nodes
                          &nt, tri, labelT,   					 // triangles and their material
                          &nb, bnd, labelB,    					 // boundary edges and their labels
                          &nc, crv, iFNC);        				 // curved edges and parametrization (dummy)
  
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
	
	char output_ps[100] = "triangulation_thickening.ps";
	
	//visualization
	ani::graph_(&nv, vrt, &nt, tri, output_ps);
	
	WriteTrInfoToFile(nv, vrt, nt, tri, nb, bnd);
	
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
