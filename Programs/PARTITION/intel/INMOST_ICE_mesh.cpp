#include "INMOST_ICE_mesh.h"
#include "coords_rotation.h"
using namespace INMOST;

// trim from start (in place)
static inline void ltrim(std::string &s) 
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), 
    [](unsigned char ch) 
    {
        return !std::isspace(ch);
    }));
};

// trim from end (in place)
static inline void rtrim(std::string &s) 
{
    s.erase(std::find_if(s.rbegin(), s.rend(), 
    [](unsigned char ch) 
    {
        return !std::isspace(ch);
    }).base(), s.end());
};

// trim from both ends (in place)
static inline void trim(std::string &s) 
{
    ltrim(s);
    rtrim(s);
};

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
	double f_x_y1 = ((x2 - x)/(x2 - x1))*f_x1_y1 + ((x - x1)/(x2 - x1))*f_x2_y1; 
	double f_x_y2 = ((x2 - x)/(x2 - x1))*f_x1_y2 + ((x - x1)/(x2 - x1))*f_x2_y2; 
	double f_x_y  = ((y2 - y)/(y2 - y1))*f_x_y1 + ((y - y1)/(y2 - y1))*f_x_y2;   
	return f_x_y;
}


INMOST_ICE_mesh::INMOST_ICE_mesh(const std::string& mesh_path)
{
    Mesh::Initialize(NULL, NULL);

#if defined(USE_PARTITIONER)
    Partitioner::Initialize(NULL, NULL);
#endif

#if defined(USE_SOLVER)
    Solver::Initialize(NULL, NULL, ""); 
#endif

    double ttt = Timer();
    bool repartition = false;
    ice_mesh = new Mesh();

#if defined(USE_MPI)
    ice_mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif

    ttt = Timer();

    if(ice_mesh->isParallelFileFormat(mesh_path))
    {
        ice_mesh->Load(mesh_path); 
        repartition = true;

        if (ice_mesh->GetProcessorRank()==0)
        {
            std::cout << "Parallel realization" << std::endl;
        }
    }
    else
    {
        if (ice_mesh->GetProcessorRank() == 0)
        {
            ice_mesh->Load(mesh_path);
            std::cout << "Serial realization" << std::endl;
        }
    }

    if (ice_mesh->GetProcessorRank()==0)
    {
        std::cout << "Mesh " << mesh_path 
                  <<  " loaded successfully: " 
                  <<  Timer() - ttt << std::endl;

        std::cout << "There are " 
                  << ice_mesh->GetProcessorsNumber() 
                  <<" processes" <<  std::endl;
    }
    BARRIER
    
    Partition();
};

void INMOST_ICE_mesh::Partition()
{
    double ttt = Timer();
#if defined(USE_PARTITIONER)
	if (ice_mesh->GetProcessorsNumber() > 1) // need repartition
	{ 
        ttt = Timer();
	    Partitioner* p = new Partitioner(ice_mesh);
#ifdef USE_PARTITIONER_PARMETIS
        p->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
#elif USE_PARTITIONER_ZOLTAN
        p->SetMethod(Partitioner::Zoltan, Partitioner::Partition);
#else
        p->SetMethod(Partitioner::INNER_KMEANS, Partitioner::Partition);
#endif
		p->Evaluate();
		delete p;
		BARRIER

	    if(ice_mesh->GetProcessorRank() == 0)
        {
            std::cout << "Mesh partition: " << Timer()-ttt << std::endl;
        } 

        BARRIER
		ttt = Timer();
		ice_mesh->Redistribute(); 
		ice_mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
        ice_mesh->ExchangeGhost(1, NODE);
		BARRIER

		if(ice_mesh->GetProcessorRank() == 0)
        {
        	std::cout << "Redistribute initial mesh data: " << Timer()-ttt << std::endl;
        } 
        BARRIER
        ice_mesh->AssignGlobalID(CELL|EDGE|FACE|NODE);
	}
#endif
};

void INMOST_ICE_mesh::PrintPMF(const std::string& filename)
{
    std::string no_spaces_filename = filename;
    rtrim(no_spaces_filename);
    if(no_spaces_filename.substr(no_spaces_filename.size()-4, 4) != ".pmf")
    {
        INMOST_ICE_ERR("Can't write mesh data to this file. Filename should ended by .pmf")
    }
    BARRIER
    double ttt = Timer();
	ice_mesh->Save(no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << " : " << Timer() - ttt << std::endl; 
    }
};


void INMOST_ICE_mesh::PrintPVTK(const std::string& filename)
{
    std::string no_spaces_filename = filename;
    rtrim(no_spaces_filename);
    if(no_spaces_filename.substr(no_spaces_filename.size()-5, 5) != ".pvtk")
    {
        INMOST_ICE_ERR("Can't write mesh to this file. Filename should ended by .pvtk")
    }
    BARRIER
    double ttt = Timer();
	ice_mesh->Save(no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << " : " << Timer() - ttt << std::endl; 
    }
};


INMOST_ICE_mesh::~INMOST_ICE_mesh()
{
#if defined(USE_PARTITIONER)
    Partitioner::Finalize();
#endif

    delete ice_mesh;

#if defined(USE_SOLVER)
    Solver::Finalize();
#endif
    
    Mesh::Finalize();
};

INMOST::Mesh* INMOST_ICE_mesh::GetMesh()
{
    return ice_mesh;
};


INMOST_ICE_nodes::INMOST_ICE_nodes(INMOST::Mesh* m)
{
    ice_mesh = m;
    double ttt = Timer();

    nc.model_coords = ice_mesh->CreateTag("model coords", DATA_REAL, NODE, NONE, 2);
    nc.geo_coords = ice_mesh->CreateTag("geo coords", DATA_REAL, NODE, NONE, 2);
    nc.topaz_coords = ice_mesh->CreateTag("TOPAZ coords", DATA_REAL, NODE, NONE, 2);    

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        // Assign model coords
		nodeit->RealArray(nc.model_coords)[0] = nodeit->Coords()[0];
        nodeit->RealArray(nc.model_coords)[1] = nodeit->Coords()[1];

        
        // Assign geo coords
        std::vector<double> geo = from_model_2_geo<double>(nodeit->RealArray(nc.model_coords)[0],
                                                           nodeit->RealArray(nc.model_coords)[1]);
        nodeit->RealArray(nc.geo_coords)[0] = geo[0];
        nodeit->RealArray(nc.geo_coords)[1] = geo[1];
        
        // Assign TOPAZ coords
        std::vector<double> topaz = from_geo_2_topaz(nodeit->RealArray(nc.geo_coords)[0],
                                                     nodeit->RealArray(nc.geo_coords)[1]);
        nodeit->RealArray(nc.topaz_coords)[0] = topaz[0]/1e5;
        nodeit->RealArray(nc.topaz_coords)[1] = topaz[1]/1e5;

        // Assign 
    }
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Assigning coords: " <<  Timer() - ttt <<std::endl;
    }
};

NodeCoords& INMOST_ICE_nodes::GetCoords()
{
    return nc;
};

INMOST::Mesh* INMOST_ICE_nodes::GetMesh()
{
    return ice_mesh;
};

void INMOST_ICE_nodes::InterpolateScalarFromNETCDF(const std::string& filename,
                                                   const std::string& xname,
                                                   const std::string& yname,
                                                   int netcdf_time_step,
                                                   const std::string& node_variable_name,
                                                   const std::string& scale_factor_name,
                                                   const std::string& invalid_value_name,
                                                   const std::string& offset_name,
                                                   double invalid_value_fill,
                                                   double no_extrapolation_fill,
                                                   const std::string& nc_variable_name,
                                                   double max_abs_value,
                                                   bool is_depth,
                                                   INMOST::Tag netcdf_coords
                                                   )
{
    if (node_data_tags.count(node_variable_name) == 0)
    {
        INMOST::Tag t;
        t = ice_mesh->CreateTag(node_variable_name, DATA_REAL, NODE, NONE, 1);
        node_data_tags[node_variable_name] = t;
    }

    //Initialize tag with no interpolation value
    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        nodeit->Real(node_data_tags[node_variable_name]) = no_extrapolation_fill;
    }

    //Find extremal coords

    double min_topaz_x;
    double min_topaz_y;
    double max_topaz_x;
    double max_topaz_y;

    double tmp_min_x = ice_mesh->BeginNode()->RealArray(netcdf_coords)[0];
    double tmp_max_x = ice_mesh->BeginNode()->RealArray(netcdf_coords)[0];
    double tmp_min_y = ice_mesh->BeginNode()->RealArray(netcdf_coords)[1];
    double tmp_max_y = ice_mesh->BeginNode()->RealArray(netcdf_coords)[1];

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double current_x = nodeit->RealArray(netcdf_coords)[0];
            double current_y = nodeit->RealArray(netcdf_coords)[1];
            if (current_x <= tmp_min_x)
            {
                tmp_min_x = current_x;
            }

            if (current_x >= tmp_max_x)
            {
                tmp_max_x = current_x;
            }

            if (current_y <= tmp_min_y)
            {
                tmp_min_y = current_y;
            }

            if (current_y >= tmp_max_y)
            {
                tmp_max_y = current_y;
            }
        }
    }
    min_topaz_x = tmp_min_x;
    min_topaz_y = tmp_min_y;
    max_topaz_x = tmp_max_x;
    max_topaz_y = tmp_max_y;

    BARRIER

    // Get coords from netcdf

    int retval;
    int fileid;

    int x_id, y_id;
    size_t x_size, y_size;

    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &fileid)))
        ERR(retval);

    if((retval = nc_inq_dimid(fileid, xname.c_str(), &x_id)))
        ERR(retval);

    if((retval = nc_inq_dimid(fileid, yname.c_str(), &y_id)))
        ERR(retval);
    
    if((retval = nc_inq_dimlen(fileid, x_id, &x_size)))
        ERR(retval);

    if((retval = nc_inq_dimlen(fileid, y_id, &y_size)))
        ERR(retval);
    
    double x_coords[x_size];
    double y_coords[y_size];

    int x_coords_id, y_coords_id;

    if ((retval = nc_inq_varid(fileid, xname.c_str(), &x_coords_id)))
        ERR(retval);
    
    if ((retval = nc_inq_varid(fileid, yname.c_str(), &y_coords_id)))
        ERR(retval);

    if ((retval = nc_get_var_double(fileid, x_coords_id, x_coords)))
        ERR(retval);
    
    if ((retval = nc_get_var_double(fileid, y_coords_id, y_coords)))
        ERR(retval);

    double dx = x_coords[1] - x_coords[0];
    double dy = y_coords[1] - y_coords[0];

    // find x start count
    size_t x_start;
    size_t x_count;

    double max_x_begin = std::max(x_coords[0], min_topaz_x);
    double min_x_end = std::min(x_coords[x_size -1], max_topaz_x);

    if (max_x_begin > min_x_end)
    {
        x_count = 0;
        x_start = 0;
    }
    else if (max_x_begin == x_coords[0])
    {
        x_start = 0;
        x_count = std::min(size_t((max_topaz_x - x_coords[0])/dx) + 1, x_size);
    }
    else if (max_x_begin == min_topaz_x)
    {
        x_start = size_t((min_topaz_x - x_coords[0])/dx);
        x_count = std::min((size_t)((max_topaz_x - (x_coords[0] + dx*x_start))/dx) + 1, (x_size - x_start));
    }
    else
    {
        INMOST_ICE_ERR("bad x_start x_count procedure");
    }

    // find y start count
    size_t y_start;
    size_t y_count;

    double max_y_begin = std::max(y_coords[0], min_topaz_y);
    double min_y_end = std::min(y_coords[y_size -1], max_topaz_y);

    if (max_y_begin > min_y_end)
    {
        y_count = 0;
        y_start = 0;
    }
    else if (max_y_begin == y_coords[0])
    {
        y_start = 0;
        y_count = std::min(size_t((max_topaz_y - y_coords[0])/dy) + 1, y_size);
    }
    else if (max_y_begin == min_topaz_y)
    {
        y_start = size_t((min_topaz_y - y_coords[0])/dy);
        y_count = std::min((size_t)((max_topaz_y - (y_coords[0] + dy*y_start))/dy) + 1, (y_size - y_start));
    }
    else
    {
        INMOST_ICE_ERR("bad y_start y_count procedure");
    }

    // Find data id
    int data_variable_id;

    if ((retval = nc_inq_varid(fileid, nc_variable_name.c_str(), &data_variable_id)))
    {
        ERR(retval);
    }

    //Get scale factor 
    double scale_factor;
    if (scale_factor_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable_id, scale_factor_name.c_str(), &scale_factor)))
        {
            ERR(retval);
        }
    }
    else
    {
        scale_factor = 1.0;
    }
    
    //Get invalid value
    nc_type atttype;
    
    int invalid_value_int;
    short invalid_value_short;
    double invalid_value_double;

    if (invalid_value_name.size() != 0)
    {
        if ((retval = nc_inq_att(fileid, data_variable_id, invalid_value_name.c_str(), &atttype, NULL)))
            ERR(retval);

        if (atttype == NC_SHORT)
        {
            if ((retval = nc_get_att(fileid, data_variable_id, invalid_value_name.c_str(), &invalid_value_short)))
                ERR(retval);
        }
        else if ((atttype == NC_INT) or (atttype == NC_LONG))
        {
            if ((retval = nc_get_att(fileid, data_variable_id, invalid_value_name.c_str(), &invalid_value_int)))
                ERR(retval);
        }
        else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
        {
            if ((retval = nc_get_att(fileid, data_variable_id, invalid_value_name.c_str(), &invalid_value_double)))
            {
                ERR(retval);
            }
        }
        else
        {
            INMOST_ICE_ERR("att type is not int or double, can't read");
        }
    }
    else
    {
        invalid_value_int = std::nan("1");
        invalid_value_double = std::nan("1");
        invalid_value_short = std::nan("1");
    }

    //Get offset value
    double offset_value;
    if (offset_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable_id, offset_name.c_str(), &offset_value)))
        {
            ERR(retval);
        }
    }
    else
    {
        offset_value = 0.0;
    }
    

    // Assemble start and count array
    size_t* start;
    size_t* count;

    if (is_depth)
    {
        start = new size_t[4];
        count = new size_t[4];
        start[0] = netcdf_time_step;
        start[1] = 0;
        start[2] = y_start;
        start[3] = x_start;
        count[0] = 1;
        count[1] = 1;
        count[2] = y_count;
        count[3] = x_count;
    }
    else
    {
        start = new size_t[3];
        count = new size_t[3];
        start[0] = netcdf_time_step;
        start[1] = y_start;
        start[2] = x_start;
        count[0] = 1;
        count[1] = y_count;
        count[2] = x_count;
    }

    // Get data in parallel
    int data_local_int[y_count][x_count];
    short data_local_short[y_count][x_count];
    double data_local_double[y_count][x_count];

    if((retval = nc_inq_var(fileid, data_variable_id, NULL, &atttype, NULL, NULL, NULL)))
        ERR(retval)

    
    if ((atttype == NC_SHORT))
    {
        if ((retval = nc_get_vara_short(fileid, data_variable_id, start, count, &data_local_short[0][0])))
            ERR(retval)
    }
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        if ((retval = nc_get_vara_int(fileid, data_variable_id, start, count, &data_local_int[0][0])))
            ERR(retval)
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        if ((retval = nc_get_vara_double(fileid, data_variable_id, start, count, &data_local_double[0][0])))
            ERR(retval)
    }
    else
    {
        INMOST_ICE_ERR("var type is not int or double, can't read");
    }

    // convert data to double
    if (atttype == NC_SHORT)
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val;
                current_val = (data_local_short[j][i] == invalid_value_short) ? invalid_value_fill 
                                 : (double)(data_local_short[j][i])*scale_factor + offset_value;
                data_local_double[j][i] = current_val;
            }
        }
    }           
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val;
                current_val = (data_local_int[j][i] == invalid_value_int) ? invalid_value_fill 
                                 : (double)(data_local_int[j][i])*scale_factor + offset_value;
                data_local_double[j][i] = current_val;
            }
        }
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val;
                current_val = (data_local_double[j][i] == invalid_value_double) ? invalid_value_fill 
                                     : (data_local_double[j][i])*scale_factor + offset_value;
                data_local_double[j][i] = current_val;
            }
        }
    }
    else
    {
        INMOST_ICE_ERR("unknown value type");
    }
    
    // get square variables     
    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
    {
        double local_x_start = x_coords[x_start];
        double local_y_start = y_coords[y_start];
        double local_x_end = x_coords[x_start + x_count];
        double local_y_end = y_coords[y_start + y_count];

        if(nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(nc.topaz_coords)[0];
            double y = nodeit->RealArray(nc.topaz_coords)[1];

            if ((x < local_x_start) or
                (x > local_x_end)   or
                (y < local_y_start) or
                (y > local_y_end))
            {
                nodeit->Real(node_data_tags[node_variable_name]) = no_extrapolation_fill;
                continue;
            }
            else
            {
                size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                size_t y_prev_pos = (size_t)((y - local_y_start)/dy);
                
                double data_ld = data_local_double[y_prev_pos][x_prev_pos];
                double data_lu = data_local_double[y_prev_pos+1][x_prev_pos];
                double data_rd = data_local_double[y_prev_pos][x_prev_pos+1];
                double data_ru = data_local_double[y_prev_pos+1][x_prev_pos+1];
                
                double xl = x_coords[x_start + x_prev_pos];
                double xr = x_coords[x_start + x_prev_pos + 1];
                double yd = y_coords[y_start + y_prev_pos];
                double yu = y_coords[y_start + y_prev_pos + 1];

                double curr_data = bilinear_interpolation(xl, xr, 
		                                                  yd, yu, 
		                                                  x , y,
		                                                  data_ld, data_lu,
                                                          data_rd, data_ru);

				nodeit->Real(node_data_tags[node_variable_name]) = 
                    (fabs(curr_data) > max_abs_value) ? invalid_value_fill : curr_data;
            }
        }
    }
    BARRIER;

    // Exchange data to ghost cells
    ice_mesh->ExchangeData(node_data_tags[node_variable_name], NODE, 0);

    // Close file
    if ((retval = nc_close(fileid)))
        ERR(retval);
};

void INMOST_ICE_nodes::InterpolateVectorFromNETCDF(
                                     const std::string& filename,
                                     const std::string& xname,
                                     const std::string& yname,
                                     const std::string& lonname,
                                     const std::string& latname,
                                     int netcdf_time_step,
                                     const std::string& node_variable_name,
                                     const std::string& scale_factor_name,
                                     const std::string& invalid_value_name,
                                     const std::string& offset_name,
                                     double invalid_value_fill,
                                     double no_extrapolation_fill,
                                     const std::string& nc_variable_name1,
                                     const std::string& nc_variable_name2,
                                     double max_abs_value,
                                     bool is_depth,
                                     INMOST::Tag netcdf_coords
                                     )
{
    if (node_data_tags.count(node_variable_name) == 0)
    {
        INMOST::Tag t;
        t = ice_mesh->CreateTag(node_variable_name, DATA_REAL, NODE, NONE, 2);
        node_data_tags[node_variable_name] = t;
    }

    //Fill tag with no interpolation value
    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        nodeit->RealArray(node_data_tags[node_variable_name])[0] = no_extrapolation_fill;
        nodeit->RealArray(node_data_tags[node_variable_name])[1] = no_extrapolation_fill;
    }

    //Find extremal coords

    double min_topaz_x;
    double min_topaz_y;
    double max_topaz_x;
    double max_topaz_y;

    double tmp_min_x = ice_mesh->BeginNode()->RealArray(netcdf_coords)[0];
    double tmp_max_x = ice_mesh->BeginNode()->RealArray(netcdf_coords)[0];
    double tmp_min_y = ice_mesh->BeginNode()->RealArray(netcdf_coords)[1];
    double tmp_max_y = ice_mesh->BeginNode()->RealArray(netcdf_coords)[1];

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double current_x = nodeit->RealArray(netcdf_coords)[0];
            double current_y = nodeit->RealArray(netcdf_coords)[1];
            if (current_x <= tmp_min_x)
            {
                tmp_min_x = current_x;
            }

            if (current_x >= tmp_max_x)
            {
                tmp_max_x = current_x;
            }

            if (current_y <= tmp_min_y)
            {
                tmp_min_y = current_y;
            }

            if (current_y >= tmp_max_y)
            {
                tmp_max_y = current_y;
            }
        }
    }
    min_topaz_x = tmp_min_x;
    min_topaz_y = tmp_min_y;
    max_topaz_x = tmp_max_x;
    max_topaz_y = tmp_max_y;

    BARRIER

    // Get coords from netcdf

    int retval;
    int fileid;

    int x_id, y_id;
    size_t x_size, y_size;

    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &fileid)))
        ERR(retval);

    if((retval = nc_inq_dimid(fileid, xname.c_str(), &x_id)))
        ERR(retval);

    if((retval = nc_inq_dimid(fileid, yname.c_str(), &y_id)))
        ERR(retval);
    
    if((retval = nc_inq_dimlen(fileid, x_id, &x_size)))
        ERR(retval);

    if((retval = nc_inq_dimlen(fileid, y_id, &y_size)))
        ERR(retval);
    
    double x_coords[x_size];
    double y_coords[y_size];

    int x_coords_id, y_coords_id;

    if ((retval = nc_inq_varid(fileid, xname.c_str(), &x_coords_id)))
        ERR(retval);
    
    if ((retval = nc_inq_varid(fileid, yname.c_str(), &y_coords_id)))
        ERR(retval);

    if ((retval = nc_get_var_double(fileid, x_coords_id, x_coords)))
        ERR(retval);
    
    if ((retval = nc_get_var_double(fileid, y_coords_id, y_coords)))
        ERR(retval);

    double dx = x_coords[1] - x_coords[0];
    double dy = y_coords[1] - y_coords[0];

    // find x start count
    size_t x_start;
    size_t x_count;

    double max_x_begin = std::max(x_coords[0], min_topaz_x);
    double min_x_end = std::min(x_coords[x_size -1], max_topaz_x);

    if (max_x_begin > min_x_end)
    {
        x_count = 0;
        x_start = 0;
    }
    else if (max_x_begin == x_coords[0])
    {
        x_start = 0;
        x_count = std::min(size_t((max_topaz_x - x_coords[0])/dx) + 1, x_size);
    }
    else if (max_x_begin == min_topaz_x)
    {
        x_start = size_t((min_topaz_x - x_coords[0])/dx);
        x_count = std::min((size_t)((max_topaz_x - (x_coords[0] + dx*x_start))/dx) + 1, (x_size - x_start));
    }
    else
    {
        INMOST_ICE_ERR("bad x_start x_count procedure");
    }

    // find y start count
    size_t y_start;
    size_t y_count;

    double max_y_begin = std::max(y_coords[0], min_topaz_y);
    double min_y_end = std::min(y_coords[y_size -1], max_topaz_y);

    if (max_y_begin > min_y_end)
    {
        y_count = 0;
        y_start = 0;
    }
    else if (max_y_begin == y_coords[0])
    {
        y_start = 0;
        y_count = std::min(size_t((max_topaz_y - y_coords[0])/dy) + 1, y_size);
    }
    else if (max_y_begin == min_topaz_y)
    {
        y_start = size_t((min_topaz_y - y_coords[0])/dy);
        y_count = std::min((size_t)((max_topaz_y - (y_coords[0] + dy*y_start))/dy) + 1, (y_size - y_start));
    }
    else
    {
        INMOST_ICE_ERR("bad y_start y_count procedure");
    }

    // Find data id
    int data_variable1_id;
    int data_variable2_id;

    int lon_id;
    int lat_id;

    if ((retval = nc_inq_varid(fileid, nc_variable_name1.c_str(), &data_variable1_id)))
    {
        ERR(retval);
    }

    if ((retval = nc_inq_varid(fileid, nc_variable_name2.c_str(), &data_variable2_id)))
    {
        ERR(retval);
    }

    if ((retval = nc_inq_varid(fileid, lonname.c_str(), &lon_id)))
    {
        ERR(retval);
    }

    if ((retval = nc_inq_varid(fileid, latname.c_str(), &lat_id)))
    {
        ERR(retval);
    }

    //Get scale factor 
    double scale_factor;
    if (scale_factor_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable1_id, scale_factor_name.c_str(), &scale_factor)))
        {
            ERR(retval);
        }
    }
    else
    {
        scale_factor = 1.0;
    }
    
    //Get invalid value
    nc_type atttype;
    
    int invalid_value_int;
    short invalid_value_short;
    double invalid_value_double;

    if (invalid_value_name.size() != 0)
    {
        if ((retval = nc_inq_att(fileid, data_variable1_id, invalid_value_name.c_str(), &atttype, NULL)))
            ERR(retval);

        if (atttype == NC_SHORT)
        {
            if ((retval = nc_get_att(fileid, data_variable1_id, invalid_value_name.c_str(), &invalid_value_short)))
                ERR(retval);
        }
        else if ((atttype == NC_INT) or (atttype == NC_LONG))
        {
            if ((retval = nc_get_att(fileid, data_variable1_id, invalid_value_name.c_str(), &invalid_value_int)))
                ERR(retval);
        }
        else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
        {
            if ((retval = nc_get_att(fileid, data_variable1_id, invalid_value_name.c_str(), &invalid_value_double)))
            {
                ERR(retval);
            }
        }
        else
        {
            INMOST_ICE_ERR("att type is not int or double, can't read");
        }
    }
    else
    {
        invalid_value_int = std::nan("1");
        invalid_value_double = std::nan("1");
        invalid_value_short = std::nan("1");
    }

    //Get offset value
    double offset_value;
    if (offset_name.size() != 0)
    {
        if ((retval = nc_get_att(fileid, data_variable1_id, offset_name.c_str(), &offset_value)))
        {
            ERR(retval);
        }
    }
    else
    {
        offset_value = 0.0;
    }

    // Assemble start and count array
    size_t* start;
    size_t* count;

    if (is_depth)
    {
        start = new size_t[4];
        count = new size_t[4];
        start[0] = netcdf_time_step;
        start[1] = 0;
        start[2] = y_start;
        start[3] = x_start;
        count[0] = 1;
        count[1] = 1;
        count[2] = y_count;
        count[3] = x_count;
    }
    else
    {
        start = new size_t[3];
        count = new size_t[3];
        start[0] = netcdf_time_step;
        start[1] = y_start;
        start[2] = x_start;
        count[0] = 1;
        count[1] = y_count;
        count[2] = x_count;
    }

    // Get data in parallel
    int vector_data_local_int1[y_count][x_count];
    int vector_data_local_int2[y_count][x_count];

    short vector_data_local_short1[y_count][x_count];
    short vector_data_local_short2[y_count][x_count];

    double vector_data_local_double1[y_count][x_count];
    double vector_data_local_double2[y_count][x_count];

    if((retval = nc_inq_var(fileid, data_variable1_id, NULL, &atttype, NULL, NULL, NULL)))
        ERR(retval)

    
    if ((atttype == NC_SHORT))
    {
        if ((retval = nc_get_vara_short(fileid, data_variable1_id, start, count, &vector_data_local_short1[0][0])))
            ERR(retval)
        
        if ((retval = nc_get_vara_short(fileid, data_variable2_id, start, count, &vector_data_local_short2[0][0])))
            ERR(retval)
    }
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        if ((retval = nc_get_vara_int(fileid, data_variable1_id, start, count, &vector_data_local_int1[0][0])))
            ERR(retval)
        
        if ((retval = nc_get_vara_int(fileid, data_variable2_id, start, count, &vector_data_local_int2[0][0])))
            ERR(retval)
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        if ((retval = nc_get_vara_double(fileid, data_variable1_id, start, count, &vector_data_local_double1[0][0])))
            ERR(retval)

        if ((retval = nc_get_vara_double(fileid, data_variable2_id, start, count, &vector_data_local_double2[0][0])))
            ERR(retval)
    }
    else
    {
        INMOST_ICE_ERR("var type is not int or double, can't read");
    }

    // get lat lon
    double lon_data_local_double[y_count][x_count];
    double lat_data_local_double[y_count][x_count];

    size_t start_coords[] = {y_start, x_start};
    size_t count_coords[] = {y_count, x_count};

    if ((retval = nc_get_vara_double(fileid, lon_id, start_coords, count_coords, &lon_data_local_double[0][0])))
        ERR(retval)

    if ((retval = nc_get_vara_double(fileid, lat_id, start_coords, count_coords, &lat_data_local_double[0][0])))
        ERR(retval)

    
    // convert velocity components to double
    if (atttype == NC_SHORT)
    {
        
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val_x, current_val_y;
                current_val_x = (vector_data_local_short1[j][i] == invalid_value_short) ? invalid_value_fill 
                                 : (double)(vector_data_local_short1[j][i])*scale_factor + offset_value;
                current_val_y = (vector_data_local_short2[j][i] == invalid_value_short) ? invalid_value_fill
                                 : (double)(vector_data_local_short2[j][i])*scale_factor + offset_value; 
                vector_data_local_double1[j][i] = current_val_x;
                vector_data_local_double2[j][i] = current_val_y;
            }
        }
    }           
    else if ((atttype == NC_INT) or (atttype == NC_LONG))
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                double current_val_x, current_val_y;
                current_val_x = (vector_data_local_int1[j][i] == invalid_value_int) ? invalid_value_fill 
                                 : (double)(vector_data_local_int1[j][i])*scale_factor + offset_value;
                current_val_y = (vector_data_local_int2[j][i] == invalid_value_int) ? invalid_value_fill
                                 : (double)(vector_data_local_int2[j][i])*scale_factor + offset_value;
                vector_data_local_double1[j][i] = current_val_x;
                vector_data_local_double2[j][i] = current_val_y;
            }
        }
    }
    else if ((atttype == NC_FLOAT) or (atttype == NC_DOUBLE))
    {
        for (size_t j = 0; j < y_count; ++j)
        {
            for(size_t i = 0; i < x_count; ++i)
            {
                 double current_val_x, current_val_y;
                current_val_x = (vector_data_local_double1[j][i] == invalid_value_double) ? invalid_value_fill 
                                     : (vector_data_local_double1[j][i])*scale_factor + offset_value;
                current_val_y = (vector_data_local_double2[j][i] == invalid_value_double) ? invalid_value_fill 
                                     : (vector_data_local_double2[j][i])*scale_factor + offset_value;
                vector_data_local_double1[j][i] = current_val_x;
                vector_data_local_double2[j][i] = current_val_y;
            }
        }
    }
    else
    {
        INMOST_ICE_ERR("unknown value type");
    }

    // rotate vector to geo directions using lat, lon info
    double rotated_vector_x[y_count][x_count];
    double rotated_vector_y[y_count][x_count];

    double radian = M_PI/180.0;
    double dlon_ip1[y_count][x_count];
    double dlat_ip1[y_count][x_count];
    double theta_ip1;
    double theta_jp1;
    
    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count - 1; ++i)
        {
            dlon_ip1[j][i] = (lon_data_local_double[j][i+1] - lon_data_local_double[j][i])*
            cos(radian* 0.5*(lat_data_local_double[j][i+1] + lat_data_local_double[j][i]));

            dlat_ip1[j][i] = (lat_data_local_double[j][i+1] - lat_data_local_double[j][i]);
        }
        dlon_ip1[j][x_count - 1] = dlon_ip1[j][x_count - 2];
        dlat_ip1[j][x_count - 1] = dlat_ip1[j][x_count - 2];
    }

    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count; ++i)
        {
            theta_ip1 = atan2(dlat_ip1[j][i], dlon_ip1[j][i]);
            theta_jp1 = theta_ip1 + radian*90.0;

            rotated_vector_x[j][i] = vector_data_local_double1[j][i]*cos(theta_ip1) + 
                                     vector_data_local_double2[j][i]*cos(theta_jp1);
            rotated_vector_y[j][i] = vector_data_local_double1[j][i]*sin(theta_ip1) + 
                                     vector_data_local_double2[j][i]*sin(theta_jp1);
        }
    } 

    // rotate vector to model directions
    for (size_t j = 0; j < y_count; ++j)
    {
        for (size_t i = 0; i < x_count; ++i)
        {
            double vec_x = rotated_vector_x[j][i];
            double vec_y = rotated_vector_y[j][i];
            std::vector<double> model_rot_vec = from_geo_2_model_vec<double>(vec_x, vec_y,
                                                                             lon_data_local_double[j][i],
                                                                             lat_data_local_double[j][i]);
            rotated_vector_x[j][i] = model_rot_vec[0];
            rotated_vector_y[j][i] = model_rot_vec[1];
       }
    }


    
    // get square variables
    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
    {
        double local_x_start = x_coords[x_start];
        double local_y_start = y_coords[y_start];
        double local_x_end = x_coords[x_start + x_count];
        double local_y_end = y_coords[y_start + y_count];

        if(nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(nc.topaz_coords)[0];
            double y = nodeit->RealArray(nc.topaz_coords)[1];

            if ((x < local_x_start) or
                (x > local_x_end) or
                (y < local_y_start) or
                (y > local_y_end))
            {
                nodeit->RealArray(node_data_tags[node_variable_name])[0] = no_extrapolation_fill;
                nodeit->RealArray(node_data_tags[node_variable_name])[1] = no_extrapolation_fill;
                continue;
            }
            else
            {
                size_t x_prev_pos = (size_t)((x - local_x_start)/dx);
                size_t y_prev_pos = (size_t)((y - local_y_start)/dy);
                
                double datax_ld = rotated_vector_x[y_prev_pos][x_prev_pos];
                double datax_lu = rotated_vector_x[y_prev_pos+1][x_prev_pos];
                double datax_rd = rotated_vector_x[y_prev_pos][x_prev_pos+1];
                double datax_ru = rotated_vector_x[y_prev_pos+1][x_prev_pos+1];

                double datay_ld = rotated_vector_y[y_prev_pos][x_prev_pos];
                double datay_lu = rotated_vector_y[y_prev_pos+1][x_prev_pos];
                double datay_rd = rotated_vector_y[y_prev_pos][x_prev_pos+1];
                double datay_ru = rotated_vector_y[y_prev_pos+1][x_prev_pos+1];

                double xl = x_coords[x_start + x_prev_pos];
                double xr = x_coords[x_start + x_prev_pos + 1];
                double yd = y_coords[y_start + y_prev_pos];
                double yu = y_coords[y_start + y_prev_pos + 1];


                double curr_data_x = bilinear_interpolation(xl, xr, 
							                                yd, yu, 
							                                x , y,
							                                datax_ld, datax_lu,
                                                            datax_rd, datax_ru);
                
                double curr_data_y = bilinear_interpolation(xl, xr, 
							                                yd, yu, 
							                                x , y,
							                                datay_ld, datay_lu,
                                                            datay_rd, datay_ru);

                double abs_val = sqrt(curr_data_x*curr_data_x + curr_data_y*curr_data_y);

				nodeit->RealArray(node_data_tags[node_variable_name])[0] = 
                (abs_val > max_abs_value) ? invalid_value_fill : curr_data_x;
                nodeit->RealArray(node_data_tags[node_variable_name])[1] = 
                (abs_val > max_abs_value) ? invalid_value_fill : curr_data_y;
            }
        }
    }
    BARRIER;
    // Exchange data to ghost cells
    ice_mesh->ExchangeData(node_data_tags[node_variable_name], NODE, 0);	
    // Close file
    if ((retval = nc_close(fileid)))
        ERR(retval);

}

//void INMOST_ICE_nodes::InterpolateVectorFromNETCDF(std::string filename,
//                                                   int time_step,
//                                                   std::string node_variable_name,
//                                                   std::string nc_variable_name_first,
//                                                   std::string nc_variable_name_second,
//                                                   INMOST::Tag netcdf_coords
//                                                   )
//{
//
//};
