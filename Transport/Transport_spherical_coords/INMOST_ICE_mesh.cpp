#include "INMOST_ICE_mesh.h"
using namespace INMOST;

// trim from start (in place)
void ltrim(std::string &s) 
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), 
    [](unsigned char ch) 
    {
        return !std::isspace(ch);
    }));
};

// trim from end (in place)
void rtrim(std::string &s) 
{
    s.erase(std::find_if(s.rbegin(), s.rend(), 
    [](unsigned char ch) 
    {
        return !std::isspace(ch);
    }).base(), s.end());
};

// trim from both ends (in place)
void trim(std::string &s) 
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

    bool repartition = false;
    ice_mesh = new Mesh();

#if defined(USE_MPI)
    ice_mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif

    if(ice_mesh->isParallelFileFormat(mesh_path))
    {
        ice_mesh->Load(mesh_path); 
        repartition = true;

        if (ice_mesh->GetProcessorRank()==0)
        {
            std::cout << "+++++++++++++ Parallel realization +++++++++++++" << std::endl;
        }
    }
    else
    {
        if (ice_mesh->GetProcessorRank() == 0)
        {
            ice_mesh->Load(mesh_path);
            std::cout << "+++++++++++++ Serial realization +++++++++++++" << std::endl;
        }
    }

    if (ice_mesh->GetProcessorRank()==0)
    {
        std::cout << "??????????? Mesh " << mesh_path 
                  <<  " loaded successfully ???????????" 
                  << std::endl;

        //std::cout << "There are " 
        //          << ice_mesh->GetProcessorsNumber() 
        //          <<" processes;" <<  std::endl;
    }
    BARRIER
    
    Partition();
};

void INMOST_ICE_mesh::Partition()
{
#if defined(USE_PARTITIONER)
	if (ice_mesh->GetProcessorsNumber() > 1) // need repartition
	{ 
	    Partitioner* p = new Partitioner(ice_mesh);
#ifdef USE_PARTITIONER_PARMETIS
        p->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
#elif USE_PARTITIONER_ZOLTAN
        p->SetMethod(Partitioner::Zoltan, Partitioner::Partition);
#else
        p->SetMethod(Partitioner::INNER_KMEANS, Partitioner::Partition);
#endif
        BARRIER
		p->Evaluate();
		delete p;
		BARRIER

		ice_mesh->Redistribute();
		ice_mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
        ice_mesh->ExchangeGhost(1, NODE); 
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
	ice_mesh->Save(no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << ";" << std::endl; 
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
	ice_mesh->Save(no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << ";" << std::endl; 
    }
};

void INMOST_ICE_mesh::PrintPVTU(const std::string& filename)
{
    std::string no_spaces_filename = filename;
    rtrim(no_spaces_filename);
    if(no_spaces_filename.substr(no_spaces_filename.size()-5, 5) != ".pvtu")
    {
        INMOST_ICE_ERR("Can't write mesh to this file. Filename should ended by .pvtu")
    }
    BARRIER
	ice_mesh->Save(no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << ";" << std::endl; 
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

    nc.model_coords = ice_mesh->CreateTag("model coords", DATA_REAL, NODE, NONE, 3);
    nc.geo_coords = ice_mesh->CreateTag("geo coords", DATA_REAL, NODE, NONE, 3);

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        // Assign model coords

        double x = nodeit->Coords()[0];
        double y = nodeit->Coords()[1];
        double z = nodeit->Coords()[2];

		nodeit->RealArray(nc.model_coords)[0] = x;
        nodeit->RealArray(nc.model_coords)[1] = y;
        nodeit->RealArray(nc.model_coords)[2] = z;

        double r = sqrt(x*x + y*y + z*z);
        double lat, lon;

    
        if ((x >= 0) and (y >= 0))
        {
            lat = std::asin(z/r);
            lon = std::atan(y/x);
        }
        else if ((x >= 0) and (y <= 0))
        {
            lat = std::asin(z/r);
            lon = 2*M_PI + std::atan(y/x);
        }
        else if ((x <= 0) and (y >= 0))
        {
            lat = std::asin(z/r);
            lon = M_PI + std::atan(y/x);
        }
        else
        {
            lat = std::asin(z/r);
            lon = M_PI + std::atan(y/x);
        }
        

        // Assign geo coords
        nodeit->RealArray(nc.geo_coords)[0] = lon;
        nodeit->RealArray(nc.geo_coords)[1] = lat;
        nodeit->RealArray(nc.geo_coords)[2] = EARTH_RADIUS;
    }
    
    // calculate global interval
    INMOST::Tag id = ice_mesh->GlobalIDTag();
    idmin = std::numeric_limits<unsigned int>::max();
    idmax = std::numeric_limits<unsigned int>::min();

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode();
        nodeit != ice_mesh->EndNode();
        ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            unsigned int pid = nodeit->Integer(id);
            if(pid < idmin)
            {
                idmin = pid;
            } 
            if((pid + 1) > idmax)
            {
                idmax = pid + 1;
            } 
        }
    }

    // calculate spatial resolution
    double proc_sp_resolution = std::numeric_limits<double>::max();

    for(Mesh::iteratorCell trianit = ice_mesh->BeginCell();
        trianit != ice_mesh->EndCell();
        ++trianit) 
    {
        INMOST::ElementArray<Node> trian = trianit->getNodes();

        double x0 = trian[0]->RealArray(nc.model_coords)[0];
        double x1 = trian[1]->RealArray(nc.model_coords)[0];

        double y0 = trian[0]->RealArray(nc.model_coords)[1];
        double y1 = trian[1]->RealArray(nc.model_coords)[1];

        double z0 = trian[0]->RealArray(nc.model_coords)[2];
        double z1 = trian[1]->RealArray(nc.model_coords)[2];


        double curr_res = EARTH_RADIUS*sqrt((x0 - x1)*(x0 - x1) +
                                            (y0 - y1)*(y0 - y1) +
                                            (z0 - z1)*(z0 - z1));
        if (curr_res < proc_sp_resolution)
        {
            proc_sp_resolution = curr_res;
        }
    }   
    BARRIER
    
    double sp_resolution = proc_sp_resolution;

#if defined(USE_MPI)
   MPI_Allreduce(&proc_sp_resolution, &sp_resolution, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
    BARRIER
    spatial_resolution = sp_resolution;
};

NodeCoords& INMOST_ICE_nodes::GetCoords()
{
    return nc;
};

INMOST::Mesh* INMOST_ICE_nodes::GetMesh()
{
    return ice_mesh;
};

std::map<std::string, INMOST::Tag>& INMOST_ICE_nodes::GetData()
{
    return node_data_tags;
}

unsigned int INMOST_ICE_nodes::GetIdMin() const
{
    return idmin;
}

unsigned int INMOST_ICE_nodes::GetIdMax() const
{
    return idmax;
}

double INMOST_ICE_nodes::GetResolution() const
{
    return spatial_resolution;
}

void INMOST_ICE_nodes::PrintPVTU(const std::string& filename)
{
    std::string no_spaces_filename = filename;
    rtrim(no_spaces_filename);
    if(no_spaces_filename.substr(no_spaces_filename.size()-5, 5) != ".pvtu")
    {
        INMOST_ICE_ERR("Can't write mesh to this file. Filename should ended by .pvtu")
    }
    BARRIER
	ice_mesh->Save(no_spaces_filename);
	BARRIER
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Mesh saved to " << no_spaces_filename << ";" << std::endl; 
    }
};