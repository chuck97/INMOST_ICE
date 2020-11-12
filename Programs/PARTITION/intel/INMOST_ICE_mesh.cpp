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
        ERR("Can't write mesh data to this file. Filename should ended by .pmf")
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
        ERR("Can't write mesh to this file. Filename should ended by .pvtk")
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
        nodeit->RealArray(nc.topaz_coords)[0] = topaz[0];
        nodeit->RealArray(nc.topaz_coords)[1] = topaz[1];
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

void INMOST_ICE_nodes::FindExtrCoords()
{
    double tmp_min_x = ice_mesh->BeginNode()->RealArray(nc.topaz_coords)[0];
    double tmp_max_x = ice_mesh->BeginNode()->RealArray(nc.topaz_coords)[0];
    double tmp_min_y = ice_mesh->BeginNode()->RealArray(nc.topaz_coords)[1];
    double tmp_max_y = ice_mesh->BeginNode()->RealArray(nc.topaz_coords)[1];

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double current_x = nodeit->RealArray(nc.topaz_coords)[0];
            double current_y = nodeit->RealArray(nc.topaz_coords)[1];
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
}

INMOST::Mesh* INMOST_ICE_nodes::GetMesh()
{
    return ice_mesh;
};
