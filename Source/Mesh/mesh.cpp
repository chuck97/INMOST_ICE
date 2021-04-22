#include "inmost_ice.h"

using namespace INMOST;

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

INMOST_ICE_mesh_data* INMOST_ICE_mesh::GetMeshData()
{
    return &data;
};