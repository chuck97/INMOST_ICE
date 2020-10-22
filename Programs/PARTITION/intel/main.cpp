#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include "inmost.h"

#define USE_MPI
#define USE_PARTITIONER
#define USE_PARTITIONER_PARMETIS

using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif


#define M_Assert(Expr, Msg) \
    if (!Expr)  \
    { \
    std::cerr << "Assert failed:\t" << Msg << "\n" \
        << "Expected:\t" << #Expr << "\n" \
        << "Source:\t\t" << __FILE__ << ", line " << __LINE__ << "\n"; \
        abort(); \
    }


int main(int argc, char **argv)
{
    Mesh::Initialize(&argc, &argv);
#if defined(USE_PARTITIONER)
    Partitioner::Initialize(&argc, &argv);
#endif

    double ttt = Timer();

    bool repartition = false;
    Mesh* m;
    m = new Mesh();

#if defined(USE_MPI)
    m -> SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif

    std::string MeshName = "/data90t/geosci/spetrov/INMOST_TESTS/GRID_TO_PMF/test.pmf";
    ttt = Timer();

    if(m -> isParallelFileFormat(MeshName))
    {
        m -> Load(MeshName); 
        repartition = true;

        if (m -> GetProcessorRank() == 0)
        {
            std::cout << "Parallel realization" << std::endl;
        }
    }
    else
    {
        if (m -> GetProcessorRank() == 0)
        {
            m -> Load(MeshName);
            std::cout << "Serial realization" << std::endl;
        }
    }


    BARRIER

    if (m -> GetProcessorRank() == 0)
    {
        std::cout << "Mesh loaded successfully: "<<  Timer() - ttt << std::endl;
        std::cout << "Processors: " << m -> GetProcessorsNumber() << std::endl;
    }

    
    Tag node_test_data = m->CreateTag("node test data", DATA_REAL, NODE, NONE, 1); 
    
#if defined(USE_PARTITIONER)
	if (m -> GetProcessorsNumber() > 1) // need repartition
	{ 
        ttt = Timer();
	    Partitioner* p = new Partitioner(m);
#ifdef USE_PARTITIONER_PARMETIS
        p -> SetMethod(Partitioner::Parmetis, Partitioner::Partition);
#elif USE_PARTITIONER_ZOLTAN
        p -> SetMethod(Partitioner::Zoltan, Partitioner::Partition);
#else
        p -> SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition);
#endif
		p -> Evaluate(); // Compute the partitioner and store new processor ID in the mesh
		delete p;
		BARRIER

	    if( m -> GetProcessorRank() == 0 )
        {
            std::cout << "Evaluate: " << Timer()-ttt << std::endl;
        } 

        BARRIER
		ttt = Timer();
		m -> Redistribute(); // Redistribute the mesh data
		m -> ReorderEmpty(CELL|FACE|EDGE|NODE); // Clean the data after reordring
        m -> ExchangeGhost(1, NODE);
		BARRIER

		if( m -> GetProcessorRank() == 0 )
        {
        	std::cout << "Redistribute: " << Timer()-ttt << std::endl;
        } 
        BARRIER
	}
#endif

    ttt = Timer();
    //assigning global id
	m -> AssignGlobalID(CELL|EDGE|FACE|NODE);
    
    //assigning local id
    int i = 0;
    for (Mesh::iteratorNode it = m->BeginNode(); it != m->EndNode(); it++)
    {
        it->Real(node_test_data) = 1.0; 
    }

    BARRIER
	if( m -> GetProcessorRank() == 0 )
    {
        std::cout << "Assign ids: " << Timer()-ttt << std::endl;
    } 


    ttt = Timer();
	m -> Save("test.pvtk");
	BARRIER

	if( m->GetProcessorRank() == 0 )
    {
        std::cout << "Save mesh:" << Timer() - ttt << std::endl;
    } 

	delete m;
    
#if defined(USE_PARTITIONER)
	Partitioner::Finalize(); // Finalize the partitioner activity
#endif

    Mesh::Finalize();
	return 0;
}