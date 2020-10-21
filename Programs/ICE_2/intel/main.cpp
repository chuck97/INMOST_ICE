#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include "inmost.h"
#define USE_MPI
#define USE_PARTITIONER
using namespace INMOST;
#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif


#define USE_PARTITIONER_PARMETIS


int main(int argc, char **argv) 
{
    int i, j, k;
    int processRank = 0, processorsCount = 1;

#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);     // Get the rank of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &processorsCount); // Get the total number of processors used
#endif

    //get parameters from arguments
    std::string parametersFileName = "";
    std::string solverName = "inner_ilu2";
    std::string MeshName = "";

    bool parametersFound = false;
    bool typeFound = false;
    bool meshFound = false;
    bool refine = false;

    //Parse argv parameters
    if (argc == 1) goto helpMessage;
    for (i = 1; i < argc; i++) 
    {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) 
        {
            helpMessage:
            if (processRank == 0) 
            {
                std::cout << "Help message: " << std::endl;
                std::cout << "Command line options: " << std::endl;
                std::cout << "-d, --database <Solver parameters file name>" << std::endl;
                std::cout << "-t, --type <Solver type name>" << std::endl;
                std::cout << "-m, --mesh <Mesh name>" << std::endl;
                std::cout << "  Available solvers:" << std::endl;
                Solver::Initialize(NULL, NULL, NULL);
                std::vector<std::string> availableSolvers = Solver::getAvailableSolvers();
                for (solvers_names_iterator_t it = availableSolvers.begin(); it != availableSolvers.end(); it++) 
                {
                    std::cout << "      " << *it << std::endl;
                }
                Solver::Finalize();
            }
            return 0;
        }
        //Parameters file name found with -d or --database options
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--database") == 0) 
        {
            if (processRank == 0) 
            {
                std::cout << "Solver parameters file found: " << argv[i + 1] << std::endl;
            }
            parametersFound = true;
            parametersFileName = std::string(argv[i + 1]);
            i++;
            continue;
        }
        //Solver type found with -t ot --type options
        if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0) 
        {
            if (processRank == 0) 
            {
                std::cout << "Solver type index found: " << argv[i + 1] << std::endl;
            }
            typeFound = true;
            solverName = std::string(argv[i + 1]);
            i++;
            continue;
        }

        if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--mesh") == 0) 
        {
            if (processRank == 0) 
            {
                std::cout << "mesh name found: " << argv[i + 1] << std::endl;
            }
            meshFound = true;
            MeshName = std::string(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "--refine") == 0) 
        {
            if (processRank == 0) 
            {
                std::cout << "mesh will be refined " << argv[i + 1] << " times" << std::endl;
            }
            refine = true;
            i++;
            continue;
        }
    }
    if (!typeFound)
        if (processRank == 0)
            std::cout
                    << "Solver type not found in command line, you can specify solver type with -t or --type option, using INNER_ILU2 solver by default."
                    << std::endl;

    if (!parametersFound)
        if (processRank == 0)
            std::cout << "Parameters not found, you can specify parameter file name with -d or --database option, "
                      << std::endl;

    if (!meshFound) {
        if (processRank == 0)
            std::cout
                    << "Mesh name not found , you can specify mesh name  with -am or --animesh option."
                    << std::endl;
        exit(-1);
    }

    // prepare  INMOST 
    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);
    Solver::Initialize(&argc, &argv, parametersFound ? parametersFileName.c_str() : NULL);

    if (!Solver::isSolverAvailable(solverName)) {
        if (processRank == 0) {
            std::cout << "Solver " << solverName << " is not available" << std::endl;
        }
        Mesh::Finalize();
        Partitioner::Finalize();
        Solver::Finalize();
        exit(-1);
    }


    if (processRank == 0) 
    {
        std::cout << "Solving with " << solverName << std::endl;
    }
    double t0, t_init, t_refine, t_prepare, t_begin_it, t_prepare_it, t_assmble, t_solve, t_exch;

    double sum_t_prepare, sum_t_assmble, sum_t_solve, sum_t_exch;
    BARRIER
    t0 = Timer();
    // load and repartition mesh 
    Mesh* mesh_init;
    mesh_init = new Mesh();
#if defined(USE_MPI)
    mesh_init->SetCommunicator(INMOST_MPI_COMM_WORLD);
#endif
    mesh_init->Load(MeshName);

#ifdef USE_PARTITIONER
    Partitioner* p = new Partitioner(mesh_init);
#ifdef USE_PARTITIONER_PARMETIS
    p->SetMethod(Partitioner::Parmetis, Partitioner::Partition);
#elif USE_PARTITIONER_ZOLTAN
    p->SetMethod(Partitioner::Zoltan, Partitioner::Partition);
#else
    p->SetMethod(Partitioner::Inner_RCM, Partitioner::Partition);
#endif
    p->Evaluate();
    delete p;
    mesh_init->Redistribute();
    mesh_init->ReorderEmpty(CELL | FACE | EDGE | NODE);
    mesh_init->ExchangeGhost(1, NODE);
#endif

    if (processRank == 0) 
    {
        std::cout << " repartitioning success" << std::endl;
    }
    BARRIER

    t_init = Timer();
    Mesh *mesh_main;

    mesh_main = mesh_init;
    
    BARRIER

    t_refine = Timer();

    //set Boundary conditions
    mesh_main->AssignGlobalID(NODE);
    long int min_index, max_index;

    long int tmp_var = mesh_main->TotalNumberOf(NODE|EDGE|FACE|CELL)+1;

    min_index = tmp_var;
    max_index = -1;
    long int tot_nodes = mesh_main->TotalNumberOf(NODE);
    if (processRank == 0)
    {
        std::cout<<"num nodes " << tot_nodes;
    }

    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++) 
    {
        tmp_var = it->GlobalID();
        if(it->GetStatus() != Element::Ghost)
        {
            if(tmp_var < min_index)
            {
                min_index = tmp_var;
            }

            if(tmp_var > max_index)
            {
                max_index = tmp_var;
            }
        }
    }

    std::cout<<processRank << ": min_index "<<min_index<<" max_index "<<max_index<<std::endl;

    double min_x = mesh_main->BeginNode()->Coords()[0];
    double max_x = mesh_main->BeginNode()->Coords()[1];
    double min_y = mesh_main->BeginNode()->Coords()[0];
    double max_y = mesh_main->BeginNode()->Coords()[1];
    for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++) 
    {
        if (it->Coords()[0] < min_x)
        {
            min_x = it->Coords()[0];
        }
        if (it->Coords()[0] > max_x)
        {
            max_x = it->Coords()[0];
        }
        if (it->Coords()[1] < min_y)
        {
            min_y = it->Coords()[1];
        }
        if (it->Coords()[1] > max_y)
        {
            max_y = it->Coords()[1];
        }
    }
    std::cout << processRank << ": min_x = " << min_x << " max_x = " << max_x << "; min_y = " << min_y << " max_y = " << max_y << std::endl; 
    BARRIER
    Tag Sol_tag = mesh_main->CreateTag("U", DATA_REAL, NODE, NONE, 1);

    double DeltaT = 1;
    int Itime = 0;
    double Time = 0;
    double FinalTime = 10;
    do 
    {
        BARRIER
        t_begin_it = Timer();
        Solver solv = Solver(solverName, "test");

        Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
        Sparse::Vector b("rhs"); // Declare the right-hand side vector
        Sparse::Vector x("sol"); // Declare the solution vector
        x.SetInterval(min_index, max_index + 1);
        b.SetInterval(min_index, max_index + 1);
        mat.SetInterval(min_index, max_index + 1);

        
        Itime++;
        solv.SetMatrix(mat);

        //    BARRIER
        if (processRank == 0) 
        {
            std::cout << "solving" << std::endl;
        }

        if (!solv.Solve(b, x)) 
        {
            std::cout << "solution failed, residual  " << solv.Residual() << " iterations " << solv.Iterations()
                      << " name "
                      << std::endl;
            return 0;
        }

        for (Mesh::iteratorNode it = mesh_main->BeginNode(); it != mesh_main->EndNode(); it++)
        {
            if (it->GetStatus() != Element::Ghost) 
            {
                it->Real(Sol_tag) = x[it->GlobalID()];
            }
        }
        mesh_main->ExchangeData(Sol_tag, NODE, 0);

        solv.Clear();
        Time += DeltaT;
    } while (Time < FinalTime);



    if (processRank == 0) 
    {
        std::cout << std::endl;
        std::cout << "End of program"  << std::endl;
    }

    mesh_main->Save("save_sol.pvtk");

    return 0;
}
