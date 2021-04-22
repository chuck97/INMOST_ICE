#include "config.h"
#include "INMOST_ICE_mesh.h"
#include "Initialization.h"
#include "Assembling.h"
#include "Errors.h"
#include "MatVec.h"
#include "FCT.h"
#include "Advection.h"

using namespace INMOST;

double sum_vec(const std::vector<double>& vec)
{
    double s = 0.0;
    for (auto item: vec)
    {
        s += item;
    }
    return s;
}


int main(int argc, char* argv[])
{
#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
#endif

    // read json
    std::string current_exec_name = argv[0]; 
    std::vector<std::string> all_args;

    std::string json_path;

    if (argc > 1) 
    {
        all_args.assign(argv + 1, argv + argc);
    }

    if (all_args.size() > 1)
    {
        INMOST_ICE_ERR("should be only *.json file on input");
    }
    else
    {
        json_path = all_args.back();
    }

    // all input parameters
    std::string mesh_path;
    double total_time_hours;
    InitialMassParams init_mass_params;
    VelocityParams velocity_params;
    AdvectionSolverParams AdvSolParams;
    OutputParameters OutpParams;

    ParseJson(json_path,
              mesh_path,
              total_time_hours,
              init_mass_params,
              velocity_params,
              AdvSolParams,
              OutpParams
              );

    // mesh initialization
    INMOST_ICE_mesh m(mesh_path);

    // nodes initialization
    INMOST_ICE_nodes n(m.GetMesh());

    // Create mass and velocity tags
    INMOST::Tag m_tag = n.GetMesh()->CreateTag("m", DATA_REAL, NODE, NONE, 1);
    INMOST::Tag init_m_tag = n.GetMesh()->CreateTag("m_init", DATA_REAL, NODE, NONE, 1);
    INMOST::Tag u_tag = n.GetMesh()->CreateTag("u", DATA_REAL, NODE, NONE, 3);

    // velocity parameters

    // setup advection test    
    AdvectionTest AdvTest(n,
                          init_mass_params,
                          velocity_params,
                          m_tag,
                          init_m_tag,
                          u_tag,
                          total_time_hours);

    // Fill AdvSolParams
    AdvSolParams.time_step_sec = AdvTest.GetTimeStepSeconds();
    AdvSolParams.total_step_number = AdvTest.GetTotalStepNumber(); 

    // Assign initial mass distribution
    AdvTest.AssignInitialMass();
    
    // setup solver of linear systems
    Solver linear_solver("inner_ilu2"); 
    linear_solver.SetParameter("absolute_tolerance", "1e-9"); 
    
    BARRIER

    // setup advection solver
    AdvectionSolver AdvSol(n,
                           AdvSolParams,
                           m_tag,
                           init_m_tag,
                           u_tag,
                           linear_solver,
                           AdvTest.GetTotalStepNumber(),
                           OutpParams
                           );
    BARRIER

    // Assemble LHS
    AdvSol.AssembleLHS();

    // Time stepping
    double current_time_hours = 0.0;

    // Times
    std::vector<double> RHS_profiler;
    std::vector<double> SOL_profiler;

    for (size_t stepn = 0; stepn < AdvTest.GetTotalStepNumber(); ++stepn)
    {
        AdvTest.UpdateVelocity(current_time_hours);
        BARRIER
        double ttt = Timer();
        AdvSol.AssembleRHS();
        RHS_profiler.push_back(Timer() - ttt);
        BARRIER
        ttt = Timer();
        AdvSol.Evaluate();
        SOL_profiler.push_back( Timer() - ttt);
        BARRIER
        current_time_hours += AdvTest.GetTimeStepHours();
        BARRIER
    }
    double sum_RHS_time = sum_vec(RHS_profiler);
    double sum_SOL_time = sum_vec(SOL_profiler);
    double max_RHS_time;
    double max_SOL_time;
    BARRIER
#if defined(USE_MPI)
   MPI_Allreduce(&sum_RHS_time, &max_RHS_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

#if defined(USE_MPI)
   MPI_Allreduce(&sum_SOL_time, &max_SOL_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    BARRIER

    if (n.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "avg RHS assembling time: " << sum_RHS_time/AdvTest.GetTotalStepNumber() << std::endl;
        std::cout << "avg SOL time: " << sum_SOL_time/AdvTest.GetTotalStepNumber() << std::endl;
    }

}