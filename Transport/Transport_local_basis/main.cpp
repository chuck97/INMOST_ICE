#include "config.h"
#include "INMOST_ICE_mesh.h"
#include "Initialization.h"
#include "Assembling.h"
#include "Errors.h"
#include "MatVec.h"
#include "FCT.h"
#include "Advection.h"
#include "Function.h"
#include "Numerical.h"

using namespace INMOST;

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


    // calculate transition matrix to local basis 
    INMOST::Tag u_transition_matricies = n.GetMesh()->CreateTag("u transition matricies", DATA_REAL, CELL, NONE, 12);
    INMOST::Tag local_node_coords = n.GetMesh()->CreateTag("local node coords", DATA_REAL, CELL, NONE, 6);

    BARRIER

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

    for (size_t stepn = 0; stepn < AdvTest.GetTotalStepNumber(); ++stepn)
    {
        AdvTest.UpdateVelocity(current_time_hours);
        BARRIER
        AdvSol.AssembleRHS();
        BARRIER
        AdvSol.Evaluate();
        BARRIER
        current_time_hours += AdvTest.GetTimeStepHours();
        BARRIER
    }

}