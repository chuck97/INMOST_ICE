#include "Advection.h"

using namespace INMOST;

AdvectionSolver::AdvectionSolver(
                    INMOST_ICE_nodes& n_,
                    const AdvectionSolverParams& params_,
                    INMOST::Tag m_tag_,
                    INMOST::Tag m_init_tag_,
                    INMOST::Tag u_tag_,
                    INMOST::Solver& Sol_,
                    size_t total_step_num_,
                    OutputParameters& output_params_)
    : n(n_),
      Sol(Sol_),
      total_step_num(total_step_num_),
      output_params(output_params_)

{
    params = params_;
    m_tag = m_tag_;
    m_init_tag = m_init_tag_;
    u_tag = u_tag_;
    int idmin = n.GetIdMin();
    int idmax = n.GetIdMax();
    BARRIER
    LHS.SetInterval(idmin, idmax);
    LHS_low.SetInterval(idmin, idmax);
    RHS.SetInterval(idmin, idmax);
    RHS_low.SetInterval(idmin, idmax);
    init_mass_integral = IntegralMass(n, m_init_tag);
    total_time = Timer();
    output_frequency = total_step_num/output_params.number_of_screenshots;
    LogInit();
}

void AdvectionSolver::AssembleLHS()
{
    INMOST::Tag id = n.GetMesh()->GlobalIDTag();
    INMOST::Tag M_L_entire = n.GetMesh()->CreateTag("M_L_entire", DATA_REAL, NODE, NONE, 1);

    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
            trianit != n.GetMesh()->EndCell();
            ++trianit) 
    {
        // get node coords of trian
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        double v0_x = trianit->RealArray(n.GetLocalNodeCoords())[0];
        double v0_y = trianit->RealArray(n.GetLocalNodeCoords())[1];
        double v1_x = trianit->RealArray(n.GetLocalNodeCoords())[2];
        double v1_y = trianit->RealArray(n.GetLocalNodeCoords())[3];
        double v2_x = trianit->RealArray(n.GetLocalNodeCoords())[4];
        double v2_y = trianit->RealArray(n.GetLocalNodeCoords())[5];
        std::vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                              {v1_x, v1_y},
                                                              {v2_x, v2_y}};
        // find global ids of local nodes
        std::vector<int> nums_global = { local_nodes[0].Integer(id),
                                         local_nodes[1].Integer(id),
                                         local_nodes[2].Integer(id) }; 
        // assemble local 3x3 stiffness matrix
        std::vector<std::vector<double>> localLHS = LocalStiffnessMatrixAssembling(local_node_coords);
        if (params.is_fct)
        {
            if (trianit->GetStatus()!= Element::Ghost)
            {
                M_C_minus_M_L.push_back(localLHS);
            } 
        }
        // assemble global LHS mass matrix
        for (int i = 0; i < 3; ++i)
        {
            if(local_nodes[i].GetStatus() != Element::Ghost)
            {
                for (int j = 0; j < 3; ++j)
                {
                    LHS[nums_global[i]][nums_global[j]] += localLHS[i][j];
                }
            }
        }
        if (params.is_fct)
        {
            for (unsigned int i = 0; i < 3; ++i)
            {
                if(local_nodes[i].GetStatus() != Element::Ghost)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        LHS_low[nums_global[i]][nums_global[i]] += localLHS[i][j];
                        local_nodes[i]->Real(M_L_entire) += localLHS[i][j];
                    }
                }
            }
            if (trianit->GetStatus() != Element::Ghost)
            {
                for (unsigned int i = 0; i < 3; ++i)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        M_C_minus_M_L.back()[i][i] -= localLHS[i][j];
                    }
                }
            }
        } 
    }   
    BARRIER
    if (params.is_fct)
    {
        n.GetMesh()->ExchangeData(M_L_entire, NODE, 0);
        BARRIER
    }

    // multiply MC_minus_ML by ML^-1 for FCT
    if (params.is_fct)
    {
        int t_num = 0;
        for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
                trianit != n.GetMesh()->EndCell();
                ++trianit) 
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                INMOST::ElementArray<Node> local_nodes = trianit->getNodes(); 
                for (unsigned int i = 0; i < 3; ++i)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        M_C_minus_M_L[t_num][i][j] /= - local_nodes[i]->Real(M_L_entire);
                    }
                }
                ++t_num;
            }
        }
    }
    BARRIER
    n.GetMesh()->DeleteTag(M_L_entire, NODE);
}

void AdvectionSolver::AssembleSingleStepRHS()
{
    INMOST::Tag id = n.GetMesh()->GlobalIDTag();
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
        trianit != n.GetMesh()->EndCell();
        ++trianit) 
    {
        // get node coords of trian
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        double v0_x = trianit->RealArray(n.GetLocalNodeCoords())[0];
        double v0_y = trianit->RealArray(n.GetLocalNodeCoords())[1];
        double v1_x = trianit->RealArray(n.GetLocalNodeCoords())[2];
        double v1_y = trianit->RealArray(n.GetLocalNodeCoords())[3];
        double v2_x = trianit->RealArray(n.GetLocalNodeCoords())[4];
        double v2_y = trianit->RealArray(n.GetLocalNodeCoords())[5];
        std::vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                              {v1_x, v1_y},
                                                              {v2_x, v2_y}};
        // get u values of trian
        std::vector<std::vector<double>> local_u_geo = 
        {{local_nodes[0].RealArray(u_tag)[0], local_nodes[0].RealArray(u_tag)[1]},
         {local_nodes[1].RealArray(u_tag)[0], local_nodes[1].RealArray(u_tag)[1]},
         {local_nodes[2].RealArray(u_tag)[0], local_nodes[2].RealArray(u_tag)[1]}};
        std::vector<double> u0_local = 
        {
            trianit->RealArray(n.GetTransitionMatricies())[0]*local_u_geo[0][0] +
            trianit->RealArray(n.GetTransitionMatricies())[1]*local_u_geo[0][1],
            trianit->RealArray(n.GetTransitionMatricies())[2]*local_u_geo[0][0] +
            trianit->RealArray(n.GetTransitionMatricies())[3]*local_u_geo[0][1]
        };
        std::vector<double> u1_local = 
        {
            trianit->RealArray(n.GetTransitionMatricies())[4]*local_u_geo[1][0] +
            trianit->RealArray(n.GetTransitionMatricies())[5]*local_u_geo[1][1],
            trianit->RealArray(n.GetTransitionMatricies())[6]*local_u_geo[1][0] +
            trianit->RealArray(n.GetTransitionMatricies())[7]*local_u_geo[1][1]
        };
        std::vector<double> u2_local = 
        {
            trianit->RealArray(n.GetTransitionMatricies())[8]*local_u_geo[2][0] +
            trianit->RealArray(n.GetTransitionMatricies())[9]*local_u_geo[2][1],
            trianit->RealArray(n.GetTransitionMatricies())[10]*local_u_geo[2][0] +
            trianit->RealArray(n.GetTransitionMatricies())[11]*local_u_geo[2][1]
        };
        std::vector<std::vector<double>> local_u = {u0_local, u1_local, u2_local};
        // get m values of trian
        std::vector<double> local_m = 
        {local_nodes[0].Real(m_tag),
         local_nodes[1].Real(m_tag),
         local_nodes[2].Real(m_tag)};
        std::vector<double> localRHS;
        // assemble local load vector of size 3
        if (params.AdvSolType == AdvectionSolverType::TG2)
        {
            localRHS = LocalTG2RhsAssembling(local_node_coords,
                                             local_u,
                                             local_m,
                                             params.time_step_sec);
        }
        else if (params.AdvSolType == AdvectionSolverType::CG2)
        {
            localRHS = LocalCG2RhsAssembling(local_node_coords,
                                             local_u,
                                             local_m,
                                             params.time_step_sec);
        }
        else
        {
            INMOST_ICE_ERR("Only avalible single step solvers: TG2, CG2");
        }
        // assemble global load vector
        for (unsigned int i = 0; i < 3; ++i)
        {
            if(local_nodes[i].GetStatus() != Element::Ghost)
            {
                RHS[local_nodes[i].Integer(id)] += localRHS[i];
            }
        }
    }   
}   

void AdvectionSolver::AssembleDoubleStepRHS(StepNumber step_num)
{
    INMOST::Tag id = n.GetMesh()->GlobalIDTag();
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
        trianit != n.GetMesh()->EndCell();
        ++trianit) 
    {
        // get node coords of trian
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        double v0_x = trianit->RealArray(n.GetLocalNodeCoords())[0];
        double v0_y = trianit->RealArray(n.GetLocalNodeCoords())[1];
        double v1_x = trianit->RealArray(n.GetLocalNodeCoords())[2];
        double v1_y = trianit->RealArray(n.GetLocalNodeCoords())[3];
        double v2_x = trianit->RealArray(n.GetLocalNodeCoords())[4];
        double v2_y = trianit->RealArray(n.GetLocalNodeCoords())[5];
        std::vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                              {v1_x, v1_y},
                                                              {v2_x, v2_y}};


        // get u values of trian
        std::vector<std::vector<double>> local_u_geo = 
        {{local_nodes[0].RealArray(u_tag)[0], local_nodes[0].RealArray(u_tag)[1]},
         {local_nodes[1].RealArray(u_tag)[0], local_nodes[1].RealArray(u_tag)[1]},
         {local_nodes[2].RealArray(u_tag)[0], local_nodes[2].RealArray(u_tag)[1]}};

        std::vector<double> u0_local = 
        {
            trianit->RealArray(n.GetTransitionMatricies())[0]*local_u_geo[0][0] +
            trianit->RealArray(n.GetTransitionMatricies())[1]*local_u_geo[0][1],
            trianit->RealArray(n.GetTransitionMatricies())[2]*local_u_geo[0][0] +
            trianit->RealArray(n.GetTransitionMatricies())[3]*local_u_geo[0][1]
        };

        std::vector<double> u1_local = 
        {
            trianit->RealArray(n.GetTransitionMatricies())[4]*local_u_geo[1][0] +
            trianit->RealArray(n.GetTransitionMatricies())[5]*local_u_geo[1][1],
            trianit->RealArray(n.GetTransitionMatricies())[6]*local_u_geo[1][0] +
            trianit->RealArray(n.GetTransitionMatricies())[7]*local_u_geo[1][1]
        };

        std::vector<double> u2_local = 
        {
            trianit->RealArray(n.GetTransitionMatricies())[8]*local_u_geo[2][0] +
            trianit->RealArray(n.GetTransitionMatricies())[9]*local_u_geo[2][1],
            trianit->RealArray(n.GetTransitionMatricies())[10]*local_u_geo[2][0] +
            trianit->RealArray(n.GetTransitionMatricies())[11]*local_u_geo[2][1]
        };

        std::vector<std::vector<double>> local_u = {u0_local, u1_local, u2_local};

        // get m values of trian
        std::vector<double> local_m = {local_nodes[0].Real(m_tag),
                                       local_nodes[1].Real(m_tag),
                                       local_nodes[2].Real(m_tag)};

        std::vector<double> local_m_half;
        if (step_num == StepNumber::second)
        {
            local_m_half = 
            {local_nodes[0].Real(m_half_tag),
             local_nodes[1].Real(m_half_tag),
             local_nodes[2].Real(m_half_tag)};
        }

        std::vector<double> localRHS;
        
        // assemble local load vector of size 3
        if (params.AdvSolType == AdvectionSolverType::TTG2)
        {
            if (step_num == StepNumber::first)
            {
                localRHS = LocalTTG2RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m,
                                                  params.time_step_sec,
                                                  1);
            }
            else
            {
                localRHS = LocalTTG2RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m_half,
                                                  params.time_step_sec,
                                                  2);
            }
        }
        else if (params.AdvSolType == AdvectionSolverType::TTG3)
        {
            if (step_num == StepNumber::first)
            {
                localRHS = LocalTTG3RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m,
                                                  params.time_step_sec,
                                                  1);
            }
            else
            {
                localRHS = LocalTTG3RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m_half,
                                                  params.time_step_sec,
                                                  2);
            }
        }
        else if (params.AdvSolType == AdvectionSolverType::TTG4)
        {
            if (step_num == StepNumber::first)
            {
                localRHS = LocalTTG4RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m,
                                                  params.time_step_sec,
                                                  1);
            }
            else
            {
                localRHS = LocalTTG4RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m_half,
                                                  params.time_step_sec,
                                                  2);
            }
        }
        else
        {
            INMOST_ICE_ERR("Only avalible double step solvers: TTG2, TTG3, TTG4");
        }
            
        // assemble global load vector
        for (unsigned int i = 0; i < 3; ++i)
        {
            if(local_nodes[i].GetStatus() != Element::Ghost)
            {
                RHS[local_nodes[i].Integer(id)] += localRHS[i];
            }
        }
    }
    BARRIER
}

void AdvectionSolver::AssembleRHS()
{
    if ((params.AdvSolType == AdvectionSolverType::CG2) or 
        (params.AdvSolType == AdvectionSolverType::TG2))
    {
        AssembleSingleStepRHS();
    }
    else if ((params.AdvSolType == AdvectionSolverType::TTG2) or 
             (params.AdvSolType == AdvectionSolverType::TTG3) or
             (params.AdvSolType == AdvectionSolverType::TTG4))
    {
        AssembleDoubleStepRHS(StepNumber::first);
    }
    else
    {
        INMOST_ICE_ERR("unknown type of solver: caan't assemble RHS");
    }
}

void AdvectionSolver::Evaluate()
{
    m_high_tag = n.GetMesh()->CreateTag("m_high", DATA_REAL, NODE, NONE, 1);
    if (params.is_fct)
    {
        m_low_tag = n.GetMesh()->CreateTag("m_low", DATA_REAL, NODE, NONE, 1);
    }
    if ((params.AdvSolType != AdvectionSolverType::TG2) and (params.AdvSolType != AdvectionSolverType::CG2))
    {
        m_half_tag = n.GetMesh()->CreateTag("m_half", DATA_REAL, NODE, NONE, 1);
    }


    INMOST::Tag id = n.GetMesh()->GlobalIDTag();
    INMOST::Sparse::Vector SOL;
    unsigned int idmin = n.GetIdMin();
    unsigned int idmax = n.GetIdMax();
    SOL.SetInterval(idmin, idmax);
    Sol.SetMatrix(LHS);
    BARRIER
    Sol.Solve(RHS, SOL);
    // update m high value
    for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            if ((params.AdvSolType != AdvectionSolverType::TG2) and
                (params.AdvSolType != AdvectionSolverType::CG2))
            {
                nodeit->Real(m_half_tag) = SOL[nodeit->Integer(id)];
            }
            else
            {
                nodeit->Real(m_high_tag) = SOL[nodeit->Integer(id)];
            }
        }
    }
    BARRIER
    if ((params.AdvSolType != AdvectionSolverType::TG2) and
        (params.AdvSolType != AdvectionSolverType::CG2))
    {
        n.GetMesh()->ExchangeData(m_half_tag, NODE, 0);       
    }
    else
    {
        n.GetMesh()->ExchangeData(m_high_tag, NODE, 0);
    }
    BARRIER

    // make second step if TTG2, TTG3 or TTG4 is used

    if ((params.AdvSolType != AdvectionSolverType::TG2) and
        (params.AdvSolType != AdvectionSolverType::CG2))
    {
        // Reset RHS and SOL vector
        SOL.Clear();
        RHS.Clear();
        SOL.SetInterval(idmin, idmax);
        RHS.SetInterval(idmin, idmax);
        BARRIER
        AssembleDoubleStepRHS(StepNumber::second);
        BARRIER
        Sol.Solve(RHS, SOL);

        // finally update m_high value
        for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_high_tag) = SOL[nodeit->Integer(id)];
            }
        }
        BARRIER
        n.GetMesh()->ExchangeData(m_high_tag, NODE, 0);
    }

    // Reset SOL vector
    SOL.Clear();
    SOL.SetInterval(idmin, idmax);
    BARRIER
    
    // Calculate low order solution if FCT is used
    if (params.is_fct)
    {
        // make old MASS vector and tmp vectors
        Sparse::Vector MASS, tmp_vec_1, tmp_vec_2, tmp_vec_3;
        MASS.SetInterval(idmin, idmax);
        tmp_vec_1.SetInterval(idmin, idmax);
        tmp_vec_2.SetInterval(idmin, idmax);
        tmp_vec_3.SetInterval(idmin, idmax);

        // fill mass vector
        for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                MASS[nodeit->Integer(id)] = nodeit->Real(m_tag);
            }
        }
        BARRIER
        mat_mult_vec(LHS, MASS, tmp_vec_1);
        BARRIER
        mat_mult_vec(LHS_low, MASS, tmp_vec_2);
        BARRIER
        vec_mult_num(tmp_vec_1, (params.fct_cd - 1.0), idmin, idmax);
        BARRIER
        vec_mult_num(tmp_vec_2, (1.0 - params.fct_cd), idmin, idmax);
        BARRIER
        vec_plus_vec(tmp_vec_1, tmp_vec_2, tmp_vec_3, idmin, idmax);
        BARRIER
        vec_plus_vec(tmp_vec_3, RHS, RHS_low, idmin, idmax);
        BARRIER


        // solve low order system 
        Sol.SetMatrix(LHS_low);
        Sol.Solve(RHS_low, SOL);

        // update m_low value
        for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_low_tag) = SOL[nodeit->Integer(id)];
            }
        }
        BARRIER
        n.GetMesh()->ExchangeData(m_low_tag, NODE, 0);

        // clear RHS_low vector
        RHS_low.Clear();
        RHS_low.SetInterval(idmin, idmax); 
    }

    // clear RHS and SOL vector
    SOL.Clear();
    RHS.Clear();
    SOL.SetInterval(idmin, idmax);
    RHS.SetInterval(idmin, idmax);

    // apply FCT (if used) and update m finally
    if (params.is_fct)
    {
        FCT_Procedure(n, M_C_minus_M_L, m_tag, m_low_tag, m_high_tag, params.fct_cd);
    }
    else
    {
       for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_tag) = nodeit->Real(m_high_tag);
            }
        }
        n.GetMesh()->ExchangeData(m_tag, NODE, 0); 
    }
    BARRIER

    // clear all temporal data
    n.GetMesh()->DeleteTag(m_high_tag);
    if (params.is_fct)
    {
        n.GetMesh()->DeleteTag(m_low_tag);
    }
    if ((params.AdvSolType != AdvectionSolverType::TG2) and
        (params.AdvSolType != AdvectionSolverType::CG2))
    {
        n.GetMesh()->DeleteTag(m_half_tag);
    }

    BARRIER

    if (current_step_number == (total_step_num - 1))
    {
        INMOST::Tag m_diff_tag = n.GetMesh()->CreateTag("m_diff", DATA_REAL, NODE, NONE, 1);
        for( Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_diff_tag) = nodeit->Real(m_tag) - nodeit->Real(m_init_tag);
            }
        }
        n.GetMesh()->ExchangeData(m_diff_tag, NODE, 0);
    }
    BARRIER

    // log if verbose
    if (output_params.is_verbose)
    {
        LogStep();   
    }

    BARRIER
    // increment step num
    ++current_step_number;
}

void AdvectionSolver::LogInit()
{
    if (n.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "========= Advection Solver Info =========" << std::endl;
        std::string mass_test;
        if (params.AdvSolType == AdvectionSolverType::TG2)
        {
            std::cout << "Advection solver type: Taylor-Galerkin (2 order);" << std::endl;
            if (params.is_fct)
            {
                std::cout << "With Flux Correction applied (c_d = " << params.fct_cd << ");" << std::endl;
            }
            std::cout << "Numerical integration time step = " << params.time_step_sec/3600.0 << " h.;" << std::endl;
            std::cout << "Every " << output_frequency << " screenshot will be printed;" << std::endl;
        }
        else if (params.AdvSolType == AdvectionSolverType::CG2)
        {
            std::cout << "Advection solver type: characteristics Galerkin (2 order);" << std::endl;
            if (params.is_fct)
            {
                std::cout << "With Flux Correction applied (c_d = " << params.fct_cd << ");" << std::endl;
            }
            std::cout << "Numerical integration time step = " << params.time_step_sec/3600.0 << " h.;" << std::endl;
            std::cout << "Every " << output_frequency << " screenshot will be printed;" << std::endl;
        }
        else if (params.AdvSolType == AdvectionSolverType::TTG2)
        {
            std::cout << "Advection solver type: Two-step Taylor-Galerkin (2 order);" << std::endl;
            if (params.is_fct)
            {
                std::cout << "With Flux Correction applied (c_d = " << params.fct_cd << ");" << std::endl;
            }
            std::cout << "Numerical integration time step = " << params.time_step_sec/3600.0 << " h.;" << std::endl;
            std::cout << "Every " << output_frequency << " screenshot will be printed;" << std::endl;
        }
        else if (params.AdvSolType == AdvectionSolverType::TTG3)
        {
            std::cout << "Advection solver type: Two-step Taylor-Galerkin (3 order);" << std::endl;
            if (params.is_fct)
            {
                std::cout << "With Flux Correction applied (c_d = " << params.fct_cd << ");" << std::endl;
            }
            std::cout << "Numerical integration time step = " << params.time_step_sec/3600.0 << " h.;" << std::endl;
            std::cout << "Every " << output_frequency << " screenshot will be printed;" << std::endl;
        }
        else if (params.AdvSolType == AdvectionSolverType::TTG4)
        {
            std::cout << "Advection solver type: Two-step Taylor-Galerkin (4 order);" << std::endl;
            if (params.is_fct)
            {
                std::cout << "With Flux Correction applied (c_d = " << params.fct_cd << ");" << std::endl;
            }
            std::cout << "Numerical integration time step = " << params.time_step_sec/3600.0 << " h.;" << std::endl;
            std::cout << "Every " << output_frequency << " screenshot will be printed;" << std::endl;
        }
        else 
        {
            INMOST_ICE_ERR("unknown type of advection solver");
        }
        std::cout << "=========================================" << std::endl << std::endl;
    }
    BARRIER
}

void AdvectionSolver:: LogStep()
{
    if (current_step_number%output_frequency == 0)
    {
        double relative_mass_integral = IntegralMass(n, m_tag)/init_mass_integral;
        double minimal_mass = min_f(n, m_tag);
        double maximal_mass = max_f(n, m_tag);

        BARRIER

        if (n.GetMesh()->GetProcessorRank() == 0)
        {
            std::cout << "============ Step " << (current_step_number + 1) << " of " << params.total_step_number << "============" <<std::endl;
            std::cout << "Current mass integral/Initial mass integral = " << relative_mass_integral << ";" << std::endl; 
            mass_integral_vector.push_back(relative_mass_integral);
            minimal_mass_vector.push_back(minimal_mass);
            maximal_mass_vector.push_back(maximal_mass);
        }

        std::stringstream ss;
        ss << std::setfill('0') << std::setw(5) << current_step_number;
        n.PrintPVTU(output_params.output_dir + "test" + ss.str() + ".pvtu");
    }

    if (current_step_number == 0)
    {
        n.PrintPVTU(output_params.keymoments_dir + "start.pvtu");
        std::stringstream ss;
    }

    if (current_step_number == (total_step_num/2))
    {
        n.PrintPVTU(output_params.keymoments_dir + "half.pvtu");
    }

    if (current_step_number == (total_step_num - 1))
    {
        n.PrintPVTU(output_params.keymoments_dir + "end.pvtu");
    }
    
    if (current_step_number%output_frequency == 0)
    {
        if (n.GetMesh()->GetProcessorRank() == 0)
        {
            std::cout << "=========================================" << std::endl << std::endl;
        }
        BARRIER
    }
}

void PrintMassIntegralToFile(const std::string& mass_output_path,
                             const std::vector<double>& ms,
                             const std::vector<double>& minn,
                             const std::vector<double>& maxx)
{
    std::ofstream outp(mass_output_path);
    outp << "mass/init_mass" << std::endl;
    for (const auto& item : ms)
    {
        outp  << item << std::endl; 
    }
    outp << std::endl;
    outp << "min(mass)" << std::endl;
    for (const auto& item : minn)
    {
        outp  << item << std::endl; 
    }
    outp << std::endl;
    outp << "max(mass)" << std::endl;
    for (const auto& item : maxx)
    {
        outp  << item << std::endl; 
    }
    outp <<std::endl;
}

AdvectionSolver::~AdvectionSolver()
{
    total_time = Timer() - total_time;
    if (n.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "=============== Final Info ===============" <<std::endl;
        std::cout << "Total time = " << total_time << " seconds;" << std::endl;
        std::cout << "Average iteration time = " << total_time/current_step_number << " seconds;" << std::endl;
        PrintMassIntegralToFile(output_params.relative_mass_file,
                                mass_integral_vector,
                                minimal_mass_vector,
                                maximal_mass_vector);
    }

    if (output_params.print_errors)
    {      
        INMOST_ICE::PrintErrors(n, m_tag, m_init_tag);
    }

    if (n.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "=========================================" << std::endl << std::endl;
    }

}