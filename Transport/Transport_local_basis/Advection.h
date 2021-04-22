#pragma once
#include "config.h"
#include "INMOST_ICE_mesh.h"
#include "Assembling.h"
#include "Errors.h"
#include "MatVec.h"
#include "FCT.h"

enum class AdvectionSolverType
{
    TG2,
    CG2,
    TTG2,
    TTG3,
    TTG4
};

enum class StepNumber
{
    first,
    second
};

struct AdvectionSolverParams
{
    AdvectionSolverType AdvSolType;
    double time_step_sec;
    bool is_fct;
    double fct_cd;
    size_t total_step_number;
};


struct OutputParameters
{   
    bool is_verbose;
    bool print_errors;
    size_t number_of_screenshots;
    std::string output_dir;
    std::string keymoments_dir;
    std::string relative_mass_file;
};

class AdvectionSolver
{

public:

    AdvectionSolver(INMOST_ICE_nodes& n_,
                    const AdvectionSolverParams& params_,
                    INMOST::Tag m_tag_,
                    INMOST::Tag m_init_tag_,
                    INMOST::Tag u_tag_,
                    INMOST::Solver& Sol_,
                    size_t total_step_num_,
                    OutputParameters& output_params_);
    void AssembleLHS();
    void AssembleRHS();
    void Evaluate();
    ~AdvectionSolver();
private:

    void AssembleSingleStepRHS();
    void AssembleDoubleStepRHS(StepNumber step_num);

    void LogInit();
    void LogStep();

    void PrintErrors();


    INMOST_ICE_nodes& n;
    AdvectionSolverParams params;
    INMOST::Tag m_tag;
    INMOST::Tag m_init_tag;
    INMOST::Tag u_tag;
    INMOST::Tag m_high_tag;
    INMOST::Tag m_low_tag;
    INMOST::Tag m_half_tag;
    INMOST::Sparse::Matrix LHS;
    INMOST::Sparse::Matrix LHS_low;
    INMOST::Sparse::Vector RHS;
    INMOST::Sparse::Vector RHS_low;
    INMOST::Solver& Sol;
    size_t total_step_num;
    std::vector<std::vector<std::vector<double>>> M_C_minus_M_L;
    size_t current_step_number = 0;
    double init_mass_integral;
    double total_time;
    std::vector<double> mass_integral_vector;
    std::vector<double> minimal_mass_vector;
    std::vector<double> maximal_mass_vector;
    OutputParameters output_params;
    size_t output_frequency;
};