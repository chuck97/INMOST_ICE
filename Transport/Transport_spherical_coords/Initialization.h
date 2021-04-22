#pragma once
#include "config.h"
#include "INMOST_ICE_mesh.h"
#include "Advection.h"

// velocity update function
enum class InitialMassType
{
    cosine_bells,
    Gaussian_hills,
    slotted_cylinders
};

enum class VelocityFieldType
{
    NonDiv1,
    NonDiv2,
    NonDiv3,
    NonDiv4,
    Div
};

struct InitialMassParams
{
    InitialMassType mass_type;
    double amplitude;
    double scale_factor;
    double radius;
    double background;
    double width;
    double first_center_lon;
    double first_center_lat;
    double second_center_lon;
    double second_center_lat;
};

struct VelocityParams
{
    VelocityFieldType velocity_type;
    double Courant_number;
    double scale_factor;
};

std::vector<double> VelocityFieldUpdate(double t,
                                        double lon,
                                        double lat,
                                        double final_time,
                                        double velocity_scale_factor,
                                        VelocityFieldType velocity_type);

// mass assignment function
double InitialMassAssignment(double lon,
                             double lat,
                             InitialMassType mass_type,
                             InitialMassParams mass_params);

class AdvectionTest
{
public:
    AdvectionTest(INMOST_ICE_nodes& n_,
                  const InitialMassParams& m_params_,
                  const VelocityParams& u_params_,
                  INMOST::Tag m_tag_,
                  INMOST::Tag init_m_tag_,
                  INMOST::Tag u_tag_,
                  double total_time_hours_);

    void Log();

    double GetTimeStepHours() const;
    double GetTimeStepSeconds() const;
    double GetTotalTimeHours() const;
    size_t GetTotalStepNumber() const;

    void AssignInitialMass();

    void UpdateVelocity(double current_time_hours);

private:
    INMOST_ICE_nodes& n; 
    InitialMassParams m_params;
    VelocityParams u_params;
    INMOST::Tag m_tag;
    INMOST::Tag init_m_tag;
    INMOST::Tag u_tag;
    double total_time_hours;
    double time_step_seconds;
    double time_step_hours;
    size_t ntotsteps;
};

void ParseJson(const std::string& json_path,
               std::string& mesh_path,
               double& total_time,
               InitialMassParams& init_mass_params,
               VelocityParams& velocity_params,
               AdvectionSolverParams& AdvSolParams,
               OutputParameters& OutpParams);