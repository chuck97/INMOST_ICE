#include "Initialization.h"

using namespace INMOST;

std::vector<double> VelocityFieldUpdate(double t,
                                        double lon,
                                        double lat,
                                        double final_time,
                                        double velocity_scale_factor,
                                        VelocityFieldType velocity_type)
{
    if (velocity_type == VelocityFieldType::NonDiv1)
    {
        double u = velocity_scale_factor*
                   std::sin(lon/2.0)*std::sin(lon/2.0)* 
                   std::sin(lat*2.0)*
                   std::cos(M_PI*t/final_time);
        double v = (velocity_scale_factor/2.0)*
                   std::sin(lon)* 
                   std::cos(lat)*
                   std::cos(M_PI*t/final_time);
        return {u, v};
    }
    else if (velocity_type == VelocityFieldType::NonDiv2)
    {
        double u = velocity_scale_factor*
                   std::sin(lon)*std::sin(lon)* 
                   std::sin(lat*2.0)*
                   std::cos(M_PI*t/final_time);
        double v = velocity_scale_factor*
                   std::sin(2.0*lon)* 
                   std::cos(lat)*
                   std::cos(M_PI*t/final_time);
        return {u, v};
    }
    else if (velocity_type == VelocityFieldType::Div)
    {
        double u = -velocity_scale_factor*
                   std::sin(lon/2.0)*std::sin(lon/2.0)* 
                   std::sin(lat*2.0)*
                   std::cos(lat)*std::cos(lat)*
                   std::cos(M_PI*t/final_time);
        double v = (velocity_scale_factor/2.0)*
                   std::sin(lon)* 
                   std::cos(lat)*std::cos(lat)*std::cos(lat)*
                   std::cos(M_PI*t/final_time);
        return {u, v};
    }
    else if (velocity_type == VelocityFieldType::NonDiv3)
    {
        double lonn = lon - 2.0*M_PI*t/final_time;
        double u = velocity_scale_factor*
                   std::sin(lonn)*std::sin(lonn)* 
                   std::sin(lat*2.0)*
                   std::cos(M_PI*t/final_time) +
                   2*M_PI*EARTH_RADIUS*std::cos(lat)/(final_time*3600.0);
        double v = velocity_scale_factor*
                   std::sin(2.0*lonn)* 
                   std::cos(lat)*
                   std::cos(M_PI*t/final_time);
        return {u, v};
    }
    else if (velocity_type == VelocityFieldType::NonDiv4)
    {
        double u = 0.0;
        double v = 2*M_PI*EARTH_RADIUS*std::sin(lat/2.0)/(final_time*3600.0);
        return {u, v};
    }
    else
    {
        INMOST_ICE_ERR("unknown velocity case");
    }
}

// mass assignment function
double InitialMassAssignment(double lon,
                             double lat,
                             InitialMassType mass_type,
                             InitialMassParams mass_params)
{
    if (mass_type == InitialMassType::cosine_bells)
    {
        double r1 = std::acos(std::sin(mass_params.first_center_lat)*std::sin(lat) +
                         std::cos(mass_params.first_center_lat)*std::cos(lat)*
                         std::cos(lon - mass_params.first_center_lon));
        
        double r2 = std::acos(std::sin(mass_params.second_center_lat)*std::sin(lat) +
                         std::cos(mass_params.second_center_lat)*std::cos(lat)*
                         std::cos(lon - mass_params.second_center_lon));

        if (r1 < mass_params.radius)
        {
            double h1 = (mass_params.amplitude/2.0)*(1.0 + std::cos(M_PI*r1/mass_params.radius));
            return (mass_params.background + mass_params.scale_factor*h1);
        }
        else if (r2 < mass_params.radius)
        {
            double h2 = (mass_params.amplitude/2.0)*(1.0 + std::cos(M_PI*r2/mass_params.radius));
            return (mass_params.background + mass_params.scale_factor*h2);
        }
        else
        {
            return mass_params.background;
        }
    }
    else if (mass_type == InitialMassType::Gaussian_hills)
    {
        double x = std::cos(lat)*std::cos(lon);
        double y = std::cos(lat)*std::sin(lon);
        double z = std::sin(lat);

        double x1 = std::cos(mass_params.first_center_lat)*std::cos(mass_params.first_center_lon);
        double y1 = std::cos(mass_params.first_center_lat)*std::sin(mass_params.first_center_lon);
        double z1 = std::sin(mass_params.first_center_lat);

        double x2 = std::cos(mass_params.second_center_lat)*std::cos(mass_params.second_center_lon);
        double y2 = std::cos(mass_params.second_center_lat)*std::sin(mass_params.second_center_lon);
        double z2 = std::sin(mass_params.second_center_lat);

        double h1 = mass_params.amplitude*std::exp( -mass_params.width*
                                                    ((x - x1)*(x - x1) +
                                                     (y - y1)*(y - y1) +
                                                     (z - z1)*(z - z1)));

        double h2 = mass_params.amplitude*std::exp(-mass_params.width*
                                                    ((x - x2)*(x - x2) +
                                                     (y - y2)*(y - y2) +
                                                     (z - z2)*(z - z2)));
        return (h1 + h2);
    }
    else if (mass_type == InitialMassType::slotted_cylinders)
    {
        double r1 = std::acos(std::sin(mass_params.first_center_lat)*std::sin(lat) +
                         std::cos(mass_params.first_center_lat)*std::cos(lat)*
                         std::cos(lon - mass_params.first_center_lon));
        
        double r2 = std::acos(std::sin(mass_params.second_center_lat)*std::sin(lat) +
                         std::cos(mass_params.second_center_lat)*std::cos(lat)*
                         std::cos(lon - mass_params.second_center_lon));

        if (((r1 <= mass_params.radius) and
             (std::fabs(lon - mass_params.first_center_lon) >= mass_params.radius/6.0)) or
            ((r2 <= mass_params.radius) and
             (std::fabs(lon - mass_params.second_center_lon) >= mass_params.radius/6.0)))
        {
            return mass_params.scale_factor;
        }
        else if (((r1 <= mass_params.radius) and
                  ((std::fabs(lon - mass_params.first_center_lon) < mass_params.radius/6.0))
                  and ((lat - mass_params.first_center_lat) < -(5.0/12.0)*mass_params.radius)) or
                  ((r2 <= mass_params.radius) and
                  ((std::fabs(lon - mass_params.second_center_lon) < mass_params.radius/6.0))
                  and ((lat - mass_params.second_center_lat) > (5.0/12.0)*mass_params.radius)))
        {
            return mass_params.scale_factor;
        }
        else
        {
            return mass_params.background;
        }
    }
    else
    {
        INMOST_ICE_ERR("unknown initial mass case");
    }
}

AdvectionTest::AdvectionTest(INMOST_ICE_nodes& n_,
                             const InitialMassParams& m_params_,
                             const VelocityParams& u_params_,
                             INMOST::Tag m_tag_,
                             INMOST::Tag init_m_tag_,
                             INMOST::Tag u_tag_,
                             double total_time_hours_):
    n(n_),
    m_params(m_params_),
    u_params(u_params_),
    m_tag(m_tag_),
    init_m_tag(init_m_tag_),
    u_tag(u_tag_),
    total_time_hours(total_time_hours_)
{
    if (u_params.velocity_type == VelocityFieldType::NonDiv3)
    {
        double vel_sc = u_params.scale_factor + 2.0*M_PI*EARTH_RADIUS/(total_time_hours*3600.0); 
        time_step_seconds = u_params.Courant_number*n.GetResolution()/vel_sc;
    }
    else if (u_params.velocity_type == VelocityFieldType::NonDiv4)
    {
        double vel_sc = 2.0*M_PI*EARTH_RADIUS/(total_time_hours*3600.0); 
        time_step_seconds = u_params.Courant_number*n.GetResolution()/vel_sc;
    }
    else
    {
        time_step_seconds = u_params.Courant_number*n.GetResolution()/u_params.scale_factor;    
    }
    
    time_step_hours = time_step_seconds/3600.0;
    ntotsteps = (size_t)(total_time_hours/time_step_hours); 

    Log();
}

void AdvectionTest::Log()
{
    if (n.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "=========== Initial Mass Info ===========" << std::endl;
        std::string mass_test;
        if (m_params.mass_type == InitialMassType::cosine_bells)
        {
            std::cout << "Initial Mass: cosine bells;" << std::endl;
            std::cout << "Amplitude = " << m_params.amplitude << ";" << std::endl;
            std::cout << "Scale factor = " << m_params.scale_factor << ";" << std::endl;
            std::cout << "Radius = " << m_params.radius << ";" << std::endl;
            std::cout << "Background = " << m_params.background << ";" << std::endl;
            std::cout << "First center coords = (" << m_params.first_center_lon << ", "
                                                   << m_params.first_center_lat << ");" 
                                                   << std::endl;
            std::cout << "Second center coords = (" << m_params.second_center_lon << ", "
                                                   << m_params.second_center_lat << ");" 
                                                   << std::endl;
        }
        else if (m_params.mass_type == InitialMassType::Gaussian_hills)
        {
            std::cout << "Initial Mass: Gaussian hills" << std::endl;
            std::cout << "Amplitude = " << m_params.amplitude << ";" << std::endl;
            std::cout << "Width = " << m_params.width << ";" << std::endl;
            std::cout << "First center coords = (" << m_params.first_center_lon << ", "
                                                   << m_params.first_center_lat << ");" 
                                                   << std::endl;
            std::cout << "Second center coords = (" << m_params.second_center_lon << ", "
                                                   << m_params.second_center_lat << ");" 
                                                   << std::endl;
        }
        else if (m_params.mass_type == InitialMassType::slotted_cylinders)
        {
            std::cout << "Initial Mass: slotted cylinders" << std::endl;
            std::cout << "Scale factor = " << m_params.scale_factor << ";" << std::endl;
            std::cout << "Radius = " << m_params.radius << ";" << std::endl;
            std::cout << "Background = " << m_params.background << ";" << std::endl;
            std::cout << "First center coords = (" << m_params.first_center_lon << ", "
                                                   << m_params.first_center_lat << ");" 
                                                   << std::endl;
            std::cout << "Second center coords = (" << m_params.second_center_lon << ", "
                                                   << m_params.second_center_lat << ");" 
                                                   << std::endl;
        }
        else
        {
            INMOST_ICE_ERR("unknown type of initial mass distribution");
        }
        std::cout << "=========================================" << std::endl << std::endl;

        std::cout << "=========== Velocity Field Info =========" << std::endl;
        if (u_params.velocity_type == VelocityFieldType::NonDiv1)
        {
            std::cout << "Velocity Field: non-divirgent №1;" << std::endl;
            std::cout << "Scale factor = " << u_params.scale_factor << " m/s;" << std::endl;
            std::cout << "Courant number = " << u_params.Courant_number << ";" << std::endl;
        }
        else if (u_params.velocity_type == VelocityFieldType::NonDiv2)
        {
            std::cout << "Velocity Field: non-divirgent №2;" << std::endl;
            std::cout << "Scale factor = " << u_params.scale_factor << " m/s;" << std::endl;
            std::cout << "Courant number = " << u_params.Courant_number << ";" << std::endl;
        }
        else if (u_params.velocity_type == VelocityFieldType::NonDiv3)
        {
            std::cout << "Velocity Field: non-divirgent №3;" << std::endl;
            std::cout << "Scale factor = " << u_params.scale_factor << " m/s;" << std::endl;
            std::cout << "Courant number = " << u_params.Courant_number << ";" << std::endl;
        }
        else if (u_params.velocity_type == VelocityFieldType::Div)
        {
            std::cout << "Velocity Field: divirgent;" << std::endl;
            std::cout << "Scale factor = " << u_params.scale_factor << " m/s;" << std::endl;
            std::cout << "Courant number = " << u_params.Courant_number << ";" << std::endl;
        }
        else if (u_params.velocity_type == VelocityFieldType::NonDiv4)
        {
            std::cout << "Velocity Field: non-divergent №4;" << std::endl;
            std::cout << "Courant number = " << u_params.Courant_number << ";" << std::endl;
        }
        else
        {
            INMOST_ICE_ERR("unknown type of velocity field");
        }
        std::cout << "=========================================" << std::endl << std::endl;

        std::cout << "=========== General Info ================" << std::endl;
        std::cout << "Processors number = " << n.GetMesh()->GetProcessorsNumber() << ";" <<std::endl;
        std::cout << "Total time = " << total_time_hours << " hours;" << std::endl;
        std::cout << "Time step = " << time_step_hours << " hours;" << std::endl;
        std::cout << "Number of steps = " << ntotsteps << ";" << std::endl;
        std::cout << "Spacial resolution = " << n.GetResolution()/1000.0 << " km.;" <<std::endl;
        std::cout << "=========================================" << std::endl << std::endl;
    }
    BARRIER
}

double AdvectionTest::GetTimeStepHours() const
{
    return time_step_hours;
}

double AdvectionTest::GetTimeStepSeconds() const
{
    return time_step_seconds;
}

double AdvectionTest::GetTotalTimeHours() const
{
    return total_time_hours;
}

size_t AdvectionTest::GetTotalStepNumber() const
{
    return ntotsteps;
}

void AdvectionTest::AssignInitialMass()
{
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
            nodeit != n.GetMesh()->EndNode();
            ++nodeit)
	{
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double lon = nodeit->RealArray(n.GetCoords().geo_coords)[0];
            double lat = nodeit->RealArray(n.GetCoords().geo_coords)[1];
            double cur_m = InitialMassAssignment(lon, lat, m_params.mass_type, m_params);
            nodeit->Real(init_m_tag) = cur_m;
            nodeit->Real(m_tag) = cur_m;
        }
    }
    BARRIER
    // exchange data
    n.GetMesh()->ExchangeData(init_m_tag, NODE, 0);
    n.GetMesh()->ExchangeData(m_tag, NODE, 0);
    BARRIER
}

void AdvectionTest::UpdateVelocity(double current_time_hours)
{
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
        nodeit != n.GetMesh()->EndNode();
        ++nodeit) 
	{
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double lon = nodeit->RealArray(n.GetCoords().geo_coords)[0];
            double lat = nodeit->RealArray(n.GetCoords().geo_coords)[1];
            auto cur = VelocityFieldUpdate(current_time_hours,
                                           lon,
                                           lat,
                                           total_time_hours,
                                           u_params.scale_factor,
                                           u_params.velocity_type);

            nodeit->RealArray(u_tag)[0] = cur[0];
            nodeit->RealArray(u_tag)[1] = cur[1];
            nodeit->RealArray(u_tag)[2] = 0.0;
        }
    }
    BARRIER
    // exchange data
    n.GetMesh()->ExchangeData(u_tag, NODE, 0);
    BARRIER
}

void ParseJson(const std::string& json_path,
               std::string& mesh_path,
               double& total_time,
               InitialMassParams& init_mass_params,
               VelocityParams& velocity_params,
               AdvectionSolverParams& AdvSolParams,
               OutputParameters& OutpParams)
{
    std::string no_spaces_json = json_path;
    rtrim(no_spaces_json);
    if(no_spaces_json.substr(no_spaces_json.size()-5, 5) != ".json")
    {
        INMOST_ICE_ERR("input file shoud be ended by .json");
    }

    std::ifstream ifs(no_spaces_json);
    nlohmann::json j_input = nlohmann::json::parse(ifs);
    BARRIER

    // Parse mesh path
    if (j_input["Mesh"].empty())
    {
        INMOST_ICE_ERR("Mesh path should be given (\"Mesh\": \"path_to_mesh.pmf\")");
    }
    else
    {
        std::string m_path = j_input["Mesh"];
        rtrim(m_path);
        mesh_path = m_path;
    }
    BARRIER
    
    //Parse total time
    if (j_input["Total time (h)"].empty())
    {
        INMOST_ICE_ERR("Total time should be given (\"Total time (h)\": \"time_step\")");
    }
    else
    {
        total_time = j_input["Total time (h)"];
    }

    //Parse initial mass data
    if (j_input["Initial Mass"].empty())
    {
        INMOST_ICE_ERR("Initial mass params should be given in object \" Initial Mass \"")
    }
    else
    {
        // m type
        if(j_input["Initial Mass"]["type"].empty())
        {
            INMOST_ICE_ERR("Initial mass type should be given in \" Initial Mass::type \"")
        }
        else
        {
            std::string type = j_input["Initial Mass"]["type"];
            if (type == "cosine bells")
            {
                init_mass_params.mass_type = InitialMassType::cosine_bells;
            }
            else if (type == "Gaussian hills")
            {
                init_mass_params.mass_type = InitialMassType::Gaussian_hills;
            }
            else if (type == "slotted cylinders")
            {
                init_mass_params.mass_type = InitialMassType::slotted_cylinders;
            }
            else
            {
                INMOST_ICE_ERR("Initial mass type should be given (Possibilities: cosine bells, Gaussian hills, slotted cylinders)");
            }
        }

        // amplitude
        if(j_input["Initial Mass"]["amplitude"].empty())
        {
            INMOST_ICE_ERR("Amplitude should be given in \" Initial Mass::amplitude \"")
        }
        else
        {
            init_mass_params.amplitude = j_input["Initial Mass"]["amplitude"];
        }

        // scale factor
        if(j_input["Initial Mass"]["scale factor"].empty())
        {
            INMOST_ICE_ERR("Scale type should be given in \" Initial Mass::scale factor \"")
        }
        else
        {
            init_mass_params.scale_factor = j_input["Initial Mass"]["scale factor"];
        }

        // radius
        if(j_input["Initial Mass"]["radius"].empty())
        {
            INMOST_ICE_ERR("Radius should be given in \" Initial Mass::radius \"")
        }
        else
        {
            init_mass_params.radius = j_input["Initial Mass"]["radius"];
        }

        // background
        if(j_input["Initial Mass"]["background"].empty())
        {
            INMOST_ICE_ERR("Background should be given in \" Initial Mass::radius \"")
        }
        else
        {
            init_mass_params.background = j_input["Initial Mass"]["background"];
        }

        // width
        if(j_input["Initial Mass"]["width"].empty())
        {
            INMOST_ICE_ERR("Width should be given in \" Initial Mass::width \"")
        }
        else
        {
            init_mass_params.width = j_input["Initial Mass"]["width"];
        }

        // first center coords
        if(j_input["Initial Mass"]["first center coords"].empty())
        {
            INMOST_ICE_ERR("First center coords should be given in \" Initial Mass::first center coords \"")
        }
        else if (j_input["Initial Mass"]["first center coords"].size() != 2)
        {
            INMOST_ICE_ERR("First center coords should have size 2")
        }
        else
        {
            init_mass_params.first_center_lon = j_input["Initial Mass"]["first center coords"][0];
            init_mass_params.first_center_lat = j_input["Initial Mass"]["first center coords"][1];
        }

        // second center coords
        if(j_input["Initial Mass"]["second center coords"].empty())
        {
            INMOST_ICE_ERR("Second center coords should be given in \" Initial Mass::second center coords \"")
        }
        else if (j_input["Initial Mass"]["second center coords"].size() != 2)
        {
            INMOST_ICE_ERR("Second center coords should have size 2")
        }
        else
        {
            init_mass_params.second_center_lon = j_input["Initial Mass"]["second center coords"][0];
            init_mass_params.second_center_lat = j_input["Initial Mass"]["second center coords"][1];
        }    
    }

    //Parse velocity params
    if (j_input["Velocity Field"].empty())
    {
        INMOST_ICE_ERR("Velocity Fields params should be given in object \" Velocity Field \"")
    }
    else
    {
        // velocity type
        if(j_input["Velocity Field"]["type"].empty())
        {
            INMOST_ICE_ERR("Initial mass type should be given in \" Velocity Field::type \"")
        }
        else
        {
            std::string type = j_input["Velocity Field"]["type"];
            if (type == "non-div1")
            {
                velocity_params.velocity_type = VelocityFieldType::NonDiv1;
            }
            else if (type == "non-div2")
            {
                velocity_params.velocity_type = VelocityFieldType::NonDiv2;
            }
            else if (type == "non-div3")
            {
                velocity_params.velocity_type = VelocityFieldType::NonDiv3;
            }
            else if (type == "non-div4")
            {
                velocity_params.velocity_type = VelocityFieldType::NonDiv4;
            }
            else if (type == "div")
            {
                velocity_params.velocity_type = VelocityFieldType::Div;
            }
            else
            {
                INMOST_ICE_ERR("Velocity Field type should be given (Possibilities: non-div1, non-div2, non-div3, div)");
            }
        }

        // Courant number
        if(j_input["Velocity Field"]["Courant number"].empty())
        {
            INMOST_ICE_ERR("Courant number should be given in \" Velocity Field::Courant number \"")
        }
        else
        {
            velocity_params.Courant_number = j_input["Velocity Field"]["Courant number"];
        }

        // scale factor
        if(j_input["Velocity Field"]["scale factor (m/s)"].empty())
        {
            INMOST_ICE_ERR("Scale factor should be given in \" Velocity Field::scale factor (m/s) \"")
        }
        else
        {
            velocity_params.scale_factor = j_input["Velocity Field"]["scale factor (m/s)"];
        }
    }
    
    // Parse Advection Solver params
    if (j_input["Advection Solver"].empty())
    {
        INMOST_ICE_ERR("Advection Solver params should be given in object \" Advection Solver \"")
    }
    else
    {
        //type
        if(j_input["Advection Solver"]["type"].empty())
        {
            INMOST_ICE_ERR("Advection Solver type should be given in \" Advection Solver::type \"")
        }
        else
        {
            std::string type = j_input["Advection Solver"]["type"];
            if (type == "TG2")
            {
                AdvSolParams.AdvSolType = AdvectionSolverType::TG2;
            }
            else if (type == "CG2")
            {
                AdvSolParams.AdvSolType = AdvectionSolverType::CG2;
            }
            else if (type == "TTG2")
            {
                AdvSolParams.AdvSolType = AdvectionSolverType::TTG2;
            }
            else if (type == "TTG3")
            {
                AdvSolParams.AdvSolType = AdvectionSolverType::TTG3;
            }
            else if (type == "TTG4")
            {
                AdvSolParams.AdvSolType = AdvectionSolverType::TTG4;
            }
            else
            {
                INMOST_ICE_ERR("Advection solver type should be given (Possibilities: TG2, CG2, TTG2, TTG3, TTG4)");
            }
        }
        //FCT
        if(j_input["Advection Solver"]["FCT"].empty())
        {
            INMOST_ICE_ERR("FCT should be given in \" Advection Solver::FCT \"");
        }
        else
        {
            AdvSolParams.is_fct = j_input["Advection Solver"]["FCT"];
        }

        //cd
        if(j_input["Advection Solver"]["cd"].empty())
        {
            INMOST_ICE_ERR("FCT cd should be given in \" Advection Solver::cd \"");
        }
        else
        {
            AdvSolParams.fct_cd = j_input["Advection Solver"]["cd"];
        }
    }

    // Parse output parameters
    if (j_input["Output Parameters"].empty())
    {
        INMOST_ICE_ERR("Output Parameters should be given in object \" Output Parameters \"")
    }
    else
    {
        // verbose
        if(j_input["Output Parameters"]["verbose"].empty())
        {
            INMOST_ICE_ERR("verbose should be given in \" Output Parameters::verbose \"");
        }
        else
        {
            OutpParams.is_verbose = j_input["Output Parameters"]["verbose"];
        }

        // print errors
        if(j_input["Output Parameters"]["print errors"].empty())
        {
            INMOST_ICE_ERR("print errors should be given in \" Output Parameters::print errors \"");
        }
        else
        {
            OutpParams.print_errors = j_input["Output Parameters"]["print errors"];
        }

        // nouts
        if(j_input["Output Parameters"]["number of screenshots"].empty())
        {
            INMOST_ICE_ERR("number of screenshots should be given in \" Output Parameters::number of screenshots \"");
        }
        else
        {
            OutpParams.number_of_screenshots = j_input["Output Parameters"]["number of screenshots"];
        }

        // output directory
        if(j_input["Output Parameters"]["output directory"].empty())
        {
            INMOST_ICE_ERR("output directory should be given in \" Output Parameters::output directory \"");
        }
        else
        {
            std::string dir = j_input["Output Parameters"]["output directory"];
            rtrim(dir);
            if(dir.substr(dir.size()-1, 1) != "/")
            {
                INMOST_ICE_ERR("dir shoud be ended by /");
            }
            else
            {
                OutpParams.output_dir = j_input["Output Parameters"]["output directory"];
            }
        }

        // keymoments directory
        if(j_input["Output Parameters"]["keymoments directory"].empty())
        {
            INMOST_ICE_ERR("keymoments directory should be given in \" Output Parameters::keymoments directory \"");
        }
        else
        {
            std::string dir = j_input["Output Parameters"]["keymoments directory"];
            rtrim(dir);
            if(dir.substr(dir.size()-1, 1) != "/")
            {
                INMOST_ICE_ERR("dir shoud be ended by /");
            }
            else
            {
                OutpParams.keymoments_dir = j_input["Output Parameters"]["keymoments directory"];
            }
        }

        // mass integral file
        if(j_input["Output Parameters"]["mass file"].empty())
        {
            INMOST_ICE_ERR("mass file should be given in \" Output Parameters::mass file \"");
        }
        else
        {
            OutpParams.relative_mass_file = j_input["Output Parameters"]["mass file"];
        }
    }
}