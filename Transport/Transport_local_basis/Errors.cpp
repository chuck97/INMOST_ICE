#include "Errors.h"

using namespace INMOST;

std::vector<std::vector<double>> GeoTrCorr(const std::vector<std::vector<double>>& tr)
{
    double lon0 = tr[0][0];
    double lat0 = tr[0][1];
    double lon1 = tr[1][0];
    double lat1 = tr[1][1];
    double lon2 = tr[2][0];
    double lat2 = tr[2][1];

    // correct triangle if needed
    std::vector<std::vector<double>> new_trnodes;

    if ((fabs(lon0 - lon1) > 2.0*M_PI/3.0) and (fabs(lon0 - lon2) > 2.0*M_PI/3.0))
    {
        if (lon0 < lon1)
        {
            lon0 += 2.0*M_PI;
        }
        else
        {
            lon0 -= 2.0*M_PI;
        }
        new_trnodes = {{lon0,lat0}, {lon1, lat1}, {lon2, lat2}};
    }
    else if ((fabs(lon1 - lon0) > 2.0*M_PI/3.0) and (fabs(lon1 - lon2) > 2.0*M_PI/3.0))
    {
        if (lon1 < lon2)
        {
            lon1 += 2.0*M_PI;
        }
        else
        {
            lon1 -= 2.0*M_PI;
        }
        new_trnodes = {{lon0,lat0}, {lon1, lat1}, {lon2, lat2}};
    }
    else if ((fabs(lon2 - lon0) > 2.0*M_PI/3.0) and (fabs(lon2 - lon1) > 2.0*M_PI/3.0))
    {
        if (lon2 < lon0)
        {
            lon2 += 2.0*M_PI;
        }
        else
        {
            lon2 -= 2.0*M_PI;
        }
        new_trnodes = {{lon0,lat0}, {lon1, lat1}, {lon2, lat2}};
    }
    else
    {
        new_trnodes = {{lon0,lat0}, {lon1, lat1}, {lon2, lat2}};
    }

    return new_trnodes;
}

double integral_over_subdomain_f(INMOST_ICE_nodes& n, INMOST::Tag& m)
{
    double res = 0.0;
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
        trianit != n.GetMesh()->EndCell();
        ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
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

            std::vector<double> fvalues = {local_nodes[0]->Real(m),
                                           local_nodes[1]->Real(m),
                                           local_nodes[2]->Real(m)};

            std::vector<double> local_params = Solve3x3({
                                                            {local_node_coords[0][0], local_node_coords[1][0], local_node_coords[2][0]},
                                                            {local_node_coords[0][1], local_node_coords[1][1], local_node_coords[2][1]},
                                                            {1.0                    , 1.0                    , 1.0                    }
                                                        }, fvalues);

            ScalarFunction f(FuncType::linear, local_params, local_node_coords);

            // calculate integral over triangle and add to res
            res += integral_over_triangle(local_node_coords, f);
        }
    }
    return res;
}

double integral_over_subdomain_abs_f1_minus_f2(INMOST_ICE_nodes& n, INMOST::Tag& m1, INMOST::Tag& m2)
{
    double res = 0.0;
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
        trianit != n.GetMesh()->EndCell();
        ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
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

            std::vector<double> f1values = {local_nodes[0]->Real(m1),
                                            local_nodes[1]->Real(m1),
                                            local_nodes[2]->Real(m1)};

            std::vector<double> f2values = {local_nodes[0]->Real(m2),
                                            local_nodes[1]->Real(m2),
                                            local_nodes[2]->Real(m2)};

            std::vector<double> local_params = Solve3x3({
                                                            {local_node_coords[0][0], local_node_coords[1][0], local_node_coords[2][0]},
                                                            {local_node_coords[0][1], local_node_coords[1][1], local_node_coords[2][1]},
                                                            {1.0                    , 1.0                    , 1.0                    }
                                                         }, {fabs(f1values[0] - f2values[0]),
                                                             fabs(f1values[1] - f2values[1]),
                                                             fabs(f1values[2] - f2values[2])});

            ScalarFunction f(FuncType::linear, local_params, local_node_coords);

            // calculate integral over triangle and add to res
            res += integral_over_triangle(local_node_coords, f);
        }
    }
    return res;
}

double integral_over_subdomain_squared_f1_minus_f2(INMOST_ICE_nodes& n, INMOST::Tag& m1, INMOST::Tag& m2)
{
    double res = 0.0;
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
        trianit != n.GetMesh()->EndCell();
        ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
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

            std::vector<double> f1values = {local_nodes[0]->Real(m1),
                                            local_nodes[1]->Real(m1),
                                            local_nodes[2]->Real(m1)};

            std::vector<double> f2values = {local_nodes[0]->Real(m2),
                                            local_nodes[1]->Real(m2),
                                            local_nodes[2]->Real(m2)};

            std::vector<double> local_params = Solve3x3({
                                                            {local_node_coords[0][0], local_node_coords[1][0], local_node_coords[2][0]},
                                                            {local_node_coords[0][1], local_node_coords[1][1], local_node_coords[2][1]},
                                                            {1.0                    , 1.0                    , 1.0                    }
                                                         }, {fabs(f1values[0] - f2values[0]),
                                                             fabs(f1values[1] - f2values[1]),
                                                             fabs(f1values[2] - f2values[2])});

            ScalarFunction f(FuncType::linear, local_params, local_node_coords);

            // calculate integral over triangle and add to res
            res += integral_over_triangle(local_node_coords, f*f);
        }
    }
    return res;
}

double subdomain_max_f(INMOST_ICE_nodes& n, INMOST::Tag& m)
{
    double res = std::numeric_limits<double>::min();
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
        nodeit != n.GetMesh()->EndNode();
        ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double val = nodeit->Real(m);
            if (val > res)
            {
                res = val;
            }
        }
    }
    return res;
}

double max_f(INMOST_ICE_nodes& n, INMOST::Tag& m)
{
    // init_m max
    double init_m_max = 0.0;
    double init_m_max_local = subdomain_max_f(n, m);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&init_m_max_local, &init_m_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    return init_m_max;
}

double subdomain_min_f(INMOST_ICE_nodes& n, INMOST::Tag& m)
{
    double res = std::numeric_limits<double>::max();
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
        nodeit != n.GetMesh()->EndNode();
        ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double val = nodeit->Real(m);
            if (val < res)
            {
                res = val;
            }
        }
    }
    return res;
}

double min_f(INMOST_ICE_nodes& n, INMOST::Tag& m)
{
    // init_m min
    double init_m_min = 0.0;
    double init_m_min_local = subdomain_min_f(n, m);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&init_m_min_local, &init_m_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
    return init_m_min;
}

double subdomain_max_abs_f(INMOST_ICE_nodes& n, INMOST::Tag& m)
{
    double res = std::numeric_limits<double>::min();
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
        nodeit != n.GetMesh()->EndNode();
        ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double val = fabs(nodeit->Real(m));
            if (val > res)
            {
                res = val;
            }
        }
    }
    return res;
}

double subdomain_max_abs_f1_minus_f2(INMOST_ICE_nodes& n, INMOST::Tag& m1, INMOST::Tag& m2)
{
    double res = std::numeric_limits<double>::min();
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
        nodeit != n.GetMesh()->EndNode();
        ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double val = fabs(nodeit->Real(m1) - nodeit->Real(m2));
            if (val > res)
            {
                res = val;
            }
        }
    }
    return res;
}

double l1_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag)
{
    // numerator
    double numerator = 0.0;
    double numerator_local = integral_over_subdomain_abs_f1_minus_f2(n,
                                                                     m_tag,
                                                                     m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&numerator_local, &numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif               
    
    // denominator
    INMOST::Tag zero;
    zero = n.GetMesh()->CreateTag("zero", DATA_REAL, NODE, NONE, 1);

    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
        nodeit != n.GetMesh()->EndNode();
        ++nodeit)
    {
        nodeit->Real(zero) = 0.0;
    }


    double denominator = 0.0;
    double denominator_local = integral_over_subdomain_abs_f1_minus_f2(n,
                                                                       m_init_tag,
                                                                       zero);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&denominator_local, &denominator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    n.GetMesh()->DeleteTag(zero);
    return (numerator/denominator);
}

double l2_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag)
{
    // numerator
    double numerator = 0.0;
    double numerator_local = integral_over_subdomain_squared_f1_minus_f2(n,
                                                                         m_tag,
                                                                         m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&numerator_local, &numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif               
    
    // denominator
    INMOST::Tag zero;
    zero = n.GetMesh()->CreateTag("zero", DATA_REAL, NODE, NONE, 1);
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
        nodeit != n.GetMesh()->EndNode();
        ++nodeit)
    {
        nodeit->Real(zero) = 0.0;
    }
    double denominator = 0.0;
    double denominator_local = integral_over_subdomain_squared_f1_minus_f2(n,
                                                                           m_init_tag,
                                                                           zero);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&denominator_local, &denominator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    n.GetMesh()->DeleteTag(zero);
    return (numerator/denominator);
}

double linf_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag)
{
    // numerator
    double numerator = 0.0;
    double numerator_local = subdomain_max_abs_f1_minus_f2(n, m_tag, m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&numerator_local, &numerator, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    // denominator
    double denominator = 0.0;
    double denominator_local = subdomain_max_abs_f(n, m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&denominator_local, &denominator, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    return (numerator/denominator);
}

double phi_max_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag)
{
    // init_m max
    double init_m_max = 0.0;
    double init_m_max_local = subdomain_max_f(n, m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&init_m_max_local, &init_m_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    // init_m min
    double init_m_min = 0.0;
    double init_m_min_local = subdomain_min_f(n, m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&init_m_min_local, &init_m_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

    // m max
    double m_max = 0.0;
    double m_max_local = subdomain_max_f(n, m_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&m_max_local, &m_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    return ((m_max - init_m_max)/(init_m_max - init_m_min));
}

double phi_min_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag)
{
    // init_m max
    double init_m_max = 0.0;
    double init_m_max_local = subdomain_max_f(n, m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&init_m_max_local, &init_m_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    // init_m min
    double init_m_min = 0.0;
    double init_m_min_local = subdomain_min_f(n, m_init_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&init_m_min_local, &init_m_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

    // m min
    double m_min = 0.0;
    double m_min_local = subdomain_min_f(n, m_tag);
    BARRIER

#if defined(USE_MPI)
   MPI_Allreduce(&m_min_local, &m_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

    return ((m_min - init_m_min)/(init_m_max - init_m_min));
}

double IntegralMass(INMOST_ICE_nodes& n, INMOST::Tag& m_tag)
{
    double local_mass = integral_over_subdomain_f(n, m_tag);
    double global_mass = 0.0;
    BARRIER 
#if defined(USE_MPI)
   MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    return global_mass;
}

void INMOST_ICE::PrintErrors(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag)
{
    BARRIER
    double l1_err = l1_error(n, m_tag, m_init_tag);
    BARRIER
    double l2_err = l2_error(n, m_tag, m_init_tag);
    BARRIER
    double linf_err = linf_error(n, m_tag, m_init_tag);
    BARRIER
    double phi_max_err = phi_max_error(n, m_tag, m_init_tag);
    BARRIER
    double phi_min_err = phi_min_error(n, m_tag, m_init_tag);
    BARRIER

    if (n.GetMesh()->GetProcessorRank() == 0)
    {
       std::cout << "l1 error is: " << l1_err << ";" << std::endl; 
       std::cout << "l2 error is: " << l2_err << ";" << std::endl;
       std::cout << "linf error is: " << linf_err << ";" << std::endl;
       std::cout << "phi max error is: " << phi_max_err << ";" << std::endl;
       std::cout << "phi min error is: " << phi_min_err << ";" << std::endl;
    }
}