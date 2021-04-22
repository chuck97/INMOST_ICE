#pragma once
#include "config.h"
#include "INMOST_ICE_mesh.h"
#include "Numerical.h"
#include "Function.h"

std::vector<std::vector<double>> GeoTrCorr(const std::vector<std::vector<double>>& tr);

double integral_over_subdomain_f(INMOST_ICE_nodes& n, INMOST::Tag& m);
double integral_over_subdomain_abs_f1_minus_f2(INMOST_ICE_nodes& n, INMOST::Tag& m1, INMOST::Tag& m2);
double integral_over_subdomain_squared_f1_minus_f2(INMOST_ICE_nodes& n, INMOST::Tag& m1, INMOST::Tag& m2);

double subdomain_max_f(INMOST_ICE_nodes& n, INMOST::Tag& m);
double subdomain_min_f(INMOST_ICE_nodes& n, INMOST::Tag& m);
double subdomain_max_abs_f(INMOST_ICE_nodes& n, INMOST::Tag& m);
double subdomain_max_abs_f1_minus_f2(INMOST_ICE_nodes& n, INMOST::Tag& m1, INMOST::Tag& m2);
double max_f(INMOST_ICE_nodes& n, INMOST::Tag& m);
double min_f(INMOST_ICE_nodes& n, INMOST::Tag& m);

double l1_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double l2_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double linf_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double phi_max_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double phi_min_error(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);

double IntegralMass(INMOST_ICE_nodes& n, INMOST::Tag& m_tag);
namespace INMOST_ICE
{
    void PrintErrors(INMOST_ICE_nodes& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
}