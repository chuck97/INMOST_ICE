#pragma once
#include "config.h"
#include "INMOST_ICE_mesh.h"
#include "Assembling.h"

void FCT_Procedure(INMOST_ICE_nodes& n,
                   const std::vector<std::vector<std::vector<double>>>& M_C_minus_M_L,
                   INMOST::Tag m_tag,
                   INMOST::Tag m_low_tag,
                   INMOST::Tag m_high_tag,
                   double fct_cd);

void CalculateAlpha(INMOST_ICE_nodes& n,
                    INMOST::Tag& Fe_trian_tag,
                    INMOST::Tag& alpha_trian_tag,
                    INMOST::Tag m_tag,
                    INMOST::Tag m_low_tag);