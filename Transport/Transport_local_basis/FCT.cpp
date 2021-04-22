#include "FCT.h"

using namespace INMOST;

void FCT_Procedure(INMOST_ICE_nodes& n,
                   const std::vector<std::vector<std::vector<double>>>& M_C_minus_M_L,
                   INMOST::Tag m_tag,
                   INMOST::Tag m_low_tag,
                   INMOST::Tag m_high_tag,
                   double fct_cd)
{
    // calculate fluxes
    std::vector<std::vector<double>> Fe_vector(M_C_minus_M_L.size());
    int tr_num = 0;

    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
             trianit != n.GetMesh()->EndCell();
             ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
            std::vector<double> m_n_vec = {local_nodes[0]->Real(m_tag),
                                           local_nodes[1]->Real(m_tag),
                                           local_nodes[2]->Real(m_tag)};
            std::vector<double> m_H_vec = {local_nodes[0]->Real(m_high_tag),
                                           local_nodes[1]->Real(m_high_tag),
                                           local_nodes[2]->Real(m_high_tag)};
            std::vector<double> vec = (fct_cd - 1.0)*m_n_vec + m_H_vec;
            Fe_vector[tr_num] = M_C_minus_M_L[tr_num]*vec;
            ++tr_num;
        }
    }
    BARRIER

    // make fluxes tags on triangles
    INMOST::Tag Fe_trian_tag;
    tr_num = 0;
    Fe_trian_tag = n.GetMesh()->CreateTag("Fe_trian", DATA_REAL, CELL, NONE, 3);
    
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
             trianit != n.GetMesh()->EndCell();
             ++trianit)
    {
        if (trianit->GetStatus() != Element::Ghost)
        {
            trianit->RealArray(Fe_trian_tag)[0] = Fe_vector[tr_num][0];
            trianit->RealArray(Fe_trian_tag)[1] = Fe_vector[tr_num][1];
            trianit->RealArray(Fe_trian_tag)[2] = Fe_vector[tr_num][2];
            ++tr_num;
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(Fe_trian_tag, CELL, 0);
    BARRIER

    // make alpha tags on triangles
    INMOST::Tag alpha_trian_tag;
    tr_num = 0;
    alpha_trian_tag = n.GetMesh()->CreateTag("alpha_trian", DATA_REAL, CELL, NONE, 1);
    BARRIER

    // calculate coefficients
    CalculateAlpha(n, Fe_trian_tag, alpha_trian_tag, m_tag, m_low_tag);
    BARRIER

    // make m to be m_low
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            nodeit->Real(m_tag) = nodeit->Real(m_low_tag);
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(m_tag, NODE, 0);
    BARRIER

    // update solution
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
             trianit != n.GetMesh()->EndCell();
             ++trianit)
    {
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        if (local_nodes[0]->GetStatus() != Element::Ghost)
        {
            local_nodes[0]->Real(m_tag) += trianit->RealArray(Fe_trian_tag)[0]*
                                                trianit->Real(alpha_trian_tag);
        }

        if (local_nodes[1]->GetStatus() != Element::Ghost)
        {
            local_nodes[1]->Real(m_tag) += trianit->RealArray(Fe_trian_tag)[1]*
                                           trianit->Real(alpha_trian_tag);
        }

        if (local_nodes[2]->GetStatus() != Element::Ghost)
        {
            local_nodes[2]->Real(m_tag) += trianit->RealArray(Fe_trian_tag)[2]*
                                           trianit->Real(alpha_trian_tag);
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(m_tag, NODE, 0);
    BARRIER

    // Delete temporal Tags
    n.GetMesh()->DeleteTag(Fe_trian_tag, CELL);
    n.GetMesh()->DeleteTag(alpha_trian_tag, CELL);
}

void CalculateAlpha(INMOST_ICE_nodes& n, INMOST::Tag& Fe_trian_tag, INMOST::Tag& alpha_trian_tag, INMOST::Tag m_tag, INMOST::Tag m_low_tag)
{
    // calculate m_node_str: 0 - min, 1 - max
    INMOST::Tag m_node_str_tag;
    m_node_str_tag = n.GetMesh()->CreateTag("m_node_str", DATA_REAL, NODE, NONE, 2);
    
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double m_low = nodeit->Real(m_low_tag);
            double m_n = nodeit->Real(m_tag);
            nodeit->RealArray(m_node_str_tag)[0] = std::min(m_low, m_n);
            nodeit->RealArray(m_node_str_tag)[1] = std::max(m_low, m_n);
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(m_node_str_tag, NODE, 0);
    BARRIER

    // calculate m_trian_str_str: 0 - min, 1 - max
    INMOST::Tag m_trian_str_str_tag;
    m_trian_str_str_tag = n.GetMesh()->CreateTag("m_trian_str_str", DATA_REAL, CELL, NONE, 2);
    
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
             trianit != n.GetMesh()->EndCell();
             ++trianit)
    {
        if(trianit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
            double min1 = local_nodes[0]->RealArray(m_node_str_tag)[0];
            double min2 = local_nodes[1]->RealArray(m_node_str_tag)[0];
            double min3 = local_nodes[2]->RealArray(m_node_str_tag)[0];

            double max1 = local_nodes[0]->RealArray(m_node_str_tag)[1];
            double max2 = local_nodes[1]->RealArray(m_node_str_tag)[1];
            double max3 = local_nodes[2]->RealArray(m_node_str_tag)[1];

            double min_trian = std::min(std::min(min1, min2), min3);
            double max_trian = std::max(std::max(max1, max2), max3);

            trianit->RealArray(m_trian_str_str_tag)[0] = min_trian;
            trianit->RealArray(m_trian_str_str_tag)[1] = max_trian;
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(m_trian_str_str_tag, CELL, 0);
    BARRIER

    // calculate m_node_min_max: 0 - min, 1 - max
    INMOST::Tag m_node_min_max_tag;
    m_node_min_max_tag = n.GetMesh()->CreateTag("m_node_min_max", DATA_REAL, NODE, NONE, 2);
    
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Cell> local_trians = nodeit->getCells();
            double node_min = std::numeric_limits<double>::max();;
            double node_max = std::numeric_limits<double>::min();;
            for (size_t i = 0; i < local_trians.size(); ++i)
            {
                if (local_trians[i]->RealArray(m_trian_str_str_tag)[0] < node_min)
                {
                    node_min = local_trians[i]->RealArray(m_trian_str_str_tag)[0];
                }

                if (local_trians[i]->RealArray(m_trian_str_str_tag)[1] > node_max)
                {
                    node_max = local_trians[i]->RealArray(m_trian_str_str_tag)[1];
                }
            }
            nodeit->RealArray(m_node_min_max_tag)[0] = node_min;
            nodeit->RealArray(m_node_min_max_tag)[1] = node_max;
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(m_node_min_max_tag, NODE, 0);
    BARRIER

    // calculate q_node_min_max: 0 - min, 1 - max
    INMOST::Tag q_node_min_max_tag;
    q_node_min_max_tag = n.GetMesh()->CreateTag("q_node_min_max", DATA_REAL, NODE, NONE, 2);
    
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double node_m_min = nodeit->RealArray(m_node_min_max_tag)[0];
            double node_m_max = nodeit->RealArray(m_node_min_max_tag)[1];

            double m_low = nodeit->Real(m_low_tag);

            nodeit->RealArray(q_node_min_max_tag)[0] = node_m_min - m_low;
            nodeit->RealArray(q_node_min_max_tag)[1] = node_m_max - m_low;
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(q_node_min_max_tag, NODE, 0);
    BARRIER

    // calculate p_node_neg_pos: 0 - negative, 1 - positive
    INMOST::Tag p_node_neg_pos_tag;
    p_node_neg_pos_tag = n.GetMesh()->CreateTag("p_node_neg_pos", DATA_REAL, NODE, NONE, 2);
    
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double sum_neg = 0.0;
            double sum_pos = 0.0;
            INMOST::ElementArray<Cell> local_trians = nodeit->getCells();
            for (size_t i = 0; i < local_trians.size(); ++i)
            {
                INMOST::ElementArray<Node> local_nodes = local_trians[i]->getNodes();
                if ((fabs(local_nodes[0].Coords()[0] - nodeit->Coords()[0]) < VAREPS) and
                    (fabs(local_nodes[0].Coords()[1] - nodeit->Coords()[1]) < VAREPS))
                {
                    sum_neg += std::min(local_trians[i]->RealArray(Fe_trian_tag)[0], 0.0);
                    sum_pos += std::max(local_trians[i]->RealArray(Fe_trian_tag)[0], 0.0);
                }
                else if ((fabs(local_nodes[1].Coords()[0] - nodeit->Coords()[0]) < VAREPS) and
                         (fabs(local_nodes[1].Coords()[1] - nodeit->Coords()[1]) < VAREPS))
                {
                    sum_neg += std::min(local_trians[i]->RealArray(Fe_trian_tag)[1], 0.0);
                    sum_pos += std::max(local_trians[i]->RealArray(Fe_trian_tag)[1], 0.0);
                }
                else
                {
                    sum_neg += std::min(local_trians[i]->RealArray(Fe_trian_tag)[2], 0.0);
                    sum_pos += std::max(local_trians[i]->RealArray(Fe_trian_tag)[2], 0.0);
                }
            }

            nodeit->RealArray(p_node_neg_pos_tag)[0] = sum_neg;
            nodeit->RealArray(p_node_neg_pos_tag)[1] = sum_pos;
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(p_node_neg_pos_tag, NODE, 0);
    BARRIER

    // calculate R_node_neg_pos: 0 - negative, 1 - positive
    INMOST::Tag R_node_neg_pos_tag;
    R_node_neg_pos_tag = n.GetMesh()->CreateTag("R_node_neg_pos", DATA_REAL, NODE, NONE, 2);
    
    for(Mesh::iteratorNode nodeit = n.GetMesh()->BeginNode();
             nodeit != n.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            double P_neg = nodeit->RealArray(p_node_neg_pos_tag)[0];
            double P_pos = nodeit->RealArray(p_node_neg_pos_tag)[1];

            double Q_neg = nodeit->RealArray(q_node_min_max_tag)[0];
            double Q_pos = nodeit->RealArray(q_node_min_max_tag)[1];

            if (fabs(P_neg) < VAREPS)
            {
                nodeit->RealArray(R_node_neg_pos_tag)[0] = 0.0;
            }
            else
            {
                nodeit->RealArray(R_node_neg_pos_tag)[0] = std::min(1.0, Q_neg/P_neg);
            }

            if (fabs(P_pos) < VAREPS)
            {
                nodeit->RealArray(R_node_neg_pos_tag)[1] = 0.0;
            }
            else
            {
                nodeit->RealArray(R_node_neg_pos_tag)[1] = std::min(1.0, Q_pos/P_pos);
            }   
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(R_node_neg_pos_tag, NODE, 0);
    BARRIER

    // calculate trian_alpha
    
    for(Mesh::iteratorCell trianit = n.GetMesh()->BeginCell();
             trianit != n.GetMesh()->EndCell();
             ++trianit)
    {
        if(trianit->GetStatus() != Element::Ghost)
        {
            INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
            double Fe1 = trianit->RealArray(Fe_trian_tag)[0];
            double Fe2 = trianit->RealArray(Fe_trian_tag)[1];
            double Fe3 = trianit->RealArray(Fe_trian_tag)[2];

            double R1, R2, R3;
            
            if (Fe1 < 0)
            {
                R1 = local_nodes[0]->RealArray(R_node_neg_pos_tag)[0];
            }
            else
            {
                R1 = local_nodes[0]->RealArray(R_node_neg_pos_tag)[1];
            }
            
            if (Fe2 < 0)
            {
                R2 = local_nodes[1]->RealArray(R_node_neg_pos_tag)[0];
            }
            else
            {
                R2 = local_nodes[1]->RealArray(R_node_neg_pos_tag)[1];
            }

            if (Fe3 < 0)
            {
                R3 = local_nodes[2]->RealArray(R_node_neg_pos_tag)[0];
            }
            else
            {
                R3 = local_nodes[2]->RealArray(R_node_neg_pos_tag)[1];
            }

            trianit->Real(alpha_trian_tag) = std::min(std::min(R1, R2), R3); 
        }
    }
    BARRIER
    n.GetMesh()->ExchangeData(alpha_trian_tag, CELL, 0);
    BARRIER

    // delete all temporal tags
    n.GetMesh()->DeleteTag(m_node_str_tag, NODE);
    n.GetMesh()->DeleteTag(m_trian_str_str_tag, CELL);
    n.GetMesh()->DeleteTag(m_node_min_max_tag, NODE);
    n.GetMesh()->DeleteTag(q_node_min_max_tag, NODE);
    n.GetMesh()->DeleteTag(p_node_neg_pos_tag, NODE);
    n.GetMesh()->DeleteTag(R_node_neg_pos_tag, NODE);
    BARRIER
}