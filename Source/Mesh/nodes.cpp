#include "inmost_ice.h"

using namespace INMOST;

INMOST_ICE_nodes::INMOST_ICE_nodes(INMOST_ICE_mesh& m)
{
    ice_mesh = m.GetMesh();
    data = m.GetMeshData();
    double ttt = Timer();

    nc.model_coords = ice_mesh->CreateTag("model coords", DATA_REAL, NODE, NONE, 2);
    nc.geo_coords = ice_mesh->CreateTag("geo coords", DATA_REAL, NODE, NONE, 2);
    nc.topaz_coords = ice_mesh->CreateTag("TOPAZ coords", DATA_REAL, NODE, NONE, 2);    

    for(Mesh::iteratorNode nodeit = ice_mesh->BeginNode(); nodeit != ice_mesh->EndNode(); ++nodeit) 
	{
        // Assign model coords
		nodeit->RealArray(nc.model_coords)[0] = nodeit->Coords()[0];
        nodeit->RealArray(nc.model_coords)[1] = nodeit->Coords()[1];

        
        // Assign geo coords
        std::vector<double> geo = from_model_2_geo<double>(nodeit->RealArray(nc.model_coords)[0],
                                                           nodeit->RealArray(nc.model_coords)[1]);
        nodeit->RealArray(nc.geo_coords)[0] = geo[0];
        nodeit->RealArray(nc.geo_coords)[1] = geo[1];
        
        // Assign TOPAZ coords
        std::vector<double> topaz = from_geo_2_topaz(nodeit->RealArray(nc.geo_coords)[0],
                                                     nodeit->RealArray(nc.geo_coords)[1]);
        nodeit->RealArray(nc.topaz_coords)[0] = topaz[0]/1e5;
        nodeit->RealArray(nc.topaz_coords)[1] = topaz[1]/1e5;
    }
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Initializing nodes: " <<  Timer() - ttt <<std::endl;
    }
};

NodeCoords& INMOST_ICE_nodes::GetCoords()
{
    return nc;
};

INMOST::Mesh* INMOST_ICE_nodes::GetMesh()
{
    return ice_mesh;
};

INMOST_ICE_mesh_data* INMOST_ICE_nodes::GetMeshData()
{
    return data;
}