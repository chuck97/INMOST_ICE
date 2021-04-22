#include "inmost_ice.h"

using namespace INMOST;

INMOST_ICE_triangles::INMOST_ICE_triangles(INMOST_ICE_mesh& m)
{
    ice_mesh = m.GetMesh();
    data = m.GetMeshData();

    double ttt = Timer();
    for(Mesh::iteratorCell trianit = ice_mesh->BeginCell(); trianit != ice_mesh->EndCell(); ++trianit) 
	{
        data = m.GetMeshData();
        ice_mesh = m.GetMesh();

        ElementArray<Node> trnodes = trianit->getNodes();

        double x0 = trnodes[0].RealArray(ice_mesh->Coords[0]);
        double y0 = trnodes[0].RealArray(ice_mesh->Coords[1]);
        double x1 = trnodes[1].RealArray(ice_mesh->Coords[0]);
        double y1 = trnodes[1].RealArray(ice_mesh->Coords[1]);
        double x2 = trnodes[2].RealArray(ice_mesh->Coords[0]);
        double y2 = trnodes[2].RealArray(ice_mesh->Coords[1]);

        INMOST::Tag t;
        t = ice_mesh->CreateTag("trsize", DATA_REAL, CELL, NONE, 1);
        data.trian_data["trsize"] = trsize({x0, y0},
                                           {x1, y1},
                                           {x2, y2}); 
    }
    if (ice_mesh->GetProcessorRank() == 0)
    {
        std::cout << "Initializing triangles: " << Timer() - ttt << std::endl;
    }
}

INMOST::Mesh* INMOST_ICE_triangles::GetMesh()
{
    return ice_mesh;
}

INMOST_ICE_mesh_data* INMOST_ICE_triangles::GetMeshData()
{
    return data;
}