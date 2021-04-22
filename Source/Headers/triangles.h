#pragma once
#include "inmost_ice.h"

class INMOST_ICE_triangles
{
public:
    INMOST_ICE_triangles(INMOST_ICE_mesh& m);
    INMOST::Mesh* GetMesh();
    INMOST_ICE_mesh_data* GetMeshData();
private:
    INMOST::Mesh* ice_mesh;
    INMOST_ICE_mesh_data* data;
};