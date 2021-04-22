#pragma once
#include "config.h"

struct NodeCoords
{
    INMOST::Tag model_coords;
    INMOST::Tag geo_coords;
    INMOST::Tag topaz_coords;
};

struct INMOST_ICE_mesh_data
{
    std::map<std::string, INMOST::Tag> node_data;
    std::map<std::string, INMOST::Tag> trian_data;
    NodeCoords node_coords;
};

class INMOST_ICE_mesh
{
public:
    INMOST_ICE_mesh(const std::string& mesh_path);  
    ~INMOST_ICE_mesh();                             
    void Partition();                               
    INMOST::Mesh* GetMesh();           
    INMOST_ICE_mesh_data* GetMeshData();             
    void PrintPMF(const std::string& filename);     
    void PrintPVTK(const std::string& filename);    
    void PrintNETCDF(const std::string& filename);  
private:
    INMOST::Mesh* ice_mesh;
    INMOST_ICE_mesh_data data;
};