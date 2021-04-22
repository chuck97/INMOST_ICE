#pragma once
#include "config.h"
#include "INMOST_ICE_mesh.h"
#include "MatVec.h"

void ltrim(std::string &s);
void rtrim(std::string &s);
void trim(std::string &s);

struct NodeCoords
{
    INMOST::Tag model_coords;
    INMOST::Tag geo_coords;
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
    void PrintPVTU(const std::string& filename);    
    void PrintNETCDF(const std::string& filename);  
private:
    INMOST::Mesh* ice_mesh;
    INMOST_ICE_mesh_data data;
};

class INMOST_ICE_nodes
{
public:
    INMOST_ICE_nodes(INMOST::Mesh* m); 
    NodeCoords& GetCoords();
    INMOST::Mesh* GetMesh();
    unsigned int GetIdMin() const;
    unsigned int GetIdMax() const;
    double GetResolution() const;
    void PrintPVTU(const std::string& filename); 
    std::map<std::string, INMOST::Tag>& GetData();
    INMOST::Tag GetLocalNodeCoords();
    INMOST::Tag GetTransitionMatricies();
private:
    INMOST::Mesh* ice_mesh;
    NodeCoords nc;
    std::map<std::string, INMOST::Tag> node_data_tags;
    unsigned int idmin;
    unsigned int idmax;
    double spatial_resolution;
    INMOST::Tag u_transition_matricies;
    INMOST::Tag local_node_coords;
};