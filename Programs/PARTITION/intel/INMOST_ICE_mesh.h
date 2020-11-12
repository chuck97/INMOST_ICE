#pragma once
#include "config.h"

static inline void ltrim(std::string &s);
static inline void rtrim(std::string &s);
static inline void trim(std::string &s);

class INMOST_ICE_mesh
{
public:
    INMOST_ICE_mesh(const std::string& mesh_path);  // +
    ~INMOST_ICE_mesh();                             // +
    void Partition();                               // +
    INMOST::Mesh* GetMesh();                        // +
    void PrintPMF(const std::string& filename);     // +
    void PrintPVTK(const std::string& filename);    // +
    void PrintNETCDF(const std::string& filename);  // -
private:
    INMOST::Mesh* ice_mesh;
};

struct NodeCoords
{
    INMOST::Tag model_coords;
    INMOST::Tag geo_coords;
    INMOST::Tag topaz_coords;
};

class INMOST_ICE_nodes
{
public:
    INMOST_ICE_nodes(INMOST::Mesh* m);
    INMOST::Mesh* GetMesh();
    void FindExtrCoords();
    NodeCoords& GetCoords();
    void ReadNodeData();                                    
private:
    INMOST::Mesh* ice_mesh;
    NodeCoords nc;
    INMOST::Tag a;
public:
    double min_topaz_x, max_topaz_x, min_topaz_y, max_topaz_y;
};
