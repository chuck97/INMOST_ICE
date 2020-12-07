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
    NodeCoords& GetCoords();     
    void TopazScalarInterpolation(const std::string& filename,
                                     const std::string& xname,
                                     const std::string& yname,
                                     int netcdf_time_step,
                                     const std::string& node_variable_name,
                                     const std::string& scale_factor_name,
                                     const std::string& invalid_value_name,
                                     const std::string& offset_name,
                                     double invalid_value_fill,
                                     double no_extrapolation_fill,
                                     const std::string& nc_variable_name,
                                     double max_abs_value,
                                     bool is_depth,
                                     INMOST::Tag netcdf_coords
                                     );

    void TopazVectorInterpolation(const std::string& filename,
                                     const std::string& xname,
                                     const std::string& yname,
                                     const std::string& lonname,
                                     const std::string& latname,
                                     int netcdf_time_step,
                                     const std::string& node_variable_name,
                                     const std::string& scale_factor_name,
                                     const std::string& invalid_value_name,
                                     const std::string& offset_name,
                                     double invalid_value_fill,
                                     double no_extrapolation_fill,
                                     const std::string& nc_variable_name1,
                                     const std::string& nc_variable_name2,
                                     double max_abs_value,
                                     bool is_depth,
                                     INMOST::Tag netcdf_coords
                                     );
    
    void CamsScalarInterpolation(const std::string& filename,
                                     const std::string& lonname,
                                     const std::string& latname,
                                     int netcdf_time_step,
                                     const std::string& node_variable_name,
                                     const std::string& scale_factor_name,
                                     const std::string& invalid_value_name,
                                     const std::string& offset_name,
                                     double invalid_value_fill,
                                     double no_extrapolation_fill,
                                     const std::string& nc_variable_name,
                                     double max_abs_value,
                                     bool is_depth,
                                     INMOST::Tag netcdf_coords
                                     );

    void CamsVectorInterpolation(const std::string& filename,
                                     const std::string& lonname,
                                     const std::string& latname,
                                     int netcdf_time_step,
                                     const std::string& node_variable_name,
                                     const std::string& scale_factor_name,
                                     const std::string& invalid_value_name,
                                     const std::string& offset_name,
                                     double invalid_value_fill,
                                     double no_extrapolation_fill,
                                     const std::string& nc_variable_name1,
                                     const std::string& nc_variable_name2,
                                     double max_abs_value,
                                     bool is_depth,
                                     INMOST::Tag netcdf_coords
                                     );
private:
    std::map<std::string, INMOST::Tag> node_data_tags;
    INMOST::Mesh* ice_mesh;
    NodeCoords nc;
    INMOST::Tag a;
};