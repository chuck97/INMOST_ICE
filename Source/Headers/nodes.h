#pragma once
#include "inmost_ice.h"

class INMOST_ICE_nodes
{
public:
    INMOST_ICE_nodes(INMOST_ICE_mesh& m);
    INMOST::Mesh* GetMesh();
    NodeCoords* GetCoords();
    INMOST_ICE_mesh_data* GetMeshData();
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
    INMOST::Mesh* ice_mesh;
    INMOST_ICE_mesh_data* data;
    NodeCoords nc;
};