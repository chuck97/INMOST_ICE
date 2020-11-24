#include "config.h"
#include "INMOST_ICE_mesh.h"

int main(int argc, char **argv)
{
    double ttt = Timer();
#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
#endif
    // mesh initialization
    INMOST_ICE_mesh m((std::string)PMF_PATH);

    // nodes initialization
    INMOST_ICE_nodes n(m.GetMesh());

    // read data from .netcdf in parallel
    //n.ReadNodeData();

    // Interpolate a ice
    ttt = Timer();
    n.InterpolateScalarFromNETCDF(TOPAZ_NC_FILENAME,
                                  X_NAME_TOPAZ,
                                  Y_NAME_TOPAZ,
                                  0,
                                  "a",
                                  SCALE_FACTOR_NAME_TOPAZ,
                                  INVALID_VALUE_NAME_TOPAZ,
                                  OFFSET_NAME_TOPAZ,
                                  0.0,
                                  0.0,
                                  A_NAME_TOPAZ,
                                  1.1,
                                  false,
                                  n.GetCoords().topaz_coords
                                  );
    BARRIER;

    if (m.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "a ice interpolation: " << Timer() - ttt << std::endl;
    }

    // Interpolate h water
    ttt = Timer();
    n.InterpolateScalarFromNETCDF(TOPAZ_NC_FILENAME,
                                  X_NAME_TOPAZ,
                                  Y_NAME_TOPAZ,
                                  0,
                                  "hw",
                                  "",
                                  INVALID_VALUE_NAME_TOPAZ,
                                  "",
                                  0.0,
                                  0.0,
                                  HW_NAME_TOPAZ,
                                  3.0,
                                  false,
                                  n.GetCoords().topaz_coords
                                  );
    
    BARRIER;
    if (m.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "h water interpolation: " << Timer() - ttt << std::endl;
    }


    // Interpolate h ice
    ttt = Timer();
    n.InterpolateScalarFromNETCDF(TOPAZ_NC_FILENAME,
                                  X_NAME_TOPAZ,
                                  Y_NAME_TOPAZ,
                                  0,
                                  "h",
                                  SCALE_FACTOR_NAME_TOPAZ,
                                  INVALID_VALUE_NAME_TOPAZ,
                                  OFFSET_NAME_TOPAZ,
                                  0.0,
                                  0.0,
                                  H_NAME_TOPAZ,
                                  20.0,
                                  false,
                                  n.GetCoords().topaz_coords
                                  );
    
    BARRIER;  
    if (m.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "h ice interpolation: " << Timer() - ttt << std::endl;
    }             

    // interpolate u ice
    ttt = Timer();
    n.InterpolateVectorFromNETCDF(TOPAZ_NC_FILENAME,
                                  X_NAME_TOPAZ,
                                  Y_NAME_TOPAZ,
                                  "longitude",
                                  "latitude",
                                  0,
                                  "u",
                                  SCALE_FACTOR_NAME_TOPAZ,
                                  INVALID_VALUE_NAME_TOPAZ,
                                  OFFSET_NAME_TOPAZ,
                                  0.0,
                                  0.0,
                                  U_NAME_TOPAZ,
                                  V_NAME_TOPAZ,
                                  2.0,
                                  false,
                                  n.GetCoords().topaz_coords
                                  );
    BARRIER;
    if (m.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "u ice interpolation: " << Timer() - ttt << std::endl;
    } 

    // interpolate u water
    ttt = Timer();
    n.InterpolateVectorFromNETCDF(TOPAZ_NC_FILENAME,
                                  X_NAME_TOPAZ,
                                  Y_NAME_TOPAZ,
                                  "longitude",
                                  "latitude",
                                  0,
                                  "uw",
                                  SCALE_FACTOR_NAME_TOPAZ,
                                  INVALID_VALUE_NAME_TOPAZ,
                                  OFFSET_NAME_TOPAZ,
                                  0.0,
                                  0.0,
                                  UW_NAME_TOPAZ,
                                  VW_NAME_TOPAZ,
                                  10.0,
                                  true,
                                  n.GetCoords().topaz_coords
                                  );
    BARRIER;
    if (m.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "u water interpolation: " << Timer() - ttt << std::endl;
    } 
                                  
    m.PrintPVTK("test.pvtk");
   
    return 0;
}