#include "config.h"
#include "INMOST_ICE_mesh.h"

void FindExtremalTopaz(INMOST_ICE_nodes& n)
{
    double ttt = Timer();
    n.FindExtrCoords();
    if (n.GetMesh()->GetProcessorRank() == 0)
    {
        std::cout << "Finding extremal topaz coords: " << Timer() - ttt << std::endl;
    }
}

int main(int argc, char **argv)
{
#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
#endif
    // mesh initialization
    INMOST_ICE_mesh m((std::string)PMF_PATH);

    // nodes initialization
    INMOST_ICE_nodes n(m.GetMesh());
    
    // finding extremal TOPAZ coords
    FindExtremalTopaz(n);

    m.PrintPVTK("./test.pvtk");
    return 0;
}