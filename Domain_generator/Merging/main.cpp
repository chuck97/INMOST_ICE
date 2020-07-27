#include "contours_merging.h"

const std::string RESOLUTION = "crude";
const bool ISLANDS           = false;

int main()
{
  Euler_rotation_info<double> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF);

  std::string input_path = "../Coastline/arctic_" + RESOLUTION + ".txt";
  std::string output_path_islands = "../Domain/data/arctic_" + RESOLUTION + "_domain_islands.txt";
  std::string output_path_external = "../Domain/arctic_" + RESOLUTION + "_domain_external.txt";
  if (ISLANDS)
  {
    Contour_Database<double> db(input_path, output_path_islands, rotation.Get_REVERSE());
    db.MakeClosedDomain();
    db.WriteToFile(ISLANDS);
  }
  else
  {
    Contour_Database<double> db(input_path, output_path_external, rotation.Get_REVERSE());
    db.MakeClosedDomain();
    db.WriteToFile(ISLANDS);
  }
}
