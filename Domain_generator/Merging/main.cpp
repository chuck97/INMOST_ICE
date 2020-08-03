#include "contours_merging.h"
#include <iostream>

const std::string RESOLUTION = "crude";

int main()
{
  Euler_rotation_info<double> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF);

  std::string input_path = "../../Coastline/arctic_" + RESOLUTION + ".txt";
  std::string output_path_islands = "../../Domain/data/arctic_" + RESOLUTION + "_domain_islands.txt";
  std::string output_path_external = "../../Domain/data/arctic_" + RESOLUTION + "_domain_external_parts.txt";
  
  //Islands
  {
    Contour_Database<double> db(input_path, output_path_islands, rotation.Get_REVERSE());
    db.MakeClosedDomain();
    db.WriteToFile(true);
    std::system("python ./plot_islands.py");
  }
  
  //External boundary parts
  {
	Contour_Database<double> db(input_path, output_path_external, rotation.Get_REVERSE());
	db.MakeClosedDomain();
    db.WriteToFile(false);
    std::system("python ./plot_external_parts.py");
  }
}
