#include <iostream>
#include <string>
#include "contours_coarsening.h"

const std::string RESOLUTION   = "crude";
const int COARSENING_DEPTH = 2;

int main()
{
  std::string input_path = "../../Domain/data/arctic_" + RESOLUTION +"_external.txt";
  std::string output_path = "../../Domain/data/arctic_" + RESOLUTION +"_coarsened_external.txt";


  Contour_Database<double> db(input_path, output_path);
  db.RemoveZeoroEdges();
  std::cout << std::endl;
  db.RemoveCuts();
  std::cout << std::endl;
  db.RemoveZeoroEdges();
  std::cout << std::endl;
  db.RemoveCuts();
  std::cout << std::endl;
  db.RemoveZeoroEdges();
  std::cout << std::endl;
  
  std::cout << std::endl;
  for (int i = 0 ; i < COARSENING_DEPTH; ++i)
  {
	std::cout << "#### Coarsening iteration " << i << " ####" << std::endl;
    db.Coarsening();
    db.RemoveZeoroEdges();
    db.RemoveCuts();
    db.RemoveZeoroEdges();
    db.RemoveCuts();
    db.RemoveZeoroEdges();
    std::cout << std::endl;
  }
  std::cout << "countour size: " << db.Size() << std::endl;
  db.WriteToFile();
  std::system("python ./plot_external_coarsened.py");
}
