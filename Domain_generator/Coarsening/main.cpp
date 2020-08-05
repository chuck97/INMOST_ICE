#include <iostream>
#include <string>
#include "contours_coarsening.h"

static const std::string RESOLUTION   = "crude";
static const int COARSENING_EXTERNAL_DEPTH = 2;
static const int COARSENING_ISLANDS_DEPTH = 6;

void external_coarsening()
{
  std::cout << std::endl;
  std::cout << "===========================" << std::endl;
  std::cout << "=== external coarsening ===" << std::endl;
  std::cout << "===========================" << std::endl;
  std::cout << std::endl;
  
  
  std::string input_path = "../../Domain/data/arctic_" + RESOLUTION +"_external.txt";
  std::string output_path = "../../Domain/data/arctic_" + RESOLUTION +"_coarsened_external.txt";


  Contour_Database<double> db(input_path, output_path, true);
  
  std::cout << "Total points before procedure: " << db.Size() << std::endl;
  std::cout << std::endl;
  
  db.RemoveZeroEdges();
  std::cout << std::endl;
  db.RemoveCuts();
  std::cout << std::endl;
  db.RemoveZeroEdges();
  std::cout << std::endl;
  db.RemoveCuts();
  std::cout << std::endl;
  db.RemoveZeroEdges();
  std::cout << std::endl;
  
  for (int i = 0 ; i < COARSENING_EXTERNAL_DEPTH; ++i)
  {
	std::cout << "#### Coarsening iteration " << i << " ####" << std::endl;
    db.Coarsening();
    db.RemoveZeroEdges();
    db.RemoveCuts();
    db.RemoveZeroEdges();
    db.RemoveCuts();
    db.RemoveZeroEdges();
    std::cout << std::endl;
  }
  
  std::cout << "Total points after procedure: " << db.Size() << std::endl;
  db.WriteToFile();
  std::system("python ./plot_external_coarsened.py");
}

void islands_coarsening()
{
	
  std::cout << std::endl;
  std::cout << "===========================" << std::endl;
  std::cout << "===  islands coarsening ===" << std::endl;
  std::cout << "===========================" << std::endl;
  std::cout << std::endl;
   
  std::string input_path = "../../Domain/data/arctic_" + RESOLUTION +"_domain_islands.txt";
  std::string output_path = "../../Domain/data/arctic_" + RESOLUTION +"_coarsened_islands.txt";
  Contour_Database<double> db(input_path, output_path, false);
  
  std::cout << "Total points before procedure: " << db.Size() << std::endl;
  std::cout << std::endl;
  
  std::cout << "Total islands before procedure: " << db.N_islands() << std::endl;
  std::cout << std::endl;
  
  db.RemoveZeroEdges();
  std::cout << std::endl;
  db.RemoveCuts();
  std::cout << std::endl;
  db.RemoveZeroEdges();
  std::cout << std::endl;
  db.RemoveCuts();
  std::cout << std::endl;
  db.RemoveZeroEdges();
  std::cout << std::endl;
  
  for (int i = 0 ; i < COARSENING_ISLANDS_DEPTH; ++i)
  {
	std::cout << "#### Coarsening iteration " << i << " ####" << std::endl;
    db.Coarsening();
    db.RemoveZeroEdges();
    db.RemoveCuts();
    db.RemoveZeroEdges();
    db.RemoveCuts();
    db.RemoveZeroEdges();
    std::cout << std::endl;
  }
 
  std::cout << "Total points after procedure: " << db.Size() << std::endl;
  db.WriteToFile();
  std::system("python ./plot_islands_coarsened.py");
}

int main()
{
  external_coarsening();
  islands_coarsening();
  // plot external boundary with islands
  std::system("python ./plot_full_coarsened.py");
}
