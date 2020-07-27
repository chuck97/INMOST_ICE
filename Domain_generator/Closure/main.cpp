#include "coords_rotation.h"
#include <string>
#include <sstream>
#include <math.h>
#include <utility>
#include <iostream>
#include <fstream>
#include <tuple>

const std::string RESOLUTION             = "crude";
const double      DESIRED_TRIANGLE_SIZE  = 1000.0;  // in km
const double      EARTH_RADIUS           = 6371.0;  // in km
const double      EXTREME_LATITUDE       = 45.0;    // extreme latitude while cropping
const double      EXTREME_LATITUDE_OFFSET = 5.0;   // in degrees


template<typename RealType>
Spherical_Coords<RealType> ReadCoords(std::string& line)
{
  std::stringstream ss(line);
  double first, second;
  ss >> first >> second;
  return Spherical_Coords<RealType>{first, second};
}

int main()
{
  Euler_rotation_info<double> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF);

  std::string input_path = "../../Domain/data/arctic_" + RESOLUTION + "_domain_external_parts.txt";
  std::string output_path = "../../Domain/data/arctic_" + RESOLUTION + "_external.txt";

  std::ifstream input(input_path);
  std::ofstream output(output_path);
  std::string line;

  output << "#external boundary" << std::endl;
  getline(input, line);
  Spherical_Coords<double> sph_rotated_east;
  Spherical_Coords<double> sph_rotated_east_first;
  Spherical_Coords<double> sph_east_first;

  //store first east point in geographical coordinates
  getline(input, line);
  sph_rotated_east = ReadCoords<double>(line);
  sph_rotated_east_first = sph_rotated_east;
  output << sph_rotated_east << std::endl;
  sph_east_first = Rotate_Spherical(sph_rotated_east_first, rotation.Get_FORWARD());

  // write east coast
  while (true)
  {
    getline(input, line);
    if (line[0] == '#')
    {
      break;
    }
    sph_rotated_east = ReadCoords<double>(line);
    output << sph_rotated_east << std::endl;
  }

  // store last east point in geographical coordinates
  Spherical_Coords<double> sph_east_last = Rotate_Spherical(sph_rotated_east, rotation.Get_FORWARD());

  // close Atlantic region
  getline(input, line);
  Spherical_Coords<double> sph_rotated_west = ReadCoords<double>(line);

  //store first west point in geographical coordinates
  Spherical_Coords<double> sph_west_first =  Rotate_Spherical(sph_rotated_west, rotation.Get_FORWARD());

  double east_longitude =  sph_east_last.Get_x();
  double west_longitude =  sph_west_first.Get_x();
  double atlantic_latitude = std::max(sph_east_first.Get_y(), sph_east_last.Get_y()) - EXTREME_LATITUDE_OFFSET;
  double degree_step_size = DESIRED_TRIANGLE_SIZE/(EARTH_RADIUS*cos(atlantic_latitude*M_PI/180.0));

  double current_longitude = east_longitude;
  while (current_longitude > west_longitude)
  {
    current_longitude -= degree_step_size;
    Spherical_Coords<double> sph_rotated_atlantic(current_longitude, atlantic_latitude);
    output << Rotate_Spherical(sph_rotated_atlantic, rotation.Get_REVERSE()) << std::endl;
  }

  // write west coast
  output << sph_rotated_west << std::endl;

  while (getline(input, line))
  {
    sph_rotated_west = ReadCoords<double>(line);
    output << sph_rotated_west << std::endl;
  }

  // close Pacific region
  east_longitude =  sph_east_first.Get_x();
  west_longitude =  Rotate_Spherical(sph_rotated_west, rotation.Get_FORWARD()).Get_x();
  double pacific_latitude = std::max(Rotate_Spherical(sph_rotated_west, rotation.Get_FORWARD()).Get_y(),
                                     sph_east_first.Get_y()) - EXTREME_LATITUDE_OFFSET;
  degree_step_size = DESIRED_TRIANGLE_SIZE/(EARTH_RADIUS*cos(pacific_latitude*M_PI/180.0));
  current_longitude = west_longitude;
  double current_longitude_tmp = current_longitude;

  while (current_longitude > east_longitude - 360.0)
  {
    current_longitude -= degree_step_size;
    if (current_longitude < -180.0)
    {
      current_longitude_tmp = current_longitude + 360.0;
    }
    else
    {
      current_longitude_tmp = current_longitude;
    }

    Spherical_Coords<double> sph_rotated_pacific(current_longitude_tmp, pacific_latitude);
    output << Rotate_Spherical(sph_rotated_pacific, rotation.Get_REVERSE()) << std::endl;
  }
  output << Rotate_Spherical(sph_east_first, rotation.Get_REVERSE()) << std::endl;
}
