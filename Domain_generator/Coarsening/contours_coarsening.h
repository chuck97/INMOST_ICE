#pragma once
#include "coords_rotation.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <list>
#include <stdexcept>
#include <algorithm>
#include <math.h>
#include <iterator>

const double MINIMAL_ANGLE     = 20.0; //in degrees   // 15.0
const double ZERO_EDGE_SIZE    = 26e-2;
const double MINIMAL_DISTANCE  = 61364e-6;         	  // 61364e-6
const int    FORWARD_COUNTER   = 6;                   // 6

template<typename RealType>
RealType dist(Spherical_Coords<RealType> first,
              Spherical_Coords<RealType> second)
{
  return sqrt((second.Get_x() - first.Get_x())*(second.Get_x() - first.Get_x()) +
              (second.Get_y() - first.Get_y())*(second.Get_y() - first.Get_y()));
}

template<typename RealType>
bool check_collinearity(Spherical_Coords<RealType>& first,
                        Spherical_Coords<RealType>& second,
                        Spherical_Coords<RealType>& third)
{
  RealType a = dist(first, second);
  RealType b = dist(second, third);
  RealType c = dist(first, third);
  RealType cos_alpha = (a*a + b*b - c*c)/(2.0*a*b);
  RealType alpha = acos(cos_alpha)*180.0/M_PI;
  if (alpha < MINIMAL_ANGLE)
  {
    return true;
  }
  return false;
}

template<typename RealType>
bool check_neighbourhood(Spherical_Coords<RealType>& target_point,
                         Spherical_Coords<RealType>& s_begin,
                         Spherical_Coords<RealType>& s_end)
{
  RealType a = dist(s_begin, s_end);
  RealType b = dist(target_point, s_begin);
  RealType c = dist(target_point, s_end);
  RealType p = (a + b + c)/2.0;
  RealType S = sqrt(p*(p-a)*(p-b)*(p-c));
  RealType h = 2.0*S/a;
  if (h < MINIMAL_DISTANCE)
  {
    return true;
  }
  return false;
}


template<typename RealType>
class Contour_Database
{
public:
  Contour_Database(const std::string& input_file_path,
                   const std::string& output_file_path)
    : input_file_path_(input_file_path)
    , output_file_path_(output_file_path)
  {
    ReadFromFile();
  }

  void ReadFromFile()
  {
    std::ifstream input(input_file_path_);
    std::string line;

    //skip first line
    getline(input, line);

    // Read all spherical Points of external_boundary_
    while (getline(input, line))
    {
      std::stringstream ss(line);
      RealType first, second;
      ss >> first >> second;
      external_boundary_.push_back({first, second});
    }
  }

  void WriteToFile()
  {
    std::ofstream output(output_file_path_);

    output << "# external boundary " << std::endl;
    for (auto& point: external_boundary_)
    {
      output << point << std::endl;
    }
  }

  void RemoveCuts()
  {
    size_t before_remove_cuts = external_boundary_.size();
    auto current_point_iter = external_boundary_.begin();
    while(*current_point_iter != *(++external_boundary_.rbegin()))
    {
      auto next = ++current_point_iter;
      --current_point_iter;
      auto next_next = ++next;
      --next;

      if (check_collinearity<RealType>(*current_point_iter, *next, *next_next))
      {
        external_boundary_.erase(next);
      }
      current_point_iter++;
    }
    size_t after_remove_cuts = external_boundary_.size();
    std::cout <<"Number of removed points due to cut suspicion: " <<
      before_remove_cuts - after_remove_cuts << std::endl;
  }

  void Coarsening()
  {
    size_t before_coarsening = external_boundary_.size();

    auto final_point_iter = external_boundary_.rbegin();
    for (int i = 0; i < FORWARD_COUNTER; ++i)
    {
      ++final_point_iter;
    }

    auto current_point_iter = external_boundary_.begin();

    while(*current_point_iter != *(final_point_iter))
    {
      auto last_erase_it = current_point_iter;
      auto segment_begin = ++current_point_iter;
      --current_point_iter;
      for(int i = 0; i < FORWARD_COUNTER; ++i)
      {
        auto segment_end = ++segment_begin;
        --segment_begin;

        if (check_neighbourhood(*current_point_iter, *segment_begin, *segment_end))
        {
          last_erase_it = segment_begin;
        }
        ++segment_begin;
      }

      if (last_erase_it == current_point_iter)
      {
        current_point_iter++;
      }
      else
      {
        auto first_erase_iter = ++current_point_iter;
        current_point_iter = ++last_erase_it;
        external_boundary_.erase(first_erase_iter, last_erase_it);
      }
    }

    size_t after_coarsening = external_boundary_.size();
    std::cout <<"Number of removed points due to coarsening: " <<
        before_coarsening - after_coarsening << std::endl;
  }
	
  void RemoveZeoroEdges()
  {
	size_t before_zero = external_boundary_.size();
	auto current_point_iter = external_boundary_.begin();
	auto last_erase_it = (++(++external_boundary_.rbegin()));
	
	while (*current_point_iter != *last_erase_it)
	{
	  auto next = ++current_point_iter;
	  --current_point_iter;
	  if (dist(*current_point_iter, *next) <  ZERO_EDGE_SIZE)
	  {
	    current_point_iter = ++next;
	    --next;
	    external_boundary_.erase(next);
	  }
	  else
	  {
	    current_point_iter++; 
      }
	}
	
	size_t after_zero = external_boundary_.size();
	std::cout <<"Number of removed points due to zero size suspicion: " <<
        before_zero - after_zero << std::endl;
  }
  
  size_t Size()
  {
    return external_boundary_.size();
  }

private:
  const std::string      input_file_path_;
  const std::string      output_file_path_;
  std::list<Spherical_Coords<RealType>> external_boundary_;
};
