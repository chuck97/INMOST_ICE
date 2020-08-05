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

const double MINIMAL_ANGLE_EXTERNAL    = 20.0; //in degrees   // 15.0
const double MINIMAL_ANGLE_ISLANDS     = 30.0;
const double ZERO_EDGE_SIZE_EXTERNAL   = 26e-2;
const double ZERO_EDGE_SIZE_ISLANDS    = 26e-2;
const double MINIMAL_DISTANCE_EXTERNAL = 61364e-6;            // 61364e-6
const double MINIMAL_DISTANCE_ISLANDS  = 61364e-6;            // 61364e-6
const size_t FORWARD_EXTERNAL_COUNTER  = 6;                   // 6
const size_t FORWARD_ISLANDS_COUNTER   = 6;                   // 6

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
                        Spherical_Coords<RealType>& third,
                        bool is_external)
{
  double MINIMAL_ANGLE;
  if (is_external)
  {
    MINIMAL_ANGLE = MINIMAL_ANGLE_EXTERNAL;
  }
  else
  {
    MINIMAL_ANGLE = MINIMAL_ANGLE_ISLANDS;
  }
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
                         Spherical_Coords<RealType>& s_end,
                         bool is_external)
{
  double MINIMAL_DISTANCE;
  if (is_external)
  {
	  MINIMAL_DISTANCE = MINIMAL_DISTANCE_EXTERNAL;
  }
  else
  {
	  MINIMAL_DISTANCE = MINIMAL_DISTANCE_ISLANDS;
  }
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
                   const std::string& output_file_path,
                   bool is_external)
    : input_file_path_(input_file_path)
    , output_file_path_(output_file_path)
    , is_external_(is_external)
  {
    ReadFromFile();
  }

  void ReadFromFile()
  {
    std::ifstream input(input_file_path_);
    std::string line;

    std::list<Spherical_Coords<RealType>> current_island;

    //skip first line
    getline(input, line);

    // Read all spherical Points of boundary_
    while (getline(input, line))
    {
      if (line[0] == '#')
      {
        boundary_.push_back(current_island);
        current_island.clear();
      }
      else
      {
        std::stringstream ss(line);
        RealType first, second;
        ss >> first >> second;
        current_island.push_back({first, second});
      }
    }

    if (current_island.size() > 0)
    {
        boundary_.push_back(current_island);
    }
  }

  void WriteToFile()
  {
    std::ofstream output(output_file_path_);

    if (boundary_.size() == 1)
    {
        output << "# external boundary " << std::endl;
        for (auto& point: boundary_.back())
        {
            output << point << std::endl;
        }
    }
    else
    {
        size_t i = 0;
        for (auto& island: boundary_)
        {
            output << "# island " << i << std::endl;
            for (auto& point: island)
            {
                output << point << std::endl;
            }
            ++i;
        }
    }
  }

  void RemoveCuts()
  {
    size_t before_remove_cuts = Size();
    auto current_island_iter = boundary_.begin();
    while (current_island_iter != boundary_.end())
    {
        auto current_point_iter = (*current_island_iter).begin();

        if ((*current_island_iter).size() > 4)
        {
            while((*current_point_iter != *(++((*current_island_iter).rbegin()))) and
                  ((*current_point_iter != *((*current_island_iter).rbegin()))))
            {
                auto next = ++current_point_iter;
                --current_point_iter;
                auto next_next = ++next;
                --next;
                if (check_collinearity<RealType>(*current_point_iter, *next, *next_next, is_external_))
                {
                  current_point_iter = next_next;
                  (*current_island_iter).erase(next);
                }
                else
                {
                    current_point_iter = next;
                }
            }
        }

        if ((*current_island_iter).size() <= 4)
        {
            auto remove_island_iter = current_island_iter;
            ++current_island_iter;
            boundary_.erase(remove_island_iter);
        }
        else
        {
            ++current_island_iter;
        }
    }
    size_t after_remove_cuts = Size();
    std::cout <<"Number of removed points due to cut suspicion: " <<
    before_remove_cuts - after_remove_cuts << std::endl;
  }

  void Coarsening()
  {
    size_t before_coarsening = Size();
    auto current_island_iter = boundary_.begin();
    size_t island_num = 0;
    size_t FORWARD_COUNTER;
		  
	if (is_external_)
	{
		FORWARD_COUNTER = FORWARD_EXTERNAL_COUNTER;
	}
	else
	{
		FORWARD_COUNTER = FORWARD_ISLANDS_COUNTER;
	}
    
    while (current_island_iter != boundary_.end())
    {
		size_t current_island_size = (*current_island_iter).size();
		if ((*current_island_iter).size() < FORWARD_ISLANDS_COUNTER + 3)
		{
			++current_island_iter;
			++island_num;
			continue;
		} 
		   
        auto current_point_iter = (*current_island_iter).begin();
		size_t current_point_num = 0;
		size_t max_point_num = current_island_size - FORWARD_COUNTER;
		size_t shift = 0;
		
        while(current_point_num < max_point_num) 
        {
            auto last_erase_it = current_point_iter;
            auto segment_begin = ++current_point_iter;
            --current_point_iter;
            size_t i = 0;
            for(int i = 0; i < FORWARD_COUNTER; ++i)
            {
                auto segment_end = ++segment_begin;
                --segment_begin;
                ++i;

                if (check_neighbourhood(*current_point_iter, *segment_begin, *segment_end, is_external_))
                {
                    last_erase_it = segment_begin;
                    shift = i;
                }
                ++segment_begin;
            }

            if (last_erase_it == current_point_iter)
            {
                current_point_iter++;
                ++current_point_num;
            }
            else
            {
				current_point_num += shift;
                auto first_erase_iter = ++current_point_iter;
				current_point_iter = ++last_erase_it; 
				(*current_island_iter).erase(first_erase_iter, last_erase_it);
            }
        }

        if ((*current_island_iter).size() <= 4)
        {
            auto remove_island_iter = current_island_iter++;
            boundary_.erase(remove_island_iter);
        }
        else
        {
			++current_island_iter;
		}
    }
    
    size_t after_coarsening = Size();
    std::cout <<"Number of removed points due to coarsening: " <<
        before_coarsening - after_coarsening << std::endl;
  }

  void RemoveZeroEdges()
  {
    size_t before_zero = Size();
    double ZERO_EDGE_SIZE;
    if (is_external_)
    {
		ZERO_EDGE_SIZE = ZERO_EDGE_SIZE_EXTERNAL;
	}
	else
	{
		ZERO_EDGE_SIZE = ZERO_EDGE_SIZE_ISLANDS;
	}
	
    auto current_island_iter = boundary_.begin();
    while (current_island_iter != boundary_.end())
    {
        auto current_point_iter = (*current_island_iter).begin();
        auto last_erase_it = (++(++(*current_island_iter).rbegin()));

        while (*current_point_iter != *last_erase_it)
        {
            auto next = ++current_point_iter;
            --current_point_iter;
            if (dist(*current_point_iter, *next) <  ZERO_EDGE_SIZE)
            {
                current_point_iter = ++next;
                --next;
                (*current_island_iter).erase(next);
            }
            else
            {
                current_point_iter++;
            }
        }

        if ((*current_island_iter).size() <= 4)
        {
            auto remove_island_iter = current_island_iter;
            ++current_island_iter;
            boundary_.erase(remove_island_iter);
        }
        else
        {
            ++current_island_iter;
        }
    }

    size_t after_zero = Size();
    std::cout <<"Number of removed points due to zero size suspicion: " <<
        before_zero - after_zero << std::endl;
  }

  size_t Size()
  {
    size_t total = 0;
    for (auto& island: boundary_)
    {
        total += island.size();
    }
    return total;
  }

  size_t N_islands()
  {
    return boundary_.size();
  }

private:
  const std::string                                input_file_path_;
  const std::string                                output_file_path_;
  bool                                             is_external_;
  std::list<std::list<Spherical_Coords<RealType>>> boundary_;
};
