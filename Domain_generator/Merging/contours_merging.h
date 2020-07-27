// This header file contains functions and classes for making closed contours from .txt files generated fom GMT:
// this code uses "coords_rotation.h"
//
//  1)
//  2)
//  3)
//  4)
//  5)
//  6)
//  7)
//  8)
//  9)
// 10)

#pragma once
#include "coords_rotation.h"

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <list>
#include <stdexcept>
#include <algorithm>


const double EXTREME_LATITUDE = 45.0;             // extreme latitude while cropping
const double DESIRED_TRIANGULATION_STEP = 10.0;   // desired triangle size in km

template<typename RealType>
class ContourPart
{
public:

  void AddPoint(const Spherical_Coords<RealType>& new_point)
  {
    if (contour_part_.empty())
    {
      contour_part_.push_back(new_point);
    }
    else if (contour_part_.front() == new_point)
    {
      contour_part_.push_back(new_point);
      is_closed_ = true;
    }
    else
    {
      contour_part_.push_back(new_point);
    }
  }

  std::list<Spherical_Coords<RealType>>& GetContour()
  {
    return contour_part_;
  }

  bool is_closed() const
  {
    return is_closed_;
  }

  Spherical_Coords<RealType> First() const
  {
    return contour_part_.front();
  }

  Spherical_Coords<RealType> Last() const
  {
    return contour_part_.back();
  }

  int Size() const
  {
    return contour_part_.size();
  }

  void Clear()
  {
    contour_part_.clear();
    is_closed_ = false;
  }

  void EraseLast()
  {
    contour_part_.pop_back();
  }

  const typename std::list<Spherical_Coords<RealType>>::iterator begin()
  {
    return contour_part_.begin();
  }

  const typename std::list<Spherical_Coords<RealType>>::iterator end()
  {
    return contour_part_.end();
  }

  void Merge(ContourPart<RealType>& other)
  {
    if (contour_part_.back() == other.First())
    {
      contour_part_.pop_back();
      std::copy
      (
          other.GetContour().begin(),
          other.GetContour().end(),
          std::back_inserter(contour_part_)
      );
    }
    else
    {
      throw std::runtime_error("can't merge contours");
    }
  }



private:
  std::list<Spherical_Coords<RealType>> contour_part_;
  bool is_closed_ = false;
};

template<typename RealType>
class Contour_Database
{
public:
  Contour_Database(const std::string& input_file_path,
                   const std::string& output_file_path,
                   const std::vector<std::vector<RealType>>& rot)
    : input_file_path_(input_file_path)
    , output_file_path_(output_file_path)
    , rot_(rot)
  {
    ReadFromFile();
  }

  void ReadFromFile()
  {
    std::ifstream input(input_file_path_);
    std::string line;

    ContourPart<RealType> CurrentPart;

    getline(input, line);

    // Read all spherical Points and fill lists closed_contours and not_closed_contours_
    while (getline(input, line))
    {
      if (line[0] == '>')
      {
        if (CurrentPart.Size() > 0)
        {
          if (CurrentPart.is_closed())
          {
            islands_.push_back(CurrentPart);
          }
          else
          {
            not_closed_contours_.push_back(CurrentPart);
          }
          CurrentPart.Clear();
        }
      }
      else
      {
        std::stringstream ss(line);
        RealType first, second;
        ss >> first >> second;
        CurrentPart.AddPoint({first, second});
      }
    }
  }

  void MakeClosedDomain()
  {
    // information after reading
    std::cout <<"After reading" << std::endl;
    std::cout <<"Number of islands: " << islands_.size() << std::endl;
    std::cout <<"Number of not closed contour parts: " << not_closed_contours_.size() << std::endl;
    std::cout << std::endl;

    // making islands
    ContourPart<RealType> CurrentPart;

    while (not_closed_contours_.size() > 0)
    {
      CurrentPart = not_closed_contours_.front();
      not_closed_contours_.pop_front();
      while (MergeNextContour(CurrentPart, not_closed_contours_))
      {}

      if (CurrentPart.First() == CurrentPart.Last())
      {
        islands_.push_back(CurrentPart);
      }
      else
      {
        external_boundary_parts_.push_back(CurrentPart);
      }
    }


    //delete contours with only 2 points from external boundary parts
    auto iter =
        std::find_if
        (
            external_boundary_parts_.begin(),
            external_boundary_parts_.end(),
            [](const ContourPart<RealType>& cont)
            {
              return (cont.Size() <= 2);
            }
        );
    while (iter != external_boundary_parts_.end())
    {
      external_boundary_parts_.erase(iter);
      iter =
        std::find_if
        (
            external_boundary_parts_.begin(),
            external_boundary_parts_.end(),
            [](const ContourPart<RealType>& cont)
            {
              return (cont.Size() <= 3);
            }
        );
    }

    //delete contours with only 2 points from islands
    iter =
         std::find_if
         (
             islands_.begin(),
             islands_.end(),
             [](const ContourPart<RealType>& cont)
             {
               return (cont.Size() <= 3);
             }
         );
    while (iter != islands_.end())
    {
      islands_.erase(iter);
      iter =
        std::find_if
        (
          islands_.begin(),
          islands_.end(),
          [](const ContourPart<RealType>& cont)
          {
            return (cont.Size() <= 2);
          }
        );
     }


    std::cout <<"After making islands" << std::endl;
    std::cout <<"Number of islands: " << islands_.size() << std::endl;
    std::cout <<"Number of external boundary parts: " << external_boundary_parts_.size() << std::endl;
    std::cout << std::endl;

    // making common external border

    // 1) remove cropped islands and lakes

    iter =
        std::find_if
        (
            external_boundary_parts_.begin(),
            external_boundary_parts_.end(),
            [](const ContourPart<RealType>& cont)
            {
              return((((cont.First().Get_x() >= 1.0) or (cont.First().Get_x() <= 135.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)) or      // inside Russian coast
                     (((cont.First().Get_x() >= 140.0) or (cont.First().Get_x() <= -135.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)) or   // inside Pacific ocean
                     (((cont.First().Get_x() <= -65.0) or (cont.First().Get_x() >= -120.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)) or  // inside US/Canada coast
                     (((cont.First().Get_x() >= -60.0) or (cont.First().Get_x() <= -5.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)));      // inside Atlantic ocean
            }
        );

    while (iter != external_boundary_parts_.end())
    {
      external_boundary_parts_.erase(iter);

      iter =
        std::find_if
        (
            external_boundary_parts_.begin(),
            external_boundary_parts_.end(),
            [](const ContourPart<RealType>& cont)
            {
              return((((cont.First().Get_x() >= 1.0) or (cont.First().Get_x() <= 135.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)) or     // inside Russian coast
                     (((cont.First().Get_x() >= 140.0) or (cont.First().Get_x() <= -135.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)) or  // inside Pacific ocean
                     (((cont.First().Get_x() <= -65.0) or (cont.First().Get_x() >= -120.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)) or // inside US/Canada coast
                     (((cont.First().Get_x() >= -60.0) or (cont.First().Get_x() <= -5.0)) and (cont.First().Get_y() == EXTREME_LATITUDE)));     // inside Atlantic ocean
            }
        );
    }


    // 2) making rest of islands

    while (external_boundary_parts_.size() > 7)
    {
      CurrentPart = external_boundary_parts_.front();
      external_boundary_parts_.pop_front();
      while (MergeNextContour(CurrentPart, external_boundary_parts_))
      {}
      if (CurrentPart.First() == CurrentPart.Last())
      {
        islands_.push_back(CurrentPart);
        std::cout << "ISLAND" << std::endl;
      }
      else
      {
        external_boundary_parts_.push_back(CurrentPart);
      }
    }

  }

  void WriteToFile(bool islands)
  {
    if (not_closed_contours_.size() > 0)
    {
      throw std::runtime_error("can't write not closed contours");
    }

    std::ofstream output(output_file_path_);

    if (islands)
    {
      int counter = 0;
      for (auto& island: islands_)
      {
        output << "# island " << counter++ << std::endl;
        for (auto& point: island)
        {
          auto sph_coord = Rotate_Spherical(point, rot_);
          output << sph_coord << std::endl;
        }
      }
    }
    else
    {
      int counter = 0;
      for (auto& p: external_boundary_parts_)
      {
        output << "# external boundary part # " << counter++ << std::endl;
        for (auto& point: p)
        {
          auto sph_coord = Rotate_Spherical(point, rot_);
          output << sph_coord << std::endl;
        }
      }
    }
  }

private:
  const std::string                  input_file_path_;
  const std::string                  output_file_path_;
  std::vector<std::vector<RealType>> rot_;

  std::list<ContourPart<RealType>>   not_closed_contours_;
  std::list<ContourPart<RealType>>   islands_;
  std::list<ContourPart<RealType>>   external_boundary_parts_;
  ContourPart<RealType>              external_boundary_;
};


template<typename RealType>
bool MergeNextContour(ContourPart<RealType>& CurrentContour, std::list<ContourPart<RealType>>& CountourParts)
{
  auto next_contour_it =
      std::find_if
      (
          CountourParts.begin(),
          CountourParts.end(),
          [CurrentContour](const ContourPart<RealType>& cont)
          {
            return (cont.First() == CurrentContour.Last());
          }
      );
  if (next_contour_it == CountourParts.end())
  {
    return false;
  }
  else
  {
    CurrentContour.Merge(*next_contour_it);
    CountourParts.erase(next_contour_it);
    return true;
  }
  return true;
}
