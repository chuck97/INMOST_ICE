//  This header file contains functions and classes for working with spherical and Cartesian coordinates/
//
//  1) ALPHA_DEF, BETA_DEF, GAMMA_DEF -- Euler angles in degrees
//  2) eps                            -- offset value
//  3) Euler_rotation_info            -- class that stores information about rotation and performs it
//  4) Spherical_Coords               -- class that stores spherical coords
//  5) Cartesian_Coords               -- class that stores Cartesian coords
//  6) Overloaded operators <<        -- functions that prints spherical and Cartesian coordinates to output stream
//  7) Overloaded operator ==         -- checks equality of spherical coordinates
//  8) From_spherical_to_Cartesian    -- function that performs transformation from spherical to Cartesian coordinates (with unit radius)
//  9) From_Cartesian_to_spherical    -- function that performs transformation from Cartesian to spherical coordinates (with unit radius)
// 10) Rotate_Cartesian               -- function that performs rotation of point with Cartesian coordinates
// 11) Rotate_Spherical               -- function that performs rotation of point with spherical coordinates

#pragma once
#include <vector>
#include <cmath>
#include <cmath>
#include <ostream>

namespace CORDROT
{

const double ALPHA_DEF               = -30.0;  // default alpha for Euler rotation in degrees
const double BETA_DEF                = -90.0;  // default beta for Euler rotation in degrees
const double GAMMA_DEF               = 0.0;    // default gamma for Euler rotation in degrees
const double eps                     = 1e-10;  // offset to avoid zero values
const int    OUTPUT_DOUBLE_PRECISION = 13;     // total number of digits in output double


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  #
// #    !!! Purpose of class - store the information about rotation !!!                                  #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                  #
// #                                                                                                     #
// #    Constructor input : Euler angles: alpha, beta and gamma                                          #
// #    Contains:                                                                                        #
// #    ALPHA   - real rotation angle about z-axis in degrees can be obtained by call .Get_ALPHA()       #
// #    BETA    - real rotation angle about new y-axis in degrees can be obtained by call .Get_BETA()    #
// #    GAMMA   - real rotation angle about new z-axis in degrees can be obtained by call .Get_GAMMA()   #
// #    FORWARD - 3x3 real matrix that perform rotation from model to geographical Cartesian coordinates #
// #          (appears after a call constructor) can be obtained by call .Get_FORWARD()                  #
// #    REVERSE - 3x3 real matrix that perform rotation from geographical to model Cartesian coordinates #
// #          (appears after a call constructor) can be obtained by call .Get_REVERSE()                  #
// #                                                                                                     #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Euler_rotation_info
{
public:
  Euler_rotation_info(double alpha, double beta, double gamma)
  : ALPHA(alpha)
  , BETA (beta)
  , GAMMA(gamma)
  {
    GenerateRotationMatricies();
  }

  const std::vector<std::vector<double> >& Get_FORWARD() const
   {
     return FORWARD;
   }

  const std::vector<std::vector<double> >& Get_REVERSE() const
   {
     return REVERSE;
   }

   double Get_ALPHA() const
   {
     return ALPHA;
   }

   double Get_BETA() const
   {
     return BETA;
   }

   double Get_GAMMA() const
   {
     return BETA;
   }
private:
 const double ALPHA;                        // rotation about z-axis
 const double BETA;                         // rotation about new y-axis
 const double GAMMA;                        // rotation about new z-axis
 std::vector<std::vector<double> > FORWARD;  // rotation matricies to geographic coordinates
 std::vector<std::vector<double> > REVERSE;  // rotation matricies to model coordinates

 void GenerateRotationMatricies()
 {
   double ALPHA_RAD = ALPHA*M_PI/180.0;
   double BETA_RAD = BETA*M_PI/180.0;
   double GAMMA_RAD = GAMMA*M_PI/180.0;

   FORWARD.push_back(std::vector<double> {
     cos(GAMMA_RAD)*cos(BETA_RAD)*cos(ALPHA_RAD) - sin(GAMMA_RAD)*sin(ALPHA_RAD),
     -sin(GAMMA_RAD)*cos(BETA_RAD)*cos(ALPHA_RAD) - cos(GAMMA_RAD)*sin(ALPHA_RAD),
     sin(BETA_RAD)*cos(ALPHA_RAD)
   });

   FORWARD.push_back(std::vector<double>{
     cos(GAMMA_RAD)*cos(BETA_RAD)*sin(ALPHA_RAD) + sin(GAMMA_RAD)*cos(ALPHA_RAD),
     -sin(GAMMA_RAD)*cos(BETA_RAD)*sin(ALPHA_RAD) + cos(GAMMA_RAD)*cos(ALPHA_RAD),
     sin(BETA_RAD)*sin(ALPHA_RAD)
   });

   FORWARD.push_back(std::vector<double>{
     -cos(GAMMA_RAD)*sin(BETA_RAD),
     sin(GAMMA_RAD)*sin(BETA_RAD),
     cos(BETA_RAD)
   });

   REVERSE.push_back(std::vector<double>{
     cos(ALPHA_RAD)*cos(BETA_RAD)*cos(GAMMA_RAD) - sin(ALPHA_RAD)*sin(GAMMA_RAD),
     sin(ALPHA_RAD)*cos(BETA_RAD)*cos(GAMMA_RAD) + cos(ALPHA_RAD)*sin(GAMMA_RAD),
     -sin(BETA_RAD)*cos(GAMMA_RAD)
   });

   REVERSE.push_back(std::vector<double>{
     -cos(ALPHA_RAD)*cos(BETA_RAD)*sin(GAMMA_RAD) - sin(ALPHA_RAD)*cos(GAMMA_RAD),
     -sin(ALPHA_RAD)*cos(BETA_RAD)*sin(GAMMA_RAD) + cos(ALPHA_RAD)*cos(GAMMA_RAD),
     sin(BETA_RAD)*sin(GAMMA_RAD)
   });

   REVERSE.push_back(std::vector<double>{
     cos(ALPHA_RAD)*sin(BETA_RAD),
     sin(ALPHA_RAD)*sin(BETA_RAD),
     cos(BETA_RAD)
   });

   FORWARD.shrink_to_fit();
   REVERSE.shrink_to_fit();
 }
};


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                            #
// #    !!! Purpose of class - store spherical coordinates of point !!!                            #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                            #
// #                                                                                               #
// #     Contains:                                                                                 #
// #     x_      - real first coordinate  (longitude) in degrees can be obtained by call .Get_x()  #
// #     y_      - real second coordinate (latitude) in degrees can be obtained by call .Get_y()   #
// #     +=      - increment real number to both components                                        #
// #     -=      - decrement real number to both components                                        #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Spherical_Coords
{
public:
  Spherical_Coords()
  {}

  Spherical_Coords(double x, double y)
  : x_(x)
  , y_(y)
  {}

  double Get_x() const
  {
    return x_;
  }

  double Get_y() const
  {
    return y_;
  }

  void operator += (double increment)
  {
    x_ += increment;
    y_ += increment;
  }

  void operator -= (double decrement)
  {
    x_ -= decrement;
    y_ -= decrement;
  }


private:
  double x_;
  double y_;
};



// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #    !!! Purpose of class - store Cartesian coordinates of point !!!      #
// #    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                         #
// #     Contains:                                                           #
// #     x_      - real first coordinate can be obtained by call .Get_x()    #
// #     y_      - real second coordinate can be obtained by call .Get_y()   #
// #     z_      - real third coordinate can be obtained by call .Get_z()    #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Cartesian_Coords
{
public:
  Cartesian_Coords()
  {}

  Cartesian_Coords(double x, double y, double z)
  : x_(x)
  , y_(y)
  , z_(z)
  {}

  double Get_x() const
  {
    return x_;
  }

  double Get_y() const
  {
    return y_;
  }

  double Get_z() const
  {
    return z_;
  }

private:
  double x_;
  double y_;
  double z_;
};

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - print spherical coords to output stream (cout or file) !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : point with spherical coordinates                                                   #
// #    output: reference to output stream                                                         #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
std::ostream& operator<<(std::ostream& stream, const Spherical_Coords& cords)
{
  stream.precision(OUTPUT_DOUBLE_PRECISION);
  stream << cords.Get_x() << ' ' << cords.Get_y(); //<< cords.IsWater() << ' ';
  return stream;
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - print Cartesian coords to output stream (cout or file) !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : point with Cartesian coordinates                                                   #
// #    output: reference to output stream                                                         #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
std::ostream& operator<<(std::ostream& stream, const Cartesian_Coords& cords)
{
  stream.precision(OUTPUT_DOUBLE_PRECISION);
  stream << cords.Get_x() << ' ' << cords.Get_y() << ' ' << cords.Get_z();
  return stream;
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - identify equality of two spherical points              !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input1 : lhs - first spherical point                                                       #
// #    input2 : rhs - second spherical point                                                      #
// #    output: boolean (true for =, false for !=)                                                 #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
bool operator==(const Spherical_Coords& lhs, const Spherical_Coords& rhs)
{
  return ((lhs.Get_x() == rhs.Get_x()) and (lhs.Get_y() == rhs.Get_y()));
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - transform spherical point to Cartesian one             !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : Sph_coords - point with spherical coords                                           #
// #    output: point with Cartesian coords                                                        #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Cartesian_Coords From_spherical_to_Cartesian(const Spherical_Coords& Sph_coords)
{
  double x_rad = Sph_coords.Get_x()*M_PI/180.0;
  double y_rad = Sph_coords.Get_y()*M_PI/180.0;

  return {cos(y_rad)*cos(x_rad),
          cos(y_rad)*sin(x_rad),
          sin(y_rad)};
}


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!S!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - transform Cartesian point to spherical one             !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input : Crt_coords - point with Cartesian coords                                           #
// #    output: point with spherical coords                                                        #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Spherical_Coords From_Cartesian_to_spherical(const Cartesian_Coords& Crt_coords)
{
  double xx = Crt_coords.Get_x();
  double yy = Crt_coords.Get_y();
  double zz = Crt_coords.Get_z();

  double phi, theta, costn;

  theta = asin(zz)*180/M_PI;
  costn = sqrt(1.0 - zz*zz);

  if ((xx > 0.0) and (yy > 0.0))
    {
      if (xx < yy)
      {
        phi = acos(xx/costn)*180.0/M_PI;
      }
      else
      {
        phi = asin(yy/costn)*180.0/M_PI;
      }
    }
    else if ((xx < 0.0) and (yy > 0.0))
    {
      if (abs(xx) < yy)
      {
        phi = 180.0 - acos(abs(xx)/costn)*180.0/M_PI;
      }
      else
      {
        phi = 180.0 - asin(yy/costn)*180.0/M_PI;
      }
    }
    else if ((xx < 0.0) and (yy < 0.0))
    {
      if (abs(xx) < abs(yy))
      {
        phi = -180.0 + acos(abs(xx)/costn)*180/M_PI;
      }
      else
      {
        phi =-180.0 + asin(abs(yy)/costn)*180/M_PI;
      }
    }
    else if ((xx > 0.0) and (yy < 0.0))
    {
      if (xx < abs(yy))
      {
        phi = -acos(abs(xx)/costn)*180.0/M_PI;
      }
      else
      {
        phi = -asin(abs(yy)/costn)*180/M_PI;
      }
    }

    //if (phi < 0.0)
    //{
    //  phi += 360.0;
    //}

    return {phi, theta};
}


// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - calculate new Cartesian coordinates of point after     !!!      #
// #     !!!                       rotation with given rotation matrix                    !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input1: cartesian_old   - Cartesian coordinates of point before rotation                   #
// #    input2: rot             - 3x3 rotation matrix                                              #
// #    output: Cartesian coordinate of point after rotation                                       #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Cartesian_Coords Rotate_Cartesian(const Cartesian_Coords&  cartesian_old, const std::vector<std::vector<double> >& rot)
{
  return {rot[0][0]*cartesian_old.Get_x() + rot[0][1]*cartesian_old.Get_y() + rot[0][2]*cartesian_old.Get_z(),
          rot[1][0]*cartesian_old.Get_x() + rot[1][1]*cartesian_old.Get_y() + rot[1][2]*cartesian_old.Get_z(),
          rot[2][0]*cartesian_old.Get_x() + rot[2][1]*cartesian_old.Get_y() + rot[2][2]*cartesian_old.Get_z()
  };
}

// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #     !!! Purpose of function - calculate new spherical coordinates of point after     !!!      #
// #     !!!                       rotation with given rotation matrix                    !!!      #
// #     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      #
// #                                                                                               #
// #    input1: spherical_old   - spherical coordinates of point before rotation                   #
// #    input2: rot             - 3x3 rotation matrix                                              #
// #    output: spherical coordinate of point after rotation                                       #
// # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
Spherical_Coords Rotate_Spherical(const Spherical_Coords& spherical_old, const std::vector<std::vector<double> >& rot)
{

  Spherical_Coords  spherical_new;   //output
  Spherical_Coords  spherical_tmp;   //temporary
  Cartesian_Coords  cartesian_old;   //old Cartesian coordinate
  Cartesian_Coords  cartesian_new;   //new Cartesian coordinate

  spherical_tmp = spherical_old;

  //avoid trouble of an exactly zero angle by adding offset
  spherical_tmp += eps;

  //spherical to Cartesian
  cartesian_old = From_spherical_to_Cartesian(spherical_tmp);

  //new Cartesian coordinates after rotation
  cartesian_new = Rotate_Cartesian(cartesian_old, rot);

  //Cartesian to spherical
  spherical_new = From_Cartesian_to_spherical(cartesian_new);

   //avoid trouble of an exactly zero angle by subtracting offset
  spherical_new -= eps;

  return spherical_new;
}

}