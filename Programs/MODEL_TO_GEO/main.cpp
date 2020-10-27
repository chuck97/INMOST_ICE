#include<iostream>
#include "coords_rotation.h"

int main()
{
   CORDROT::Euler_rotation_info rotation(CORDROT::ALPHA_DEF, CORDROT::BETA_DEF, CORDROT::GAMMA_DEF); 
   

   double lat, lon, model_lat, model_lon;   

   // GEO TO MODEL
   
   std::cout << "geo (lon, lat) in degrees:";   
   std::cin >>  lon >> lat;
   
   CORDROT::Spherical_Coords geo_rev(lon, lat);
   CORDROT::Spherical_Coords model_rev = CORDROT::Rotate_Spherical(geo_rev, rotation.Get_REVERSE());
   
   model_lon = model_rev.Get_x();
   model_lat = model_rev.Get_y();
   
   std::cout << "calculated model (lon, lat) in degrees:" << model_lon << " " << model_lat << std::endl;


   // MODEL TO GEO

   std::cout << "model (lon, lat) in degrees:";   
   std::cin >>  model_lon >> model_lat;
   
   CORDROT::Spherical_Coords model_forw(model_lon, model_lat);
   CORDROT::Spherical_Coords geo_forw = CORDROT::Rotate_Spherical(model_forw, rotation.Get_FORWARD());
   
   lon = geo_forw.Get_x();
   lat = geo_forw.Get_y();
   
   std::cout << "calculated geo (lon, lat):" << lon << " " << lat << std::endl;

   return 0;
   
}
