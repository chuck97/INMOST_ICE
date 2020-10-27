#include<iostream>
#include "coords_rotation.h"

int main()
{
   Euler_rotation_info<double> rotation(ALPHA_DEF, BETA_DEF, GAMMA_DEF); 
   

   double lat, lon, model_lat, model_lon;   

   // GEO TO MODEL
   
   std::cout << "geo (lon, lat) in degrees:";   
   std::cin >>  lon >> lat;
   
   Spherical_Coords<double> geo_rev(lon, lat);
   Spherical_Coords<double> model_rev = Rotate_Spherical<double>(geo_rev, rotation.Get_REVERSE());
   
   model_lon = model_rev.Get_x();
   model_lat = model_rev.Get_y();
   
   std::cout << "calculated model (lon, lat) in degrees:" << model_lon << " " << model_lat << std::endl;


   // MODEL TO GEO

   std::cout << "model (lon, lat) in degrees:";   
   std::cin >>  model_lon >> model_lat;
   
   Spherical_Coords<double> model_forw(model_lon, model_lat);
   Spherical_Coords<double> geo_forw = Rotate_Spherical<double>(model_forw, rotation.Get_FORWARD());
   
   lon = geo_forw.Get_x();
   lat = geo_forw.Get_y();
   
   std::cout << "calculated geo (lon, lat):" << lon << " " << lat << std::endl;

   return 0;
   
}
