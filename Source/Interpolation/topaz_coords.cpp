#include "inmost_ice.h"

std::vector<double> from_geo_2_topaz(double lon, double lat)
{
  projPJ pj_longlat, pj_topaz;
  int p;

  if (!(pj_longlat = pj_init_plus("+proj=longlat +ellps=clrk66")) )
  {
    INMOST_ICE_ERR("error pj_longlat");
  }

  if (!(pj_topaz = pj_init_plus("+units=m +proj=stere +a=6378273.0 +b=6378273.0 +lon_0=-45.0 +lat_0=90.0 +lat_ts=90.0")) )
  {
    INMOST_ICE_ERR("error pj_topaz");
  }
   
  double lon_rad = lon*DEG_TO_RAD;
  double lat_rad = lat*DEG_TO_RAD; 

  p = pj_transform(pj_longlat, pj_topaz, 1, 1, &lon_rad, &lat_rad, NULL);
  
  pj_free(pj_longlat);
  pj_free(pj_topaz);
  
  return{lon_rad, lat_rad};
}

std::vector<double> from_topaz_2_geo(double x, double y)
{
  projPJ pj_longlat, pj_topaz;
  int p;

  if (!(pj_longlat = pj_init_plus("+proj=longlat +ellps=clrk66")) )
  {
    INMOST_ICE_ERR("error pj_longlat");
  }

  if (!(pj_topaz = pj_init_plus("+units=m +proj=stere +a=6378273.0 +b=6378273.0 +lon_0=-45.0 +lat_0=90.0 +lat_ts=90.0")) )
  {
    INMOST_ICE_ERR("error pj_topaz");
  }
   
  p = pj_transform(pj_topaz, pj_longlat, 1, 1, &x, &y, NULL);

  pj_free(pj_longlat);
  pj_free(pj_topaz);

  return {x*RAD_TO_DEG, y*RAD_TO_DEG};  
}
