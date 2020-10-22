#include<iostream>
#include <proj_api.h>

int main(int argc, char **argv) 
{
    projPJ pj_longlat, pj_topaz;
    double x, y;
    int p;

    if (!(pj_longlat = pj_init_plus("+proj=longlat +ellps=clrk66")) )
    {
        std::cout << "error pj_longlat" << std::endl;
        return 1;
    }
    if (!(pj_topaz = pj_init_plus("+units=m +proj=stere +a=6378273.0 +b=6378273.0 +lon_0=-45.0 +lat_0=90.0 +lat_ts=90.0")) )
    {
        std::cout << "error pj_topaz" << std::endl;
        return 1;
    }
   
    std::cout << "geo to topaz (long, lat):" << std::endl;
    std::cin >> x >> y;
    x *= DEG_TO_RAD;
    y *= DEG_TO_RAD;
    p = pj_transform(pj_longlat, pj_topaz, 1, 1, &x, &y, NULL);
    std::cout << x << " " << y << std::endl;

    std::cout << "topaz to geo (x, y):" << std::endl;    
    std::cin >> x >> y;
    p = pj_transform(pj_topaz, pj_longlat, 1, 1, &x, &y, NULL);
    std::cout << x*RAD_TO_DEG << " " << y*RAD_TO_DEG << std::endl;

    pj_free(pj_longlat);
    pj_free(pj_topaz);

    return 0;
}
