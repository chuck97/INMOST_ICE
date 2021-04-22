#include "inmost_ice.h"


double trsize(const 2dPoint& v1, 
              const 2dPoint& v2,
              const 2dPoint& v3)
{
    double a = sqrt((v1.x - v2.x)*(v1.x - v2.x) +   
                    (v1.y - v2.y)*(v1.y - v2.y));
    
    double b = sqrt((v1.x - v3.x)*(v1.x - v3.x) +   
                    (v1.y - v3.y)*(v1.y - v3.y));

    double c = sqrt((v3.x - v2.x)*(v3.x - v2.x) +   
                    (v3.y - v2.y)*(v3.y - v2.y));

    double p = (a+b+c)/2.0;
    return sqrt(p*(p-a)*(p-b)*(p-c)); 
}