#pragma once
  
namespace mapx
{
struct mapx_class;
mapx_class* init_mapx(char* filename);
int forward_xy_mapx(mapx_class* this_, double lat, double lon, double* x, double* y);
void close_mapx(mapx_class* this_);
}
