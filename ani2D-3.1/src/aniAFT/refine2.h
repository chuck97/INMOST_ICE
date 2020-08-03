#ifndef H_REFINE2_MESH2D
#define H_REFINE2_MESH2D

#include <math.h>
#include <stdio.h>

#include "struct2.h"

double angle( int v1, int v2, int v );
void test_quality( void );
void smoothing( void );
void regularity( void );

#endif

