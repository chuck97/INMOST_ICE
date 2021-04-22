#pragma once

// line trim
void ltrim(std::string& s);
void rtrim(std::string& s);
void trim(std::string& s);

// 2d point
struct _2dPoint
{
    double x = 0.0;
    double y = 0.0;
};

// triangle square from verticies
double trsize(const _2dPoint& v1, 
              const _2dPoint& v2,
              const _2dPoint& v3);