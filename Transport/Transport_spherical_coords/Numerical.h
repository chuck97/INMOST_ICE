#pragma once
#include "config.h"
#include "Function.h"

std::vector<double> FromBariToOrdinary(const std::vector<std::vector<double>>& trnodes,
                                       const std::vector<double>& bcoords);

double triansize(double x0, double y0, double x1, double y1, double x2, double y2);

double integral_over_triangle(const std::vector<std::vector<double>>& trnodes,
                              ScalarFunction f);