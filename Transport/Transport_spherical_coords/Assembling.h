#pragma once

#include "Function.h"
#include "Numerical.h"

std::vector<double> operator*(double num, const std::vector<double>& v);

std::vector<std::vector<double>> operator*(double num, const std::vector<std::vector<double>>& m);

std::vector<double> operator*(const std::vector<std::vector<double>>& m, const std::vector<double>& v);

std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& m1,
                                           const std::vector<std::vector<double>>& m2);

std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& m1,
                                           const std::vector<std::vector<double>>& m2);

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);


std::vector<std::vector<double>> LocalStiffnessMatrixAssembling(const std::vector<std::vector<double>>& trnodes);

double L_ij_entire(const std::vector<std::vector<double>>& trnodes,
                   const ScalarFunction& phi_i,
                   const ScalarFunction& phi_j,
                   const ScalarFunction& u0,
                   const ScalarFunction& u1);


std::vector<double> LocalTG2RhsAssembling(const std::vector<std::vector<double>>& trnodes,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    double time_step);

std::vector<double> LocalCG2RhsAssembling(const std::vector<std::vector<double>>& trnodes,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    double time_step);

std::vector<double> LocalTTG2RhsAssembling(const std::vector<std::vector<double>>& trnodes,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    const std::vector<double>& localmass_half,
                                                    double time_step,
                                                    int step);

std::vector<double> LocalTTG3RhsAssembling(const std::vector<std::vector<double>>& trnodes,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    const std::vector<double>& localmass_half,
                                                    double time_step,
                                                    int step);

std::vector<double> LocalTTG4RhsAssembling(const std::vector<std::vector<double>>& trnodes,
                                                    const std::vector<std::vector<double>>& uvalues,
                                                    const std::vector<double>& localmass,
                                                    const std::vector<double>& localmass_half,
                                                    double time_step,
                                                    int step);