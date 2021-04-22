#pragma once
#include "config.h"


void vec_plus_vec(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax);
void vec_minus_vec(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax);
void mat_mult_vec(INMOST::Sparse::Matrix& M, INMOST::Sparse::Vector& b, INMOST::Sparse::Vector& res);
void vec_mult_num(INMOST::Sparse::Vector& b, double scale, unsigned int idmin, unsigned int idmax);

std::vector<double> operator*(double num, const std::vector<double>& v);

std::vector<std::vector<double>> operator*(double num, const std::vector<std::vector<double>>& m);

std::vector<double> operator*(const std::vector<std::vector<double>>& m, const std::vector<double>& v);

std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& m1,
                                           const std::vector<std::vector<double>>& m2);

std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& m1,
                                           const std::vector<std::vector<double>>& m2);

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);


double Determinant3x3(const std::vector<std::vector<double>>& m);

// matrix should be stored as vector of columns!!!
std::vector<double> Solve3x3(const std::vector<std::vector<double>>& m, 
                             const std::vector<double>& rhs);

double VectorsScalarProduct(std::vector<double> v1,
                            std::vector<double> v2);