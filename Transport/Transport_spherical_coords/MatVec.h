#pragma once
#include "config.h"


void vec_plus_vec(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax);
void vec_minus_vec(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax);
void mat_mult_vec(INMOST::Sparse::Matrix& M, INMOST::Sparse::Vector& b, INMOST::Sparse::Vector& res);
void vec_mult_num(INMOST::Sparse::Vector& b, double scale, unsigned int idmin, unsigned int idmax);