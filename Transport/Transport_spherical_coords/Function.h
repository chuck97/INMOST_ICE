#pragma once

#include "config.h"

double Determinant3x3(const std::vector<std::vector<double>>& m);

// matrix should be stored as vector of columns!!!
std::vector<double> Solve3x3(const std::vector<std::vector<double>>& m, 
                             const std::vector<double>& rhs);

enum class FuncType 
{
     constant,
     linear,
     quadratic
};

double VectorsScalarProduct(std::vector<double> v1,
                            std::vector<double> v2);

class ScalarFunction
{
public:
    ScalarFunction(FuncType ft_,
                   std::vector<double> params_,
                   std::vector<std::vector<double>> trian_);
    double Evaluate(double x, double y) const;
    FuncType GetType() const;
    std::vector<std::vector<double>> GetTrian() const;
    std::vector<double> GetParams() const;
private:
    FuncType ft;
    std::vector<double> params;
    std::vector<std::vector<double>> trian;
};

class VectorFunction
{
public:
    VectorFunction(FuncType ft_,
                   std::vector<double> params1_,
                   std::vector<double> params2_,
                   std::vector<std::vector<double>> trian_);
    FuncType GetType() const;
    std::vector<std::vector<double>> GetTrian() const;
    std::vector<std::vector<double>>  GetParams() const;
private:
    FuncType ft;
    std::vector<double> params1;
    std::vector<double> params2;
    std::vector<std::vector<double>> trian;
};

ScalarFunction operator*(const VectorFunction& v1, const VectorFunction& v2);

ScalarFunction mult_linear_by_cos(const ScalarFunction& s);

VectorFunction Gradient(ScalarFunction s);

ScalarFunction Divirgence(VectorFunction v);

ScalarFunction d_dx(ScalarFunction s);

ScalarFunction d_dy(ScalarFunction s);

ScalarFunction operator*(const ScalarFunction& s1, const ScalarFunction& s2);

VectorFunction operator*(const ScalarFunction& s, const VectorFunction& v);

ScalarFunction operator+(const ScalarFunction& s1, const ScalarFunction& s2);