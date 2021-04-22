#include "MatVec.h"

void vec_plus_vec(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax)
{
    for (unsigned int i = idmin; i < idmax; ++i)
    {
        res[i] = v1[i] + v2[i];
    }
}

void vec_minus_vec(const INMOST::Sparse::Vector& v1, const INMOST::Sparse::Vector& v2, INMOST::Sparse::Vector& res, unsigned int idmin, unsigned int idmax)
{
    for (unsigned int i = idmin; i < idmax; ++i)
    {
        res[i] = v1[i] - v2[i];
    }
}

void vec_mult_num(INMOST::Sparse::Vector& b, double scale, unsigned int idmin, unsigned int idmax)
{
    for (unsigned int i = idmin; i < idmax; ++i)
    {
        b[i]*=scale;
    }
}

void mat_mult_vec(INMOST::Sparse::Matrix& M, INMOST::Sparse::Vector& b, INMOST::Sparse::Vector& res)
{
    INMOST::Solver::OrderInfo info;
    info.PrepareMatrix(M, 0);
    info.PrepareVector(b);
    info.Update(b);
    BARRIER
    M.MatVec(1.0, b, 0.0, res);
    BARRIER
    info.RestoreVector(b);
    info.RestoreMatrix(M);
}

std::vector<double> operator*(double num, const std::vector<double>& v)
{
    std::vector<double> res = v;

    for (size_t i = 0; i < v.size(); ++i)
    {
        res[i] *= num; 
    }
    return res;
}

std::vector<std::vector<double>> operator*(double num, const std::vector<std::vector<double>>& m)
{
    std::vector<std::vector<double>> res = m;
    for (size_t i = 0; i < m[0].size(); ++i)
    {
        res[i] = num*m[i]; 
    }
    return res;
}

std::vector<double> operator*(const std::vector<std::vector<double>>& m, const std::vector<double>& v)
{
    if (m.size() != v.size())
    {
        INMOST_ICE_ERR("cant mult matr by vec - wrong sizes");
    }
    else
    {
        std::vector<double> res(v.size());
        for(size_t i = 0; i < v.size(); ++i)
        {
            res[i] = 0.0;
            for (size_t j = 0; j < v.size(); ++j)
            {
                res[i]+= v[j]*m[i][j];
            }
        }
        return res;
    }
}

std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& m1,
                                           const std::vector<std::vector<double>>& m2)
{
    std::vector<std::vector<double>> res = m1;
    if ((m1.size() != m2.size()) or (m1[0].size() != m2[0].size()))
    {
        throw(std::invalid_argument("wrong matricies sizes for summation"));
    }

    for (size_t i = 0; i < m2.size(); ++i)
    {
        for (size_t j = 0; j < m2[0].size(); ++j)
        {
            res[i][j]+=m2[i][j];
        }
    }
    return res;
}


std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& m1,
                                           const std::vector<std::vector<double>>& m2)
{
    std::vector<std::vector<double>> res = m1;
    if ((m1.size() != m2.size()) or (m1[0].size() != m2[0].size()))
    {
        throw(std::invalid_argument("wrong matricies sizes for difference"));
    }

    for (size_t i = 0; i < m2.size(); ++i)
    {
        for (size_t j = 0; j < m2[0].size(); ++j)
        {
            res[i][j]-=m2[i][j];
        }
    }
    return res;
}

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2)
{
    std::vector<double> res = v1;
    if (v1.size() != v2.size())
    {
        throw(std::invalid_argument("wrong vectors sizes for summation"));
    }

    for (size_t i = 0; i < v1.size(); ++i)
    {
        res[i]+=v2[i];
    }
    return res;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2)
{
    std::vector<double> res = v1;
    if (v1.size() != v2.size())
    {
        throw(std::invalid_argument("wrong vectors sizes for difference"));
    }

    for (size_t i = 0; i < v1.size(); ++i)
    {
        res[i]-=v2[i];
    }
    return res;
}

double Determinant3x3(const std::vector<std::vector<double>>& m)
{
    double a11 = m[0][0];
    double a12 = m[0][1];
    double a13 = m[0][2];
    double a21 = m[1][0];
    double a22 = m[1][1];
    double a23 = m[1][2];
    double a31 = m[2][0];
    double a32 = m[2][1];
    double a33 = m[2][2];

    return {a11*a22*a33 + a12*a23*a31 + a21*a32*a13 -
            a13*a22*a31 - a12*a21*a33 - a23*a32*a11};
}

// matrix should be stored as vector of columns!!!
std::vector<double> Solve3x3(const std::vector<std::vector<double>>& m, 
                             const std::vector<double>& rhs)
{
    double d = Determinant3x3(m);
    if (fabs(d) < 1e-15)
    {
        INMOST_ICE_ERR("almost singular 3x3 system, can't solve");
    }
    std::vector<std::vector<double>> m1 = m;
    std::vector<std::vector<double>> m2 = m;
    std::vector<std::vector<double>> m3 = m;

    m1[0] = rhs;
    m2[1] = rhs;
    m3[2] = rhs;

    double d1 = Determinant3x3(m1);
    double d2 = Determinant3x3(m2);
    double d3 = Determinant3x3(m3);

    return {d1/d, d2/d, d3/d};
}


double VectorsScalarProduct(std::vector<double> v1,
                            std::vector<double> v2)
{
    if (v1.size() != v2.size())
    {
        INMOST_ICE_ERR("scalar product of vector with different lengths");
    }
    else
    {
        double res = 0.0;
        for (size_t i = 0; i < v1.size(); ++i)
        {
            res += v1[i]*v2[i];
        }
        return res;
    }
}