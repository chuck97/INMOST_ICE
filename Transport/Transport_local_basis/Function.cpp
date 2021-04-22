#include "Function.h"

ScalarFunction::ScalarFunction(FuncType ft_,
                               std::vector<double> params_, 
                               std::vector<std::vector<double>> trian_)
{
    if (ft_ == FuncType::constant)
    {
        ft = FuncType::constant;
        if (params_.size() != 1)
        {
            INMOST_ICE_ERR("const function must have 1 param");
        }
    }
    else if (ft_ == FuncType::linear)
    {
        ft = FuncType::linear;
        if (params_.size() != 3)
        {
            INMOST_ICE_ERR("linear function must have 3 param");
        }
    }
    else if (ft_ == FuncType::quadratic)
    {
        ft = FuncType::quadratic;
        if (params_.size() != 6)
        {
            INMOST_ICE_ERR("quadratic function must have 6 param");
        }
    }
    else
    {
        INMOST_ICE_ERR("unknown type of scalar function");
    }
    params = params_;
    trian = trian_;
}

double ScalarFunction::Evaluate(double x, double y) const
{
    if (ft == FuncType::constant)
    {
        return params[0];
    }
    else if (ft == FuncType::linear)
    {
        return x*params[0] + y*params[1] + params[2];
    }
    else if (ft == FuncType::quadratic)
    {
        return x*x*params[0] + x*y*params[1] + y*y*params[2] +
               x*params[3] + y*params[4] + params[5];
    }
    else 
    {
        INMOST_ICE_ERR("unknown type of function");
    }
}

FuncType ScalarFunction::GetType() const
{
    return ft;
}

std::vector<std::vector<double>> ScalarFunction::GetTrian() const
{
    return trian;
}

std::vector<double> ScalarFunction::GetParams() const
{
    return params;
}

VectorFunction::VectorFunction(FuncType ft_,
                               std::vector<double> params1_,
                               std::vector<double> params2_,
                               std::vector<std::vector<double>> trian_)
{
    if (ft_ == FuncType::constant)
    {
        ft = FuncType::constant;
        if (params1_.size() != 1)
        {
            INMOST_ICE_ERR("const function must have 1 param");
        }
        if (params2_.size() != 1)
        {
            INMOST_ICE_ERR("const function must have 1 param");
        }
    }
    else if (ft_ == FuncType::linear)
    {
        ft = FuncType::linear;
        if (params1_.size() != 3)
        {
            INMOST_ICE_ERR("linear function must have 3 param");
        }
        if (params2_.size() != 3)
        {
            INMOST_ICE_ERR("linear function must have 3 param");
        }
    }
    else if (ft_ == FuncType::quadratic)
    {
        ft = FuncType::quadratic;
        if (params1_.size() != 6)
        {
            INMOST_ICE_ERR("quadratic function must have 6 param");
        }
        if (params2_.size() != 6)
        {
            INMOST_ICE_ERR("quadratic function must have 6 param");
        }
    }
    else
    {
        INMOST_ICE_ERR("unknown type of scalar function");
    }
    params1 = params1_;
    params2 = params2_;
    trian = trian_;
}

FuncType VectorFunction::GetType() const  
{
    return ft;
}

std::vector<std::vector<double>> VectorFunction::GetTrian() const
{
    return trian;
}

std::vector<std::vector<double>>  VectorFunction::GetParams() const
{
    return {params1, params2};
}

ScalarFunction operator*(const VectorFunction& v1, const VectorFunction& v2)
{
    std::vector<std::vector<double>> tr = v1.GetTrian();
    if (v1.GetType() == FuncType::constant)
    {
        std::vector<double> lhs = {v1.GetParams()[0][0], v1.GetParams()[1][0]}; 
        if (v2.GetType() == FuncType::constant)
        {
            std::vector<double> rhs = {v2.GetParams()[0][0], v2.GetParams()[1][0]};
            double res = VectorsScalarProduct(lhs, rhs);
            return {FuncType::constant, {res}, tr};
        }
        else if (v2.GetType() == FuncType::linear)
        {
            std::vector<double> rhs1 = {v2.GetParams()[0][0], v2.GetParams()[1][0]};
            std::vector<double> rhs2 = {v2.GetParams()[0][1], v2.GetParams()[1][1]};
            std::vector<double> rhs3 = {v2.GetParams()[0][2], v2.GetParams()[1][2]};

            double res1 = VectorsScalarProduct(lhs, rhs1);
            double res2 = VectorsScalarProduct(lhs, rhs2);
            double res3 = VectorsScalarProduct(lhs, rhs3);

            return {FuncType::linear, {res1, res2, res3}, tr};
        }
        else if (v2.GetType() == FuncType::quadratic)
        {
            std::vector<double> rhs1 = {v2.GetParams()[0][0], v2.GetParams()[1][0]};
            std::vector<double> rhs2 = {v2.GetParams()[0][1], v2.GetParams()[1][1]};
            std::vector<double> rhs3 = {v2.GetParams()[0][2], v2.GetParams()[1][2]};
            std::vector<double> rhs4 = {v2.GetParams()[0][3], v2.GetParams()[1][3]};
            std::vector<double> rhs5 = {v2.GetParams()[0][4], v2.GetParams()[1][4]};
            std::vector<double> rhs6 = {v2.GetParams()[0][5], v2.GetParams()[1][5]};

            double res1 = VectorsScalarProduct(lhs, rhs1);
            double res2 = VectorsScalarProduct(lhs, rhs2);
            double res3 = VectorsScalarProduct(lhs, rhs3);
            double res4 = VectorsScalarProduct(lhs, rhs4);
            double res5 = VectorsScalarProduct(lhs, rhs5);
            double res6 = VectorsScalarProduct(lhs, rhs6);

            return {FuncType::quadratic, {res1, res2, res3, res4, res5, res6}, tr};
        }
        else
        {
            INMOST_ICE_ERR("max degree is 3");
        }
    }
    else if (v1.GetType() == FuncType::linear)
    {
        std::vector<double> lhs1 = {v1.GetParams()[0][0], v1.GetParams()[1][0]};
        std::vector<double> lhs2 = {v1.GetParams()[0][1], v1.GetParams()[1][1]};
        std::vector<double> lhs3 = {v1.GetParams()[0][2], v1.GetParams()[1][2]};
        if (v2.GetType() == FuncType::constant)
        {
            std::vector<double> rhs = {v2.GetParams()[0][0], v2.GetParams()[1][0]};
            double res1 = VectorsScalarProduct(lhs1, rhs);
            double res2 = VectorsScalarProduct(lhs2, rhs);
            double res3 = VectorsScalarProduct(lhs3, rhs);
            return {FuncType::linear, {res1, res2, res3}, tr};
        }
        else if (v2.GetType() == FuncType::linear)
        {
            std::vector<double> rhs1 = {v2.GetParams()[0][0], v2.GetParams()[1][0]};
            std::vector<double> rhs2 = {v2.GetParams()[0][1], v2.GetParams()[1][1]};
            std::vector<double> rhs3 = {v2.GetParams()[0][2], v2.GetParams()[1][2]};

            double res1 = VectorsScalarProduct(lhs1, rhs1);
            double res2 = VectorsScalarProduct(lhs1, rhs2) + VectorsScalarProduct(lhs2, rhs1);
            double res3 = VectorsScalarProduct(lhs2, rhs2);
            double res4 = VectorsScalarProduct(lhs1, rhs3) + VectorsScalarProduct(lhs3, rhs1);
            double res5 = VectorsScalarProduct(lhs2, rhs3) + VectorsScalarProduct(lhs3, rhs2);
            double res6 = VectorsScalarProduct(lhs3, rhs3);

            return {FuncType::quadratic, {res1, res2, res3, res4, res5, res6}, tr};
        }
        else
        {
            INMOST_ICE_ERR("max degree is 3");
        }
    }
    else if (v1.GetType() == FuncType::quadratic)
    {
        std::vector<double> lhs1 = {v1.GetParams()[0][0], v1.GetParams()[1][0]};
        std::vector<double> lhs2 = {v1.GetParams()[0][1], v1.GetParams()[1][1]};
        std::vector<double> lhs3 = {v1.GetParams()[0][2], v1.GetParams()[1][2]};
        std::vector<double> lhs4 = {v1.GetParams()[0][3], v1.GetParams()[1][3]};
        std::vector<double> lhs5 = {v1.GetParams()[0][4], v1.GetParams()[1][4]};
        std::vector<double> lhs6 = {v1.GetParams()[0][5], v1.GetParams()[1][5]};
        if (v2.GetType() == FuncType::constant)
        {
            std::vector<double> rhs = {v2.GetParams()[0][0], v2.GetParams()[1][0]};

            double res1 = VectorsScalarProduct(lhs1, rhs);
            double res2 = VectorsScalarProduct(lhs2, rhs);
            double res3 = VectorsScalarProduct(lhs3, rhs);
            double res4 = VectorsScalarProduct(lhs4, rhs);
            double res5 = VectorsScalarProduct(lhs5, rhs);
            double res6 = VectorsScalarProduct(lhs6, rhs);

            return {FuncType::quadratic, {res1, res2, res3, res4, res5, res6}, tr};
        }
        else
        {
            INMOST_ICE_ERR("max degree is 3");
        }
    }
    else
    {
        INMOST_ICE_ERR("max degree is 3");
    }
}

VectorFunction Gradient(ScalarFunction s)
{
    if (s.GetType() == FuncType::constant)
    {
        INMOST_ICE_ERR("grad of const is 0");
    }
    else if (s.GetType() == FuncType::linear)
    {
        std::vector<std::vector<double>> tr = s.GetTrian();
        return{FuncType::constant, {s.GetParams()[0]},
                                   {s.GetParams()[1]}, tr};
    }
    else if (s.GetType() == FuncType::quadratic)
    {
        INMOST_ICE_ERR("can't grad quadratic function");
    }
    else
    {
        INMOST_ICE_ERR("max degree for grad is 3");
    } 
}

ScalarFunction Divirgence(VectorFunction v)
{
    if (v.GetType() == FuncType::constant)
    {
        INMOST_ICE_ERR("divergence of const is 0");
    }
    else if (v.GetType() == FuncType::linear)
    {
        std::vector<std::vector<double>> tr = v.GetTrian();        
        return{FuncType::constant, {v.GetParams()[0][0]+
                                    v.GetParams()[1][1]}, tr};
    }
    else if (v.GetType() == FuncType::quadratic)
    {
        INMOST_ICE_ERR("can't div quadratic function");
    }
    else
    {
        INMOST_ICE_ERR("max degree for div is 2");
    }
}

ScalarFunction d_dx(ScalarFunction s)
{
    if (s.GetType() == FuncType::constant)
    {
        INMOST_ICE_ERR("d_dx of const is 0");
    }
    else if (s.GetType() == FuncType::linear)
    {
        std::vector<std::vector<double>> tr = s.GetTrian();

        return{FuncType::constant, {s.GetParams()[0]}, tr};
    }
    else if (s.GetType() == FuncType::quadratic)
    {
        INMOST_ICE_ERR("can't d/dx quadratic function");
    }
    else
    {
        INMOST_ICE_ERR("max degree for d_dx is 2");
    }
}

ScalarFunction d_dy(ScalarFunction s)
{
    if (s.GetType() == FuncType::constant)
    {
        INMOST_ICE_ERR("d_dy of const is 0");
    }
    else if (s.GetType() == FuncType::linear)
    {
        std::vector<std::vector<double>> tr = s.GetTrian();
        return{FuncType::constant, {s.GetParams()[1]}, tr};
    }
    else if (s.GetType() == FuncType::quadratic)
    {
        INMOST_ICE_ERR("can't d/dy quadratic function");
    }
    else
    {
        INMOST_ICE_ERR("max degree for d_dy is 2");
    }
}

ScalarFunction operator*(const ScalarFunction& s1, const ScalarFunction& s2)
{
    std::vector<std::vector<double>> tr = s1.GetTrian();
    if (s1.GetType() == FuncType::constant)
    {
        double k = s1.GetParams()[0];
        if (s2.GetType() == FuncType::constant)
        {
            return {FuncType::constant, {k*s2.GetParams()[0]}, tr};
        }
        else if (s2.GetType() == FuncType::linear)
        {
            double a = s2.GetParams()[0];
            double b = s2.GetParams()[1];
            double c = s2.GetParams()[2];

            return {FuncType::linear, {k*a, k*b, k*c}, tr};
        }
        else if (s2.GetType() == FuncType::quadratic)
        {
            double a0 = s2.GetParams()[0];
            double a1 = s2.GetParams()[1];
            double a2 = s2.GetParams()[2];
            double a3 = s2.GetParams()[3];
            double a4 = s2.GetParams()[4];
            double a5 = s2.GetParams()[5];

            return {FuncType::quadratic, {k*a0, k*a1, k*a2, k*a3, k*a4, k*a5}, tr};
        }
        else
        {
            INMOST_ICE_ERR("max result degree for mult is 2");
        }
    }
    else if (s1.GetType() == FuncType::linear)
    {
        double a0 = s1.GetParams()[0];
        double a1 = s1.GetParams()[1];
        double a2 = s1.GetParams()[2];
            
        if (s2.GetType() == FuncType::constant)
        {
            double b = s2.GetParams()[0];
            return {FuncType::linear, {a0*b, a1*b, a2*b}, tr};
        }
        else if (s2.GetType() == FuncType::linear)
        {
            double b0 = s2.GetParams()[0];
            double b1 = s2.GetParams()[1];
            double b2 = s2.GetParams()[2];

            return {FuncType::quadratic, {a0*b0, a0*b1 + a1*b0, a1*b1,
                                       a0*b2 + b0*a2, a1*b2 + b1*a2, a2*b2}, tr};
        }
        else
        {
            INMOST_ICE_ERR("max result degree for mult is 2");
        }
    }
    else if (s1.GetType() == FuncType::quadratic)
    {
        double a0 = s1.GetParams()[0];
        double a1 = s1.GetParams()[1];
        double a2 = s1.GetParams()[2];
        double a3 = s1.GetParams()[3];
        double a4 = s1.GetParams()[4];
        double a5 = s1.GetParams()[5];
        if (s2.GetType() == FuncType::constant)
        {
            double b = s2.GetParams()[0];
            return {FuncType::quadratic, {a0*b, a1*b, a2*b, a3*b, a4*b, a5*b}, tr};
        }
    }
    else
    {
        INMOST_ICE_ERR("max result degree for mult is 2");
    }
}

VectorFunction operator*(const ScalarFunction& s, const VectorFunction& v)
{
    std::vector<std::vector<double>> tr = s.GetTrian();
    if (s.GetType() == FuncType::constant)
    {
        double a = s.GetParams()[0];
        if (v.GetType() == FuncType::constant)
        {
            double b = v.GetParams()[0][0];
            double c = v.GetParams()[1][0];

            return {FuncType::constant, {a*b}, {a*c}, tr};
        }
        else if (v.GetType() == FuncType::linear)
        {
            double b0 = v.GetParams()[0][0];
            double b1 = v.GetParams()[0][1];
            double b2 = v.GetParams()[0][2];

            double c0 = v.GetParams()[1][0];
            double c1 = v.GetParams()[1][1];
            double c2 = v.GetParams()[1][2];

            return {FuncType::linear, {a*b0, a*b1, a*b2},
                                      {a*c0, a*c1, a*c2}, tr};
        }
        else if (v.GetType() == FuncType::quadratic)
        {
            double b0 = v.GetParams()[0][0];
            double b1 = v.GetParams()[0][1];
            double b2 = v.GetParams()[0][2];
            double b3 = v.GetParams()[0][3];
            double b4 = v.GetParams()[0][4];
            double b5 = v.GetParams()[0][5];

            double c0 = v.GetParams()[1][0];
            double c1 = v.GetParams()[1][1];
            double c2 = v.GetParams()[1][2];
            double c3 = v.GetParams()[1][3];
            double c4 = v.GetParams()[1][4];
            double c5 = v.GetParams()[1][5];

            return {FuncType::quadratic, {a*b0, a*b1, a*b2, a*b3, a*b4, a*b5},
                                         {a*c0, a*c1, a*c2, a*c3, a*c4, a*c5}, tr};
        }
        else
        {
            INMOST_ICE_ERR("max result degree for mult is 2");
        }
    }
    else if (s.GetType() == FuncType::linear)
    {
        double a0 = s.GetParams()[0];
        double a1 = s.GetParams()[1];
        double a2 = s.GetParams()[2];
            
        if (v.GetType() == FuncType::constant)
        {
            double b = v.GetParams()[0][0];
            double c = v.GetParams()[1][0];
            return {FuncType::linear, {a0*b, a1*b, a2*b},
                                      {a0*c, a1*c, a2*c}, tr};
        }
        else if (v.GetType() == FuncType::linear)
        {
            double b0 = v.GetParams()[0][0];
            double b1 = v.GetParams()[0][1];
            double b2 = v.GetParams()[0][2];

            double c0 = v.GetParams()[1][0];
            double c1 = v.GetParams()[1][1];
            double c2 = v.GetParams()[1][2];

            return {FuncType::quadratic, {a0*b0, a0*b1 + b0*a1, a1*b1, a0*b2 + b0*a2, a1*b2 + b1*a2, a2*b2}, 
                                         {a0*c0, a0*c1 + c0*a1, a1*c1, a0*c2 + c0*a2, a1*c2 + c1*a2, a2*c2}, tr};
        }
        else
        {
            INMOST_ICE_ERR("max result degree for mult is 2");
        }
    }
    else if (s.GetType() == FuncType::quadratic)
    {
        double a0 = s.GetParams()[0];
        double a1 = s.GetParams()[1];
        double a2 = s.GetParams()[2];
        double a3 = s.GetParams()[3];
        double a4 = s.GetParams()[4];
        double a5 = s.GetParams()[5];
        if (v.GetType() == FuncType::constant)
        {
            double b = v.GetParams()[0][0];
            double c = v.GetParams()[1][0];
            return {FuncType::quadratic, {a0*b, a1*b, a2*b, a3*b, a4*b, a5*b},
                                         {a0*c, a1*c, a2*c, a3*c, a4*c, a5*c}, tr};
        }
    }
    else
    {
        INMOST_ICE_ERR("max result degree for mult is 2");
    }
}

ScalarFunction operator+(const ScalarFunction& s1, const ScalarFunction& s2)
{
    std::vector<std::vector<double>> tr = s1.GetTrian();
    if (s1.GetType() == FuncType::constant)
    {
        double a = s1.GetParams()[0];
        if (s2.GetType() == FuncType::constant)
        {
            double b = s2.GetParams()[0];
            return {FuncType::constant, {a+b}, tr};
        }
        else
        {
            INMOST_ICE_ERR("no sense to sum different degree function");
        }
    }
    else if (s1.GetType() == FuncType::linear)
    {
        double a0 = s1.GetParams()[0];
        double a1 = s1.GetParams()[1];
        double a2 = s1.GetParams()[2];
            
        if (s2.GetType() == FuncType::linear)
        {
            double b0 = s2.GetParams()[0];
            double b1 = s2.GetParams()[1];
            double b2 = s2.GetParams()[2];
            return {FuncType::linear, {a0 + b0, a1 + b1, a2 + b2}, tr};
        }
        else
        {
            INMOST_ICE_ERR("no sense to sum different degree function");
        }
    }
    else if (s1.GetType() == FuncType::quadratic)
    {
        double a0 = s1.GetParams()[0];
        double a1 = s1.GetParams()[1];
        double a2 = s1.GetParams()[2];
        double a3 = s1.GetParams()[3];
        double a4 = s1.GetParams()[4];
        double a5 = s1.GetParams()[5];
        if (s2.GetType() == FuncType::quadratic)
        {
            double b0 = s2.GetParams()[0];
            double b1 = s2.GetParams()[1];
            double b2 = s2.GetParams()[2];
            double b3 = s2.GetParams()[3];
            double b4 = s2.GetParams()[4];
            double b5 = s2.GetParams()[5];

            return {FuncType::quadratic, {a0 + b0, a1 + b1, a2 + b2, 
                                          a3 + b3, a4 + b4, a5 + b5}, tr};
        }
        else
        {
            INMOST_ICE_ERR("no sense to sum different degree function");
        }
        
    }
    else
    {
        INMOST_ICE_ERR("unknown type of function for summation");
    }
}