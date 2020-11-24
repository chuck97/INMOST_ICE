#include<any>
#include<vector>
#include<iostream>

auto Vec()
{
    std::vector<int> a = {1};
    return a;
}

int main()
{
    auto b = Vec();
    std::cout << b[0] <<std::endl;
    std::any a = 1;
    std::cout << a.type().name() << ": " << std::any_cast<int>(a) << std::endl;
    a = 3.14;
    std::cout << a.type().name() << ": " << std::any_cast<double>(a) << std::endl;
    a = true;
    std::cout << a.type().name() << ": " << std::any_cast<bool>(a) << std::endl;
}