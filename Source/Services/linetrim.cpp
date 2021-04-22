#include "inmost_ice.h"

// trim from start (in place)
void ltrim(std::string& s) 
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), 
    [](unsigned char ch) 
    {
        return !std::isspace(ch);
    }));
};

// trim from end (in place)
void rtrim(std::string& s) 
{
    s.erase(std::find_if(s.rbegin(), s.rend(), 
    [](unsigned char ch) 
    {
        return !std::isspace(ch);
    }).base(), s.end());
};

// trim from both ends (in place)
void trim(std::string& s) 
{
    ltrim(s);
    rtrim(s);
};
