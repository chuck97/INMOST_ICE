#include "netcdf_INMOST.h"

int main()
{
    int retval;
    Dataset d(FILE_NAME);
    d.ReadVars(TIME_STEP);
    d.VarsInfo();
    //d.ReadDims();
    //d.ReadAtts();
    //d.AttsInfo();
    //d.DimsInfo();
}