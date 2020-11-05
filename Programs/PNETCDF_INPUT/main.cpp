#include <stdio.h>
#include <stdlib.h>

#include <iostream>
using namespace std;

#include <string.h> /* strcpy(), strncpy() */
#include <unistd.h> /* getopt() */
#include <pnetcdf>

using namespace PnetCDF;
using namespace PnetCDF::exceptions;

static void
usage(char *argv0)
{
    cerr <<
    "Usage: %s [-h] | [-q] [file_name]\n"
    "       [-h] Print help\n"
    "       [-q] Quiet mode (reports when fail)\n"
    "       [filename] input netCDF file name\n"
    << argv0;
}

int main(int argc, char** argv)
{
    extern int optind;
    // .nc filename
    char filename[256] = "/data90t/geosci/spetrov/nc_topaz/dataset-ice.nc";
    char str_att[NC_MAX_NAME];
    int i, rank, nprocs, err, verbose=1;
    MPI_Offset len, global_ny, global_nx, local_ny, local_nx;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    try 
    {
        //open file
        NcmpiFile ncFile(MPI_COMM_WORLD, filename, NcmpiFile::read);

        //check file format
        if (ncFile.getFormat() != NcmpiFile::classic) 
        {
            cout << "unexpected file format" << endl;
            throw NcmpiException("read Error ",__FILE__,__LINE__);
        }

        // get global attribute named "CDI" 
        NcmpiGroupAtt att = ncFile.getAtt("CDI");
        att.getValues(str_att);
        len = att.getAttLength();
        str_att[len] = '\0'; /* add a NULL char at the end */
        if (rank == 0)
        {
            std::cout << "global attribute \"CDI\" :" << str_att << std::endl;
        }

        // get dimensio  IDs
        NcmpiDim xDim = ncFile.getDim("x");
        NcmpiDim yDim = ncFile.getDim("y");

        // get dimension lengths 
        global_nx = xDim.getSize();
        global_ny = yDim.getSize();

        if (rank == 0)
        {
            cout << "xdim: " << global_nx << std::endl;
            cout << "ydim: " << global_ny << std::endl;
        }

        // get the variable ID of a 2D variable of double type "fice" 
        NcmpiVar var = ncFile.getVar("fice");

        // get variable's attribute named "standard_name"
        NcmpiVarAtt vatt = var.getAtt("standard_name");
        vatt.getValues(str_att);
        len = vatt.getAttLength();
        str_att[len] = '\0'; /* add a NULL char at the end */
        if (rank == 0)
        {
            std::cout << "variable \"fice\" attribute \"standart_name\" :" << str_att << std::endl;
        }


        //assign local sizes
        local_nx = global_nx/nprocs;
        local_ny = global_ny;

        // get "fice" data

        vector<MPI_Offset> start(3), count(3);
        start[0] = 0;
        start[1] = local_nx*rank;
        start[3] = 0;

        count[0] = local_ny;
        count[1] = local_nx;
        count[3] = 1;

        double* buf = new double[local_nx*local_ny];

        var.getVar_all(start, count, buf);

        delete[] buf;

    }
    catch(NcmpiException& e) 
    {
       cout << e.what() << " error code=" << e.errorCode() << " Error!\n";
       return 1;
    }

    MPI_Finalize();
    return 0;
}
