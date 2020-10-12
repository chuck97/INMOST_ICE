#include <utility>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include "inmost.h"

#define USE_MPI
#define USE_PARTITIONER
using namespace INMOST;
#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD);
#else
#define BARRIER
#endif

#define M_Assert(Expr, Msg) \
    if (!Expr)  \
    { \
    std::cerr << "Assert failed:\t" << Msg << "\n" \
        << "Expected:\t" << #Expr << "\n" \
        << "Source:\t\t" << __FILE__ << ", line " << __LINE__ << "\n"; \
        abort(); \
    }



int main(int argc, char **argv) 
{
    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);

    Mesh *m;
    m = new Mesh();
    const char * name_vrt = "/data90t/geosci/spetrov/INMOST_ICE/Triangulation/Triangulation_thickening/triangulation_output/vrt.txt";
    const char * name_bnd = "/data90t/geosci/spetrov/INMOST_ICE/Triangulation/Triangulation_thickening/triangulation_output/bnd.txt";
    const char * name_tri = "/data90t/geosci/spetrov/INMOST_ICE/Triangulation/Triangulation_thickening/triangulation_output/tri.txt";
    std::ifstream file;

    Tag Label_tag = m->CreateTag("LabelTag", DATA_INTEGER, CELL | FACE | EDGE | NODE, NONE, 1);
    Tag Node_id = m->CreateTag("Node_id", DATA_INTEGER, NODE, NONE, 1);
    Tag Tri_id = m->CreateTag("Tri_id", DATA_INTEGER, CELL, NONE, 1);
    Tag Bnd_id = m->CreateTag("Bnd_id", DATA_INTEGER, FACE, FACE, 1);
    Tag Is_node_bnd = m->CreateTag("Is_node_bnd", DATA_INTEGER, NODE, NODE, 1);

    std::cout << "Start to read mesh!" <<std::endl;

    // Read verticies
    file.open(name_vrt);
    M_Assert((file.is_open()), "can not open vrt file");

    std::string tmp_STR;
    int nv,nt,nf;
    file >> tmp_STR;
    file >> tmp_STR;
    file >> tmp_STR;
    file >> nv;
    ElementArray<Node> newverts(m);
    newverts.reserve(nv);
    std::cout << "num of verts: " << nv << std::endl;
    for (int i = 0; i < nv; i++) 
    {
        Storage::real xyz[3];

        file>>xyz[0];
        file>>xyz[1];
        xyz[2] = 0.0;

        Node nod = m->CreateNode(xyz);
        nod->Integer(Node_id) = i+1;
        nod->Integer(Label_tag) = 1;
        nod->Integer(Is_node_bnd) = 0;
        newverts.push_back(nod);
    }
    file.close();

    // Read triangles
    file.open(name_tri);
    M_Assert((file.is_open()), "can not open tri file");
    file >> tmp_STR;
    file >> tmp_STR;
    file >> tmp_STR;
    file >> nt;
    std::cout << "num of triagles: " << nt << std::endl;
    int num_node_tmp;
    for (int i = 0; i < nt; i++) 
    {
        ElementArray<Node> verts(m);
        for (int j = 0; j < 3; j++) 
        {
            file >> num_node_tmp;
            verts.push_back(newverts[num_node_tmp - 1]);
        }

        const INMOST_DATA_INTEGER_TYPE ne_face_nodes[6] = {0, 1, 1, 2, 2, 0};
        const INMOST_DATA_INTEGER_TYPE ne_num_nodes[3] = {2, 2, 2};

        std::pair<Cell, bool> pair = m -> CreateCell(verts, ne_face_nodes, ne_num_nodes, 3);

        pair.first.Integer(Tri_id) = i+1;

        pair.first.Integer(Label_tag) = 1;
    }
    file.close();

    // Read boundary edges
    file.open(name_bnd);
    M_Assert((file.is_open()),"can not open bnd file");
    file >> tmp_STR;
    file >> tmp_STR;
    file >> tmp_STR;
    file >> tmp_STR;
    file >> nf;
    std::cout << "num of edges: " << nf << std::endl;
    int nums[2];
    for(int i = 0; i < nf; ++i)
    {
        ElementArray<Node> verts(m);
        for (int j = 0; j < 2; j++) 
        {
            file >> nums[j];
            newverts[nums[j] + 1].Integer(Is_node_bnd) = 1; 
            verts.push_back(newverts[nums[j]-1]);
        }

        std::pair<Face, bool> pair = m -> CreateFace(verts);
        M_Assert((!pair.second), "we didn't find face!!!");
        pair.first.Integer(Bnd_id) = i+1;
    }
    ////////////////////////////////////// load and repartition mesh //////////////////////////
    m->Save("test.pmf");
    m->Save("test.vtk");
    //////////////////////////////////create discretization/////////////////////////////

    return 0;
}