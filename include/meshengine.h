#pragma once

/* Handles everything on mesh from reading .mphtxt file */

#define _HAS_STD_BYTE 0
#include<string>
#include<fstream>
#include <algorithm>
#include <vector>
#include <functional> 
#include <cctype>
#include <locale>
#include <filesystem>
#include "utilities.h"

#include <map>


using namespace std;
struct Mesh2D{
        dmatrix <double> nodes;
        dmatrix <int> triangles;
        dmatrix <int> tritype;
        dmatrix <int> tri2edges;
        dmatrix <int> pecnodes;
        dmatrix <int> pecedges;
        dmatrix <int> signs;
        int nedges;

    };

class MeshEngine{

private:
	Util ut;
    struct TagSet { 
        string name;
        vector <int> tags;
        int tagId;
    }; 

    string meshfilename;

    //Number of nodes, edges, tetrahedrons and triangles


    int domcount, boundcount;
    vector<TagSet> Boundaries;
    vector<TagSet> Domains;
    map<int, int> boundDictionary;
    map<int, int> domDictionary;

public: 
    int nnodes;
    int nedges;
    int ntets;
    int ntris;
    
    dmatrix<double> nodes;
    dmatrix<int> tetrahedrons;
    dmatrix<int> triangles;
    int **tet2edge;
    dmatrix<int> tri2edge;
 
    int *tettype;
    dmatrix <int> tritype;
    int *tetoftris;
    int **signs_tet;
    int **signs_tri;

    MeshEngine();
    ~MeshEngine();

    Mesh2D getLocal2DMesh(int fno);


    void displayTagInfo();
    void dsplayMeshInfo();
    void  readMesh(string str);


    dmatrix<int> getTriangles(int t);
    dmatrix<int> getTetrahedrons(int t);


	void getVolume() ;
private:
    int getTetofTri();
    void signs();
    void getEdgeData();
    inline int getId(int a, int b);
    void getMeshInfo();
    void getTagInfo();


};