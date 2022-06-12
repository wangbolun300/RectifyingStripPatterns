#pragma once
#include <lsc/MeshProcessing.h>

// Efunc represent a elementary value, which is the linear combination of
// some function values on their corresponding vertices, of this vertex.
// i.e. Efunc[0]*f_0+Efunc[1]*f_1+...+Efunc[n]*f_n.
typedef Eigen::SparseVector<double> Efunc;
// // this class is the wrapper of the linear combination
// class lsElem{
//  public:
//  Eigen::SparseVector<double> ele;
// };

// vector of the linear combination of function values (LSC).
class Vectorlf
{
private:
    std::vector<Efunc> mat; // the vector of LCF, for each vertex there is a LCF

public:
    Vectorlf(){};
    int size()
    {
        return mat.size();
    }

    // since for each mesh vertex there is a LSC, and for each LSC, it is the combination of all function values of
    // the vertices, which can be represented as a vector of size sz (some cofficients are 0, so we use sparse matrix)
    void resize(int sz)
    {
        mat.resize(sz);
        for (int i = 0; i < sz; i++)
        {
            mat[i].resize(sz);
        }
    }

    Efunc &operator()(int id)
    {
        assert(id >= 0 && id < mat.size());
        return mat[id];
    }
};
// The basic tool of LSC. Please initialize it with a mesh
class lsTools
{
private:
    MeshProcessing MP;
    CGMesh lsmesh;                       // the input mesh
    Eigen::MatrixXd norm_f;              // normal per face
    Eigen::MatrixXd norm_v;              // normal perf vertex
    Eigen::MatrixXd angF;                // angel per vertex of each F. nx3, n is the number of faces
    Eigen::VectorXd areaF;               // the area of each face.
    std::vector<Eigen::Matrix3d> Rotate; // the 90 degree rotation matrices for each face.

public:
    lsTools(CGMesh &mesh);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi bnd;   // boundary loop
    Eigen::MatrixXd paras; // parameters of the mesh vertices, nx2

    // convert the parameters into a mesh
    void convert_paras_as_meshes(CGMesh &output);
    void get_mesh_angles();
    void get_mesh_normals_per_face();
    void get_mesh_normals_per_ver();
    void get_face_rotation_matices();
    // template <typename Tp>
    // void gradient(Tp&)
};