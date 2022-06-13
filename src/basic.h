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

// // vector of 
// class VecScalar{

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

    
    void resize(int sz)
    {
        mat.resize(sz);
    }
    void resize(int sz1, int sz2)
    {
        mat.resize(sz1);
        for (int i = 0; i < sz1; i++)
        {
            mat[i].resize(sz2);
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
    // the rotated 3 half edges for each face, in the plane of this face
    // the order of the 3 edges are: v0-v2, v1-v0, v2-v1
    std::vector<std::array<Eigen::Vector3d, 3>> Erotate; 
    Vectorlf gradVF; //gradient of function in each face
    Vectorlf gradV; //gradient of function in each vertex
    void get_mesh_angles();
    void get_mesh_normals_per_face();
    void get_mesh_normals_per_ver();

    // This is to calculate the 90 degree rotation matrices in
    // each triangle face. It can be used to calculate gradients.
    void get_face_rotation_matices();
    void get_rotated_edges_for_each_face();
    template <class Cls>
    void gradient_v2f( Cls& values,std::array<Cls,3> &output);// calculate gradient in each face from vertex values 
    template <class Cls>
    void gradient_f2v( std::array<Cls, 3>& values,std::array<Cls,3> &output);// calculate gradient in each vertex by averging face values 
public:
    // parametrization and find the boundary loop 
    lsTools(CGMesh &mesh);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXi E;
    Eigen::VectorXi bnd;   // boundary loop
    Eigen::MatrixXd paras; // parameters of the mesh vertices, nx2
    
    // this function should be calculated first once the class get constructed
    //  1. get face normals; 
    //  2. get all the angles;
    //  3. get vertex normals using angles and face normals;
    //  4. get the face rotation matrices.
    //  5. get the rotated edges in each face.
    void initialize_mesh_properties(){
        get_mesh_normals_per_face();
        get_mesh_angles();
        get_mesh_normals_per_ver();
        get_face_rotation_matices();
        get_rotated_edges_for_each_face();
    }
     // convert the parameters into a mesh, for visulization purpose
    void convert_paras_as_meshes(CGMesh &output);

    void debug_tool(int id=0, double value=0){
        // std::cout<<"the rotation matrix of face id "<<id<<":\n"<<Rotate[id]<<"\n";
        // std::cout<<"points 0 and 1:\np0: "<<V.row(F(id,0))<<", \np1: "<<V.row(F(id,1))<<""<<std::endl;
        // Eigen::Vector3d dirc=Eigen::Vector3d(V.row(F(id,1))-V.row(F(id,0)));
        // std::cout<<"direction v01: "<<dirc<<std::endl;

        // Eigen::Vector3d rotated=Rotate[id]*dirc;
        // std::cout<<"Rotated "<<rotated<<std::endl;
        // std::cout<<"the inner product: "<<dirc.dot(rotated)<<std::endl;
        // std::cout<<"the difference of the lengths "<<rotated.norm()-dirc.norm()<<std::endl;
        Efunc vec, vec1;
        vec.resize(3);
        vec1.resize(3);
        std::cout<<"print the test sparse vector\n"<<vec<<std::endl;
        std::cout<<"size "<<vec.size()<<std::endl;
        std::cout<<"values "<<vec.coeffRef(0)<<std::endl;
        vec.coeffRef(0)=1.0;
        vec.coeffRef(2)=2.0;
        vec=vec*2;
        std::cout<<"print the test sparse vector\n"<<vec<<std::endl;
        std::cout<<"size "<<vec.size()<<std::endl;
        std::cout<<"values "<<vec.coeffRef(0)<<std::endl;
    }

    
};